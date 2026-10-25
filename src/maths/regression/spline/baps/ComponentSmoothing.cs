using System;
using Maths.Rng;
using Maths.Samplers.Blocks;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// How the anisotropic schemes update their log smoothing variances (<see cref="ComponentSmoothing"/>). Both leave the full
/// conditional invariant. Slice is the default: in the pre-registered comparison of 2026-10-24
/// (<c>SmoothingUpdateBenchmark</c>) it gave 1.28–1.37 times Newton's worst-component ESS per second on four BAPS regimes
/// and 1.05 times on a hierarchical one.
/// </summary>
public enum SmoothingVarianceUpdate
{
    /// <summary>One slice-sampling update on each ln τⱼ² in turn (<see cref="SliceSampler"/>, width 1, at most 32
    /// steps out): values only, no derivatives, and no rejections. The default.</summary>
    Slice,

    /// <summary>Bach &amp; Klein's Newton-proposal Metropolis step on every ln τⱼ² at once (<see cref="NewtonMetropolis"/>),
    /// then one random-walk Metropolis step on each, which reach the tails the Newton step alone visits too
    /// rarely.</summary>
    Newton,
}

/// <summary>
/// The full conditional of an anisotropic penalty's smoothing variances, on <c>ρⱼ = ln τⱼ²</c>, and its sampler. The
/// coefficients' prior is <c>N(0, K⁻)</c> with <c>K = Σⱼ Pⱼ/τⱼ²</c>, and each τⱼ² has the prior q
/// (<see cref="ISmoothingVariancePrior"/>: the inverse gamma by default, or the Weibull). With <c>wⱼ = 1/τⱼ²</c> and
/// <c>Rⱼ = βᵀPⱼβ</c>, <c>ln p(ρ | β) = ½ ln Det K − ½ Σⱼ Rⱼwⱼ + Σⱼ ln q(ρⱼ) + const</c> (Bach &amp; Klein 2025, eq. 12). Its
/// gradient and Hessian follow from the spectrum's derivatives in <c>ln w = −ρ</c> (<see cref="PenaltySpectrum"/>) and
/// the prior's. With one component and the IG prior it is the IG conditional of the isotropic Gibbs scheme.
/// <para>
/// The update is <see cref="SmoothingVarianceUpdate.Newton"/> or <see cref="SmoothingVarianceUpdate.Slice"/>. The Newton
/// update is a Newton-proposal Metropolis step on every ρⱼ at once (<see cref="NewtonMetropolis"/>), Bach &amp;
/// Klein's update, then one random-walk Metropolis step on each ρⱼ in turn with standard deviation
/// <see cref="EscapeStep"/>. Both leave the conditional invariant. The random-walk steps are there for the tails: ρⱼ's
/// upper tail decays only exponentially, more slowly than the Gaussian proposal, and from deep in it the capped Newton
/// step lands far beyond the mode, so a chain of Newton steps alone stays there for thousands of steps and visits the
/// tail too rarely. Measured on a fixed β: 4·10⁶ Newton steps never reached a region the conditional gives
/// probability 4.4·10⁻⁴, with rejection runs of 10⁴ steps; with the random-walk steps they hit it at 4.47·10⁻⁴, and
/// the longest rejection run was 8 steps.
/// </para>
/// </summary>
internal sealed class ComponentSmoothing
{
    /// <summary>The random-walk steps' standard deviation, in ln τ².</summary>
    public const double EscapeStep = 1.0;

    /// <summary>The slice update's initial width, in ln τ², and its stepping-out limit.</summary>
    public const double SliceWidth = 1.0;
    public const int SliceSteps = 32;

    private readonly IAnisotropicPenalty _penalty;
    private readonly ISmoothingVariancePrior _prior;
    private readonly SmoothingVarianceUpdate _update;
    private readonly NewtonMetropolis _step;
    private readonly NewtonMetropolis.LogDensity _density;
    private readonly Func<double, double>[] _coordinate;   // ρ ↦ ln p with only ρⱼ replaced, for the random-walk steps
    private readonly double[] _roughness;
    private readonly double[] _weights;
    private readonly double[] _gradient;   // of ln Det K in ln w
    private readonly double[] _hessian;
    private double[] _current = Array.Empty<double>();

    public ComponentSmoothing(IAnisotropicPenalty penalty, ISmoothingVariancePrior prior,
                              SmoothingVarianceUpdate update = SmoothingVarianceUpdate.Slice)
    {
        if (update is not (SmoothingVarianceUpdate.Newton or SmoothingVarianceUpdate.Slice)) throw new ArgumentOutOfRangeException(nameof(update));
        _penalty = penalty;
        _prior = prior;
        _update = update;
        int p = penalty.Components;
        _step = new NewtonMetropolis(p);
        _density = LogDensity;
        _roughness = new double[p];
        _weights = new double[p];
        _gradient = new double[p];
        _hessian = new double[p * p];
        _coordinate = new Func<double, double>[p];
        for (int j = 0; j < p; j++)
        {
            int jj = j;
            _coordinate[j] = u =>
            {
                double keep = _current[jj];
                _current[jj] = u;
                double value = LogDensity(_current);
                _current[jj] = keep;
                return value;
            };
        }
    }

    /// <summary>The log-density evaluations spent so far, with or without derivatives.</summary>
    public long Evaluations { get; private set; }

    /// <summary>One update of <paramref name="logTau2"/> given the coefficients. The Newton update returns whether its
    /// Newton step was accepted; the slice update always moves.</summary>
    public bool Update(double[] logTau2, double[] beta, Xoshiro256PlusPlus rng)
    {
        for (int j = 0; j < _roughness.Length; j++) _roughness[j] = _penalty.ComponentRoughness(beta, j);
        bool accepted = _update == SmoothingVarianceUpdate.Slice || _step.Step(logTau2, _density, rng);
        _current = logTau2;
        double value = LogDensity(logTau2);
        for (int j = 0; j < logTau2.Length; j++)
        {
            if (_update == SmoothingVarianceUpdate.Slice)
                SliceSampler.Step(ref logTau2[j], ref value, SliceWidth, SliceSteps, _coordinate[j], rng);
            else
                RandomWalkMetropolis.Step(ref logTau2[j], ref value, EscapeStep, _coordinate[j], rng);
        }
        return accepted;
    }

    /// <summary>The log full conditional at ρ.</summary>
    public double LogDensity(ReadOnlySpan<double> rho)
    {
        Evaluations++;
        if (!Weights(rho)) return double.NegativeInfinity;
        double value = 0.5 * _penalty.Spectrum.LogPseudoDeterminant(_weights);
        for (int j = 0; j < rho.Length; j++) value += _prior.LogDensityOfLog(rho[j], out _, out _) - 0.5 * _roughness[j] * _weights[j];
        return value;
    }

    /// <summary>The log full conditional at ρ, with its gradient and Hessian in ρ.</summary>
    public double LogDensity(ReadOnlySpan<double> rho, Span<double> gradient, Span<double> hessian)
    {
        Evaluations++;
        if (!Weights(rho)) return double.NegativeInfinity;
        int p = rho.Length;
        double value = 0.5 * _penalty.Spectrum.LogPseudoDeterminant(_weights, _gradient, _hessian);
        for (int j = 0; j < p; j++)
        {
            double pull = 0.5 * _roughness[j] * _weights[j];
            value += _prior.LogDensityOfLog(rho[j], out double first, out double second) - pull;
            gradient[j] = -0.5 * _gradient[j] + pull + first;
            for (int k = 0; k < p; k++) hessian[j * p + k] = 0.5 * _hessian[j * p + k] + (j == k ? second - pull : 0.0);
        }
        return value;
    }

    // wⱼ = e^(−ρⱼ); false beyond the range where every weight is positive and finite.
    private bool Weights(ReadOnlySpan<double> rho)
    {
        for (int j = 0; j < rho.Length; j++)
        {
            _weights[j] = Math.Exp(-rho[j]);
            if (!(_weights[j] > 0.0) || double.IsInfinity(_weights[j])) return false;
        }
        return true;
    }
}
