using System;
using Maths.Rng;
using Maths.Samplers.Blocks;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// The full conditional of an anisotropic penalty's smoothing variances, on <c>ρⱼ = ln τⱼ²</c>, and its sampler. The
/// coefficients' prior is <c>N(0, K⁻)</c> with <c>K = Σⱼ Pⱼ/τⱼ²</c>, and each τⱼ² has the prior q
/// (<see cref="ISmoothingVariancePrior"/>: the inverse gamma by default, or the Weibull). With <c>wⱼ = 1/τⱼ²</c> and
/// <c>Rⱼ = βᵀPⱼβ</c>, <c>ln p(ρ | β) = ½ ln Det K − ½ Σⱼ Rⱼwⱼ + Σⱼ ln q(ρⱼ) + const</c> (Bach &amp; Klein 2025, eq. 12). Its
/// gradient and Hessian follow from the spectrum's derivatives in <c>ln w = −ρ</c> (<see cref="PenaltySpectrum"/>) and
/// the prior's. With one component and the IG prior it is the IG conditional of the isotropic Gibbs scheme.
/// <para>
/// Each update is a Newton-proposal Metropolis step on every ρⱼ at once (<see cref="NewtonMetropolis"/>), Bach &amp;
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

    private readonly IAnisotropicPenalty _penalty;
    private readonly ISmoothingVariancePrior _prior;
    private readonly NewtonMetropolis _step;
    private readonly NewtonMetropolis.LogDensity _density;
    private readonly Func<double, double>[] _coordinate;   // ρ ↦ ln p with only ρⱼ replaced, for the random-walk steps
    private readonly double[] _roughness;
    private readonly double[] _weights;
    private readonly double[] _gradient;   // of ln Det K in ln w
    private readonly double[] _hessian;
    private double[] _current = Array.Empty<double>();

    public ComponentSmoothing(IAnisotropicPenalty penalty, ISmoothingVariancePrior prior)
    {
        _penalty = penalty;
        _prior = prior;
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

    /// <summary>One update of <paramref name="logTau2"/> given the coefficients: the Newton step, then the random-walk
    /// steps. Returns whether the Newton step was accepted.</summary>
    public bool Update(double[] logTau2, double[] beta, Xoshiro256PlusPlus rng)
    {
        for (int j = 0; j < _roughness.Length; j++) _roughness[j] = _penalty.ComponentRoughness(beta, j);
        bool accepted = _step.Step(logTau2, _density, rng);
        _current = logTau2;
        double value = LogDensity(logTau2);
        for (int j = 0; j < logTau2.Length; j++)
            RandomWalkMetropolis.Step(ref logTau2[j], ref value, EscapeStep, _coordinate[j], rng);
        return accepted;
    }

    /// <summary>The log full conditional at ρ.</summary>
    public double LogDensity(ReadOnlySpan<double> rho)
    {
        if (!Weights(rho)) return double.NegativeInfinity;
        double value = 0.5 * _penalty.Spectrum.LogPseudoDeterminant(_weights);
        for (int j = 0; j < rho.Length; j++) value += _prior.LogDensityOfLog(rho[j], out _, out _) - 0.5 * _roughness[j] * _weights[j];
        return value;
    }

    /// <summary>The log full conditional at ρ, with its gradient and Hessian in ρ.</summary>
    public double LogDensity(ReadOnlySpan<double> rho, Span<double> gradient, Span<double> hessian)
    {
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
