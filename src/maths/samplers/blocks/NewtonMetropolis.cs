using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Samplers.Blocks;

/// <summary>
/// One Metropolis–Hastings step with a Newton proposal, also called a Taylored or IWLS proposal (Geweke &amp; Tanizaki
/// 2003). At x, with the log target's gradient g and Hessian H, it proposes <c>x′ ~ N(x + M⁻¹g, M⁻¹)</c> with
/// <c>M = −H̃</c>, where H̃ is H with its eigenvalues capped at −δ (Nocedal &amp; Wright 2006, §3.4), so the proposal is
/// defined everywhere and its steps are bounded. The step accepts on the target ratio times the reverse proposal's
/// density over the forward one's, the reverse proposal built the same way at x′. On a Gaussian target whose curvature
/// exceeds δ the proposal is the target itself and every step is accepted. Bach &amp; Klein (2025) use it for
/// anisotropic smoothing variances with δ = 1/π, the default here. Each step evaluates the target with its derivatives
/// at the current point and at the proposal, and draws the uniform whatever the proposal; a non-finite target at the
/// proposal rejects it.
/// </summary>
public sealed class NewtonMetropolis
{
    /// <summary>A log density at x, with its gradient and its Hessian (row-major) written by the call.</summary>
    public delegate double LogDensity(ReadOnlySpan<double> x, Span<double> gradient, Span<double> hessian);

    private readonly int _d;
    private readonly double _delta;
    private readonly double[] _gradient;
    private readonly double[] _hessian;
    private readonly double[] _proposal;

    /// <param name="dimension">The dimension of x.</param>
    /// <param name="delta">The cap on the Hessian's eigenvalues, −δ: the proposal's precision is at least δ in every
    /// direction.</param>
    public NewtonMetropolis(int dimension, double delta = 1.0 / Math.PI)
    {
        if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
        if (!(delta > 0.0) || double.IsInfinity(delta)) throw new ArgumentOutOfRangeException(nameof(delta), "δ must be positive and finite.");
        _d = dimension;
        _delta = delta;
        _gradient = new double[dimension];
        _hessian = new double[dimension * dimension];
        _proposal = new double[dimension];
    }

    /// <summary>The dimension of x.</summary>
    public int Dimension => _d;

    /// <summary>The cap δ on the Hessian's eigenvalues.</summary>
    public double Delta => _delta;

    /// <summary>Take one step from <paramref name="x"/>; on acceptance x moves to the proposal. Returns whether the
    /// proposal was accepted.</summary>
    public bool Step(Span<double> x, LogDensity logDensity, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(logDensity);
        ArgumentNullException.ThrowIfNull(rng);
        if (x.Length != _d) throw new ArgumentException($"x needs {_d} entries.", nameof(x));

        double current = logDensity(x, _gradient, _hessian);
        if (!double.IsFinite(current)) throw new InvalidOperationException("The log density is not finite at the current point.");
        Proposal forward = Build(x);
        for (int i = 0; i < _d; i++)
        {
            double z = Normal.SampleStandard(rng) / Math.Sqrt(forward.Precision[i]);   // along eigenvector i
            for (int k = 0; k < _d; k++) _proposal[k] = (i == 0 ? forward.Mean[k] : _proposal[k]) + forward.Basis[i][k] * z;
        }
        double logU = Math.Log(1.0 - rng.NextDouble());

        double proposed = logDensity(_proposal, _gradient, _hessian);
        if (!double.IsFinite(proposed)) return false;
        Proposal reverse = Build(_proposal);
        double logRatio = proposed - current + reverse.LogDensityAt(x) - forward.LogDensityAt(_proposal);
        if (!(logU < logRatio)) return false;
        _proposal.CopyTo(x);
        return true;
    }

    // The proposal at x from the gradient and Hessian last written: M = −H̃ = Σᵢ mᵢqᵢqᵢᵀ with mᵢ = max(−λᵢ, δ), and the
    // mean x + M⁻¹g.
    private Proposal Build(ReadOnlySpan<double> x)
    {
        var h = new double[_d, _d];
        for (int i = 0; i < _d; i++)
            for (int k = 0; k < _d; k++) h[i, k] = 0.5 * (_hessian[i * _d + k] + _hessian[k * _d + i]);
        EigenResult eigen = DenseEigen.DecomposeSymmetric(h);
        var precision = new double[_d];
        var mean = x.ToArray();
        for (int i = 0; i < _d; i++)
        {
            precision[i] = Math.Max(-eigen.Eigenvalues[i], _delta);
            double[] q = eigen.Eigenvectors[i];
            double projection = 0.0;
            for (int k = 0; k < _d; k++) projection += q[k] * _gradient[k];
            for (int k = 0; k < _d; k++) mean[k] += q[k] * projection / precision[i];
        }
        return new Proposal(mean, eigen.Eigenvectors, precision);
    }

    private sealed record Proposal(double[] Mean, double[][] Basis, double[] Precision)
    {
        // ln N(p; mean, M⁻¹) up to the constant −(d/2) ln 2π, which cancels in the ratio.
        public double LogDensityAt(ReadOnlySpan<double> p)
        {
            double sum = 0.0;
            for (int i = 0; i < Precision.Length; i++)
            {
                double c = 0.0;
                for (int k = 0; k < Mean.Length; k++) c += Basis[i][k] * (p[k] - Mean[k]);
                sum += 0.5 * Math.Log(Precision[i]) - 0.5 * Precision[i] * c * c;
            }
            return sum;
        }
    }
}
