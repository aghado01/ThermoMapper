using System;
using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// A <see cref="LocalBetaKernel"/> whose spread τ is tuned online toward a target acceptance rate — adaptive
/// Metropolis (Haario; Roberts &amp; Rosenthal) for the free-knot relocate/birth proposal. A fixed τ can't be
/// right across datasets of different knot density; <see cref="Adapt"/> nudges <c>log τ</c> by a Robbins–Monro
/// step from the observed acceptance (higher acceptance ⇒ proposals too timid ⇒ loosen, i.e. <i>smaller</i> τ).
/// To stay ergodic, adapt only during burn-in and freeze for the sampling phase (diminishing/zero adaptation).
/// Delegates the Beta math to an immutable inner kernel, rebuilt on each adapt — no duplicated proposal code.
/// </summary>
public sealed class AdaptiveLocalBetaKernel : IKnotKernel
{
    private const double MinLogTau = 0.0;          // τ ≥ 1
    private const double MaxLogTau = 11.512925;    // τ ≤ 1e5
    private LocalBetaKernel _inner;
    private double _logTau;

    public AdaptiveLocalBetaKernel(double initialTau = 50.0)
    {
        if (!(initialTau > 0.0)) throw new ArgumentOutOfRangeException(nameof(initialTau), "Spread τ must be positive.");
        _logTau = Math.Clamp(Math.Log(initialTau), MinLogTau, MaxLogTau);
        _inner = new LocalBetaKernel(Math.Exp(_logTau));
    }

    /// <summary>The current spread τ.</summary>
    public double Tau => Math.Exp(_logTau);

    public double Sample(double center, Xoshiro256PlusPlus rng) => _inner.Sample(center, rng);

    public double LogDensity(double x, double center) => _inner.LogDensity(x, center);

    /// <summary>
    /// One Robbins–Monro tuning step toward <paramref name="target"/> acceptance:
    /// <c>log τ −= stepSize·(acceptanceRate − target)</c> (clamped to a sane τ range), then rebuild the kernel.
    /// Use a vanishing <paramref name="stepSize"/> (or freeze after burn-in) to preserve ergodicity. Default
    /// target 0.44 is the 1-D random-walk optimum (the relocate move perturbs a single knot).
    /// </summary>
    public void Adapt(double acceptanceRate, double stepSize, double target = 0.44)
    {
        if (!(stepSize >= 0.0)) throw new ArgumentOutOfRangeException(nameof(stepSize));
        _logTau = Math.Clamp(_logTau - stepSize * (acceptanceRate - target), MinLogTau, MaxLogTau);
        _inner = new LocalBetaKernel(Math.Exp(_logTau));
    }
}
