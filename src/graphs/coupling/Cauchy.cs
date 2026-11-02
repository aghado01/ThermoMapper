using System;

namespace Graphs.Coupling;

/// <summary>
/// Math evaluator for the Cauchy coupling kernel.
/// </summary>
public static class CauchyKernel
{
    /// <summary>
    /// J = 1 / (1 + d² / δ²). Heavy-tailed — decays much slower than Gaussian.
    /// Retains weak long-range coupling that Gaussian would zero out.
    /// </summary>
    public static double Evaluate(double distance, double delta)
        => 1.0 / (1.0 + (distance * distance) / (delta * delta));

    /// <summary>
    /// ln J = −ln(1 + d²/δ²), by log1p so a small ratio keeps its digits; past d/δ = 1e150, where d² overflows,
    /// the tail ln J = −2 ln(d/δ).
    /// </summary>
    public static double LogEvaluate(double distance, double delta)
    {
        double ratio = distance / delta;
        double squared = ratio * ratio;
        return double.IsInfinity(squared) ? -2.0 * Math.Log(ratio) : -LogMath.Log1P(squared);
    }
}

/// <summary>
/// Cauchy-kernel descriptor for the fluent compiler API. Pass via
/// <c>GraphCompiler.CouplingStrategy(new Cauchy(0.5))</c> or
/// <c>Kernels.Cauchy(0.5)</c>. A <see cref="Bandwidth"/> of 0.0 is the
/// auto-estimate sentinel.
/// </summary>
public readonly record struct Cauchy(double Bandwidth = 0.0) : IKernelDescriptor;
