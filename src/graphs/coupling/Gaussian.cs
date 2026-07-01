using System;

namespace Graphs.Coupling;

/// <summary>
/// Math evaluator for the Gaussian coupling kernel. Static — the
/// fluent API consumes the <see cref="Gaussian"/> descriptor record
/// below; the scaler dispatches here via direct call (no virtual
/// indirection on the hot path).
/// </summary>
public static class GaussianKernel
{
    /// <summary>
    /// J = exp(-d² / 2δ²). Smooth, differentiable decay.
    /// Natural choice for Euclidean-geometry metrics.
    /// </summary>
    public static double Evaluate(double distance, double delta)
    {
        double twoSigmaSq = 2.0 * delta * delta;
        return Math.Exp(-(distance * distance) / twoSigmaSq);
    }
}

/// <summary>
/// Gaussian-kernel descriptor for the fluent compiler API. Pass via
/// <c>GraphCompiler.CouplingStrategy(new Gaussian(0.5))</c> or
/// <c>Kernels.Gaussian(0.5)</c>. A <see cref="Bandwidth"/> of 0.0 is
/// the auto-estimate sentinel.
/// </summary>
public readonly record struct Gaussian(double Bandwidth = 0.0) : IKernelDescriptor;
