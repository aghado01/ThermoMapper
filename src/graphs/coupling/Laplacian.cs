using System;

namespace Graphs.Coupling;

/// <summary>
/// Math evaluator for the Laplacian coupling kernel.
/// </summary>
public static class LaplacianKernel
{
    /// <summary>
    /// J = exp(-d / δ). Decays on distance, not distance². Sharper local contrast
    /// than Gaussian — near neighbors are favoured more aggressively.
    /// </summary>
    public static double Evaluate(double distance, double delta)
        => Math.Exp(-distance / delta);
}

/// <summary>
/// Laplacian-kernel descriptor for the fluent compiler API. Pass via
/// <c>GraphCompiler.CouplingStrategy(new Laplacian(0.5))</c> or
/// <c>Kernels.Laplacian(0.5)</c>. A <see cref="Bandwidth"/> of 0.0 is
/// the auto-estimate sentinel.
/// </summary>
public readonly record struct Laplacian(double Bandwidth = 0.0) : IKernelDescriptor;
