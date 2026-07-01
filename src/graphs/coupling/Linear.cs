using System;

namespace Graphs.Coupling;

/// <summary>
/// Math evaluator for the Linear (triangular) coupling kernel.
/// </summary>
public static class LinearKernel
{
    /// <summary>
    /// J = max(0, 1 - d / δ). Compact support — exactly zero beyond one scale length.
    /// Hard locality boundary; no long-range coupling whatsoever.
    /// </summary>
    public static double Evaluate(double distance, double delta)
        => Math.Max(0.0, 1.0 - distance / delta);
}

/// <summary>
/// Linear (triangular) kernel descriptor for the fluent compiler API.
/// Pass via <c>GraphCompiler.CouplingStrategy(new Linear(0.5))</c> or
/// <c>Kernels.Linear(0.5)</c>. A <see cref="Bandwidth"/> of 0.0 is the
/// auto-estimate sentinel.
/// </summary>
public readonly record struct Linear(double Bandwidth = 0.0) : IKernelDescriptor;
