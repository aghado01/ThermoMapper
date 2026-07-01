namespace Graphs.Coupling;

/// <summary>
/// Math evaluator for the kernel-mixture coupling. Linear combination
/// of Gaussian + Cauchy + Laplacian components, each with its own
/// weight and bandwidth.
/// </summary>
public static class MixtureKernel
{
    public static double Evaluate(double distance, MixtureBandwidth b, MixtureWeights w)
        => w.Gaussian  * GaussianKernel.Evaluate(distance, b.Gaussian)
         + w.Cauchy    * CauchyKernel.Evaluate(distance, b.Cauchy)
         + w.Laplacian * LaplacianKernel.Evaluate(distance, b.Laplacian);
}

/// <summary>
/// Mixture-kernel descriptor for the fluent compiler API. Carries
/// per-component weight and per-component bandwidth in one record.
/// Pass via <c>GraphCompiler.CouplingStrategy(new Mixture(...))</c> or
/// <c>Kernels.Mixture(...)</c>.
/// </summary>
/// <remarks>
/// <para>This descriptor converts to the <see cref="MixtureWeights"/> +
/// <see cref="MixtureBandwidth"/> pair the scaler math consumes (see
/// <see cref="ToLegacy"/>). Phase 8 of the graphs-maturity refactor folds
/// <see cref="MixtureWeights"/> into this record once the REPL persistence
/// layer is migrated off it; <see cref="MixtureBandwidth"/> stays — it is the
/// resolved-bandwidth result type carried by <c>ScalerResult</c> and
/// <c>GraphBuildResult</c>, not a deletable legacy record.</para>
/// </remarks>
public readonly record struct Mixture(
    double GaussianWeight     = 0.0,
    double CauchyWeight       = 0.0,
    double LaplacianWeight    = 0.0,
    double GaussianBandwidth  = 0.0,
    double CauchyBandwidth    = 0.0,
    double LaplacianBandwidth = 0.0) : IKernelDescriptor
{
    /// <summary>Split into the legacy weight/bandwidth pair the
    /// existing scaler math consumes.</summary>
    public (MixtureWeights Weights, MixtureBandwidth? Bandwidth) ToLegacy()
    {
        var weights = new MixtureWeights(GaussianWeight, CauchyWeight, LaplacianWeight);
        bool anyExplicitBandwidth =
            GaussianBandwidth > 0.0 || CauchyBandwidth > 0.0 || LaplacianBandwidth > 0.0;
        MixtureBandwidth? bandwidth = anyExplicitBandwidth
            ? new MixtureBandwidth(GaussianBandwidth, CauchyBandwidth, LaplacianBandwidth)
            : null;
        return (weights, bandwidth);
    }
}
