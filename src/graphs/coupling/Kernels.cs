namespace Graphs.Coupling;

/// <summary>
/// Factory entry points for the kernel descriptors that the fluent
/// <c>GraphCompiler.CouplingStrategy(...)</c> consumes. Equivalent to
/// constructing the descriptors directly (<c>new Gaussian(0.5)</c>) —
/// import with <c>using static Graphs.Coupling.Kernels;</c> for the
/// unqualified <c>Gaussian(0.5)</c> syntax at call sites.
/// </summary>
public static class Kernels
{
    public static Gaussian  Gaussian(double bandwidth = 0.0)  => new(bandwidth);
    public static Cauchy    Cauchy(double bandwidth = 0.0)    => new(bandwidth);
    public static Laplacian Laplacian(double bandwidth = 0.0) => new(bandwidth);
    public static Linear    Linear(double bandwidth = 0.0)    => new(bandwidth);

    /// <summary>
    /// Construct a kernel-mixture descriptor. All weights default to
    /// 0.0; pass per-component weights (and optionally per-component
    /// bandwidths) by name. A component with weight 0 contributes
    /// nothing; a component with a non-zero weight but bandwidth 0
    /// will be auto-estimated.
    /// </summary>
    public static Mixture Mixture(
        double gaussianWeight     = 0.0,
        double cauchyWeight       = 0.0,
        double laplacianWeight    = 0.0,
        double gaussianBandwidth  = 0.0,
        double cauchyBandwidth    = 0.0,
        double laplacianBandwidth = 0.0)
        => new(gaussianWeight, cauchyWeight, laplacianWeight,
               gaussianBandwidth, cauchyBandwidth, laplacianBandwidth);
}
