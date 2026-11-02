using System;

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

    /// <summary>
    /// ln J of the mixture: the log-sum-exp of each weighted component's log, so a mixture whose every
    /// component has underflowed keeps its value. A component of weight 0 is left out; no component gives −∞.
    /// </summary>
    public static double LogEvaluate(double distance, MixtureBandwidth b, MixtureWeights w)
    {
        Span<double> terms = stackalloc double[3];
        int count = 0;
        if (w.Gaussian > 0.0)  terms[count++] = Math.Log(w.Gaussian)  + GaussianKernel.LogEvaluate(distance, b.Gaussian);
        if (w.Cauchy > 0.0)    terms[count++] = Math.Log(w.Cauchy)    + CauchyKernel.LogEvaluate(distance, b.Cauchy);
        if (w.Laplacian > 0.0) terms[count++] = Math.Log(w.Laplacian) + LaplacianKernel.LogEvaluate(distance, b.Laplacian);
        if (count == 0) return double.NegativeInfinity;

        double max = double.NegativeInfinity;
        for (int i = 0; i < count; i++) max = Math.Max(max, terms[i]);
        if (double.IsNegativeInfinity(max)) return double.NegativeInfinity;

        double sum = 0.0;
        for (int i = 0; i < count; i++) sum += Math.Exp(terms[i] - max);
        return max + Math.Log(sum);
    }
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
