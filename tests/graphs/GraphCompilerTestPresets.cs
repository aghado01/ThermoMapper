using System;
using Graphs.Coupling;
using Graphs.Distance;

namespace Graphs.TestSupport;

/// <summary>
/// Test presets: legacy graph-builder semantics with explicit
/// <see cref="GraphCompilerConfig"/> construction.
/// </summary>
public static class GraphCompilerTestPresets
{
    public static GraphBuildResult BuildResult(
        int n,
        Func<int, int, double> dist,
        TopologyKind topologyKind = TopologyKind.Knn,
        FilterKind filterKind = FilterKind.OrRule,
        int k = 10,
        double epsilon = 0.0,
        KernelType kernel = KernelType.Gaussian,
        double bandwidth = 0.0,
        bool ensureConnected = false,
        BandwidthStrategy bandwidthStrategy = BandwidthStrategy.MadConsistencyFactor,
        MutualBandwidthSource mutualBandwidthSource = MutualBandwidthSource.DirectedKth)
    {
        IKernelDescriptor kernelDescriptor = kernel switch
        {
            KernelType.Gaussian  => new Gaussian(bandwidth),
            KernelType.Cauchy    => new Cauchy(bandwidth),
            KernelType.Laplacian => new Laplacian(bandwidth),
            KernelType.Linear    => new Linear(bandwidth),
            _ => new Gaussian(bandwidth),
        };

        var config = BuildConfig(
            topologyKind,
            filterKind,
            k,
            epsilon,
            kernelDescriptor,
            ensureConnected,
            RefinementKind.Auto,
            false,
            bandwidthStrategy,
            mutualBandwidthSource);

        return GraphCompiler.Build(config, n, new GraphMetric(dist));
    }

    public static GraphBuildResult BuildWithMixtureResult(
        int n,
        Func<int, int, double> dist,
        MixtureWeights weights,
        TopologyKind topologyKind = TopologyKind.Knn,
        FilterKind filterKind = FilterKind.OrRule,
        int k = 10,
        double epsilon = 0.0,
        MixtureBandwidth? bandwidth = null,
        bool ensureConnected = false,
        BandwidthStrategy bandwidthStrategy = BandwidthStrategy.MadConsistencyFactor,
        MutualBandwidthSource mutualBandwidthSource = MutualBandwidthSource.DirectedKth)
    {
        var mixture = new Mixture(
            weights.Gaussian, weights.Cauchy, weights.Laplacian,
            bandwidth?.Gaussian ?? 0.0,
            bandwidth?.Cauchy ?? 0.0,
            bandwidth?.Laplacian ?? 0.0);

        var config = BuildConfig(
            topologyKind,
            filterKind,
            k,
            epsilon,
            mixture,
            ensureConnected,
            RefinementKind.Auto,
            false,
            bandwidthStrategy,
            mutualBandwidthSource);

        return GraphCompiler.Build(config, n, new GraphMetric(dist));
    }

    private static GraphCompilerConfig BuildConfig(
        TopologyKind topologyKind,
        FilterKind filterKind,
        int k,
        double epsilon,
        IKernelDescriptor kernel,
        bool ensureConnected,
        RefinementKind refinement,
        bool? lmpRescale,
        BandwidthStrategy bandwidthStrategy,
        MutualBandwidthSource mutualBandwidthSource)
    {
        return new GraphCompilerConfig
        {
            Topology = topologyKind == TopologyKind.EpsilonBall
                ? new TopologyConfig { Kind = TopologyKind.EpsilonBall, Epsilon = epsilon }
                : new TopologyConfig { Kind = TopologyKind.Knn, K = k },
            Filter = new FilterConfig
            {
                Kind = filterKind,
                MutualBandwidthSource = mutualBandwidthSource,
            },
            Repair = new RepairConfig
            {
                Kind = ensureConnected ? RepairKind.MstMin : RepairKind.NoRepair,
            },
            Refinement = new RefinementConfig
            {
                Kind = refinement,
            },
            Projection = new CouplingProjection
            {
                Kernel = kernel,
                LmpRescale = lmpRescale,
                BandwidthOverride = bandwidthStrategy,
            },
        };
    }
}
