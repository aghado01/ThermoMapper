using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using Graphs.Primitives;
using Graphs.TestSupport;
using Maths.Geometry;
using Xunit;

namespace VizCore.Tests;

public sealed class GraphBuilderBandwidthTests
{
    [Fact]
    public void Build_EnsureConnected_PreservesWeightsOnExistingEdges()
    {
        DisconnectedControl.Fixture control = DisconnectedControl.Generate(pointsPerComponent: 8, separation: 100.0);
        double[][] points = control.Points;

        CsrGraph withoutRepair = GraphCompilerTestPresets.BuildResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: 1,
            kernel: KernelType.Gaussian,
            ensureConnected: false).Graph;

        CsrGraph withRepair = GraphCompilerTestPresets.BuildResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: 1,
            kernel: KernelType.Gaussian,
            ensureConnected: true).Graph;

        Dictionary<long, double> withoutWeights = BuildUndirectedWeightMap(withoutRepair);
        Dictionary<long, double> withWeights = BuildUndirectedWeightMap(withRepair);

        Assert.NotEmpty(withoutWeights);
        Assert.True(withWeights.Count > withoutWeights.Count, "Fixture must add at least one MST bridge edge.");

        foreach ((long edgeKey, double expectedWeight) in withoutWeights)
        {
            Assert.True(withWeights.TryGetValue(edgeKey, out double repairedWeight), $"Missing original edge {edgeKey} after EnsureConnected.");
            Assert.InRange(Math.Abs(repairedWeight - expectedWeight), 0.0, 1e-12);
        }
    }

    [Fact]
    public void Build_WithPoincareMetric_UsesHyperbolicGeometryAwareBandwidth()
    {
        double[][] features =
        {
            new[] { 0.00, 0.00 },
            new[] { 0.18, 0.05 },
            new[] { 0.34, 0.11 },
            new[] { 0.49, 0.16 },
        };

        var metric = new PoincareMetric();
        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, CreateAutoBandwidthConfig(), metric);

        double expected = EstimateGaussianBandwidth(build.DirectedSelection.KthNeighborDistances, metric.Properties);
        double legacy = EstimateLegacyHyperbolicPlaceholderBandwidth(build.DirectedSelection.KthNeighborDistances);

        Assert.Equal(SpaceGeometry.Hyperbolic, build.Metric?.Geometry);
        Assert.NotNull(build.SingleBandwidth);
        Assert.Equal(expected, build.SingleBandwidth!.Value, 12);
        Assert.True(Math.Abs(build.SingleBandwidth.Value - legacy) > 1e-9,
            $"Expected geometry-aware Poincare bandwidth to differ from the old hyperbolic placeholder factor path. New={build.SingleBandwidth.Value}, old={legacy}.");
    }

    [Fact]
    public void Build_WithFisherRaoHalfPlaneMetric_UsesHyperbolicLogScaleBandwidth()
    {
        double[][] features =
        {
            new[] { -1.0, Math.Log(0.7) },
            new[] {  0.0, Math.Log(1.0) },
            new[] {  1.2, Math.Log(1.8) },
            new[] {  2.1, Math.Log(3.2) },
        };

        var metric = new FisherRaoHalfPlaneMetric();
        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, CreateAutoBandwidthConfig(), metric);

        double expected = EstimateGaussianBandwidth(build.DirectedSelection.KthNeighborDistances, metric.Properties);
        double legacy = EstimateLegacyHyperbolicPlaceholderBandwidth(build.DirectedSelection.KthNeighborDistances);

        Assert.Equal(SpaceGeometry.Hyperbolic, build.Metric?.Geometry);
        Assert.Equal(BandwidthStrategy.LogScaleHyperbolic, build.Metric?.BandwidthStrategy);
        Assert.NotNull(build.SingleBandwidth);
        Assert.Equal(expected, build.SingleBandwidth!.Value, 12);
        Assert.True(Math.Abs(build.SingleBandwidth.Value - legacy) > 1e-9,
            $"Expected FisherRaoHalfPlane bandwidth to differ from the old hyperbolic placeholder factor path. New={build.SingleBandwidth.Value}, old={legacy}.");
    }

    [Fact]
    public void BandwidthForHyperbolicMixture_UsesPerKernelFactors()
    {
        double[] sample = { 0.09, 0.18, 0.36, 0.72, 1.44 };
        double[] scratch = new double[sample.Length];

        MixtureBandwidth bandwidth = BandwidthEstimation.ForMixture(
            sample,
            scratch,
            BandwidthStrategy.LogScaleHyperbolic,
            SpaceGeometry.Hyperbolic);

        Assert.Equal(
            BandwidthEstimation.LogScaleBandwidth(sample, new double[sample.Length], BandwidthEstimation.GaussianHyperbolicFactor),
            bandwidth.Gaussian,
            12);
        Assert.Equal(
            BandwidthEstimation.LogScaleBandwidth(sample, new double[sample.Length], BandwidthEstimation.CauchyHyperbolicFactor),
            bandwidth.Cauchy,
            12);
        Assert.Equal(
            BandwidthEstimation.LogScaleBandwidth(sample, new double[sample.Length], BandwidthEstimation.LaplacianHyperbolicFactor),
            bandwidth.Laplacian,
            12);

        Assert.True(bandwidth.Cauchy < bandwidth.Gaussian,
            $"Expected hyperbolic Cauchy fallback to stay below Gaussian. Cauchy={bandwidth.Cauchy}, Gaussian={bandwidth.Gaussian}");
        Assert.True(bandwidth.Gaussian < bandwidth.Laplacian,
            $"Expected Laplacian hyperbolic factor to exceed Gaussian. Gaussian={bandwidth.Gaussian}, Laplacian={bandwidth.Laplacian}");
    }

    [Fact]
    public void Build_WithPoincareMetric_AutoFidelityMatchesExplicitGeodesicLinear()
    {
        double[][] features = CreatePoincareBallFixture();
        var metric = new PoincareMetric();

        GraphBuildResult autoBuild = SpcGraphBuilder.BuildResult(features, CreatePoincareConfig(), metric);
        GraphBuildResult explicitBuild = SpcGraphBuilder.BuildResult(
            features,
            CreatePoincareConfig(fidelity: CouplingFidelity.GeodesicLinear),
            metric);

        Assert.Equal(autoBuild.SingleBandwidth, explicitBuild.SingleBandwidth);
        Assert.Equal(BuildUndirectedWeightMap(autoBuild.Graph), BuildUndirectedWeightMap(explicitBuild.Graph));
        AssertContainsScalingMessage(autoBuild, "Fidelity=GeodesicLinear");
        AssertContainsScalingMessage(explicitBuild, "Fidelity=GeodesicLinear");
    }

    [Fact]
    public void Build_WithPoincareMetric_IntrinsicSuppressesLongerEdges()
    {
        double[][] features = CreatePoincareBallFixture();
        var metric = new PoincareMetric();

        GraphBuildResult linearBuild = SpcGraphBuilder.BuildResult(
            features,
            CreatePoincareConfig(kernel: new Gaussian(0.0), fidelity: CouplingFidelity.GeodesicLinear),
            metric);
        GraphBuildResult intrinsicBuild = SpcGraphBuilder.BuildResult(
            features,
            CreatePoincareConfig(kernel: new Gaussian(0.0), fidelity: CouplingFidelity.Intrinsic),
            metric);

        Dictionary<long, double> linearWeights = BuildUndirectedWeightMap(linearBuild.Graph);
        Dictionary<long, double> intrinsicWeights = BuildUndirectedWeightMap(intrinsicBuild.Graph);

        double shortDistance = metric.Distance(features[0], features[1]);
        double mediumDistance = metric.Distance(features[1], features[3]);
        double longDistance = metric.Distance(features[0], features[4]);

        long shortKey = EdgeKey(0, 1);
        long mediumKey = EdgeKey(1, 3);
        long longKey = EdgeKey(0, 4);

        Assert.NotNull(linearBuild.SingleBandwidth);
        Assert.NotNull(intrinsicBuild.SingleBandwidth);
        Assert.Equal(
            EstimateGaussianBandwidth(linearBuild.DirectedSelection.KthNeighborDistances, metric.Properties),
            linearBuild.SingleBandwidth!.Value,
            12);
        Assert.Equal(
            EstimateIntrinsicGaussianBandwidth(intrinsicBuild.DirectedSelection.KthNeighborDistances, features[0].Length),
            intrinsicBuild.SingleBandwidth!.Value,
            12);
        Assert.True(Math.Abs(intrinsicWeights[shortKey] - linearWeights[shortKey]) < 1e-3,
            $"Expected near-origin edge to stay close under Intrinsic fidelity. r={shortDistance}, linear={linearWeights[shortKey]}, intrinsic={intrinsicWeights[shortKey]}");
        Assert.True(intrinsicWeights[mediumKey] < linearWeights[mediumKey],
            $"Expected Intrinsic fidelity to suppress medium hyperbolic edge. r={mediumDistance}, linear={linearWeights[mediumKey]}, intrinsic={intrinsicWeights[mediumKey]}");
        Assert.True(intrinsicWeights[longKey] < linearWeights[longKey],
            $"Expected Intrinsic fidelity to suppress long hyperbolic edge. r={longDistance}, linear={linearWeights[longKey]}, intrinsic={intrinsicWeights[longKey]}");

        AssertContainsScalingMessage(linearBuild, "Fidelity=GeodesicLinear");
        AssertContainsScalingMessage(intrinsicBuild, "Fidelity=Intrinsic");
        AssertContainsScalingWarning(linearBuild, "GeodesicLinear coupling on Hyperbolic geometry");
        Assert.DoesNotContain(
            intrinsicBuild.Diagnostics?.Warnings ?? Array.Empty<Graphs.Diagnostics.DiagnosticMessage>(),
            message => message.Stage == "Scaling" && message.Text.Contains("GeodesicLinear coupling", StringComparison.Ordinal));
    }

    [Fact]
    public void Build_WithSphericalMetric_UsesIntrinsicGaussianSphericalBandwidth()
    {
        double[][] features = CreateSphericalFixture();
        var metric = new SphericalGeodesicMetric();

        // Gaussian(0.0) requests auto-estimation: the intrinsic spherical heat
        // bandwidth is derived from the NN sample (an explicit bandwidth > 0 is
        // respected verbatim, as the hyperbolic A/B harness relies on). Mirrors
        // the Poincaré intrinsic sibling and the spherical monotonicity harness.
        GraphBuildResult build = SpcGraphBuilder.BuildResult(
            features,
            CreateSphericalConfig(kernel: new Gaussian(0.0)),
            metric);

        Assert.Equal(SpaceGeometry.Spherical, build.Metric?.Geometry);
        Assert.NotNull(build.SingleBandwidth);
        Assert.Equal(
            EstimateIntrinsicSphericalBandwidth(build.DirectedSelection.KthNeighborDistances, features[0].Length),
            build.SingleBandwidth!.Value,
            12);
        AssertContainsScalingMessage(build, "Fidelity=Intrinsic");
    }

    [Fact]
    public void IntrinsicSphericalVolumeCorrection_UsesSinCorrection()
    {
        // S^3 (ambient n=4, intrinsic m=3) so the Van Vleck exponent (m-1)/2 == 1
        // and the correction is the clean base factor r/sin(r) — the positive-curvature
        // mirror of IntrinsicHyperbolicVolumeCorrection_InThreeDimensions (H^3, exponent 1).
        // DO NOT drop this to SphericalManifold(3): that is S^2 (m=2), whose correct
        // correction is (r/sin r)^(1/2), per the (m-1)/2 spec (off-by-one ward).
        var scaler = new Graphs.Pipeline.Scalers.GlobalBandwidthScaler(
            KernelType.Gaussian,
            bandwidth: 1.0,
            strategy: BandwidthStrategy.QuantileNormalized,
            geometry: SpaceGeometry.Spherical,
            fidelity: CouplingFidelity.Intrinsic,
            ambientDimension: 4,
            manifold: new SphericalManifold(4));

        System.Reflection.MethodInfo? volumeCorrection = typeof(Graphs.Pipeline.Scalers.GlobalBandwidthScaler)
            .GetMethod("VolumeCorrection", System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic);

        Assert.NotNull(volumeCorrection);

        foreach (double distance in new[] { 1e-9, 0.03, 0.15, 0.75, 1.15 })
        {
            double actual = (double)volumeCorrection.Invoke(scaler, new object[] { distance })!;
            double expected = distance / Math.Sin(distance);
            Assert.InRange(Math.Abs(actual - expected), 0.0, 1e-12);
            Assert.True(actual >= 1.0);
        }
    }

    [Fact]
    public void IntrinsicHyperbolicGaussianBandwidth_RecoversSyntheticScale()
    {
        const int ambientDimension = 3;
        const double expectedBandwidth = 0.60;

        double[] sample = SimulateIntrinsicNearestNeighborSample(ambientDimension, expectedBandwidth, count: 1024, seed: 1729);
        double recoveredBandwidth = BandwidthEstimation.ForIntrinsicGaussianHyperbolic(sample, ambientDimension);

        Assert.InRange(RelativeError(recoveredBandwidth, expectedBandwidth), 0.0, 0.20);
        Assert.InRange(
            RelativeError(
                BandwidthEstimation.IntrinsicHeatTimeFromBandwidth(recoveredBandwidth),
                BandwidthEstimation.IntrinsicHeatTimeFromBandwidth(expectedBandwidth)),
            0.0,
            0.36);
    }

    [Fact]
    public void IntrinsicHyperbolicVolumeCorrection_InThreeDimensions_MatchesLegacyPocExactly()
    {
        var scaler = new Graphs.Pipeline.Scalers.GlobalBandwidthScaler(
            KernelType.Gaussian,
            bandwidth: 1.0,
            strategy: BandwidthStrategy.LogScaleHyperbolic,
            geometry: SpaceGeometry.Hyperbolic,
            fidelity: CouplingFidelity.Intrinsic,
            ambientDimension: 3);

        System.Reflection.MethodInfo? volumeCorrection = typeof(Graphs.Pipeline.Scalers.GlobalBandwidthScaler)
            .GetMethod("VolumeCorrection", System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic);

        Assert.NotNull(volumeCorrection);

        foreach (double distance in new[] { 1e-9, 0.03, 0.15, 0.75, 2.5, 4.0 })
        {
            double actual = (double)volumeCorrection.Invoke(scaler, new object[] { distance })!;
            Assert.Equal(LegacyIntrinsicPocCorrection(distance), actual);
        }
    }

    [Fact]
    public void Build_WithHyperbolicIntrinsic_AndMissingAmbientDimension_Throws()
    {
        double[][] features = CreatePoincareBallFixture();
        var metric = new PoincareMetric();
        var substrate = new GraphMetric(
            Distance: (i, j) => metric.Distance(features[i], features[j]),
            Properties: metric.Properties,
            AmbientDimension: null);

        InvalidOperationException ex = Assert.Throws<InvalidOperationException>(() =>
            GraphCompiler.Build(
                CreatePoincareConfig(fidelity: CouplingFidelity.Intrinsic),
                features.Length,
                substrate));

        Assert.Contains("requires the ambient dimension", ex.Message, StringComparison.Ordinal);
    }

    [Fact]
    public void Build_WithIntrinsicFidelity_AndNonGaussianKernel_Throws()
    {
        double[][] features = CreatePoincareBallFixture();
        var metric = new PoincareMetric();

        NotSupportedException ex = Assert.Throws<NotSupportedException>(() =>
            SpcGraphBuilder.BuildResult(
                features,
                CreatePoincareConfig(kernel: new Cauchy(0.0), fidelity: CouplingFidelity.Intrinsic),
                metric));

        Assert.Contains("Gaussian kernel only", ex.Message, StringComparison.Ordinal);
    }


    private static Dictionary<long, double> BuildUndirectedWeightMap(CsrGraph graph)
    {
        var map = new Dictionary<long, double>();

        for (int source = 0; source < graph.NodeCount; source++)
        {
            int rowStart = graph.RowPointers[source];
            int rowEnd = graph.RowPointers[source + 1];
            for (int edge = rowStart; edge < rowEnd; edge++)
            {
                int target = graph.Targets[edge];
                if (target <= source)
                    continue;

                long key = (((long)source) << 32) | (uint)target;
                map[key] = graph.Weights[edge];
            }
        }

        return map;
    }

    private static long EdgeKey(int left, int right)
    {
        int source = Math.Min(left, right);
        int target = Math.Max(left, right);
        return (((long)source) << 32) | (uint)target;
    }

    private static double EuclideanDistance(double[] left, double[] right)
    {
        double sum = 0.0;
        for (int i = 0; i < left.Length; i++)
        {
            double diff = left[i] - right[i];
            sum += diff * diff;
        }

        return Math.Sqrt(sum);
    }

    private static GraphCompilerConfig CreateAutoBandwidthConfig() => new()
    {
        Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 2 },
        Filter = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
        Repair = new RepairConfig { Kind = RepairKind.NoRepair },
        Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
        Projection = new CouplingProjection { Kernel = new Gaussian(0.0), LmpRescale = false },
    };

    private static double EstimateGaussianBandwidth(double[] sample, MetricProperties props)
    {
        double[] scratch = new double[sample.Length];
        return BandwidthEstimation.ForGaussian(sample, scratch, props.BandwidthStrategy, props.Geometry);
    }

    private static double EstimateIntrinsicGaussianBandwidth(double[] sample, int ambientDimension)
        => BandwidthEstimation.ForIntrinsicGaussianHyperbolic(sample, ambientDimension);

    private static double EstimateIntrinsicSphericalBandwidth(double[] sample, int ambientDimension)
        => BandwidthEstimation.ForIntrinsicGaussianSpherical(sample, ambientDimension - 1);

    private static double GaussianWeight(double radius, double bandwidth)
        => Math.Exp(-(radius * radius) / (2.0 * bandwidth * bandwidth));

    private static double EstimateLegacyHyperbolicPlaceholderBandwidth(double[] sample)
    {
        double[] scratch = new double[sample.Length];
        return BandwidthEstimation.LogScaleBandwidth(sample, scratch, 1.0);
    }

    private static double LegacyIntrinsicPocCorrection(double distance)
    {
        if (distance < 1e-12)
            return 1.0;

        return distance / Math.Sinh(distance);
    }

    private static double[] SimulateIntrinsicNearestNeighborSample(int ambientDimension, double bandwidth, int count, int seed)
    {
        const int candidatePoolSize = 32;
        double radiusMax = Math.Max(8.0, ((ambientDimension - 1) * bandwidth * bandwidth) / 2.0 + 12.0 * bandwidth + ambientDimension);
        double logMax = EstimateIntrinsicLogDensityMax(ambientDimension, bandwidth, radiusMax);
        var rng = new Random(seed);
        var samples = new double[count];
        var pool = new double[candidatePoolSize];
        int accepted = 0;
        int attempts = 0;
        int maxAttempts = count * candidatePoolSize * 20000;

        while (accepted < count)
        {
            int poolFill = 0;
            while (poolFill < candidatePoolSize)
            {
                if (attempts++ >= maxAttempts)
                {
                    throw new InvalidOperationException(
                        $"Intrinsic nearest-neighbor sampler stalled for d={ambientDimension} bandwidth={bandwidth:G4}.");
                }

                double radius = rng.NextDouble() * radiusMax;
                double logAccept = IntrinsicLogDensity(ambientDimension, bandwidth, radius) - logMax;
                if (Math.Log(rng.NextDouble()) <= logAccept)
                    pool[poolFill++] = radius;
            }

            Array.Sort(pool);
            samples[accepted++] = pool[0];
        }

        return samples;
    }

    private static double EstimateIntrinsicLogDensityMax(int ambientDimension, double bandwidth, double radiusMax)
    {
        const int gridCount = 4096;
        double logMax = double.NegativeInfinity;

        for (int index = 1; index <= gridCount; index++)
        {
            double radius = radiusMax * index / gridCount;
            double logDensity = IntrinsicLogDensity(ambientDimension, bandwidth, radius);
            if (logDensity > logMax)
                logMax = logDensity;
        }

        return logMax;
    }

    private static double IntrinsicLogDensity(int ambientDimension, double bandwidth, double radius)
    {
        double kernelTerm = -(radius * radius) / (2.0 * bandwidth * bandwidth);
        double exponent = (ambientDimension - 1) / 2.0;
        if (exponent <= 0.0)
            return kernelTerm;
        if (radius <= 0.0)
            return double.NegativeInfinity;

        return kernelTerm + exponent * (Math.Log(radius) + LogSinh(radius));
    }

    private static double LogSinh(double radius)
    {
        if (radius < 1e-8)
            return Math.Log(Math.Max(radius, 1e-12));

        if (radius < 20.0)
            return Math.Log(Math.Sinh(radius));

        return radius - Math.Log(2.0);
    }

    private static double RelativeError(double estimate, double truth)
        => Math.Abs(estimate - truth) / Math.Max(1e-12, truth);

    private static double EuclideanNorm(ReadOnlySpan<double> vector)
    {
        double sumSq = 0.0;
        for (int i = 0; i < vector.Length; i++)
            sumSq += vector[i] * vector[i];

        return Math.Sqrt(sumSq);
    }

    private static double[][] CreatePoincareBallFixture() =>
    [
        new[] { 0.00, 0.00, 0.00 },
        new[] { 0.03, 0.00, 0.00 },
        new[] { 0.15, 0.00, 0.00 },
        new[] { 0.45, 0.00, 0.00 },
        new[] { 0.75, 0.00, 0.00 },
    ];

    private static double[][] CreateSphericalFixture() =>
    [
        new[] { 1.0, 0.0, 0.0 },
        new[] { 0.0, 1.0, 0.0 },
        new[] { 0.0, 0.0, 1.0 },
        new[] { 1.0 / Math.Sqrt(2.0), 1.0 / Math.Sqrt(2.0), 0.0 },
        new[] { 1.0 / Math.Sqrt(2.0), 0.0, 1.0 / Math.Sqrt(2.0) },
    ];

    private static GraphCompilerConfig CreateSphericalConfig(
        IKernelDescriptor? kernel = null,
        CouplingFidelity fidelity = CouplingFidelity.Auto) => new()
    {
        Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
        Filter = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
        Repair = new RepairConfig { Kind = RepairKind.NoRepair },
        Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
        Projection = new CouplingProjection
        {
            Kernel = kernel ?? new Gaussian(1.0),
            LmpRescale = false,
            Fidelity = fidelity,
        },
    };

    private static GraphCompilerConfig CreatePoincareConfig(
        IKernelDescriptor? kernel = null,
        CouplingFidelity fidelity = CouplingFidelity.Auto) => new()
    {
        Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
        Filter = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
        Repair = new RepairConfig { Kind = RepairKind.NoRepair },
        Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
        Projection = new CouplingProjection
        {
            Kernel = kernel ?? new Gaussian(1.0),
            LmpRescale = false,
            Fidelity = fidelity,
        },
    };

    private static void AssertContainsScalingMessage(GraphBuildResult build, string snippet)
    {
        IEnumerable<Graphs.Diagnostics.DiagnosticMessage> messages =
            build.Diagnostics?.Messages ?? System.Linq.Enumerable.Empty<Graphs.Diagnostics.DiagnosticMessage>();
        Assert.Contains(
            messages,
            message => message.Stage == "Scaling"
                && message.Text.Contains(snippet, StringComparison.Ordinal));
    }

    private static void AssertContainsScalingWarning(GraphBuildResult build, string snippet)
    {
        IEnumerable<Graphs.Diagnostics.DiagnosticMessage> warnings =
            build.Diagnostics?.Warnings ?? System.Linq.Enumerable.Empty<Graphs.Diagnostics.DiagnosticMessage>();
        Assert.Contains(
            warnings,
            message => message.Stage == "Scaling"
                && message.Text.Contains(snippet, StringComparison.Ordinal));
    }
}
