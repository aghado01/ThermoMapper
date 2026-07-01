using System;
using System.Collections.Generic;
using System.Linq;
using Graphs;
using Graphs.Coupling;
using Graphs.TestSupport;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

public sealed class MixtureKernelTests
{
    [Fact]
    public void BandwidthEstimation_ForMixture_EmptyDistances_ReturnsFallbackTuple()
    {
        MixtureBandwidth bandwidth = BandwidthEstimation.ForMixture(Array.Empty<double>(), Array.Empty<double>(), fallback: 7.5);

        Assert.Equal(7.5, bandwidth.Gaussian);
        Assert.Equal(7.5, bandwidth.Cauchy);
        Assert.Equal(7.5, bandwidth.Laplacian);
    }

    [Fact]
    public void BandwidthEstimation_ForMixture_AllEqualDistances_ReturnsZeroTuple()
    {
        double[] distances = { 3.0, 3.0, 3.0 };
        double[] scratch = new double[distances.Length];

        MixtureBandwidth bandwidth = BandwidthEstimation.ForMixture(distances, scratch);

        Assert.Equal(0.0, bandwidth.Gaussian);
        Assert.Equal(0.0, bandwidth.Cauchy);
        Assert.Equal(0.0, bandwidth.Laplacian);
    }

    [Fact]
    public void BandwidthEstimation_ForMixture_TranslatesSharedMadThroughFamilyFactors()
    {
        double[] distances = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        double[] scratch = new double[distances.Length];

        MixtureBandwidth bandwidth = BandwidthEstimation.ForMixture(distances, scratch);

        Assert.InRange(Math.Abs(bandwidth.Gaussian - 1.4826), 0.0, 1e-12);
        Assert.InRange(Math.Abs(bandwidth.Cauchy - 1.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(bandwidth.Laplacian - 1.4427), 0.0, 1e-12);
    }

    [Theory]
    [InlineData(1.0, 0.0, 0.0)]
    [InlineData(0.0, 1.0, 0.0)]
    [InlineData(0.0, 0.0, 1.0)]
    public void Mixture_Evaluate_PureKernelWeights_RoundTripExistingKernels(
        double gaussianWeight,
        double cauchyWeight,
        double laplacianWeight)
    {
        const double distance = 2.5;
        MixtureBandwidth bandwidth = new(1.2, 2.3, 3.4);
        MixtureWeights weights = new(gaussianWeight, cauchyWeight, laplacianWeight);

        double actual = MixtureKernel.Evaluate(distance, bandwidth, weights);
        double expected = gaussianWeight == 1.0
            ? GaussianKernel.Evaluate(distance, bandwidth.Gaussian)
            : cauchyWeight == 1.0
                ? CauchyKernel.Evaluate(distance, bandwidth.Cauchy)
                : LaplacianKernel.Evaluate(distance, bandwidth.Laplacian);

        Assert.Equal(expected, actual);
    }

    [Fact]
    public void Mixture_Evaluate_IsLinearInWeights()
    {
        const double alpha = 2.75;
        const double distance = 1.5;
        MixtureBandwidth bandwidth = new(2.0, 3.0, 4.0);
        MixtureWeights weights = new(0.5, 1.25, 2.0);
        MixtureWeights scaledWeights = new(alpha * weights.Gaussian, alpha * weights.Cauchy, alpha * weights.Laplacian);

        double baseline = MixtureKernel.Evaluate(distance, bandwidth, weights);
        double scaled = MixtureKernel.Evaluate(distance, bandwidth, scaledWeights);

        Assert.InRange(Math.Abs(scaled - alpha * baseline), 0.0, 1e-12);
    }

    [Fact]
    public void Mixture_Evaluate_AtZeroDistance_ReturnsWeightSum()
    {
        MixtureBandwidth bandwidth = new(1.0, 2.0, 3.0);
        MixtureWeights weights = new(0.25, 1.5, 2.25);

        double value = MixtureKernel.Evaluate(0.0, bandwidth, weights);

        Assert.Equal(weights.Gaussian + weights.Cauchy + weights.Laplacian, value);
    }

    [Theory]
    [InlineData(KernelType.Gaussian, 1.0, 0.0, 0.0)]
    [InlineData(KernelType.Cauchy, 0.0, 1.0, 0.0)]
    [InlineData(KernelType.Laplacian, 0.0, 0.0, 1.0)]
    public void BuildWithMixture_PureKernelWeights_MatchStandardBuild(
        KernelType kernel,
        double gaussianWeight,
        double cauchyWeight,
        double laplacianWeight)
    {
        double[][] points =
        {
            new[] { 0.0 },
            new[] { 1.0 },
            new[] { 3.0 },
            new[] { 6.0 }
        };

        CsrGraph standard = GraphCompilerTestPresets.BuildResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: points.Length - 1,
            kernel: kernel).Graph;

        CsrGraph mixture = GraphCompilerTestPresets.BuildWithMixtureResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            weights: new MixtureWeights(gaussianWeight, cauchyWeight, laplacianWeight),
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: points.Length - 1).Graph;

        Dictionary<long, double> standardWeights = BuildUndirectedWeightMap(standard);
        Dictionary<long, double> mixtureWeights = BuildUndirectedWeightMap(mixture);

        Assert.Equal(standardWeights.Count, mixtureWeights.Count);

        foreach ((long edgeKey, double expectedWeight) in standardWeights)
        {
            Assert.True(mixtureWeights.TryGetValue(edgeKey, out double actualWeight), $"Missing edge {edgeKey} in mixture build.");
            Assert.InRange(Math.Abs(actualWeight - expectedWeight), 0.0, 1e-12);
        }
    }

    [Fact]
    public void BuildWithMixture_EnsureConnected_PreservesWeightsOnExistingEdges()
    {
        DisconnectedControl.Fixture control = DisconnectedControl.Generate(pointsPerComponent: 8, separation: 100.0);
        double[][] points = control.Points;
        MixtureWeights weights = new(1.0, 0.5, 0.25);

        CsrGraph withoutRepair = GraphCompilerTestPresets.BuildWithMixtureResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            weights: weights,
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: 1,
            ensureConnected: false).Graph;

        CsrGraph withRepair = GraphCompilerTestPresets.BuildWithMixtureResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            weights: weights,
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: 1,
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
    public void BuildWithMixture_ExplicitBandwidthOverridesEstimation()
    {
        double[][] points =
        {
            new[] { 0.0 },
            new[] { 1.0 },
            new[] { 3.0 },
            new[] { 6.0 }
        };

        MixtureBandwidth explicitBandwidth = new(7.0, 11.0, 13.0);
        MixtureWeights weights = new(1.0, 0.0, 0.0);

        CsrGraph explicitGraph = GraphCompilerTestPresets.BuildWithMixtureResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            weights: weights,
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: points.Length - 1,
            bandwidth: explicitBandwidth).Graph;

        CsrGraph autoGraph = GraphCompilerTestPresets.BuildWithMixtureResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            weights: weights,
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: points.Length - 1).Graph;

        Dictionary<long, double> explicitWeights = BuildUndirectedWeightMap(explicitGraph);
        Dictionary<long, double> autoWeights = BuildUndirectedWeightMap(autoGraph);

        Assert.Equal(MixtureKernel.Evaluate(1.0, explicitBandwidth, weights), explicitWeights[EdgeKey(0, 1)]);
        Assert.Equal(MixtureKernel.Evaluate(2.0, explicitBandwidth, weights), explicitWeights[EdgeKey(1, 2)]);
        Assert.Equal(MixtureKernel.Evaluate(6.0, explicitBandwidth, weights), explicitWeights[EdgeKey(0, 3)]);
        Assert.Contains(explicitWeights, kvp => Math.Abs(kvp.Value - autoWeights[kvp.Key]) > 1e-6);
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

                map[EdgeKey(source, target)] = graph.Weights[edge];
            }
        }

        return map;
    }

    private static long EdgeKey(int source, int target)
        => (((long)Math.Min(source, target)) << 32) | (uint)Math.Max(source, target);

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
}
