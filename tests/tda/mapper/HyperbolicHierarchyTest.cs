// ============================================================================
// TDA.Mapper.Tests — HyperbolicHierarchyMapperTest.cs
// ============================================================================
// Regression test fixture for the MAPPER pipeline against hierarchical data
// embedded in the Poincaré ball B^3.
//
// Exercises three diagnostic levels in sequence (per the discussion):
//   1. Point-cloud MAPPER with Poincaré-radial filter
//   2. Graph MAPPER with the lifted Poincaré-radial filter
//      (FromPointFilter adapter; same lens, different clustering)
//   3. Graph MAPPER with the graph-native geodesic-distance filter
//      (BFS from the node nearest to origin)
//
// All three should produce tree-like nerves on hierarchical data. The cross-
// check at the end asserts the three results agree on tree-likeness — if they
// disagree, the pipeline is producing inconsistent topology at different
// stages, which is exactly the regression class this fixture is designed to
// catch.
//
// Self-contained: generates its own hierarchical-in-B^3 data and Poincaré
// distance function inline. No dependency on the in-flux HyperbolicHierarchy
// generator or future Poincaré manifold adapter — when those stabilize this
// test can be migrated to use them, but for now the test owns its inputs.
//
// Wiring: the helper fixture remains callable via the static Run() entry
// point, while a small wrapper fact class exposes a stable end-to-end smoke
// for xUnit and the parallel test harness.
// ============================================================================

#nullable enable
using System;
using System.Collections.Generic;
using Graphs;
using Graphs.Coupling;
using Graphs.Primitives;
using Graphs.TestSupport;
using TDA.Mapper.Clusterers;
using TDA.Mapper.Cover;
using TDA.Mapper.Diagnostics;
using TDA.Mapper.Filters;
using Xunit;

namespace TDA.Mapper.Tests;

public sealed class HyperbolicHierarchyMapperTest
{
    // ── Test parameters ─────────────────────────────────────────────────────
    //
    // Tuned for a moderately-sized test that runs in ~1 second and produces
    // a clean tree nerve when the pipeline is healthy. Adjust if extending
    // the test to stress different scales.

    private const int HierarchyLevels = 3;
    private const int BranchingFactor = 2;
    private const int PointsPerLeaf = 18;
    private const double RadiusStep = 0.20;
    private const double LeafSpread = 0.04;
    private const int DataSeed = 42;

    private const int KNearestNeighbors = 8;
    private const double KernelBandwidthAuto = 0.0;  // 0 = auto via MAD

    private const int CoverIntervals = 8;
    private const double CoverOverlap = 0.35;
    private const int KMeansK = BranchingFactor;
    private const int MaxLoopsTolerance = 2;    // small loops from cover overlap are OK

    private const double BoundaryEps = 1e-12;

    // ── Entry point ─────────────────────────────────────────────────────────

    /// <summary>
    /// Runs the full regression fixture in sequence. Throws on the first
    /// failure with a diagnostic message identifying which level failed.
    /// Returns normally on success.
    /// </summary>
    public static void Run()
    {
        var fixture = new HyperbolicHierarchyMapperTest();

        Console.WriteLine("HyperbolicHierarchyMapperTest — preparing fixture...");
        fixture.SetUp();
        Console.WriteLine($"  generated {fixture._data.Length} points across {fixture._expectedLeafCount} leaf clusters");
        Console.WriteLine($"  built CsrGraph: {fixture._graph.NodeCount} nodes, {fixture._graph.Targets.Length / 2} undirected edges");

        Console.WriteLine("\n[1/4] Level 1: point-cloud MAPPER with PoincareRadial...");
        fixture.PointCloudMapper_PoincareRadial_ShouldProduceTreeStructure();

        Console.WriteLine("[2/4] Level 2: graph MAPPER with lifted PoincareRadial...");
        fixture.GraphMapper_LiftedPoincareRadial_ShouldProduceTreeStructure();

        Console.WriteLine("[3/4] Level 3: graph MAPPER with native GeodesicDistance from origin-nearest...");
        fixture.GraphMapper_GeodesicDistance_ShouldProduceTreeStructure();

        Console.WriteLine("[4/4] Cross-check: three diagnostic levels should agree on tree-likeness...");
        fixture.AllThreeMapperVariants_ShouldAgreeOnTreeStructure();

        Console.WriteLine("\nAll HyperbolicHierarchyMapperTest assertions passed.");
    }

    // ── Fixture state ───────────────────────────────────────────────────────

    private double[][] _data = null!;
    private CsrGraph _graph;
    private int _originNearestNode;
    private int _expectedLeafCount;

    public HyperbolicHierarchyMapperTest()
    {
        SetUp();
    }

    private void SetUp()
    {
        _data = GenerateHierarchicalBallData(
            levels: HierarchyLevels,
            branchingFactor: BranchingFactor,
            pointsPerLeaf: PointsPerLeaf,
            radiusStep: RadiusStep,
            leafSpread: LeafSpread,
            seed: DataSeed);

        _expectedLeafCount = 1;
        for (int level = 0; level < HierarchyLevels; level++) _expectedLeafCount *= BranchingFactor;

        // Build a Poincaré-distance / mutual-kNN / Cauchy-kernel proximity graph.
        // Cauchy is the right kernel for hyperbolic distances (heavy-tailed
        // preserves multi-scale couplings); per spc-maturity discussion.
        Func<int, int, double> dist = (i, j) => PoincareDistance(_data[i], _data[j]);

        _graph = GraphCompilerTestPresets.BuildResult(
            n: _data.Length,
            dist: dist,
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.MutualKnn,
            k: KNearestNeighbors,
            kernel: KernelType.Cauchy,
            bandwidth: KernelBandwidthAuto,
            ensureConnected: true).Graph;

        _originNearestNode = FindOriginNearestNode(_data);
    }

    // ── Level 1: point-cloud MAPPER ─────────────────────────────────────────

    public void PointCloudMapper_PoincareRadial_ShouldProduceTreeStructure()
    {
        var result = Mapper.Build(
            data: _data,
            filter: HyperbolicFilters.PoincareRadial,
            cover: new BalancedHistogramCover(CoverIntervals, CoverOverlap),
            clusterer: new KMeansClusterer(k: KMeansK));

        ReportResult("PointCloudMapper(PoincareRadial)", result);
        AssertHealthyMapperResult(result, "PointCloud PoincareRadial");
    }

    // ── Level 2: graph MAPPER with lifted point filter ──────────────────────

    public void GraphMapper_LiftedPoincareRadial_ShouldProduceTreeStructure()
    {
        var result = Mapper.Build(
            graph: _graph,
            features: _data,
            filter: GraphFilters.FromPointFilter(HyperbolicFilters.PoincareRadial),
            cover: new BalancedHistogramCover(CoverIntervals, CoverOverlap),
            clusterer: new ConnectedComponentsClusterer());

        ReportResult("GraphMapper(lifted PoincareRadial)", result);
        AssertHealthyMapperResult(result, "Graph lifted-PoincareRadial");
    }

    // ── Level 3: graph MAPPER with native geodesic distance ─────────────────

    public void GraphMapper_GeodesicDistance_ShouldProduceTreeStructure()
    {
        var result = Mapper.Build(
            graph: _graph,
            features: null,
            filter: GraphFilters.GeodesicDistance(_originNearestNode),
            cover: new UniformCover(CoverIntervals, CoverOverlap),
            clusterer: new ConnectedComponentsClusterer());

        ReportResult("GraphMapper(GeodesicDistance)", result);
        AssertHealthyMapperResult(result, "Graph GeodesicDistance");
    }

    // ── Cross-check ─────────────────────────────────────────────────────────

    public void AllThreeMapperVariants_ShouldAgreeOnTreeStructure()
    {
        // Re-run all three. Since SetUp built the data + graph deterministically
        // (fixed seed), and KMeans uses a fixed seed too, results are reproducible.
        var pointResult = Mapper.Build(
            _data, HyperbolicFilters.PoincareRadial,
            new BalancedHistogramCover(CoverIntervals, CoverOverlap),
            new KMeansClusterer(k: KMeansK));

        var liftedResult = Mapper.Build(
            _graph, _data,
            GraphFilters.FromPointFilter(HyperbolicFilters.PoincareRadial),
            new BalancedHistogramCover(CoverIntervals, CoverOverlap),
            new ConnectedComponentsClusterer());

        var geodesicResult = Mapper.Build(
            _graph, null,
            GraphFilters.GeodesicDistance(_originNearestNode),
            new UniformCover(CoverIntervals, CoverOverlap),
            new ConnectedComponentsClusterer());

        var pointTopology = NerveTopology.From(pointResult);
        var liftedTopology = NerveTopology.From(liftedResult);
        var geodesicTopology = NerveTopology.From(geodesicResult);

        // Strong consistency: all three should classify the nerve as tree-like
        // (loop count below tolerance) OR all three should not.
        // Disagreement signals a pipeline anomaly worth investigating.
        bool p = pointTopology.LoopCount <= MaxLoopsTolerance;
        bool l = liftedTopology.LoopCount <= MaxLoopsTolerance;
        bool g = geodesicTopology.LoopCount <= MaxLoopsTolerance;

        if (!(p && l && g) && (p || l || g))
            throw new Exception(
                $"MAPPER diagnostic levels disagree on tree-likeness:\n" +
                $"  point-cloud:     {p} ({pointTopology.LoopCount} loops)\n" +
                $"  graph (lifted):  {l} ({liftedTopology.LoopCount} loops)\n" +
                $"  graph (geo):     {g} ({geodesicTopology.LoopCount} loops)\n" +
                "Pipeline anomaly — one stage is producing different topology than the others.");

        Console.WriteLine("  ✓ all three levels agree on tree-likeness");
    }

    // ── Shared assertion helper ─────────────────────────────────────────────

    private void AssertHealthyMapperResult(MapperResult result, string label)
    {
        var topology = NerveTopology.From(result);
        var nodeStats = NodeStats.From(result);
        var coverage = Coverage.From(result, _data.Length);
        var warnings = MapperWarnings.From(topology, nodeStats, coverage);

        // Sanity: produced any nodes at all
        if (topology.NerveNodeCount == 0)
            throw new Exception($"[{label}] MAPPER produced no nodes — pipeline broken upstream");

        // Connectivity: nerve should be one component (or close to it)
        // We allow up to 2 components since cover overlap can occasionally leave
        // gaps; more than that signals broken topology.
        if (topology.ConnectedComponents > 2)
            throw new Exception(
            $"[{label}] Nerve has {topology.ConnectedComponents} connected components, " +
                $"expected ≤ 2. Pipeline is fragmenting the hierarchy.");

        // Tree-likeness: loop count should be small for hierarchical data.
        if (topology.LoopCount > MaxLoopsTolerance)
            throw new Exception(
            $"[{label}] Nerve has {topology.LoopCount} loops, expected ≤ {MaxLoopsTolerance}. " +
                "Hierarchical data should produce a tree-like nerve.");

        // Node count sanity: shouldn't collapse to a single node (over-merging)
        // or explode wildly (under-merging / cover too fine).
        // Expected ballpark: somewhere between hierarchy depth and total leaf count × 2.
        int minExpected = HierarchyLevels;
        int maxExpected = _expectedLeafCount * 4;
        if (topology.NerveNodeCount < minExpected || topology.NerveNodeCount > maxExpected)
            throw new Exception(
            $"[{label}] Nerve has {topology.NerveNodeCount} nodes, expected in [{minExpected}, {maxExpected}]. " +
                "Cover or clusterer parameters may be miscalibrated.");

        // Bin populations should not have any size-1 nodes (those indicate
        // noise or cover/clusterer mismatch).
        if (nodeStats.MinSize < 2)
            throw new Exception(
            $"[{label}] Minimum bin population is {nodeStats.MinSize}; expected ≥ 2. " +
                "Singleton nodes indicate noise — increase cover overlap or reduce clusterer K.");

        Console.WriteLine($"  ✓ [{label}] {topology.NerveNodeCount} nodes, {topology.NerveEdgeCount} edges, " +
                  $"{topology.ConnectedComponents} component(s), {topology.LoopCount} loops, " +
                  $"empty bins: {result.EmptyBinCount}");
        foreach (var warning in warnings)
            Console.WriteLine($"    (warn) {warning}");
    }

    private static void ReportResult(string label, MapperResult result)
    {
        Console.WriteLine($"  {label}:");
        Console.WriteLine($"    filter:    {result.FilterName}");
        Console.WriteLine($"    cover:     {result.CoverName}");
        Console.WriteLine($"    clusterer: {result.ClustererName}");
    }

    // ── Synthetic data generation ───────────────────────────────────────────
    //
    // Hierarchical placement: root at origin; level-d nodes placed at distance
    // RadiusStep from their parent in random directions. Leaves get
    // PointsPerLeaf samples drawn from a Gaussian around the leaf center.
    //
    // Coordinates are clamped to ||x|| < 0.95 to keep all points safely inside
    // the unit ball (the Poincaré boundary causes distance blowup, which we
    // want to avoid in the test). The placement is "approximately hyperbolic"
    // — sufficient for MAPPER pipeline validation; the production
    // HyperbolicHierarchy generator (when stabilized) will use proper exp map
    // tangent transport.

    private static double[][] GenerateHierarchicalBallData(
        int levels,
        int branchingFactor,
        int pointsPerLeaf,
        double radiusStep,
        double leafSpread,
        int seed)
    {
        var rng = new Random(seed);
        var points = new List<double[]>();

        void Recurse(double[] center, int depth)
        {
            if (depth == levels)
            {
                for (int p = 0; p < pointsPerLeaf; p++)
                {
                    var pt = new double[3];
                    for (int d = 0; d < 3; d++)
                        pt[d] = center[d] + leafSpread * SampleStandardNormal(rng);
                    ClampToBall(pt, maxNorm: 0.95);
                    points.Add(pt);
                }
                return;
            }

            for (int b = 0; b < branchingFactor; b++)
            {
                var direction = SampleUnitVector(rng, dim: 3);
                var childCenter = new double[3];
                for (int d = 0; d < 3; d++)
                    childCenter[d] = center[d] + radiusStep * direction[d];
                ClampToBall(childCenter, maxNorm: 0.9);
                Recurse(childCenter, depth + 1);
            }
        }

        Recurse(new double[3], depth: 0);
        return points.ToArray();
    }

    private static int FindOriginNearestNode(double[][] data)
    {
        int best = 0;
        double bestNormSq = double.PositiveInfinity;
        for (int i = 0; i < data.Length; i++)
        {
            double normSq = 0;
            for (int d = 0; d < data[i].Length; d++) normSq += data[i][d] * data[i][d];
            if (normSq < bestNormSq) { bestNormSq = normSq; best = i; }
        }
        return best;
    }

    // ── Geometric utilities ─────────────────────────────────────────────────

    /// <summary>
    /// Poincaré distance on B^n: <c>d(x, y) = arcosh(1 + 2||x − y||² / ((1 − ||x||²)(1 − ||y||²)))</c>.
    /// </summary>
    private static double PoincareDistance(double[] x, double[] y)
    {
        double normXSq = 0, normYSq = 0, diffSq = 0;
        int dim = Math.Min(x.Length, y.Length);
        for (int d = 0; d < dim; d++)
        {
            normXSq += x[d] * x[d];
            normYSq += y[d] * y[d];
            double diff = x[d] - y[d];
            diffSq += diff * diff;
        }
        double denom = (1.0 - normXSq) * (1.0 - normYSq);
        if (denom <= BoundaryEps) return double.PositiveInfinity;
        double arg = 1.0 + 2.0 * diffSq / denom;
        return Math.Acosh(arg);
    }

    private static double SampleStandardNormal(Random rng)
    {
        // Box-Muller.
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double[] SampleUnitVector(Random rng, int dim)
    {
        var v = new double[dim];
        double normSq = 0;
        for (int d = 0; d < dim; d++)
        {
            v[d] = SampleStandardNormal(rng);
            normSq += v[d] * v[d];
        }
        double norm = Math.Sqrt(normSq);
        if (norm < 1e-12)
        {
            v[0] = 1.0;
            return v;
        }
        for (int d = 0; d < dim; d++) v[d] /= norm;
        return v;
    }

    private static void ClampToBall(double[] point, double maxNorm)
    {
        double normSq = 0;
        for (int d = 0; d < point.Length; d++) normSq += point[d] * point[d];
        double norm = Math.Sqrt(normSq);
        if (norm >= maxNorm)
        {
            double scale = maxNorm / norm;
            for (int d = 0; d < point.Length; d++) point[d] *= scale;
        }
    }
}

public sealed class HyperbolicHierarchyMapperFacts
{
    [Fact]
    public void AllThreeMapperVariants_ShouldAgreeOnTreeStructure()
        => new HyperbolicHierarchyMapperTest().AllThreeMapperVariants_ShouldAgreeOnTreeStructure();
}
