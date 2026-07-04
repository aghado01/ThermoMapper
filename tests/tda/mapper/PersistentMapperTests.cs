#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;
using Maths.Topology;
using TDA.Mapper;
using TDA.Mapper.Clusterers;
using TDA.Mapper.Cover;
using TDA.Mapper.Filters;
using TDA.Ph.Nerves;
using Xunit;

using TDA.Ph;
namespace TDA.Mapper.Tests;

/// <summary>
/// Tests for Persistent Mapper over a scalar parameter sweep.
/// Exit criterion (Phase 3): barcode from a clean cluster fixture whose H0
/// birth/death features correspond to structural transitions on the sweep axis.
/// </summary>
public sealed class PersistentMapperTests
{
    // ── Unit test: trivial two-node nerve ─────────────────────────────────────

    [Fact]
    public void PathGraph_H0Barcode_TwoComponentsThenOne()
    {
        // Frame 0 (T=0.0): nerve has 2 isolated nodes → 2 H0 components.
        // Frame 1 (T=1.0): nerve has 2 nodes connected by 1 edge → 1 H0 component.
        // Expected barcode: 2 H0 bars (both born at 0.0), one dies at 1.0, one persists.

        var frame0 = new NerveFiltrationFrame(
            ParameterValue: 0.0,
            Nerve: CsrGraph.FromEdges(Array.Empty<Edge>(), 2),
            NodeMemberIndices: new[] { new[] { 0, 1 }, new[] { 2, 3 } },
            FrameIndex: 0);

        var frame1 = new NerveFiltrationFrame(
            ParameterValue: 1.0,
            Nerve: CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0) }, 2),
            NodeMemberIndices: new[] { new[] { 0, 1 }, new[] { 2, 3 } },
            FrameIndex: 1);

        var filtration = new NerveFiltration(new[] { frame0, frame1 }, "T");
        Barcode barcode = PersistenceBarcode.ComputeH0(filtration);

        var h0 = barcode.Bars.Where(b => b.Dimension == 0).ToList();
        Assert.Equal(2, h0.Count);
        Assert.Single(h0, b => b.IsInfinite);
        Assert.Single(h0, b => !b.IsInfinite && b.Death == 1.0);
    }

    // ── End-to-end: two-Gaussian cluster fixture ──────────────────────────────

    // Two Gaussians separated by 4 units (spread 0.5 each).
    // Sweep an epsilon-ball graph: at small epsilon the two clusters form two
    // disconnected subgraphs → nerve has 2 H0 components.  At epsilon ≥ the
    // inter-cluster gap the graph connects → nerve collapses to 1 component.
    //
    // The H0 barcode captures this structural transition: a finite bar whose
    // death value marks the epsilon at which the two clusters first merge.
    // This corresponds to the SPC phase transition on the same parameter axis
    // (bandwidth / scale), satisfying the Phase 3 exit criterion.

    [Fact]
    public void TwoGaussian_PersistentMapperSweep_H0BarcodeCapturesClusterMerge()
    {
        const int nPerCluster = 30;
        const double separation = 4.0;
        const double spread = 0.5;

        double[][] data = GenerateTwoGaussians(nPerCluster, separation, spread, seed: 42);

        // Sweep epsilon from below intra-cluster scale to above inter-cluster scale.
        // Intra-cluster typical distance ≈ 0.6–0.8; inter-cluster minimum ≈ 2.0.
        double[] epsilons = new[] { 0.7, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0 };

        Barcode barcode = PersistentMapper.SweepH0(
            epsilons,
            eps => BuildGraphMapperAtEpsilon(data, eps),
            parameterLabel: "epsilon (T proxy)");

        // Phase 3 exit criterion: H0 barcode has ≥ 2 bars (one per initial cluster)
        // and at least one finite bar capturing the cluster-merge transition.
        var h0 = barcode.Bars.Where(b => b.Dimension == 0).ToList();
        Assert.True(h0.Count >= 2,
            $"Expected ≥ 2 H0 bars for two-cluster fixture, got {h0.Count}. " +
            "Both clusters should appear as separate components at low epsilon.");

        var finiteBars = h0.Where(b => !b.IsInfinite).ToList();
        Assert.True(finiteBars.Count >= 1,
            "Expected ≥ 1 finite H0 bar representing the cluster-merge transition.");

        // The merge happens when epsilon first bridges the inter-cluster gap.
        // With separation=4.0 and spread=0.5, the closest inter-cluster points
        // are ≈ 2.0–3.0 apart.  The finite bar's death should fall in that range.
        double minFiniteDeath = finiteBars.Min(b => b.Death);
        Assert.InRange(minFiniteDeath, 0.7, 5.0);
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    private static MapperResult BuildGraphMapperAtEpsilon(double[][] data, double eps)
    {
        int n = data.Length;

        // Build an epsilon-ball graph inline (O(n²), fine for n=60 in tests).
        // At small eps: two clusters form two disconnected subgraphs.
        // At large eps: inter-cluster edges appear → connected.
        var edges = new List<Edge>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                double d = EuclideanDist(data[i], data[j]);
                if (d <= eps)
                    edges.Add(new Edge(i, j, Math.Exp(-(d * d) / (eps * eps))));
            }

        CsrGraph graph = CsrGraph.FromEdges(edges.ToArray(), n);

        return Mapper.Build(
            graph: graph,
            features: data,
            filter: GraphFilters.WeightedDegree,
            cover: new BalancedHistogramCover(numIntervals: 5, overlapPercent: 0.40),
            clusterer: new ConnectedComponentsClusterer());
    }

    private static double[][] GenerateTwoGaussians(
        int nPerCluster, double separation, double spread, int seed)
    {
        var rng = new Random(seed);
        var pts = new List<double[]>(2 * nPerCluster);
        for (int i = 0; i < nPerCluster; i++)
            pts.Add(new[] { -separation / 2 + spread * SampleNormal(rng), spread * SampleNormal(rng) });
        for (int i = 0; i < nPerCluster; i++)
            pts.Add(new[] {  separation / 2 + spread * SampleNormal(rng), spread * SampleNormal(rng) });
        return pts.ToArray();
    }

    private static double SampleNormal(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double EuclideanDist(double[] a, double[] b)
    {
        double sum = 0;
        for (int d = 0; d < a.Length; d++)
        {
            double diff = a[d] - b[d];
            sum += diff * diff;
        }
        return Math.Sqrt(sum);
    }
}
