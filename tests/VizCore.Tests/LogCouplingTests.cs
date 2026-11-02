using System;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using Graphs.Observables;
using Graphs.Pipeline;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The compiler's log couplings (<see cref="CsrGraph.LogWeights"/>): ln J exact where the linear coupling has
/// underflowed, a pure shift under the gauge, and consistent with the linear weights wherever those are
/// representable. The oracle for the tail is the kernel's closed form on the bridge's own distance.
/// </summary>
public sealed class LogCouplingTests
{
    private static double[][] Lines(int lines, int pointsPerLine, double separation)
    {
        var points = new double[lines * pointsPerLine][];
        for (int line = 0; line < lines; line++)
            for (int i = 0; i < pointsPerLine; i++)
                points[line * pointsPerLine + i] = new[] { line * separation + 0.1 * i + 0.01 * i * i, 0.0 };
        return points;
    }

    private static GraphMetric Metric(double[][] p) => new((i, j) =>
    {
        double dx = p[i][0] - p[j][0], dy = p[i][1] - p[j][1];
        return Math.Sqrt(dx * dx + dy * dy);
    });

    private static GraphBuildResult Build(
        double[][] points, IKernelDescriptor kernel, bool? lmp = false, double floor = 1e-10,
        GaugeKind gauge = GaugeKind.None, int k = 3, RepairKind repair = RepairKind.MstMin) =>
        GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = k },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = repair },
                Projection = new CouplingProjection { Kernel = kernel, LmpRescale = lmp, CouplingFloor = floor, Gauge = gauge },
            },
            points.Length, Metric(points));

    /// <summary>The slot of the weakest edge (by log) crossing from the first half of the nodes to the second, and its endpoints.</summary>
    private static (int Slot, int From, int To) WeakestCrossing(CsrGraph g, int half)
    {
        int best = -1, from = -1;
        for (int i = 0; i < half; i++)
            for (int s = g.RowPointers[i]; s < g.RowPointers[i + 1]; s++)
                if (g.Targets[s] >= half && (best < 0 || g.LogWeights![s] < g.LogWeights[best])) { best = s; from = i; }
        return (best, from, g.Targets[best]);
    }

    [Fact]
    public void AFarGaussianBridge_KeepsItsExactLogWhereTheLinearCouplingHasUnderflowed()
    {
        double[][] points = Lines(lines: 2, pointsPerLine: 10, separation: 10.0);
        GraphBuildResult r = Build(points, new Gaussian());
        CsrGraph g = r.Graph;

        var (slot, a, b) = WeakestCrossing(g, points.Length / 2);
        double d = Metric(points).Distance(a, b);
        double bandwidth = r.SingleBandwidth!.Value;
        double expected = -(d * d) / (2.0 * bandwidth * bandwidth);

        Assert.True(expected < -1000.0, "The fixture's bridge must lie past the double range.");
        Assert.InRange(g.LogWeights![slot] / expected, 1 - 1e-12, 1 + 1e-12);
        Assert.Equal(ScalingContext.RepairEdgeMinimum, g.Weights[slot]);   // the linear weight floors at the repair minimum
    }

    [Theory]
    [InlineData("gaussian")]
    [InlineData("cauchy")]
    [InlineData("laplacian")]
    [InlineData("linear")]
    [InlineData("mixture")]
    public void LogWeights_AgreeWithTheLinearWeights_WhereTheyAreRepresentable(string name)
    {
        IKernelDescriptor kernel = name switch
        {
            "gaussian" => new Gaussian(2.0),
            "cauchy" => new Cauchy(2.0),
            "laplacian" => new Laplacian(2.0),
            "linear" => new Linear(40.0),
            _ => new Mixture(GaussianWeight: 0.3, CauchyWeight: 0.4, LaplacianWeight: 0.3,
                             GaussianBandwidth: 2.0, CauchyBandwidth: 2.0, LaplacianBandwidth: 2.0),
        };
        double[][] points = Lines(lines: 1, pointsPerLine: 14, separation: 0.0);

        CsrGraph g = Build(points, kernel, floor: 0.0).Graph;

        Assert.NotNull(g.LogWeights);
        Assert.Equal(g.Weights.Length, g.LogWeights!.Length);
        for (int s = 0; s < g.Weights.Length; s++)
        {
            Assert.True(g.Weights[s] > 1e-200, "Every edge of the fixture is representable.");
            double ln = Math.Log(g.Weights[s]);
            double tolerance = 1e-12 * Math.Max(1.0, Math.Abs(ln));
            Assert.InRange(g.LogWeights[s], ln - tolerance, ln + tolerance);
        }
    }

    [Fact]
    public void LogWeights_AreSymmetric()
    {
        double[][] points = Lines(lines: 2, pointsPerLine: 8, separation: 6.0);
        CsrGraph g = Build(points, new Gaussian()).Graph;

        for (int i = 0; i < g.NodeCount; i++)
            for (int s = g.RowPointers[i]; s < g.RowPointers[i + 1]; s++)
                Assert.Equal(g.LogWeights![s], g.LogWeights[g.FindSlot(g.Targets[s], i)]);
    }

    [Fact]
    public void TheGauge_IsAnAdditiveShiftOfLnG_InLogSpace()
    {
        double[][] points = Lines(lines: 2, pointsPerLine: 10, separation: 10.0);
        CsrGraph plain = Build(points, new Gaussian()).Graph;
        CsrGraph gauged = Build(points, new Gaussian(), gauge: GaugeKind.InverseMeanDegree).Graph;
        double shift = Math.Log(plain.NodeCount / (double)plain.Targets.Length);

        for (int s = 0; s < plain.LogWeights!.Length; s++)
            Assert.InRange(gauged.LogWeights![s] - plain.LogWeights[s], shift - 1e-12, shift + 1e-12);
    }

    [Theory]
    [InlineData(LmpRankRule.Strict)]
    [InlineData(LmpRankRule.Inclusive)]
    public void Lmp_StoresTheLogOfItsRankProducts_AndKeepsTheRepairEdgesKernelLog(LmpRankRule rank)
    {
        double[][] points = Lines(lines: 2, pointsPerLine: 10, separation: 10.0);
        GraphBuildResult r = GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 3 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection { Kernel = new Gaussian(), LmpRescale = true, LmpRank = rank },
            },
            points.Length, Metric(points));
        CsrGraph g = r.Graph;
        var (bridge, from, to) = WeakestCrossing(g, points.Length / 2);
        int mirror = g.FindSlot(to, from);

        for (int s = 0; s < g.Weights.Length; s++)
        {
            if (s == bridge || s == mirror) continue;
            if (g.Weights[s] > 0.0)
            {
                double ln = Math.Log(g.Weights[s]);
                Assert.InRange(g.LogWeights![s], ln - 1e-12, ln + 1e-12);
            }
            else
            {
                Assert.Equal(double.NegativeInfinity, g.LogWeights![s]);
            }
        }

        // The protected bridge keeps the kernel's own log, far below any rank product.
        Assert.True(g.LogWeights![bridge] < -1000.0);
        Assert.Equal(ScalingContext.RepairEdgeMinimum, g.Weights[bridge]);
    }

    [Fact]
    public void ACouplingFloorOfZero_KeepsAnUnderflowedTopologyEdgeAtItsExactLog()
    {
        // The outlier's three nominations (d near 14) underflow in the linear weights; with no floor they stay,
        // at zero linear coupling and an exact finite log. (MST-min repair symmetrizes the rows, so the edges the
        // OR rule leaves directed reach the graph.)
        double[][] line = Lines(lines: 1, pointsPerLine: 10, separation: 0.0);
        double[][] points = line.Append(new[] { 15.0, 0.0 }).ToArray();
        int outlier = points.Length - 1;

        GraphBuildResult r = Build(points, new Gaussian(), floor: 0.0);
        CsrGraph g = r.Graph;

        Assert.Equal(3, g.Degree(outlier));
        double bandwidth = r.SingleBandwidth!.Value;
        for (int s = g.RowPointers[outlier]; s < g.RowPointers[outlier + 1]; s++)
        {
            double d = Metric(points).Distance(outlier, g.Targets[s]);
            double expected = -(d * d) / (2.0 * bandwidth * bandwidth);
            Assert.Equal(0.0, g.Weights[s]);
            Assert.InRange(g.LogWeights![s] / expected, 1 - 1e-12, 1 + 1e-12);
        }

        // With the historical 1e-10 floor two of the three are cut. The strongest stays, or the outlier would be
        // cut off, with its linear weight raised to the repair edges' minimum.
        CsrGraph thinned = Build(points, new Gaussian(), floor: 1e-10).Graph;
        Assert.Equal(1, thinned.Degree(outlier));
        Assert.Equal(ScalingContext.RepairEdgeMinimum, thinned.Weights[thinned.RowPointers[outlier]]);
    }

    [Fact]
    public void CoupledConnectivity_ReadsTheLog_WhenTheGraphHasIt()
    {
        // Edge 0-1 has zero linear weight and a finite log; edge 1-2 the reverse. The log decides.
        var edges = new[] { new Edge(0, 1, 0.0, -2000.0), new Edge(1, 2, 0.5, double.NegativeInfinity) };
        CsrGraph g = CsrGraph.FromEdges(edges, 3);

        ConnectivityReport coupled = Connectivity.Validate(g, positiveWeightsOnly: true);

        Assert.Equal(2, coupled.ComponentCount);          // {0,1} and {2}
        Assert.Equal(1, coupled.EdgeCount);
        Assert.Equal(1, Connectivity.Validate(g).ComponentCount);   // structure ignores weights
    }

    [Fact]
    public void BinaryRoundTrip_CarriesTheLogs_AndTheirAbsence()
    {
        var withLog = CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.25, Math.Log(0.25)), new Edge(1, 2, 0.0, double.NegativeInfinity) }, 3);
        var without = CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.25), new Edge(1, 2, 0.5) }, 3);

        foreach (CsrGraph g in new[] { withLog, without })
        {
            using var stream = new MemoryStream();
            using (var w = new BinaryWriter(stream, System.Text.Encoding.UTF8, leaveOpen: true)) g.WriteTo(w);
            stream.Position = 0;
            CsrGraph read = CsrGraph.FromBinary(new BinaryReader(stream));

            Assert.Equal(g.Weights, read.Weights);
            Assert.Equal(g.LogWeights, read.LogWeights);
        }
        Assert.NotNull(withLog.LogWeights);
        Assert.Null(without.LogWeights);
    }

    [Fact]
    public void FromEdges_KeepsTheLogOnlyWhenEveryEdgeBringsOne()
    {
        var mixed = CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.5, Math.Log(0.5)), new Edge(1, 2, 0.5) }, 3);

        Assert.Null(mixed.LogWeights);
    }

    [Fact]
    public void InducedSubgraph_CarriesTheLogs()
    {
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 0.5, Math.Log(0.5)), new Edge(1, 2, 0.25, Math.Log(0.25)), new Edge(2, 3, 0.0, -3000.0),
        }, 4);

        CsrGraph sub = g.InducedSubgraph(new[] { false, true, true, true }, out _, out _);

        Assert.NotNull(sub.LogWeights);
        for (int s = 0; s < sub.Weights.Length; s++)
        {
            if (sub.Weights[s] == 0.0) Assert.Equal(-3000.0, sub.LogWeights![s]);
            else Assert.InRange(sub.LogWeights![s], Math.Log(sub.Weights[s]) - 1e-12, Math.Log(sub.Weights[s]) + 1e-12);
        }
    }

    [Fact]
    public void PreRepairGraph_KeepsTheLogsOfItsEdges()
    {
        double[][] points = Lines(lines: 2, pointsPerLine: 10, separation: 10.0);
        GraphBuildResult r = Build(points, new Cauchy());

        CsrGraph pre = r.PreRepairGraph!.Value;

        Assert.NotNull(pre.LogWeights);
        for (int i = 0; i < pre.NodeCount; i++)
            for (int s = pre.RowPointers[i]; s < pre.RowPointers[i + 1]; s++)
                Assert.Equal(r.Graph.LogWeights![r.Graph.FindSlot(i, pre.Targets[s])], pre.LogWeights![s]);
    }

    [Fact]
    public void IntrinsicHyperbolicGaussian_LogWeightsMatchTheLinearVanVleckCorrectedValues()
    {
        double[][] features =
        {
            new[] { 0.00, 0.00 }, new[] { 0.18, 0.05 }, new[] { 0.34, 0.11 }, new[] { 0.49, 0.16 }, new[] { 0.62, 0.30 },
        };
        var config = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 2 },
            Projection = new CouplingProjection { Kernel = new Gaussian(), LmpRescale = false, Fidelity = CouplingFidelity.Intrinsic },
        };

        GraphBuildResult r = SpcGraphBuilder.BuildResult(features, config, new PoincareMetric());
        CsrGraph g = r.Graph;

        Assert.NotNull(g.LogWeights);
        for (int s = 0; s < g.Weights.Length; s++)
        {
            Assert.True(g.Weights[s] > 0.0);
            Assert.InRange(g.LogWeights![s], Math.Log(g.Weights[s]) - 1e-9, Math.Log(g.Weights[s]) + 1e-9);
        }
    }
}
