using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Primitives;
using Maths.Rng;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The bracket spanning a graph from its first break to its last melt, read from the couplings alone
/// (primer §13.21), against hand-computed ends, a brute-force merge tree and a run of the sampler.
/// </summary>
public sealed class SpanBracketTests
{
    private const int Q = 20;

    private static readonly double Melt = Math.Log(1.0 + Math.Sqrt(Q));
    private static readonly double Bridge = 4.0 * Math.Log(1.0 + Q);
    private static readonly double Rung = Math.Pow(80.0, 1.0 / 47.0);

    /// <summary>Two cliques of <paramref name="side"/> nodes at unit coupling, joined by one edge of coupling <paramref name="bridge"/>.</summary>
    private static CsrGraph TwoCliques(int side, double bridge)
    {
        var edges = new List<Edge>();
        for (int offset = 0; offset <= side; offset += side)
            for (int i = 0; i < side; i++)
                for (int j = i + 1; j < side; j++)
                    edges.Add(new Edge(offset + i, offset + j, 1.0));
        edges.Add(new Edge(0, side, bridge));
        return CsrGraph.FromEdges(edges.ToArray(), 2 * side);
    }

    private static GraphBuildResult Build(double[][] points, IKernelDescriptor kernel, GaugeKind gauge) =>
        GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection { Kernel = kernel, LmpRescale = false, BandwidthOverride = BandwidthStrategy.MeanEdgeDistance, Gauge = gauge },
            },
            points.Length,
            new GraphMetric((i, j) =>
            {
                double dx = points[i][0] - points[j][0], dy = points[i][1] - points[j][1];
                return Math.Sqrt(dx * dx + dy * dy);
            }));

    /// <summary>Two runs of twelve points a <paramref name="gap"/> apart.</summary>
    private static double[][] TwoRuns(double gap)
    {
        var points = new double[24][];
        for (int i = 0; i < 12; i++)
        {
            points[i] = new[] { 0.1 * i + 0.01 * i * i, 0.0 };
            points[12 + i] = new[] { gap + 0.12 * i + 0.012 * i * i, 0.5 };
        }
        return points;
    }

    private static double ColdestSplitByBruteForce(int n, IReadOnlyList<Edge> edges, int m, double floor)
    {
        foreach (double j in edges.Select(e => e.J).Where(w => w > floor).Distinct().OrderBy(w => w))
        {
            int[] above = Components(n, edges.Where(e => e.J > j));
            int[] at = Components(n, edges.Where(e => e.J >= j));
            var large = Enumerable.Range(0, n).GroupBy(i => above[i]).Where(g => g.Count() >= m).ToList();
            if (large.GroupBy(g => at[g.First()]).Any(g => g.Count() >= 2))
                return j;
        }
        return double.NaN;
    }

    private static int[] Components(int n, IEnumerable<Edge> edges)
    {
        int[] label = Enumerable.Range(0, n).ToArray();
        Edge[] list = edges.ToArray();
        for (bool changed = true; changed;)
        {
            changed = false;
            foreach (Edge e in list)
            {
                int low = Math.Min(label[e.Source], label[e.Target]);
                if (label[e.Source] != low || label[e.Target] != low)
                {
                    label[e.Source] = label[e.Target] = low;
                    changed = true;
                }
            }
        }
        return label;
    }

    [Fact]
    public void TwoCliquesAndABridge_TakeTheHandComputedEnds()
    {
        (double lo, double hi, int levels, string basis) = SpcScheduleHelpers.SpanBracket(TwoCliques(12, 0.01), Q, null, 5);

        Assert.Equal(11.0 / Melt, hi, 12);
        Assert.Equal(0.01 / Bridge, lo, 12);
        Assert.Equal((int)Math.Ceiling(Math.Log(hi / lo) / Math.Log(Rung)) + 1, levels);
        Assert.Equal("span(m=5,J_first=0.01)", basis);
    }

    [Fact]
    public void WithNoSplitIntoPartsThatLarge_TheBottomIsTheTopOver80()
    {
        (double lo, double hi, int levels, string basis) = SpcScheduleHelpers.SpanBracket(TwoCliques(12, 0.01), Q, null, 13);

        Assert.Equal(hi / 80.0, lo);
        Assert.Equal(48, levels);
        Assert.Equal("span(m=13,no-split)", basis);
    }

    [Fact]
    public void ACouplingAtOrUnderTheFloor_IsNoSplit_AndTheGaugeCarriesTheFloorWithIt()
    {
        // Far enough apart that the only bridge is the repair edge, a Gaussian coupling far under the floor.
        foreach (GaugeKind gauge in new[] { GaugeKind.None, GaugeKind.InverseMeanDegree })
        {
            GraphBuildResult build = Build(TwoRuns(100.0), new Gaussian(), gauge);
            Edge[] edges = EdgesOf(build.Graph);
            double floor = build.Resolved!.Coupling!.CouplingFloor * build.Resolved.Coupling.GaugeFactor;
            Assert.Contains(edges, e => e.J <= floor);

            var declared = SpcScheduleHelpers.SpanBracket(build.Graph, Q, build.Resolved, 5);
            double undeclared = SpcScheduleHelpers.SpanBracket(build.Graph, Q, null, 5).Lo;

            double split = ColdestSplitByBruteForce(24, edges, 5, floor);
            Assert.Equal(double.IsNaN(split) ? declared.Hi / 80.0 : split / Bridge, declared.Lo, 1e-9 * declared.Lo);
            Assert.Equal(ColdestSplitByBruteForce(24, edges, 5, 0.0) / Bridge, undeclared, 1e-9 * undeclared);
            Assert.True(declared.Lo > 1e6 * undeclared, $"{declared.Lo} against {undeclared}");
        }
    }

    [Theory]
    [InlineData(8.0)]
    [InlineData(100.0)]
    public void TheGaugeMovesTheBracketByExactlyG(double gap)
    {
        foreach (IKernelDescriptor kernel in new IKernelDescriptor[] { new Gaussian(), new Cauchy() })
        {
            GraphBuildResult plain = Build(TwoRuns(gap), kernel, GaugeKind.None);
            GraphBuildResult gauged = Build(TwoRuns(gap), kernel, GaugeKind.InverseMeanDegree);
            double g = gauged.Resolved!.Coupling!.GaugeFactor;

            var atPlain = SpcScheduleHelpers.SpanBracket(plain.Graph, Q, plain.Resolved, 5);
            var atGauged = SpcScheduleHelpers.SpanBracket(gauged.Graph, Q, gauged.Resolved, 5);

            Assert.InRange(atGauged.Lo / (g * atPlain.Lo), 1 - 1e-9, 1 + 1e-9);
            Assert.InRange(atGauged.Hi / (g * atPlain.Hi), 1 - 1e-9, 1 + 1e-9);
            Assert.Equal(atPlain.Levels, atGauged.Levels);
        }
    }

    [Fact]
    public void EdgesTiedAtACoupling_MergeAtOnce_SoTheirOrderDoesNotMatter()
    {
        // A path of equal couplings: walked from one end it never joins two pairs, walked from both it would.
        Edge[] fromOneEnd = { new(0, 1, 1.0), new(1, 2, 1.0), new(2, 3, 1.0) };
        Edge[] fromBoth = { new(0, 1, 1.0), new(2, 3, 1.0), new(1, 2, 1.0) };

        Assert.Equal("span(m=2,no-split)", SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(fromOneEnd, 4), Q, null, 2).Basis);
        Assert.Equal("span(m=2,no-split)", SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(fromBoth, 4), Q, null, 2).Basis);
    }

    [Fact]
    public void ASpanTooWideForTheLevels_KeepsItsEnds_AndWidensItsRungs()
    {
        (double lo, double hi, int levels, _) = SpcScheduleHelpers.SpanBracket(
            CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0), new Edge(1, 2, 1e-290) }, 3), Q, null, 1);

        Assert.Equal(240, levels);
        Assert.Equal(1e-290 / Bridge, lo, 1e-9 * lo);
        Assert.True(hi > 0.5);
    }

    [Fact]
    public void WhenFewNodesHoldTheCouplings_TheStrongestNodeSetsTheTop()
    {
        // Most nodes have no coupling, so the percentile is zero.
        var lone = SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0) }, 30), Q, null, 1);
        Assert.Equal(1.0 / Melt, lone.Hi, 12);
        Assert.Equal(1.0 / Bridge, lone.Lo, 12);

        // Weak pairs hold the percentile well under a strong split between two strong pairs.
        var edges = new List<Edge> { new(0, 1, 2.0), new(2, 3, 2.0), new(1, 2, 1.0) };
        for (int i = 4; i < 76; i += 2)
            edges.Add(new Edge(i, i + 1, 1e-6));
        var strong = SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(edges.ToArray(), 76), Q, null, 2);
        Assert.Equal(3.0 / Melt, strong.Hi, 12);
        Assert.Equal(1.0 / Bridge, strong.Lo, 12);
    }

    [Fact]
    public void AGraphWithNoCouplings_StillGivesAGrid_AndBadArgumentsThrow()
    {
        foreach (CsrGraph empty in new[]
        {
            CsrGraph.FromEdges(Array.Empty<Edge>(), 5),
            CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.0) }, 2),
            CsrGraph.FromEdges(new[] { new Edge(0, 1, double.Epsilon) }, 2),
        })
        {
            (double lo, double hi, int levels, string basis) = SpcScheduleHelpers.SpanBracket(empty, Q, null, 1);
            Assert.True(lo > 0.0 && hi > lo);
            Assert.Equal(48, levels);
            Assert.Equal("span(no-couplings)", basis);
        }

        // A subnormal coupling beside a normal one counts as none.
        var beside = SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0), new Edge(1, 2, 1e-320) }, 3), Q, null, 1);
        Assert.Equal(1.0 / Bridge, beside.Lo, 12);

        CsrGraph graph = TwoCliques(4, 0.1);
        Assert.Throws<ArgumentOutOfRangeException>(() => SpcScheduleHelpers.SpanBracket(graph, 1, null, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpcScheduleHelpers.SpanBracket(graph, Q, null, 0));
    }

    [Fact]
    public void ASeededSweep_AgreesWithABruteForceMergeTree_AndHoldsItsShape()
    {
        var rng = new Xoshiro256PlusPlus(seed: 20261104);
        int[] qs = { 2, 3, 5, 10, 20, 40, 100 };
        int capped = 0, noSplit = 0, fellBack = 0;

        for (int trial = 0; trial < 600; trial++)
        {
            int n = 2 + rng.NextInt(39);
            int q = qs[rng.NextInt(qs.Length)];
            int m = 1 + rng.NextInt(n + 1);
            double[] exponents = (trial % 3) switch
            {
                0 => new[] { -6.0, 0.0 },
                1 => new[] { -300.0, 300.0 },
                _ => new[] { -1.0, 1.0 },
            };

            Edge[] edges = RandomEdges(rng, n, 1 + rng.NextInt(3 * n), trial % 3 == 2, exponents);
            var (lo, hi, levels, _) = SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(edges, n), q, null, m);

            var strength = new double[n];
            foreach (Edge e in edges)
            {
                strength[e.Source] += e.J;
                strength[e.Target] += e.J;
            }
            Array.Sort(strength);

            double melt = Math.Log(1.0 + Math.Sqrt(q));
            double split = ColdestSplitByBruteForce(n, edges, m, 0.0);
            double expectedLo = split / (4.0 * Math.Log(1.0 + q));
            double expectedHi = strength[(int)Math.Floor(0.9 * (n - 1))] / melt;
            if (expectedHi <= 0.0 || expectedLo >= expectedHi)
            {
                expectedHi = strength[^1] / melt;
                fellBack++;
            }
            if (double.IsNaN(split))
            {
                expectedLo = expectedHi / 80.0;
                noSplit++;
            }

            Assert.True(Math.Abs(lo / expectedLo - 1.0) < 1e-12, $"trial {trial}: lo {lo} against {expectedLo}");
            Assert.True(Math.Abs(hi / expectedHi - 1.0) < 1e-12, $"trial {trial}: hi {hi} against {expectedHi}");
            Assert.True(lo > 0.0 && lo < hi && double.IsFinite(hi), $"trial {trial}: [{lo}, {hi}]");
            Assert.InRange(levels, 2, 240);

            double[] grid = SpcScheduleHelpers.LogSpaceGrid(lo, hi, levels);
            Assert.InRange(grid[0] / lo, 1 - 1e-9, 1 + 1e-9);
            Assert.InRange(grid[^1] / hi, 1 - 1e-9, 1 + 1e-9);
            if (levels < 240)
                Assert.True(grid[1] / grid[0] <= Rung * (1 + 1e-9), $"trial {trial}: rung {grid[1] / grid[0]}");
            else
                capped++;

            // The graph's labels and its edges' order say nothing about the temperatures.
            int[] relabel = Enumerable.Range(0, n).OrderBy(_ => rng.NextDouble()).ToArray();
            Edge[] renamed = edges.OrderBy(_ => rng.NextDouble()).Select(e => new Edge(relabel[e.Source], relabel[e.Target], e.J)).ToArray();
            var other = SpcScheduleHelpers.SpanBracket(CsrGraph.FromEdges(renamed, n), q, null, m);
            Assert.True(Math.Abs(other.Lo / lo - 1.0) < 1e-12, $"trial {trial}: relabelled lo {other.Lo} against {lo}");
            Assert.True(Math.Abs(other.Hi / hi - 1.0) < 1e-12, $"trial {trial}: relabelled hi {other.Hi} against {hi}");
        }

        // The sweep reached the branches it's here for.
        Assert.True(capped > 0 && noSplit > 0 && fellBack > 0, $"capped {capped}, no split {noSplit}, fell back {fellBack}");
    }

    private static Edge[] RandomEdges(Xoshiro256PlusPlus rng, int n, int count, bool tied, double[] exponents)
    {
        double[] ties = { 0.25, 0.5, 1.0 };
        var seen = new HashSet<(int, int)>();
        var edges = new List<Edge>();
        for (int k = 0; k < count; k++)
        {
            int a = rng.NextInt(n), b = rng.NextInt(n);
            if (a == b || !seen.Add((Math.Min(a, b), Math.Max(a, b))))
                continue;
            double exponent = exponents[0] + (exponents[1] - exponents[0]) * rng.NextDouble();
            edges.Add(new Edge(a, b, tied ? ties[rng.NextInt(ties.Length)] : Math.Pow(10.0, exponent)));
        }
        if (edges.Count == 0)
            edges.Add(new Edge(0, 1, 1.0));
        return edges.ToArray();
    }

    private static Edge[] EdgesOf(CsrGraph graph)
    {
        var edges = new List<Edge>();
        foreach (UndirectedEdge e in graph.UndirectedEdges())
            edges.Add(new Edge(e.Source, e.Target, graph.Weights[e.Slot]));
        return edges.ToArray();
    }

    [Fact]
    public void TheSampler_ReadsTheGraphJoinedAtTheBottom_MeltedAtTheTop_AndSplitBetween()
    {
        CsrGraph graph = TwoCliques(12, 0.01);
        (double lo, double hi, int levels, _) = SpcScheduleHelpers.SpanBracket(graph, Q, null, 5);

        (double Chi, double Mean)[] curve = SpcScheduleHelpers.LogSpaceGrid(lo, hi, levels).Select(t =>
        {
            SwRunResult run = SwRunner.Run(new SwRunSpec
            {
                Graph = graph, Temperature = t, Q = Q, ReplicaIndex = 0,
                Seed = SpcSeedHelper.Derive(20261104, t, 0),
                Budget = new RunBudget(50, 100),
            });
            return (run.Accumulator.RunningSumSqClusterSizesExcl / (run.Accumulator.DrawCount * (double)graph.NodeCount), run.FkSusceptibility);
        }).ToArray();

        Assert.True(curve[0].Chi < 25.0 / graph.NodeCount, $"bottom's χ {curve[0].Chi}");
        Assert.True(curve[^1].Mean < 2.0, $"top's mean cluster {curve[^1].Mean}");
        Assert.True(curve.Max(c => c.Chi) > 4.0, $"largest χ {curve.Max(c => c.Chi)}");
    }

    [Fact]
    public void TheCommand_ResolvesSpan_ByTheReadingsClusterSize_AndRecordsTheGridExactly()
    {
        GraphBuildResult build = Build(TwoRuns(8.0), new Cauchy(), GaugeKind.InverseMeanDegree);
        var options = new SpcCommand.Options { TemperaturesSpec = " SPAN ", MinClusterSize = 5 };

        (double[] grid, string resolved, string basis) = SpcCommand.ResolveTemperatureGrid(options, build.Graph, build.Resolved);

        (double lo, double hi, int levels, string expectedBasis) = SpcScheduleHelpers.SpanBracket(build.Graph, options.Q, build.Resolved, 5);
        Assert.Equal(expectedBasis, basis);
        Assert.Equal(levels, grid.Length);
        Assert.InRange(grid[0] / lo, 1 - 1e-12, 1 + 1e-12);
        Assert.InRange(grid[^1] / hi, 1 - 1e-12, 1 + 1e-12);
        Assert.Equal(grid, TemperatureGridSpec.Parse(resolved));

        // A size under 1 counts as 1, as the lineage reading takes it.
        options.MinClusterSize = 0;
        Assert.Equal(SpcScheduleHelpers.SpanBracket(build.Graph, options.Q, build.Resolved, 1).Basis,
            SpcCommand.ResolveTemperatureGrid(options, build.Graph, build.Resolved).Basis);
    }

    [Fact]
    public void TheCommand_WarnsOnlyWhenALoneNodeWouldSetTheBottom()
    {
        Assert.NotNull(SpcCommand.SpanClusterSizeWarning(new SpcCommand.Options { TemperaturesSpec = "span" }));
        Assert.NotNull(SpcCommand.SpanClusterSizeWarning(new SpcCommand.Options { TemperaturesSpec = " Span ", MinClusterSize = 0 }));
        Assert.Null(SpcCommand.SpanClusterSizeWarning(new SpcCommand.Options { TemperaturesSpec = "span", MinClusterSize = 2 }));
        Assert.Null(SpcCommand.SpanClusterSizeWarning(new SpcCommand.Options { TemperaturesSpec = "auto" }));
        Assert.Null(SpcCommand.SpanClusterSizeWarning(new SpcCommand.Options()));
    }

    [Fact]
    public void APreset_StatesTheClusterSizeASpanGridReads()
    {
        var options = new SpcCommand.Options();

        new SpcPreset { Temperatures = "span", MinClusterSize = 5 }.ApplyTo(options);

        Assert.Equal(5, options.MinClusterSize);
        Assert.Null(SpcCommand.SpanClusterSizeWarning(options));
    }

    [Fact]
    public void AutoStaysAsItWas()
    {
        GraphBuildResult build = Build(TwoRuns(8.0), new Cauchy(), GaugeKind.InverseMeanDegree);
        var options = new SpcCommand.Options { TemperaturesSpec = "auto" };

        (double[] grid, _, string basis) = SpcCommand.ResolveTemperatureGrid(options, build.Graph, build.Resolved);

        (double lo, double hi, string expectedBasis) = SpcScheduleHelpers.AutoBracket(build.Graph, options.Q, build.Resolved);
        Assert.Equal(expectedBasis, basis);
        Assert.Equal(SpcScheduleHelpers.LogSpaceGrid(lo, hi, 48), grid);
    }
}
