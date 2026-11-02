using System;
using System.Linq;
using Clustering.Graphical.CumulativeField;
using Graphs;
using Graphs.Coupling;
using Graphs.Models.Potts;
using Graphs.Observables;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The cumulative-field solver on the log ladder: it agrees with the linear solver wherever both are representable, and
/// it reaches the heights of a Gaussian gap the linear solver cannot (a coupling of e^-1250 cuts at ln T = ln Hcum - ln 2).
/// The oracles are a hand-summed ladder and the exact closed form 1 - exp(-Hcum/T).
/// </summary>
public sealed class CumulativeFieldLogTests
{
    internal static double[][] LinePoints(double separation = 10.0)
    {
        var points = new double[20][];
        for (int line = 0; line < 2; line++)
            for (int i = 0; i < 10; i++)
                points[line * 10 + i] = new[] { line * separation + 0.1 * i + 0.01 * i * i, 0.0 };
        return points;
    }

    internal static GraphBuildResult TwoLines(IKernelDescriptor kernel, double separation = 10.0)
    {
        double[][] points = LinePoints(separation);
        return GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 3 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection { Kernel = kernel, LmpRescale = false },
            },
            points.Length, new GraphMetric((i, j) => Math.Abs(points[i][0] - points[j][0])));
    }

    [Theory]
    [InlineData(Field.Mean)]
    [InlineData(Field.Local)]
    public void TheLogSolve_AgreesWithTheLinearSolve_WhereBothAreRepresentable(Field field)
    {
        CsrGraph graph = TwoLines(new Cauchy()).Graph;
        CumulativeFieldContext context = CumulativeFieldSolver.Prepare(graph, EdgeWeightKind.Coupling, field);

        foreach (double t in new[] { 1e-6, 1e-3, 0.05, 0.4, 3.0 })
        {
            double[] linear = CumulativeFieldSolver.Solve(context, t).G;
            double[] log = CumulativeFieldSolver.SolveLog(context, Math.Log(t)).G;
            for (int e = 0; e < linear.Length; e++)
                Assert.InRange(log[e], linear[e] - 1e-10, linear[e] + 1e-10);
        }
    }

    [Fact]
    public void AHandSummedLadder_GivesTheClosedFormAffinities()
    {
        // Path 0-1-2-3 with J = 0.5, 0.2, 0.05. Ascending: 0.05, 0.2, 0.5, so Hcum = 0.05, 0.25, 0.75.
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 0.5, Math.Log(0.5)), new Edge(1, 2, 0.2, Math.Log(0.2)), new Edge(2, 3, 0.05, Math.Log(0.05)),
        }, 4);
        CumulativeFieldContext context = CumulativeFieldSolver.Prepare(g, EdgeWeightKind.Coupling, Field.Mean);

        double t = 0.3;
        double[] affinity = CumulativeFieldSolver.SolveLog(context, Math.Log(t)).G;

        Assert.Equal(1.0 - Math.Exp(-0.75 / t), affinity[g.FindSlot(0, 1)], 12);
        Assert.Equal(1.0 - Math.Exp(-0.25 / t), affinity[g.FindSlot(1, 2)], 12);
        Assert.Equal(1.0 - Math.Exp(-0.05 / t), affinity[g.FindSlot(2, 3)], 12);
    }

    [Fact]
    public void AGaussianGapPastTheDoubleRange_CutsWhereTheClosedFormSaysAndTheLinearSolverCannotGo()
    {
        GraphBuildResult build = TwoLines(new Gaussian());
        CsrGraph graph = build.Graph;
        double gap = MergeHeights.Compute(graph)[^1].LogCoupling;     // the bridge's ln J, below -1000
        Assert.True(gap < -1000.0);
        CumulativeFieldContext context = CumulativeFieldSolver.Prepare(graph, EdgeWeightKind.Coupling, Field.Mean);

        // The bridge is the weakest rung, so its Hcum is its own J: it cuts at ln T = gap - ln(ln 2).
        double cut = gap - Math.Log(Math.Log(2.0));
        Assert.Equal(1, CumulativeFieldSolver.ClusterLog(context, cut - 5.0).Count);   // colder: the bridge holds
        Assert.Equal(2, CumulativeFieldSolver.ClusterLog(context, cut + 5.0).Count);   // warmer: the gap is cut, the lines hold

        // The linear solver has no temperature there.
        Assert.Throws<ArgumentOutOfRangeException>(() => CumulativeFieldSolver.Solve(context, Math.Exp(cut)));
    }

    [Fact]
    public void ALogSweep_CountsMoreClustersAsTheTemperatureRises_AcrossTheGap()
    {
        CsrGraph graph = TwoLines(new Gaussian()).Graph;
        double gap = MergeHeights.Compute(graph)[^1].LogCoupling;
        double[] grid = Enumerable.Range(0, 60).Select(i => gap - 10.0 + i * ((3.0 - (gap - 10.0)) / 59.0)).ToArray();

        CumulativeFieldSweepResult sweep = CumulativeFieldSweep.RunLog(graph, EdgeWeightKind.Coupling, grid, Field.Mean);

        int[] counts = sweep.ClusterCounts();
        Assert.Equal(1, counts[0]);
        Assert.True(counts[^1] > 2);
        for (int i = 1; i < counts.Length; i++) Assert.True(counts[i] >= counts[i - 1]);
        Assert.Equal(grid, sweep.LogTemperatures);
        Assert.Equal(0.0, sweep.Temperatures[0]);                // exp(-1260) underflows; the log grid is the record
    }

    [Fact]
    public void TheLogBracket_IsTheLogOfTheLinearBracket_WhereBothAreRepresentable()
    {
        CsrGraph graph = TwoLines(new Cauchy()).Graph;

        var (lo, hi) = CumulativeFieldSolver.EstimateBracket(graph, EdgeWeightKind.Coupling, Field.Mean);
        var (lnLo, lnHi) = CumulativeFieldSolver.EstimateLogBracket(graph, EdgeWeightKind.Coupling, Field.Mean);

        Assert.InRange(lnLo, Math.Log(lo) - 1e-9, Math.Log(lo) + 1e-9);
        Assert.InRange(lnHi, Math.Log(hi) - 1e-9, Math.Log(hi) + 1e-9);
    }

    [Fact]
    public void AGraphWithoutLogWeights_StillPreparesLinearly_AndRefusesTheLogSolve()
    {
        var g = CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.5), new Edge(1, 2, 0.25) }, 3);

        CumulativeFieldContext context = CumulativeFieldSolver.Prepare(g, EdgeWeightKind.Coupling, Field.Mean);

        Assert.Equal(1, CumulativeFieldSolver.Cluster(context, 0.01).Count);
        Assert.Throws<InvalidOperationException>(() => CumulativeFieldSolver.SolveLog(context, -3.0));
    }

    [Fact]
    public void LogBondProbability_MatchesTheLinearOne_AndKeepsWhatTheLinearOneLoses()
    {
        Assert.Equal(FkKernel.BondProbability(0.4, 0.7), FkKernel.LogBondProbability(Math.Log(0.4), Math.Log(0.7)), 12);

        // A tiny ratio: 1 - exp(-1e-30) rounds to 0; the log form keeps it.
        Assert.Equal(0.0, FkKernel.BondProbability(1e-30, 1.0));
        Assert.InRange(FkKernel.LogBondProbability(Math.Log(1e-30), 0.0) / 1e-30, 1 - 1e-12, 1 + 1e-12);

        // Past the double range on both sides: a coupling of e^-1250 at a temperature of e^-1240.
        double expected = -(Math.Exp(-Math.Exp(-10.0)) - 1.0);
        Assert.InRange(FkKernel.LogBondProbability(-1250.0, -1240.0) / expected, 1 - 1e-9, 1 + 1e-9);

        Assert.Equal(1.0, FkKernel.LogBondProbability(0.0, -800.0));                       // frozen
        Assert.Equal(0.0, FkKernel.LogBondProbability(double.NegativeInfinity, 0.0));      // no coupling
    }
}
