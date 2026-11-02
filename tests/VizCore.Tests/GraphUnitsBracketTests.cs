using System;
using Clustering.Graphical.CumulativeField;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Models.Potts;
using Graphs.Primitives;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The temperature bracket a graph declares through its resolved construction
/// (graph-compiler-brief §4.1 finding 5, decision 5): one function, read by the REPL
/// and the library alike, and equivariant under the gauge.
/// </summary>
public sealed class GraphUnitsBracketTests
{
    private static readonly double[][] Points = MakePoints();

    private static double[][] MakePoints()
    {
        var points = new double[24][];
        for (int i = 0; i < 12; i++)
        {
            points[i] = new[] { 0.1 * i + 0.01 * i * i, 0.0 };
            points[12 + i] = new[] { 8.0 + 0.12 * i + 0.012 * i * i, 0.5 };
        }
        return points;
    }

    private static GraphBuildResult Build(IKernelDescriptor kernel, GaugeKind gauge, BandwidthStrategy? strategy = null, bool? lmp = false) =>
        GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection { Kernel = kernel, LmpRescale = lmp, BandwidthOverride = strategy, Gauge = gauge },
            },
            Points.Length,
            new GraphMetric((i, j) =>
            {
                double dx = Points[i][0] - Points[j][0], dy = Points[i][1] - Points[j][1];
                return Math.Sqrt(dx * dx + dy * dy);
            }));

    [Fact]
    public void ClassicUnits_TakeTheQOnlyBracket()
    {
        GraphBuildResult classic = Build(new Gaussian(), GaugeKind.InverseMeanDegree, BandwidthStrategy.MeanEdgeDistance);

        (double lo, double hi, string basis) = SpcScheduleHelpers.AutoBracket(classic.Graph, 20, classic.Resolved);

        (double expectedLo, double expectedHi) = BwdPottsCriticalEstimate.Bracket(20);
        Assert.Equal(expectedLo, lo);
        Assert.Equal(expectedHi, hi);
        Assert.Equal("bwd-t_ps(q=20)", basis);
    }

    [Fact]
    public void AnyOtherConstruction_TakesTheHeuristicDividedByTheGauge_SoTheGaugeMovesTheBracketByExactlyG()
    {
        // The same Cauchy graph with and without the gauge: J → gJ must move the bracket T → gT.
        GraphBuildResult plain = Build(new Cauchy(), GaugeKind.None);
        GraphBuildResult gauged = Build(new Cauchy(), GaugeKind.InverseMeanDegree);
        double g = plain.Graph.NodeCount / (double)plain.Graph.Targets.Length;

        var atPlain = SpcScheduleHelpers.AutoBracket(plain.Graph, 20, plain.Resolved);
        var atGauged = SpcScheduleHelpers.AutoBracket(gauged.Graph, 20, gauged.Resolved);

        Assert.InRange(atGauged.Lo / (g * atPlain.Lo), 1 - 1e-9, 1 + 1e-9);
        Assert.InRange(atGauged.Hi / (g * atPlain.Hi), 1 - 1e-9, 1 + 1e-9);
        Assert.Equal("edge-statistics", atPlain.Basis);
        Assert.StartsWith("edge-statistics/gauge=", atGauged.Basis, StringComparison.Ordinal);
    }

    [Fact]
    public void CauchyOrLmpGraphs_NoLongerTakeTheQOnlyBracketJustForTheMeanEdgeBandwidth()
    {
        GraphBuildResult cauchy = Build(new Cauchy(), GaugeKind.InverseMeanDegree, BandwidthStrategy.MeanEdgeDistance);
        GraphBuildResult lmp = Build(new Gaussian(), GaugeKind.InverseMeanDegree, BandwidthStrategy.MeanEdgeDistance, lmp: true);

        Assert.StartsWith("edge-statistics", SpcScheduleHelpers.AutoBracket(cauchy.Graph, 20, cauchy.Resolved).Basis, StringComparison.Ordinal);
        Assert.StartsWith("edge-statistics", SpcScheduleHelpers.AutoBracket(lmp.Graph, 20, lmp.Resolved).Basis, StringComparison.Ordinal);
    }

    [Fact]
    public void ALibraryGridAndTheReplsAutoGridAreTheSameBracket()
    {
        foreach (GraphBuildResult build in new[]
        {
            Build(new Gaussian(), GaugeKind.InverseMeanDegree, BandwidthStrategy.MeanEdgeDistance),
            Build(new Cauchy(), GaugeKind.InverseMeanDegree),
            Build(new Cauchy(), GaugeKind.None),
        })
        {
            (double lo, double hi, _) = SpcScheduleHelpers.AutoBracket(build.Graph, 20, build.Resolved);

            var grid = AutoGridFixedSweep.BuildConfig(build.Graph, construction: build.Resolved).Temperatures;

            Assert.InRange(grid[0] / lo, 1 - 1e-9, 1 + 1e-9);
            Assert.InRange(grid[^1] / hi, 1 - 1e-9, 1 + 1e-9);
        }
    }

    [Fact]
    public void AGraphWithNoDeclaration_TakesTheHeuristicAsItStands()
    {
        GraphBuildResult build = Build(new Cauchy(), GaugeKind.None);

        var undeclared = SpcScheduleHelpers.AutoBracket(build.Graph, 20, construction: null);

        (double lo, double hi) = SpcScheduleHelpers.EstimateBracket(build.Graph, 20, 0.05, 5.0);
        Assert.Equal(lo, undeclared.Lo);
        Assert.Equal(hi, undeclared.Hi);
    }

    [Fact]
    public void FieldLocalOnAnLmpGraph_WarnsThatTheNormalizationsStack()
    {
        GraphBuildResult lmp = Build(new Cauchy(), GaugeKind.None, lmp: true);
        GraphBuildResult plain = Build(new Cauchy(), GaugeKind.None, lmp: false);

        Assert.NotNull(CumulativeFieldCommand.StackedNormalizationWarning(Field.Local, lmp.Resolved));
        Assert.Null(CumulativeFieldCommand.StackedNormalizationWarning(Field.Mean, lmp.Resolved));
        Assert.Null(CumulativeFieldCommand.StackedNormalizationWarning(Field.Local, plain.Resolved));
    }
}
