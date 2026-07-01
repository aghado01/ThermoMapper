using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Profiling.Signals;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

public sealed class SpcTabularProjectionsTests
{
    /// <summary>
    /// Regression: a rich-accumulation sweep (bond-survival counts present, so
    /// <see cref="SweepProfile.BondEntropy"/> is non-null) once threw
    /// <c>ArgumentException: Duplicate column name 'BondEntropy'</c> from
    /// <see cref="SpcTabularProjections.CreateSessionSummaryProjection"/>: the
    /// explicit chosen-T <c>BondEntropy</c> column collided with the
    /// <c>BondEntropy</c> key <see cref="ChiPeakSignalAnalyzer"/> always stuffs
    /// into <see cref="ProfileCriteria.Enrichments"/>. The two are distinct
    /// readouts (chosen-T vs the analyzer's plateau TClus), so the fix namespaces
    /// the chosen-T column as <c>BondEntropyAtChosenT</c> rather than dropping
    /// either. Asserts no throw, exactly one literal <c>BondEntropy</c> column,
    /// the disambiguated chosen-T column, and globally-unique column names.
    /// </summary>
    [Fact]
    public void CreateSessionSummaryProjection_RichAccumulation_DeDuplicatesBondEntropyColumn()
    {
        // Two replicas at one T, each concentrating its bond activity on a
        // different edge: SweepProfile.From mints a non-null BondEntropy curve
        // (the ColumnIf precondition) — mirrors SweepProfile_From_BondEntropy_*.
        var graph = BuildGraph(4, (0, 1, 1.0), (2, 3, 1.0));
        var runs = new List<SpcRunResult>
        {
            new() { Graph = graph, Accumulator = BondFrame(temperature: 1.0, draws: 10, new[] { 10, 0 }) },
            new() { Graph = graph, Accumulator = BondFrame(temperature: 1.0, draws: 10, new[] { 0, 10 }) },
        };

        SweepProfile profile = SweepProfile.From(runs);
        Assert.NotNull(profile.BondEntropy);  // the explicit-column precondition

        // The real analyzer the pipeline uses — its criteria ALWAYS carries a
        // "BondEntropy" enrichment, which is exactly what collided.
        ProfileCriteria criteria = new ChiPeakSignalAnalyzer().Analyze(profile);
        Assert.Contains("BondEntropy", criteria.Enrichments.Keys);

        var result = new SpcSessionResult(
            Partition:          Assignment.FromLabels(new[] { 0, 0, 1, 1 }),
            Profile:            profile,
            ProfileCriteria:    criteria,
            ScheduleSummary:    new SweepSummary { SubgraphNodes = 4, SubgraphEdges = 2, ChosenTemperature = 1.0 },
            Graph:              graph,
            ChosenAffinities:   new Affinities { Temperature = 1.0, G = new double[graph.Targets.Length] },
            ChosenAlignments:   null,
            ChosenCoMembership: null,
            SweepRuns:          runs,
            EvaluatorScores:    new Dictionary<string, double>());

        // Pre-fix this threw ArgumentException: Duplicate column name 'BondEntropy'.
        var projection = SpcTabularProjections.CreateSessionSummaryProjection(result);

        Assert.Equal(1, projection.Columns.Count(c => c == "BondEntropy"));       // analyzer plateau readout
        Assert.Contains("BondEntropyAtChosenT", projection.Columns);              // chosen-T readout, disambiguated
        Assert.Equal(projection.Columns.Count, projection.Columns.Distinct().Count());  // the invariant the bug broke
    }

    private static Accumulator BondFrame(double temperature, int draws, int[] bondFormedCount) => new()
    {
        Temperature = temperature, Q = 4, DrawCount = draws,
        Spins = new int[4], ClusterSizeHistogram = new int[4],
        RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
        RunningSumSqClusterSizes = 0.0, RunningSumSqClusterSizesExcl = 0.0,
        RunningSumEnergy = 0.0, RunningSumEnergySq = 0.0,
        RunningSumMag = 0.0, RunningSumMagSq = 0.0,
        BondFormedCount = bondFormedCount,
    };

    private static CsrGraph BuildGraph(int nodeCount, params (int Source, int Target, double Weight)[] edges)
    {
        var graphEdges = new Edge[edges.Length];
        for (int i = 0; i < edges.Length; i++)
            graphEdges[i] = new Edge(edges[i].Source, edges[i].Target, edges[i].Weight);
        return CsrGraph.FromEdges(graphEdges, nodeCount);
    }
}
