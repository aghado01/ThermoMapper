using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json.Nodes;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;
using Clustering.Evaluation.Internal;
using Clustering.Export;
using Clustering.Graphical.CumulativeField;
using Clustering.Graphical.CumulativeField.Export;
using Clustering.Primitives;
using Graphs;
using Graphs.Observables;
using Graphs.Primitives;
using UserRepl.Commands;

namespace UserRepl;

/// <summary>
/// Cumulative-field solver session — the analytic sibling of <see cref="SpcUserSession"/>.
/// Consumes the SAME shared envelope the SW path does (a Coupling-weighted
/// <see cref="CsrGraph"/> built by the graph engine + the shared temperature
/// grid) and reduces to fewer knobs by construction. Built against the stable
/// <see cref="CumulativeFieldSolver.Prepare"/>/<see cref="CumulativeFieldSolver.Solve"/>/<see cref="CumulativeFieldSolver.Cluster"/>
/// lip only.
/// </summary>
/// <remarks>
/// <para><b>Resolution.</b> The single partition is the thermodynamic-EOM
/// resolution — the EXACT twin of SW's <c>--resolver thermal-eom</c>. The solver
/// mints the two producer-agnostic inputs from its closed form (the thermal
/// single-linkage dendrogram from per-T affinity columns, the bond-mass
/// landscape) and feeds the shared <see cref="ThermalEom.Resolve(CsrGraph,
/// Dendrogram, Landscape, int, ThermalPeripheryCompletion)"/> core. Abstained
/// leaves are completed by modal ascent when <paramref name="completion"/> is
/// <see cref="ThermalPeripheryCompletion.Ascend"/>. There is no susceptibility
/// peak to key a single chosen-T on — the resolution is cross-temperature, which
/// is the solver's strength.</para>
///
/// <para><b>Export.</b> The run writes its own bundle (<see cref="CumulativeFieldRunPaths"/>).
/// The currency tables use the same shared writer as SW:
/// <see cref="CurrencyCsvWriter.WritePartition"/> over the resolved
/// <see cref="Assignment"/>, <see cref="CurrencyCsvWriter.WriteEdgeCurrencies"/> over a representative single-T <see cref="Affinities"/> slice (the coldest /
/// most-coupled column, or <c>Run</c>'s <c>cutTemperature</c>). The MC-only
/// exports (susceptibility profile, criteria, replica traces) have no closed-form
/// analog and are absent; in their place the solver writes a thin per-temperature
/// cluster-count sweep.</para>
///
/// <para><b>A grid in ln T.</b> A <see cref="TemperatureGrid"/> of <c>ln T</c> runs the whole session on the log
/// ladder (<see cref="CumulativeFieldSolver.SolveLog"/>), so a Gaussian gap below the double range is a stage with a
/// temperature. The thermal dendrogram and the landscape then carry the
/// <see cref="ThermalDendrogram.LogTemperatureAxis"/> label, and the EOM's lifetimes and masses integrate over ln T:
/// persistence is measured in e-folds of T, where a temperature grid measures it in T. The two weigh the same stages
/// differently (a stage from T to 10T is 9T wide on a temperature grid and 2.3 e-folds on a log grid), so they can
/// resolve different partitions; the sweep gains a <c>log_temperature</c> column and the edge snapshot a
/// <c>log_weight</c> column.</para>
/// </remarks>
public static class CumulativeFieldUserSession
{
    public static CumulativeFieldUserRunResult Run(
        SpcUserDataset dataset,
        CumulativeFieldRunPaths paths,
        CsrGraph graph,
        EdgeWeightKind weightKind,
        Field field,
        SymmetrizationRule symmetrization,
        TemperatureGrid grid,
        double theta,
        int minClusterSize = 1,
        ThermalPeripheryCompletion completion = ThermalPeripheryCompletion.None,
        double? cutTemperature = null,
        IEnumerable<IExternalClusterEvaluator>? externalEvaluators = null,
        IEnumerable<IGraphPartitionEvaluator>? spcEvaluators = null,
        int[]? referenceLabels = null,
        double? cutLogTemperature = null)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        ArgumentNullException.ThrowIfNull(paths);
        ArgumentNullException.ThrowIfNull(grid);
        if (grid.Values.Length == 0)
            throw new ArgumentException("Temperature grid must be non-empty.", nameof(grid));
        if (cutTemperature is not null && cutLogTemperature is not null)
            throw new ArgumentException("Give the edge snapshot's temperature as T or as ln T, not both.", nameof(cutLogTemperature));
        double[] temperatures = grid.Values;     // T, or ln T on a log grid; the solver call below is where the two part

        string runDirectory = paths.Scope.EnsureDirectory().Dir;
        var csv = paths.Csv.EnsureDirectory();
        string sweepPath     = csv.File(CumulativeFieldCsvWriter.SweepFileName);
        string partitionPath = csv.File(CumulativeFieldCsvWriter.PartitionFileName);
        string edgesPath     = csv.File(CumulativeFieldCsvWriter.EdgeCurrenciesFileName);

        // Prepare the energy ladder once; one closed-form pass yields both the
        // per-T affinity (the edge currency) and the per-T cluster count.
        CumulativeFieldContext context = CumulativeFieldSolver.Prepare(graph, weightKind, field, symmetrization);

        int tCount = temperatures.Length;
        var affinities = new Affinities[tCount];
        var counts = new int[tCount];
        for (int i = 0; i < tCount; i++)
        {
            affinities[i] = grid.IsLog
                ? CumulativeFieldSolver.SolveLog(context, temperatures[i])
                : CumulativeFieldSolver.Solve(context, temperatures[i]);
            counts[i] = AffinityThreshold.Connect(graph, affinities[i].G, theta).Count;
        }

        // Solver-native sweep diagnostic (caller order): per-T cluster count.
        if (grid.IsLog) CumulativeFieldCsvWriter.WriteLogSweep(temperatures, counts, sweepPath);
        else            CumulativeFieldCsvWriter.WriteSweep(temperatures, counts, sweepPath);

        // Ascending view for the thermal resolution (FromEdgeCurves requires a
        // strictly ascending grid; column 0 is then the coldest/most-coupled).
        int[] order = AscendingOrder(temperatures);
        var ascTemps   = new double[tCount];
        var edgeColumns = new double[tCount][];
        var bondColumns = new double[tCount][];
        for (int k = 0; k < tCount; k++)
        {
            int s = order[k];
            ascTemps[k]    = temperatures[s];
            edgeColumns[k] = affinities[s].G;
            bondColumns[k] = AffinityNodeMarginals.BondMass(graph, affinities[s]);
        }

        // Structure (thermal single-linkage) + height (bond-mass landscape),
        // both EXACT from the closed form — the producer-agnostic inputs the
        // shared resolver eats. SW builds the same two from sampled frames.
        string axis = grid.IsLog ? ThermalDendrogram.LogTemperatureAxis : ThermalDendrogram.TemperatureAxis;
        Dendrogram dendrogram = ThermalDendrogram.FromEdgeCurves(graph, ascTemps, edgeColumns, theta, axis);
        Landscape landscape = Landscape.Create(
            axis: axis,
            grid: ascTemps,
            valuesByGridPoint: bondColumns,
            provenance: new LandscapeProvenance("BondMass", GraphId: "cumulative-field", GaugeNote: "cumulative-field:closed-form"));

        ThermalEomResult eom = ThermalEom.Resolve(graph, dendrogram, landscape, minClusterSize, completion);
        Assignment partition = eom.Assignment;
        int abstained = CountAbstained(partition.Labels);

        // Representative single-T slice for the edge snapshot: the coldest
        // (most-coupled) column the periphery ascent keys on, unless overridden.
        int repIdx = RepresentativeIndex(grid, cutTemperature, cutLogTemperature, order[0]);
        double repTemperature = grid.Temperature(repIdx);
        double repLogTemperature = grid.LogTemperature(repIdx);
        Affinities repAffinities = affinities[repIdx];

        // Currency-typed exports — identical writers/formats as the SW path.
        CurrencyCsvWriter.WritePartition(
            partition, partitionPath, features: dataset.Features, trueLabels: dataset.Labels);
        string? edgesWritten = CurrencyCsvWriter.WriteEdgeCurrencies(
            graph, repAffinities, alignments: null, edgesPath, logWeight: grid.IsLog);

        // Evaluator scores on the resolved cut (method-agnostic indices —
        // IGraphPartitionEvaluator explicitly scores the cumulative-field affinity currency).
        var evaluatorScores = new Dictionary<string, double>(StringComparer.Ordinal);
        if (referenceLabels is not null && referenceLabels.Length == partition.Labels.Length)
            foreach (var ev in externalEvaluators ?? Array.Empty<IExternalClusterEvaluator>())
                evaluatorScores[ev.Name] = ev.Evaluate(partition.Labels, referenceLabels);
        foreach (var ev in spcEvaluators ?? Array.Empty<IGraphPartitionEvaluator>())
            evaluatorScores[ev.Name] = ev.Evaluate(
                graph, repAffinities.G, partition.Labels, partition.Count);

        // Summary JSON — cumulative-field-native shape. JsonObject keeps it reflection- and
        // source-gen-free (no dependency on the SW summary's serializer context).
        var scoresNode = new JsonObject();
        foreach (var kv in evaluatorScores) scoresNode[kv.Key] = kv.Value;
        var summary = new JsonObject
        {
            ["algorithm"]              = "cumulative-field",
            ["resolver"]               = "thermal-eom",
            ["field"]                  = field.ToString(),
            ["symmetrization"]         = symmetrization.ToString(),
            ["theta"]                  = theta,
            ["minClusterSize"]         = minClusterSize,
            ["periphery"]              = completion.ToString(),
            ["nodeCount"]              = graph.NodeCount,
            ["temperatureAxis"]        = axis,
            ["temperatureCount"]       = tCount,
            ["clusterCount"]           = partition.Count,
            ["abstained"]              = abstained,
            ["representativeTemperature"] = repTemperature,
            ["representativeLogTemperature"] = repLogTemperature,
            ["representativeBy"]       = cutTemperature is null && cutLogTemperature is null ? "coldest-column" : "explicit",
            ["evaluatorScores"]        = scoresNode,
            ["runDirectory"]           = runDirectory,
            ["sweepCsv"]               = sweepPath,
            ["partitionCsv"]           = partitionPath,
            ["edgeCurrenciesCsv"]      = edgesWritten,
        };
        string summaryPath = paths.Summary;
        File.WriteAllText(summaryPath, summary.ToJsonString(s_jsonOptions));

        return new CumulativeFieldUserRunResult(
            RunDirectory:            runDirectory,
            SweepCsvPath:            sweepPath,
            PartitionCsvPath:        partitionPath,
            EdgeCurrenciesCsvPath:   edgesWritten,
            SummaryJsonPath:         summaryPath,
            ClusterCount:            partition.Count,
            Abstained:               abstained,
            RepresentativeTemperature: repTemperature,
            RepresentativeLogTemperature: repLogTemperature,
            EvaluatorScores:         evaluatorScores);
    }

    private static readonly System.Text.Json.JsonSerializerOptions s_jsonOptions = new() { WriteIndented = true };

    private static int CountAbstained(int[] labels)
    {
        int n = 0;
        for (int i = 0; i < labels.Length; i++)
            if (labels[i] == Assignment.Unassigned) n++;
        return n;
    }

    private static int[] AscendingOrder(IReadOnlyList<double> temps)
    {
        var order = new int[temps.Count];
        for (int i = 0; i < order.Length; i++) order[i] = i;
        Array.Sort(order, (a, b) => temps[a].CompareTo(temps[b]));
        return order;
    }

    /// <summary>
    /// The grid point of the edge snapshot: the one nearest the requested temperature, measured in the grid's own
    /// coordinate (T on a temperature grid, ln T on a log grid); <paramref name="coldest"/> when none was requested.
    /// </summary>
    private static int RepresentativeIndex(TemperatureGrid grid, double? cutTemperature, double? cutLogTemperature, int coldest)
    {
        if (cutTemperature is null && cutLogTemperature is null) return coldest;
        if (cutTemperature is <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(cutTemperature), "A temperature must be positive; use the ln T form for one below the double range.");
        double target = grid.IsLog
            ? cutLogTemperature ?? Math.Log(cutTemperature!.Value)
            : cutTemperature ?? Math.Exp(cutLogTemperature!.Value);
        return NearestIndex(grid.Values, target);
    }

    private static int NearestIndex(IReadOnlyList<double> temps, double t)
    {
        int best = 0;
        double bestDist = Math.Abs(temps[0] - t);
        for (int i = 1; i < temps.Count; i++)
        {
            double d = Math.Abs(temps[i] - t);
            if (d < bestDist) { bestDist = d; best = i; }
        }
        return best;
    }
}

/// <summary>Paths and headline figures from one <see cref="CumulativeFieldUserSession.Run"/> call.</summary>
public sealed record CumulativeFieldUserRunResult(
    string RunDirectory,
    string SweepCsvPath,
    string PartitionCsvPath,
    string? EdgeCurrenciesCsvPath,
    string SummaryJsonPath,
    int ClusterCount,
    int Abstained,
    double RepresentativeTemperature,
    double RepresentativeLogTemperature,
    IReadOnlyDictionary<string, double> EvaluatorScores);
