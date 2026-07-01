using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json.Nodes;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;
using Clustering.Evaluation.Internal;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Partitions.Thermal;
using Clustering.Graphical.SPC.Runtime.Core.Solver;
using Clustering.Primitives;
using Graphs;
using Graphs.Observables;
using Graphs.Primitives;

namespace UserRepl;

/// <summary>
/// PKWang solver session — the analytic sibling of <see cref="SpcUserSession"/>.
/// Consumes the SAME shared envelope the SW path does (a Coupling-weighted
/// <see cref="CsrGraph"/> built by the graph engine + the shared temperature
/// grid) and reduces to fewer knobs by construction. Built against the stable
/// <see cref="PKWang.Prepare"/>/<see cref="PKWang.Solve"/>/<see cref="PKWang.Cluster"/>
/// lip only.
/// </summary>
/// <remarks>
/// <para><b>Resolution.</b> The single partition is the thermodynamic-EOM
/// resolution — the EXACT twin of SW's <c>--resolver thermal-eom</c>. PKWang
/// mints the two producer-agnostic inputs from its closed form (the thermal
/// single-linkage dendrogram from per-T affinity columns, the bond-mass
/// landscape) and feeds the shared <see cref="ThermalEom.Resolve(CsrGraph,
/// Dendrogram, Landscape, int, ThermalPeripheryCompletion)"/> core. Abstained
/// leaves are completed by modal ascent when <paramref name="completion"/> is
/// <see cref="ThermalPeripheryCompletion.Ascend"/>. There is no susceptibility
/// peak to key a single chosen-T on — the resolution is cross-temperature, which
/// is the solver's strength.</para>
///
/// <para><b>Export parity.</b> The currency-typed exports are written with the
/// IDENTICAL writers as SW: <see cref="SpcCsvWriter.WritePartition"/> over the
/// resolved <see cref="Assignment"/>, <see cref="SpcCsvWriter.WriteEquilibriumEdges"/>
/// over a representative single-T <see cref="Affinities"/> slice (the coldest /
/// most-coupled column, or <paramref name="cutTemperature"/>). The MC-only
/// exports (susceptibility profile, criteria, replica traces) have no closed-form
/// analog and are absent; in their place the solver writes a thin per-temperature
/// cluster-count sweep.</para>
/// </remarks>
public static class PKWangUserSession
{
    public static PKWangUserRunResult Run(
        SpcUserDataset dataset,
        SpcRunPaths paths,
        CsrGraph graph,
        EdgeWeightKind weightKind,
        Field field,
        SymmetrizationRule symmetrization,
        IReadOnlyList<double> temperatures,
        double theta,
        int minClusterSize = 1,
        ThermalPeripheryCompletion completion = ThermalPeripheryCompletion.None,
        double? cutTemperature = null,
        IEnumerable<IExternalClusterEvaluator>? externalEvaluators = null,
        IEnumerable<IGraphPartitionEvaluator>? spcEvaluators = null,
        int[]? referenceLabels = null)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        ArgumentNullException.ThrowIfNull(paths);
        ArgumentNullException.ThrowIfNull(temperatures);
        if (temperatures.Count == 0)
            throw new ArgumentException("Temperature grid must be non-empty.", nameof(temperatures));

        string runDirectory = paths.Scope.EnsureDirectory().Dir;
        var csv = paths.Csv.EnsureDirectory();
        string sweepPath     = csv.File(SpcCsvWriter.SweepFileName);
        string partitionPath = csv.File(SpcCsvWriter.PartitionFileName);
        string edgesPath     = csv.File(SpcCsvWriter.EquilibriumEdgesFileName);

        // Prepare the energy ladder once; one closed-form pass yields both the
        // per-T affinity (the edge currency) and the per-T cluster count.
        PKWangContext context = PKWang.Prepare(graph, weightKind, field, symmetrization);

        int tCount = temperatures.Count;
        var affinities = new Affinities[tCount];
        var counts = new int[tCount];
        for (int i = 0; i < tCount; i++)
        {
            affinities[i] = PKWang.Solve(context, temperatures[i]);
            counts[i] = PKWang.Cluster(context, temperatures[i], theta).Count;
        }

        // Solver-native sweep diagnostic (caller order): per-T cluster count.
        WriteSweepCounts(sweepPath, temperatures, counts);

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
        Dendrogram dendrogram = ThermalDendrogram.FromEdgeCurves(graph, ascTemps, edgeColumns, theta);
        Landscape landscape = Landscape.Create(
            axis: "temperature",
            grid: ascTemps,
            valuesByGridPoint: bondColumns,
            provenance: new LandscapeProvenance("BondMass", GraphId: "pkwang", GaugeNote: "pkwang:closed-form"));

        ThermalEomResult eom = ThermalEom.Resolve(graph, dendrogram, landscape, minClusterSize, completion);
        Assignment partition = eom.Assignment;
        int abstained = CountAbstained(partition.Labels);

        // Representative single-T slice for the edge snapshot: the coldest
        // (most-coupled) column the periphery ascent keys on, unless overridden.
        int repIdx = cutTemperature is double cut ? NearestIndex(temperatures, cut) : order[0];
        double repTemperature = temperatures[repIdx];
        Affinities repAffinities = affinities[repIdx];

        // Currency-typed exports — identical writers/formats as the SW path.
        SpcCsvWriter.WritePartition(
            partition, partitionPath, features: dataset.Features, trueLabels: dataset.Labels);
        string? edgesWritten = SpcCsvWriter.WriteEquilibriumEdges(
            graph, repAffinities, alignments: null, edgesPath);

        // Evaluator scores on the resolved cut (method-agnostic indices —
        // IGraphPartitionEvaluator explicitly scores PKWang's affinity currency).
        var evaluatorScores = new Dictionary<string, double>(StringComparer.Ordinal);
        if (referenceLabels is not null && referenceLabels.Length == partition.Labels.Length)
            foreach (var ev in externalEvaluators ?? Array.Empty<IExternalClusterEvaluator>())
                evaluatorScores[ev.Name] = ev.Evaluate(partition.Labels, referenceLabels);
        foreach (var ev in spcEvaluators ?? Array.Empty<IGraphPartitionEvaluator>())
            evaluatorScores[ev.Name] = ev.Evaluate(
                graph, repAffinities.G, partition.Labels, partition.Count);

        // Summary JSON — PKWang-native shape. JsonObject keeps it reflection- and
        // source-gen-free (no dependency on the SW summary's serializer context).
        var scoresNode = new JsonObject();
        foreach (var kv in evaluatorScores) scoresNode[kv.Key] = kv.Value;
        var summary = new JsonObject
        {
            ["algorithm"]              = "spc-pkwang",
            ["resolver"]               = "thermal-eom",
            ["field"]                  = field.ToString(),
            ["symmetrization"]         = symmetrization.ToString(),
            ["theta"]                  = theta,
            ["minClusterSize"]         = minClusterSize,
            ["periphery"]              = completion.ToString(),
            ["nodeCount"]              = graph.NodeCount,
            ["temperatureCount"]       = temperatures.Count,
            ["clusterCount"]           = partition.Count,
            ["abstained"]              = abstained,
            ["representativeTemperature"] = repTemperature,
            ["representativeBy"]       = cutTemperature is null ? "coldest-column" : "explicit",
            ["evaluatorScores"]        = scoresNode,
            ["runDirectory"]           = runDirectory,
            ["sweepCsv"]               = sweepPath,
            ["partitionCsv"]           = partitionPath,
            ["equilibriumEdgesCsv"]    = edgesWritten,
        };
        string summaryPath = paths.Summary;
        File.WriteAllText(summaryPath, summary.ToJsonString(s_jsonOptions));

        return new PKWangUserRunResult(
            RunDirectory:            runDirectory,
            SweepCsvPath:            sweepPath,
            PartitionCsvPath:        partitionPath,
            EquilibriumEdgesCsvPath: edgesWritten,
            SummaryJsonPath:         summaryPath,
            ClusterCount:            partition.Count,
            Abstained:               abstained,
            RepresentativeTemperature: repTemperature,
            EvaluatorScores:         evaluatorScores);
    }

    private static readonly System.Text.Json.JsonSerializerOptions s_jsonOptions = new() { WriteIndented = true };

    private static void WriteSweepCounts(string path, IReadOnlyList<double> temps, int[] counts)
    {
        var lines = new List<string>(temps.Count + 1) { "temperature,cluster_count" };
        for (int i = 0; i < temps.Count; i++)
            lines.Add(FormattableString.Invariant($"{temps[i]:G9},{counts[i]}"));
        File.WriteAllLines(path, lines);
    }

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

/// <summary>Paths and headline figures from one <see cref="PKWangUserSession.Run"/> call.</summary>
public sealed record PKWangUserRunResult(
    string RunDirectory,
    string SweepCsvPath,
    string PartitionCsvPath,
    string? EquilibriumEdgesCsvPath,
    string SummaryJsonPath,
    int ClusterCount,
    int Abstained,
    double RepresentativeTemperature,
    IReadOnlyDictionary<string, double> EvaluatorScores);
