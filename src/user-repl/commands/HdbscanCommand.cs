using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using Clustering.Evaluation.External;
using Clustering.Graphical.HdbScan;

namespace UserRepl.Commands;

/// <summary>
/// The <c>hdbscan</c> subcommand. Runs HDBSCAN on the same dataset
/// shapes the SPC subcommand accepts (synthetic generators or CSV) and
/// writes a labels CSV + evaluator-score summary so output can sit
/// alongside SPC runs for like-for-like comparison.
/// </summary>
public static class HdbscanCommand
{
    private static readonly IReadOnlyList<string> DatasetKinds =
        SpcUserSession.ListAvailableSyntheticGeneratorNames();

    public static int Run(string[] args)
    {
        try
        {
            var options = ParseArguments(args);
            if (options.ShowHelp) { PrintHelp(); return 0; }

            string baseDirectory = options.BaseDirectory ?? Path.Combine("artifacts", "hdbscan-user");
            Directory.CreateDirectory(baseDirectory);

            SpcUserDataset dataset = LoadDataset(options);

            string runName = options.RunName ?? options.DatasetKind;
            string runDirectory = CreateRunDirectory(baseDirectory, runName, includeGuid: !options.NoGuid);

            // Manifest is the single source of truth for run provenance —
            // landed before the algorithm runs so even an interrupted
            // session leaves a self-describing directory.
            BuildManifest(options, runDirectory, args).WriteTo(runDirectory);

            Console.WriteLine($"Run directory      : {runDirectory}");
            Console.WriteLine(
                $"HDBSCAN            : minPts={options.MinPts}, " +
                $"minClusterSize={(options.MinClusterSize?.ToString() ?? "(=minPts)")}, " +
                $"allowSingleCluster={options.AllowSingleCluster}, " +
                $"metric={options.DistanceMetricSpec}");

            // HdbscanSession owns flatten → core-distance clamp → struct-metric
            // dispatch → noise-remap + evaluator scoring — the orchestration that
            // used to be hand-rolled here, in MapperHdbScan, and in the smoke
            // driver. The 6 external evaluators run only when GT labels exist.
            var settings = new HdbscanSettings
            {
                MinPts             = options.MinPts,
                MinClusterSize     = options.MinClusterSize,
                AllowSingleCluster = options.AllowSingleCluster,
                Metric             = options.DistanceMetricSpec,
            };

            IExternalClusterEvaluator[]? evaluators = dataset.Labels is { Length: > 0 }
                ? new IExternalClusterEvaluator[]
                {
                    new Purity(),
                    new NormalizedMutualInformation(),
                    new AdjustedRandIndex(),
                    new Homogeneity(),
                    new Completeness(),
                    new VMeasure(),
                }
                : null;

            HdbscanSessionResult session = HdbscanSession.Run(
                dataset.Features, settings, evaluators, dataset.Labels);

            HdbscanResult result = session.Result;
            int noiseCount = session.NoiseCount;
            IReadOnlyDictionary<string, double> evaluatorScores = session.EvaluatorScores;

            Console.WriteLine($"Result             : K={result.ClusterCount}, noise={noiseCount}");

            string partitionPath  = WritePartitionCsv(runDirectory, dataset, result);
            string dendrogramPath = WriteDendrogramJson(runDirectory, result);
            var clusterStats = ComputeClusterStats(result);
            string summaryPath = WriteSummaryJson(
                runDirectory, dataset, options, result, noiseCount, evaluatorScores,
                clusterStats, partitionPath, dendrogramPath);

            Console.WriteLine($"Partition CSV      : {partitionPath}");
            Console.WriteLine($"Dendrogram JSON    : {dendrogramPath}");
            Console.WriteLine($"Summary JSON       : {summaryPath}");
            foreach (var kvp in evaluatorScores)
                Console.WriteLine($"  {kvp.Key,-30} {kvp.Value.ToString("F6", CultureInfo.InvariantCulture)}");

            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }
    }

    // ── Manifest construction ────────────────────────────────────────────
    private static RunManifest BuildManifest(Options options, string runDirectory, string[] args)
    {
        DatasetSpec dataset = !string.IsNullOrWhiteSpace(options.DatasetFile)
            ? new DatasetSpec(
                Source:              "csv",
                GeneratorName:       null,
                GeneratorParameters: null,
                Seed:                null,
                CsvPath:             Path.GetFullPath(options.DatasetFile),
                LabelColumn:         options.LabelColumn,
                HasHeader:           options.HasHeader,
                Delimiter:           options.Delimiter == '\t' ? "tab" : options.Delimiter.ToString())
            : new DatasetSpec(
                Source:              "synthetic",
                GeneratorName:       options.DatasetKind,
                GeneratorParameters: new Dictionary<string, object?>(options.GeneratorParameters, StringComparer.OrdinalIgnoreCase),
                Seed:                options.Seed,
                CsvPath:             null,
                LabelColumn:         null,
                HasHeader:           null,
                Delimiter:           null);

        var hdbscan = new HdbscanSpec(
            MinPts:             options.MinPts,
            MinClusterSize:     options.MinClusterSize,
            AllowSingleCluster: options.AllowSingleCluster,
            DistanceMetric:     options.DistanceMetricSpec);

        return new RunManifest(
            SchemaVersion: RunManifest.CurrentSchemaVersion,
            CreatedUtc:    DateTime.UtcNow,
            Algorithm:     "hdbscan",
            CommandLine:   string.Join(" ", args),
            Dataset:       dataset,
            Graph:         null,
            Sweep:         null,
            Hdbscan:       hdbscan,
            Output:        new OutputSpec(RunDirectory: runDirectory, CheckpointDirectory: null));
    }

    // ── Dataset loading (shared shape with SpcCommand) ───────────────────
    private static SpcUserDataset LoadDataset(Options options)
    {
        if (!string.IsNullOrWhiteSpace(options.DatasetFile))
        {
            Console.WriteLine($"Loading dataset    : CSV {options.DatasetFile}");
            return SpcUserSession.FromCsv(
                options.DatasetFile, options.LabelColumn, options.HasHeader, options.Delimiter);
        }

        Console.WriteLine($"Loading dataset    : synthetic '{options.DatasetKind}'");
        var parameters = new Dictionary<string, object?>(StringComparer.OrdinalIgnoreCase)
        {
            ["seed"] = options.Seed,
        };
        foreach (var kvp in options.GeneratorParameters)
            parameters[kvp.Key] = kvp.Value;
        return SpcUserSession.GenerateDataset(options.DatasetKind, parameters);
    }

    private static string WritePartitionCsv(string runDirectory, SpcUserDataset dataset, HdbscanResult result)
    {
        string path = Path.Combine(runDirectory, "hdbscan_partition.csv");
        int n   = dataset.Features.Length;
        int dim = dataset.Features[0].Length;
        bool hasGt = dataset.Labels is { Length: > 0 } && dataset.Labels.Length == n;

        using var sw = new StreamWriter(path);
        // Header: feature_0..feature_{dim-1}, label, membership_probability, true_label?
        var headerCols = new List<string>(dim + 3);
        for (int d = 0; d < dim; d++) headerCols.Add($"feature_{d}");
        headerCols.Add("label");
        headerCols.Add("membership_probability");
        if (hasGt) headerCols.Add("true_label");
        sw.WriteLine(string.Join(",", headerCols));

        var ci = CultureInfo.InvariantCulture;
        for (int i = 0; i < n; i++)
        {
            var row = new List<string>(dim + 3);
            for (int d = 0; d < dim; d++)
                row.Add(dataset.Features[i][d].ToString("G17", ci));
            row.Add(result.Labels[i].ToString(ci));
            row.Add(result.MembershipProbabilities[i].ToString("G17", ci));
            if (hasGt) row.Add(dataset.Labels[i].ToString(ci));
            sw.WriteLine(string.Join(",", row));
        }
        return path;
    }

    private static string WriteSummaryJson(
        string runDirectory, SpcUserDataset dataset, Options options,
        HdbscanResult result, int noiseCount,
        IReadOnlyDictionary<string, double> evaluatorScores,
        IReadOnlyList<HdbscanClusterStat> clusterStats,
        string partitionPath, string dendrogramPath)
    {
        string path = Path.Combine(runDirectory, "summary.json");
        var payload = new HdbscanSummaryPayload(
            Algorithm: "HDBSCAN",
            Dataset: dataset.Metadata,
            Hdbscan: new HdbscanRunPayload(
                MinPts: options.MinPts,
                MinClusterSize: options.MinClusterSize,
                AllowSingleCluster: options.AllowSingleCluster,
                Metric: options.DistanceMetricSpec),
            ReferenceLabels: dataset.Labels,
            Result: new HdbscanResultPayload(
                ClusterCount: result.ClusterCount,
                NoiseCount: noiseCount,
                EvaluatorScores: evaluatorScores,
                Clusters: clusterStats),
            Run: new HdbscanRunPaths(
                RunDirectory: runDirectory,
                PartitionCsv: partitionPath,
                DendrogramJson: dendrogramPath));

        UserReplJson.Writer.WriteDocumentToFile(payload, path);
        return path;
    }

    /// <summary>
    /// Persist the single-linkage dendrogram as a JSON sidecar so
    /// downstream plotting (Python dendrogram renders, matplotlib /
    /// seaborn) can read the merge tree directly. Mutual-reachability
    /// distance is the cost axis; λ = 1/distance is the persistence
    /// scalar HDBSCAN's condensation pass consumed to produce the
    /// labels.
    /// </summary>
    private static string WriteDendrogramJson(string runDirectory, HdbscanResult result)
    {
        string path = Path.Combine(runDirectory, "hdbscan_dendrogram.json");
        var payload = new HdbscanDendrogramPayload(
            LeafCount: result.Dendrogram.LeafCount,
            CostAxis: result.Dendrogram.CostAxis,
            Merges: result.Dendrogram.Merges.Select(m => new HdbscanDendrogramMerge(
                LeftChild: m.LeftChild,
                RightChild: m.RightChild,
                Distance: m.Distance,
                Size: m.Size,
                Lambda: m.Distance > 0.0 ? 1.0 / m.Distance : double.PositiveInfinity)).ToArray());

        UserReplJson.Writer.WriteDocumentToFile(payload, path);
        return path;
    }

    /// <summary>
    /// Per-cluster size + mean membership probability. Cheap to derive
    /// from <see cref="HdbscanResult.Labels"/> + <see cref="HdbscanResult.MembershipProbabilities"/>;
    /// gives the user a quick "is this cluster well-formed?" read in
    /// summary.json without opening the per-point CSV.
    /// </summary>
    private static IReadOnlyList<HdbscanClusterStat> ComputeClusterStats(HdbscanResult result)
    {
        int k = result.ClusterCount;
        if (k <= 0) return Array.Empty<HdbscanClusterStat>();

        var sizes = new int[k];
        var probSums = new double[k];
        for (int i = 0; i < result.Labels.Length; i++)
        {
            int label = result.Labels[i];
            if (label < 0) continue;   // noise
            sizes[label]++;
            probSums[label] += result.MembershipProbabilities[i];
        }

        var stats = new HdbscanClusterStat[k];
        for (int c = 0; c < k; c++)
        {
            double mean = sizes[c] > 0 ? probSums[c] / sizes[c] : 0.0;
            stats[c] = new HdbscanClusterStat(
                ClusterId:              c,
                Size:                   sizes[c],
                MeanMembershipProbability: mean);
        }
        return stats;
    }

    private static string CreateRunDirectory(string baseDirectory, string runName, bool includeGuid)
    {
        string folder = includeGuid
            ? $"{runName}-{Guid.NewGuid():N}"
            : runName;
        string runDirectory = Path.Combine(baseDirectory, folder);
        Directory.CreateDirectory(runDirectory);
        return runDirectory;
    }

    // ── Argv parsing ─────────────────────────────────────────────────────
    private static Options ParseArguments(string[] args)
    {
        var options = new Options();

        // Pre-pass: load JSON presets first so CLI flags in the same argv override them.
        for (int i = 0; i < args.Length; i++)
        {
            if (args[i].Equals("--config", StringComparison.OrdinalIgnoreCase) && i + 1 < args.Length)
            {
                HdbscanPreset.Load(args[i + 1]).ApplyTo(options);
            }
        }

        for (int i = 0; i < args.Length; i++)
        {
            string arg = args[i];
            string? next = i + 1 < args.Length ? args[i + 1] : null;
            switch (arg.ToLowerInvariant())
            {
                case "--help": case "-h": options.ShowHelp = true; return options;
                case "--config": i++; break;   // already consumed in pre-pass

                case "--dataset": options.DatasetKind = ParseDatasetKind(RequireValue(next, arg)); i++; break;
                case "--param": ParseGeneratorParameter(RequireValue(next, arg), options.GeneratorParameters); i++; break;
                case "--dataset-file": options.DatasetFile = RequireValue(next, arg); i++; break;
                case "--label-column": options.LabelColumn = RequireValue(next, arg); i++; break;
                case "--delimiter": options.Delimiter = ParseDelimiter(RequireValue(next, arg)); i++; break;
                case "--no-header": options.HasHeader = false; break;

                case "--base-dir": options.BaseDirectory = RequireValue(next, arg); i++; break;
                case "--run-name": options.RunName = RequireValue(next, arg); i++; break;
                case "--no-guid": options.NoGuid = true; break;
                case "--seed": options.Seed = int.Parse(RequireValue(next, arg), CultureInfo.InvariantCulture); i++; break;

                case "--min-pts": options.MinPts = int.Parse(RequireValue(next, arg), CultureInfo.InvariantCulture); i++; break;
                case "--min-cluster-size":
                    options.MinClusterSize = int.Parse(RequireValue(next, arg), CultureInfo.InvariantCulture); i++; break;
                case "--allow-single-cluster": options.AllowSingleCluster = true; break;
                case "--no-allow-single-cluster": options.AllowSingleCluster = false; break;
                case "--distance-metric": options.DistanceMetricSpec = RequireValue(next, arg); i++; break;

                default:
                    throw new ArgumentException($"Unknown argument '{arg}'. Use --help for usage.");
            }
        }
        return options;

        static string RequireValue(string? value, string option)
        {
            if (string.IsNullOrWhiteSpace(value))
                throw new ArgumentException($"Option '{option}' requires a value.");
            return value;
        }
    }

    internal static string ParseDatasetKind(string value)
    {
        var match = DatasetKinds.FirstOrDefault(name => string.Equals(name, value, StringComparison.OrdinalIgnoreCase));
        if (match is not null) return match;
        throw new ArgumentException($"Unknown dataset kind '{value}'. Valid: {string.Join(", ", DatasetKinds)}.");
    }

    internal static char ParseDelimiter(string value) => value.Trim() switch
    {
        "tab" or "\\t" => '\t',
        var s when s.Length == 1 => s[0],
        _ => throw new ArgumentException($"Invalid delimiter '{value}'. Use a single character or 'tab'."),
    };

    private static void ParseGeneratorParameter(string token, IDictionary<string, object?> parameters)
    {
        int eq = token.IndexOf('=');
        if (eq <= 0) throw new ArgumentException($"Generator parameter must be name=value: '{token}'");
        parameters[token[..eq].Trim()] = token[(eq + 1)..].Trim();
    }

    private static void PrintHelp()
    {
        Console.WriteLine("Usage: userrepl hdbscan --dataset <generator> [options]");
        Console.WriteLine("       userrepl hdbscan --dataset-file <path> [options]");
        Console.WriteLine("       userrepl hdbscan --config <preset.json> [options]");
        Console.WriteLine();
        Console.WriteLine("Presets:");
        Console.WriteLine("  --config <path>               Load a JSON preset. CLI flags passed alongside override it.");
        Console.WriteLine();
        Console.WriteLine("Dataset:");
        Console.WriteLine("  --dataset <name>              Synthetic generator (use 'userrepl spc --list-generators' to see options)");
        Console.WriteLine("  --param <name>=<value>        Per-generator parameter; may repeat");
        Console.WriteLine("  --dataset-file <path>         CSV input");
        Console.WriteLine("  --label-column <name|idx>     Label column (default: last)");
        Console.WriteLine("  --delimiter <char|tab>        CSV separator (default: ,)");
        Console.WriteLine("  --no-header                   Treat CSV as having no header row");
        Console.WriteLine("  --seed <int>                  Synthetic generator seed (default: 42)");
        Console.WriteLine();
        Console.WriteLine("Output:");
        Console.WriteLine("  --base-dir <path>             Run-root base directory (default: artifacts/hdbscan-user)");
        Console.WriteLine("  --run-name <name>             Subdirectory name (default: dataset kind)");
        Console.WriteLine("  --no-guid                     Don't append a GUID to the run folder");
        Console.WriteLine();
        Console.WriteLine("HDBSCAN:");
        Console.WriteLine("  --min-pts <int>               Core-distance neighbor count (default: 5)");
        Console.WriteLine("  --min-cluster-size <int>      Min subtree size; defaults to --min-pts when absent");
        Console.WriteLine("  --allow-single-cluster        Permit selecting the root cluster (default: on)");
        Console.WriteLine("  --no-allow-single-cluster     Disallow single-cluster selection (sklearn-default behavior)");
        Console.WriteLine("  --distance-metric <spec>      euclidean | manhattan | minkowski:p=N | hamming | poincare | cosine");
        Console.WriteLine("                                Default: euclidean. Minkowski exponent folds into the spec, e.g. minkowski:p=1.5");
        Console.WriteLine();
        Console.WriteLine("  --help, -h                    Show this help text");
    }

    // ── Options bag ──────────────────────────────────────────────────────
    internal sealed class Options
    {
        public bool ShowHelp { get; set; }

        // Dataset
        public string DatasetKind { get; set; } = DatasetKinds.FirstOrDefault() ?? "unknown";
        public string? DatasetFile { get; set; }
        public string? LabelColumn { get; set; }
        public bool HasHeader { get; set; } = true;
        public char Delimiter { get; set; } = ',';
        public Dictionary<string, object?> GeneratorParameters { get; } = new(StringComparer.OrdinalIgnoreCase);

        // Output
        public string? BaseDirectory { get; set; }
        public string? RunName { get; set; }
        public bool NoGuid { get; set; }
        public int Seed { get; set; } = 42;

        // HDBSCAN
        public int MinPts { get; set; } = 5;
        public int? MinClusterSize { get; set; }
        public bool AllowSingleCluster { get; set; } = true;
        public string DistanceMetricSpec { get; set; } = "euclidean";
    }
}
