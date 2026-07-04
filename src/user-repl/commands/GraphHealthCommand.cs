using System;
using System.Globalization;
using System.IO;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;

namespace UserRepl.Commands;

/// <summary>
/// The <c>graph-health</c> subcommand. Re-evaluates a run's graph
/// health from its manifest (without touching the SPC checkpoints) and
/// refreshes <c>&lt;run-dir&gt;/graph_health.json</c>. Useful for
/// ad-hoc inspection of an existing run without re-rebuilding sweep
/// CSVs or paying the extract command's hydration cost.
/// </summary>
/// <remarks>
/// <para>Subset of <see cref="ExtractCommand"/>: same manifest read +
/// graph rebuild path, but stops after the diagnostic — no checkpoint
/// walk, no profile reconstruction. Output is the same
/// <c>graph_health.json</c> file the other commands write, so
/// downstream consumers can rely on a single file location regardless
/// of which command produced it.</para>
/// </remarks>
public static class GraphHealthCommand
{
    public static int Run(string[] args)
    {
        try
        {
            var options = ParseArguments(args);
            if (options.ShowHelp) { PrintHelp(); return 0; }

            string runDirectory = options.RunDirectory
                ?? throw new ArgumentException("--run-dir is required.");

            RunManifest? manifest = TryLoadManifest(runDirectory);
            SpcUserDataset dataset = LoadDataset(options, manifest);
            GraphBuildResult build = BuildGraphResult(options, manifest, dataset);

            int k = manifest?.Graph?.K ?? options.K;
            GraphHealthReport health = GraphHealth.Evaluate(build, k);
            string path = GraphHealthFile.WriteTo(runDirectory, health);

            Console.WriteLine($"Run directory      : {runDirectory}");
            Console.WriteLine($"Manifest           : {(manifest is null ? "(not found — using CLI flags)" : "loaded")}");
            Console.WriteLine($"Graph health JSON  : {path}");
            Console.WriteLine();
            if (!string.IsNullOrEmpty(health.Verdict.PrimaryRecommendation))
                Console.WriteLine($"Recommendation     : {health.Verdict.PrimaryRecommendation}");
            else
                Console.WriteLine("Recommendation     : (none — verdict flags are all clear)");

            PrintFlagSummary(health.Verdict);
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }
    }

    private static void PrintFlagSummary(GraphHealthVerdict v)
    {
        Console.WriteLine();
        Console.WriteLine("Flags:");
        Console.WriteLine($"  bandwidth_too_large   : {v.BandwidthTooLarge}");
        Console.WriteLine($"  bandwidth_too_small   : {v.BandwidthTooSmall}");
        Console.WriteLine($"  hubness_concern       : {v.HubnessConcern}");
        Console.WriteLine($"  connectivity_concern  : {v.ConnectivityConcern}");
        Console.WriteLine($"  forced_edges_concern  : {v.ForcedEdgesConcern}");
        Console.WriteLine($"  underconnected_nodes  : {v.UnderconnectedNodes}");
    }

    private static RunManifest? TryLoadManifest(string runDirectory)
    {
        string path = RunManifest.PathFor(runDirectory);
        return File.Exists(path) ? RunManifest.ReadFrom(runDirectory) : null;
    }

    private static SpcUserDataset LoadDataset(Options options, RunManifest? manifest)
    {
        if (!string.IsNullOrWhiteSpace(options.DatasetFile))
        {
            return SpcUserSession.FromCsv(
                options.DatasetFile, options.LabelColumn, options.HasHeader, options.Delimiter);
        }

        if (manifest is null)
            throw new ArgumentException(
                "No manifest.json in run dir, so dataset reconstruction needs --dataset-file.");

        return manifest.Dataset.Materialize();
    }

    private static GraphBuildResult BuildGraphResult(Options options, RunManifest? manifest, SpcUserDataset dataset)
    {
        if (manifest?.Graph is { } graphSpec && !options.HasGraphOverride)
        {
            IDistanceMetric? metric = string.IsNullOrWhiteSpace(graphSpec.DistanceMetric)
                ? null
                : DistanceMetricFactory.Create(graphSpec.DistanceMetric);
            return SpcGraphBuilder.BuildResult(
                dataset.Features, graphSpec.ToConfig(), metric,
                protectedEdges: TDA.Ph.H1CycleEdges.FromDistanceGraph);
        }

        var config = new GraphCompilerConfig
        {
            Topology = new TopologyConfig
            {
                Kind = TopologyKind.Knn,
                K = options.K,
            },
            Filter = new FilterConfig
            {
                Kind = FilterKind.OrRule,
            },
            Repair = new RepairConfig
            {
                Kind = options.EnsureConnected ? RepairKind.MstMin : RepairKind.NoRepair,
            },
            Projection = new CouplingProjection
            {
                Kernel = new Gaussian(1.0),
            },
        };
        return SpcGraphBuilder.BuildResult(
            dataset.Features, config,
            protectedEdges: TDA.Ph.H1CycleEdges.FromDistanceGraph);
    }

    private static Options ParseArguments(string[] args)
    {
        var options = new Options();
        for (int i = 0; i < args.Length; i++)
        {
            string arg = args[i];
            string? next = i + 1 < args.Length ? args[i + 1] : null;
            switch (arg.ToLowerInvariant())
            {
                case "--help": case "-h": options.ShowHelp = true; return options;
                case "--run-dir": options.RunDirectory = RequireValue(next, arg); i++; break;
                case "--dataset-file": options.DatasetFile = RequireValue(next, arg); i++; break;
                case "--label-column": options.LabelColumn = RequireValue(next, arg); i++; break;
                case "--delimiter": options.Delimiter = ParseDelimiter(RequireValue(next, arg)); i++; break;
                case "--no-header": options.HasHeader = false; break;
                case "--k": options.K = int.Parse(RequireValue(next, arg), CultureInfo.InvariantCulture); options.HasGraphOverride = true; i++; break;
                case "--ensure-connected": options.EnsureConnected = true; options.HasGraphOverride = true; break;
                default: throw new ArgumentException($"Unknown argument '{arg}'. Use --help for usage.");
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

    private static char ParseDelimiter(string value) => value.Trim() switch
    {
        "tab" or "\\t" => '\t',
        var s when s.Length == 1 => s[0],
        _ => throw new ArgumentException($"Invalid delimiter '{value}'. Use a single character or 'tab'."),
    };

    private static void PrintHelp()
    {
        Console.WriteLine("Usage: userrepl graph-health --run-dir <path> [options]");
        Console.WriteLine();
        Console.WriteLine("Re-evaluates graph health from <run-dir>/manifest.json (or CLI overrides)");
        Console.WriteLine("and refreshes <run-dir>/graph_health.json. Subset of 'userrepl extract' —");
        Console.WriteLine("does NOT walk the checkpoint directory or rebuild sweep CSVs.");
        Console.WriteLine();
        Console.WriteLine("Required:");
        Console.WriteLine("  --run-dir <path>              The SPC run directory");
        Console.WriteLine();
        Console.WriteLine("Optional (overrides manifest.json):");
        Console.WriteLine("  --dataset-file <csv>          Use a different dataset");
        Console.WriteLine("  --label-column <name|idx>     CSV label column");
        Console.WriteLine("  --delimiter <char|tab>        CSV separator (default: ,)");
        Console.WriteLine("  --no-header                   CSV has no header");
        Console.WriteLine("  --k <int>                     KNN k for graph rebuild");
        Console.WriteLine("  --ensure-connected            Apply MST-repair on graph rebuild");
        Console.WriteLine();
        Console.WriteLine("  --help, -h                    Show this help");
    }

    internal sealed class Options
    {
        public bool ShowHelp { get; set; }
        public string? RunDirectory { get; set; }
        public string? DatasetFile { get; set; }
        public string? LabelColumn { get; set; }
        public bool HasHeader { get; set; } = true;
        public char Delimiter { get; set; } = ',';
        public int K { get; set; } = 10;
        public bool EnsureConnected { get; set; }
        public bool HasGraphOverride { get; set; }
    }
}
