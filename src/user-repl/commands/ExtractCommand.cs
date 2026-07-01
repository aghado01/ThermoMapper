using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Profiling.Signals;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Execution.Sinks;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;
using Graphs.Primitives;

namespace UserRepl.Commands;

/// <summary>
/// The <c>extract</c> subcommand. Walks an existing checkpoint
/// directory (the output of one or more interrupted/partial SPC sweeps),
/// hydrates the available <c>.spcx</c>/<c>.spce</c> files into
/// <see cref="SpcRunResult"/>s, and rebuilds the standard SPC CSVs
/// (<c>spc_sweep.csv</c>, <c>spc_criteria.csv</c>) without re-running
/// the sampler.
/// </summary>
/// <remarks>
/// <para><b>Manifest-driven.</b> Reads <c>&lt;run-dir&gt;/manifest.json</c>
/// to reconstruct the original dataset + graph + checkpoint location.
/// No flags are needed beyond <c>--run-dir</c> for a typical case;
/// dataset/graph overrides exist for manifests that pre-date this
/// feature or for ad-hoc reanalysis with a different graph.</para>
///
/// <para><b>Partial sweeps are fine.</b> Missing checkpoints just shrink
/// the assembled profile — there's no requirement that every (T, replica)
/// originally requested be present. The output is "the best SweepProfile
/// I can build from what's on disk."</para>
/// </remarks>
public static class ExtractCommand
{
    public static int Run(string[] args)
    {
        try
        {
            var options = ParseArguments(args);
            if (options.ShowHelp) { PrintHelp(); return 0; }

            string runDirectory = options.RunDirectory
                ?? throw new ArgumentException("--run-dir is required.");

            // Prefer manifest.json for dataset + graph reconstruction; fall
            // back to CLI overrides if the manifest is absent or any flag
            // explicitly overrides it.
            RunManifest? manifest = TryLoadManifest(runDirectory);

            SpcUserDataset   dataset = LoadDataset(options, manifest);
            GraphBuildResult build   = BuildGraphResult(options, manifest, dataset);
            CsrGraph         graph   = build.Graph;

            string checkpointDirectory = options.CheckpointDir
                ?? manifest?.Output.CheckpointDirectory
                ?? Path.Combine(runDirectory, "checkpoints");

            if (!Directory.Exists(checkpointDirectory))
                throw new DirectoryNotFoundException(
                    $"Checkpoint directory not found: {checkpointDirectory}");

            Console.WriteLine($"Run directory      : {runDirectory}");
            Console.WriteLine($"Manifest           : {(manifest is null ? "(not found — using CLI flags)" : "loaded")}");
            Console.WriteLine($"Checkpoint dir     : {checkpointDirectory}");

            List<SpcRunResult> hydrated = HydrateAllCheckpoints(graph, checkpointDirectory);
            Console.WriteLine($"Hydrated runs      : {hydrated.Count}");
            if (hydrated.Count == 0)
            {
                Console.Error.WriteLine("No .spcx files found under checkpoint dir; nothing to extract.");
                return 1;
            }

            SweepProfile profile = SweepProfile.From(hydrated);
            ISignalAnalyzer analyzer = BuildAnalyzer(manifest?.Sweep?.Analyzer);
            ProfileCriteria criteria = analyzer.Analyze(profile);

            string sweepPath    = SpcOutputPathHelper.GetSweepCsvPath(runDirectory);
            string criteriaPath = SpcOutputPathHelper.GetCriteriaCsvPath(runDirectory);

            SpcCsvWriter.WriteSweepProfile(profile, sweepPath);
            SpcCsvWriter.WriteCriteria(criteria, criteriaPath);

            // Persist the per-(T, replica) trace alongside the averaged sweep
            // CSV — partial extracts get partial traces, which is the right
            // shape for inspecting in-flight runs in Python.
            string replicaTracesPath = SpcOutputPathHelper.GetReplicaTracesCsvPath(runDirectory);
            SpcCsvWriter.WriteReplicaTraces(hydrated, replicaTracesPath);

            // Edge dump anchored on the PARTITION TEMPERATURE — re-derived from the
            // rebuilt profile exactly as the live run picks it (PickPeakTemperature),
            // not inferred from a cycle-count proxy. The result the partition is taken
            // at is the edge-bearing one at that T (cycle count breaks ties among
            // co-located runs). If no edge-bearing result exists at the partition T —
            // a pure scan, or that T was never run at the edge tier — final-edges is
            // inapplicable by construction, so we skip it rather than guess a "final".
            double partitionT = SpcProfileAnalysis.PickPeakTemperature(profile);
            SpcRunResult? partitionResult = hydrated
                .Where(r => r.Accumulator.BondFormedCount is not null
                         && Math.Abs(r.Accumulator.Temperature - partitionT) < 1e-9)
                .OrderByDescending(r => r.Accumulator.DrawCount)
                .FirstOrDefault();
            string equilibriumEdgesPath = SpcOutputPathHelper.GetEquilibriumEdgesCsvPath(runDirectory);
            string? equilibriumEdgesWritten = null;
            if (partitionResult is not null)
            {
                Affinities   affinities    = SwCurrencies.ToAffinities(partitionResult.Accumulator);
                Alignments?  alignments    = partitionResult.Accumulator.SpinAgreementCount is null
                    ? null : SwCurrencies.ToAlignments(partitionResult.Accumulator);
                CoMembership? coMembership = partitionResult.Accumulator.CoMembershipCount is null
                    ? null : SwCurrencies.ToCoMembership(partitionResult.Accumulator);
                equilibriumEdgesWritten = SpcCsvWriter.WriteEquilibriumEdges(
                    partitionResult.Graph, affinities, alignments, equilibriumEdgesPath,
                    coMembership: coMembership);
            }

            // Refresh graph_health.json from the reconstructed build so the
            // run directory carries a self-consistent diagnostic snapshot
            // alongside the rebuilt sweep CSVs. Same rule used by SpcCommand
            // at run start; extract just re-runs it now that the graph has
            // been rebuilt deterministically from the manifest.
            int k = manifest?.Graph?.K ?? options.K;
            GraphHealthReport health = GraphHealth.Evaluate(build, k);
            string graphHealthPath = GraphHealthFile.WriteTo(runDirectory, health);

            Console.WriteLine($"Sweep CSV          : {sweepPath}");
            Console.WriteLine($"Criteria CSV       : {criteriaPath}");
            Console.WriteLine($"Replica traces CSV : {replicaTracesPath}");
            if (equilibriumEdgesWritten is not null)
                Console.WriteLine($"Equilibrium edges  : {equilibriumEdgesWritten}");
            Console.WriteLine($"Graph health JSON  : {graphHealthPath}");
            if (!string.IsNullOrEmpty(health.Verdict.PrimaryRecommendation))
                Console.WriteLine($"Graph health       : {health.Verdict.PrimaryRecommendation}");
            Console.WriteLine($"Profile coverage   : {profile.Temperatures.Count} unique T's");
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }
    }

    // ── Manifest / dataset / graph reconstruction ────────────────────────
    private static RunManifest? TryLoadManifest(string runDirectory)
    {
        string path = RunManifest.PathFor(runDirectory);
        return File.Exists(path) ? RunManifest.ReadFrom(runDirectory) : null;
    }

    private static SpcUserDataset LoadDataset(Options options, RunManifest? manifest)
    {
        // Explicit --dataset-file always wins (allows ad-hoc reanalysis with a
        // different dataset against the same checkpoint shape).
        if (!string.IsNullOrWhiteSpace(options.DatasetFile))
        {
            return SpcUserSession.FromCsv(
                options.DatasetFile, options.LabelColumn, options.HasHeader, options.Delimiter);
        }

        if (manifest is null)
            throw new ArgumentException(
                "No manifest.json in run dir, so dataset reconstruction needs --dataset-file " +
                "(or re-run the original command after this version which writes a manifest).");

        return manifest.Dataset.Materialize();
    }

    private static GraphBuildResult BuildGraphResult(Options options, RunManifest? manifest, SpcUserDataset dataset)
    {
        // Manifest-only path when no overrides are passed.
        if (manifest?.Graph is { } graphSpec && !options.HasGraphOverride)
        {
            IDistanceMetric? metric = string.IsNullOrWhiteSpace(graphSpec.DistanceMetric)
                ? null
                : DistanceMetricFactory.Create(graphSpec.DistanceMetric);
            return SpcGraphBuilder.BuildResult(dataset.Features, graphSpec.ToConfig(), metric);
        }

        // Override path or pre-manifest fallback.
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
        return SpcGraphBuilder.BuildResult(dataset.Features, config);
    }

    // Rebuild the signal analyzer the original run recorded (manifest.Sweep.Analyzer
    // is the AnalyzerKind enum name) so the reconstructed criteria match the live run.
    // Falls back to χ-peak when there is no manifest or the value is unknown — including
    // old "MultiSignalConsensus" manifests, since that analyzer is parked (parking-lot/);
    // their reconstructed criteria are the χ-peak reading of the same profile.
    private static ISignalAnalyzer BuildAnalyzer(string? recorded) => new ChiPeakSignalAnalyzer();

    /// <summary>
    /// Walks the checkpoint directory recursively, loading every
    /// <c>.spcx</c> file (and its paired <c>.spce</c> when present)
    /// into an in-memory <see cref="SpcRunResult"/>. The sink's
    /// <see cref="IFrameSink.TryLoad"/> implementation handles the
    /// per-file serialization; we just enumerate and synthesize task
    /// specs to feed it.
    /// </summary>
    private static List<SpcRunResult> HydrateAllCheckpoints(CsrGraph graph, string checkpointDirectory)
    {
        var sink = SpcxDiskFrameSink.Instance;
        var results = new List<SpcRunResult>();
        var spcxFiles = Directory.GetFiles(checkpointDirectory, "*.spcx", SearchOption.AllDirectories);

        foreach (string spcxPath in spcxFiles)
        {
            // Probe whether a paired .spce exists; pick the accumulation that
            // lets TryLoad succeed — Currencies tier signals per-edge arrays.
            string spcePath = Path.ChangeExtension(spcxPath, ".spce");
            AccumulationSpec accumulation = File.Exists(spcePath)
                ? AccumulationSpec.Currencies
                : AccumulationSpec.None;

            // Most fields of the task spec are unused by TryLoad on the
            // disk sink; the only field that matters is CheckpointPath
            // (path on disk) and Accumulation (presence of per-edge arrays).
            var task = new SpcTaskSpec
            {
                Temperature    = 0.0,
                ReplicaIndex   = 0,
                Q              = 0,
                Accumulation   = accumulation,
                CheckpointPath = spcxPath,
            };

            SpcRunResult? r = sink.TryLoad(task, graph);
            if (r is not null) results.Add(r);
        }
        return results;
    }

    // ── Argv parsing ─────────────────────────────────────────────────────
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
                case "--checkpoint-dir": options.CheckpointDir = RequireValue(next, arg); i++; break;
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
        Console.WriteLine("Usage: userrepl extract --run-dir <path> [options]");
        Console.WriteLine();
        Console.WriteLine("Walks the checkpoint files (.spcx + .spce) in a run directory, hydrates");
        Console.WriteLine("them, and rebuilds spc_sweep.csv + spc_criteria.csv without re-running");
        Console.WriteLine("the sampler. Partial sweeps are fine — the output reflects only what's");
        Console.WriteLine("on disk.");
        Console.WriteLine();
        Console.WriteLine("By default, the run's manifest.json (written by 'userrepl spc' at run start)");
        Console.WriteLine("supplies the dataset + graph rebuild parameters. Pre-manifest runs or ad-hoc");
        Console.WriteLine("reanalysis can override via the flags below.");
        Console.WriteLine();
        Console.WriteLine("Required:");
        Console.WriteLine("  --run-dir <path>              The original SPC run directory");
        Console.WriteLine();
        Console.WriteLine("Optional (overrides manifest.json):");
        Console.WriteLine("  --checkpoint-dir <path>       Override checkpoint location");
        Console.WriteLine("  --dataset-file <csv>          Use a different dataset for graph rebuild");
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
        public string? CheckpointDir { get; set; }
        public string? DatasetFile { get; set; }
        public string? LabelColumn { get; set; }
        public bool HasHeader { get; set; } = true;
        public char Delimiter { get; set; } = ',';
        public int K { get; set; } = 10;
        public bool EnsureConnected { get; set; }

        /// <summary>True when the user passed any graph-override flag —
        /// suppresses the manifest-driven graph reconstruction so the
        /// override actually takes effect.</summary>
        public bool HasGraphOverride { get; set; }
    }
}
