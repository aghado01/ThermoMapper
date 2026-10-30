using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Dendrograms;
using Archivory;
using Clustering.Evaluation.External;
using Clustering.Evaluation.Internal;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Profiling.Signals;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;
using Graphs.Primitives;

using Graphs.Models.Potts;


namespace UserRepl.Commands;

/// <summary>
/// The <c>spc</c> subcommand. Runs an end-to-end SPC clustering session
/// against a synthetic-generator dataset or a CSV input, with adaptive
/// or fixed-grid temperature sweeps and optional checkpoint resume.
/// </summary>
public static class SpcCommand
{
    public static int Run(string[] args)
    {
        try
        {
            var options = ParseArguments(args);
            if (options.ShowHelp) { SpcCommandHelp.PrintHelp(); return 0; }
            if (options.ListGenerators) { SpcCommandHelp.PrintGeneratorList(); return 0; }
            if (!string.IsNullOrWhiteSpace(options.GeneratorSchema))
            {
                SpcCommandHelp.PrintGeneratorSchema(options.GeneratorSchema);
                return 0;
            }

            string baseDirectory = options.BaseDirectory ?? "artifacts";

            SpcUserDataset dataset = GraphRun.LoadDataset(options);
            var graphConfig = GraphRun.BuildGraphConfig(options);
            IDistanceMetric? distanceMetric = GraphRun.DistanceMetric(options);

            ISignalAnalyzer analyzer = BuildAnalyzer(options.AnalyzerKind);
            IPartitionStrategy partitionStrategy = BuildPartitionStrategy(
                options.PartitionStrategyKind, options.Theta, options.PeripheralCapture);
            IHierarchicalPartitionStrategy? hierarchicalStrategy = BuildHierarchicalStrategy(
                options.HierarchicalStrategyKind);

            // One scope owns the whole run tree (artifacts/{family}/{stamp}/), created
            // once and threaded down — no second owner, no double run-directory.
            // family: explicit --run-name, else the "run_name" placeholder until
            // RunIdentity auto-resolution lands.
            var runId = RunIdentity.Resolve(options.RunName, callerStub: "spc");
            var runPaths = new SpcRunPaths(
                ArtifactScope.Root(baseDirectory, runId.Family, RunStamp.Now()).EnsureDirectory());
            string runDirectory = runPaths.Scope.Dir;

            // Honor --checkpoint-dir as-is (shared resume dirs across runs); else the
            // run's own checkpoints sub-scope.
            if (options.Resolver != ResolverKind.None && !options.Accumulation.CoMembership)
            {
                throw new ArgumentException(
                    $"--resolver {ResolverToken(options.Resolver)} needs the co-membership currency: " +
                    "add --accumulation comembership (+ cluster-size-landscape for thermal-eom/hierarchy).");
            }
            // thermal-eom / hierarchy walk the per-node landscape; lineage does not.
            bool resolverNeedsLandscape =
                options.Resolver is ResolverKind.ThermalEom or ResolverKind.Hierarchy;
            if (resolverNeedsLandscape && !options.Accumulation.ClusterSizeLandscape)
            {
                throw new ArgumentException(
                    $"--resolver {ResolverToken(options.Resolver)} needs a rich sweep: " +
                    "add --accumulation comembership,cluster-size-landscape.");
            }

            string checkpointDirectory = options.CheckpointDir ?? runPaths.Checkpoints.Dir;

            // Persist the manifest before any compute happens so a crash mid-run
            // still leaves a self-describing run directory that ExtractCommand
            // (or any reanalysis) can pick up.
            var manifest = BuildManifest(options, graphConfig, runDirectory, checkpointDirectory, args, runId, null);
            manifest.WriteTo(runDirectory);

            Console.WriteLine($"Run directory      : {runDirectory}");
            Console.WriteLine($"Checkpoint dir     : {checkpointDirectory}");
            GraphRun.PrintGraphSummary(options, graphConfig);
            Console.WriteLine($"Schedule           : {options.Schedule}");
            Console.WriteLine($"Analyzer           : {analyzer.GetType().Name}");
            Console.WriteLine($"PartitionStrategy  : {partitionStrategy.GetType().Name} (theta={options.Theta})");
            if (hierarchicalStrategy is not null)
                Console.WriteLine($"HierarchicalStrategy: {hierarchicalStrategy.GetType().Name}");

            IExternalClusterEvaluator[] externalEvaluators =
            {
                new Purity(),
                new NormalizedMutualInformation(),
                new AdjustedRandIndex(),
                new Homogeneity(),
                new Completeness(),
                new VMeasure(),
            };
            IGraphPartitionEvaluator[] spcEvaluators =
            {
                new BondModularity(),
                new BondCoverage(),
                new BondConductance(),
            };

            // Load the persisted graph artifact when available, else compile and cache it.
            CsrGraph prebuiltGraph = GraphRun.LoadOrBuildGraph(options, dataset, graphConfig, baseDirectory, runDirectory);

            (double[] temperatures, string temperaturesResolved) = ResolveTemperatureGrid(options, prebuiltGraph);

            // Rewrite the manifest with the resolved temperature grid provenance
            manifest = manifest with
            {
                Sweep = manifest.Sweep! with { TemperaturesResolved = temperaturesResolved }
            };
            manifest.WriteTo(runDirectory);

            ISweepStrategy sweepStrategy = BuildSweepStrategy(options, checkpointDirectory, temperatures);
            Console.WriteLine($"Sweep strategy     : {sweepStrategy.GetType().Name} (temperatures: {temperaturesResolved})");

            var result = SpcUserSession.Run(
                dataset,
                graphConfig,
                distanceMetric,
                runPaths,
                externalEvaluators: externalEvaluators,
                spcEvaluators: spcEvaluators,
                referenceLabels: dataset.Labels,
                partitionStrategy: partitionStrategy,
                analyzer: analyzer,
                sweepStrategy: sweepStrategy,
                hierarchicalStrategy: hierarchicalStrategy,
                prebuiltGraph: prebuiltGraph);

            Console.WriteLine("Session finished.");
            Console.WriteLine($"Sweep CSV          : {result.SweepCsvPath}");
            Console.WriteLine($"Partition CSV      : {result.PartitionCsvPath}");
            Console.WriteLine($"Replica traces CSV : {result.ReplicaTracesCsvPath}");
            if (result.EquilibriumEdgesCsvPath is not null)
                Console.WriteLine($"Equilibrium edges  : {result.EquilibriumEdgesCsvPath}");
            if (result.PartitionHierarchyJsonPath is not null)
                Console.WriteLine($"Partition hierarchy: {result.PartitionHierarchyJsonPath}");
            Console.WriteLine($"Criteria CSV       : {result.CriteriaCsvPath}");
            Console.WriteLine($"Session CSV        : {result.SessionCsvPath}");
            Console.WriteLine($"Summary JSON       : {result.SummaryJsonPath}");

            if (options.Resolver != ResolverKind.None)
            {
                var frames = result.SessionResult.SweepRuns.Select(r => r.Accumulator).ToArray();
                CsrGraph resolverGraph = result.SessionResult.Graph;

                (Clustering.Primitives.Assignment assignment, string label, string csvStub) = options.Resolver switch
                {
                    ResolverKind.ThermalEom => (
                        Clustering.Graphical.SPC.Partitions.Thermal.SweepThermalEom.Resolve(
                            resolverGraph, frames, options.Theta,
                            minClusterSize: options.MinClusterSize,
                            completion: options.PeripheryCompletion).Assignment,
                        "thermal-eom", "spc_thermal_eom"),
                    ResolverKind.Hierarchy => ResolveHierarchy(resolverGraph, frames, options),
                    ResolverKind.LineagePersistence  => ResolveLineagePersistence(resolverGraph, frames, result.SessionResult, options),
                    _ => throw new ArgumentOutOfRangeException(nameof(options.Resolver)),
                };

                int[] labels = assignment.Labels;
                int abstained = labels.Count(l => l == Clustering.Primitives.Assignment.Unassigned);
                var sizes = labels
                    .Where(l => l >= 0)
                    .GroupBy(l => l)
                    .Select(g => g.Count())
                    .OrderByDescending(size => size)
                    .ToArray();
                Console.WriteLine(
                    $"Resolver ({label,-11}): {assignment.Count} clusters; " +
                    $"top sizes: {string.Join(",", sizes.Take(8))}; abstained: {abstained}");

                string resolverCsvPath = Path.Combine(
                    Path.GetDirectoryName(result.PartitionCsvPath)!, $"{csvStub}.csv");
                File.WriteAllLines(resolverCsvPath, new[] { "node_index,resolver_label" }
                    .Concat(labels.Select((l, i) => FormattableString.Invariant($"{i},{l}"))));
                Console.WriteLine($"Resolver CSV       : {resolverCsvPath}");
            }

            string tabularPath = result.WriteTabularExports(runPaths.Tabular.Dir);
            Console.WriteLine($"Tabular exports    : {tabularPath}");
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }
    }

    // ── Manifest construction ────────────────────────────────────────────
    private static RunManifest BuildManifest(
        Options options, GraphCompilerConfig graphConfig,
        string runDirectory, string checkpointDirectory, string[] args, RunIdentity identity,
        string? temperaturesResolved)
    {
        DatasetSpec dataset = GraphRun.DatasetSpec(options);
        GraphSpec graph = GraphRun.GraphSpec(options, graphConfig);

        var sweep = new SweepSpec(
            Schedule:          options.Schedule.ToString(),
            TemperaturesSpec:  options.TemperaturesSpec,
            Replicas:          options.Replicas,
            SweepBudget:       options.SweepBudget,
            EquilibriumBudget: options.EquilibriumBudget,
            Q:                 options.Q,
            Analyzer:          options.AnalyzerKind.ToString(),
            PartitionStrategy: options.PartitionStrategyKind.ToString(),
            Theta:             options.Theta,
            TemperaturesResolved: temperaturesResolved);

        return new RunManifest(
            SchemaVersion: RunManifest.CurrentSchemaVersion,
            CreatedUtc:    DateTime.UtcNow,
            Algorithm:     "spc-sw",
            CommandLine:   string.Join(" ", args),
            Dataset:       dataset,
            Graph:         graph,
            Sweep:         sweep,
            Hdbscan:       null,
            CumulativeField: null,
            Output:        new OutputSpec(RunDirectory: runDirectory, CheckpointDirectory: checkpointDirectory),
            Identity:      new RunIdentitySpec(identity.Family, identity.Source, identity.Requested));
    }

    // ── Sweep strategy assembly ──────────────────────────────────────────
    /// <summary>
    /// Resolve the temperature grid: an explicit spec, or the auto bracket the
    /// coupling construction declares (T_ps(q) under 1/K̂, else the data heuristic).
    /// </summary>
    private static (double[] Temperatures, string Resolved) ResolveTemperatureGrid(Options options, CsrGraph graph)
    {
        if (string.IsNullOrWhiteSpace(options.TemperaturesSpec))
            throw new ArgumentException(
                "--schedule fixed-grid requires --temperatures. Examples: " +
                "auto | linspace:0.01,0.5,100 | logspace:0.001,1.0,80 | 0.01,0.05,0.1,0.5");

        string spec = options.TemperaturesSpec.Trim();
        if (spec.Equals("auto", StringComparison.OrdinalIgnoreCase))
        {
            double lo, hi;

            // The Potts sampler's coupling construction declares its T-estimate:
            // under the 1/K̂-normalized replication kernel (MeanEdgeDistance) the
            // q-only physics bracket is valid (T_ps(q) anchor, graphs/models/potts); the
            // data-dependent heuristic is dimensionally J² and lands ~K̂× too cold
            // there — it remains the estimator for the un-normalized default kernel.
            if (options.BandwidthStrategy == Graphs.Distance.BandwidthStrategy.MeanEdgeDistance)
            {
                double tps = BwdPottsCriticalEstimate.TpsUpperBound(options.Q);
                (lo, hi) = (0.05 * tps, 4.0 * tps);
            }
            else
            {
                (lo, hi) = SpcScheduleHelpers.EstimateBracket(graph, options.Q, coldOvershoot: 0.05, hotOvershoot: 5.0);
            }
            return (SpcScheduleHelpers.LogSpaceGrid(lo, hi, 48), $"logspace:{lo:G6},{hi:G6},48");
        }

        return (TemperatureGridSpec.Parse(options.TemperaturesSpec), options.TemperaturesSpec);
    }

    private static ISweepStrategy BuildSweepStrategy(Options options, string checkpointDirectory, double[] temperatures)
    {
        var cfg = new FixedGridSweepConfig
        {
            Temperatures        = temperatures,
            Replicas            = options.Replicas,
            SweepBudget         = options.SweepBudget,
            Sampler             = new PottsModelConfig { Q = options.Q },
            EquilibriumBudget   = options.EquilibriumBudget,
            SusceptibilityKind  = options.SusceptibilityKind,
            Accumulation        = options.Accumulation,
            BaseSeed            = options.Seed,
            CheckpointDirectory = checkpointDirectory,
        };
        return new FixedGridSweepStrategy(cfg);
    }

    // ── Argv parsing ─────────────────────────────────────────────────────
    private static Options ParseArguments(string[] args)
    {
        var options = new Options();

        // Pre-pass: load any --config JSON presets first so CLI flags
        // passed in the same argv override the preset's values. Last
        // --config in argv wins if multiple are passed (each one stacks
        // on top of the previous).
        for (int i = 0; i < args.Length; i++)
        {
            if (args[i].Equals("--config", StringComparison.OrdinalIgnoreCase) && i + 1 < args.Length)
            {
                SpcPreset.Load(args[i + 1]).ApplyTo(options);
            }
        }

        for (int i = 0; i < args.Length; i++)
        {
            if (GraphRun.TryParse(options, args, ref i)) continue;
            string arg = args[i];
            string? next = i + 1 < args.Length ? args[i + 1] : null;
            switch (arg.ToLowerInvariant())
            {
                case "--help": case "-h": options.ShowHelp = true; return options;
                case "--config": i++; break;   // already consumed in pre-pass

                // Output
                case "--no-guid": options.NoGuid = true; break;
                case "--checkpoint-dir": options.CheckpointDir = RequireValue(next, arg); i++; break;

                // Schedule
                case "--schedule": options.Schedule = ParseScheduleMode(RequireValue(next, arg)); i++; break;
                case "--temperatures": options.TemperaturesSpec = RequireValue(next, arg); i++; break;
                case "--replicas": options.Replicas = GraphRun.ParseInt(RequireValue(next, arg)); i++; break;
                case "--sweep": options.SweepBudget = ParseRunBudget(RequireValue(next, arg), arg); i++; break;
                case "--equilibrium": options.EquilibriumBudget = ParseRunBudget(RequireValue(next, arg), arg); i++; break;
                case "--accumulation": options.Accumulation = ParseAccumulationSpec(RequireValue(next, arg)); i++; break;
                case "--q": options.Q = GraphRun.ParseInt(RequireValue(next, arg)); i++; break;

                // Profile analyzer + cut policy
                case "--analyzer": options.AnalyzerKind = ParseAnalyzerKind(RequireValue(next, arg)); i++; break;
                case "--partition-strategy":
                    options.PartitionStrategyKind = ParsePartitionStrategyKind(RequireValue(next, arg));
                    i++; break;
                case "--peripheral-capture": options.PeripheralCapture = true; break;
                case "--hierarchical-strategy":
                    options.HierarchicalStrategyKind = ParseHierarchicalStrategyKind(RequireValue(next, arg));
                    i++; break;
                case "--theta": options.Theta = GraphRun.ParseDouble(RequireValue(next, arg)); i++; break;
                case "--susceptibility":
                    options.SusceptibilityKind = ParseSusceptibilityKind(RequireValue(next, arg)); i++; break;

                case "--resolver":
                    options.Resolver = RequireValue(next, arg).Trim().ToLowerInvariant() switch
                    {
                        "thermal-eom" => ResolverKind.ThermalEom,
                        "hierarchy"   => ResolverKind.Hierarchy,
                        "lineage"   => ResolverKind.LineagePersistence,
                        "none"        => ResolverKind.None,
                        var v => throw new ArgumentException(
                            $"Unknown resolver '{v}'. Valid: none, thermal-eom, hierarchy, lineage."),
                    };
                    i++; break;
                case "--min-cluster-size":
                    options.MinClusterSize = GraphRun.ParseInt(RequireValue(next, arg)); i++; break;
                case "--periphery":
                    options.PeripheryCompletion = RequireValue(next, arg).Trim().ToLowerInvariant() switch
                    {
                        "ascend" => Clustering.Dendrograms.ThermalPeripheryCompletion.Ascend,
                        "none"   => Clustering.Dendrograms.ThermalPeripheryCompletion.None,
                        var v => throw new ArgumentException($"Unknown periphery completion '{v}'. Valid: none, ascend."),
                    };
                    i++; break;

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

    // ── Parsing primitives ───────────────────────────────────────────────
    internal static RunBudget ParseRunBudget(string s, string argName)
    {
        int burnIn = 0, cycles = 0;
        foreach (string part in s.Split(','))
        {
            int eq = part.IndexOf('=');
            if (eq < 0)
                throw new ArgumentException($"'{argName}' expects burnin=N,cycles=N; got '{s}'.");
            string key = part[..eq].Trim();
            string val = part[(eq + 1)..].Trim();
            if (key.Equals("burnin", StringComparison.OrdinalIgnoreCase))
                burnIn = GraphRun.ParseInt(val);
            else if (key.Equals("cycles", StringComparison.OrdinalIgnoreCase))
                cycles = GraphRun.ParseInt(val);
            else
                throw new ArgumentException($"'{argName}': unknown key '{key}'. Expected burnin=N,cycles=N.");
        }
        return new RunBudget(burnIn, cycles);
    }

    internal static AccumulationSpec ParseAccumulationSpec(string s)
    {
        bool affinities = false, alignments = false, coMembership = false;
        bool clusterSizeLandscape = false, orderLandscape = false;
        bool burnInTrace = false, exactClusterSizeHistogram = false;
        foreach (string token in s.Split(','))
        {
            switch (token.Trim().ToLowerInvariant())
            {
                case "affinities":                   affinities                = true; break;
                case "alignments":                   alignments                = true; break;
                case "comembership":                 coMembership              = true; break;
                case "cluster-size-landscape":       clusterSizeLandscape      = true; break;
                case "order-landscape":              orderLandscape            = true; break;
                case "burnin-trace":                 burnInTrace               = true; break;
                case "exact-cluster-size-histogram": exactClusterSizeHistogram = true; break;
                case "none": break;
                default:
                    throw new ArgumentException(
                        $"Unknown accumulation token '{token}'. Valid: none, affinities, alignments, " +
                        "comembership, cluster-size-landscape, order-landscape, burnin-trace, " +
                        "exact-cluster-size-histogram.");
            }
        }
        return new AccumulationSpec
        {
            Affinities = affinities,
            Alignments = alignments,
            CoMembership = coMembership,
            ClusterSizeLandscape = clusterSizeLandscape,
            OrderLandscape = orderLandscape,
            BurnInTrace = burnInTrace,
            ExactClusterSizeHistogram = exactClusterSizeHistogram,
        };
    }

    internal static AnalyzerKind ParseAnalyzerKind(string value) => value.ToLowerInvariant() switch
    {
        "chi-peak"               => AnalyzerKind.ChiPeak,
        "multi-signal-consensus" => throw new ArgumentException(
            "The multi-signal-consensus analyzer is parked (parking-lot/) pending the analysis " +
            "rewrite. Use --analyzer chi-peak."),
        _ => throw new ArgumentException($"Unknown analyzer '{value}'. Valid: chi-peak."),
    };

    internal static SusceptibilityKind ParseSusceptibilityKind(string value) => value.ToLowerInvariant() switch
    {
        "fk-cluster"             => SusceptibilityKind.FkCluster,
        "fk-reduced"             => SusceptibilityKind.FkReduced,
        "magnetization"          => SusceptibilityKind.Magnetization,
        "magnetization-variance" => SusceptibilityKind.MagnetizationVariance,
        _ => throw new ArgumentException(
            $"Unknown susceptibility '{value}'. Valid: fk-cluster, fk-reduced, magnetization, magnetization-variance."),
    };

    internal static PartitionStrategyKind ParsePartitionStrategyKind(string value) => value.ToLowerInvariant() switch
    {
        "co-membership"  => PartitionStrategyKind.CoMembership,
        "spin-agreement" => PartitionStrategyKind.SpinAgreement,
        "bond-frequency" => PartitionStrategyKind.BondFrequency,
        _ => throw new ArgumentException($"Unknown partition strategy '{value}'. Valid: co-membership, spin-agreement, bond-frequency."),
    };

    internal static HierarchicalStrategyKind ParseHierarchicalStrategyKind(string value) => value.ToLowerInvariant() switch
    {
        "none"  or ""      => HierarchicalStrategyKind.None,
        "blatt"            => HierarchicalStrategyKind.Blatt,
        _ => throw new ArgumentException($"Unknown hierarchical strategy '{value}'. Valid: none, blatt."),
    };

    internal static ScheduleMode ParseScheduleMode(string value) => value.ToLowerInvariant() switch
    {
        "fixed-grid" => ScheduleMode.FixedGrid,
        "adaptive"   => throw new ArgumentException(
            "Adaptive scheduling is parked. Use --schedule fixed-grid with --temperatures " +
            "(e.g. logspace:0.001,1.0,80 | linspace:0.01,0.5,100)."),
        _ => throw new ArgumentException($"Unknown schedule '{value}'. Valid: fixed-grid."),
    };

    // ── Analyzer / partition factory ─────────────────────────────────────
    private static ISignalAnalyzer BuildAnalyzer(AnalyzerKind kind) => kind switch
    {
        AnalyzerKind.ChiPeak => new ChiPeakSignalAnalyzer(),
        _ => throw new ArgumentOutOfRangeException(nameof(kind)),
    };

    private static IPartitionStrategy BuildPartitionStrategy(PartitionStrategyKind kind, double theta, bool peripheralCapture) => kind switch
    {
        PartitionStrategyKind.CoMembership  => new ThresholdCoMembership  { Theta = theta, PeripheralCapture = peripheralCapture },
        PartitionStrategyKind.SpinAgreement => new ThresholdSpinAgreement { Theta = theta, PeripheralCapture = peripheralCapture },
        PartitionStrategyKind.BondFrequency => new ThresholdBondFrequency { Theta = theta, PeripheralCapture = peripheralCapture },
        _ => throw new ArgumentOutOfRangeException(nameof(kind)),
    };

    private static IHierarchicalPartitionStrategy? BuildHierarchicalStrategy(
        HierarchicalStrategyKind? kind) => kind switch
    {
        null                                => null,
        HierarchicalStrategyKind.None       => null,
        HierarchicalStrategyKind.Blatt      => new BlattPartitionStrategy(),
        _ => throw new ArgumentOutOfRangeException(nameof(kind)),
    };

    /// <summary>
    /// Track 1 — the Blatt/Domany hierarchical resolver: dense T-stack →
    /// nested-degenerate dendrogram bridge → excess-of-mass walk. Reports the
    /// FKG nesting diagnostic (did the raw sampled stack nest; was it restored).
    /// </summary>
    private static (Clustering.Primitives.Assignment, string, string) ResolveHierarchy(
        CsrGraph graph, Accumulator[] frames, Options options)
    {
        HierarchyEomResult result = SweepHierarchyEom.Resolve(
            graph, frames, options.Theta,
            minClusterSize: options.MinClusterSize,
            completion: options.PeripheryCompletion);
        Console.WriteLine(
            $"Hierarchy nesting  : raw-nested={result.RawNestingHeld}, restored={result.Restored}, " +
            $"levels={result.Stack.Count}");
        return (result.Assignment, "hierarchy", "spc_hierarchy_eom");
    }

    /// <summary>
    /// Track 2 — the lineage-persistence resolver: select cluster lineages by
    /// persistence over the T-stack, bounded by the SP-plateau (T_fs/T_ps).
    /// Reports the data-driven split-share and the persistence span of each
    /// selected lineage.
    /// </summary>
    private static (Clustering.Primitives.Assignment, string, string) ResolveLineagePersistence(
        CsrGraph graph, Accumulator[] frames, Clustering.Graphical.SPC.SpcSessionResult session, Options options)
    {
        // Regime border = the SP-plateau the analyzer already located.
        var plateau = SpcProfileAnalysis.SpPlateau(session.Profile);
        (double, double)? window = plateau.CliffFound ? (plateau.TFs, plateau.TPs) : null;

        LineagePersistenceResult result = SweepLineagePersistence.Resolve(
            graph, frames, options.Theta,
            minClusterSize: Math.Max(options.MinClusterSize, 1),
            temperatureWindow: window);

        Console.WriteLine(
            $"Lineage persistence: {result.AllLineages.Count} tracked, {result.Selected.Count} selected; " +
            $"split-share={result.SplitShare:F2}" +
            (window is { } w ? $"; SP-window=[{w.Item1:G4},{w.Item2:G4}]" : "; SP-window=full grid"));
        return (result.Assignment, "lineage", "spc_lineage_persistence");
    }

    private static string ResolverToken(ResolverKind kind) => kind switch
    {
        ResolverKind.ThermalEom => "thermal-eom",
        ResolverKind.Hierarchy  => "hierarchy",
        ResolverKind.LineagePersistence   => "lineage",
        _                       => "none",
    };

    // ── Enums ────────────────────────────────────────────────────────────
    public enum ScheduleMode { FixedGrid }
    public enum AnalyzerKind { ChiPeak }
    public enum PartitionStrategyKind { CoMembership, SpinAgreement, BondFrequency }
    public enum HierarchicalStrategyKind { None, Blatt }

    /// <summary>The post-sweep T-stack resolver run over the rich sweep's frames.
    /// <c>ThermalEom</c> builds the merge tree from per-edge co-membership
    /// curves; <c>Hierarchy</c> reads the canonical dendrogram-across-T off the
    /// partition stack (Blatt/Domany) — both end in the excess-of-mass walk.
    /// <c>LineagePersistence</c> selects cluster lineages by persistence over the stack
    /// (no dendrogram, no landscape).</summary>
    public enum ResolverKind { None, ThermalEom, Hierarchy, LineagePersistence }

    // ── Options bag ──────────────────────────────────────────────────────
    internal sealed class Options : GraphRunOptions
    {
        // Output (SPC-only: run-folder GUID, checkpoint resume)
        public bool NoGuid { get; set; }
        public string? CheckpointDir { get; set; }

        // Schedule
        public ScheduleMode Schedule { get; set; } = ScheduleMode.FixedGrid;
        public string? TemperaturesSpec { get; set; }
        public int Replicas { get; set; } = 1;
        public RunBudget SweepBudget { get; set; } = new(200, 1000);
        public RunBudget EquilibriumBudget { get; set; } = new(1000, 5000);
        public AccumulationSpec Accumulation { get; set; } = AccumulationSpec.None;
        public int Q { get; set; } = 20;

        // Profile analyzer + cut policy
        public AnalyzerKind AnalyzerKind { get; set; } = AnalyzerKind.ChiPeak;
        public PartitionStrategyKind PartitionStrategyKind { get; set; } = PartitionStrategyKind.CoMembership;
        public HierarchicalStrategyKind HierarchicalStrategyKind { get; set; } = HierarchicalStrategyKind.None;
        public double Theta { get; set; } = 0.5;
        public bool PeripheralCapture { get; set; } = false;
        public SusceptibilityKind SusceptibilityKind { get; set; } = SusceptibilityKind.FkCluster;

        /// <summary>Post-sweep T-stack resolver run over the rich sweep frames
        /// after the flat session. Non-<c>None</c> requires comembership +
        /// cluster-size-landscape accumulation.</summary>
        public ResolverKind Resolver { get; set; } = ResolverKind.None;

        /// <summary>Selection eligibility floor for the thermal resolver.</summary>
        public int MinClusterSize { get; set; } = 1;

        /// <summary>Periphery completion for the thermal resolver's abstains.</summary>
        public Clustering.Dendrograms.ThermalPeripheryCompletion PeripheryCompletion { get; set; }
            = Clustering.Dendrograms.ThermalPeripheryCompletion.None;
    }
}
