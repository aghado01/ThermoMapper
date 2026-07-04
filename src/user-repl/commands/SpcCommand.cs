using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
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

using Clustering.Graphical.SPC.Runtime.Core.Solver;
using SolverField = Clustering.Graphical.SPC.Runtime.Core.Solver.Field;

namespace UserRepl.Commands;

/// <summary>
/// The <c>spc</c> subcommand. Runs an end-to-end SPC clustering session
/// against a synthetic-generator dataset or a CSV input, with adaptive
/// or fixed-grid temperature sweeps and optional checkpoint resume.
/// </summary>
public static class SpcCommand
{
    private static readonly IReadOnlyList<string> DatasetKinds =
        SpcUserSession.ListAvailableSyntheticGeneratorNames();

    public static int Run(string[] args)
    {
        try
        {
            var options = ParseArguments(args);
            if (options.ShowHelp) { SpcCommandHelp.PrintHelp(DatasetKinds); return 0; }
            if (options.ListGenerators) { SpcCommandHelp.PrintGeneratorList(); return 0; }
            if (!string.IsNullOrWhiteSpace(options.GeneratorSchema))
            {
                SpcCommandHelp.PrintGeneratorSchema(options.GeneratorSchema);
                return 0;
            }

            string baseDirectory = options.BaseDirectory ?? "artifacts";

            SpcUserDataset dataset = LoadDataset(options);

            var graphConfig = BuildGraphConfig(options);
            IDistanceMetric? distanceMetric = string.IsNullOrWhiteSpace(options.DistanceMetricSpec)
                ? null
                : DistanceMetricFactory.Create(options.DistanceMetricSpec);

            string datasetFingerprint = GraphConstructionPersistence.ComputeDatasetFingerprint(
                dataset, distanceMetric);
            string graphConfigFingerprint = GraphConstructionPersistence.ComputeConfigFingerprint(graphConfig);
            // Graph cache is orthogonal to the run tree: a content-addressed sibling
            // under artifacts/graph_cache/, never a parent of the run.
            string datasetRoot = GraphConstructionPersistence.ChooseDatasetRoot(
                Path.Combine(baseDirectory, "graph_cache"), dataset, datasetFingerprint, graphConfigFingerprint);

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
            // Each kernel descriptor is a record that already renders itself
            // (e.g. "Gaussian { Bandwidth = 0.5 }"); no kernel-variant switch here.
            // Distance mode has no kernel — fall back to the projection's name.
            IKernelDescriptor? displayKernel = (graphConfig.Projection as CouplingProjection)?.Kernel;
            string kernelDesc = displayKernel?.ToString() ?? graphConfig.Projection.GetType().Name;

            string topologyDesc = graphConfig.Topology.Kind == TopologyKind.EpsilonBall
                ? $"epsilon (ε={graphConfig.Topology.Epsilon:G6})"
                : $"knn (k={graphConfig.Topology.K ?? 10})";

            Console.WriteLine(
                $"Graph              : {topologyDesc}, " +
                $"Kernel={kernelDesc}, " +
                $"DistanceMetric={(options.DistanceMetricSpec ?? "EuclideanDefault")}, " +
                $"MST-repair={(graphConfig.Repair.Kind == RepairKind.MstMin)}, " +
                $"LMP={((graphConfig.Projection as CouplingProjection)?.LmpRescale == true)}");
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

            // Load a persisted graph artifact when available to avoid
            // rebuilding the same graph repeatedly for SPC sweep runs.
            CsrGraph prebuiltGraph;
            GraphConstructionManifest? graphManifest = null;
            if (GraphConstructionPersistence.TryLoadGraphArtifact(datasetRoot, out var persistedGraph, out var persistedManifest))
            {
                prebuiltGraph = persistedGraph;
                graphManifest = persistedManifest;
                Console.WriteLine($"Loaded persisted graph artifact from {datasetRoot}");
                Console.WriteLine($"Graph artifact : {GraphConstructionPersistence.GetGraphPath(datasetRoot)}");
                Console.WriteLine($"Graph manifest : {GraphConstructionPersistence.GetManifestPath(datasetRoot)}");
                Console.WriteLine($"Graph manifest contains {graphManifest?.Diagnostics.Messages.Count ?? 0} diagnostic entries.");
            }
            else
            {
                GraphBuildResult buildResult = SpcGraphBuilder.BuildResult(
                    dataset.Features, graphConfig, distanceMetric,
                    protectedEdges: TDA.Ph.H1CycleEdges.FromDistanceGraph);

                graphManifest = GraphConstructionPersistence.MaterializeManifest(
                    graphConfig, buildResult, datasetFingerprint);
                GraphConstructionPersistence.WriteGraphArtifact(datasetRoot, buildResult, graphManifest);
                prebuiltGraph = buildResult.Graph;

                Console.WriteLine($"Persisted graph artifact : {GraphConstructionPersistence.GetGraphPath(datasetRoot)}");
                Console.WriteLine($"Persisted graph manifest : {GraphConstructionPersistence.GetManifestPath(datasetRoot)}");

                GraphHealthReport health = GraphHealth.Evaluate(buildResult, k: graphConfig.Topology.K ?? 0);
                string graphHealthPath = GraphHealthFile.WriteTo(runDirectory, health);
                if (!string.IsNullOrEmpty(health.Verdict.PrimaryRecommendation))
                {
                    Console.WriteLine();
                    Console.WriteLine($"Graph health       : {health.Verdict.PrimaryRecommendation}");
                    Console.WriteLine();
                }

                Console.WriteLine($"Graph health JSON  : {graphHealthPath}");
            }

            // Shared schedule envelope: both siblings ride the identical grid
            // (same bracket, same q-anchor). The fork is below — SW runs the
            // MC sweep strategy; PKWang runs the closed-form solver session.
            (double[] temperatures, string temperaturesResolved) = ResolveTemperatureGrid(options, prebuiltGraph);

            // Rewrite the manifest with the resolved temperature grid provenance
            manifest = manifest with
            {
                Sweep = manifest.Sweep! with { TemperaturesResolved = temperaturesResolved }
            };
            manifest.WriteTo(runDirectory);

            if (options.Solver == SolverKind.PKWang)
            {
                SolverField field = options.Field ?? SolverField.Mean;
                Console.WriteLine($"Solver             : PKWang (field={field}, symmetrization={options.Symmetrization})");
                Console.WriteLine($"Resolver           : thermal-eom (periphery={options.PeripheryCompletion}, min-size={options.MinClusterSize})");
                Console.WriteLine($"Temperatures       : {temperaturesResolved}");

                PKWangUserRunResult pk = PKWangUserSession.Run(
                    dataset,
                    runPaths,
                    prebuiltGraph,
                    EdgeWeightKind.Coupling,
                    field,
                    options.Symmetrization,
                    temperatures,
                    options.Theta,
                    options.MinClusterSize,
                    options.PeripheryCompletion,
                    options.CutTemperature,
                    externalEvaluators,
                    spcEvaluators,
                    dataset.Labels);

                Console.WriteLine("Session finished.");
                Console.WriteLine($"Thermal EOM        : {pk.ClusterCount} clusters; abstained: {pk.Abstained}");
                Console.WriteLine($"Sweep CSV          : {pk.SweepCsvPath}");
                Console.WriteLine($"Partition CSV      : {pk.PartitionCsvPath}");
                if (pk.EquilibriumEdgesCsvPath is not null)
                    Console.WriteLine($"Equilibrium edges  : {pk.EquilibriumEdgesCsvPath} (T={pk.RepresentativeTemperature:G6})");
                Console.WriteLine($"Summary JSON       : {pk.SummaryJsonPath}");
                foreach (var kv in pk.EvaluatorScores)
                    Console.WriteLine($"  {kv.Key,-24}: {kv.Value:F4}");
                return 0;
            }

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
                        Clustering.Graphical.SPC.Partitions.Thermal.ThermalEom.Resolve(
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

    // ── Dataset loading ──────────────────────────────────────────────────
    private static SpcUserDataset LoadDataset(Options options)
    {
        if (!string.IsNullOrWhiteSpace(options.DatasetFile))
        {
            Console.WriteLine($"Loading dataset    : CSV {options.DatasetFile}");
            return SpcUserSession.FromCsv(
                options.DatasetFile, options.LabelColumn, options.HasHeader, options.Delimiter);
        }

        Console.WriteLine($"Loading dataset    : synthetic '{options.DatasetKind}'");
        var generatorParameters = new Dictionary<string, object?>(StringComparer.OrdinalIgnoreCase)
        {
            ["seed"] = options.Seed,
        };
        foreach (var kvp in options.GeneratorParameters)
            generatorParameters[kvp.Key] = kvp.Value;
        return SpcUserSession.GenerateDataset(options.DatasetKind, generatorParameters);
    }

    // ── Graph config assembly ────────────────────────────────────────────
    private static GraphCompilerConfig BuildGraphConfig(Options options)
    {
        IDistanceMetric? distanceMetric = string.IsNullOrWhiteSpace(options.DistanceMetricSpec)
            ? null
            : DistanceMetricFactory.Create(options.DistanceMetricSpec);

        IKernelDescriptor kernel = string.IsNullOrWhiteSpace(options.MixtureSpec)
            ? options.KernelType switch
            {
                KernelType.Gaussian  => new Gaussian(options.Bandwidth),
                KernelType.Cauchy    => new Cauchy(options.Bandwidth),
                KernelType.Laplacian => new Laplacian(options.Bandwidth),
                KernelType.Linear    => new Linear(options.Bandwidth),
                _ => throw new NotSupportedException($"Unsupported kernel type {options.KernelType}.")
            }
            : new Mixture(
                GaussianWeight:  MixtureSpecParser.ParseWeights(options.MixtureSpec).Gaussian,
                CauchyWeight:    MixtureSpecParser.ParseWeights(options.MixtureSpec).Cauchy,
                LaplacianWeight: MixtureSpecParser.ParseWeights(options.MixtureSpec).Laplacian,
                GaussianBandwidth: options.MixtureBandwidthSpec is null ? 0.0 : MixtureSpecParser.ParseBandwidth(options.MixtureBandwidthSpec).Gaussian,
                CauchyBandwidth:   options.MixtureBandwidthSpec is null ? 0.0 : MixtureSpecParser.ParseBandwidth(options.MixtureBandwidthSpec).Cauchy,
                LaplacianBandwidth: options.MixtureBandwidthSpec is null ? 0.0 : MixtureSpecParser.ParseBandwidth(options.MixtureBandwidthSpec).Laplacian);

        return new GraphCompilerConfig
        {
            Topology = options.TopologyKind == TopologyKind.EpsilonBall
                ? new TopologyConfig { Kind = TopologyKind.EpsilonBall, Epsilon = options.Epsilon }
                : new TopologyConfig { Kind = TopologyKind.Knn, K = options.K },
            Filter = new FilterConfig
            {
                Kind = options.FilterKind,
            },
            Repair = new RepairConfig
            {
                Kind = options.EnsureConnected ? RepairKind.MstMin : RepairKind.NoRepair,
            },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection
            {
                Kernel = kernel,
                LmpRescale = options.ApplyLmp,
                BandwidthOverride = options.BandwidthStrategy,
            },
            Interrupts = new PathologyInterruptConfig(),
        };
    }

    // ── Manifest construction ────────────────────────────────────────────
    private static RunManifest BuildManifest(
        Options options, GraphCompilerConfig graphConfig,
        string runDirectory, string checkpointDirectory, string[] args, RunIdentity identity,
        string? temperaturesResolved)
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

            // The projection (kernel + LMP + bandwidth, or distance pass-through) is
            // embedded in the manifest as-is and round-trips polymorphically — no kernel
            // flattening switch to keep in sync with ManifestMaterialization.ToConfig.
            var graph = new GraphSpec(
                TopologyKind:    graphConfig.Topology.Kind == TopologyKind.EpsilonBall ? "epsilon" : "knn",
                FilterKind:      graphConfig.Filter.Kind == FilterKind.MutualKnn ? "mutualknn" : "or_rule",
                K:               graphConfig.Topology.K ?? 0,
                Epsilon:         graphConfig.Topology.Epsilon ?? 0.0,
                DistanceMetric:  options.DistanceMetricSpec,
                EnsureConnected: graphConfig.Repair.Kind == RepairKind.MstMin,
                Projection:      graphConfig.Projection);

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
            Algorithm:     options.Solver == SolverKind.PKWang ? "spc-pkwang" : "spc-sw",
            CommandLine:   string.Join(" ", args),
            Dataset:       dataset,
            Graph:         graph,
            Sweep:         sweep,
            Hdbscan:       null,
            Output:        new OutputSpec(RunDirectory: runDirectory, CheckpointDirectory: checkpointDirectory),
            Identity:      new RunIdentitySpec(identity.Family, identity.Source, identity.Requested));
    }

    // ── Sweep strategy assembly ──────────────────────────────────────────
    /// <summary>
    /// Resolve the temperature grid — the LAST step shared by both siblings.
    /// SW and PKWang ride the identical bracket so PKWang is a true
    /// apples-to-apples (zero-variance) probe of the SW schedule.
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

            // The analytical solver brackets off its OWN cumulative-energy ladder:
            // T is a single-linkage cut height (Hcum = T·ln2), so the only
            // temperatures that change anything span [min,max] Hcum / ln2. q-free
            // and gauge-free — no borrowed Potts T_ps(q) (the solver has no
            // susceptibility transition to anchor on).
            if (options.Solver == SolverKind.PKWang)
            {
                (lo, hi) = PKWang.EstimateBracket(
                    graph, EdgeWeightKind.Coupling, options.Field ?? SolverField.Mean, options.Symmetrization);
                return (SpcScheduleHelpers.LogSpaceGrid(lo, hi, 48), $"logspace:{lo:G6},{hi:G6},48 (pkwang-hcum)");
            }

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
            string arg = args[i];
            string? next = i + 1 < args.Length ? args[i + 1] : null;
            switch (arg.ToLowerInvariant())
            {
                case "--help": case "-h": options.ShowHelp = true; return options;
                case "--config": i++; break;   // already consumed in pre-pass
                case "--list-generators": options.ListGenerators = true; break;
                case "--generator-schema": options.GeneratorSchema = RequireValue(next, arg); i++; break;

                // Dataset
                case "--dataset": options.DatasetKind = ParseDatasetKind(RequireValue(next, arg)); i++; break;
                case "--param": ParseGeneratorParameter(RequireValue(next, arg), options.GeneratorParameters); i++; break;
                case "--dataset-file": options.DatasetFile = RequireValue(next, arg); i++; break;
                case "--label-column": options.LabelColumn = RequireValue(next, arg); i++; break;
                case "--delimiter": options.Delimiter = ParseDelimiter(RequireValue(next, arg)); i++; break;
                case "--no-header": options.HasHeader = false; break;

                // Output
                case "--base-dir": options.BaseDirectory = RequireValue(next, arg); i++; break;
                case "--run-name": options.RunName = RequireValue(next, arg); i++; break;
                case "--no-guid": options.NoGuid = true; break;
                case "--checkpoint-dir": options.CheckpointDir = RequireValue(next, arg); i++; break;

                // Reproducibility
                case "--seed": options.Seed = ParseInt(RequireValue(next, arg)); i++; break;

                // Graph topology
                case "--k": options.K = ParseInt(RequireValue(next, arg)); i++; break;
                case "--epsilon": options.Epsilon = ParseDouble(RequireValue(next, arg)); i++; break;
                case "--topology": options.TopologyKind = ParseTopologyKind(RequireValue(next, arg)); i++; break;
                case "--filter": options.FilterKind = ParseFilterKind(RequireValue(next, arg)); i++; break;
                case "--distance-metric": options.DistanceMetricSpec = RequireValue(next, arg); i++; break;
                case "--ensure-connected": options.EnsureConnected = true; break;
                case "--lmp": options.ApplyLmp = true; break;

                // Coupling kernel
                case "--kernel": options.KernelType = ParseKernelType(RequireValue(next, arg)); i++; break;
                case "--bandwidth": options.Bandwidth = ParseDouble(RequireValue(next, arg)); i++; break;
                case "--mixture": options.MixtureSpec = RequireValue(next, arg); i++; break;
                case "--mixture-bandwidth": options.MixtureBandwidthSpec = RequireValue(next, arg); i++; break;

                // Schedule
                case "--schedule": options.Schedule = ParseScheduleMode(RequireValue(next, arg)); i++; break;
                case "--temperatures": options.TemperaturesSpec = RequireValue(next, arg); i++; break;
                case "--replicas": options.Replicas = ParseInt(RequireValue(next, arg)); i++; break;
                case "--sweep": options.SweepBudget = ParseRunBudget(RequireValue(next, arg), arg); i++; break;
                case "--equilibrium": options.EquilibriumBudget = ParseRunBudget(RequireValue(next, arg), arg); i++; break;
                case "--accumulation": options.Accumulation = ParseAccumulationSpec(RequireValue(next, arg)); i++; break;
                case "--q": options.Q = ParseInt(RequireValue(next, arg)); i++; break;

                // Profile analyzer + cut policy
                case "--analyzer": options.AnalyzerKind = ParseAnalyzerKind(RequireValue(next, arg)); i++; break;
                case "--partition-strategy":
                    options.PartitionStrategyKind = ParsePartitionStrategyKind(RequireValue(next, arg));
                    i++; break;
                case "--peripheral-capture": options.PeripheralCapture = true; break;
                case "--hierarchical-strategy":
                    options.HierarchicalStrategyKind = ParseHierarchicalStrategyKind(RequireValue(next, arg));
                    i++; break;
                case "--theta": options.Theta = ParseDouble(RequireValue(next, arg)); i++; break;
                case "--susceptibility":
                    options.SusceptibilityKind = ParseSusceptibilityKind(RequireValue(next, arg)); i++; break;
                case "--bandwidth-strategy":
                    options.BandwidthStrategy = ParseBandwidthStrategy(RequireValue(next, arg)); i++; break;

                // Inference method (sampler vs analytic solver)
                case "--solver": options.Solver = ParseSolverKind(RequireValue(next, arg)); i++; break;
                case "--field": options.Field = ParseSolverField(RequireValue(next, arg)); i++; break;
                case "--symmetrization":
                    options.Symmetrization = ParseSymmetrization(RequireValue(next, arg)); i++; break;
                case "--cut-temperature":
                    options.CutTemperature = ParseDouble(RequireValue(next, arg)); i++; break;
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
                    options.MinClusterSize = ParseInt(RequireValue(next, arg)); i++; break;
                case "--periphery":
                    options.PeripheryCompletion = RequireValue(next, arg).Trim().ToLowerInvariant() switch
                    {
                        "ascend" => Clustering.Graphical.SPC.Partitions.Thermal.ThermalPeripheryCompletion.Ascend,
                        "none"   => Clustering.Graphical.SPC.Partitions.Thermal.ThermalPeripheryCompletion.None,
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
    private static int ParseInt(string s)
        => int.Parse(s, CultureInfo.InvariantCulture);

    private static double ParseDouble(string s)
        => double.Parse(s, CultureInfo.InvariantCulture);

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
                burnIn = ParseInt(val);
            else if (key.Equals("cycles", StringComparison.OrdinalIgnoreCase))
                cycles = ParseInt(val);
            else
                throw new ArgumentException($"'{argName}': unknown key '{key}'. Expected burnin=N,cycles=N.");
        }
        return new RunBudget(burnIn, cycles);
    }

    internal static AccumulationSpec ParseAccumulationSpec(string s)
    {
        bool affinities = false, alignments = false, coMembership = false;
        bool clusterSizeLandscape = false, orderLandscape = false;
        foreach (string token in s.Split(','))
        {
            switch (token.Trim().ToLowerInvariant())
            {
                case "affinities":             affinities           = true; break;
                case "alignments":             alignments           = true; break;
                case "comembership":           coMembership         = true; break;
                case "cluster-size-landscape": clusterSizeLandscape = true; break;
                case "order-landscape":        orderLandscape       = true; break;
                case "none": break;
                default:
                    throw new ArgumentException(
                        $"Unknown accumulation token '{token}'. Valid: none, affinities, alignments, " +
                        "comembership, cluster-size-landscape, order-landscape.");
            }
        }
        return new AccumulationSpec
        {
            Affinities = affinities,
            Alignments = alignments,
            CoMembership = coMembership,
            ClusterSizeLandscape = clusterSizeLandscape,
            OrderLandscape = orderLandscape,
        };
    }

    internal static string ParseDatasetKind(string value)
    {
        var match = DatasetKinds.FirstOrDefault(name => string.Equals(name, value, StringComparison.OrdinalIgnoreCase));
        if (match is not null) return match;
        throw new ArgumentException($"Unknown dataset kind '{value}'. Valid: {string.Join(", ", DatasetKinds)}.");
    }

    internal static KernelType ParseKernelType(string value) => value.ToLowerInvariant() switch
    {
        "gaussian"  => KernelType.Gaussian,
        "cauchy"    => KernelType.Cauchy,
        "laplacian" => KernelType.Laplacian,
        "linear"    => KernelType.Linear,
        _ => throw new ArgumentException($"Unknown kernel '{value}'. Valid: gaussian, cauchy, laplacian, linear."),
    };

    internal static TopologyKind ParseTopologyKind(string value) => value.ToLowerInvariant() switch
    {
        "knn"     => TopologyKind.Knn,
        "epsilon" => TopologyKind.EpsilonBall,
        _ => throw new ArgumentException($"Unknown topology kind '{value}'. Valid: knn, epsilon."),
    };

    internal static FilterKind ParseFilterKind(string value) => value.ToLowerInvariant() switch
    {
        "or_rule"   => FilterKind.OrRule,
        "or"        => FilterKind.OrRule,
        "mutualknn" => FilterKind.MutualKnn,
        _ => throw new ArgumentException($"Unknown filter kind '{value}'. Valid: or_rule, mutualknn."),
    };

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

    internal static Graphs.Distance.BandwidthStrategy ParseBandwidthStrategy(string value) => value.ToLowerInvariant() switch
    {
        "mad"                  => Graphs.Distance.BandwidthStrategy.MadConsistencyFactor,
        "quantile-normalized"  => Graphs.Distance.BandwidthStrategy.QuantileNormalized,
        "log-scale-hyperbolic" => Graphs.Distance.BandwidthStrategy.LogScaleHyperbolic,
        "mean-edge-distance"   => Graphs.Distance.BandwidthStrategy.MeanEdgeDistance,
        _ => throw new ArgumentException(
            $"Unknown bandwidth strategy '{value}'. Valid: mad, quantile-normalized, log-scale-hyperbolic, mean-edge-distance."),
    };

    internal static SolverKind ParseSolverKind(string value) => value.ToLowerInvariant() switch
    {
        "sw" or "swendsen-wang" or "swendsenwang" => SolverKind.SwendsenWang,
        "pkwang" or "pk-wang" or "wang"           => SolverKind.PKWang,
        _ => throw new ArgumentException($"Unknown solver '{value}'. Valid: sw, pkwang."),
    };

    internal static SolverField ParseSolverField(string value) => value.ToLowerInvariant() switch
    {
        "mean"  => SolverField.Mean,
        "local" => SolverField.Local,
        _ => throw new ArgumentException($"Unknown field '{value}'. Valid: mean, local."),
    };

    internal static SymmetrizationRule ParseSymmetrization(string value) => value.ToLowerInvariant() switch
    {
        "mutual"    => SymmetrizationRule.Mutual,
        "inclusive" => SymmetrizationRule.Inclusive,
        "mean"      => SymmetrizationRule.Mean,
        _ => throw new ArgumentException($"Unknown symmetrization '{value}'. Valid: mutual, inclusive, mean."),
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
        string key = token[..eq].Trim();
        string val = token[(eq + 1)..].Trim();
        if (string.IsNullOrWhiteSpace(key))
            throw new ArgumentException($"Generator parameter name cannot be empty: '{token}'");
        parameters[key] = val;
    }

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
        HierarchyEomResult result = HierarchyEom.Resolve(
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

        LineagePersistenceResult result = LineagePersistence.Resolve(
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
    public enum SolverKind { SwendsenWang, PKWang }
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
    internal sealed class Options
    {
        public bool ShowHelp { get; set; }
        public bool ListGenerators { get; set; }
        public string? GeneratorSchema { get; set; }

        // Dataset
        public string DatasetKind { get; set; } = DatasetKinds.FirstOrDefault() ?? "unknown";
        public string? DatasetFile { get; set; }
        public string? LabelColumn { get; set; }
        public bool HasHeader { get; set; } = true;
        public char Delimiter { get; set; } = ',';
        public Dictionary<string, object?> GeneratorParameters { get; } = new(StringComparer.OrdinalIgnoreCase);

        // Output / reproducibility
        public string? BaseDirectory { get; set; }
        public string? RunName { get; set; }
        public bool NoGuid { get; set; }
        public string? CheckpointDir { get; set; }
        public int Seed { get; set; } = 42;

        // Graph topology
        public int K { get; set; } = 10;
        public double Epsilon { get; set; } = 0.0;
        public TopologyKind TopologyKind { get; set; } = TopologyKind.Knn;
        public FilterKind FilterKind { get; set; } = FilterKind.OrRule;
        public string? DistanceMetricSpec { get; set; }
        public bool EnsureConnected { get; set; }
        public bool ApplyLmp { get; set; }

        // Coupling kernel
        public KernelType KernelType { get; set; } = KernelType.Gaussian;
        public double Bandwidth { get; set; } = 0.0;
        public string? MixtureSpec { get; set; }
        public string? MixtureBandwidthSpec { get; set; }

        // Inference method (sampler vs analytic solver)
        public SolverKind Solver { get; set; } = SolverKind.SwendsenWang;

        /// <summary>PKWang energy-ladder field (null defaults to Mean). Inert for SW.</summary>
        public SolverField? Field { get; set; }

        /// <summary>Directed-field reconciliation for PKWang LocalField. Inert for SW / MeanField.</summary>
        public SymmetrizationRule Symmetrization { get; set; } = SymmetrizationRule.Mutual;

        /// <summary>Explicit single-cut temperature for the PKWang solver; null = longest
        /// cluster-count plateau (the solver's native stability signal).</summary>
        public double? CutTemperature { get; set; }

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

        /// <summary>Null = resolve from metric properties / the MAD default;
        /// <c>MeanEdgeDistance</c> selects the BWD replication kernel
        /// (mean-pair <c>a</c> + 1/K̂) and the q-only auto bracket.</summary>
        public Graphs.Distance.BandwidthStrategy? BandwidthStrategy { get; set; }

        /// <summary>Post-sweep T-stack resolver run over the rich sweep frames
        /// after the flat session. Non-<c>None</c> requires comembership +
        /// cluster-size-landscape accumulation.</summary>
        public ResolverKind Resolver { get; set; } = ResolverKind.None;

        /// <summary>Selection eligibility floor for the thermal resolver.</summary>
        public int MinClusterSize { get; set; } = 1;

        /// <summary>Periphery completion for the thermal resolver's abstains.</summary>
        public Clustering.Graphical.SPC.Partitions.Thermal.ThermalPeripheryCompletion PeripheryCompletion { get; set; }
            = Clustering.Graphical.SPC.Partitions.Thermal.ThermalPeripheryCompletion.None;
    }
}
