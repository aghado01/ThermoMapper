using System;
using System.Globalization;
using Archivory;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;
using Clustering.Evaluation.Internal;
using Clustering.Graphical.CumulativeField;
using Clustering.Graphical.CumulativeField.Export;
using Graphs;
using Graphs.Primitives;

namespace UserRepl.Commands;

/// <summary>
/// The <c>cumulative-field</c> subcommand: compiles the coupling graph through the
/// shared <see cref="GraphRun"/> surface, evaluates the closed-form
/// <see cref="CumulativeFieldSolver"/> over a temperature grid, and resolves one
/// partition by thermal excess of mass into its own run bundle
/// (<see cref="CumulativeFieldRunPaths"/>).
/// </summary>
public static class CumulativeFieldCommand
{
    public static int Run(string[] args)
    {
        try
        {
            var options = ParseArguments(args);
            if (options.ShowHelp) { PrintHelp(); return 0; }
            if (options.ListGenerators) { SpcCommandHelp.PrintGeneratorList(); return 0; }
            if (!string.IsNullOrWhiteSpace(options.GeneratorSchema))
            {
                SpcCommandHelp.PrintGeneratorSchema(options.GeneratorSchema);
                return 0;
            }

            string baseDirectory = options.BaseDirectory ?? "artifacts";
            SpcUserDataset dataset = GraphRun.LoadDataset(options);
            GraphCompilerConfig graphConfig = GraphRun.BuildGraphConfig(options);

            // One scope owns the run tree (artifacts/{family}/{stamp}/).
            var runId = RunIdentity.Resolve(options.RunName, callerStub: "cumulative_field");
            var paths = new CumulativeFieldRunPaths(
                ArtifactScope.Root(baseDirectory, runId.Family, RunStamp.Now()).EnsureDirectory());
            string runDirectory = paths.Scope.Dir;

            // Persist the manifest before any compute so a crash still leaves a
            // self-describing run directory.
            var manifest = BuildManifest(options, graphConfig, runDirectory, args, runId);
            manifest.WriteTo(runDirectory);

            Console.WriteLine($"Run directory      : {runDirectory}");
            GraphRun.PrintGraphSummary(options, graphConfig);
            Console.WriteLine($"Field              : {options.Field} (symmetrization={options.Symmetrization}, theta={options.Theta})");
            Console.WriteLine($"Resolver           : thermal-eom (periphery={options.PeripheryCompletion}, min-size={options.MinClusterSize})");

            CsrGraph graph = GraphRun.LoadOrBuildGraph(options, dataset, graphConfig, baseDirectory, runDirectory);

            (double[] temperatures, string temperaturesResolved) = ResolveTemperatureGrid(options, graph);
            manifest = manifest with
            {
                CumulativeField = manifest.CumulativeField! with { TemperaturesResolved = temperaturesResolved },
            };
            manifest.WriteTo(runDirectory);
            Console.WriteLine($"Temperatures       : {temperaturesResolved}");

            IExternalClusterEvaluator[] externalEvaluators =
            {
                new Purity(),
                new NormalizedMutualInformation(),
                new AdjustedRandIndex(),
                new Homogeneity(),
                new Completeness(),
                new VMeasure(),
            };
            IGraphPartitionEvaluator[] graphEvaluators =
            {
                new BondModularity(),
                new BondCoverage(),
                new BondConductance(),
            };

            CumulativeFieldUserRunResult result = CumulativeFieldUserSession.Run(
                dataset,
                paths,
                graph,
                EdgeWeightKind.Coupling,
                options.Field,
                options.Symmetrization,
                temperatures,
                options.Theta,
                options.MinClusterSize,
                options.PeripheryCompletion,
                options.CutTemperature,
                externalEvaluators,
                graphEvaluators,
                dataset.Labels);

            Console.WriteLine("Session finished.");
            Console.WriteLine($"Thermal EOM        : {result.ClusterCount} clusters; abstained: {result.Abstained}");
            Console.WriteLine($"Sweep CSV          : {result.SweepCsvPath}");
            Console.WriteLine($"Partition CSV      : {result.PartitionCsvPath}");
            if (result.EdgeCurrenciesCsvPath is not null)
                Console.WriteLine($"Edge currencies    : {result.EdgeCurrenciesCsvPath} (T={result.RepresentativeTemperature:G6})");
            Console.WriteLine($"Summary JSON       : {result.SummaryJsonPath}");
            foreach (var kv in result.EvaluatorScores)
                Console.WriteLine($"  {kv.Key,-24}: {kv.Value:F4}");
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }
    }

    // ── Manifest ─────────────────────────────────────────────────────────
    private static RunManifest BuildManifest(
        Options options, GraphCompilerConfig graphConfig, string runDirectory, string[] args, RunIdentity identity)
        => new RunManifest(
            SchemaVersion:   RunManifest.CurrentSchemaVersion,
            CreatedUtc:      DateTime.UtcNow,
            Algorithm:       "cumulative-field",
            CommandLine:     string.Join(" ", args),
            Dataset:         GraphRun.DatasetSpec(options),
            Graph:           GraphRun.GraphSpec(options, graphConfig),
            Sweep:           null,
            Hdbscan:         null,
            CumulativeField: new CumulativeFieldSpec(
                Field:            options.Field.ToString(),
                Symmetrization:   options.Symmetrization.ToString(),
                Theta:            options.Theta,
                TemperaturesSpec: options.TemperaturesSpec,
                MinClusterSize:   options.MinClusterSize,
                Periphery:        options.PeripheryCompletion.ToString(),
                CutTemperature:   options.CutTemperature),
            Output:          new OutputSpec(RunDirectory: runDirectory, CheckpointDirectory: null),
            Identity:        new RunIdentitySpec(identity.Family, identity.Source, identity.Requested));

    // ── Temperature grid ─────────────────────────────────────────────────
    /// <summary>
    /// An explicit spec, or <c>auto</c>: the solver brackets off its own
    /// cumulative-energy ladder (q-free, no Potts anchor) and the grid is 48
    /// log-spaced points across it. The resolved spec is recorded exactly
    /// (round-trip formatting), so re-parsing it reproduces the grid.
    /// </summary>
    private static (double[] Temperatures, string Resolved) ResolveTemperatureGrid(Options options, CsrGraph graph)
    {
        string spec = options.TemperaturesSpec.Trim();
        if (!spec.Equals("auto", StringComparison.OrdinalIgnoreCase))
            return (TemperatureGridSpec.Parse(spec), spec);

        (double lo, double hi) = CumulativeFieldSolver.EstimateBracket(
            graph, EdgeWeightKind.Coupling, options.Field, options.Symmetrization);
        string resolved = string.Create(CultureInfo.InvariantCulture, $"logspace:{lo:R},{hi:R},48");
        return (TemperatureGridSpec.Parse(resolved), resolved);
    }

    // ── Argv ─────────────────────────────────────────────────────────────
    private static Options ParseArguments(string[] args)
    {
        var options = new Options();
        for (int i = 0; i < args.Length; i++)
        {
            if (GraphRun.TryParse(options, args, ref i)) continue;
            string arg = args[i];
            string? next = i + 1 < args.Length ? args[i + 1] : null;
            switch (arg.ToLowerInvariant())
            {
                case "--help": case "-h": options.ShowHelp = true; return options;

                // Cumulative field
                case "--field": options.Field = ParseField(GraphRun.RequireValue(next, arg)); i++; break;
                case "--symmetrization":
                    options.Symmetrization = ParseSymmetrization(GraphRun.RequireValue(next, arg)); i++; break;
                case "--theta": options.Theta = GraphRun.ParseDouble(GraphRun.RequireValue(next, arg)); i++; break;
                case "--temperatures": options.TemperaturesSpec = GraphRun.RequireValue(next, arg); i++; break;

                // Resolution
                case "--min-cluster-size":
                    options.MinClusterSize = GraphRun.ParseInt(GraphRun.RequireValue(next, arg)); i++; break;
                case "--periphery":
                    options.PeripheryCompletion = GraphRun.RequireValue(next, arg).Trim().ToLowerInvariant() switch
                    {
                        "ascend" => ThermalPeripheryCompletion.Ascend,
                        "none"   => ThermalPeripheryCompletion.None,
                        var v => throw new ArgumentException($"Unknown periphery completion '{v}'. Valid: none, ascend."),
                    };
                    i++; break;
                case "--cut-temperature":
                    options.CutTemperature = GraphRun.ParseDouble(GraphRun.RequireValue(next, arg)); i++; break;

                default:
                    throw new ArgumentException($"Unknown argument '{arg}'. Use --help for usage.");
            }
        }
        return options;
    }

    internal static Field ParseField(string value) => value.ToLowerInvariant() switch
    {
        "mean"  => Field.Mean,
        "local" => Field.Local,
        _ => throw new ArgumentException($"Unknown field '{value}'. Valid: mean, local."),
    };

    internal static SymmetrizationRule ParseSymmetrization(string value) => value.ToLowerInvariant() switch
    {
        "mutual"    => SymmetrizationRule.Mutual,
        "inclusive" => SymmetrizationRule.Inclusive,
        "mean"      => SymmetrizationRule.Mean,
        _ => throw new ArgumentException($"Unknown symmetrization '{value}'. Valid: mutual, inclusive, mean."),
    };

    // ── Help ─────────────────────────────────────────────────────────────
    private static void PrintHelp()
    {
        Console.WriteLine("Usage: userrepl cumulative-field --dataset <generator> [options]");
        Console.WriteLine("       userrepl cumulative-field --dataset-file <path> [options]");
        Console.WriteLine();
        Console.WriteLine("Builds the coupling graph, evaluates the closed-form cumulative-field affinity");
        Console.WriteLine("over a temperature grid, and resolves one partition by thermal excess of mass.");
        Console.WriteLine();
        GraphRun.PrintHelp();
        Console.WriteLine("Cumulative field:");
        Console.WriteLine("  --field <mean|local>                     Ladder pooling (default: mean)");
        Console.WriteLine("  --symmetrization <mutual|inclusive|mean> Reconcile a local edge's two directed values");
        Console.WriteLine("                                           (default: mutual)");
        Console.WriteLine("  --theta <double>                         Affinity threshold (default: 0.5)");
        Console.WriteLine("  --temperatures <spec>                    auto (default: bracket from the ladder)");
        Console.WriteLine("                                           | linspace:Tmin,Tmax,N | logspace:Tmin,Tmax,N");
        Console.WriteLine("                                           | 0.01,0.05,0.1  (explicit list)");
        Console.WriteLine();
        Console.WriteLine("Resolution:");
        Console.WriteLine("  --min-cluster-size <int>                 Selection floor (default: 1)");
        Console.WriteLine("  --periphery <none|ascend>                Complete abstains by modal ascent (default: none)");
        Console.WriteLine("  --cut-temperature <double>               Temperature of the exported edge-currency slice");
        Console.WriteLine("                                           (default: the coldest grid point)");
        Console.WriteLine();
        Console.WriteLine("  --help, -h                    Show this help text");
    }

    // ── Options bag ──────────────────────────────────────────────────────
    internal sealed class Options : GraphRunOptions
    {
        public Field Field { get; set; } = Field.Mean;
        public SymmetrizationRule Symmetrization { get; set; } = SymmetrizationRule.Mutual;
        public double Theta { get; set; } = 0.5;
        public string TemperaturesSpec { get; set; } = "auto";
        public int MinClusterSize { get; set; } = 1;
        public ThermalPeripheryCompletion PeripheryCompletion { get; set; } = ThermalPeripheryCompletion.None;
        public double? CutTemperature { get; set; }
    }
}
