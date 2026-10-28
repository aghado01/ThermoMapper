using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Graphs;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Distance;
using Graphs.Primitives;

namespace UserRepl.Commands;

/// <summary>
/// The shared front half of every graph-based clustering command: parse the
/// dataset + coupling-graph flags of <see cref="GraphRunOptions"/>, load the
/// dataset, assemble the <see cref="GraphCompilerConfig"/>, load or compile the
/// graph through the content-addressed cache (writing the graph-health report on
/// a fresh build), and describe both in the run manifest. Each command adds its
/// method on top of the compiled <see cref="CsrGraph"/>.
/// </summary>
internal static class GraphRun
{
    internal static readonly IReadOnlyList<string> DatasetKinds =
        SpcUserSession.ListAvailableSyntheticGeneratorNames();

    // ── Argv ─────────────────────────────────────────────────────────────
    /// <summary>
    /// Consume the shared flag at <c>args[i]</c> (advancing <paramref name="i"/>
    /// past its value). Returns <see langword="false"/> when the flag is not part
    /// of the shared surface, leaving it to the command's own switch.
    /// </summary>
    internal static bool TryParse(GraphRunOptions options, string[] args, ref int i)
    {
        string arg = args[i];
        string? next = i + 1 < args.Length ? args[i + 1] : null;
        switch (arg.ToLowerInvariant())
        {
            case "--list-generators": options.ListGenerators = true; return true;
            case "--generator-schema": options.GeneratorSchema = RequireValue(next, arg); i++; return true;

            // Dataset
            case "--dataset": options.DatasetKind = ParseDatasetKind(RequireValue(next, arg)); i++; return true;
            case "--param": ParseGeneratorParameter(RequireValue(next, arg), options.GeneratorParameters); i++; return true;
            case "--dataset-file": options.DatasetFile = RequireValue(next, arg); i++; return true;
            case "--label-column": options.LabelColumn = RequireValue(next, arg); i++; return true;
            case "--delimiter": options.Delimiter = ParseDelimiter(RequireValue(next, arg)); i++; return true;
            case "--no-header": options.HasHeader = false; return true;

            // Output / reproducibility
            case "--base-dir": options.BaseDirectory = RequireValue(next, arg); i++; return true;
            case "--run-name": options.RunName = RequireValue(next, arg); i++; return true;
            case "--seed": options.Seed = ParseInt(RequireValue(next, arg)); i++; return true;

            // Graph topology
            case "--k": options.K = ParseInt(RequireValue(next, arg)); i++; return true;
            case "--epsilon": options.Epsilon = ParseDouble(RequireValue(next, arg)); i++; return true;
            case "--topology": options.TopologyKind = ParseTopologyKind(RequireValue(next, arg)); i++; return true;
            case "--filter": options.FilterKind = ParseFilterKind(RequireValue(next, arg)); i++; return true;
            case "--distance-metric": options.DistanceMetricSpec = RequireValue(next, arg); i++; return true;
            case "--ensure-connected": options.EnsureConnected = true; return true;
            case "--lmp": options.ApplyLmp = true; return true;

            // Coupling kernel
            case "--kernel": options.KernelType = ParseKernelType(RequireValue(next, arg)); i++; return true;
            case "--bandwidth": options.Bandwidth = ParseDouble(RequireValue(next, arg)); i++; return true;
            case "--mixture": options.MixtureSpec = RequireValue(next, arg); i++; return true;
            case "--mixture-bandwidth": options.MixtureBandwidthSpec = RequireValue(next, arg); i++; return true;
            case "--bandwidth-strategy": options.BandwidthStrategy = ParseBandwidthStrategy(RequireValue(next, arg)); i++; return true;

            default: return false;
        }
    }

    internal static string RequireValue(string? value, string option)
    {
        if (string.IsNullOrWhiteSpace(value))
            throw new ArgumentException($"Option '{option}' requires a value.");
        return value;
    }

    // ── Dataset ──────────────────────────────────────────────────────────
    internal static SpcUserDataset LoadDataset(GraphRunOptions options)
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

    // ── Graph ────────────────────────────────────────────────────────────
    internal static IDistanceMetric? DistanceMetric(GraphRunOptions options)
        => string.IsNullOrWhiteSpace(options.DistanceMetricSpec)
            ? null
            : DistanceMetricFactory.Create(options.DistanceMetricSpec);

    internal static GraphCompilerConfig BuildGraphConfig(GraphRunOptions options)
    {
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

    /// <summary>The one-line graph description every command prints.</summary>
    internal static void PrintGraphSummary(GraphRunOptions options, GraphCompilerConfig config)
    {
        // Each kernel descriptor is a record that already renders itself
        // (e.g. "Gaussian { Bandwidth = 0.5 }"); no kernel-variant switch here.
        // Distance mode has no kernel — fall back to the projection's name.
        IKernelDescriptor? displayKernel = (config.Projection as CouplingProjection)?.Kernel;
        string kernelDesc = displayKernel?.ToString() ?? config.Projection.GetType().Name;

        string topologyDesc = config.Topology.Kind == TopologyKind.EpsilonBall
            ? $"epsilon (ε={config.Topology.Epsilon:G6})"
            : $"knn (k={config.Topology.K ?? 10})";

        Console.WriteLine(
            $"Graph              : {topologyDesc}, " +
            $"Kernel={kernelDesc}, " +
            $"DistanceMetric={(options.DistanceMetricSpec ?? "EuclideanDefault")}, " +
            $"MST-repair={(config.Repair.Kind == RepairKind.MstMin)}, " +
            $"LMP={((config.Projection as CouplingProjection)?.LmpRescale == true)}");
    }

    /// <summary>
    /// Load the compiled graph from the content-addressed cache under
    /// <c>{baseDirectory}/graph_cache/</c> when present; otherwise compile it,
    /// persist the artifact, and write the graph-health report into the run.
    /// </summary>
    internal static CsrGraph LoadOrBuildGraph(
        GraphRunOptions options, SpcUserDataset dataset, GraphCompilerConfig config,
        string baseDirectory, string runDirectory)
    {
        IDistanceMetric? distanceMetric = DistanceMetric(options);
        string datasetFingerprint = GraphConstructionPersistence.ComputeDatasetFingerprint(dataset, distanceMetric);
        string configFingerprint = GraphConstructionPersistence.ComputeConfigFingerprint(config);
        // Graph cache is orthogonal to the run tree: a content-addressed sibling
        // under {base}/graph_cache/, never a parent of the run.
        string datasetRoot = GraphConstructionPersistence.ChooseDatasetRoot(
            Path.Combine(baseDirectory, "graph_cache"), dataset, datasetFingerprint, configFingerprint);

        if (GraphConstructionPersistence.TryLoadGraphArtifact(datasetRoot, out var persistedGraph, out var persistedManifest))
        {
            Console.WriteLine($"Loaded persisted graph artifact from {datasetRoot}");
            Console.WriteLine($"Graph artifact : {GraphConstructionPersistence.GetGraphPath(datasetRoot)}");
            Console.WriteLine($"Graph manifest : {GraphConstructionPersistence.GetManifestPath(datasetRoot)}");
            Console.WriteLine($"Graph manifest contains {persistedManifest?.Diagnostics.Messages.Count ?? 0} diagnostic entries.");
            return persistedGraph;
        }

        if (dataset.Features.Length == 0)
            throw new ArgumentException("Dataset must contain at least one observation.");

        GraphBuildResult buildResult = GraphCompiler.Build(
            config, dataset.Features.Length, GraphMetric.FromFeatures(dataset.Features, distanceMetric),
            protectedEdges: TDA.Ph.H1CycleEdges.FromDistanceGraph);

        GraphConstructionManifest graphManifest = GraphConstructionPersistence.MaterializeManifest(
            config, buildResult, datasetFingerprint);
        GraphConstructionPersistence.WriteGraphArtifact(datasetRoot, buildResult, graphManifest);

        Console.WriteLine($"Persisted graph artifact : {GraphConstructionPersistence.GetGraphPath(datasetRoot)}");
        Console.WriteLine($"Persisted graph manifest : {GraphConstructionPersistence.GetManifestPath(datasetRoot)}");

        GraphHealthReport health = GraphHealth.Evaluate(buildResult, k: config.Topology.K ?? 0);
        string graphHealthPath = GraphHealthFile.WriteTo(runDirectory, health);
        if (!string.IsNullOrEmpty(health.Verdict.PrimaryRecommendation))
        {
            Console.WriteLine();
            Console.WriteLine($"Graph health       : {health.Verdict.PrimaryRecommendation}");
            Console.WriteLine();
        }
        Console.WriteLine($"Graph health JSON  : {graphHealthPath}");

        return buildResult.Graph;
    }

    // ── Manifest ─────────────────────────────────────────────────────────
    internal static DatasetSpec DatasetSpec(GraphRunOptions options)
        => !string.IsNullOrWhiteSpace(options.DatasetFile)
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

    /// <summary>
    /// The projection (kernel + LMP + bandwidth, or distance pass-through) is
    /// embedded as-is and round-trips polymorphically — no kernel flattening
    /// switch to keep in sync with <c>ManifestMaterialization.ToConfig</c>.
    /// </summary>
    internal static GraphSpec GraphSpec(GraphRunOptions options, GraphCompilerConfig config)
        => new GraphSpec(
            TopologyKind:    config.Topology.Kind == TopologyKind.EpsilonBall ? "epsilon" : "knn",
            FilterKind:      config.Filter.Kind == FilterKind.MutualKnn ? "mutualknn" : "or_rule",
            K:               config.Topology.K ?? 0,
            Epsilon:         config.Topology.Epsilon ?? 0.0,
            DistanceMetric:  options.DistanceMetricSpec,
            EnsureConnected: config.Repair.Kind == RepairKind.MstMin,
            Projection:      config.Projection);

    // ── Help ─────────────────────────────────────────────────────────────
    /// <summary>The dataset, output and graph sections of every graph command's help.</summary>
    internal static void PrintHelp()
    {
        Console.WriteLine("Available synthetic generators:");
        foreach (var kind in DatasetKinds) Console.WriteLine($"  {kind}");
        Console.WriteLine();
        Console.WriteLine("Dataset:");
        Console.WriteLine("  --dataset <name>              Synthetic generator (see list above)");
        Console.WriteLine("  --param <name>=<value>        Per-generator parameter; may repeat");
        Console.WriteLine("  --list-generators             Print generators and exit");
        Console.WriteLine("  --generator-schema <name>     Print a generator's parameter schema and exit");
        Console.WriteLine("  --dataset-file <path>         CSV input instead of synthetic data");
        Console.WriteLine("  --label-column <name|idx>     Label column (default: last)");
        Console.WriteLine("  --delimiter <char|tab>        CSV separator (default: ,)");
        Console.WriteLine("  --no-header                   Treat CSV as having no header row");
        Console.WriteLine("  --seed <int>                  RNG seed (default: 42)");
        Console.WriteLine();
        Console.WriteLine("Output:");
        Console.WriteLine("  --base-dir <path>             Run-root base directory (default: artifacts)");
        Console.WriteLine("  --run-name <name>             Run family directory (default: the command)");
        Console.WriteLine();
        Console.WriteLine("Graph topology:");
        Console.WriteLine("  --topology <knn|epsilon>            Graph topology kind (default: knn)");
        Console.WriteLine("  --filter <or_rule|mutualknn>        Graph filter kind (default: or_rule)");
        Console.WriteLine("  --k <int>                             KNN k (default: 10)");
        Console.WriteLine("  --epsilon <double>                    Epsilon-ball radius (epsilon mode only)");
        Console.WriteLine("  --distance-metric <spec>              euclidean | manhattan | minkowski:p=N");
        Console.WriteLine("                                        | poincare | cosine | hamming.");
        Console.WriteLine("                                        Default: euclidean (inline, no metric obj).");
        Console.WriteLine("  --ensure-connected                    MST-repair the graph (Blatt: typically on)");
        Console.WriteLine("  --lmp                                 Apply Local Mutual Proximity scaling");
        Console.WriteLine();
        Console.WriteLine("Coupling kernel:");
        Console.WriteLine("  --kernel <gaussian|cauchy|laplacian|linear>   Single-kernel coupling (default: gaussian)");
        Console.WriteLine("  --bandwidth <double>                          Kernel bandwidth (0=auto)");
        Console.WriteLine("  --bandwidth-strategy <mad|quantile-normalized|log-scale-hyperbolic|mean-edge-distance>");
        Console.WriteLine("                                                Auto-bandwidth estimator (default: from the metric)");
        Console.WriteLine("  --mixture <spec>                              Kernel mixture weights, e.g.");
        Console.WriteLine("                                                'gauss=0.5,cauchy=0.3,laplace=0.2'.");
        Console.WriteLine("                                                When set, overrides --kernel.");
        Console.WriteLine("  --mixture-bandwidth <spec>                    Per-component bandwidths, same format.");
        Console.WriteLine("                                                Omit to auto-estimate per component.");
        Console.WriteLine();
    }

    // ── Parsing primitives ───────────────────────────────────────────────
    internal static int ParseInt(string s) => int.Parse(s, CultureInfo.InvariantCulture);

    internal static double ParseDouble(string s) => double.Parse(s, CultureInfo.InvariantCulture);

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

    internal static Graphs.Distance.BandwidthStrategy ParseBandwidthStrategy(string value) => value.ToLowerInvariant() switch
    {
        "mad"                  => Graphs.Distance.BandwidthStrategy.MadConsistencyFactor,
        "quantile-normalized"  => Graphs.Distance.BandwidthStrategy.QuantileNormalized,
        "log-scale-hyperbolic" => Graphs.Distance.BandwidthStrategy.LogScaleHyperbolic,
        "mean-edge-distance"   => Graphs.Distance.BandwidthStrategy.MeanEdgeDistance,
        _ => throw new ArgumentException(
            $"Unknown bandwidth strategy '{value}'. Valid: mad, quantile-normalized, log-scale-hyperbolic, mean-edge-distance."),
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
}
