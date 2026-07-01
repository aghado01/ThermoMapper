using System;
using System.Collections.Generic;
using System.IO;
using Graphs;
using Graphs.Coupling;
using Graphs.Primitives;
using Graphs.Proximity;
using Graphs.Spectral;
using Maths.Geometry;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Euclidean;
using Synthetic.Manifolds;
using Viz;
using Viz.Adapters.Synthetic;
using Viz.Renderers;

ArtifactRun run = args.Length > 0
    ? HarnessArtifacts.Attach(
        runKind: "smoke-runs",
        suiteName: "VizCoreSmoke",
        runName: "Main",
        runDirectory: args[0],
        metadata: new Dictionary<string, object?>
        {
            ["OutputDirectoryMode"] = "user-specified",
        })
    : HarnessArtifacts.Create(
        runKind: "smoke-runs",
        suiteName: "VizCoreSmoke",
        runName: "Main",
        metadata: new Dictionary<string, object?>
        {
            ["OutputDirectoryMode"] = "default",
        });
string outDir = run.RunDirectory;
var sceneOutputs = new List<SceneArtifact>();

// ── Scene 1: Crescent + Ellipsoid (baseline) ──────────────────────────────────
RunScene(
    CrescentAndEllipsoid.Generate(),
    "Crescent + Ellipsoid (smoke test)",
    Path.Combine(outDir, "viz-smoke-crescent.html"),
    sceneOutputs);

// ── Scene 2–5: Möbius + Ellipsoid across placement modes ─────────────────────
foreach (var placement in new[] {
    MobiusPlacement.NearSeam,
    MobiusPlacement.CenterCrossOrtho,
    MobiusPlacement.PeripheralElbow })
{
    RunScene(
        MobiusAndEllipsoid.Generate(placement: placement, dimensions: 3),
        $"Möbius + Ellipsoid — {placement} (3D)",
        Path.Combine(outDir, $"viz-smoke-mobius-{placement.ToString().ToLowerInvariant()}-3d.html"),
        sceneOutputs);
}

// ── Scene 6: Möbius 4D ────────────────────────────────────────────────────────
RunScene(
    MobiusAndEllipsoid.Generate(placement: MobiusPlacement.CenterCrossOrtho, dimensions: 4),
    "Möbius + Ellipsoid — CenterCrossOrtho (4D projected)",
    Path.Combine(outDir, "viz-smoke-mobius-orthogonal-4d.html"),
    sceneOutputs);

// ── Scene 7: Möbius cross-section variants ────────────────────────────────────
foreach (var xsec in new[] {
    TubeCrossSection.UniformDisk,
    TubeCrossSection.Annular,
    TubeCrossSection.GaussianAnisotropic })
{
    RunScene(
        MobiusAndEllipsoid.Generate(crossSection: xsec, placement: MobiusPlacement.NearSeam),
        $"Möbius + Ellipsoid — {xsec}",
        Path.Combine(outDir, $"viz-smoke-mobius-{xsec.ToString().ToLowerInvariant()}.html"),
        sceneOutputs);
}

string summaryPath = run.WriteJson(
    "summary.json",
    new
    {
        SceneCount = sceneOutputs.Count,
        Scenes = sceneOutputs,
    });

Console.WriteLine($"RunRoot\t{run.RunDirectory}");
Console.WriteLine($"Manifest\t{run.ManifestPath}");
Console.WriteLine($"Summary\t{summaryPath}");

static void RunScene(
    SyntheticDataset dataset,
    string title,
    string outPath,
    List<SceneArtifact> sceneOutputs)
{
    var adapted = new SyntheticDatasetAdapter().Adapt(dataset);

    double[] features = adapted.Points.Features.ToArray();
    int n = adapted.Points.N;
    int d = adapted.Points.D;
    int[] gtLabels = ExtractGtLabels(adapted);

    // Project 4D to 3D for edge-building so distances are in 3D ambient space
    double[] edgeFeatures = d == 4
        ? FlattenRows(MobiusAndEllipsoid.Project4DTo3D(UnflattenRows(features, n, d)), n, 3)
        : features;
    int edgeDim = Math.Min(d, 3);
    double[][] displayPoints = d == 4
        ? MobiusAndEllipsoid.Project4DTo3D(UnflattenRows(features, n, d))
        : UnflattenRows(features, n, d);

    var edgeLayers = new List<EdgeLayer>
    {
        BuildEdgeLayer(edgeFeatures, n, edgeDim, gtLabels, "euclidean", new KnnSpec(7)),
        BuildEdgeLayer(edgeFeatures, n, edgeDim, gtLabels, "euclidean", new MutualKnnSpec(7)),
    };
    IReadOnlyList<LineFieldLayer> lineFieldLayers = BuildIntrinsicLineFieldLayers(displayPoints, k: 7);

    var vizDataset = new VizDataset(
        adapted.Points,
        adapted.LabelLayers,
        adapted.NodeSignalLayers,
        edgeLayers,
        adapted.GaussianLayers,
        adapted.TemporalSequences,
        adapted.SpineLayers,
        generatorParams: adapted.GeneratorParams,
        generatorParamSchema: adapted.GeneratorParamSchema,
        lineFieldLayers: lineFieldLayers);

    var package = SceneBuilder.Build(vizDataset, new SceneDescriptor
    {
        Title = title,
        Hints = SceneRenderHints.Default,
    });

    using var stream = File.Create(outPath);
    new ThreeJsHtmlRenderTarget().Render(package, stream);
    sceneOutputs.Add(new SceneArtifact(title, outPath, adapted.Points.N, adapted.Points.D, edgeLayers.Count));
}

static double[][] UnflattenRows(double[] flat, int n, int d)
{
    var rows = new double[n][];
    for (int i = 0; i < n; i++)
    {
        rows[i] = new double[d];
        Array.Copy(flat, i * d, rows[i], 0, d);
    }
    return rows;
}

static double[] FlattenRows(double[][] rows, int n, int d)
{
    var flat = new double[n * d];
    for (int i = 0; i < n; i++)
        Array.Copy(rows[i], 0, flat, i * d, d);
    return flat;
}

static IReadOnlyList<LineFieldLayer> BuildIntrinsicLineFieldLayers(double[][] points, int k)
{
    int n = points.Length;
    int d = points[0].Length;
    CsrGraph graph = GraphCompiler.Build(
        new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = k },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.MstMin },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection { Kernel = new Gaussian(0.0), LmpRescale = false },
        },
        n,
        new GraphMetric((i, j) => EuclideanDistPair(points[i], points[j]))).Graph;

    int[][] adjacency = BuildAdjacency(graph);
    double[] localTangents = LocalTangent.Compute(points, adjacency);

    return new[]
    {
        new LineFieldLayer(
            "SpectralBridge",
            SpectralBridge.LineFieldFromFiedler(graph, points),
            n,
            d,
            LineFieldSource.SpectralGradient),
        new LineFieldLayer(
            "LocalTangent",
            localTangents,
            n,
            d,
            LineFieldSource.LocalPca),
    };
}

static int[][] BuildAdjacency(CsrGraph graph)
{
    var adjacency = new int[graph.NodeCount][];
    for (int i = 0; i < graph.NodeCount; i++)
    {
        int start = graph.RowPointers[i];
        int length = graph.RowPointers[i + 1] - start;
        var neighbors = new int[length];
        Array.Copy(graph.Targets, start, neighbors, 0, length);
        adjacency[i] = neighbors;
    }

    return adjacency;
}

// ── Helpers ───────────────────────────────────────────────────────────────────

static int[] ExtractGtLabels(VizDataset ds)
{
    // Prefer the GroundTruth label layer; fall back to first available.
    var layer = ds.LabelLayers.Count > 0 ? ds.LabelLayers[0] : null;
    return layer?.Labels.ToArray() ?? Array.Empty<int>();
}

static EdgeLayer BuildEdgeLayer(
    double[] features, int n, int d, int[] gtLabels,
    string metric, ProximitySpec proximity)
{
    if (!string.Equals(metric, "euclidean", StringComparison.OrdinalIgnoreCase))
        throw new NotSupportedException($"Metric {metric} not yet wired in harness");

    Func<int, int, double> dist = (i, j) => EuclideanDist(features, d, i, j);

    NeighborSelection sel = proximity switch
    {
        KnnSpec s => ProximityGraph.SelectKnn(n, s.K, dist),
        MutualKnnSpec s => ProximityGraph.SelectMutualKnn(n, s.K, dist),
        _ => throw new NotSupportedException($"NeighborRule {proximity.Kind} not yet wired in harness"),
    };

    // Flatten per-node neighbor lists → parallel edge arrays
    int edgeCount = 0;
    foreach (var row in sel.AllNeighbors) edgeCount += row.Length;

    var src = new int[edgeCount];
    var dst = new int[edgeCount];
    var weight = new double[edgeCount];
    int idx = 0;

    for (int i = 0; i < n; i++)
        foreach (var nb in sel.AllNeighbors[i])
        {
            src[idx] = i;
            dst[idx] = nb.Index;
            weight[idx] = nb.Distance;
            idx++;
        }

    // False-bridge arrays: GT cluster of each edge endpoint.
    // Renderer colors edge magenta when edgeClusterSrc[e] != edgeClusterDst[e].
    int[]? edgeClusterSrc = null;
    int[]? edgeClusterDst = null;
    if (gtLabels.Length > 0)
    {
        edgeClusterSrc = new int[edgeCount];
        edgeClusterDst = new int[edgeCount];
        for (int e = 0; e < edgeCount; e++)
        {
            edgeClusterSrc[e] = gtLabels[src[e]];
            edgeClusterDst[e] = gtLabels[dst[e]];
        }
    }

    string name = LayerName(metric, proximity);
    return new EdgeLayer(name, src, dst, weight, metric, proximity, edgeClusterSrc, edgeClusterDst);
}

static double EuclideanDist(double[] features, int d, int i, int j)
{
    double sum = 0;
    for (int k = 0; k < d; k++)
    {
        double diff = features[i * d + k] - features[j * d + k];
        sum += diff * diff;
    }
    return Math.Sqrt(sum);
}

static double EuclideanDistPair(double[] left, double[] right)
{
    double sum = 0;
    for (int i = 0; i < left.Length; i++)
    {
        double diff = left[i] - right[i];
        sum += diff * diff;
    }

    return Math.Sqrt(sum);
}

static string LayerName(string metric, ProximitySpec proximity)
{
    string m = metric.ToLowerInvariant();
    string p = proximity switch
    {
        KnnSpec s => $"knn:k={s.K}",
        MutualKnnSpec s => $"mutual_knn:k={s.K}",
        EpsilonBallSpec s => $"epsilon_ball:eps={s.Epsilon:G4}",
        MstAugmentedSpec s => $"mst_aug:k={s.K}",
        _ => proximity.Kind.ToString().ToLowerInvariant(),
    };
    return $"{m}:{p}";
}

internal sealed record SceneArtifact(
    string Title,
    string OutputPath,
    int PointCount,
    int Dimension,
    int EdgeLayerCount);
