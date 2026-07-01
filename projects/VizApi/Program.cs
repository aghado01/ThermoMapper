using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Clustering.Statistical.GMM;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Primitives;
using Graphs.Proximity;
using Graphs.Spectral;
using Maths.Geometry;
using Microsoft.AspNetCore.Builder;
using Microsoft.AspNetCore.Http;
using Microsoft.Extensions.Hosting;
using Synthetic;
using Synthetic.Euclidean;
using Synthetic.Manifolds;
using Viz;
using VizApi;
using Viz.Adapters.Gmm;
using Viz.Adapters.Synthetic;
using Viz.Renderers;

var builder = WebApplication.CreateBuilder(args);
var app = builder.Build();

// ── GET / → initial viewer HTML with default scene ────────────────────────────
app.MapGet("/", () =>
{
    var package = BuildPackage(new RegenRequest());
    using var ms = new MemoryStream();
    new ThreeJsHtmlRenderTarget().Render(package, ms);
    string html = Encoding.UTF8.GetString(ms.ToArray());
    return Results.Content(html, "text/html; charset=utf-8");
});

// ── POST /api/regen → ScenePackage JSON ───────────────────────────────────────
// Receives generator params from viewer.html; returns full ScenePackageJson.
// Dispatches on req.Generator; add new generators to the switch below.
app.MapPost("/api/regen", (RegenRequest req) =>
{
    var package = BuildPackage(req);
    using var ms = new MemoryStream();
    new JsonExportRenderTarget(compact: true).Render(package, ms);
    string json = Encoding.UTF8.GetString(ms.ToArray());
    return Results.Content(json, "application/json");
});

app.MapGet("/api/metrics", () => Results.Json(
    MetricRegistry.Available.Select(metric => new
    {
        id = metric.Id,
        displayName = metric.DisplayName,
        requiresParameter = metric.RequiresParameter,
        parameterName = metric.ParameterName,
    })));

app.Run();

// ── Scene construction ────────────────────────────────────────────────────────

// Phase 1 dispatch: generator-specific synthesis + adaptation only.
// Phase 2 (shared): edges, GMM overlays, mergedParams, scene build — all in BuildPackage.
static ScenePackage BuildPackage(RegenRequest req)
{
    var (adapted, title) = req.Generator switch
    {
        "MobiusAndEllipsoid" => BuildMobiusAdapted(req),
        "HyperbolicBlobs" => BuildHyperbolicBlobsAdapted(req),
        "HyperbolicBlattHierarchy" => BuildHyperbolicBlattHierarchyAdapted(req),
        "Simplex" => BuildSimplexAdapted(req),
        "GaussianManifold" => BuildGaussianManifoldAdapted(req),
        "TwoMoons" => BuildTwoMoonsAdapted(req),
        "BlattHierarchy" => BuildBlattHierarchyAdapted(req),
        "BlattThreeCluster" => BuildBlattThreeClusterAdapted(req),
        _ => BuildCrescentAdapted(req),
    };

    // ── Shared: proximity graph ───────────────────────────────────────────────
    double[] features = adapted.Points.Features.ToArray();
    int n = adapted.Points.N;
    int d = adapted.Points.D;
    int[] gtLabels = ExtractGtLabels(adapted);

    string metricSpec = ResolveCompatibleMetric(req.Generator, req.Metric);
    KernelType kernelType = Enum.TryParse<KernelType>(req.Kernel, out var kt) ? kt : KernelType.Gaussian;
    TopologyKind topologyKind = req.NeighborRule == "EpsilonBall"
        ? TopologyKind.EpsilonBall
        : TopologyKind.Knn;
    FilterKind filterKind = req.NeighborRule == "MutualKnn"
        ? FilterKind.MutualKnn
        : FilterKind.OrRule;
    // MST is a connectivity modifier — it bridges disconnected components on
    // top of any rule. Used to be a fake 4th rule; now it's an orthogonal flag.
    bool ensureConnected = req.EnsureConnected;

    // ETL/manifold adaptation is quarantined for now. The active synthetic
    // generators already emit graph-ready feature coordinates for the metrics
    // currently exposed by the viewer, so build distances directly from them.
    Func<int, int, double> dist = DistanceFactory.Create(features, d, metricSpec);

    IKernelDescriptor kernelDescriptor = kernelType switch
    {
        KernelType.Gaussian => new Gaussian(req.Bandwidth),
        KernelType.Cauchy => new Cauchy(req.Bandwidth),
        KernelType.Laplacian => new Laplacian(req.Bandwidth),
        KernelType.Linear => new Linear(req.Bandwidth),
        _ => new Gaussian(req.Bandwidth),
    };

    var graphConfig = new GraphCompilerConfig
    {
        Topology = topologyKind == TopologyKind.EpsilonBall
            ? new TopologyConfig { Kind = TopologyKind.EpsilonBall, Epsilon = req.EpsilonBallEpsilon }
            : new TopologyConfig { Kind = TopologyKind.Knn, K = req.KnnK },
        Filter = new FilterConfig
        {
            Kind = filterKind,
            MutualBandwidthSource = MutualBandwidthSource.DirectedKth,
        },
        Repair = new RepairConfig
        {
            Kind = ensureConnected ? RepairKind.MstMin : RepairKind.NoRepair,
        },
        Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
        Projection = new CouplingProjection { Kernel = kernelDescriptor, LmpRescale = false },
    };

    CsrGraph csrGraph = GraphCompiler.Build(graphConfig, n, new GraphMetric(dist)).Graph;

    // ── Edge layer: serialize CsrGraph — weights are coupling strengths in (0,1] ──
    // MST is a connectivity modifier (req.EnsureConnected), not its own rule —
    // the spec captures the primary rule only. The applied modifiers travel
    // back to the viewer via generator_params so they can be shown in the UI.
    ProximitySpec spec = req.NeighborRule switch
    {
        "MutualKnn" => new MutualKnnSpec(req.KnnK),
        "EpsilonBall" => new EpsilonBallSpec(req.EpsilonBallEpsilon),
        _ => new KnnSpec(req.KnnK),
    };
    var edgeLayers = new List<EdgeLayer> { BuildEdgeLayerFromCsr(csrGraph, n, gtLabels, metricSpec, spec) };

    // ── Wing-2: local tangent flow — reuses same CsrGraph adjacency ──────────
    var points2d = new double[n][];
    for (int i = 0; i < n; i++)
    {
        points2d[i] = new double[d];
        Array.Copy(features, i * d, points2d[i], 0, d);
    }
    double[][] displayPoints = d == 4 ? MobiusAndEllipsoid.Project4DTo3D(points2d) : points2d;
    int fieldDim = displayPoints[0].Length;

    var adjacency = new int[n][];
    for (int i = 0; i < n; i++)
    {
        int rowStart = csrGraph.RowPointers[i];
        int rowEnd = csrGraph.RowPointers[i + 1];
        int len = rowEnd - rowStart;
        adjacency[i] = new int[len];
        for (int j = 0; j < len; j++)
            adjacency[i][j] = csrGraph.Targets[rowStart + j];
    }
    var tangentVectors = LocalTangent.Compute(displayPoints, adjacency);
    var localTangentLineField = (double[])tangentVectors.Clone();
    var lineFieldLayers = new List<LineFieldLayer>(adapted.LineFieldLayers)
    {
        new LineFieldLayer(
            "SpectralBridge",
            SpectralBridge.LineFieldFromFiedler(csrGraph, displayPoints),
            n,
            fieldDim,
            LineFieldSource.SpectralGradient),
        new LineFieldLayer(
            "LocalTangent",
            localTangentLineField,
            n,
            fieldDim,
            LineFieldSource.LocalPca),
    };

    // BFS orientation propagation: aligns tangent signs across the graph so
    // incoherence is a data signal (manifold discontinuity, false bridges) not
    // an artefact of independent power-iteration sign choices.
    LocalTangent.PropagateOrientation(
        tangentVectors, n, fieldDim, csrGraph.RowPointers, csrGraph.Targets);

    // Per-point coherence: mean dot(t_i, t_j) over CSR neighbours.
    // High = intrinsic edge; near-zero = ambient shortcut or degenerate neighbourhood.
    double[] coherence = LocalTangent.ComputeCoherence(
        tangentVectors, n, fieldDim, csrGraph.RowPointers, csrGraph.Targets);

    var vectorFieldLayers = new List<VectorFieldLayer>
    {
        new VectorFieldLayer("Local PCA Flow", tangentVectors, n, fieldDim),
    };

    // ── Shared: GMM overlay (single layer, gated by GmmMode) ─────────────────
    // The viewer picks Oracle (spine-aware, analytic σ) or EM (real fit on the
    // point cloud) via a dropdown. We emit exactly one GaussianLayer so the two
    // modes don't visually pile up. ComponentToClusterMap is populated so each
    // ellipsoid is coloured by the cluster its points belong to, matching the
    // point-cloud palette.
    int gmmPos = Math.Clamp((int)Math.Round(req.GmmComponents), 0, 4);
    int gmmK = 1 << gmmPos;
    var gaussianLayers = new List<GaussianLayer>();

    bool useEm = string.Equals(req.GmmMode, "EM", StringComparison.OrdinalIgnoreCase);
    if (useEm)
    {
        var fittedGmm = new GaussianMixtureModel(k: gmmK, dimension: d);
        fittedGmm.RobustInitialize(points2d);
        fittedGmm.Fit(points2d);
        int[] ctcMap = ComponentToClusterMap(fittedGmm, points2d, gtLabels);
        gaussianLayers.Add(BuildGaussianLayerFromGmm(fittedGmm, $"GMM (EM) K={gmmK}", ctcMap));
    }
    else
    {
        foreach (var gmmLayer in BuildGmmOverlaysFromSpines(adapted.SpineLayers, gmmK))
            gaussianLayers.Add(gmmLayer);
    }

    var nodeSignalLayers = new List<NodeSignalLayer>(adapted.NodeSignalLayers)
    {
        new NodeSignalLayer("Flow Coherence", coherence, ScalarSource.CoherenceScore),
    };

    // ── Shared: merge params so they round-trip without resetting on regen ────
    var mergedParams = adapted.GeneratorParams != null
        ? new Dictionary<string, object>(adapted.GeneratorParams)
        : new Dictionary<string, object>();
    mergedParams["generator"] = req.Generator;
    mergedParams["seed"] = req.Seed;
    mergedParams["knnK"] = req.KnnK;
    mergedParams["metric"] = metricSpec;
    mergedParams["neighborRule"] = req.NeighborRule;
    mergedParams["epsilonBallEpsilon"] = req.EpsilonBallEpsilon;
    mergedParams["kernel"] = req.Kernel;
    mergedParams["bandwidth"] = req.Bandwidth;
    mergedParams["ensureConnected"] = req.EnsureConnected;
    mergedParams["gmmComponents"] = req.GmmComponents;
    mergedParams["gmmMode"] = req.GmmMode;
    mergedParams["showFlow"] = req.ShowFlow;
    mergedParams["flowMode"] = req.FlowMode;

    var vizDataset = new VizDataset(
        adapted.Points,
        adapted.LabelLayers,
        nodeSignalLayers,
        edgeLayers,
        gaussianLayers,
        adapted.TemporalSequences,
        adapted.SpineLayers,
        generatorParams: mergedParams,
        generatorParamSchema: adapted.GeneratorParamSchema ?? SchemaCatalog.ForGenerator(req.Generator),
        vectorFieldLayers: vectorFieldLayers,
        lineFieldLayers: lineFieldLayers);

    return SceneBuilder.Build(vizDataset, new SceneDescriptor
    {
        Title = title,
        Hints = new SceneRenderHints
        {
            ShowVectorField = req.ShowFlow,
        },
    });
}

static string ResolveCompatibleMetric(string generator, string? requested)
{
    const string fallback = "euclidean";

    if (TryNormalizeIntegratedMetricSpec(requested, out string integratedMetric))
    {
        return integratedMetric switch
        {
            "poincare" when !string.Equals(generator, "HyperbolicBlobs", StringComparison.OrdinalIgnoreCase)
                && !string.Equals(generator, "HyperbolicBlattHierarchy", StringComparison.OrdinalIgnoreCase)
                => fallback,
            _ => integratedMetric,
        };
    }

    if (!Enum.TryParse(requested, true, out VizApiExtraMetric requestedMetric))
        return fallback;

    return requestedMetric switch
    {
        VizApiExtraMetric.FisherRaoSimplex when !string.Equals(generator, "Simplex", StringComparison.OrdinalIgnoreCase)
            => fallback,
        VizApiExtraMetric.FisherRaoSimplex => nameof(VizApiExtraMetric.FisherRaoSimplex),
        VizApiExtraMetric.FisherRaoHalfPlane when !string.Equals(generator, "GaussianManifold", StringComparison.OrdinalIgnoreCase)
            => fallback,
        VizApiExtraMetric.FisherRaoHalfPlane => nameof(VizApiExtraMetric.FisherRaoHalfPlane),
        VizApiExtraMetric.Wasserstein when !string.Equals(generator, "Simplex", StringComparison.OrdinalIgnoreCase)
            => fallback,
        VizApiExtraMetric.Wasserstein => nameof(VizApiExtraMetric.Wasserstein),
        VizApiExtraMetric.Jaccard => nameof(VizApiExtraMetric.Jaccard),
        _ => fallback,
    };
}

static bool TryNormalizeIntegratedMetricSpec(string? requested, out string metricSpec)
{
    if (string.IsNullOrWhiteSpace(requested))
    {
        metricSpec = string.Empty;
        return false;
    }

    if (string.Equals(requested, "Minkowski", StringComparison.OrdinalIgnoreCase))
    {
        metricSpec = "minkowski:p=2";
        return true;
    }

    try
    {
        metricSpec = FormatMetricSpec(MetricRegistry.Parse(requested));
        return true;
    }
    catch (ArgumentException)
    {
        metricSpec = string.Empty;
        return false;
    }
}

static string FormatMetricSpec(MetricSpec spec)
{
    return spec.Kind switch
    {
        MetricKind.Euclidean => "euclidean",
        MetricKind.Manhattan => "manhattan",
        MetricKind.Minkowski => $"minkowski:p={spec.MinkowskiP.ToString(CultureInfo.InvariantCulture)}",
        MetricKind.Hamming => "hamming",
        MetricKind.Poincare => "poincare",
        MetricKind.Cosine => "cosine",
        _ => throw new ArgumentOutOfRangeException(nameof(spec), spec.Kind, "Unhandled metric kind."),
    };
}

// ── GMM overlay: generator-agnostic, driven purely by spine geometry ──────────
// For each SpineLayer, auto-derives:
//   sigmaLong = arcLength / K           (one arc-segment width)
//   sigmaPerp = TypicalScale × 0.9      (from generator) or stepSize × 2 (fallback)
static IEnumerable<GaussianLayer> BuildGmmOverlaysFromSpines(
    IReadOnlyList<SpineLayer> spineLayers, int K)
{
    foreach (var spine in spineLayers)
    {
        var ss = spine.SpineSamples;
        int M = ss.Length;
        if (M < 2) continue;

        double arcLength = 0;
        for (int i = 1; i < M; i++)
        {
            double dx = ss[i][0] - ss[i - 1][0];
            double dy = ss[i][1] - ss[i - 1][1];
            double dz = ss[i].Length > 2 ? ss[i][2] - ss[i - 1][2] : 0.0;
            arcLength += Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        double stepSize = arcLength / (M - 1);
        double sigmaLong = arcLength / K;
        double sigmaPerp = spine.TypicalScale > 0 ? spine.TypicalScale * 0.9 : stepSize * 2.0;

        int stride = Math.Max(1, M / K);
        int offset = stride / 2;
        yield return BuildGmmOverlayLayer(
            spine, sigmaLong, sigmaPerp, stride, offset,
            $"GMM (Oracle) K={K}", spine.ClusterIdx);
    }
}

// ── EM → cluster mapping ──────────────────────────────────────────────────────
// Per-component dominant GT cluster, used to colour ellipsoids in lockstep with
// the underlying point cloud. Hard-assigns points via gmm.Predict, then takes
// the mode of GT labels per component. Empty components fall back to component
// index so legacy palette behaviour is preserved.
static int[] ComponentToClusterMap(GaussianMixtureModel gmm, double[][] points, int[] gtLabels)
{
    int k = gmm.K;
    int[] map = new int[k];
    if (gtLabels.Length == 0)
    {
        for (int i = 0; i < k; i++) map[i] = i;
        return map;
    }

    int[] preds = gmm.Predict(points);

    int maxLabel = 0;
    for (int i = 0; i < gtLabels.Length; i++)
        if (gtLabels[i] > maxLabel) maxLabel = gtLabels[i];
    int labelSpan = maxLabel + 1;

    var counts = new int[k * labelSpan];
    for (int i = 0; i < preds.Length; i++)
    {
        int c = preds[i];
        int lbl = gtLabels[i];
        if (lbl < 0) continue;
        counts[c * labelSpan + lbl]++;
    }

    for (int c = 0; c < k; c++)
    {
        int bestLabel = c;  // fallback: empty component keeps component index
        int bestCount = 0;
        for (int lbl = 0; lbl < labelSpan; lbl++)
        {
            int v = counts[c * labelSpan + lbl];
            if (v > bestCount) { bestCount = v; bestLabel = lbl; }
        }
        map[c] = bestLabel;
    }
    return map;
}

static GaussianLayer BuildGaussianLayerFromGmm(GaussianMixtureModel gmm, string name, int[]? ctcMap)
{
    int k = gmm.K;
    int d = gmm.Dimension;
    var means = new double[k * d];
    var covs = new double[k * d * d];
    var weights = new double[k];
    for (int ki = 0; ki < k; ki++)
    {
        var comp = gmm.Components[ki];
        weights[ki] = comp.Weight;
        for (int dim = 0; dim < d; dim++) means[ki * d + dim] = comp.Mean[dim];
        for (int r = 0; r < d; r++)
            for (int c = 0; c < d; c++)
                covs[ki * d * d + r * d + c] = comp.Covariance[r, c];
    }
    return new GaussianLayer(name, means, covs, weights, k, d, ctcMap);
}

// ── Phase 1: Crescent synthesis + adaptation (no edges, no overlay) ──────────
static (VizDataset adapted, string title) BuildCrescentAdapted(RegenRequest req)
{
    EllipsoidPlacement placement = Enum.TryParse<EllipsoidPlacement>(req.Placement, out var p)
        ? p
        : EllipsoidPlacement.OrthogonalElbowIntersect;

    EllipsoidShellMode shellMode = Enum.TryParse<EllipsoidShellMode>(req.EllipsoidShellMode, out var sm)
        ? sm : EllipsoidShellMode.Solid;

    var dataset = CrescentAndEllipsoid.Generate(
        crescentPoints: req.CrescentPoints,
        crescentRadius: req.CrescentRadius,
        crescentWidth: req.CrescentWidth,
        arcHalfAngle: req.ArcHalfAngle,
        ellipsoidPoints: req.EllipsoidPoints,
        ellipsoidAxes: req.EllipsoidAxes,
        ellipsoidShellMode: shellMode,
        placement: placement,
        intersectDepth: req.IntersectDepth,
        intersectRadialShift: req.IntersectRadialShift,
        gapScale: req.GapScale,
        seed: req.Seed);

    return (new SyntheticDatasetAdapter().Adapt(dataset), "Crescent + Ellipsoid");
}

// ── Phase 1: Möbius synthesis + adaptation (no edges, no overlay) ────────────
static (VizDataset adapted, string title) BuildMobiusAdapted(RegenRequest req)
{
    MobiusPlacement placement = Enum.TryParse<MobiusPlacement>(req.Placement, out var mp)
        ? mp
        : MobiusPlacement.CenterCrossOrtho;

    TubeCrossSection crossSection = Enum.TryParse<TubeCrossSection>(req.CrossSection, out var cs)
        ? cs
        : TubeCrossSection.GaussianIsotropic;

    MobiusSpineShape spineShape = Enum.TryParse<MobiusSpineShape>(req.SpineShape, out var ss)
        ? ss
        : MobiusSpineShape.Circle;

    EllipsoidShellMode shellMode = Enum.TryParse<EllipsoidShellMode>(req.EllipsoidShellMode, out var sm)
        ? sm : EllipsoidShellMode.Solid;

    var dataset = MobiusAndEllipsoid.Generate(
        mobiusPoints: req.MobiusPoints,
        spineRadius: req.SpineRadius,
        halfWidth: req.HalfWidth,
        halfThickness: req.HalfThickness,
        noiseSigma: req.NoiseSigma,
        twistCount: req.TwistCount,
        crossSection: crossSection,
        radialBias: req.RadialBias,
        spineShape: spineShape,
        splayFactor: req.SplayFactor,
        ellipsoidPoints: req.EllipsoidPoints,
        ellipsoidAxes: req.EllipsoidAxes,
        ellipsoidShellMode: shellMode,
        placement: placement,
        intersectDepth: req.IntersectDepth,
        intersectRadialShift: req.IntersectRadialShift,
        gapScale: req.GapScale,
        dimensions: req.Dimensions,
        seed: req.Seed);

    return (new SyntheticDatasetAdapter().Adapt(dataset), "Möbius + Ellipsoid");
}

// ── Phase 1: HyperbolicBlobs synthesis + adaptation ─────────────────────────
// Pairs with the Poincaré metric — points live strictly inside the open unit
// ball so the hyperbolic geodesic is exact rather than saturated to 1e9.
static (VizDataset adapted, string title) BuildHyperbolicBlobsAdapted(RegenRequest req)
{
    var dataset = HyperbolicBlobs.Generate(
        clusterCount: req.ClusterCount,
        pointsPerCluster: req.PointsPerCluster,
        dimensions: req.Dimensions,
        separation: req.Separation,
        spread: req.Spread,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Hyperbolic Blobs");
}

// ── Phase 1: HyperbolicBlattHierarchy synthesis + adaptation ───────────────
// Hierarchical data directly in the Poincare ball. Meant for testing whether
// graph phase transitions peel off nested semantic or structural scales.
static (VizDataset adapted, string title) BuildHyperbolicBlattHierarchyAdapted(RegenRequest req)
{
    var dataset = HyperbolicHierarchy.Generate(
        nPoints: req.HierarchyPoints,
        hierarchyDepth: req.HierarchyDepth,
        branchesPerNode: req.BranchesPerNode,
        basePointsPerLeaf: req.BasePointsPerLeaf,
        dimensions: req.Dimensions,
        baseSeparation: req.Separation,
        radiusDecay: req.RadiusDecay,
        noiseScale: req.Spread,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Hyperbolic Blatt Hierarchy");
}

// ── Phase 1: Simplex synthesis + adaptation ─────────────────────────────────
// Pairs with Fisher–Rao (simplex). Each point is a probability vector. With
// disjointSupports = true, clusters concentrate on different category ranges;
// the first three categories project as XYZ for the 3D viewer.
static (VizDataset adapted, string title) BuildSimplexAdapted(RegenRequest req)
{
    var dataset = Simplex.Generate(
        clusterCount: req.ClusterCount,
        pointsPerCluster: req.PointsPerCluster,
        categories: req.Categories,
        disjointSupports: req.DisjointSupports,
        concentration: req.Concentration,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Simplex (probability vectors)");
}

// ── Phase 1: GaussianManifold synthesis + adaptation ────────────────────────
// Pairs with Fisher–Rao (half-plane). Each point is (μ, log σ) on the 1D
// Gaussian statistical manifold. Euclidean and Fisher–Rao topologies disagree —
// the storytelling case for why information geometry matters.
static (VizDataset adapted, string title) BuildGaussianManifoldAdapted(RegenRequest req)
{
    var dataset = GaussianManifold.Generate(
        clusterCount: req.ClusterCount,
        pointsPerCluster: req.PointsPerCluster,
        clusterRadius: req.ClusterRadius,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Gaussian Manifold (μ, log σ)");
}

// ── Phase 1: TwoMoons synthesis + adaptation ───────────────────────────────
// Canonical non-convex 2-cluster toy for graph / topology diagnostics.
static (VizDataset adapted, string title) BuildTwoMoonsAdapted(RegenRequest req)
{
    var dataset = TwoMoons.Generate(
        pointsPerMoon: req.PointsPerMoon,
        noise: req.Noise,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Two Moons");
}

// ── Phase 1: BlattHierarchy synthesis + adaptation ─────────────────────────
// Three-level Gaussian hierarchy from Blatt, Wiseman, Domany (1997).
static (VizDataset adapted, string title) BuildBlattHierarchyAdapted(RegenRequest req)
{
    var dataset = BlattHierarchy.Generate(
        coarseClusters: req.CoarseClusters,
        mediumPerCoarse: req.MediumPerCoarse,
        finePerMedium: req.FinePerMedium,
        pointsPerFine: req.PointsPerFine,
        coarseSeparation: req.CoarseSeparation,
        mediumSeparation: req.MediumSeparation,
        fineSeparation: req.FineSeparation,
        leafSpread: req.LeafSpread,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Blatt Hierarchy");
}

// ── Phase 1: BlattThreeCluster synthesis + adaptation ──────────────────────
// Canonical 3-Gaussian SPC validation toy from Blatt, Wiseman, Domany (1996).
static (VizDataset adapted, string title) BuildBlattThreeClusterAdapted(RegenRequest req)
{
    var dataset = BlattThreeCluster.Generate(
        pointsPerCluster: req.PointsPerCluster,
        stdDev: req.StdDev,
        seed: req.Seed);
    return (new SyntheticDatasetAdapter().Adapt(dataset), "Blatt Three Cluster");
}

static int[] ExtractGtLabels(VizDataset ds)
{
    var layer = ds.LabelLayers.Count > 0 ? ds.LabelLayers[0] : null;
    return layer?.Labels.ToArray() ?? Array.Empty<int>();
}

// ── Proximity graph helpers ───────────────────────────────────────────────────

// Serialize a CsrGraph into an EdgeLayer.
// Weights are coupling strengths in (0,1] — high weight = high similarity = bright.
// Iterates upper triangle only (j > i) to emit each undirected edge once.
static EdgeLayer BuildEdgeLayerFromCsr(
    CsrGraph graph, int n, int[] gtLabels,
    string metricSpec, ProximitySpec proximity)
{
    // Count edges in upper triangle
    int edgeCount = 0;
    for (int i = 0; i < n; i++)
        for (int idx = graph.RowPointers[i]; idx < graph.RowPointers[i + 1]; idx++)
            if (graph.Targets[idx] > i) edgeCount++;

    var src = new int[edgeCount];
    var dst = new int[edgeCount];
    var weight = new double[edgeCount];
    int e = 0;

    for (int i = 0; i < n; i++)
        for (int idx = graph.RowPointers[i]; idx < graph.RowPointers[i + 1]; idx++)
        {
            int j = graph.Targets[idx];
            if (j <= i) continue;
            src[e] = i;
            dst[e] = j;
            weight[e] = graph.Weights[idx];  // coupling strength, not distance
            e++;
        }

    if (e != edgeCount)
        throw new InvalidOperationException($"Edge layer fill mismatch: expected {edgeCount}, filled {e}.");

    int[]? edgeClusterSrc = null;
    int[]? edgeClusterDst = null;
    if (gtLabels.Length > 0)
    {
        edgeClusterSrc = new int[edgeCount];
        edgeClusterDst = new int[edgeCount];
        for (int ei = 0; ei < edgeCount; ei++)
        {
            edgeClusterSrc[ei] = gtLabels[src[ei]];
            edgeClusterDst[ei] = gtLabels[dst[ei]];
        }
    }

    string proximityName = proximity switch
    {
        KnnSpec s => $"knn:k={s.K}",
        MutualKnnSpec s => $"mutual_knn:k={s.K}",
        EpsilonBallSpec s => $"epsilon_ball:eps={s.Epsilon:G4}",
        MstAugmentedSpec s => $"mst_aug:k={s.K}",
        _ => proximity.Kind.ToString().ToLowerInvariant(),
    };

    string name = metricSpec.ToLowerInvariant() + ":" + proximityName;
    return new EdgeLayer(name, src, dst, weight, metricSpec, proximity, edgeClusterSrc, edgeClusterDst);
}

// ── GMM overlay helper ───────────────────────────────────────────────────────
// Places K ellipsoidal Gaussians evenly along a spine, with the long axis
// aligned to the spine tangent and short axes equal to sigmaPerp.
// offset centres the components within their arc segments (pass stride/2).
static GaussianLayer BuildGmmOverlayLayer(
    SpineLayer spineLayer,
    double sigmaLong,
    double sigmaPerp,
    int stride = 4,
    int offset = 0,
    string name = "GMM K=1",
    int clusterIdx = 0)
{
    var ss = spineLayer.SpineSamples;   // M × D (D=3 for crescent in XY)
    int M = ss.Length;
    var indices = new List<int>();
    for (int i = offset; i < M; i += stride)
        indices.Add(i);

    int K = indices.Count;
    const int geomDim = 3;
    var means = new double[K * geomDim];
    var covs = new double[K * geomDim * geomDim];
    var weights = new double[K];
    double w = 1.0 / K;

    // Oracle ellipsoids all belong to the same physical cluster (the spine's),
    // so the ctc map is constant. This makes the viewer colour every Oracle
    // ellipsoid with that cluster's palette entry.
    var ctcMap = new int[K];
    for (int i = 0; i < K; i++) ctcMap[i] = clusterIdx;

    for (int ci = 0; ci < K; ci++)
    {
        int i = indices[ci];
        int j = Math.Min(i + stride, M - 1);

        // Tangent: direction along spine.
        double tx = ss[j][0] - ss[i][0];
        double ty = ss[j][1] - ss[i][1];
        double tz = ss[j].Length > 2 ? ss[j][2] - ss[i][2] : 0.0;
        double tLen = Math.Sqrt(tx * tx + ty * ty + tz * tz) + 1e-12;
        tx /= tLen; ty /= tLen; tz /= tLen;

        // Binormal approximation: crescent lives in XY plane, so B ≈ Z.
        // N = B × T, then re-orthogonalise B = T × N.
        double nx = -ty * 1.0;  // cross([0,0,1], T) = [-Ty, Tx, 0]
        double ny = tx * 1.0;
        double nz = 0.0;
        double nLen = Math.Sqrt(nx * nx + ny * ny) + 1e-12;
        nx /= nLen; ny /= nLen;

        double bx = ty * nz - tz * ny;
        double by = tz * nx - tx * nz;
        double bz = tx * ny - ty * nx;  // should be ~1.0

        // Covariance: Σ = R * diag(σl², σp², σp²) * Rᵀ  where R = [T | N | B]
        double sl2 = sigmaLong * sigmaLong;
        double sp2 = sigmaPerp * sigmaPerp;
        double[] T3 = { tx, ty, tz };
        double[] N3 = { nx, ny, nz };
        double[] B3 = { bx, by, bz };

        int covBase = ci * geomDim * geomDim;
        for (int r = 0; r < geomDim; r++)
            for (int c = 0; c < geomDim; c++)
                covs[covBase + r * geomDim + c] =
                    T3[r] * T3[c] * sl2 + N3[r] * N3[c] * sp2 + B3[r] * B3[c] * sp2;

        means[ci * geomDim] = ss[i][0];
        means[ci * geomDim + 1] = ss[i].Length > 1 ? ss[i][1] : 0.0;
        means[ci * geomDim + 2] = ss[i].Length > 2 ? ss[i][2] : 0.0;
        weights[ci] = w;
    }

    return new GaussianLayer(name, means, covs, weights, K, geomDim, ctcMap);
}

// ── Request model ─────────────────────────────────────────────────────────────
// JSON property names are camelCase (default System.Text.Json / ASP.NET Core).
// Generator-specific fields default to sensible values; the discriminator field
// "generator" determines which BuildXxx function is called.

record RegenRequest(
    // ── Discriminator ──────────────────────────────────────────────────────────
    string Generator = "CrescentAndEllipsoid",
    // ── CrescentAndEllipsoid ───────────────────────────────────────────────────
    int CrescentPoints = 10000,
    double CrescentRadius = 3.0,
    double CrescentWidth = 0.40,
    double ArcHalfAngle = 2.04,
    double GmmComponents = 1,   // slider pos: 0→K=1, 1→K=2, 2→K=4, 3→K=8, 4→K=16
                                // ── MobiusAndEllipsoid ─────────────────────────────────────────────────────
    int MobiusPoints = 10000,
    double SpineRadius = 2.5,
    double HalfWidth = 1.1,
    double HalfThickness = 0.12,
    double NoiseSigma = 0.06,
    int TwistCount = 1,
    double RadialBias = 1.0,
    string CrossSection = "Annular",
    string SpineShape = "FigureEight",
    double SplayFactor = 0.7,
    // ── HyperbolicBlobs / HyperbolicBlattHierarchy / Simplex / GaussianManifold shared knobs ────────────
    int ClusterCount = 3,
    int PointsPerCluster = 200,
    double Separation = 3.0,
    double Spread = 0.5,
    // ── HyperbolicBlattHierarchy-specific ────────────────────────────────────
    int HierarchyPoints = 1200,
    int HierarchyDepth = 3,
    int BranchesPerNode = 3,
    int BasePointsPerLeaf = 50,
    double RadiusDecay = 0.65,
    // ── Simplex-specific ───────────────────────────────────────────────────────
    int Categories = 12,
    bool DisjointSupports = true,
    double Concentration = 40.0,
    // ── GaussianManifold-specific ─────────────────────────────────────────────
    double ClusterRadius = 0.3,
    // ── TwoMoons ───────────────────────────────────────────────────────────────
    int PointsPerMoon = 100,
    double Noise = 0.1,
    // ── BlattHierarchy ─────────────────────────────────────────────────────────
    int CoarseClusters = 2,
    int MediumPerCoarse = 3,
    int FinePerMedium = 4,
    int PointsPerFine = 25,
    double CoarseSeparation = 20.0,
    double MediumSeparation = 4.0,
    double FineSeparation = 0.8,
    double LeafSpread = 0.15,
    // ── BlattThreeCluster ─────────────────────────────────────────────────────
    double StdDev = 1.0,
    // ── Shared ellipsoid ───────────────────────────────────────────────────────
    int EllipsoidPoints = 5000,
    double[]? EllipsoidAxes = null,
    string EllipsoidShellMode = "Solid",
    // ── Embedding (Möbius only) ────────────────────────────────────────────
    int Dimensions = 3,
    // ── Placement (shared shape) ───────────────────────────────────────────────
    string Placement = "OrthogonalElbowIntersect",
    double IntersectDepth = 0.0,
    double IntersectRadialShift = 0.0,
    double GapScale = 1.0,
    // ── Graph ──────────────────────────────────────────────────────────────────
    int Seed = 42,
    int KnnK = 5,
    string Metric = "euclidean",
    string NeighborRule = "Knn",
    double EpsilonBallEpsilon = 0.5,
    string Kernel = "Gaussian",
    double Bandwidth = 0.0,    // 0 = auto-estimate via BandwidthEstimation
    bool EnsureConnected = false,  // MST modifier: bridge disconnected components
                                   // ── Overlay ────────────────────────────────────────────────────────────────
    string GmmMode = "Oracle", // "Oracle" or "EM"
    bool ShowFlow = false,
    string FlowMode = "Beam");
