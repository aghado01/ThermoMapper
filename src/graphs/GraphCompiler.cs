using System;
using System.Collections.Generic;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;
using Graphs.Pipeline;
using Graphs.Pipeline.Filters;
using Graphs.Pipeline.Generators;
using Graphs.Pipeline.Refinement;
using Graphs.Pipeline.Repair;
using Graphs.Pipeline.Scalers;
using Graphs.Primitives;
using Maths.Geometry;

namespace Graphs;

/// <summary>
/// Caller-supplied source of protected edges for demoting passes (today: LMP
/// rescaling). Receives the conditioned distance graph of the selection being
/// scaled; returns normalized <c>(Lo, Hi)</c> vertex pairs to shield, or
/// <c>null</c> for none. Keeps construction consumer-neutral: topological
/// criteria (e.g. involuted-persistence H1 cycle edges) are computed above
/// this seam and injected — never imported by a build stage.
/// </summary>
public delegate IReadOnlySet<(int Lo, int Hi)>? ProtectedEdgeSource(CsrGraph distanceGraph);

/// <summary>
/// Pure declarative graph-construction engine. Takes a single immutable
/// <see cref="GraphCompilerConfig"/> + the data plumbing it needs;
/// returns a <see cref="GraphBuildResult"/> carrying the assembled
/// <see cref="CsrGraph"/> and the full forensic trail.
/// </summary>
/// <remarks>
/// <para>No instance state, no fluent chain. The engine evaluates the
/// whole config holistically before executing — order-dependent
/// concerns (LMP interrupts global bandwidth, post-repair-timing
/// only meaningful when repair ran, etc.) are resolved internally.</para>
///
/// <para>Fluent / chained APIs that translate user intent into a
/// <see cref="GraphCompilerConfig"/> live at the CLI / REPL boundary
/// (see <c>UserRepl.Commands.GraphCommandBuilder</c> — Task #14) and
/// are not exposed by the backend.</para>
///
/// <para><b>Two entry points:</b> the static <see cref="Build"/> for
/// one-shot construction, and the <see cref="GraphSession"/> class
/// for the build-many-against-the-same-data sweep pattern.</para>
/// </remarks>
public static class GraphCompiler
{
    /// <summary>
    /// Build a graph from a config + data substrate. Returns the
    /// assembled graph plus the diagnostic trail. Throws
    /// <see cref="GraphPathologyException"/> if any configured
    /// pathology threshold is crossed; throws
    /// <see cref="ArgumentException"/> on a structurally invalid
    /// config.
    /// </summary>
    public static GraphBuildResult Build(
        GraphCompilerConfig config,
        int n,
        GraphMetric metric,
        ProtectedEdgeSource? protectedEdges = null)
    {
        if (config is null) throw new ArgumentNullException(nameof(config));
        if (metric is null) throw new ArgumentNullException(nameof(metric));
        if (metric.Distance is null) throw new ArgumentNullException(nameof(metric));
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

        const double HubnessWarningThreshold = 2.0;
        const double CoverageWarningThreshold = 0.98;
        const double MedianWeightWarningBand = 0.20; // 20% above/below configured limits
        const double DefaultAutoLmpNearZeroRatio = 0.05;

        var log = new DiagnosticsLog();
        Func<int, int, double> distance = metric.Distance;

        // Stage 1 — Topology (required, no auto)
        ITopologyGenerator generator = ResolveTopology(config.Topology, log);
        NeighborSelection directed = generator.Generate(n, distance);
        log.Info("Topology", $"Directed neighbor pass complete (n={n}).");

        // Diagnostic gate 1: hubness skewness
        // (Auto-pick for Filter stage when applicable; pathology check otherwise.)
        var hubness = Hubness.Analyze(directed, ResolveK(config.Topology));
        log.Info("Topology", $"InDegreeSkewness={hubness.InDegreeSkewness:F3}.");
        if (hubness.InDegreeSkewness > HubnessWarningThreshold)
        {
            log.Warning(
                "Topology",
                $"Hubness skewness is elevated ({hubness.InDegreeSkewness:F3}); " +
                "consider MutualKnn filtering or reviewing the topology choice.");
        }
        if (config.Interrupts.MaxHubnessSkewness is double skewLimit
            && hubness.InDegreeSkewness > skewLimit)
        {
            throw new GraphPathologyException(
                $"InDegreeSkewness {hubness.InDegreeSkewness:F3} exceeds configured " +
                $"MaxHubnessSkewness {skewLimit:F3}.");
        }

        // Stage 2 — Filter (Auto resolved here)
        ITopologyFilter filter = ResolveFilter(config.Filter, hubness, log);
        NeighborSelection filtered = filter.Filter(directed, n, distance);

        // Diagnostic gate 2: connectivity coverage
        var connectivityPre = Connectivity.Validate(BuildIndicatorGraph(filtered, n));
        double coveragePre = connectivityPre.NodeCount == 0
            ? 1.0
            : (double)connectivityPre.LargestComponent / connectivityPre.NodeCount;
        log.Info(
            "Filter",
            $"Largest component covers {connectivityPre.LargestComponent}/{connectivityPre.NodeCount} " +
            $"after filtering ({coveragePre:F3}).");

        if (coveragePre < CoverageWarningThreshold)
        {
            log.Warning(
                "Filter",
                $"Coverage after filtering is {coveragePre:P1}; repair may be advisable to avoid a fragmented graph.");
        }

        // Stage 3 — Repair (Auto resolved here)
        ITopologyRepair repair = ResolveRepair(config.Repair, connectivityPre, log);
        int edgesBeforeRepair = CountUndirectedEdges(filtered, n);
        NeighborSelection repaired = repair.Repair(filtered, n, distance);
        bool repairApplied = CountUndirectedEdges(repaired, n) > edgesBeforeRepair;
        RepairEdgeSet repairEdges = RepairEdgeSet.Between(filtered, repaired, n);

        if (config.Repair.Kind == RepairKind.Auto && repairApplied)
        {
            log.Warning(
                "Repair",
                "MST-based repair was applied because filtered coverage was below the auto-repair threshold.");
        }

        // Re-check connectivity post-repair (NoRepair leaves it the same;
        // MstMinRepair should fix it). The returned graph is compared with this.
        var connectivityPost = Connectivity.Validate(BuildIndicatorGraph(repaired, n));

        if (config.Interrupts.MinLargestCoverage is double coverageLimit)
        {
            double coverage = connectivityPost.NodeCount == 0
                ? 1.0
                : (double)connectivityPost.LargestComponent / connectivityPost.NodeCount;
            if (coverage < coverageLimit)
            {
                throw new GraphPathologyException(
                    $"Largest-component coverage {coverage:F3} after Stage 3 falls below " +
                    $"configured MinLargestCoverage {coverageLimit:F3}.");
            }
            else if (coverage < CoverageWarningThreshold)
            {
                log.Warning(
                    "Repair",
                    $"Coverage after repair is {coverage:P1}; graph is formally connected but still sparse.");
            }
        }

        // Stage 4 — Refinement (Euclidean is the only impl until PathNeighbor lands)
        IMetricRefiner refiner = ResolveRefinement(config.Refinement, log);
        NeighborSelection refined = refiner.Refine(repaired, n);

        // Stage 5 — Projection (L3): emit exactly the weight kind requested.
        if (config.Projection is DistanceProjection)
        {
            log.Info("Projection", "Distance pass-through — conditioned distance graph (L2).");
            CsrGraph distanceGraph = BuildDistanceGraph(refined, n);
            LogReturnedGraph(log, distanceGraph, n, repairEdges, connectivityPost, cutEdges: 0, couplingFloor: 0.0, weightsAreCouplings: false);
            return new GraphBuildResult(
                Graph:                  distanceGraph,
                DirectedSelection:      directed,
                SingleBandwidth:        null,
                MixtureBandwidth:       null,
                EnsureConnectedApplied: repairApplied,
                PreRepairGraph:         repairApplied ? WithoutEdges(distanceGraph, repairEdges, n) : null,
                Metric:                 metric.Properties,
                Diagnostics:            log,
                WeightKind:             EdgeWeightKind.Distance);
        }

        if (config.Projection is not CouplingProjection coupling)
            throw new NotSupportedException(
                $"Projection '{config.Projection.GetType().Name}' is not yet wired in " +
                "GraphCompiler.Build. Supported: DistanceProjection, CouplingProjection.");

        // Coupling projection — kernel → coupling weights (the former Stage 5).
        bool? lmpSetting = coupling.LmpRescale;
        log.Info("Scaling", lmpSetting switch
        {
            true  => "LMP rescaling enabled (explicit).",
            false => "LMP rescaling disabled (explicit).",
            null  => "LMP rescaling auto-decided from tentative edge statistics."
        });

        IEdgeScaler scaler = BuildScaler(
            coupling,
            repairEdges,
            protectedEdges,
            lmpSetting,
            metric.Properties,
            metric.AmbientDimension,
            metric.Features,
            metric.Manifold,
            log,
            DefaultAutoLmpNearZeroRatio);
        ScalerResult scaled = scaler.Scale(refined, n);
        LogReturnedGraph(log, scaled.Graph, n, repairEdges, connectivityPost, scaled.CutEdges, coupling.CouplingFloor, weightsAreCouplings: true);

        if (scaled.SingleBandwidth is double delta && delta > 0.0)
        {
            double sampleMedian = ComputeMedian(refined.KthNeighborDistances);
            if (sampleMedian > 0.0)
            {
                if (delta < sampleMedian * 0.1)
                {
                    log.Warning(
                        "Scaling",
                        $"Estimated bandwidth {delta:G3} is much smaller than median k-NN " +
                        $"distance {sampleMedian:G3}; edge weights may decay too quickly.");
                }
                else if (delta > sampleMedian * 10.0)
                {
                    log.Warning(
                        "Scaling",
                        $"Estimated bandwidth {delta:G3} is much larger than median k-NN " +
                        $"distance {sampleMedian:G3}; edge weights may be overly flat.");
                }
            }
        }

        // Diagnostic gate 4: edge-weight summary (median / near-zero ratio)
        if (config.Interrupts.MinMedianWeight is double minMed
            || config.Interrupts.MaxMedianWeight is double maxMed
            || config.Interrupts.MaxNearZeroEdgeRatio is double maxNZ)
        {
            var ew = EdgeWeights.Summary(scaled.Graph);
            log.Info("Scaling", $"MedianWeight={ew.MedianWeight:G3} over {ew.EdgeCount} edges.");

            if (config.Interrupts.MinMedianWeight is double minM)
            {
                if (ew.MedianWeight < minM)
                {
                    throw new GraphPathologyException(
                        $"MedianWeight {ew.MedianWeight:G3} below configured MinMedianWeight {minM:G3} " +
                        "— bandwidth likely too large.");
                }
                else if (ew.MedianWeight < minM * (1.0 + MedianWeightWarningBand))
                {
                    log.Warning(
                        "Scaling",
                        $"Median edge weight {ew.MedianWeight:G3} is near the configured minimum {minM:G3}; " +
                        "consider increasing bandwidth or inspecting edge decay.");
                }
            }
            if (config.Interrupts.MaxMedianWeight is double maxM)
            {
                if (ew.MedianWeight > maxM)
                {
                    throw new GraphPathologyException(
                        $"MedianWeight {ew.MedianWeight:G3} above configured MaxMedianWeight {maxM:G3} " +
                        "— bandwidth likely too small.");
                }
                else if (ew.MedianWeight > maxM * (1.0 - MedianWeightWarningBand))
                {
                    log.Warning(
                        "Scaling",
                        $"Median edge weight {ew.MedianWeight:G3} is near the configured maximum {maxM:G3}; " +
                        "consider decreasing bandwidth or reviewing the kernel choice.");
                }
            }
            if (config.Interrupts.MaxNearZeroEdgeRatio is double maxRatio)
            {
                double nzRatio = EdgeWeightThresholds.NearZeroRatio(scaled.Graph, scaled.CutEdges);
                log.Info("Scaling", $"NearZeroEdgeRatio={nzRatio:F4} (threshold={maxRatio:F4}).");
                if (nzRatio > maxRatio)
                {
                    throw new GraphPathologyException(
                        $"NearZeroEdgeRatio {nzRatio:F4} exceeds configured " +
                        $"MaxNearZeroEdgeRatio {maxRatio:F4} — delta-collapse pathology.");
                }
            }
        }

        // The pre-repair view is the returned graph without the repair edges: the same units and
        // weights, so MstBridge.Compare sees the topology change alone.
        CsrGraph? preRepairGraph = repairApplied ? WithoutEdges(scaled.Graph, repairEdges, n) : null;

        return new GraphBuildResult(
            Graph:                  scaled.Graph,
            DirectedSelection:      directed,
            SingleBandwidth:        scaled.SingleBandwidth,
            MixtureBandwidth:       scaled.MixtureBandwidth,
            EnsureConnectedApplied: repairApplied,
            PreRepairGraph:         preRepairGraph,
            Metric:                 metric.Properties,
            Diagnostics:            log,
            WeightKind:             EdgeWeightKind.Coupling);
    }

    // ── Stage resolution ─────────────────────────────────────────────────

    private static ITopologyGenerator ResolveTopology(TopologyConfig cfg, DiagnosticsLog log)
    {
        switch (cfg.Kind)
        {
            case TopologyKind.Knn:
                int k = cfg.K ?? 10;
                log.Info("Topology", $"Knn (k={k}).");
                return new KnnGenerator(k);
            case TopologyKind.EpsilonBall:
                if (cfg.Epsilon is not double eps || eps <= 0.0)
                    throw new ArgumentException(
                        "TopologyConfig.Epsilon is required and must be positive when Kind == EpsilonBall.");
                log.Info("Topology", $"EpsilonBall (ε={eps:G3}).");
                return new EpsilonBallGenerator(eps);
            default:
                throw new ArgumentOutOfRangeException(nameof(cfg.Kind), cfg.Kind, "Unknown TopologyKind.");
        }
    }

    private static int ResolveK(TopologyConfig cfg) =>
        cfg.Kind == TopologyKind.Knn ? (cfg.K ?? 10) : 0;

    private static ITopologyFilter ResolveFilter(
        FilterConfig cfg, HubnessReport hubness, DiagnosticsLog log)
    {
        FilterKind kind = cfg.Kind;
        if (kind == FilterKind.Auto)
        {
            const double SkewnessAutoThreshold = 3.0;
            kind = hubness.InDegreeSkewness > SkewnessAutoThreshold
                ? FilterKind.MutualKnn
                : FilterKind.OrRule;
            if (kind == FilterKind.MutualKnn)
            {
                log.Warning(
                    "Filter",
                    $"Auto-picked MutualKnn because skewness={hubness.InDegreeSkewness:F3} " +
                    "exceeded the safe threshold (3.0). This is a conservative topology choice.");
            }
            else
            {
                log.Info(
                    "Filter",
                    $"Auto: skewness={hubness.InDegreeSkewness:F3} " +
                    "→ OrRule (≤ 3.0)." );
            }
        }
        else
        {
            log.Info("Filter", $"Explicit: {kind}.");
        }

        return kind switch
        {
            FilterKind.OrRule    => new PassThroughFilter(),
            FilterKind.MutualKnn => new MutualKnnFilter(
                cfg.MutualBandwidthSource ?? Graphs.MutualBandwidthSource.DirectedKth),
            _ => throw new ArgumentOutOfRangeException(nameof(kind), kind, "Unresolvable FilterKind."),
        };
    }

    private static ITopologyRepair ResolveRepair(
        RepairConfig cfg, ConnectivityReport connectivity, DiagnosticsLog log)
    {
        RepairKind kind = cfg.Kind;
        if (kind == RepairKind.Auto)
        {
            const double CoverageAutoThreshold = 0.95;
            double coverage = connectivity.NodeCount == 0
                ? 1.0
                : (double)connectivity.LargestComponent / connectivity.NodeCount;
            kind = coverage < CoverageAutoThreshold ? RepairKind.MstMin : RepairKind.NoRepair;
            if (kind == RepairKind.MstMin)
            {
                log.Warning(
                    "Repair",
                    $"Auto-picked MstMin because coverage={coverage:F3} fell below the safe threshold (0.95)." );
            }
            else
            {
                log.Info(
                    "Repair",
                    $"Auto: coverage={coverage:F3} → NoRepair (≥ 0.95)." );
            }
        }
        else
        {
            log.Info("Repair", $"Explicit: {kind}.");
        }

        return kind switch
        {
            RepairKind.NoRepair => new NoRepair(),
            RepairKind.MstMin   => new MstMinRepair(),
            _ => throw new ArgumentOutOfRangeException(nameof(kind), kind, "Unresolvable RepairKind."),
        };
    }

    private static IMetricRefiner ResolveRefinement(RefinementConfig cfg, DiagnosticsLog log)
    {
        if (cfg.Kind == RefinementKind.Auto || cfg.Kind == RefinementKind.PathNeighbor)
        {
            log.Info("Refinement", cfg.Kind == RefinementKind.Auto
                ? "Auto → PathNeighbor."
                : $"Explicit: {cfg.Kind}.");

            return new PathNeighborRefiner(cfg.MaxDistance);
        }

        throw new NotSupportedException(
            $"Unsupported refinement kind {cfg.Kind}. Euclidean refinement is retired.");
    }

    private static IEdgeScaler BuildScaler(
        CouplingProjection coupling,
        RepairEdgeSet repairEdges,
        ProtectedEdgeSource? protectedEdges,
        bool? lmpSetting,
        MetricProperties? metricProps,
        int? ambientDimension,
        double[][]? features,
        IRiemannianManifold? manifold,
        DiagnosticsLog log,
        double defaultAutoLmpThreshold)
    {
        if (double.IsNaN(coupling.CouplingFloor) || coupling.CouplingFloor < 0.0)
            throw new ArgumentException("CouplingProjection.CouplingFloor must be zero or positive.");
        var context = new ScalingContext(repairEdges, coupling.CouplingFloor);

        IKernelDescriptor kernel = coupling.Kernel;
        BandwidthStrategy strategy =
            coupling.BandwidthOverride
            ?? metricProps?.BandwidthStrategy
            ?? BandwidthStrategy.MadConsistencyFactor;
        SpaceGeometry geometry = metricProps?.Geometry ?? SpaceGeometry.Euclidean;
        CouplingFidelity fidelity = coupling.Fidelity switch
        {
            // Spherical's geometrically-exact path is Intrinsic (the heat-kernel
            // parametrix). GeodesicLinear on the sphere is only a first-order
            // fallback and warns below, so Auto defaults to Intrinsic there;
            // every other geometry keeps the GeodesicLinear default.
            CouplingFidelity.Auto => geometry == SpaceGeometry.Spherical
                ? CouplingFidelity.Intrinsic
                : CouplingFidelity.GeodesicLinear,
            var explicitChoice => explicitChoice,
        };

        SphericalIntrinsicMode sphericalMode = coupling.SphericalMode switch
        {
            SphericalIntrinsicMode.Auto => SphericalIntrinsicMode.LocalParametrix,
            var explicitChoice => explicitChoice,
        };

        if (fidelity == CouplingFidelity.Intrinsic
            && geometry == SpaceGeometry.Hyperbolic
            && ambientDimension is null)
        {
            throw new InvalidOperationException(
                "Intrinsic coupling on hyperbolic geometry requires the ambient dimension; build via GraphMetric.FromFeatures.");
        }

        if (coupling.BandwidthOverride is null && metricProps is null)
        {
            log.Warning(
                "Scaling",
                "Kernel projection has no declared metric provenance; bandwidth strategy fell back to Euclidean MAD.");
        }

        if (metricProps is { } mp && NaturalStrategyFor(mp.Geometry) != mp.BandwidthStrategy)
        {
            log.Warning(
                "Scaling",
                $"Metric '{mp.Name}' declares geometry {mp.Geometry} but strategy {mp.BandwidthStrategy}; expected {NaturalStrategyFor(mp.Geometry)}.");
        }

        if (fidelity == CouplingFidelity.GeodesicLinear
            && (geometry == SpaceGeometry.Hyperbolic || geometry == SpaceGeometry.Spherical))
        {
            log.Warning(
                "Scaling",
                $"GeodesicLinear coupling on {geometry} geometry — first-order pathway; Intrinsic is the geometrically-exact alternative.");
        }

        log.Info("Scaling",
            $"Kernel={kernel.GetType().Name}, Strategy={strategy}, Geometry={geometry}, Fidelity={fidelity}, SphericalMode={sphericalMode}, CouplingFloor={coupling.CouplingFloor:G}, Gauge={coupling.Gauge}.");

        // LMP leaves the injected protect set (H1 cycle edges) and the repair bridges their kernel coupling.
        Func<NeighborSelection, int, IReadOnlySet<(int Lo, int Hi)>?>? protectedProvider = null;
        if (protectedEdges is not null || repairEdges.Count > 0)
        {
            protectedProvider = (sel, nodeCount) =>
            {
                IReadOnlySet<(int Lo, int Hi)>? injected =
                    protectedEdges?.Invoke(BuildDistanceGraph(sel, nodeCount));
                if (repairEdges.Count == 0) return injected;
                if (injected is null) return repairEdges.Pairs;
                var union = new HashSet<(int Lo, int Hi)>(injected);
                union.UnionWith(repairEdges.Pairs);
                return union;
            };
            if (protectedEdges is not null)
                log.Info("Scaling", "Protected-edge source supplied; LMP will preserve the designated edges.");
            if (repairEdges.Count > 0)
                log.Info("Scaling", $"LMP will keep the {repairEdges.Count} repair edge(s) at their kernel coupling.");
        }

        IEdgeScaler inner = kernel switch
        {
            Gaussian g  => new GlobalBandwidthScaler(KernelType.Gaussian,  g.Bandwidth,   strategy, geometry, fidelity, ambientDimension, features, manifold, sphericalMode, context),
            Cauchy c    => new GlobalBandwidthScaler(KernelType.Cauchy,    c.Bandwidth,   strategy, geometry, fidelity, ambientDimension, features, manifold, sphericalMode, context),
            Laplacian l => new GlobalBandwidthScaler(KernelType.Laplacian, l.Bandwidth,   strategy, geometry, fidelity, ambientDimension, features, manifold, sphericalMode, context),
            Linear lin  => new GlobalBandwidthScaler(KernelType.Linear,    lin.Bandwidth, strategy, geometry, fidelity, ambientDimension, features, manifold, sphericalMode, context),
            Mixture m   => BuildMixtureScaler(m, strategy, geometry, fidelity, ambientDimension, features, manifold, sphericalMode, context),
            _ => throw new NotSupportedException($"Unknown kernel descriptor: {kernel.GetType().Name}"),
        };

        IEdgeScaler measured = lmpSetting switch
        {
            true => new LocalMutualProximityScaler(inner, protectedEdgeProvider: protectedProvider, rankRule: coupling.LmpRank),
            false => inner,
            null => new AutoLmpScaler(inner, log, defaultAutoLmpThreshold, protectedProvider, repairEdges, coupling.LmpRank),
        };

        // The gauge goes last, after the floor and LMP: a pure rescaling of temperature.
        return coupling.Gauge == GaugeKind.None ? measured : new GaugeScaler(measured, coupling.Gauge);
    }

    private static GlobalBandwidthScaler BuildMixtureScaler(
        Mixture m,
        BandwidthStrategy strategy,
        SpaceGeometry geometry,
        CouplingFidelity fidelity,
        int? ambientDimension,
        double[][]? features,
        IRiemannianManifold? manifold,
        SphericalIntrinsicMode sphericalMode,
        ScalingContext context)
    {
        var (weights, bandwidth) = m.ToLegacy();
        return new GlobalBandwidthScaler(weights, bandwidth, strategy, geometry, fidelity, ambientDimension, features, manifold, sphericalMode, context);
    }

    private static BandwidthStrategy NaturalStrategyFor(SpaceGeometry geometry) => geometry switch
    {
        SpaceGeometry.Hyperbolic => BandwidthStrategy.LogScaleHyperbolic,
        SpaceGeometry.Spherical  => BandwidthStrategy.QuantileNormalized,
        _                        => BandwidthStrategy.MadConsistencyFactor,
    };

    // ── Helpers ──────────────────────────────────────────────────────────

    /// <summary>
    /// Build a minimal unweighted CSR view from a NeighborSelection for
    /// pre-scaling diagnostic passes (connectivity). Weight=1.0 — only
    /// the topology matters.
    /// </summary>
    private static CsrGraph BuildIndicatorGraph(NeighborSelection sel, int n)
    {
        var edges = new System.Collections.Generic.List<Edge>();
        for (int i = 0; i < n; i++)
        {
            foreach (var nb in sel.AllNeighbors[i])
            {
                if (i > nb.Index) continue;
                edges.Add(new Edge(i, nb.Index, 1.0));
            }
        }
        return CsrGraph.FromEdges(edges.ToArray(), n);
    }

    /// <summary>
    /// Build the L2 conditioned distance graph: a symmetric CSR whose weights
    /// are the refined edge <em>distances</em> (not couplings). Sibling of
    /// <see cref="BuildIndicatorGraph"/>; this is the <see cref="DistanceProjection"/>
    /// output and the distance-valued graph HDBSCAN-sparse / Mapper consume.
    /// </summary>
    private static CsrGraph BuildDistanceGraph(NeighborSelection sel, int n)
    {
        var edges = new System.Collections.Generic.List<Edge>();
        for (int i = 0; i < n; i++)
        {
            foreach (var nb in sel.AllNeighbors[i])
            {
                if (i > nb.Index) continue;
                edges.Add(new Edge(i, nb.Index, nb.Distance));
            }
        }
        return CsrGraph.FromEdges(edges.ToArray(), n);
    }

    /// <summary>
    /// Log what the build returned: the component count and largest-component
    /// coverage of the returned graph, over all its edges and over the edges that
    /// carry coupling; the edges the coupling floor dropped; and what became of the
    /// repair stage's bridges. Repair counts its bridges as connecting; a scaler
    /// that drops them, or Local Mutual Proximity writing them an exact zero, leaves
    /// a graph the repair stage never reported.
    /// </summary>
    private static void LogReturnedGraph(
        DiagnosticsLog log,
        CsrGraph returned,
        int n,
        RepairEdgeSet repairEdges,
        ConnectivityReport repairedConnectivity,
        int cutEdges,
        double couplingFloor,
        bool weightsAreCouplings)
    {
        ConnectivityReport report = Connectivity.Validate(returned);
        log.Info(
            "Connectivity",
            $"Returned graph: {report.ComponentCount} component(s); largest covers " +
            $"{report.LargestComponent}/{report.NodeCount} ({report.LargestCoverage:F3}); {report.EdgeCount} edges.");

        if (cutEdges > 0)
        {
            log.Info(
                "Scaling",
                $"{cutEdges} edge(s) fell under the {couplingFloor:G} coupling floor and were dropped (repair edges are exempt).");
        }

        if (weightsAreCouplings)
        {
            // An edge of weight zero is in the graph but joins nothing: the kernel can underflow
            // to it, and Local Mutual Proximity writes it for a node's weakest edge.
            ConnectivityReport coupled = Connectivity.Validate(returned, positiveWeightsOnly: true);
            int zeroEdges = report.EdgeCount - coupled.EdgeCount;
            if (zeroEdges > 0)
            {
                string text =
                    $"{zeroEdges} of {report.EdgeCount} edge(s) carry zero coupling; over the {coupled.EdgeCount} " +
                    $"coupled edges the graph has {coupled.ComponentCount} component(s), largest " +
                    $"{coupled.LargestComponent}/{coupled.NodeCount} ({coupled.LargestCoverage:F3}).";
                if (coupled.ComponentCount > report.ComponentCount) log.Warning("Connectivity", text);
                else log.Info("Connectivity", text);
            }
        }

        if (repairEdges.Count > 0)
        {
            Dictionary<long, double> returnedWeights = EdgeWeightsByKey(returned, n);
            int missing = 0, zeroed = 0;
            foreach (long key in repairEdges.Keys)
            {
                if (!returnedWeights.TryGetValue(key, out double weight)) missing++;
                else if (weightsAreCouplings && !(weight > 0.0)) zeroed++;
            }

            if (missing > 0)
            {
                log.Warning(
                    "Connectivity",
                    $"{missing} of {repairEdges.Count} repair edge(s) are absent from the returned graph: " +
                    "the repair stage counted them as bridges.");
            }
            if (zeroed > 0)
            {
                log.Warning(
                    "Connectivity",
                    $"{zeroed} of {repairEdges.Count} repair edge(s) carry zero coupling in the returned graph: " +
                    "the repair stage counted them as bridges.");
            }
        }

        if (report.ComponentCount > repairedConnectivity.ComponentCount)
        {
            log.Warning(
                "Connectivity",
                $"The repaired topology had {repairedConnectivity.ComponentCount} component(s); the returned graph has " +
                $"{report.ComponentCount}. The graph is less connected than the repair stage reported.");
        }
    }

    /// <summary>A copy of <paramref name="graph"/> without the edges in <paramref name="removed"/>; weights unchanged.</summary>
    private static CsrGraph WithoutEdges(CsrGraph graph, RepairEdgeSet removed, int n)
    {
        var edges = new List<Edge>(graph.Targets.Length / 2);
        for (int i = 0; i < n; i++)
            for (int slot = graph.RowPointers[i]; slot < graph.RowPointers[i + 1]; slot++)
            {
                int j = graph.Targets[slot];
                if (j > i && !removed.Contains(i, j)) edges.Add(new Edge(i, j, graph.Weights[slot]));
            }
        return CsrGraph.FromEdges(edges.ToArray(), n);
    }

    private static Dictionary<long, double> EdgeWeightsByKey(CsrGraph graph, int n)
    {
        var weights = new Dictionary<long, double>(graph.Targets.Length / 2);
        for (int i = 0; i < n; i++)
            for (int slot = graph.RowPointers[i]; slot < graph.RowPointers[i + 1]; slot++)
                if (graph.Targets[slot] > i) weights[RepairEdgeSet.KeyOf(i, graph.Targets[slot], n)] = graph.Weights[slot];
        return weights;
    }

    private static int CountUndirectedEdges(NeighborSelection selection, int n)
    {
        int count = 0;
        for (int i = 0; i < n; i++)
        {
            foreach (var nb in selection.AllNeighbors[i])
            {
                if (i < nb.Index)
                    count++;
            }
        }

        return count;
    }

    private static double ComputeMedian(double[] values)
    {
        if (values is null || values.Length == 0)
            return 0.0;

        double[] copy = new double[values.Length];
        Array.Copy(values, copy, values.Length);
        Array.Sort(copy);

        int mid = copy.Length / 2;
        return copy.Length % 2 == 0
            ? (copy[mid - 1] + copy[mid]) / 2.0
            : copy[mid];
    }
}

/// <summary>
/// Thrown by <see cref="GraphCompiler.Build"/> when a configured
/// <see cref="PathologyInterruptConfig"/> threshold is crossed. The
/// message identifies the offending metric and threshold.
/// </summary>
public sealed class GraphPathologyException : Exception
{
    public GraphPathologyException(string message) : base(message) { }
}

/// <summary>
/// Stateful **data substrate** holder for the build-many-against-the-
/// same-data sweep pattern. Holds (<c>n</c>, <c>metric</c>); each
/// <see cref="Build"/> call takes a fresh
/// <see cref="GraphCompilerConfig"/> and produces a result. Future:
/// memoize intermediate stages keyed by config sub-trees so a
/// parameter sweep reuses the directed-KNN computation across runs.
/// </summary>
/// <remarks>
/// <para>This is **not** a fluent builder. The session has no
/// configuration state — only data substrate. Each Build call is fully
/// declarative.</para>
/// </remarks>
public sealed class GraphSession
{
    private readonly int _n;
    private readonly GraphMetric _metric;

    public GraphSession(int n, GraphMetric metric)
    {
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
        _n = n;
        _metric = metric ?? throw new ArgumentNullException(nameof(metric));
    }

    /// <summary>Build a graph for the given config against this session's data.</summary>
    public GraphBuildResult Build(GraphCompilerConfig config, ProtectedEdgeSource? protectedEdges = null) =>
        GraphCompiler.Build(config, _n, _metric, protectedEdges);
}
