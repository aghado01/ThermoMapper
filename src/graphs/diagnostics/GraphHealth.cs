using System.Collections.Generic;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Primitives;
using Graphs.Observables;

namespace Graphs.Diagnostics
{
    /// <summary>
    /// Verdict bits derived from a <see cref="GraphHealthReport"/> — the
    /// "did the graph come out OK?" boolean fan-out. Each flag has a
    /// documented detection rule; <see cref="PrimaryRecommendation"/> is
    /// the single actionable string the REPL surfaces by default.
    /// </summary>
    public readonly record struct GraphHealthVerdict(
        bool    BandwidthTooLarge,
        bool    BandwidthTooSmall,
        bool    HubnessConcern,
        bool    ConnectivityConcern,
        bool    ForcedEdgesConcern,
        bool    UnderconnectedNodes,
        string? PrimaryRecommendation);

    /// <summary>
    /// Composite diagnostic over the rich
    /// <see cref="GraphBuildResult"/> produced by
    /// <see cref="GraphCompiler.Build"/>. Runs the affordable diagnostic
    /// subset, packages each report alongside its peers, and distils
    /// the lot into a <see cref="GraphHealthVerdict"/> with one
    /// actionable recommendation. The expensive
    /// <see cref="AlgebraicConnectivity"/> pass is included only when
    /// the graph is small enough to compute it dense; for larger graphs
    /// the cheaper <see cref="Connectivity"/> result drives the verdict.
    /// </summary>
    public sealed record GraphHealthReport(
        ConnectivityReport           Connectivity,
        DegreeReport                 Degree,
        EdgeWeightSummary            EdgeWeights,
        HubnessReport                Hubness,
        NeighborhoodScaleReport      NeighborhoodScale,
        CycleReport                  Cycles,
        AlgebraicConnectivityReport? AlgebraicConnectivity,
        MstBridgeReport?             MstBridge,
        // ── Graph initialization fingerprint (PS.1) ──────────────────────
        // Captured alongside the diagnostic battery so the persisted
        // report carries enough provenance to answer "what bandwidth /
        // metric was actually used?" without re-deriving them.
        double?                      ResolvedBandwidth,
        MixtureBandwidth?            ResolvedMixtureBandwidth,
        MetricProperties?            Metric,
        GraphHealthVerdict           Verdict);

    public static class GraphHealth
    {
        /// <summary>
        /// Run the standard diagnostic battery on a fresh
        /// <see cref="GraphBuildResult"/>. Returns a packaged report plus
        /// a single-line recommendation suitable for REPL display.
        /// </summary>
        /// <param name="build">Rich-shape build result.</param>
        /// <param name="k">Effective k used during construction (drives
        /// Hubness expectations + NeighborhoodScale window).</param>
        /// <param name="maxNodesForAlgebraic">Skip the dense
        /// <see cref="AlgebraicConnectivity"/> pass above this size and
        /// fall back to <see cref="Connectivity"/>-based signals.</param>
        public static GraphHealthReport Evaluate(
            GraphBuildResult build, int k, int maxNodesForAlgebraic = 2000)
        {
            var graph        = build.Graph;
            var connectivity = Connectivity.Validate(graph);
            var degree       = Degree.Distribution(graph);
            var edgeWeights  = EdgeWeights.Summary(graph);
            var hubness      = Hubness.Analyze(build.DirectedSelection, k);
            var scale        = NeighborhoodScale.Compute(build.DirectedSelection, k);

            // Cycles consume the already-computed Connectivity for the
            // cyclomatic-complexity arithmetic — no extra component pass.
            // The triangle / girth statistics skip themselves for graphs
            // above CycleReport's internal threshold (default 5000 nodes).
            var cycles = Cycles.Compute(graph, connectivity);

            AlgebraicConnectivityReport? algebraic = graph.NodeCount <= maxNodesForAlgebraic
                ? AlgebraicConnectivity.Compute(graph, maxNodesForDense: maxNodesForAlgebraic)
                : null;

            // GraphCompiler stashes the pre-repair graph in
            // PreRepairGraph when ensureConnected was used; if present,
            // compare the two so MstBridge can surface bridge-edge weight
            // statistics. Skipped otherwise (no comparison possible).
            MstBridgeReport? mstBridge = build.PreRepairGraph is CsrGraph preRepair
                ? MstBridge.Compare(preRepair, graph)
                : null;

            var verdict = DeriveVerdict(
                connectivity, degree, edgeWeights, hubness, build.EnsureConnectedApplied);

            return new GraphHealthReport(
                Connectivity:             connectivity,
                Degree:                   degree,
                EdgeWeights:              edgeWeights,
                Hubness:                  hubness,
                NeighborhoodScale:        scale,
                Cycles:                   cycles,
                AlgebraicConnectivity:    algebraic,
                MstBridge:                mstBridge,
                ResolvedBandwidth:        build.SingleBandwidth,
                ResolvedMixtureBandwidth: build.MixtureBandwidth,
                Metric:                   build.Metric,
                Verdict:                  verdict);
        }

        private static GraphHealthVerdict DeriveVerdict(
            ConnectivityReport connectivity,
            DegreeReport       degree,
            EdgeWeightSummary  edgeWeights,
            HubnessReport      hubness,
            bool               ensureConnectedApplied)
        {
            // Thresholds picked from the perplexity-notes design discussion;
            // tuned for the proximity-graph regime SPC operates in (~50-5000
            // nodes, bounded k, kernel-transformed weights in (0, 1]).
            const double BandwidthLowMedian   = 0.05;   // collapsed weights → bandwidth too large
            const double BandwidthHighMedian  = 0.95;   // saturated weights → bandwidth too small
            const double HubnessSkewLimit     = 3.0;    // in-degree skewness above this is concerning
            const double LargestCoverageLimit = 0.95;   // less than this → meaningful disconnection
            const double UndersampledShare    = 0.10;   // fraction of nodes at degree==1

            bool bandwidthTooLarge =
                edgeWeights.EdgeCount > 0 && edgeWeights.MedianWeight < BandwidthLowMedian;
            bool bandwidthTooSmall =
                edgeWeights.EdgeCount > 0 && edgeWeights.MedianWeight > BandwidthHighMedian;
            bool hubnessConcern    = hubness.InDegreeSkewness > HubnessSkewLimit;
            bool connectivityConcern =
                connectivity.LargestComponent < (int)(LargestCoverageLimit * connectivity.NodeCount);
            // Forced-edges concern is meaningful only when MST repair was
            // actually requested — otherwise there are no "forced" edges
            // to be concerned about.
            bool forcedEdgesConcern   = ensureConnectedApplied && connectivityConcern;
            bool underconnectedNodes  =
                degree.NodeCount > 0
                && degree.UndersampledCount > (int)(UndersampledShare * degree.NodeCount);

            string? recommendation = PickPrimaryRecommendation(
                bandwidthTooLarge, bandwidthTooSmall, hubnessConcern,
                connectivityConcern, forcedEdgesConcern, underconnectedNodes,
                edgeWeights, hubness, connectivity, degree);

            return new GraphHealthVerdict(
                BandwidthTooLarge:     bandwidthTooLarge,
                BandwidthTooSmall:     bandwidthTooSmall,
                HubnessConcern:        hubnessConcern,
                ConnectivityConcern:   connectivityConcern,
                ForcedEdgesConcern:    forcedEdgesConcern,
                UnderconnectedNodes:   underconnectedNodes,
                PrimaryRecommendation: recommendation);
        }

        /// <summary>
        /// Picks the single most-actionable recommendation by walking the
        /// signal flags in priority order: bandwidth issues first (most
        /// likely root cause, easiest to fix), then connectivity, then
        /// hubness, then degree distribution.
        /// </summary>
        private static string? PickPrimaryRecommendation(
            bool bandwidthTooLarge, bool bandwidthTooSmall, bool hubnessConcern,
            bool connectivityConcern, bool forcedEdgesConcern, bool underconnectedNodes,
            EdgeWeightSummary edgeWeights, HubnessReport hubness,
            ConnectivityReport connectivity, DegreeReport degree)
        {
            if (bandwidthTooLarge)
                return $"Median edge weight {edgeWeights.MedianWeight:G3} is below 0.05 — " +
                       "kernel bandwidth is likely too large (weights collapsed). Try halving --bandwidth, " +
                       "or omit it to let MAD auto-estimate from NN distances.";

            if (bandwidthTooSmall)
                return $"Median edge weight {edgeWeights.MedianWeight:G3} is above 0.95 — " +
                       "kernel bandwidth is likely too small (weights saturated). Try doubling --bandwidth.";

            if (forcedEdgesConcern)
                return $"--ensure-connected papered over a real disconnection " +
                       $"({connectivity.ComponentCount} components before repair). " +
                       "Consider raising --k, switching --proximity from mutualknn to knn, or relaxing the kernel.";

            if (connectivityConcern)
                return $"Largest component covers only {connectivity.LargestComponent}/{connectivity.NodeCount} nodes. " +
                       "The graph is disconnected — add --ensure-connected, or raise --k.";

            if (hubnessConcern)
                return $"In-degree skewness {hubness.InDegreeSkewness:F2} is high. " +
                       "Hub concentration may indicate the curse of dimensionality; " +
                       "consider switching --proximity to mutualknn or pre-projecting to lower dimensions.";

            if (underconnectedNodes)
                return $"{degree.UndersampledCount}/{degree.NodeCount} nodes have degree 1. " +
                       "Local structure is sparse; consider raising --k.";

            return null;
        }
    }
}
