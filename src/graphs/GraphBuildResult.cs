using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;
using Graphs.Primitives;

namespace Graphs
{
    /// <summary>
    /// Rich-shape return for <see cref="GraphCompiler.Build"/>. Carries
    /// the assembled <see cref="CsrGraph"/> plus the intermediate artifacts
    /// that downstream diagnostics need: the directed (pre-symmetrization)
    /// <see cref="NeighborSelection"/>, the resolved bandwidth, and the
    /// nearest-neighbor distance sample that bandwidth was derived from.
    /// </summary>
    /// <remarks>
    /// <para><b>Why this exists.</b> The bare <c>CsrGraph</c> return loses
    /// information that diagnostics need — <c>Hubness.Analyze</c> wants the
    /// directed selection (counting in-degree pre-symmetrization),
    /// <c>NeighborhoodScale</c> wants raw NN distances, and a bandwidth
    /// retroactive check wants to know what bandwidth was actually used. The
    /// historical <c>Build</c> path discarded all of this and forced
    /// consumers to re-run KNN to recover it.</para>
    ///
    /// <para><b>Bandwidth discriminator.</b> Exactly one of
    /// <see cref="SingleBandwidth"/> and <see cref="MixtureBandwidth"/> is
    /// non-null: single-kernel paths populate the former, mixture paths
    /// populate the latter. Consumers that need a representative scale
    /// without caring which side of the union it came from can read
    /// <see cref="RepresentativeBandwidth"/>.</para>
    ///
    /// <para><b>Metric provenance.</b> <see cref="Metric"/> is non-null
    /// when the caller supplied a custom <see cref="IDistanceMetric"/> —
    /// it carries the metric's name, bounded/domain properties, and the
    /// bandwidth strategy that was selected. Null when the inline default
    /// Euclidean path was used. Graph-health persistence records this so
    /// a reviewer can confirm "Poincaré → LogScaleHyperbolic" without
    /// re-reading the manifest's metric spec string.</para>
    /// </remarks>
    public sealed record GraphBuildResult(
        CsrGraph              Graph,
        NeighborSelection     DirectedSelection,
        double?               SingleBandwidth,
        MixtureBandwidth?     MixtureBandwidth,
        bool                  EnsureConnectedApplied,
        CsrGraph?             PreRepairGraph = null,
        MetricProperties?     Metric         = null,
        DiagnosticsLog?       Diagnostics    = null,
        EdgeWeightKind        WeightKind     = EdgeWeightKind.Coupling)
    {
        /// <summary>
        /// Convenience accessor: returns <see cref="SingleBandwidth"/> for
        /// single-kernel builds, or the Gaussian component of the mixture
        /// bandwidth otherwise. Diagnostic code that just wants "a number
        /// in distance units" can use this without branching on which build
        /// variant produced the result.
        /// </summary>
        public double RepresentativeBandwidth =>
            SingleBandwidth ?? MixtureBandwidth?.Gaussian ?? 0.0;
    }
}
