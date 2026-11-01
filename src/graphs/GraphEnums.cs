namespace Graphs;

/// <summary>Internal kernel-family identifier. The fluent
/// <see cref="Coupling.IKernelDescriptor"/> descriptors resolve down to this at the
/// <see cref="Pipeline.Scalers.GlobalBandwidthScaler"/> boundary — it is the live
/// currency between the compiler and the scaler, not a public-facing API.</summary>
public enum KernelType { Gaussian, Cauchy, Laplacian, Linear }

/// <summary>Bandwidth sample source when mutual-KNN filtering is active.</summary>
public enum MutualBandwidthSource { DirectedKth, MutualKth }

/// <summary>
/// Semantics of a graph's edge weights — what the numbers in
/// <see cref="Primitives.CsrGraph.Weights"/> mean. Echoed onto
/// <see cref="GraphBuildResult.WeightKind"/> from the requested
/// <see cref="IEdgeProjection"/> so consumers can assert they were handed the
/// weight kind they expect (Potts requires Coupling; HDBSCAN-sparse requires
/// Distance).
/// </summary>
public enum EdgeWeightKind { Unweighted, Distance, Affinity, Coupling }

/// <summary>
/// How faithfully the metric's curvature is honored when a kernel turns a
/// distance into a coupling. Orthogonal to kernel shape
/// (<see cref="Coupling.IKernelDescriptor"/>) and metric geometry
/// (<see cref="Distance.SpaceGeometry"/>).
/// </summary>
public enum CouplingFidelity
{
	/// <summary>Resolve from geometry at build time (the delegated default).</summary>
	Auto,

	/// <summary>
	/// Geodesic distance into a flat kernel shape (first-order; historical
	/// behavior). Scale-absolute, sparse-graph-safe, not a Mercer kernel.
	/// </summary>
	GeodesicLinear,

	/// <summary>
	/// The manifold's canonical kernel (heat kernel) — PD by construction.
	/// Implemented for hyperbolic geometry (Gaussian; Van Vleck correction
	/// (r/sinh r)^((d-1)/2), exact for H^3).
	/// </summary>
	Intrinsic,

}

/// <summary>
/// Execution mode for positive-curvature Intrinsic spherical kernels,
/// dispatched at builder level.
/// </summary>
public enum SphericalIntrinsicMode
{
	Auto,
	LocalParametrix,
	GlobalSchoenberg
}

/// <summary>
/// How Local Mutual Proximity ranks an edge among a node's incident edges.
/// </summary>
public enum LmpRankRule
{
    /// <summary>Counts the strictly weaker edges (the historical rule). A node's weakest edge ranks 0
    /// at that node, so LMP writes it exactly 0; a degree-1 node loses its only edge.</summary>
    Strict,

    /// <summary>Counts the edges at most as strong, the edge itself included: the empirical CDF.
    /// The weakest edge ranks 1/degree, so no edge is written 0.</summary>
    Inclusive,
}

/// <summary>
/// A global factor applied to the finished couplings, last in the pipeline after every threshold,
/// so it is exactly a rescaling of temperature: <c>J → gJ</c> is <c>T → gT</c>.
/// </summary>
public enum GaugeKind
{
    /// <summary>No factor: couplings are the kernel's (or LMP's) own units.</summary>
    None,

    /// <summary>Every coupling multiplied by 1/K-hat, with K-hat the mean degree of the returned graph:
    /// the classical SPC normalization J = K-hat⁻¹·exp(−d²/2a²), which makes T comparable across data sets.</summary>
    InverseMeanDegree,
}

/// <summary>
/// Which edges the mean-edge-distance bandwidth averages.
/// </summary>
public enum ScaleEdgeSet
{
    /// <summary>The edges the topology and filter produced, before repair: the data's own neighbourhoods.
    /// The repair bridges, the longest edges by construction, do not move every coupling's scale.</summary>
    BeforeRepair,

    /// <summary>Every edge of the repaired graph, bridges included: the classical SPC reading, where the
    /// papers superimpose the minimum spanning tree before averaging.</summary>
    AfterRepair,
}
