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

