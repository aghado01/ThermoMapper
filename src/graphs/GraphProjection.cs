using System.Text.Json.Serialization;
using Graphs.Coupling;
using Graphs.Distance;

namespace Graphs;

// =========================================================================
// IEdgeProjection — the L3 stage of the graph spine
// =========================================================================
// The caller DECLARES the projection it wants; GraphCompiler.Build emits
// exactly that and echoes the kind onto GraphBuildResult.WeightKind, so a
// consumer can assert it was handed the weight semantics it expects (Potts
// requires Coupling, HDBSCAN-sparse requires Distance).
//
// This replaces the former mandatory/terminal ScalingConfig: "couple via a
// kernel" is now one projection among several, not the only way to finish a
// build. The conditioned distance graph (L2 output) is reachable directly via
// DistanceProjection.
// =========================================================================

/// <summary>
/// Declares how the conditioned (post-refinement) edge distances become the
/// output edge weights a given consumer needs. See <see cref="EdgeWeightKind"/>.
/// </summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "kind")]
[JsonDerivedType(typeof(DistanceProjection), "distance")]
[JsonDerivedType(typeof(CouplingProjection), "coupling")]
[JsonDerivedType(typeof(AffinityProjection), "affinity")]
public interface IEdgeProjection
{
    [JsonIgnore] EdgeWeightKind Kind { get; }
}

/// <summary>
/// Pass-through: output weights are the conditioned distances themselves —
/// the L2 artifact. Consumed by HDBSCAN-sparse (which lifts to mutual
/// reachability) and Mapper (which wants the structural distance graph and no
/// kernel). Carries no parameters.
/// </summary>
public sealed record DistanceProjection : IEdgeProjection
{
    [JsonIgnore] public EdgeWeightKind Kind => EdgeWeightKind.Distance;
}

/// <summary>
/// Coupling kernel: <c>J = kernel(distance)</c>. The SPC / Potts terminal
/// projection and the role the former <c>ScalingConfig</c> played. The kernel,
/// LMP rescaling, and bandwidth-strategy override live here because they are
/// meaningful only for a kernel projection — they no longer clutter
/// <see cref="GraphCompilerConfig"/> as top-level fields that go inert in
/// distance mode.
/// </summary>
public sealed record CouplingProjection : IEdgeProjection
{
    /// <summary><b>Required.</b> Kernel descriptor driving coupling
    /// computation (Gaussian, Cauchy, Laplacian, Linear, Mixture).</summary>
    public required IKernelDescriptor Kernel { get; init; }

    /// <summary>Local Mutual Proximity post-rescaling. <c>null</c> =
    /// diagnostic-driven; <c>true</c>/<c>false</c> = explicit user assertion.</summary>
    public bool? LmpRescale { get; init; }

    /// <summary>
    /// When LMP runs, preserve couplings on H1 load-bearing edges discovered from the
    /// conditioned distance graph via involuted persistence (default <c>true</c>).
    /// </summary>
    public bool PreserveH1Cycles { get; init; } = true;

    /// <summary>Override for the metric's declared
    /// <see cref="MetricProperties.BandwidthStrategy"/>. <c>null</c> = use the
    /// metric's preference; non-null is an explicit override. (Was
    /// <c>BandwidthConfig.StrategyOverride</c>.)</summary>
    public BandwidthStrategy? BandwidthOverride { get; init; }

    /// <summary>Coupling fidelity — how the metric's curvature enters the kernel.
    /// <c>Auto</c> resolves from geometry at build time.</summary>
    public CouplingFidelity Fidelity { get; init; } = CouplingFidelity.Auto;

    /// <summary>Spherical intrinsic mode — whether to use a local heat kernel parametrix
    /// or a global Schoenberg Gegenbauer kernel. Only active when geometry is Spherical
    /// and fidelity is Intrinsic.</summary>
    public SphericalIntrinsicMode SphericalMode { get; init; } = SphericalIntrinsicMode.Auto;

    [JsonIgnore] public EdgeWeightKind Kind => EdgeWeightKind.Coupling;
}

/// <summary>
/// Affinity (similarity in [0,1], distance-decreasing) for spectral / diffusion
/// consumers. Declared here as the third L3 sibling; <b>not yet wired</b> in
/// <see cref="GraphCompiler.Build"/> — it throws until a consumer earns it.
/// </summary>
public sealed record AffinityProjection : IEdgeProjection
{
    public required IKernelDescriptor Kernel { get; init; }
    [JsonIgnore] public EdgeWeightKind Kind => EdgeWeightKind.Affinity;
}
