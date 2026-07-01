using System.Text.Json.Serialization;

namespace Graphs.Coupling;

/// <summary>
/// Marker interface for kernel-coupling descriptor value types. Each
/// concrete descriptor (<see cref="Gaussian"/>, <see cref="Cauchy"/>,
/// <see cref="Laplacian"/>, <see cref="Linear"/>, <see cref="Mixture"/>)
/// carries its kernel parameters (bandwidth, mixture weights) and is
/// consumed by the fluent <c>GraphCompiler.CouplingStrategy(...)</c>
/// entry point.
/// </summary>
/// <remarks>
/// <para>The descriptor is decoupled from the math evaluator. Math
/// lives in the <c>*Kernel</c> static classes (<see cref="GaussianKernel"/>,
/// <see cref="CauchyKernel"/>, etc.) so the JIT can dispatch directly
/// at the call site without virtual indirection. The descriptor is a
/// pure value object — the scaler pattern-matches on the descriptor
/// type to select the right evaluator at type-check time.</para>
///
/// <para>A descriptor with <c>Bandwidth = 0.0</c> (the default)
/// signals "auto-estimate via the configured BandwidthEstimation
/// strategy"; any positive value is treated as an explicit override.</para>
///
/// <para>The discriminated-union shape is described <i>once</i>, here, via
/// <see cref="JsonPolymorphicAttribute"/> — STJ source-gen handles the
/// <c>"kind"</c> tag and per-variant fields for both read and write. The former
/// hand-rolled <c>KernelDescriptorJsonConverter</c> (one of several
/// per-output-target re-descriptions of this union) is gone; JSON, the config
/// fingerprint, and the manifest all derive from this single declaration.</para>
/// </remarks>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "kind")]
[JsonDerivedType(typeof(Gaussian), "Gaussian")]
[JsonDerivedType(typeof(Cauchy), "Cauchy")]
[JsonDerivedType(typeof(Laplacian), "Laplacian")]
[JsonDerivedType(typeof(Linear), "Linear")]
[JsonDerivedType(typeof(Mixture), "Mixture")]
public interface IKernelDescriptor { }
