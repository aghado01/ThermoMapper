using Clustering.Graphical.SPC.Profiling;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Detects pseudo-transitions in an SPC sweep — the temperature points
/// where the system reorganizes between super-paramagnetic phases in the
/// Blatt 1996 / Blatt-Wiseman-Domany 1997 picture.
/// </summary>
/// <remarks>
/// <para><b>Classical signal.</b> The canonical Blatt detector watches
/// the magnetization susceptibility <c>χ_m = N · (⟨m²⟩ − ⟨m⟩²) / T</c>
/// for peaks. Each peak marks a phase boundary where a cluster of size
/// <c>K_i</c> breaks into <c>K_{i+1} &gt; K_i</c> sub-clusters; between
/// consecutive peaks the system is in a stable super-paramagnetic phase
/// and the canonical friends-of-friends cut yields a clean partition.
/// The default implementation
/// (<see cref="MagnetizationPeakDetector"/>) consumes
/// <see cref="SweepProfile.AdditionalChannels"/>'s
/// <c>"MagnetizationVariance"</c> channel — N/T scaling is a strict
/// monotone transform and does not change peak locations on a fixed
/// graph, so the raw variance suffices for relative peak detection.</para>
///
/// <para><b>Pluggable.</b> Concrete detectors may swap the signal
/// (specific heat <c>Cv</c> for energy-axis transitions, label entropy
/// for partition-shape jumps) or the discrimination policy (prominence
/// threshold, peak persistence, multi-signal consensus). All return
/// the same shape: a sorted set of temperatures at which a
/// pseudo-transition was detected.</para>
/// </remarks>
public interface IPseudoTransitionDetector
{
    /// <summary>
    /// Returns the temperatures at which pseudo-transitions were
    /// detected, sorted ascending. Empty when the profile shows no
    /// resolvable phase structure (e.g. a single-component system or a
    /// sweep that didn't bracket any transition).
    /// </summary>
    double[] Detect(SweepProfile profile);
}
