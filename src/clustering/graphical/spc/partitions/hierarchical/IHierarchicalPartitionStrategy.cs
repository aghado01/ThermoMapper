using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Hierarchical analogue of <see cref="IPartitionStrategy"/>: consumes
/// the full cross-temperature output of an
/// <see cref="ISweepStrategy"/> rather than a single equilibrium frame,
/// and returns a multi-level <see cref="PartitionHierarchy"/> indexed
/// by temperature (or whichever axis the strategy's detector chose).
/// </summary>
/// <remarks>
/// <para><b>Why a separate seam.</b> The flat
/// <see cref="IPartitionStrategy.Apply"/> contract operates on a single
/// <c>SpcRunResult</c> — the chosen-T equilibrium picked by the sweep.
/// A hierarchical strategy needs the whole sweep (to detect phase
/// boundaries from signal trajectories) <i>and</i> the graph (to
/// re-sample equilibria at phase-representative temperatures with the
/// per-edge observables a cut policy needs). Different inputs, different
/// shape — separate interface.</para>
///
/// <para><b>Typical implementations.</b>
/// <see cref="BlattPartitionStrategy"/> implements the classical Blatt
/// 1996 / Blatt-Wiseman-Domany 1997 picture: detect pseudo-transitions
/// on the χ_m trajectory, run a Tier-1 equilibrium per stable phase,
/// apply the friends-of-friends cut at each, emit the nested
/// sequence.</para>
/// </remarks>
public interface IHierarchicalPartitionStrategy
{
    /// <summary>
    /// Build the hierarchy. <paramref name="sweep"/> is typically the
    /// output of a <see cref="FixedGridSweepStrategy"/> with a
    /// dense T-grid; <paramref name="graph"/> is the same graph that
    /// produced the sweep (used to run per-phase equilibria).
    /// </summary>
    PartitionHierarchy Apply(SweepResult sweep, CsrGraph graph);
}
