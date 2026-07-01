using Clustering.Graphical.SPC.Profiling.Signals;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Strategy seam for SPC sweep schedules. A sweep strategy explores the
/// temperature axis (and any auxiliary axes), assembles a
/// <see cref="Profiling.SweepProfile"/>, picks an equilibrium point, mints
/// the <see cref="SweepResult.ChosenAffinities"/> (and optionally
/// <see cref="SweepResult.ChosenAlignments"/>) currencies, and returns the
/// chosen state plus the probe traces that justified the choice. The
/// end-to-end <see cref="SpcClusteringSession"/> passes the currencies to
/// a partition strategy.
/// </summary>
/// <remarks>
/// <para><b>Implementations.</b> <see cref="FixedGridSweepStrategy"/> is the
/// current bedrock — it runs a user-supplied (or auto-bracketed, via
/// AutoGridFixedSweep) temperature grid. Signal-driven adaptive refinement
/// is parked (see <c>.depr</c>); future strategies may include BARS-style
/// sampling or replica-exchange pre-warmups; all share this contract.</para>
///
/// <para><b>Result type.</b> The return type is
/// <see cref="SweepResult"/> — the shape (summary + sweep runs + chosen-T
/// currencies + profile criteria) is exactly what any sweep strategy must
/// produce. It is the canonical, strategy-agnostic sweep result.</para>
///
/// <para><b>Stateless contract.</b> Implementations should be safe to
/// reuse across multiple <see cref="Run"/> calls on different graphs —
/// any per-run scratch lives inside the call, configuration on the
/// instance.</para>
/// </remarks>
public interface ISweepStrategy
{
    /// <summary>
    /// Execute the sweep and return the chosen equilibrium plus the
    /// probe traces. <paramref name="analyzer"/> defaults to
    /// <see cref="ChiPeakSignalAnalyzer"/> when null and annotates the
    /// result's <see cref="SweepResult.ProfileCriteria"/> — it does not
    /// drive the chosen temperature (the multi-signal consensus variant
    /// is parked in <c>parking-lot/</c> pending the analysis rewrite).
    /// </summary>
    SweepResult Run(CsrGraph graph, ISignalAnalyzer? analyzer = null);
}
