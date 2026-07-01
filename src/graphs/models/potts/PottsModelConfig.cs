namespace Graphs.Models.Potts;

/// <summary>
/// Sampler-level configuration for the Potts model — knobs that are
/// chosen once per run and don't vary across the temperature sweep
/// or the chosen scheduling strategy. Lives separately from
/// scheduler/sweep configs so the same sampler settings can be
/// reused regardless of which <c>ISweepStrategy</c> drives them.
/// </summary>
/// <remarks>
/// <para><b>Why separate from the sweep configs.</b> The number of
/// Potts colors <see cref="Q"/> is a property of the Potts model
/// itself, not of how the temperature axis is explored. Holding it
/// inside each sweep config (e.g. <c>FixedGridSweepConfig</c>)
/// duplicated the field and conflated two layers of concern. The
/// sweep configs compose this record instead.</para>
///
/// <para><b>What stays on the sweep configs.</b>
/// <c>BaseSeed</c> remains on the sweep configs because it seeds the
/// schedule's per-task seed derivation, which is a scheduler-level
/// concern (the sampler only ever sees an already-derived per-task
/// seed). The worker-budget policy (<c>Parallelism</c>) likewise lives on
/// the sweep configs and execution options, not the sampler.</para>
/// </remarks>
public sealed class PottsModelConfig
{
    /// <summary>
    /// Number of Potts colors (q). Default 20 — enough that the
    /// paramagnetic <c>⟨1_{s_i = s_j}⟩ ≈ 1/q</c> baseline sits well
    /// below the Blatt cut threshold (typically 0.5), while keeping
    /// per-cycle color-bucket scratch small (O(Q) memory, single
    /// digit ns per pass at q=20).
    /// </summary>
    public int Q { get; init; } = 20;
}
