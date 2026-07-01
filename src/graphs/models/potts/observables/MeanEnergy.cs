namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Mean energy ⟨E⟩ of the q-state Potts model. Model vocabulary; inference- and
/// sweep-agnostic (the ⟨E⟩(T) curve is assembled in profiling).
/// </summary>
public static class MeanEnergy
{
    /// <summary>
    /// Per-draw-set reduction: ⟨E⟩ = <paramref name="runningSumEnergy"/> /
    /// <paramref name="draws"/>. Returns 0 when there are no draws.
    /// </summary>
    /// <param name="runningSumEnergy">Σ over retained draws of the energy.</param>
    /// <param name="draws">Number of retained MC draws averaged (0 ⇒ 0).</param>
    public static double Reduce(double runningSumEnergy, long draws)
        => draws > 0 ? runningSumEnergy / draws : 0.0;
}
