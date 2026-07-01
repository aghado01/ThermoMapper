using System;

namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Specific heat C_v of the q-state Potts model — the energy fluctuation scaled by
/// temperature, C_v = (⟨E²⟩ − ⟨E⟩²) / T². Model vocabulary; inference- and
/// sweep-agnostic (the C_v(T) curve assembly lives in profiling).
/// </summary>
public static class SpecificHeat
{
    /// <summary>
    /// Per-draw-set reduction: C_v = max(0, ⟨E²⟩ − ⟨E⟩²) / T², with ⟨E⟩ =
    /// <paramref name="runningSumEnergy"/> / <paramref name="draws"/> and ⟨E²⟩ =
    /// <paramref name="runningSumEnergySq"/> / <paramref name="draws"/>. Returns 0
    /// when there are no draws or T ≤ 0.
    /// </summary>
    /// <param name="runningSumEnergy">Σ over retained draws of the energy.</param>
    /// <param name="runningSumEnergySq">Σ over retained draws of the squared energy.</param>
    /// <param name="draws">Number of retained MC draws averaged (0 ⇒ 0).</param>
    /// <param name="temperature">T &gt; 0.</param>
    public static double Reduce(double runningSumEnergy, double runningSumEnergySq, long draws, double temperature)
    {
        if (draws <= 0 || temperature <= 0.0)
            return 0.0;

        double meanEnergy = runningSumEnergy / draws;
        double meanEnergySq = runningSumEnergySq / draws;
        double variance = Math.Max(0.0, meanEnergySq - meanEnergy * meanEnergy);
        return variance / (temperature * temperature);
    }
}
