using System;

namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Magnetic susceptibility χ_m of the q-state Potts model — the literal paper
/// susceptibility (BWD eq 5 / Domany): χ_m = (N/T)(⟨m²⟩ − ⟨m⟩²), the order-parameter
/// fluctuation scaled by N/T. Model vocabulary; inference- and sweep-agnostic (the
/// χ_m(T) curve assembly lives in profiling).
/// </summary>
public static class MagnetizationSusceptibility
{
    /// <summary>
    /// χ_m = (N/T)·max(0, ⟨m²⟩ − ⟨m⟩²). Returns 0 when N ≤ 0 or T ≤ 0.
    /// </summary>
    /// <remarks>
    /// Feed the <b>pooled</b> moments (averaged across the whole ensemble at this T), not
    /// a single replica's: the variance is nonlinear, so it must be reduced <i>once</i> on
    /// the pooled moments — never per-frame-then-averaged (Jensen).
    /// </remarks>
    /// <param name="meanMag">⟨m⟩, the sweep-pooled magnetization at this temperature.</param>
    /// <param name="secondMomentMag">⟨m²⟩, the sweep-pooled second moment at this temperature.</param>
    /// <param name="siteCount">N, the number of sites.</param>
    /// <param name="temperature">T &gt; 0.</param>
    public static double Reduce(double meanMag, double secondMomentMag, int siteCount, double temperature)
    {
        if (siteCount <= 0 || temperature <= 0.0)
            return 0.0;

        double variance = Math.Max(0.0, secondMomentMag - meanMag * meanMag);
        return siteCount / temperature * variance;
    }
}
