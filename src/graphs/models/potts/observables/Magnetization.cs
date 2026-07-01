using System;

namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Potts magnetization moments at one draw-set: the order parameter ⟨m⟩, its
/// second moment ⟨m²⟩, and the variance ⟨m²⟩ − ⟨m⟩² (the magnetic
/// susceptibility χ_m).
/// </summary>
public readonly record struct MagnetizationMoments(double Mean, double SecondMoment, double Variance);

/// <summary>
/// Magnetization M of the q-state Potts model — the order parameter and its
/// fluctuations. Model vocabulary; inference- and sweep-agnostic (the M(T) curves
/// are assembled in profiling).
/// </summary>
public static class Magnetization
{
    /// <summary>
    /// Per-draw-set reduction to ⟨m⟩, ⟨m²⟩, and Var(m), with the moments averaged
    /// over <paramref name="draws"/>. Returns zeros when there are no draws.
    /// </summary>
    /// <param name="runningSumMag">Σ over retained draws of the magnetization.</param>
    /// <param name="runningSumMagSq">Σ over retained draws of the squared magnetization.</param>
    /// <param name="draws">Number of retained MC draws averaged (0 ⇒ zeros).</param>
    public static MagnetizationMoments Reduce(double runningSumMag, double runningSumMagSq, long draws)
    {
        if (draws <= 0)
            return new MagnetizationMoments(0.0, 0.0, 0.0);

        double mean = runningSumMag / draws;
        double secondMoment = runningSumMagSq / draws;
        double variance = Math.Max(0.0, secondMoment - mean * mean);
        return new MagnetizationMoments(mean, secondMoment, variance);
    }
}
