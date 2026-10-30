using System;

namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Raw energy moments at one draw-set: ⟨E⟩, ⟨E²⟩, ⟨E³⟩ and ⟨E⁴⟩.
/// </summary>
public readonly record struct EnergyMoments(double Mean, double SecondMoment, double ThirdMoment, double FourthMoment);

/// <summary>
/// Energy cumulant V of the q-state Potts model, V = 1 − ⟨E⁴⟩/(3⟨E²⟩²) (Challa, Landau and Binder 1986): a reduction
/// of the per-draw energy series. Away from a transition, and at a continuous one, V tends to 2/3 as the system grows;
/// at a first-order transition its minimum stays below 2/3, by an amount the latent heat sets. Model vocabulary;
/// inference- and sweep-agnostic (the V(T) assembly lives in profiling).
/// </summary>
public static class EnergyCumulant
{
    /// <summary>
    /// ⟨E⟩, ⟨E²⟩, ⟨E³⟩ and ⟨E⁴⟩ over a per-draw energy series, each summed in draw order, so ⟨E⟩ and ⟨E²⟩ match the
    /// running sums' means bit for bit. Returns zeros for an empty series.
    /// </summary>
    /// <param name="energyPerDraw">E for each retained draw.</param>
    public static EnergyMoments Moments(ReadOnlySpan<double> energyPerDraw)
    {
        if (energyPerDraw.IsEmpty)
            return new EnergyMoments(0.0, 0.0, 0.0, 0.0);

        double sum = 0.0, sumSq = 0.0, sumCube = 0.0, sumFourth = 0.0;
        foreach (double e in energyPerDraw)
        {
            double sq = e * e;
            sum += e;
            sumSq += sq;
            sumCube += sq * e;
            sumFourth += sq * sq;
        }
        int n = energyPerDraw.Length;
        return new EnergyMoments(sum / n, sumSq / n, sumCube / n, sumFourth / n);
    }

    /// <summary>
    /// V = 1 − ⟨E⁴⟩/(3⟨E²⟩²) from the moments. NaN when ⟨E²⟩ is not positive (no draws, or every draw at E = 0),
    /// where V is undefined.
    /// </summary>
    /// <remarks>
    /// Feed the <b>pooled</b> moments across the replicas at one temperature: V is nonlinear, so it is reduced once on
    /// the pooled moments, never per replica and then averaged.
    /// </remarks>
    /// <param name="secondMoment">⟨E²⟩.</param>
    /// <param name="fourthMoment">⟨E⁴⟩.</param>
    public static double Reduce(double secondMoment, double fourthMoment)
        => secondMoment > 0.0
            ? 1.0 - fourthMoment / (3.0 * secondMoment * secondMoment)
            : double.NaN;

    /// <summary>V of one per-draw energy series.</summary>
    public static double Reduce(ReadOnlySpan<double> energyPerDraw)
    {
        EnergyMoments moments = Moments(energyPerDraw);
        return Reduce(moments.SecondMoment, moments.FourthMoment);
    }
}
