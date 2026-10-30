using System;

namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Potts magnetization moments for the Binder cumulant at one draw-set: ⟨m²⟩ and ⟨m⁴⟩.
/// </summary>
public readonly record struct BinderMoments(double SecondMoment, double FourthMoment);

/// <summary>
/// Binder cumulant U of the q-state Potts model, U = 1 − ⟨m⁴⟩/(3⟨m²⟩²) (Binder 1981): a reduction of the per-draw
/// magnetization series, whose curves for different system sizes cross at a continuous transition's T_c. Model
/// vocabulary; inference- and sweep-agnostic (the U(T) assembly lives in profiling). For q &gt; 4 the 2D transition is
/// first-order and the cumulant dips instead of crossing.
/// </summary>
public static class BinderCumulant
{
    /// <summary>
    /// ⟨m²⟩ and ⟨m⁴⟩ over a per-draw magnetization series, each summed in draw order. Returns zeros for an empty series.
    /// </summary>
    /// <param name="magPerDraw">m for each retained draw.</param>
    public static BinderMoments Moments(ReadOnlySpan<double> magPerDraw)
    {
        if (magPerDraw.IsEmpty)
            return new BinderMoments(0.0, 0.0);

        double sumSq = 0.0, sumFourth = 0.0;
        foreach (double m in magPerDraw)
        {
            double sq = m * m;
            sumSq += sq;
            sumFourth += sq * sq;
        }
        return new BinderMoments(sumSq / magPerDraw.Length, sumFourth / magPerDraw.Length);
    }

    /// <summary>
    /// U = 1 − ⟨m⁴⟩/(3⟨m²⟩²) from the moments. NaN when ⟨m²⟩ is not positive (no draws, or every draw at m = 0),
    /// where U is undefined.
    /// </summary>
    /// <remarks>
    /// Feed the <b>pooled</b> moments across the replicas at one temperature: U is nonlinear, so it is reduced once on
    /// the pooled moments, never per replica and then averaged.
    /// </remarks>
    /// <param name="secondMoment">⟨m²⟩.</param>
    /// <param name="fourthMoment">⟨m⁴⟩.</param>
    public static double Reduce(double secondMoment, double fourthMoment)
        => secondMoment > 0.0
            ? 1.0 - fourthMoment / (3.0 * secondMoment * secondMoment)
            : double.NaN;

    /// <summary>U of one per-draw magnetization series.</summary>
    public static double Reduce(ReadOnlySpan<double> magPerDraw)
    {
        BinderMoments moments = Moments(magPerDraw);
        return Reduce(moments.SecondMoment, moments.FourthMoment);
    }
}
