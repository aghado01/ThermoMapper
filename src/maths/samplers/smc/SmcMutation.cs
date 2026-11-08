using System;
using System.Collections.Generic;

namespace Maths.Samplers.Smc;

/// <summary>
/// An adaptive count of kernel moves per stage (<see cref="SmcOptions.AdaptiveMutation"/>), spent from the fixed count's
/// budget. Stage t reads stage t − 1's report: k moves left a fraction 1 − m of its particles where they were
/// (<see cref="SmcStageReport.MovedFraction"/>), so stage t wants ⌈k·ln c / ln(1 − m)⌉ moves, at least 1, the count at which
/// that rate leaves at most a fraction c (<see cref="StayProbability"/>) in place. It takes the want within its credit,
/// k̄·t − Σ_{s&lt;t} k_s with k̄ = <see cref="SmcOptions.MutationSteps"/>, so easy stages pass their unspent moves to hard
/// ones and a run never moves more than the fixed count would over the same stages. Stage 1, with no moves before it, takes
/// k̄. Setting a stage's moves to ⌈ln c / ln(1 − p)⌉ from an acceptance rate p is Drovandi and Pettitt's (2011); reading
/// the moved fraction instead lets particles that never move raise the count. The count is fixed before any of the stage's moves, so
/// every move leaves the stage invariant, and a resumed run reads the same reports.
/// </summary>
public sealed record AdaptiveMutation
{
    /// <summary>c: the fraction of particles a stage's moves may leave where they were.</summary>
    public double StayProbability { get; init; } = 0.1;

    /// <summary>Stage t's moves, t = <paramref name="reports"/>.Count, the reports of the stages before it with the initial
    /// draw first; <paramref name="rate"/> is k̄.</summary>
    public int StepsFor(int rate, IReadOnlyList<SmcStageReport> reports)
    {
        ArgumentNullException.ThrowIfNull(reports);
        ArgumentOutOfRangeException.ThrowIfLessThan(rate, 1);
        int t = reports.Count;
        ArgumentOutOfRangeException.ThrowIfLessThan(t, 1, nameof(reports));
        long credit = (long)rate * t;
        for (int s = 1; s < t; s++) credit -= reports[s].Steps;
        if (t == 1) return (int)Math.Min(rate, credit);
        SmcStageReport previous = reports[t - 1];
        double want = Want(previous.Steps, previous.MovedFraction);
        return (int)Math.Max(1, Math.Min(want, credit));
    }

    // ⌈k·ln c / ln(1 − m)⌉; 1 when every particle moved, unbounded when none did or the stage took no moves to read.
    private double Want(int k, double moved)
    {
        if (k < 1 || double.IsNaN(moved) || moved <= 0.0) return double.PositiveInfinity;
        if (moved >= 1.0) return 1.0;
        return Math.Ceiling(k * Math.Log(StayProbability) / Math.Log(1.0 - moved));
    }
}

/// <summary>
/// A run's stages side by side, the initial draw first, read from its reports: each stage's b, its moves per particle, the
/// kernel's acceptance over them, the fraction of particles they moved, the reweighted population's ESS over N, whether it
/// resampled, and the lineages left. A run whose particles stop moving shows where (<see cref="LeastMoved"/>).
/// </summary>
public sealed record SmcStageReadout(
    double[] InverseTemperature, int[] Steps, double[] AcceptanceRate, double[] MovedFraction, double[] EssFraction, bool[] Resampled,
    int[] Families)
{
    /// <summary>The readout of <paramref name="stages"/>; N is the initial draw's ESS.</summary>
    public static SmcStageReadout Of(IReadOnlyList<SmcStageReport> stages)
    {
        ArgumentNullException.ThrowIfNull(stages);
        if (stages.Count == 0) throw new ArgumentException("There are no stages.", nameof(stages));
        int count = stages.Count;
        double n = stages[0].EffectiveSampleSize;
        var b = new double[count];
        var steps = new int[count];
        var acceptance = new double[count];
        var moved = new double[count];
        var ess = new double[count];
        var resampled = new bool[count];
        var families = new int[count];
        for (int s = 0; s < count; s++)
        {
            SmcStageReport r = stages[s];
            b[s] = r.Temperature.InverseTemperature;
            steps[s] = r.Steps;
            acceptance[s] = r.AcceptanceRate;
            moved[s] = r.MovedFraction;
            ess[s] = r.EffectiveSampleSize / n;
            resampled[s] = r.Resampled;
            families[s] = r.Families;
        }
        return new SmcStageReadout(b, steps, acceptance, moved, ess, resampled, families);
    }

    /// <summary>The moves per particle over the run.</summary>
    public long Moves
    {
        get
        {
            long total = 0;
            foreach (int k in Steps) total += k;
            return total;
        }
    }

    /// <summary>The stage whose moves moved the smallest fraction of particles, and that fraction; (−1, NaN) when no stage
    /// took moves.</summary>
    public (int Stage, double Fraction) LeastMoved
    {
        get
        {
            (int stage, double fraction) = (-1, double.NaN);
            for (int s = 1; s < MovedFraction.Length; s++)
                if (!double.IsNaN(MovedFraction[s]) && (stage < 0 || MovedFraction[s] < fraction)) (stage, fraction) = (s, MovedFraction[s]);
            return (stage, fraction);
        }
    }
}
