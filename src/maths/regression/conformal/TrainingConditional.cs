using System;
using Maths.Distributions;

namespace Maths.Regression.Conformal;

/// <summary>
/// A split set's own miscoverage bound (the conformal layer's slice 18 follow-up): the split form's training-conditional law (Vovk 2012,
/// "Conditional validity of inductive conformal predictors", Proposition 2a; Bian and Barber 2023). A test unit of a stratum with n
/// calibration units is accepted when its conformity reaches S₍ₖ₎, the k-th smallest of theirs, k = ⌊α·m⌋, m = n + 1 (CB-1). Given the
/// calibration units, the sets' miscoverage over new test units is F(S₍ₖ₎−), F the conformity's law, which the probability integral
/// transform puts at or below U₍ₖ₎, the k-th smallest of n uniforms: Beta(k, m − k), that law exactly where the scores are continuous,
/// stochastically larger than the miscoverage where they tie. So with chance at least β over the calibration units the miscoverage is at
/// most Beta(k, m − k)'s β-quantile. Its mean is k/m, CB-1's marginal reading.
/// <para>Gaffke's b*_β (<see cref="LossBound"/>) bounds a rule fixed before its losses are read, such as a validation split's misses at a
/// set already chosen; a set whose threshold is the calibration's own k-th score needs this law. It is the split form's alone. The full
/// step and the cross-fitted forms rank each case among units whose scores depend on it, and carry no law of this kind (Bian and Barber);
/// the weighted rule's threshold moves with each case's weights.</para>
/// </summary>
public static class TrainingConditional
{
    /// <summary>The β-quantile of Beta(k, m − k) at <paramref name="units"/> m and <paramref name="level"/> k: with chance at least
    /// <paramref name="beta"/> over the calibration units, the miscoverage of the sets they give is at most this. 0 at k = 0, where every
    /// set is the whole space.</summary>
    /// <param name="units">m, the calibration units and the test unit.</param>
    /// <param name="level">k = ⌊α·m⌋, below m.</param>
    /// <param name="beta">The confidence β, in (0, 1), as <see cref="LossBound.Gaffke(double)"/> takes it.</param>
    public static double Miscoverage(int units, int level, double beta)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(units, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(level);
        if (level >= units) throw new ArgumentOutOfRangeException(nameof(level), $"k = {level} lies below m = {units}.");
        if (!(beta > 0.0 && beta < 1.0)) throw new ArgumentOutOfRangeException(nameof(beta), "β lies in (0, 1).");
        return level == 0 ? 0.0 : Beta.Shapes(level, units - level).Quantile(beta);
    }

    /// <summary>The (α, β) level at <paramref name="units"/> m: the largest k at most ⌊α·m⌋ whose Beta(k, m − k) β-quantile is at most
    /// α, so the sets it gives miss with chance at most α with chance at least <paramref name="beta"/> over the calibration units, and cover
    /// at least 1 − k/m ≥ 1 − α marginally (Vovk 2012's (ε, δ)-validity). 0, the whole space, where even k = 1 passes α.</summary>
    public static int Level(int units, double alpha, double beta) => Level(units, alpha, beta, 1);

    // Each of `tails` tails' k, the tails reading one score: the largest k at most ⌊α·m/tails⌋ whose Beta(tails·k, m − tails·k) β-quantile
    // is at most α, by bisection, the quantile rising in k.
    internal static int Level(int units, double alpha, double beta, int tails)
    {
        SplitConformal.CheckAlpha(alpha);
        ArgumentOutOfRangeException.ThrowIfLessThan(units, 1);
        if (!(beta > 0.0 && beta < 1.0)) throw new ArgumentOutOfRangeException(nameof(beta), "β lies in (0, 1).");
        int lo = 0, hi = (int)Math.Floor(alpha * units / tails);
        if (hi == 0 || Miscoverage(units, tails * hi, beta) <= alpha) return hi;
        while (hi - lo > 1)
        {
            int mid = lo + (hi - lo) / 2;
            if (Miscoverage(units, tails * mid, beta) <= alpha) lo = mid; else hi = mid;
        }
        return lo;
    }
}
