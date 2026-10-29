using System;

namespace Maths.Samplers.Rjmcmc;

/// <summary>How a parallel-tempering ladder chooses the adjacent pairs it attempts to swap in each scan.</summary>
public enum SwapSchedule
{
    /// <summary>
    /// Deterministic even-odd (DEO; Syed, Bouchard-Côté, Deligiannidis and Doucet 2022): pairs (0,1), (2,3), … on even
    /// scans and (1,2), (3,4), … on odd ones. Non-reversible: a state keeps moving in one direction along the ladder
    /// until a swap is rejected, and the round-trip rate stays bounded below as rungs are added.
    /// </summary>
    DeterministicEvenOdd,

    /// <summary>Stochastic even-odd (SEO): the even or the odd pairs by a fair coin each scan. The reversible baseline,
    /// whose states random-walk along the ladder, so its round-trip rate decays as rungs are added.</summary>
    StochasticEvenOdd,

    /// <summary>Every adjacent pair each scan, cold to hot, later pairs seeing earlier swaps: this engine's schedule
    /// before step 3.3, kept for comparison.</summary>
    Sequential,
}

/// <summary>One adjacent pair's swap statistics: attempts, accepts, and the mean rejection probability over every scan
/// since the ladder took its count of rungs (<see cref="ParallelTempering{TState}.PairScans"/>), attempted or not
/// (Rao–Blackwellized: 1 − α from the cached evaluations, rather than the 0/1 outcomes).</summary>
public sealed record SwapPairStat(double ColderInverseTemperature, double HotterInverseTemperature, long Attempts, long Accepts, double MeanRejection);

/// <summary>
/// Ladder arithmetic from swap rejection rates (Syed et al. 2022). With r_i the rejection rate of adjacent pair i:
/// <list type="bullet">
/// <item>the global communication barrier is Λ ≈ Σ r_i, nearly the same for every ladder (their Corollary 2);</item>
/// <item>under efficient local exploration, the round-trip rate is 1/(2 + 2E) for DEO and 1/(2N + 2E) for SEO, with
/// N the number of pairs and E = Σ r_i/(1 − r_i) (their Theorem 1);</item>
/// <item>for a fixed number of rungs, equal rejection rates minimise E, and N ≈ 2Λ pairs is near the best use of the
/// rungs (their §5.1 and §5.3).</item>
/// </list>
/// A self-sizing ladder reads Λ̂ after each burn-in round, takes the count <see cref="Size"/> gives, places that many rungs
/// on the barrier (<see cref="Place"/>) and starts each new one from the nearest current rung (<see cref="Sources"/>); the
/// shared ladder's <see cref="ParallelTempering{TState}.Resize"/> takes the placement from its caller.
/// </summary>
public static class TemperingLadder
{
    /// <summary>The rungs <see cref="Size"/> lets a reading ask for below the current count before the ladder shrinks:
    /// it shrinks only when <see cref="RecommendedLevels"/> asks for at least this many plus one fewer.</summary>
    public const int ShrinkDeadband = 1;

    /// <summary>The global communication barrier estimate Λ̂ = Σ r_i.</summary>
    public static double Barrier(ReadOnlySpan<double> rejection)
    {
        double sum = 0.0;
        foreach (double r in rejection) sum += r;
        return sum;
    }

    /// <summary>
    /// The round-trip rate per scan, summed over all states, that Theorem 1 predicts under efficient local exploration
    /// (each rung's local move an exact independent draw). Real local moves fall short of it, so the gap between it
    /// and the observed rate measures how far local exploration is from ideal. Undefined for the sequential schedule.
    /// </summary>
    public static double IdealRoundTripRate(ReadOnlySpan<double> rejection, SwapSchedule schedule)
    {
        double inefficiency = 0.0;
        foreach (double r in rejection) inefficiency += r / (1.0 - r);
        return schedule switch
        {
            SwapSchedule.DeterministicEvenOdd => 1.0 / (2.0 + 2.0 * inefficiency),
            SwapSchedule.StochasticEvenOdd => 1.0 / (2.0 * rejection.Length + 2.0 * inefficiency),
            _ => throw new NotSupportedException("Theorem 1 covers the even-odd schedules only."),
        };
    }

    /// <summary>The number of rungs whose DEO round-trip rate per rung is near its best: N ≈ 2Λ pairs, at least 2 rungs.</summary>
    public static int RecommendedLevels(double barrier)
    {
        if (!(barrier >= 0.0) || double.IsPositiveInfinity(barrier))
            throw new ArgumentOutOfRangeException(nameof(barrier));
        return Math.Max(2, (int)Math.Ceiling(2.0 * barrier) + 1);
    }

    /// <summary>
    /// The rung count a self-sizing ladder of <paramref name="current"/> rungs takes after a round that read the barrier
    /// <paramref name="barrier"/>, by <see cref="RecommendedLevels"/>: it grows at once to the reading when that asks for
    /// more, since a saturated ladder reads Λ̂ ≤ its pairs and so under-reads, and keeps growing until the reading stops
    /// asking; it shrinks to the reading only past <see cref="ShrinkDeadband"/>, since early rounds over-read; either way
    /// it is clamped to [<paramref name="minimum"/>, <paramref name="maximum"/>].
    /// </summary>
    public static int Size(double barrier, int current, int minimum, int maximum)
    {
        if (current < 2) throw new ArgumentOutOfRangeException(nameof(current), "A ladder has at least two rungs.");
        if (minimum < 2) throw new ArgumentOutOfRangeException(nameof(minimum), "A ladder has at least two rungs.");
        if (maximum < minimum) throw new ArgumentOutOfRangeException(nameof(maximum), "The maximum lies below the minimum.");
        int asked = RecommendedLevels(barrier);
        int next = asked > current || asked < current - ShrinkDeadband ? asked : current;
        return Math.Clamp(next, minimum, maximum);
    }

    /// <summary>
    /// A ladder with the same rungs and endpoints whose adjacent rejection rates are equal, to the accuracy of the
    /// estimates: <see cref="Place"/> at the current count. <paramref name="inverseTemperatures"/> is descending from the cold rung;
    /// <paramref name="rejection"/>[i] is pair (i, i+1)'s rate. A ladder with no rejections, or only two rungs, comes
    /// back unchanged.
    /// </summary>
    public static double[] Equalize(double[] inverseTemperatures, ReadOnlySpan<double> rejection)
    {
        ArgumentNullException.ThrowIfNull(inverseTemperatures);
        return Place(inverseTemperatures, rejection, inverseTemperatures.Length);
    }

    /// <summary>
    /// <paramref name="levels"/> rungs on <paramref name="inverseTemperatures"/>' endpoints whose adjacent rejection rates are equal, to
    /// the accuracy of the estimates: the cumulative barrier Λ(b), interpolated from the observed rates by a monotone cubic
    /// (Fritsch and Carlson 1980), is inverted at <paramref name="levels"/> − 1 equal steps (Syed et al.'s Algorithms 2 and
    /// 3, at any count). <paramref name="inverseTemperatures"/> is descending from the cold rung; <paramref name="rejection"/>[i] is pair
    /// (i, i+1)'s rate. Where the rates cannot place the rungs (none rejected, or estimates too coarse to keep the rungs
    /// strictly descending), the current count comes back unchanged, as <see cref="Equalize"/> does, and another count is
    /// placed as though every pair rejected alike, by the rungs' index.
    /// </summary>
    public static double[] Place(double[] inverseTemperatures, ReadOnlySpan<double> rejection, int levels)
    {
        ArgumentNullException.ThrowIfNull(inverseTemperatures);
        int n = inverseTemperatures.Length;
        CheckRates(n, rejection);
        if (levels == n && n <= 2) return (double[])inverseTemperatures.Clone();
        if (levels < 2) throw new ArgumentOutOfRangeException(nameof(levels), "A ladder has at least two rungs.");
        if (n < 2) throw new ArgumentException("Placing rungs needs a ladder with two ends.", nameof(inverseTemperatures));

        (double[] x, double[] y) = Cumulative(inverseTemperatures, rejection);
        if (y[n - 1] > 0.0 && Placed(x, y, levels) is { } placed) return placed;
        if (levels == n) return (double[])inverseTemperatures.Clone();
        for (int j = 0; j < n; j++) y[j] = j;   // no reading to place by: every pair alike
        return Placed(x, y, levels) ?? throw new ArgumentException($"The rungs lie too close to place {levels} between them.", nameof(inverseTemperatures));
    }

    /// <summary>
    /// For each rung of <paramref name="placed"/>, the rung of <paramref name="inverseTemperatures"/> whose state it should start from
    /// in a resize (<see cref="ParallelTempering{TState}.Resize"/>): the nearest on the cumulative barrier Λ(b), read from
    /// <paramref name="rejection"/> by the monotone cubic <see cref="Place"/> inverts, so the new rung starts from the state
    /// a swap would hand it most readily; on a tie the nearest in b, and on a further tie the colder. With no rejection the
    /// barrier is flat and the nearest in b decides. A ladder placed on the same ends maps each end to itself.
    /// </summary>
    public static int[] Sources(double[] inverseTemperatures, ReadOnlySpan<double> rejection, double[] placed)
    {
        ArgumentNullException.ThrowIfNull(inverseTemperatures);
        ArgumentNullException.ThrowIfNull(placed);
        int n = inverseTemperatures.Length;
        CheckRates(n, rejection);
        if (n < 1) throw new ArgumentException("A ladder needs a rung.", nameof(inverseTemperatures));
        foreach (double b in placed)
            if (!double.IsFinite(b)) throw new ArgumentException("Every placed b must be finite.", nameof(placed));

        double[]? x = null, y = null, slope = null;
        if (n >= 2)
        {
            (x, y) = Cumulative(inverseTemperatures, rejection);
            if (y[n - 1] > 0.0) slope = MonotoneSlopes(x, y);
        }
        var sources = new int[placed.Length];
        for (int k = 0; k < placed.Length; k++)
        {
            double b = placed[k], at = slope is null ? 0.0 : BarrierAt(x!, y!, slope, b);
            int best = 0;
            double bestBarrier = double.PositiveInfinity, bestDistance = double.PositiveInfinity;
            for (int i = 0; i < n; i++)
            {
                double barrier = slope is null ? 0.0 : Math.Abs(at - y![n - 1 - i]), distance = Math.Abs(b - inverseTemperatures[i]);
                if (barrier < bestBarrier || barrier == bestBarrier && distance < bestDistance)
                    (best, bestBarrier, bestDistance) = (i, barrier, distance);
            }
            sources[k] = best;
        }
        return sources;
    }

    private static void CheckRates(int rungs, ReadOnlySpan<double> rejection)
    {
        if (rejection.Length != rungs - 1)
            throw new ArgumentException("One rejection rate per adjacent pair is required.", nameof(rejection));
        foreach (double r in rejection)
            if (!(r >= 0.0 && r <= 1.0)) throw new ArgumentOutOfRangeException(nameof(rejection), "Rejection rates lie in [0, 1].");
    }

    // Ascending from the hot end: x_j = b, y_j = the cumulative barrier from the hottest rung.
    private static (double[] X, double[] Y) Cumulative(double[] inverseTemperatures, ReadOnlySpan<double> rejection)
    {
        int n = inverseTemperatures.Length;
        var x = new double[n];
        var y = new double[n];
        for (int j = 0; j < n; j++) x[j] = inverseTemperatures[n - 1 - j];
        for (int j = 0; j + 1 < n; j++) y[j + 1] = y[j] + rejection[n - 2 - j];
        return (x, y);
    }

    // `levels` rungs, descending, on the ends of x, the cubic through (x, y) inverted at equal steps of y; null when the
    // inversion does not keep them strictly descending.
    private static double[]? Placed(double[] x, double[] y, int levels)
    {
        int n = x.Length;
        double total = y[n - 1];
        double[] slope = MonotoneSlopes(x, y);
        var ascending = new double[levels];
        ascending[0] = x[0];
        ascending[levels - 1] = x[n - 1];
        int segment = 0;
        for (int k = 1; k < levels - 1; k++)
        {
            double level = total * k / (levels - 1);
            while (segment < n - 2 && !(y[segment + 1] > level)) segment++;
            ascending[k] = Invert(x, y, slope, segment, level);
        }
        var result = new double[levels];
        for (int i = 0; i < levels; i++) result[i] = ascending[levels - 1 - i];
        for (int i = 1; i < levels; i++)
            if (!(result[i] < result[i - 1])) return null;   // estimates too coarse to place the rungs
        return result;
    }

    // Λ at b, the cubic through (x, y), held at its ends outside them.
    private static double BarrierAt(double[] x, double[] y, double[] m, double b)
    {
        int n = x.Length;
        if (!(b > x[0])) return y[0];
        if (!(b < x[n - 1])) return y[n - 1];
        int j = 0;
        while (j < n - 2 && !(x[j + 1] > b)) j++;
        return Hermite(x, y, m, j, b);
    }

    // Fritsch–Carlson tangents for non-decreasing data: the interpolant never leaves [y_j, y_{j+1}] on a segment.
    private static double[] MonotoneSlopes(double[] x, double[] y)
    {
        int n = x.Length;
        var delta = new double[n - 1];
        for (int j = 0; j < n - 1; j++) delta[j] = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
        var m = new double[n];
        m[0] = delta[0];
        m[n - 1] = delta[n - 2];
        for (int j = 1; j < n - 1; j++) m[j] = delta[j - 1] * delta[j] > 0.0 ? 0.5 * (delta[j - 1] + delta[j]) : 0.0;
        for (int j = 0; j < n - 1; j++)
        {
            if (delta[j] == 0.0) { m[j] = 0.0; m[j + 1] = 0.0; continue; }
            double a = m[j] / delta[j], b = m[j + 1] / delta[j];
            double s = a * a + b * b;
            if (s > 9.0)
            {
                double t = 3.0 / Math.Sqrt(s);
                m[j] = t * a * delta[j];
                m[j + 1] = t * b * delta[j];
            }
        }
        return m;
    }

    // The b in segment j at which the monotone cubic reaches level, by bisection.
    private static double Invert(double[] x, double[] y, double[] m, int j, double level)
    {
        double lo = x[j], hi = x[j + 1];
        for (int it = 0; it < 200 && hi - lo > 1e-15 * Math.Max(1.0, Math.Abs(hi)); it++)
        {
            double mid = 0.5 * (lo + hi);
            if (Hermite(x, y, m, j, mid) < level) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    private static double Hermite(double[] x, double[] y, double[] m, int j, double at)
    {
        double h = x[j + 1] - x[j];
        double t = (at - x[j]) / h, t2 = t * t, t3 = t2 * t;
        return (2 * t3 - 3 * t2 + 1) * y[j] + (t3 - 2 * t2 + t) * h * m[j] + (-2 * t3 + 3 * t2) * y[j + 1] + (t3 - t2) * h * m[j + 1];
    }
}
