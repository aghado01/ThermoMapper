using System;
using System.Collections.Generic;
using System.Numerics;

namespace Maths.Regression.Conformal;

/// <summary>
/// A batch prediction set on the score scale (<see cref="BatchConformal"/>): [<paramref name="Lower"/>, <paramref name="Upper"/>] for a
/// function of the batch's m scores, read from the calibration scores at ranks <paramref name="LowerRank"/> − 1 and
/// <paramref name="UpperRank"/> (S_(0) the scores' floor, S_(n+1) their ceiling), and <paramref name="Coverage"/>, the share of the rank
/// sequences it holds, at least 1 − α: the set covers on at least that share of the test sets (<c>batch_coverage</c>).
/// </summary>
public sealed record BatchInterval(double Lower, double Upper, int LowerRank, int UpperRank, double Coverage)
{
    /// <summary>Whether <paramref name="value"/> lies in the closed set.</summary>
    public bool Contains(double value) => value >= Lower && value <= Upper;

    /// <summary>Upper − Lower.</summary>
    public double Width => Upper - Lower;
}

/// <summary>
/// Simultaneous bounds on order statistics of the batch's scores (<see cref="BatchConformal.Simultaneous"/>): for each order t_j, the
/// ranks [<paramref name="LowerRanks"/>_j, <paramref name="UpperRanks"/>_j] and the score bounds [S_(w_j − 1), S_(q_j)]; all hold at once on
/// at least <paramref name="Coverage"/> of the test sets (<c>batch_simultaneous</c>).
/// </summary>
public sealed record BatchBounds(IReadOnlyList<int> Orders, IReadOnlyList<int> LowerRanks, IReadOnlyList<int> UpperRanks,
                                 IReadOnlyList<double> Lower, IReadOnlyList<double> Upper, double Coverage)
{
    /// <summary>Whether every order statistic of <paramref name="batch"/>'s scores lies within its bounds.</summary>
    public bool Contains(IReadOnlyList<double> batch)
    {
        double[] sorted = [.. batch];
        Array.Sort(sorted);
        for (int j = 0; j < Orders.Count; j++)
        {
            double v = sorted[Orders[j] - 1];
            if (v < Lower[j] || v > Upper[j]) return false;
        }
        return true;
    }
}

/// <summary>
/// The selection rule with false claims controlled (<see cref="BatchConformal.Select"/>): the threshold T̂ on the predictions, and the test
/// units whose prediction exceeds it, each declared to have its outcome above c. At most <paramref name="FalseClaims"/> of them have it
/// at or below c, with chance at least 1 − α (<c>batch_exceed</c>, <c>false_claim_iff</c>).
/// </summary>
public sealed record BatchSelection(double Threshold, int Rank, IReadOnlyList<int> Selected, int FalseClaims, double Coverage);

/// <summary>
/// Batch predictive inference (Lee, Tchetgen Tchetgen and Dobriban, arXiv 2409.13990; the conformal layer's slice 16): sets for functions
/// of a batch of m test units' scores, from n calibration scores, under exchangeability of the n + m units. The test units' ranks among the
/// calibration scores, sorted, are uniform on H = {1 ≤ r_1 ≤ … ≤ r_m ≤ n + 1}, |H| = C(n + m, m) (<c>card_testRank_mem</c>), and every
/// sorted test score lies between S_(r_j − 1) and S_(r_j) at once (<c>batch_sandwich</c>, ties broken at random, which moves no calibration
/// order statistic). Each set below is that sandwich read on a set of rank sequences holding at least 1 − α of H.
/// <list type="bullet">
/// <item><see cref="Quantile"/>: the ζ-th smallest score of the batch, from the exact law of r_(ζ) (their quantile corollary).</item>
/// <item><see cref="Simultaneous"/> and <see cref="Box"/>: several order statistics at once (their joint-bounds theorem).</item>
/// <item><see cref="Function"/>: any function h of the sorted batch scores, monotone in each (their Condition 1 and Algorithm 1), by
/// enumerating H; NP-hard in general, so only for a small H.</item>
/// <item><see cref="Select"/>: units declared to have outcome above c, with at most η false claims (their selection corollary).</item>
/// </list>
/// The scores are any per-unit scores from a fit independent of the calibration and test units; the sets bound functions of the batch's
/// scores. S_(0) and S_(n+1) are the scores' floor and ceiling, ±∞ unless given (0 for absolute residuals, say). The rank counts are exact
/// integers, and each level is met exactly: a share of H at most γ is read as the count at most γ·|H| for the given γ.
/// </summary>
public static class BatchConformal
{
    /// <summary>The law of r_(ζ), the ζ-th smallest of the batch's ranks: p(k) = C(k + ζ − 2, ζ − 1)·C(n + m − k − ζ + 1, m − ζ)/C(n + m, m),
    /// k = 1, …, n + 1, at index k − 1 (the order statistic of a sample drawn without replacement).</summary>
    public static double[] RankLaw(int calibration, int batch, int order)
    {
        Check(calibration, batch, order);
        BigInteger[] counts = RankCounts(calibration, batch, order);
        BigInteger total = Binomial(calibration + batch, batch);
        var law = new double[counts.Length];
        for (int k = 0; k < counts.Length; k++) law[k] = Ratio(counts[k], total);
        return law;
    }

    /// <summary>
    /// The batch set for the ζ-th smallest of the batch's m scores (<paramref name="order"/> ζ, 1-based; ⌈(1 − δ)m⌉ for the (1 − δ) quantile):
    /// [S_(q_L − 1), S_(q_U)] with q_U the least rank whose law above it is at most γ and q_L the greatest whose law below it is at most β,
    /// so it covers with chance at least 1 − β − γ, and at most that plus max_k p(k) with distinct scores. The shortest set of its form
    /// (their optimality proposition). At m = 1 it is split conformal's.
    /// </summary>
    /// <param name="scores">The n calibration scores.</param>
    /// <param name="batch">m, the test units.</param>
    /// <param name="order">ζ in 1, …, m.</param>
    /// <param name="lowerError">β, the chance of falling below; 0 gives the floor.</param>
    /// <param name="upperError">γ, the chance of falling above; 0 gives the ceiling.</param>
    /// <param name="floor">S_(0), at most every score.</param>
    /// <param name="ceiling">S_(n+1), at least every score.</param>
    public static BatchInterval Quantile(IReadOnlyList<double> scores, int batch, int order, double lowerError, double upperError,
                                         double floor = double.NegativeInfinity, double ceiling = double.PositiveInfinity)
    {
        double[] sorted = Sorted(scores, floor, ceiling);
        int n = sorted.Length;
        Check(n, batch, order);
        CheckErrors(lowerError, upperError);
        BigInteger[] counts = RankCounts(n, batch, order);
        BigInteger total = Binomial(n + batch, batch);
        (int lower, int upper, BigInteger held) = Ranks(counts, total, lowerError, upperError);
        return new BatchInterval(At(sorted, lower - 1, floor, ceiling), At(sorted, upper, floor, ceiling), lower, upper, Ratio(held, total));
    }

    /// <summary>
    /// The share of H in the box w_j ≤ r_(t_j) ≤ q_j for every j, F/|H| (their F_{n,m}): the chance that every order statistic t_j of the
    /// batch lies in [S_(w_j − 1), S_(q_j)] at once is at least this (<c>batch_simultaneous</c>). Counted exactly over L(ρ), the solutions
    /// with r_(t_j) = ρ_j, with the product over consecutive orders (their closed form's product runs over j = 1, …, n, a typo the Lean
    /// corrects).
    /// </summary>
    /// <param name="calibration">n.</param>
    /// <param name="batch">m.</param>
    /// <param name="orders">t_1 &lt; … &lt; t_l in 1, …, m.</param>
    /// <param name="lowerRanks">w_j in 1, …, n + 1.</param>
    /// <param name="upperRanks">q_j in w_j, …, n + 1.</param>
    public static double BoxCoverage(int calibration, int batch, IReadOnlyList<int> orders, IReadOnlyList<int> lowerRanks, IReadOnlyList<int> upperRanks)
    {
        CheckBox(calibration, batch, orders, lowerRanks, upperRanks);
        return Ratio(BoxCount(calibration, batch, orders, lowerRanks, upperRanks), Binomial(calibration + batch, batch));
    }

    /// <summary>The bounds of a given box (<see cref="BoxCoverage"/>): [S_(w_j − 1), S_(q_j)] for each order t_j, holding at once on at least
    /// the box's share of the test sets.</summary>
    public static BatchBounds Box(IReadOnlyList<double> scores, int batch, IReadOnlyList<int> orders, IReadOnlyList<int> lowerRanks,
                                  IReadOnlyList<int> upperRanks, double floor = double.NegativeInfinity, double ceiling = double.PositiveInfinity)
    {
        double[] sorted = Sorted(scores, floor, ceiling);
        int n = sorted.Length;
        CheckBox(n, batch, orders, lowerRanks, upperRanks);
        double coverage = Ratio(BoxCount(n, batch, orders, lowerRanks, upperRanks), Binomial(n + batch, batch));
        return Bounds(sorted, orders, lowerRanks, upperRanks, coverage, floor, ceiling);
    }

    /// <summary>
    /// Simultaneous bounds on order statistics t_1 &lt; … &lt; t_l of the batch's scores at level 1 − α, by their centred search: each
    /// order's box centred at t̃_j = round(t_j·n/m), clamped to 1, …, n + 1, widened by the least a with F ≥ (1 − α)|H|.
    /// </summary>
    public static BatchBounds Simultaneous(IReadOnlyList<double> scores, int batch, IReadOnlyList<int> orders, double alpha,
                                           double floor = double.NegativeInfinity, double ceiling = double.PositiveInfinity)
    {
        double[] sorted = Sorted(scores, floor, ceiling);
        int n = sorted.Length;
        CheckOrders(batch, orders);
        CheckErrors(alpha, 0.0);
        int l = orders.Count;
        var centre = new int[l];
        for (int j = 0; j < l; j++) centre[j] = Math.Clamp((int)Math.Round((double)orders[j] * n / batch, MidpointRounding.AwayFromZero), 1, n + 1);
        BigInteger total = Binomial(n + batch, batch);
        (int[] W, int[] Q) BoxAt(int a)
        {
            var w = new int[l];
            var q = new int[l];
            for (int j = 0; j < l; j++)
            {
                w[j] = Math.Max(1, centre[j] - a);
                q[j] = Math.Min(n + 1, centre[j] + a);
            }
            return (w, q);
        }
        bool Holds(int a)
        {
            (int[] w, int[] q) = BoxAt(a);
            return AtMost(total - BoxCount(n, batch, orders, w, q), total, alpha);
        }
        int lo = 0, hi = n + 1;   // at a = n + 1 the box is all of H
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (Holds(mid)) hi = mid; else lo = mid + 1;
        }
        (int[] lower, int[] upper) = BoxAt(lo);
        return Bounds(sorted, orders, lower, upper, Ratio(BoxCount(n, batch, orders, lower, upper), total), floor, ceiling);
    }

    /// <summary>
    /// The batch set for any function <paramref name="function"/> h of the batch's sorted scores, monotone non-decreasing in each
    /// (their Condition 1), by their Algorithm 1: with the rank order h̃ (<paramref name="rankOrder"/>; h on the ranks by default), q_L the
    /// greatest value with at most β of H below it and q_U the least with at most γ above, B_L = min h(S_(r − 1)) over h̃(r) ≥ q_L and
    /// B_U = max h(S_(r)) over h̃(r) ≤ q_U; [B_L, B_U] covers with chance at least 1 − β − γ (<c>batch_coverage</c>). It enumerates H twice,
    /// refused above <paramref name="maxSequences"/> sequences: computing it is NP-hard in general (their proposition).
    /// </summary>
    /// <param name="scores">The n calibration scores.</param>
    /// <param name="batch">m.</param>
    /// <param name="function">h, read on a sorted vector of m scores.</param>
    /// <param name="lowerError">β.</param>
    /// <param name="upperError">γ.</param>
    /// <param name="rankOrder">h̃, read on a sorted rank sequence (1-based); h on the ranks as numbers when null.</param>
    /// <param name="floor">S_(0).</param>
    /// <param name="ceiling">S_(n+1).</param>
    /// <param name="maxSequences">The largest |H| enumerated.</param>
    public static BatchInterval Function(IReadOnlyList<double> scores, int batch, Func<double[], double> function, double lowerError,
                                         double upperError, Func<int[], double>? rankOrder = null, double floor = double.NegativeInfinity,
                                         double ceiling = double.PositiveInfinity, long maxSequences = 2_000_000)
    {
        ArgumentNullException.ThrowIfNull(function);
        double[] sorted = Sorted(scores, floor, ceiling);
        int n = sorted.Length;
        Check(n, batch, 1);
        CheckErrors(lowerError, upperError);
        BigInteger size = Binomial(n + batch, batch);
        if (size > maxSequences) throw new ArgumentException($"|H| = C({n + batch}, {batch}) = {size} rank sequences exceed {maxSequences}: batch sets for a general function are NP-hard; the quantile and simultaneous forms serve large batches.", nameof(batch));
        rankOrder ??= r => function(Array.ConvertAll(r, v => (double)v));
        int count = (int)size;
        var values = new double[count];
        int at = 0;
        Enumerate(n, batch, r => values[at++] = rankOrder(r));
        double[] order = (double[])values.Clone();
        Array.Sort(order);
        // q_U: the least value with at most γ of H above it; q_L: the greatest with at most β below it.
        double qU = order[^1], qL = order[0];
        for (int i = 0; i < count; i++)
            if (AtMost(count - CountAtMost(order, order[i]), count, upperError)) { qU = order[i]; break; }
        for (int i = count - 1; i >= 0; i--)
            if (AtMost(CountBelow(order, order[i]), count, lowerError)) { qL = order[i]; break; }
        double low = double.PositiveInfinity, high = double.NegativeInfinity;
        int held = 0;
        var buffer = new double[batch];
        at = 0;
        Enumerate(n, batch, r =>
        {
            double v = values[at++];
            if (v >= qL)
            {
                for (int j = 0; j < batch; j++) buffer[j] = At(sorted, r[j] - 1, floor, ceiling);
                low = Math.Min(low, function(buffer));
            }
            if (v <= qU)
            {
                for (int j = 0; j < batch; j++) buffer[j] = At(sorted, r[j], floor, ceiling);
                high = Math.Max(high, function(buffer));
            }
            if (v >= qL && v <= qU) held++;
        });
        return new BatchInterval(low, high, -1, -1, (double)held / count);
    }

    /// <summary>
    /// Selection with false claims controlled (their selection corollary): with calibration scores s_i = μ̂(x_i)·1{y_i ≤ c}, the threshold
    /// T̂ = S_(q_η), q_η the least rank with at most α of the law of r_(m − η) above it, and the test units whose prediction exceeds T̂ declared
    /// to have outcome above c. At most η of them have outcome at or below c, with chance at least 1 − α: a false claim is a score above T̂
    /// (<c>false_claim_iff</c>), and at most η scores exceed it there (<c>batch_exceed</c>). The predictions are nonnegative, from a fit
    /// independent of the calibration and test units.
    /// </summary>
    /// <param name="calibrationPredictions">μ̂(x_i) at the n calibration units.</param>
    /// <param name="calibrationResponses">y_i at the calibration units.</param>
    /// <param name="cutoff">c.</param>
    /// <param name="testPredictions">μ̂(x) at the m test units.</param>
    /// <param name="falseClaims">η in 0, …, m − 1.</param>
    /// <param name="alpha">α.</param>
    public static BatchSelection Select(IReadOnlyList<double> calibrationPredictions, IReadOnlyList<double> calibrationResponses, double cutoff,
                                        IReadOnlyList<double> testPredictions, int falseClaims, double alpha)
    {
        ArgumentNullException.ThrowIfNull(calibrationPredictions);
        ArgumentNullException.ThrowIfNull(calibrationResponses);
        ArgumentNullException.ThrowIfNull(testPredictions);
        if (calibrationResponses.Count != calibrationPredictions.Count)
            throw new ArgumentException("One response per calibration prediction is required.", nameof(calibrationResponses));
        int n = calibrationPredictions.Count, m = testPredictions.Count;
        if (m == 0) throw new ArgumentException("At least one test unit is required.", nameof(testPredictions));
        ArgumentOutOfRangeException.ThrowIfNegative(falseClaims);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(falseClaims, m);
        CheckErrors(alpha, 0.0);
        if (double.IsNaN(cutoff)) throw new ArgumentException("The cutoff is NaN.", nameof(cutoff));
        var scores = new double[n];
        for (int i = 0; i < n; i++)
        {
            double mu = calibrationPredictions[i], y = calibrationResponses[i];
            if (!(mu >= 0.0 && double.IsFinite(mu))) throw new ArgumentException($"Calibration prediction {i} is {mu}: predictions are finite and nonnegative.", nameof(calibrationPredictions));
            if (double.IsNaN(y)) throw new ArgumentException($"Calibration response {i} is NaN.", nameof(calibrationResponses));
            scores[i] = y <= cutoff ? mu : 0.0;
        }
        for (int j = 0; j < m; j++)
            if (!(testPredictions[j] >= 0.0 && double.IsFinite(testPredictions[j])))
                throw new ArgumentException($"Test prediction {j} is {testPredictions[j]}: predictions are finite and nonnegative.", nameof(testPredictions));
        Array.Sort(scores);
        int order = m - falseClaims;
        BigInteger[] counts = RankCounts(n, m, order);
        BigInteger total = Binomial(n + m, m);
        (_, int rank, _) = Ranks(counts, total, 0.0, alpha);
        BigInteger held = 0;
        for (int k = 1; k <= rank; k++) held += counts[k - 1];
        double threshold = At(scores, rank, 0.0, double.PositiveInfinity);
        var selected = new List<int>();
        for (int j = 0; j < m; j++) if (testPredictions[j] > threshold) selected.Add(j);
        return new BatchSelection(threshold, rank, selected, falseClaims, Ratio(held, total));
    }

    // ---- The counts ----

    // |{r ∈ H : r_(ζ) = k}| for k = 1, …, n + 1: C(k + ζ − 2, ζ − 1)·C(n + m − k − ζ + 1, m − ζ).
    private static BigInteger[] RankCounts(int n, int m, int zeta)
    {
        var counts = new BigInteger[n + 1];
        for (int k = 1; k <= n + 1; k++) counts[k - 1] = Multiset(k, zeta - 1) * Multiset(n - k + 2, m - zeta);
        return counts;
    }

    // F: the rank sequences with w_j ≤ r_(t_j) ≤ q_j for every j, by dynamic programming over ρ_j = r_(t_j) on L(ρ) =
    // ₍ρ₁₎H₍t₁₋₁₎ · Π_j ₍ρ_{j+1} − ρ_j + 1₎H₍t_{j+1} − t_j − 1₎ · ₍n − ρ_l + 2₎H₍m − t_l₎.
    private static BigInteger BoxCount(int n, int m, IReadOnlyList<int> orders, IReadOnlyList<int> lower, IReadOnlyList<int> upper)
    {
        int l = orders.Count;
        var f = new BigInteger[n + 2];
        for (int rho = lower[0]; rho <= upper[0]; rho++) f[rho] = Multiset(rho, orders[0] - 1);
        for (int j = 1; j < l; j++)
        {
            int gap = orders[j] - orders[j - 1] - 1;
            var kernel = new BigInteger[n + 1];
            for (int d = 0; d <= n; d++) kernel[d] = Multiset(d + 1, gap);
            var next = new BigInteger[n + 2];
            for (int rho = lower[j]; rho <= upper[j]; rho++)
            {
                BigInteger sum = BigInteger.Zero;
                for (int prev = lower[j - 1]; prev <= Math.Min(upper[j - 1], rho); prev++)
                    if (!f[prev].IsZero) sum += f[prev] * kernel[rho - prev];
                next[rho] = sum;
            }
            f = next;
        }
        BigInteger count = BigInteger.Zero;
        for (int rho = lower[l - 1]; rho <= upper[l - 1]; rho++)
            if (!f[rho].IsZero) count += f[rho] * Multiset(n - rho + 2, m - orders[l - 1]);
        return count;
    }

    // q_L, the greatest rank with at most β of the law below it, and q_U, the least with at most γ above it; and the count they hold.
    private static (int Lower, int Upper, BigInteger Held) Ranks(BigInteger[] counts, BigInteger total, double beta, double gamma)
    {
        int ranks = counts.Length;
        int upper = ranks;
        BigInteger above = total;
        for (int k = 1; k <= ranks; k++)
        {
            above -= counts[k - 1];
            if (AtMost(above, total, gamma)) { upper = k; break; }
        }
        int lower = 1;
        BigInteger below = total;
        for (int k = ranks; k >= 1; k--)
        {
            below -= counts[k - 1];   // the count strictly below k
            if (AtMost(below, total, beta)) { lower = k; break; }
        }
        BigInteger held = BigInteger.Zero;
        for (int k = lower; k <= upper; k++) held += counts[k - 1];
        return (lower, upper, held);
    }

    // ₍a₎H₍b₎ = C(a + b − 1, b), the multisets of size b from a kinds; 1 at b = 0, 0 for a ≤ 0 with b > 0.
    private static BigInteger Multiset(int a, int b)
    {
        if (b == 0) return BigInteger.One;
        if (a <= 0) return BigInteger.Zero;
        return Binomial(a + b - 1, b);
    }

    private static BigInteger Binomial(int n, int k)
    {
        if (k < 0 || k > n) return BigInteger.Zero;
        k = Math.Min(k, n - k);
        BigInteger result = BigInteger.One;
        for (int i = 1; i <= k; i++) result = result * (n - k + i) / i;
        return result;
    }

    // count ≤ level·total exactly, the level read as the double it is.
    private static bool AtMost(BigInteger count, BigInteger total, double level)
    {
        if (level >= 1.0) return count <= total;
        if (level <= 0.0) return count <= BigInteger.Zero;
        long bits = BitConverter.DoubleToInt64Bits(level);
        int exponent = (int)((bits >> 52) & 0x7FF);
        long fraction = bits & ((1L << 52) - 1);
        long mantissa = exponent == 0 ? fraction : fraction | (1L << 52);
        int shift = exponent == 0 ? 1074 : 1075 - exponent;   // level = mantissa·2^(−shift)
        return (count << shift) <= mantissa * total;
    }

    private static double Ratio(BigInteger count, BigInteger total)
    {
        if (count.IsZero) return 0.0;
        double log = BigInteger.Log(count) - BigInteger.Log(total);
        return Math.Exp(log);
    }

    // Every non-decreasing sequence 1 ≤ r_1 ≤ … ≤ r_m ≤ n + 1, in lexicographic order.
    private static void Enumerate(int n, int m, Action<int[]> visit)
    {
        var r = new int[m];
        void Fill(int j, int start)
        {
            if (j == m) { visit(r); return; }
            for (int v = start; v <= n + 1; v++)
            {
                r[j] = v;
                Fill(j + 1, v);
            }
        }
        Fill(0, 1);
    }

    private static int CountAtMost(double[] sorted, double v)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi) { int mid = (lo + hi) >>> 1; if (sorted[mid] <= v) lo = mid + 1; else hi = mid; }
        return lo;
    }

    private static int CountBelow(double[] sorted, double v)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi) { int mid = (lo + hi) >>> 1; if (sorted[mid] < v) lo = mid + 1; else hi = mid; }
        return lo;
    }

    private static BatchBounds Bounds(double[] sorted, IReadOnlyList<int> orders, IReadOnlyList<int> lower, IReadOnlyList<int> upper, double coverage,
                                      double floor, double ceiling)
    {
        int l = orders.Count;
        var lo = new double[l];
        var hi = new double[l];
        for (int j = 0; j < l; j++)
        {
            lo[j] = At(sorted, lower[j] - 1, floor, ceiling);
            hi[j] = At(sorted, upper[j], floor, ceiling);
        }
        return new BatchBounds([.. orders], [.. lower], [.. upper], lo, hi, coverage);
    }

    // S_(k): the floor at 0, the ceiling at n + 1.
    private static double At(double[] sorted, int k, double floor, double ceiling) => k <= 0 ? floor : k > sorted.Length ? ceiling : sorted[k - 1];

    private static double[] Sorted(IReadOnlyList<double> scores, double floor, double ceiling)
    {
        ArgumentNullException.ThrowIfNull(scores);
        if (double.IsNaN(floor) || double.IsNaN(ceiling) || floor > ceiling) throw new ArgumentException("The floor and ceiling are ordered numbers.", nameof(floor));
        var sorted = new double[scores.Count];
        for (int i = 0; i < sorted.Length; i++)
        {
            double v = scores[i];
            if (double.IsNaN(v)) throw new ArgumentException($"Score {i} is NaN.", nameof(scores));
            if (v < floor || v > ceiling) throw new ArgumentException($"Score {i} is {v}, outside the floor {floor} and ceiling {ceiling}.", nameof(scores));
            sorted[i] = v;
        }
        Array.Sort(sorted);
        return sorted;
    }

    private static void Check(int n, int m, int order)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfLessThan(m, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(order, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(order, m);
    }

    private static void CheckErrors(double lower, double upper)
    {
        if (!(lower >= 0.0 && upper >= 0.0 && lower + upper < 1.0))
            throw new ArgumentOutOfRangeException(nameof(lower), "The error levels are nonnegative with a sum below 1.");
    }

    private static void CheckOrders(int m, IReadOnlyList<int> orders)
    {
        ArgumentNullException.ThrowIfNull(orders);
        ArgumentOutOfRangeException.ThrowIfLessThan(m, 1);
        if (orders.Count == 0) throw new ArgumentException("At least one order is required.", nameof(orders));
        for (int j = 0; j < orders.Count; j++)
            if (orders[j] < 1 || orders[j] > m || (j > 0 && orders[j] <= orders[j - 1]))
                throw new ArgumentException("The orders rise strictly within 1, …, m.", nameof(orders));
    }

    private static void CheckBox(int n, int m, IReadOnlyList<int> orders, IReadOnlyList<int> lower, IReadOnlyList<int> upper)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        CheckOrders(m, orders);
        ArgumentNullException.ThrowIfNull(lower);
        ArgumentNullException.ThrowIfNull(upper);
        if (lower.Count != orders.Count || upper.Count != orders.Count) throw new ArgumentException("One rank bound pair per order is required.", nameof(lower));
        for (int j = 0; j < orders.Count; j++)
            if (lower[j] < 1 || upper[j] > n + 1 || lower[j] > upper[j])
                throw new ArgumentException($"Order {orders[j]}'s ranks [{lower[j]}, {upper[j]}] must lie within 1, …, n + 1, in order.", nameof(lower));
    }
}
