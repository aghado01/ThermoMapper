using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Sensitivity;

/// <summary>
/// The pick-freeze design over p inputs (the sampled-sensitivity brief, decision 3): base sets A and B of N rows each, drawn
/// independently from a product law, and for each block u the hybrid A_B^u, A's rows with the columns in u taken from B. The blocks
/// are every input alone, for the first-order and total pieces, and every nonempty subset of each requested subset, for its main
/// variance by Möbius inversion; first-order and total indices cost N(p + 2) values.
/// </summary>
public sealed class PickFreezeDesign
{
    private readonly Dictionary<int[], int> _index;
    private readonly int[][] _terms, _signs;

    /// <summary>The design over <paramref name="inputs"/> inputs, reporting every input alone and each subset in
    /// <paramref name="subsets"/> (input indices, any order, no repeats) with all of its own; refused above
    /// <paramref name="budget"/> blocks.</summary>
    public PickFreezeDesign(int inputs, IReadOnlyList<int[]>? subsets = null, int budget = 10_000)
    {
        if (inputs < 1) throw new ArgumentOutOfRangeException(nameof(inputs), "A design needs an input.");
        if (budget < inputs) throw new ArgumentOutOfRangeException(nameof(budget), $"The budget must hold the {inputs} inputs alone.");
        var found = new HashSet<int[]>(Comparer.Instance);
        for (int i = 0; i < inputs; i++) found.Add([i]);
        if (subsets is not null)
            for (int s = 0; s < subsets.Count; s++)
            {
                int[] u = [.. subsets[s] ?? throw new ArgumentException($"Subset {s} is null.", nameof(subsets))];
                if (u.Length == 0) throw new ArgumentException($"Subset {s} is empty.", nameof(subsets));
                Array.Sort(u);
                for (int k = 0; k < u.Length; k++)
                {
                    if (u[k] < 0 || u[k] >= inputs) throw new ArgumentException($"Subset {s} names input {u[k]}, outside the {inputs} inputs.", nameof(subsets));
                    if (k > 0 && u[k] == u[k - 1]) throw new ArgumentException($"Subset {s} names input {u[k]} twice.", nameof(subsets));
                }
                if (u.Length > 30 || (1L << u.Length) - 1 + inputs > budget)
                    throw new ArgumentException($"Subset {s}'s {u.Length} inputs need more blocks than the budget of {budget}.", nameof(subsets));
                for (int mask = 1; mask < 1 << u.Length; mask++)
                    if (found.Add([.. Enumerable.Range(0, u.Length).Where(k => (mask >> k & 1) != 0).Select(k => u[k])]) && found.Count > budget)
                        throw new ArgumentException($"The subsets need more than {budget} blocks: name fewer or raise the budget.", nameof(subsets));
            }
        var blocks = found.ToList();
        blocks.Sort(Order);
        Inputs = inputs;
        Blocks = Array.AsReadOnly(blocks.ToArray());
        _index = new Dictionary<int[], int>(blocks.Count, Comparer.Instance);
        for (int k = 0; k < blocks.Count; k++) _index[blocks[k]] = k;

        // Each block's Möbius terms: (−1)^{|u|−|v|} on the hybrid of each v ⊆ u, v = ∅ being A itself (−1 here).
        _terms = new int[blocks.Count][];
        _signs = new int[blocks.Count][];
        for (int k = 0; k < blocks.Count; k++)
        {
            int[] u = blocks[k];
            int count = 1 << u.Length;
            _terms[k] = new int[count];
            _signs[k] = new int[count];
            for (int mask = 0; mask < count; mask++)
            {
                int size = System.Numerics.BitOperations.PopCount((uint)mask);
                _terms[k][mask] = mask == 0 ? -1 : _index[[.. Enumerable.Range(0, u.Length).Where(t => (mask >> t & 1) != 0).Select(t => u[t])]];
                _signs[k][mask] = (u.Length - size) % 2 == 0 ? 1 : -1;
            }
        }
    }

    /// <summary>The inputs.</summary>
    public int Inputs { get; }

    /// <summary>The blocks, as ascending input indices, smaller first and then lexicographic, so block i is input i alone for i &lt; p.
    /// Each is also a reported subset.</summary>
    public IReadOnlyList<int[]> Blocks { get; }

    /// <summary>The values one estimate reads at <paramref name="points"/> rows: A, B and each block's hybrid.</summary>
    public long Values(int points) => (long)points * (Blocks.Count + 2);

    /// <summary>The index of block <paramref name="subset"/> (ascending input indices), or −1 where it is not a block.</summary>
    public int IndexOf(int[] subset) => _index.TryGetValue(subset, out int k) ? k : -1;

    /// <summary>The hybrid A_B^u of block <paramref name="block"/>: <paramref name="a"/>'s rows with u's columns from
    /// <paramref name="b"/>.</summary>
    public double[,] Hybrid(double[,] a, double[,] b, int block)
    {
        int n = CheckBases(a, b);
        var hybrid = (double[,])a.Clone();
        foreach (int v in Blocks[block])
            for (int r = 0; r < n; r++) hybrid[r, v] = b[r, v];
        return hybrid;
    }

    internal int CheckBases(double[,] a, double[,] b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        if (a.GetLength(1) != Inputs || b.GetLength(1) != Inputs) throw new ArgumentException($"The base sets have one column per input, {Inputs}.");
        if (a.GetLength(0) != b.GetLength(0)) throw new ArgumentException("The base sets A and B have as many rows.");
        return a.GetLength(0);
    }

    internal int[] Terms(int block) => _terms[block];

    internal int[] Signs(int block) => _signs[block];

    private static int Order(int[] x, int[] y)
    {
        if (x.Length != y.Length) return x.Length.CompareTo(y.Length);
        for (int k = 0; k < x.Length; k++)
            if (x[k] != y[k]) return x[k].CompareTo(y[k]);
        return 0;
    }

    private sealed class Comparer : IEqualityComparer<int[]>
    {
        public static readonly Comparer Instance = new();

        public bool Equals(int[]? x, int[]? y) => x is not null && y is not null && x.AsSpan().SequenceEqual(y);

        public int GetHashCode(int[] obj)
        {
            var hash = new HashCode();
            foreach (int v in obj) hash.Add(v);
            return hash.ToHashCode();
        }
    }
}

/// <summary>A pick-freeze estimate's status (decision 4).</summary>
[Flags]
public enum PickFreezeStatus
{
    /// <summary>The indices are defined.</summary>
    Defined = 0,

    /// <summary>The variance is zero within rounding: the indices and their errors are undefined (NaN), never 0 or 1.</summary>
    ZeroVariance = 1,

    /// <summary>A main or total variance is below zero by more than its error and rounding allow. Nothing is clamped.</summary>
    BelowZeroBeyondError = 2,
}

/// <summary>
/// One function's Sobol pieces from a pick-freeze design (<see cref="PickFreeze.Estimate(PickFreezeDesign, ReadOnlySpan{double},
/// ReadOnlySpan{double}, IReadOnlyList{double[]}, double, double)"/>), each with its Monte Carlo standard error: the mean and variance,
/// each block's main variance (the variance of the function's ANOVA component on it), each input's total variance, and their ratios
/// to the variance, whose errors are the delta method's. Estimates are reported as they are, below zero included.
/// </summary>
public sealed class PickFreezeEstimate
{
    internal PickFreezeEstimate(int points, double mean, double meanError, double variance, double varianceError, double[] main,
                                double[] mainError, double[] total, double[] totalError, double[] mainIndex, double[] mainIndexError,
                                double[] totalIndex, double[] totalIndexError, PickFreezeStatus status)
    {
        Points = points;
        Mean = mean;
        MeanError = meanError;
        Variance = variance;
        VarianceError = varianceError;
        Main = main;
        MainError = mainError;
        Total = total;
        TotalError = totalError;
        MainIndex = mainIndex;
        MainIndexError = mainIndexError;
        TotalIndex = totalIndex;
        TotalIndexError = totalIndexError;
        Status = status;
    }

    /// <summary>N, the rows of each base set.</summary>
    public int Points { get; }

    /// <summary>The mean over A and B.</summary>
    public double Mean { get; }

    /// <summary>The mean's standard error.</summary>
    public double MeanError { get; }

    /// <summary>The variance over A and B.</summary>
    public double Variance { get; }

    /// <summary>The variance's standard error.</summary>
    public double VarianceError { get; }

    /// <summary>Each block's main variance, in <see cref="PickFreezeDesign.Blocks"/>' order.</summary>
    public double[] Main { get; }

    /// <summary>Each main variance's standard error.</summary>
    public double[] MainError { get; }

    /// <summary>Each input's total variance.</summary>
    public double[] Total { get; }

    /// <summary>Each total variance's standard error.</summary>
    public double[] TotalError { get; }

    /// <summary>Each block's main index, the main variance over the variance; NaN where the variance is zero.</summary>
    public double[] MainIndex { get; }

    /// <summary>Each main index's standard error by the delta method; NaN where the variance is zero.</summary>
    public double[] MainIndexError { get; }

    /// <summary>Each input's total index; NaN where the variance is zero.</summary>
    public double[] TotalIndex { get; }

    /// <summary>Each total index's standard error by the delta method; NaN where the variance is zero.</summary>
    public double[] TotalIndexError { get; }

    /// <summary>The status.</summary>
    public PickFreezeStatus Status { get; }
}

/// <summary>
/// Sobol pieces of any function by pick-freeze (the sampled-sensitivity brief, decision 3): it accumulates the function's values on
/// the base sets A and B and on each block's hybrid A_B^u, and reduces them, row by row, to <see cref="PickFreezeEstimate"/>'s pieces.
/// <list type="bullet">
/// <item>A block's main variance is Saltelli's 2010 form by Möbius inversion: the sample covariance, over N − 1, of f(B) with
/// D_u = Σ_{v⊆u} (−1)^{|u|−|v|} f(A_B^v), A_B^∅ = A. f(B) and f(A_B^v) share x_v, and E[D_u] = 0, so it is unbiased for the
/// variance of the ANOVA component on u. Centring f(B) on its own mean makes it shift-invariant; the uncentred form's variance carries
/// 2μ²·T_u, μ the function's mean.</item>
/// <item>An input's total variance is Jansen's 1999 form, Σ (f(A) − f(A_B^i))²/(2N), unbiased.</item>
/// <item>The variance is the sample variance over A and B pooled, and each index is a ratio to it.</item>
/// </list>
/// The rows are independent, so each piece's error is its per-row term's sample deviation over √N, and an index's is the delta
/// method's on its piece and the variance. The forms and their errors are derived in the brief's slice 0 report; the citations are from
/// memory.
/// </summary>
public static class PickFreeze
{
    /// <summary>
    /// The pieces from the values of one function at <paramref name="design"/>'s points: <paramref name="a"/> on A,
    /// <paramref name="b"/> on B, and <paramref name="hybrids"/>, one per block in its order, on the hybrids, each of N values, N ≥ 2.
    /// The variance is zero when it is at most <paramref name="tolerance"/> times its rounding scale, the mean of |y|·|y − ȳ| over A and
    /// B, so a variation below about <paramref name="tolerance"/> times the mean's size is zero; a piece below zero by more than
    /// <paramref name="errorMultiple"/> standard errors and <paramref name="tolerance"/> times that scale is flagged.
    /// </summary>
    public static PickFreezeEstimate Estimate(PickFreezeDesign design, ReadOnlySpan<double> a, ReadOnlySpan<double> b,
                                              IReadOnlyList<double[]> hybrids, double tolerance = 1e-10, double errorMultiple = 3.0)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(hybrids);
        if (!(tolerance >= 0.0) || double.IsPositiveInfinity(tolerance)) throw new ArgumentOutOfRangeException(nameof(tolerance), "The tolerance is nonnegative and finite.");
        if (!(errorMultiple >= 0.0) || double.IsPositiveInfinity(errorMultiple)) throw new ArgumentOutOfRangeException(nameof(errorMultiple), "The error multiple is nonnegative and finite.");
        int n = a.Length, blocks = design.Blocks.Count, p = design.Inputs;
        if (n < 2) throw new ArgumentException("An estimate needs two rows in each base set.", nameof(a));
        if (b.Length != n) throw new ArgumentException($"B has {b.Length} values and A {n}.", nameof(b));
        if (hybrids.Count != blocks) throw new ArgumentException($"The design has {blocks} blocks and {hybrids.Count} hybrids were given.", nameof(hybrids));
        Finite(a, "A");
        Finite(b, "B");
        for (int k = 0; k < blocks; k++)
        {
            if (hybrids[k] is null || hybrids[k].Length != n) throw new ArgumentException($"Hybrid {k} does not hold {n} values.", nameof(hybrids));
            Finite(hybrids[k], $"hybrid {k}");
        }

        // The pooled mean and variance, by the corrected two-pass sums; each row's variance term ψ_j = ½((a_j − m)² + (b_j − m)²).
        double m = 0.0;
        for (int j = 0; j < n; j++) m += a[j] + b[j];
        m /= 2.0 * n;
        double shift = 0.0;
        for (int j = 0; j < n; j++) shift += (a[j] - m) + (b[j] - m);
        m += shift / (2.0 * n);
        double ss = 0.0, s1 = 0.0, size = 0.0;
        var psi = new double[n];
        for (int j = 0; j < n; j++)
        {
            double da = a[j] - m, db = b[j] - m;
            ss += da * da + db * db;
            s1 += da + db;
            size += Math.Abs(a[j]) * Math.Abs(da) + Math.Abs(b[j]) * Math.Abs(db);
            psi[j] = 0.5 * (da * da + db * db);
        }
        double variance = (ss - s1 * s1 / (2.0 * n)) / (2.0 * n - 1.0);
        double meanError = Math.Sqrt(Math.Max(variance, 0.0) / (2.0 * n));
        double varianceError = Deviation(psi) / Math.Sqrt(n);
        bool zero = !(variance > tolerance * size / (2.0 * n));

        // f(B) centred on its own mean, the main pieces' left factor.
        double bm = Mean(b);
        var bc = new double[n];
        double bcSum = 0.0;
        for (int j = 0; j < n; j++)
        {
            bc[j] = b[j] - bm;
            bcSum += bc[j];
        }

        var main = new double[blocks];
        var mainError = new double[blocks];
        var mainIndex = new double[blocks];
        var mainIndexError = new double[blocks];
        var total = new double[p];
        var totalError = new double[p];
        var totalIndex = new double[p];
        var totalIndexError = new double[p];
        var d = new double[n];
        var term = new double[n];
        var g = new double[n];
        for (int k = 0; k < blocks; k++)
        {
            int[] terms = design.Terms(k), signs = design.Signs(k);
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int t = 0; t < terms.Length; t++) s += signs[t] * (terms[t] < 0 ? a[j] : hybrids[terms[t]][j]);
                d[j] = s;
            }
            double dm = Mean(d), dSum = 0.0, cross = 0.0;
            for (int j = 0; j < n; j++)
            {
                double dc = d[j] - dm;
                dSum += dc;
                term[j] = bc[j] * dc;
                cross += term[j];
            }
            main[k] = (cross - bcSum * dSum / n) / (n - 1.0);
            mainError[k] = Deviation(term) / Math.Sqrt(n);
            (mainIndex[k], mainIndexError[k]) = Ratio(main[k], term, psi, variance, zero, g);
        }
        for (int i = 0; i < p; i++)
        {
            double[] h = hybrids[i];   // block i is input i alone
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                double delta = a[j] - h[j];
                term[j] = 0.5 * delta * delta;
                sum += term[j];
            }
            total[i] = sum / n;
            totalError[i] = Deviation(term) / Math.Sqrt(n);
            (totalIndex[i], totalIndexError[i]) = Ratio(total[i], term, psi, variance, zero, g);
        }

        var status = zero ? PickFreezeStatus.ZeroVariance : PickFreezeStatus.Defined;
        double rounding = tolerance * size / (2.0 * n);   // the rounding scale, so a zero-variance function's rounding is not flagged
        for (int k = 0; k < blocks; k++)
            if (main[k] < -(errorMultiple * mainError[k] + rounding)) status |= PickFreezeStatus.BelowZeroBeyondError;
        for (int i = 0; i < p; i++)
            if (total[i] < -(errorMultiple * totalError[i] + rounding)) status |= PickFreezeStatus.BelowZeroBeyondError;
        return new PickFreezeEstimate(n, m, meanError, variance, varianceError, main, mainError, total, totalError, mainIndex, mainIndexError,
                                      totalIndex, totalIndexError, status);
    }

    /// <summary>
    /// The pieces of <paramref name="f"/> at <paramref name="design"/>'s points, from base sets <paramref name="a"/> and
    /// <paramref name="b"/> (N × p, rows drawn independently from the product law): f is read on A, B and each hybrid, the values given
    /// to <see cref="Estimate(PickFreezeDesign, ReadOnlySpan{double}, ReadOnlySpan{double}, IReadOnlyList{double[]}, double, double)"/>.
    /// </summary>
    public static PickFreezeEstimate Estimate(PickFreezeDesign design, Func<ReadOnlySpan<double>, double> f, double[,] a, double[,] b,
                                              double tolerance = 1e-10, double errorMultiple = 3.0)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(f);
        int n = design.CheckBases(a, b), p = design.Inputs;
        var row = new double[p];
        double[] Read(double[,] x, int[]? fromB)
        {
            var values = new double[n];
            for (int j = 0; j < n; j++)
            {
                for (int v = 0; v < p; v++) row[v] = x[j, v];
                if (fromB is not null)
                    foreach (int v in fromB) row[v] = b[j, v];
                values[j] = f(row);
            }
            return values;
        }
        double[] fa = Read(a, null), fb = Read(b, null);
        var hybrids = new double[design.Blocks.Count][];
        for (int k = 0; k < hybrids.Length; k++) hybrids[k] = Read(a, design.Blocks[k]);
        return Estimate(design, fa, fb, hybrids, tolerance, errorMultiple);
    }

    // A piece's ratio to the variance and its delta-method error: the deviation of term_j − r·ψ_j over √N·V.
    private static (double Index, double Error) Ratio(double piece, double[] term, double[] psi, double variance, bool zero, double[] g)
    {
        if (zero) return (double.NaN, double.NaN);
        double r = piece / variance;
        for (int j = 0; j < g.Length; j++) g[j] = term[j] - r * psi[j];
        return (r, Deviation(g) / (Math.Sqrt(g.Length) * variance));
    }

    // The mean by one corrected pass.
    internal static double Mean(ReadOnlySpan<double> x)
    {
        double m = 0.0;
        foreach (double v in x) m += v;
        m /= x.Length;
        double shift = 0.0;
        foreach (double v in x) shift += v - m;
        return m + shift / x.Length;
    }

    // The sample deviation, over n − 1, by the corrected two-pass sums, scaled by the largest deviation: the terms it reads are in
    // squared units already, and squaring them again would leave the range beyond about 10^±77.
    internal static double Deviation(ReadOnlySpan<double> x)
    {
        double m = Mean(x), largest = 0.0;
        foreach (double v in x) largest = Math.Max(largest, Math.Abs(v - m));
        if (largest == 0.0) return 0.0;
        double ss = 0.0, s1 = 0.0;
        foreach (double v in x)
        {
            double c = (v - m) / largest;
            ss += c * c;
            s1 += c;
        }
        return largest * Math.Sqrt(Math.Max(ss - s1 * s1 / x.Length, 0.0) / (x.Length - 1.0));
    }

    private static void Finite(ReadOnlySpan<double> x, string where)
    {
        for (int j = 0; j < x.Length; j++)
            if (!double.IsFinite(x[j])) throw new ArgumentException($"The function's value on {where} at row {j} is not finite.");
    }
}
