using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Conformal;

/// <summary>
/// One test case's conformal predictive distribution on its grid (E16, slice 4; Vovk, Shen, Manokhin and Xie 2019). The candidate y is
/// ranked by a CDF score among its stratum's m units: Q(y, τ) = (B + τ·(T + 1))/m, with B the units scoring below the candidate and T
/// those tied with it. At the true y, with τ uniform, Q is uniform whatever the score (<c>coverage_tieBreak</c>, <c>rank_uniform</c>):
/// the PIT. Where Q rises along the grid it is also a CDF in y, Vovk et al.'s distribution, with quantiles and CRPS.
/// <list type="bullet">
/// <item>The split form's score, the training fit's predictive CDF F̂(y | x), rises by construction: the default
/// (<see cref="SplitConformal.Distribution"/>; Vovk, Petej, Nouretdinov, Manokhin and Gammerman 2020).</item>
/// <item>The full step's, the in-sample predictive CDF F̃ under the add-one-in weights, falls where the draws' predictives at x
/// disagree (<see cref="ConformalBayes.Distribution"/>, opt-in). Such a case is flagged with how far it falls. It keeps its p-values and
/// set, but has no distribution: its quantiles and CRPS are NaN.</item>
/// </list>
/// The lower and upper p-values at y are Q⁺ = π(y) = (B + T + 1)/m and 1 − Q⁻ = (m − B)/m, ties counted for the candidate in each
/// tail, so the equi-tailed set takes k = ⌊α·m/2⌋ per tail and covers at least 1 − 2k/m, exactly with distinct scores (CB-1 per tail).
/// </summary>
/// <param name="Case">The test case's index in the candidates.</param>
/// <param name="Stratum">Its stratum.</param>
/// <param name="Units">m, the stratum's units with the candidate.</param>
/// <param name="Grid">The case's grid, ascending.</param>
/// <param name="Below">At each grid value, the units whose score is below the candidate's.</param>
/// <param name="Tied">At each grid value, the units whose score equals the candidate's.</param>
/// <param name="Form">The split form or the full step.</param>
/// <param name="Fall">The largest fall of Q⁻ or Q⁺ between neighbouring grid values, in units of 1/m: 0 where Q rises, the
/// distribution's condition.</param>
/// <param name="ScoreFall">The largest fall of the candidate's own score between neighbouring grid values, as a probability: 0 for the
/// split form, and where the full step's draws agree. A fall marks an x where the draws' predictives disagree, a mixing signal.</param>
public sealed record ConformalDistribution(int Case, int Stratum, int Units, double[] Grid, int[] Below, int[] Tied, ConformalForm Form,
                                           int Fall, double ScoreFall)
{
    /// <summary>The optional full-step CDF-score error readouts, including inflated ranks; null on the original API and split form.</summary>
    public ConformalCdfMargin? Margin { get; init; }

    /// <summary>Whether set queries use the inflated tails; distribution, PIT and monotonicity queries retain the original ranks.</summary>
    public bool Inflated { get; init; }

    /// <summary>The conditional fixed/refined-grid coverage diagnostic at alpha; NaN for an uncomputed plain margin.</summary>
    public double CoverageUnderError(double alpha) => WholeSpace(alpha) ? 1.0 : Margin?.CoverageUnderError(Level(alpha), Units, Inflated) ?? double.NaN;

    /// <summary>Whether Q rises along the grid, so the case has a distribution there.</summary>
    public bool Monotone => Fall == 0;

    /// <summary>Q⁻ = Q(y_g, 0) = B/m, the lower distribution.</summary>
    public double Lower(int g) => (double)Below[g] / Units;

    /// <summary>Q⁺ = Q(y_g, 1) = (B + T + 1)/m, the upper distribution, and π(y_g) under the CDF score.</summary>
    public double Upper(int g) => (double)(Below[g] + Tied[g] + 1) / Units;

    /// <summary>Q(y_g, τ) = (B + τ·(T + 1))/m; at the true y with τ uniform, the PIT.</summary>
    public double At(int g, double tau)
    {
        if (!(tau >= 0.0 && tau <= 1.0)) throw new ArgumentOutOfRangeException(nameof(tau), "τ must lie in [0, 1].");
        return (Below[g] + tau * (Tied[g] + 1)) / Units;
    }

    /// <summary>Q's mean over τ, (B + (T + 1)/2)/m, the distribution the quantiles read.</summary>
    public double Mid(int g) => (Below[g] + 0.5 * (Tied[g] + 1)) / Units;

    /// <summary>k = ⌊α·m/2⌋, each tail's level at miscoverage α.</summary>
    public int Level(double alpha)
    {
        CheckAlpha(alpha);
        return (int)Math.Floor(0.5 * alpha * Units);
    }

    /// <summary>1 − 2k/m, the equi-tailed set's coverage.</summary>
    public double Coverage(double alpha) => 1.0 - 2.0 * Level(alpha) / Units;

    /// <summary>Whether the set at α is the whole line: m below ⌈2/α⌉ leaves each tail no level.</summary>
    public bool WholeSpace(double alpha) => Level(alpha) == 0;

    /// <summary>Whether the equi-tailed set at α holds y_g: both p-values above k/m, (B + T + 1) &gt; k and m − B &gt; k.</summary>
    public bool Accepts(int g, double alpha)
        => Inflated ? AcceptsInflated(g, alpha) : AcceptsPlain(g, alpha);

    /// <summary>The original equi-tailed set, irrespective of the inflation choice.</summary>
    public bool AcceptsPlain(int g, double alpha)
    {
        int k = Level(alpha);
        return k == 0 || (Below[g] + Tied[g] + 1 > k && Units - Below[g] > k);
    }

    /// <summary>Each tail inflated by its held CDF-score errors; requires the options overload's error readouts.</summary>
    public bool AcceptsInflated(int g, double alpha)
    {
        int k = Level(alpha);
        if (k == 0) return true;
        if (Margin is null) throw new InvalidOperationException("CDF-score inflation requires the full-step options overload.");
        return Margin.InflatedLowerRanks[g] > k && Margin.InflatedUpperRanks[g] > k;
    }

    /// <summary>The set's maximal runs of accepted grid values, one interval where Q rises; empty for the whole line.</summary>
    public IReadOnlyList<GridInterval> Intervals(double alpha)
    {
        var runs = new List<GridInterval>();
        if (WholeSpace(alpha)) return runs;
        for (int g = 0; g < Grid.Length; g++)
        {
            if (!Accepts(g, alpha)) continue;
            int first = g;
            while (g + 1 < Grid.Length && Accepts(g + 1, alpha)) g++;
            runs.Add(new GridInterval(Grid[first], Grid[g]));
        }
        return runs;
    }

    /// <summary>The set's length on the grid's cells, as <see cref="ConformalSet.Width"/> reads it; +∞ for the whole line.</summary>
    public double Width(double alpha)
    {
        if (WholeSpace(alpha)) return double.PositiveInfinity;
        if (Grid.Length == 1) return Accepts(0, alpha) ? 1.0 : 0.0;
        double width = 0.0;
        for (int g = 0; g < Grid.Length; g++)
            if (Accepts(g, alpha)) width += ConformalBayes.Cell(Grid, g);
        return width;
    }

    /// <summary>Whether the lowest grid value is accepted, so the set may reach below the grid.</summary>
    public bool OpenBelow(double alpha) => !WholeSpace(alpha) && Accepts(0, alpha);

    /// <summary>Whether the highest grid value is accepted, so the set may reach above the grid.</summary>
    public bool OpenAbove(double alpha) => !WholeSpace(alpha) && Accepts(Grid.Length - 1, alpha);

    /// <summary>
    /// The p-quantile: the first grid value where Q's mean over τ reaches p, within the grid's step. −∞ when the first grid value
    /// already reaches it, so the quantile lies at or below the grid; +∞ when none does. NaN for a case with no distribution.
    /// </summary>
    public double Quantile(double p)
    {
        if (!(p > 0.0 && p < 1.0)) throw new ArgumentOutOfRangeException(nameof(p), "p must lie in (0, 1).");
        if (!Monotone) return double.NaN;
        for (int g = 0; g < Grid.Length; g++)
            if (Mid(g) >= p) return g == 0 ? double.NegativeInfinity : Grid[g];
        return double.PositiveInfinity;
    }

    /// <summary>
    /// The CRPS at an observed <paramref name="y"/>: E_τ ∫ (Q(z, τ) − 1{z ≥ y})² dz = ∫ (Q̄ − 1{z ≥ y})² + (Q⁺ − Q⁻)²/12 dz, with
    /// Q̄ = <see cref="Mid"/>. Between grid values Q̄ and Q⁺ − Q⁻ are read linearly, and each piece's square is integrated exactly,
    /// with the cell holding y split there. Q's tails are cut at the grid's ends: the mass below the grid sits at its first value, the
    /// mass above at its last. NaN for a case with no distribution or a grid of one value.
    /// </summary>
    public double Crps(double y)
    {
        if (!double.IsFinite(y)) throw new ArgumentOutOfRangeException(nameof(y), "The observation must be finite.");
        if (!Monotone || Grid.Length < 2) return double.NaN;
        int last = Grid.Length - 1;
        double total = Math.Max(0.0, Grid[0] - y) + Math.Max(0.0, y - Grid[last]);
        for (int g = 0; g < last; g++)
        {
            double a = Grid[g], b = Grid[g + 1], qa = Mid(g), qb = Mid(g + 1), da = Spread(g), db = Spread(g + 1);
            if (y <= a) total += Piece(a, b, qa - 1.0, qb - 1.0, da, db);
            else if (y >= b) total += Piece(a, b, qa, qb, da, db);
            else
            {
                double t = (y - a) / (b - a), qy = qa + t * (qb - qa), dy = da + t * (db - da);
                total += Piece(a, y, qa, qy, da, dy) + Piece(y, b, qy - 1.0, qb - 1.0, dy, db);
            }
        }
        return total;
    }

    // Q⁺ − Q⁻ = (T + 1)/m.
    private double Spread(int g) => (Tied[g] + 1.0) / Units;

    // ∫ u² + d²/12 over [a, b], u and d linear between their end values: (b − a)(u_a² + u_a·u_b + u_b²)/3, and d's alike.
    private static double Piece(double a, double b, double ua, double ub, double da, double db)
        => (b - a) * ((ua * ua + ua * ub + ub * ub) / 3.0 + (da * da + da * db + db * db) / 36.0);

    private static void CheckAlpha(double alpha)
    {
        if (!(alpha > 0.0 && alpha < 1.0)) throw new ArgumentOutOfRangeException(nameof(alpha), "α must lie in (0, 1).");
    }
}

/// <summary>The distributions, one per test case in the candidates' order, from the split form or the full step
/// (<paramref name="Form"/>).</summary>
public sealed record ConformalDistributions(IReadOnlyList<ConformalDistribution> Distributions, ConformalForm Form)
{
    /// <summary>The cases whose score falls along their grid, so they have no distribution: none in the split form.</summary>
    public int Flagged
    {
        get
        {
            int count = 0;
            foreach (ConformalDistribution d in Distributions)
                if (!d.Monotone) count++;
            return count;
        }
    }
}

/// <summary>
/// A CDF score held by both tails: Σ w·F and Σ w·(1 − F), each with full relative precision where it is small, from each draw's
/// ln F (<see cref="PointwiseLogDensity.LogCdf"/>). Two scores compare on the tail that holds their digits: both at or below ½ by
/// their lower sums, both above by their upper sums reversed, otherwise the one at or below ½ is the smaller.
/// </summary>
internal readonly record struct CdfScore(double Lower, double Upper) : IComparable<CdfScore>
{
    public int CompareTo(CdfScore other)
    {
        bool low = Lower <= 0.5, otherLow = other.Lower <= 0.5;
        if (low && otherLow) return Lower.CompareTo(other.Lower);
        if (!low && !otherLow) return other.Upper.CompareTo(Upper);
        return low ? -1 : 1;
    }

    /// <summary>This score less <paramref name="other"/>, from the tail that holds both.</summary>
    public double Minus(CdfScore other) => Lower <= 0.5 && other.Lower <= 0.5 ? Lower - other.Lower : other.Upper - Upper;

    private const double Ln2 = 0.69314718055994530942;

    /// <summary>
    /// F = e^ℓ and 1 − F of a log CDF ℓ. At or below the median 1 − F is at least ½, so it is taken directly and loses nothing; above it
    /// −expm1(ℓ) by Kahan's correction on the same e^ℓ, which keeps 1 − F's digits near 1. One exp a value, and a log above the median.
    /// </summary>
    public static (double F, double Complement) Of(double logCdf)
    {
        double f = Math.Exp(logCdf);
        if (logCdf <= -Ln2) return (f, 1.0 - f);
        return (f, f == 1.0 ? -logCdf : -((f - 1.0) * logCdf / Math.Log(f)));
    }
}

/// <summary>The pieces both forms share: the counts' and the score's falls along a grid, and the CDF scores under fixed weights.</summary>
internal static class CdfScores
{
    /// <summary>The largest fall of B or of B + T between neighbouring grid values, in counts.</summary>
    public static int Fall(int[] below, int[] tied)
    {
        int fall = 0;
        for (int g = 1; g < below.Length; g++)
            fall = Math.Max(fall, Math.Max(below[g - 1] - below[g], below[g - 1] + tied[g - 1] - below[g] - tied[g]));
        return fall;
    }

    /// <summary>The largest fall of the candidate's own score between neighbouring grid values.</summary>
    public static double ScoreFall(CdfScore[] own)
    {
        double fall = 0.0;
        for (int g = 1; g < own.Length; g++)
            if (own[g - 1].CompareTo(own[g]) > 0) fall = Math.Max(fall, own[g - 1].Minus(own[g]));
        return fall;
    }

    /// <summary>Σ_s W̄_s F(y | θ_s) at columns <paramref name="first"/>.. of <paramref name="logCdf"/>, by both tails, under fixed
    /// normalized weights <paramref name="weights"/> on <paramref name="rows"/>: each column summed in row order, the rows walked along
    /// the matrix's own order.</summary>
    public static void Fixed(double[,] logCdf, int[] rows, double[] weights, int first, int count, Span<CdfScore> into)
    {
        var lower = new double[count];
        var upper = new double[count];
        for (int d = 0; d < rows.Length; d++)
        {
            int row = rows[d];
            double w = weights[d];
            for (int c = 0; c < count; c++)
            {
                (double f, double complement) = CdfScore.Of(logCdf[row, first + c]);
                lower[c] += w * f;
                upper[c] += w * complement;
            }
        }
        for (int c = 0; c < count; c++) into[c] = new CdfScore(lower[c], upper[c]);
    }

    /// <summary>The draws' normalized weights on <paramref name="rows"/>: e^(ln W̄) from their normalized logs, 1/S without.</summary>
    public static double[] Weights(int[] rows, double[]? logWeights)
    {
        var w = new double[rows.Length];
        for (int d = 0; d < rows.Length; d++) w[d] = logWeights is null ? 1.0 / rows.Length : Math.Exp(logWeights[d]);
        return w;
    }

    /// <summary>The units of <paramref name="sorted"/> below <paramref name="score"/>, and those equal to it.</summary>
    public static (int Below, int Tied) Rank(CdfScore[] sorted, CdfScore score)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid].CompareTo(score) < 0) lo = mid + 1;
            else hi = mid;
        }
        int below = lo;
        hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid].CompareTo(score) <= 0) lo = mid + 1;
            else hi = mid;
        }
        return (below, lo - below);
    }
}

/// <summary>
/// The full step's CDF score for one stratum against batches of candidates: each member unit's F(y_i | θ_s) and 1 − F under the draws
/// of positive weight, units × draws, the left operands; each batch's normalized add-one-in weights w̃_s ∝ W_s·p(y | x_t, θ_s), draws ×
/// candidates, the right one. Two products give every unit's score by both tails, and the candidate's own is Σ_s w̃_s F(y | θ_s), summed
/// as the products' entries are, so equal columns tie. A candidate no draw supports (every ln p at −∞) has no add-one-in posterior, and
/// is scored under the draws' own weights. Not thread-safe: one per stratum and call.
/// </summary>
internal sealed class CdfScorer
{
    internal double[] CopyWeights(int count) => [.. _weights.AsSpan(0, _s * count)];
    private readonly int[] _rows;
    private readonly double[]? _logWeights;
    private readonly double[] _fixed, _lower, _upper;
    private readonly int _n, _s;
    private double[] _weights = [], _a = [], _b = [], _top = [], _total = [], _ownLower = [], _ownUpper = [];
    private PackedColumns? _right;

    public CdfScorer(double[,] unitLogCdf, int[] rows, double[]? logWeights, int[] members)
    {
        _rows = rows;
        _logWeights = logWeights;
        _n = members.Length;
        _s = rows.Length;
        _fixed = CdfScores.Weights(rows, logWeights);
        _lower = new double[(long)_n * _s];
        _upper = new double[(long)_n * _s];
        // In tiles of 32 draws by 32 units, so the strided reads of the draws × units matrix stay in cache.
        const int tile = 32;
        for (int d0 = 0; d0 < _s; d0 += tile)
            for (int i0 = 0; i0 < _n; i0 += tile)
                for (int i = i0; i < Math.Min(_n, i0 + tile); i++)
                {
                    int unit = members[i];
                    for (int d = d0; d < Math.Min(_s, d0 + tile); d++)
                    {
                        (double f, double complement) = CdfScore.Of(unitLogCdf[rows[d], unit]);
                        _lower[(long)i * _s + d] = f;
                        _upper[(long)i * _s + d] = complement;
                    }
                }
    }

    /// <summary>Each of <paramref name="count"/> candidates from column <paramref name="first"/>: the units below and tied with it, into
    /// <paramref name="below"/> and <paramref name="tied"/> from <paramref name="at"/>, and its own score into <paramref name="own"/>.</summary>
    public void Score(double[,] logDensity, double[,] logCdf, int first, int count, int[] below, int[] tied, CdfScore[] own, int at)
    {
        int s = _s, size = s * count;
        if (_weights.Length < size) _weights = new double[size];
        if (_top.Length < count)
        {
            _top = new double[count];
            _total = new double[count];
            _ownLower = new double[count];
            _ownUpper = new double[count];
        }
        Span<double> w = _weights.AsSpan(0, size), top = _top.AsSpan(0, count), total = _total.AsSpan(0, count);
        Span<double> lower = _ownLower.AsSpan(0, count), upper = _ownUpper.AsSpan(0, count);
        // The rows walked in the matrices' own order; each candidate's sums still run in draw order.
        top.Fill(double.NegativeInfinity);
        for (int d = 0; d < s; d++)
        {
            int row = _rows[d];
            double lw = _logWeights?[d] ?? 0.0;
            for (int b = 0; b < count; b++) top[b] = Math.Max(top[b], lw + logDensity[row, first + b]);
        }
        total.Clear();
        for (int d = 0; d < s; d++)
        {
            int row = _rows[d];
            double lw = _logWeights?[d] ?? 0.0;
            Span<double> wd = w.Slice(d * count, count);
            for (int b = 0; b < count; b++)
                total[b] += wd[b] = double.IsNegativeInfinity(top[b]) ? _fixed[d] : Math.Exp(lw + logDensity[row, first + b] - top[b]);
        }
        // Fused multiply-adds in draw order, as each entry of the products is summed, so a unit whose column equals the candidate's
        // scores exactly as it does, and ties are exact.
        lower.Clear();
        upper.Clear();
        for (int d = 0; d < s; d++)
        {
            int row = _rows[d];
            Span<double> wd = w.Slice(d * count, count);
            for (int b = 0; b < count; b++)
            {
                wd[b] /= total[b];
                (double f, double complement) = CdfScore.Of(logCdf[row, first + b]);
                lower[b] = Math.FusedMultiplyAdd(wd[b], f, lower[b]);
                upper[b] = Math.FusedMultiplyAdd(wd[b], complement, upper[b]);
            }
        }
        for (int b = 0; b < count; b++) own[at + b] = new CdfScore(lower[b], upper[b]);
        Array.Clear(below, at, count);
        Array.Clear(tied, at, count);
        if (_n == 0) return;
        if (_right is null || _right.Columns != count || _right.Inner != s) _right = new PackedColumns(s, count);
        _right.Fill(w);
        if (_a.Length < _n * count)
        {
            _a = new double[_n * count];
            _b = new double[_n * count];
        }
        DenseProduct.Multiply(_lower, _n, _right, _a);
        DenseProduct.Multiply(_upper, _n, _right, _b);
        for (int i = 0; i < _n; i++)
            for (int b = 0; b < count; b++)
            {
                int order = new CdfScore(_a[i * count + b], _b[i * count + b]).CompareTo(own[at + b]);
                if (order < 0) below[at + b]++;
                else if (order == 0) tied[at + b]++;
            }
    }
}
