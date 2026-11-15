using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Conformal;

/// <summary>Which law on the grid's cells the proper scores read (<see cref="ConformalDistribution.CellMasses"/>).</summary>
public enum CellLaw
{
    /// <summary>The predictive's own law, Q̄'s rises between grid values: a step function of the calibration ranks, so a cell no calibration
    /// PIT rises across has no mass.</summary>
    Steps,

    /// <summary>The smoothed law (<see cref="SmoothedCdf"/>): positive on every cell where the fit's CDF rises, so the log score is finite.
    /// The split form's, unweighted.</summary>
    Smoothed,
}

/// <summary>
/// The split form's smoothed predictive CDF along a case's grid (slice 18's follow-up): G̃(F̂(y | x)), G̃ the calibration PITs' empirical
/// CDF interpolated linearly between its knots (0, 0), (PIT₍ᵢ₎, i/(n + 1)) and (1, 1). Between neighbouring PITs it rises where the step
/// law Q̄ is flat, and it lies inside the conformal band [Q⁻, Q⁺] at every y: it is Q(y, τ) with τ moving linearly across each gap. Where
/// the fit's CDF rises every cell has mass, so the log score is finite (<c>logLoss_proper</c>'s positivity hypothesis met). Its PIT at the
/// truth is near uniform, not exactly so: the randomized Q keeps the exact guarantee; this law is for scoring.
/// </summary>
/// <param name="Lower">G̃ at each grid value.</param>
/// <param name="Upper">1 − G̃ at each grid value, read on the tail that holds its digits.</param>
public sealed record SmoothedCdf(double[] Lower, double[] Upper);

/// <summary>
/// A case's weighted ranks along its grid (slice 14's follow-up 1; <see cref="ShiftWeights"/>): at each grid value the share of its
/// stratum's weight, the candidate's own included, that lies below the candidate, Q⁻ = A/(W + w_c), and at or below it with its own,
/// Q⁺ = (A + B)/(W + w_c), B the tied weight with the candidate's. The PIT T = Q⁻ + V·(Q⁺ − Q⁻), V uniform, is exactly uniform at the test
/// unit's value under covariate shift w(x) or response shift w(y) with the true ratio, in the split form and the full step
/// (<c>pit_weighted</c>, <c>sum_mul_wPitCdf</c>); with an estimated ratio it lies within the data set's total variation between the true and
/// the estimated normalized weights of uniform, from either side (<c>wPit_estimated</c>). With ties counted for the candidate (V = 1) it is
/// <c>weightLe</c>, super-uniform at every level (<c>coverage_weighted</c>). Where W + w_c is 0, Q⁻ is 0 and Q⁺ is 1.
/// </summary>
/// <param name="Lower">Q⁻ at each grid value.</param>
/// <param name="Upper">Q⁺ at each grid value.</param>
/// <param name="Fall">The largest fall of Q⁻ or Q⁺ between neighbouring grid values: 0 where both rise, the weighted distribution's
/// condition. Under response shift the candidate's weight moves along the grid, and the shares can fall where the counts rise.</param>
public sealed record WeightedRanks(double[] Lower, double[] Upper, double Fall);

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
/// A case weighed by unequal shift weights reads its <see cref="Weighted"/> ranks instead, in every query.
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

    /// <summary>
    /// The weighted ranks where the case was weighed by unequal shift weights (<see cref="ShiftWeights"/>); null where unweighted, or where
    /// its stratum's weights and its candidate weights are all one value, CB-1's distribution bit for bit. Where set every query reads them:
    /// Q⁻, Q⁺, the PIT, the quantiles, the CRPS and monotonicity, and the equi-tailed set, which holds y_g where Q⁺ &gt; α/2 and
    /// 1 − Q⁻ &gt; α/2, each tail refusing at most α/2 of the weight, so it covers at least 1 − α.
    /// </summary>
    public WeightedRanks? Weighted { get; init; }

    /// <summary>The smoothed predictive CDF along the grid (<see cref="SmoothedCdf"/>, slice 18's follow-up), the split form's where the case is
    /// unweighted; null for the full step and for a weighed case.</summary>
    public SmoothedCdf? Smoothed { get; init; }

    /// <summary>The conditional fixed/refined-grid coverage diagnostic at alpha; NaN for an uncomputed plain margin. Where weighed
    /// (slice 14's follow-up 3): 1 − α − ν for the inflated set, 1 − α − M̂ − ν for the plain one (<see cref="ConformalCdfMargin.WeightedBand"/>).</summary>
    public double CoverageUnderError(double alpha)
    {
        if (WholeSpace(alpha)) return 1.0;
        if (Weighted is not null) return Margin?.WeightedCoverageUnderError(Checked(alpha), Inflated) ?? double.NaN;
        return Margin?.CoverageUnderError(Level(alpha), Units, Inflated) ?? double.NaN;
    }

    /// <summary>Whether Q rises along the grid, so the case has a distribution there: the counts' <see cref="Fall"/>, or the weighted
    /// ranks' where the case was weighed.</summary>
    public bool Monotone => Weighted is { } w ? w.Fall == 0.0 : Fall == 0;

    /// <summary>Q⁻ = Q(y_g, 0) = B/m, the lower distribution; A/(W + w_c) where weighed.</summary>
    public double Lower(int g) => Weighted is { } w ? w.Lower[g] : (double)Below[g] / Units;

    /// <summary>Q⁺ = Q(y_g, 1) = (B + T + 1)/m, the upper distribution, and π(y_g) under the CDF score; (A + B)/(W + w_c) where weighed.</summary>
    public double Upper(int g) => Weighted is { } w ? w.Upper[g] : (double)(Below[g] + Tied[g] + 1) / Units;

    /// <summary>Q(y_g, τ) = (B + τ·(T + 1))/m, or Q⁻ + τ·(Q⁺ − Q⁻) where weighed; at the true y with τ uniform, the PIT.</summary>
    public double At(int g, double tau)
    {
        if (!(tau >= 0.0 && tau <= 1.0)) throw new ArgumentOutOfRangeException(nameof(tau), "τ must lie in [0, 1].");
        return Weighted is { } w ? w.Lower[g] + tau * (w.Upper[g] - w.Lower[g]) : (Below[g] + tau * (Tied[g] + 1)) / Units;
    }

    /// <summary>Q's mean over τ, (B + (T + 1)/2)/m, or (Q⁻ + Q⁺)/2 where weighed: the distribution the quantiles read.</summary>
    public double Mid(int g) => Weighted is { } w ? 0.5 * (w.Lower[g] + w.Upper[g]) : (Below[g] + 0.5 * (Tied[g] + 1)) / Units;

    /// <summary>k = ⌊α·m/2⌋, each tail's level at miscoverage α: the counts', which a weighed case's set does not read.</summary>
    public int Level(double alpha)
    {
        CheckAlpha(alpha);
        return (int)Math.Floor(0.5 * alpha * Units);
    }

    /// <summary>1 − 2k/m, the equi-tailed set's coverage; 1 − α where weighed.</summary>
    public double Coverage(double alpha) => Weighted is null ? 1.0 - 2.0 * Level(alpha) / Units : 1.0 - Checked(alpha);

    /// <summary>
    /// The equi-tailed set's own miscoverage bound at α under the split form (<see cref="TrainingConditional"/>). Both tails read one
    /// score, the PIT, so given the calibration a miss is the chance mass of 2k of the m spacings its PITs cut, Beta(2k, m − 2k), that law
    /// exactly where the PITs are continuous and conservative where they tie: with chance at least <paramref name="beta"/> over the
    /// calibration units the set's miscoverage is at most its β-quantile. Null for the full step, a weighed case and the inflated set.
    /// </summary>
    public double? TrainingConditionalMiscoverage(double alpha, double beta)
        => Form == ConformalForm.Split && Weighted is null && !Inflated ? TrainingConditional.Miscoverage(Units, 2 * Level(alpha), beta) : null;

    /// <summary>Whether the set at α is the whole line: m below ⌈2/α⌉ leaves each tail no level. A weighed case judges each grid value
    /// by its shares, so it is never the whole line by its count.</summary>
    public bool WholeSpace(double alpha) => Level(alpha) == 0 && Weighted is null;

    /// <summary>Whether the equi-tailed set at α holds y_g: both p-values above k/m, (B + T + 1) &gt; k and m − B &gt; k.</summary>
    public bool Accepts(int g, double alpha)
        => Inflated ? AcceptsInflated(g, alpha) : AcceptsPlain(g, alpha);

    /// <summary>The original equi-tailed set, irrespective of the inflation choice.</summary>
    public bool AcceptsPlain(int g, double alpha)
    {
        if (Weighted is { } w) return w.Upper[g] > 0.5 * Checked(alpha) && 1.0 - w.Lower[g] > 0.5 * alpha;
        int k = Level(alpha);
        return k == 0 || (Below[g] + Tied[g] + 1 > k && Units - Below[g] > k);
    }

    /// <summary>Each tail inflated by its held CDF-score errors; requires the options overload's error readouts.</summary>
    public bool AcceptsInflated(int g, double alpha)
    {
        if (Weighted is not null)
        {
            // The weighted inflated rule (slice 14's follow-up 3): both tails' shares above α/2.
            if (Margin?.WeightedInflatedLower is not { } lower) throw new InvalidOperationException("CDF-score inflation requires the full-step options overload.");
            double half = 0.5 * Checked(alpha);
            return lower[g] > half && Margin.WeightedInflatedUpper![g] > half;
        }
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

    /// <summary>
    /// The predictive's mass on each of the grid's G + 1 cells, (−∞, g₀], (g₀, g₁], …, (g_{G−1}, ∞): Q̄ at g₀, the rises of Q̄ between grid
    /// values, and 1 − Q̄ at the last (slice 18). A law on the cells, which the proper scores read. Null for a case with no distribution.
    /// </summary>
    public double[]? CellMasses(CellLaw law = CellLaw.Steps)
    {
        if (law == CellLaw.Smoothed) return SmoothedMasses();
        if (!Enum.IsDefined(law)) throw new ArgumentOutOfRangeException(nameof(law));
        if (!Monotone) return null;
        var masses = new double[Grid.Length + 1];
        double below = 0.0;
        for (int g = 0; g < Grid.Length; g++)
        {
            double at = Mid(g);
            masses[g] = at - below;
            below = at;
        }
        masses[^1] = 1.0 - below;
        return masses;
    }

    /// <summary>
    /// The log score at an observed <paramref name="y"/>, −ln p of the cell holding it (Gressmann et al.'s log-loss on the cells): strictly
    /// proper over the cells for laws positive where the truth is (<c>logLoss_proper</c>, <c>logLoss_strict</c>; Gneiting and Raftery's
    /// Example 3). +∞ where the cell has no mass; NaN for a case with no distribution.
    /// </summary>
    public double LogScore(double y, CellLaw law = CellLaw.Steps) => CellMasses(law) is { } p ? CellScores.Log(p, Cell(y)) : double.NaN;

    /// <summary>
    /// The quadratic score at an observed <paramref name="y"/>, Σ p² − 2p of the cell holding it (minus Gneiting and Raftery's quadratic
    /// score, their (18)): strictly proper over the cells, its gap to the truth's Σ (p − q)² (<c>gneitingLoss_gap</c>). NaN for a case with
    /// no distribution.
    /// </summary>
    public double QuadraticScore(double y, CellLaw law = CellLaw.Steps) => CellMasses(law) is { } p ? CellScores.Quadratic(p, Cell(y)) : double.NaN;

    // The smoothed law's mass on each cell: each rise read on the tail that holds both ends' digits; null without the smoothed CDF.
    private double[]? SmoothedMasses()
    {
        if (Smoothed is not { } s) return null;
        int last = Grid.Length - 1;
        var masses = new double[Grid.Length + 1];
        masses[0] = s.Lower[0];
        for (int g = 1; g <= last; g++)
            masses[g] = Math.Max(0.0, s.Lower[g] <= 0.5 ? s.Lower[g] - s.Lower[g - 1] : s.Upper[g - 1] - s.Upper[g]);
        masses[^1] = s.Upper[last];
        return masses;
    }

    // The cell holding y: the first grid value at or above it, or the cell above the grid.
    private int Cell(double y)
    {
        if (double.IsNaN(y)) throw new ArgumentException("The observation is NaN.", nameof(y));
        int lo = 0, hi = Grid.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (Grid[mid] >= y) hi = mid; else lo = mid + 1;
        }
        return lo;
    }

    // Q⁺ − Q⁻ = (T + 1)/m, or the weighted ranks' difference.
    private double Spread(int g) => Weighted is { } w ? w.Upper[g] - w.Lower[g] : (Tied[g] + 1.0) / Units;

    // ∫ u² + d²/12 over [a, b], u and d linear between their end values: (b − a)(u_a² + u_a·u_b + u_b²)/3, and d's alike.
    private static double Piece(double a, double b, double ua, double ub, double da, double db)
        => (b - a) * ((ua * ua + ua * ub + ub * ub) / 3.0 + (da * da + da * db + db * db) / 36.0);

    private static void CheckAlpha(double alpha)
    {
        if (!(alpha > 0.0 && alpha < 1.0)) throw new ArgumentOutOfRangeException(nameof(alpha), "α must lie in (0, 1).");
    }

    private static double Checked(double alpha)
    {
        CheckAlpha(alpha);
        return alpha;
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

/// <summary>Proper scores of a law on finitely many outcomes at an observed outcome (slice 18; Gressmann et al., Gneiting and Raftery).</summary>
internal static class CellScores
{
    /// <summary>The log-loss, −ln p_y; +∞ where p_y is 0.</summary>
    public static double Log(double[] p, int outcome) => p[outcome] > 0.0 ? -Math.Log(p[outcome]) : double.PositiveInfinity;

    /// <summary>The quadratic loss, Σ p² − 2p_y.</summary>
    public static double Quadratic(double[] p, int outcome)
    {
        double squares = 0.0;
        foreach (double v in p) squares += v * v;
        return squares - 2.0 * p[outcome];
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

    /// <summary>
    /// The weighted ranks (<see cref="WeightedRanks"/>) from each grid value's unit weight below the candidate, <paramref name="below"/>,
    /// and at or below it, <paramref name="atOrBelow"/>, the stratum's unit total <paramref name="units"/>, and the candidate weights from
    /// column <paramref name="offset"/>: Q⁻ = A/(W + w_c) and Q⁺ = (A + B)/(W + w_c), each read from its own sum, so sums that grow along
    /// the grid under a constant w_c give shares that do not fall.
    /// </summary>
    public static WeightedRanks Ranks(double[] below, double[] atOrBelow, double units, IReadOnlyList<double> candidates, int offset)
    {
        int count = below.Length;
        var lower = new double[count];
        var upper = new double[count];
        double fall = 0.0;
        for (int g = 0; g < count; g++)
        {
            double own = candidates[offset + g], total = units + own;
            lower[g] = total > 0.0 ? below[g] / total : 0.0;
            upper[g] = total > 0.0 ? (atOrBelow[g] + own) / total : 1.0;
            if (g > 0) fall = Math.Max(fall, Math.Max(lower[g - 1] - lower[g], upper[g - 1] - upper[g]));
        }
        return new WeightedRanks(lower, upper, fall);
    }

    /// <summary>
    /// The smoothed CDF G̃ at <paramref name="own"/>, the fit's CDF at a candidate, against the calibration PITs <paramref name="sorted"/>
    /// (<see cref="SmoothedCdf"/>): with B PITs below it and T tied, (B + T)/m at a tie; otherwise (B + f)/m with f the candidate's place
    /// between its neighbouring knots, PIT₍B₎ (or 0) and PIT₍B+1₎ (or 1), read on the scale that holds all three's digits, and the
    /// complement (m − 1 − B + (1 − f))/m alike.
    /// </summary>
    public static (double Lower, double Upper) Smooth(CdfScore[] sorted, CdfScore own)
    {
        int n = sorted.Length, m = n + 1;
        (int below, int tied) = Rank(sorted, own);
        if (tied > 0) return ((double)(below + tied) / m, (double)(m - below - tied) / m);
        double ua = below == 0 ? 0.0 : sorted[below - 1].Lower, ub = below == n ? 1.0 : sorted[below].Lower;
        double va = below == 0 ? 1.0 : sorted[below - 1].Upper, vb = below == n ? 0.0 : sorted[below].Upper;
        double f, rest;
        if (ub <= 0.5)
        {
            double span = ub - ua;
            f = span > 0.0 ? (own.Lower - ua) / span : 0.5;
            rest = span > 0.0 ? (ub - own.Lower) / span : 0.5;
        }
        else
        {
            double span = va - vb;
            f = span > 0.0 ? (va - own.Upper) / span : 0.5;
            rest = span > 0.0 ? (own.Upper - vb) / span : 0.5;
        }
        f = Math.Clamp(f, 0.0, 1.0);
        rest = Math.Clamp(rest, 0.0, 1.0);
        return ((below + f) / m, (m - 1 - below + rest) / m);
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
    /// <paramref name="below"/> and <paramref name="tied"/> from <paramref name="at"/>, and its own score into <paramref name="own"/>. With
    /// <paramref name="unitWeights"/> (the members' shift weights, in their order) the weight of the units below and at or below it too, into
    /// <paramref name="weightBelow"/> and <paramref name="weightAtOrBelow"/>, each summed in the members' order.</summary>
    public void Score(double[,] logDensity, double[,] logCdf, int first, int count, int[] below, int[] tied, CdfScore[] own, int at,
                      double[]? unitWeights = null, double[]? weightBelow = null, double[]? weightAtOrBelow = null)
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
        if (unitWeights is not null)
        {
            Array.Clear(weightBelow!, at, count);
            Array.Clear(weightAtOrBelow!, at, count);
        }
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
                if (unitWeights is not null && order <= 0)
                {
                    if (order < 0) weightBelow![at + b] += unitWeights[i];
                    weightAtOrBelow![at + b] += unitWeights[i];
                }
            }
    }
}
