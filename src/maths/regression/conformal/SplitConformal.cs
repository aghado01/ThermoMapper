using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace Maths.Regression.Conformal;

/// <summary>Which form a conformal prediction came from.</summary>
public enum ConformalForm
{
    /// <summary>The full step (<see cref="ConformalBayes.FullStep"/>): every unit scored in sample, from the retained draws.</summary>
    Full,

    /// <summary>The split form (<see cref="SplitConformal"/>): held-out calibration units scored by a fit on the training split.</summary>
    Split,
}

/// <summary>
/// The split form's calibration (CB-2): each stratum's calibration conformity scores, ascending, and the level its m takes. A test
/// unit in stratum h is ranked among the stratum's n_h calibration units, m = n_h + 1 with itself, at k = ⌊α·m⌋, ties counted for
/// it (CB-1, <c>card_rankLe_le</c>): accepted when #{calibration scores ≤ its score} + 1 &gt; k, which is its score reaching the k-th
/// smallest calibration score (<see cref="Threshold"/>). Below m = ⌈1/α⌉, k is 0 and every score is accepted. A stratum no
/// calibration unit holds has m = 1. Immutable.
/// </summary>
public sealed class SplitCalibration
{
    private readonly double[][] _sorted;   // each stratum's calibration scores, ascending

    internal SplitCalibration(IReadOnlyList<double> scores, int[]? stratumOf, double alpha)
    {
        ArgumentNullException.ThrowIfNull(scores);
        SplitConformal.CheckAlpha(alpha);
        if (stratumOf is not null && stratumOf.Length != scores.Count)
            throw new ArgumentException("One stratum per calibration unit is required.", nameof(stratumOf));
        int strata = 1;
        if (stratumOf is not null)
            foreach (int h in stratumOf)
            {
                if (h < 0) throw new ArgumentOutOfRangeException(nameof(stratumOf), "Stratum numbers start at 0.");
                strata = Math.Max(strata, h + 1);
            }
        var counts = new int[strata];
        for (int i = 0; i < scores.Count; i++)
        {
            if (double.IsNaN(scores[i])) throw new ArgumentException($"Calibration unit {i}'s score is NaN.", nameof(scores));
            counts[stratumOf?[i] ?? 0]++;
        }
        _sorted = new double[strata][];
        for (int h = 0; h < strata; h++) _sorted[h] = new double[counts[h]];
        Array.Clear(counts);
        for (int i = 0; i < scores.Count; i++)
        {
            int h = stratumOf?[i] ?? 0;
            _sorted[h][counts[h]++] = scores[i];
        }
        foreach (double[] s in _sorted) Array.Sort(s);
        Alpha = alpha;
    }

    /// <summary>The miscoverage α.</summary>
    public double Alpha { get; }

    /// <summary>The strata numbered from 0 below this; any stratum at or past it holds no calibration unit.</summary>
    public int Strata => _sorted.Length;

    /// <summary>m: the stratum's calibration units and the test unit.</summary>
    public int Units(int stratum) => Scores(stratum).Length + 1;

    /// <summary>k = ⌊α·m⌋.</summary>
    public int Level(int stratum) => (int)Math.Floor(Alpha * Units(stratum));

    /// <summary>The stratum's calibration scores, ascending.</summary>
    public ReadOnlySpan<double> Scores(int stratum)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(stratum);
        return stratum < _sorted.Length ? _sorted[stratum] : [];
    }

    /// <summary>
    /// The least score accepted: the k-th smallest calibration score, since a score reaching it has at least k calibration scores at
    /// or below it; −∞ when k is 0.
    /// </summary>
    public double Threshold(int stratum)
    {
        int k = Level(stratum);
        return k == 0 ? double.NegativeInfinity : _sorted[stratum][k - 1];
    }

    /// <summary>The stratum's calibration units scoring strictly below <paramref name="score"/>, and those equal to it.</summary>
    public (int Below, int Tied) Rank(int stratum, double score)
    {
        if (double.IsNaN(score)) throw new ArgumentException("The score is NaN.", nameof(score));
        ReadOnlySpan<double> s = Scores(stratum);
        int below = CountWhere(s, score, strict: true);
        return (below, CountWhere(s, score, strict: false) - below);
    }

    /// <summary>π = (Below + Tied + 1)/m, ties counted for the test unit.</summary>
    public double PValue(int stratum, double score)
    {
        (int below, int tied) = Rank(stratum, score);
        return (below + tied + 1) / (double)Units(stratum);
    }

    /// <summary>Whether a test unit scoring <paramref name="score"/> is accepted: π &gt; k/m, its score at least
    /// <see cref="Threshold"/>.</summary>
    public bool Accepts(int stratum, double score)
    {
        if (double.IsNaN(score)) throw new ArgumentException("The score is NaN.", nameof(score));
        return score >= Threshold(stratum);
    }

    // The sorted scores below the score (strict) or at or below it, by binary search.
    private static int CountWhere(ReadOnlySpan<double> sorted, double score, bool strict)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (strict ? sorted[mid] < score : sorted[mid] <= score) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
}

/// <summary>A fit's two quantile functions read at each case, CQR's q_lo(x) and q_hi(x), with each case's stratum (null: all in
/// stratum 0). Any fit on the training split serves: a quantile fit's surface, or a predictive's quantiles.</summary>
public sealed record QuantileReadings(IReadOnlyList<double> Lower, IReadOnlyList<double> Upper, int[]? StratumOf = null);

/// <summary>A fit's centre μ(x) and scale ρ(x) read at each case, the standardized residual's, with each case's stratum (null: all
/// in stratum 0). A null scale is 1 at every case, the plain absolute residual.</summary>
public sealed record ResidualReadings(IReadOnlyList<double> Center, IReadOnlyList<double>? Scale = null, int[]? StratumOf = null);

/// <summary>
/// One test case's split set when it is an interval in closed form (CQR's score, <c>cqr_ordConnected</c>; the standardized
/// residual's): [<see cref="Lower"/>, <see cref="Upper"/>], closed, its endpoints to rounding. Empty when Lower &gt; Upper, which
/// CQR gives when the test case's band q_hi − q_lo is less than −2Q: a negative Q, when most calibration units fall well inside
/// their bands, or quantiles that cross by more than 2Q. A stratum too small for the level (k = 0) has the whole line.
/// </summary>
/// <param name="Case">The test case, in the readings' order.</param>
/// <param name="Stratum">Its stratum.</param>
/// <param name="Units">m: the stratum's calibration units and the test case.</param>
/// <param name="Level">k = ⌊α·m⌋.</param>
/// <param name="Threshold">The least conformity score accepted, <see cref="SplitCalibration.Threshold"/>; −∞ for the whole
/// line.</param>
/// <param name="Lower">The interval's lower end; −∞ for the whole line.</param>
/// <param name="Upper">Its upper end; +∞ for the whole line.</param>
public sealed record ConformalInterval(int Case, int Stratum, int Units, int Level, double Threshold, double Lower, double Upper)
{
    /// <summary>Whether the stratum is too small for the level, so the set is the whole line.</summary>
    public bool WholeSpace => Level == 0;

    /// <summary>The guaranteed coverage 1 − k/m: at least this for any fit (<c>split_coverage</c>), exactly this with distinct
    /// scores.</summary>
    public double Coverage => 1.0 - (double)Level / Units;

    /// <summary>Whether the set holds nothing.</summary>
    public bool Empty => Lower > Upper;

    /// <summary>The interval's length: +∞ for the whole line, 0 when empty.</summary>
    public double Width => WholeSpace ? double.PositiveInfinity : Empty ? 0.0 : Upper - Lower;

    /// <summary>Whether <paramref name="y"/> is in the set.</summary>
    public bool Contains(double y) => y >= Lower && y <= Upper;
}

/// <summary>The split form's intervals, one per test case in the readings' order, at miscoverage <paramref name="Alpha"/>, with
/// the calibration they were read from.</summary>
public sealed record ConformalIntervals(IReadOnlyList<ConformalInterval> Intervals, double Alpha, SplitCalibration Calibration);

// Keep the fit's draw metadata after calibration, so the sorted PITs do not retain the draw-by-calibration matrices.
internal readonly record struct SplitDraws(int Count, IReadOnlyList<int>? ChainLengths, IReadOnlyList<double>? Weights,
                                           IReadOnlyList<LatentReport> Latents)
{
    internal SplitDraws(CandidateLogDensity calibration)
        : this(calibration.Draws, calibration.ChainLengths is { } c ? [.. c] : null,
               calibration.Weights is { } w ? [.. w] : null, [.. calibration.Latents]) { }

    internal void Check(CandidateLogDensity candidates)
    {
        ArgumentNullException.ThrowIfNull(candidates);
        candidates.CheckServes(ConformalForms.Split, nameof(candidates));
        const string same = "every unit, calibration or test, is read by one function of the unit (CB-2, splitScore_equivariant)";
        if (candidates.Draws != Count)
            throw new ArgumentException($"The candidates hold {candidates.Draws} draws, the calibration {Count}: {same}.", nameof(candidates));
        if (!ConformalBayes.SameSequence(ChainLengths, candidates.ChainLengths))
            throw new ArgumentException($"The candidates' chain lengths are not the calibration's: {same}.", nameof(candidates));
        if (!ConformalBayes.SameSequence(Weights, candidates.Weights))
            throw new ArgumentException($"The candidates' draw weights are not the calibration's: {same}.", nameof(candidates));
        if (!ConformalBayes.SameLatents(Latents, candidates.Latents))
            throw new ArgumentException($"The candidates' latents are not treated as the calibration's are: {same}.", nameof(candidates));
    }
}

/// <summary>
/// The split form (CB-2; E16): conformal sets from a fit on a training split, ranked among held-out calibration units. The fit is
/// any function of the training split, a sampler's seeds and its whole chain included, mixed or not, right model or wrong; the
/// calibration and test units are read by that one function (<c>splitScore_equivariant</c>), and given the training split they
/// are exchangeable when every unit is (<c>exchangeable_given_training</c>). So a test case is covered with chance at least
/// 1 − k/m, exactly 1 − k/m with distinct scores (<c>split_coverage</c>), and the draws' Monte Carlo error costs no coverage: it is
/// part of the fit. The price is the calibration units, which the fit does not see. Scores, higher conforming better:
/// <list type="bullet">
/// <item><b>The predictive density</b> (<see cref="Density"/>): ln p̂(y | x) = ln Σ_s W_s p(y | x, θ_s), the fit's draws' predictive
/// read through the candidate seam (<see cref="CandidateLogDensity"/>), at the calibration units' responses and on each test
/// case's grid. Its set is a union of intervals on the grid.</item>
/// <item><b>CQR's</b> (<see cref="Quantiles"/>; Romano, Patterson &amp; Candès 2019; <c>cqrScore</c>): min(y − q_lo(x),
/// q_hi(x) − y), the negative of their max(q_lo − y, y − q_hi). Its set is the interval [q_lo − Q, q_hi + Q] with Q the
/// negated threshold (<c>cqr_ordConnected</c>).</item>
/// <item><b>The standardized residual</b> (<see cref="Residuals"/>): −|y − μ(x)|/ρ(x). Its set is the interval μ ± Qρ.</item>
/// <item><b>The predictive CDF</b> (<see cref="Distribution"/>; Vovk, Petej, Nouretdinov, Manokhin &amp; Gammerman 2020): F̂(y | x) =
/// Σ_s W̄_s F(y | x, θ_s), the fit's predictive CDF, whose ranks give the conformal predictive distribution
/// (<see cref="ConformalDistribution"/>), the layer's default.</item>
/// </list>
/// Ranks are taken within each test case's stratum, among the calibration units of that stratum (<c>coverage_stratum</c>).
/// Thread-safe: it reads its inputs and writes only its own buffers.
/// </summary>
public static class SplitConformal
{
    /// <summary>The calibration from any conformity scores, each calibration unit's under the fit, with their strata (null: all in
    /// stratum 0).</summary>
    public static SplitCalibration Calibrate(IReadOnlyList<double> scores, int[]? stratumOf = null, double alpha = 0.1)
        => new(scores, stratumOf, alpha);

    /// <summary>
    /// The split form under the predictive density score: <paramref name="calibration"/> is the fit's candidate density at the
    /// calibration cases, each case's grid the one value of its observed response; <paramref name="candidates"/> the same fit's at
    /// the test cases, on their grids. Both must hold the same draws, chain lengths, weights and latent report, so that every unit is
    /// read by one function. Each case's set is π on its grid within its stratum, as the full step reports it
    /// (<see cref="ConformalSet"/>), with no add-one-in weights: no ESS, no margin, no diagnostics.
    /// </summary>
    public static ConformalPrediction Density(CandidateLogDensity calibration, CandidateLogDensity candidates, double alpha = 0.1)
    {
        CheckAlpha(alpha);
        CheckPair(calibration, candidates);

        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(calibration.Weights, calibration.Draws);
        var scores = new double[calibration.Cases];
        LogPredictive(calibration.Values, rows, logWeights, 0, calibration.Cases, scores);
        var fit = new SplitCalibration(scores, calibration.StratumOf, alpha);
        var sets = new ConformalSet[candidates.Cases];
        for (int t = 0; t < candidates.Cases; t++)
        {
            int h = candidates.StratumOf[t], m = fit.Units(h), k = fit.Level(h);
            double[] grid = candidates.Grids[t];
            if (k == 0)
            {
                sets[t] = new ConformalSet(t, h, m, k, grid, null, null, null, null, null);
                continue;
            }
            int g = grid.Length;
            var own = new double[g];
            LogPredictive(candidates.Values, rows, logWeights, candidates.Offset(t), g, own);
            var below = new int[g];
            var tied = new int[g];
            var p = new double[g];
            for (int j = 0; j < g; j++)
            {
                (below[j], tied[j]) = fit.Rank(h, own[j]);
                p[j] = (below[j] + tied[j] + 1) / (double)m;
            }
            sets[t] = new ConformalSet(t, h, m, k, grid, p, below, tied, own, null);
        }
        return new ConformalPrediction(sets, alpha, ConformalForm.Split);
    }

    /// <summary>
    /// The split form's conformal predictive distribution, the layer's default (the user's decision of 2026-11-08; Vovk, Petej,
    /// Nouretdinov, Manokhin and Gammerman 2020). The score is the fit's predictive CDF F̂(y | x) = Σ_s W̄_s F(y | x, θ_s), read by both
    /// tails from the seams' log CDF: <paramref name="calibration"/> at the calibration cases, each on a grid of its observed response,
    /// gives the calibration PITs, and <paramref name="candidates"/> the test cases on their grids. The calibration's scores do not move
    /// with the candidate and the candidate's rises with y, so both of Vovk's conditions hold by construction and every case has a
    /// distribution (its <see cref="ConformalDistribution.Fall"/> is 0, and reads the CDFs' rounding if not). It is the fit's predictive
    /// recalibrated by the calibration PITs, valid for any fit (CB-2), PIT included. Both seams must hold their log CDF, the same
    /// draws, chain lengths, weights and latent report.
    /// </summary>
    public static ConformalDistributions Distribution(CandidateLogDensity calibration, CandidateLogDensity candidates)
        => CalibrateDistribution(calibration).OnGrid(candidates);

    /// <summary>Sort the held-out PITs once per stratum for grid-free quantiles and equi-tailed sets at new cases.</summary>
    public static SplitPitCalibration CalibrateDistribution(CandidateLogDensity calibration) => new(calibration);

    // CB-2's alignment: calibration units each on their response, and every unit read by one function.
    private static void CheckPair(CandidateLogDensity calibration, CandidateLogDensity candidates)
    {
        CheckCalibration(calibration);
        new SplitDraws(calibration).Check(candidates);
    }

    internal static void CheckCalibration(CandidateLogDensity calibration)
    {
        ArgumentNullException.ThrowIfNull(calibration);
        calibration.CheckServes(ConformalForms.Split, nameof(calibration));
        for (int i = 0; i < calibration.Cases; i++)
            if (calibration.Grids[i].Length != 1)
                throw new ArgumentException($"Calibration case {i}'s grid holds {calibration.Grids[i].Length} values: each calibration unit is read at its "
                                            + "observed response, a grid of one value.", nameof(calibration));
    }

    /// <summary>
    /// The split form under CQR's score: the fit's quantiles at the calibration cases with their observed
    /// <paramref name="response"/>, and at the test cases. Each test case's set is the interval [q_lo − Q, q_hi + Q], Q the negated
    /// threshold of its stratum (<c>cqr_ordConnected</c>): the k-th largest of the calibration units' max(q_lo − y, y − q_hi).
    /// </summary>
    public static ConformalIntervals Quantiles(QuantileReadings calibration, IReadOnlyList<double> response, QuantileReadings test, double alpha = 0.1)
    {
        CheckAlpha(alpha);
        int n = Check(calibration, response, nameof(calibration));
        int cases = Check(test, null, nameof(test));
        var scores = new double[n];
        for (int i = 0; i < n; i++) scores[i] = CqrScore(calibration.Lower[i], calibration.Upper[i], response[i]);
        var fit = new SplitCalibration(scores, calibration.StratumOf, alpha);
        var intervals = new ConformalInterval[cases];
        for (int t = 0; t < cases; t++)
        {
            int h = test.StratumOf?[t] ?? 0;
            double c = fit.Threshold(h);
            intervals[t] = fit.Level(h) == 0
                ? WholeLine(t, h, fit)
                : new ConformalInterval(t, h, fit.Units(h), fit.Level(h), c, test.Lower[t] + c, test.Upper[t] - c);
        }
        return new ConformalIntervals(intervals, alpha, fit);
    }

    /// <summary>
    /// The split form under the standardized residual's score: the fit's centre and scale at the calibration cases with their
    /// observed <paramref name="response"/>, and at the test cases. Each test case's set is the interval μ ± Qρ, Q the negated
    /// threshold of its stratum: the k-th largest of the calibration units' |y − μ|/ρ.
    /// </summary>
    public static ConformalIntervals Residuals(ResidualReadings calibration, IReadOnlyList<double> response, ResidualReadings test, double alpha = 0.1)
    {
        CheckAlpha(alpha);
        int n = Check(calibration, response, nameof(calibration));
        int cases = Check(test, null, nameof(test));
        var scores = new double[n];
        for (int i = 0; i < n; i++) scores[i] = ResidualScore(calibration.Center[i], calibration.Scale?[i] ?? 1.0, response[i]);
        var fit = new SplitCalibration(scores, calibration.StratumOf, alpha);
        var intervals = new ConformalInterval[cases];
        for (int t = 0; t < cases; t++)
        {
            int h = test.StratumOf?[t] ?? 0;
            double c = fit.Threshold(h), half = -c * (test.Scale?[t] ?? 1.0);
            intervals[t] = fit.Level(h) == 0
                ? WholeLine(t, h, fit)
                : new ConformalInterval(t, h, fit.Units(h), fit.Level(h), c, test.Center[t] - half, test.Center[t] + half);
        }
        return new ConformalIntervals(intervals, alpha, fit);
    }

    /// <summary>CQR's conformity score at (x, y), the Lean <c>cqrScore</c>: min(y − q_lo, q_hi − y).</summary>
    internal static double CqrScore(double lower, double upper, double y) => Math.Min(y - lower, upper - y);

    /// <summary>The standardized residual's conformity score at (x, y): −|y − μ|/ρ.</summary>
    internal static double ResidualScore(double center, double scale, double y) => -Math.Abs(y - center) / scale;

    /// <summary>
    /// ln Σ_s W_s p(y | θ_s) at columns [<paramref name="first"/>, first + count): each column's terms, with the draws' normalized
    /// log weights when weighted, shifted by their largest, exponentiated and summed in draw order, less ln S at equal weights. It is
    /// the full step's arithmetic for the draws' predictive, bit for bit, so equal columns score equal.
    /// </summary>
    internal static void LogPredictive(double[,] values, int[] rows, double[]? logWeights, int first, int count, Span<double> into)
    {
        if (count == 0) return;
        var top = new double[count];
        var total = new double[count];
        top.AsSpan().Fill(double.NegativeInfinity);
        for (int d = 0; d < rows.Length; d++)
        {
            ReadOnlySpan<double> row = MemoryMarshal.CreateReadOnlySpan(ref values[rows[d], first], count);
            if (logWeights is null)
                for (int b = 0; b < count; b++) top[b] = Math.Max(top[b], row[b]);
            else
                for (int b = 0; b < count; b++) top[b] = Math.Max(top[b], logWeights[d] + row[b]);
        }
        for (int d = 0; d < rows.Length; d++)
        {
            ReadOnlySpan<double> row = MemoryMarshal.CreateReadOnlySpan(ref values[rows[d], first], count);
            double w = logWeights is null ? 0.0 : logWeights[d];
            for (int b = 0; b < count; b++)
                if (!double.IsNegativeInfinity(top[b])) total[b] += Math.Exp((logWeights is null ? row[b] : w + row[b]) - top[b]);
        }
        double shift = logWeights is null ? Math.Log(rows.Length) : 0.0;
        for (int b = 0; b < count; b++)
            into[b] = total[b] > 0.0 ? (logWeights is null ? top[b] + Math.Log(total[b]) - shift : top[b] + Math.Log(total[b])) : double.NegativeInfinity;
    }

    internal static void CheckAlpha(double alpha)
    {
        if (!(alpha > 0.0 && alpha < 1.0)) throw new ArgumentOutOfRangeException(nameof(alpha), "α must lie in (0, 1).");
    }

    private static ConformalInterval WholeLine(int t, int h, SplitCalibration fit)
        => new(t, h, fit.Units(h), 0, double.NegativeInfinity, double.NegativeInfinity, double.PositiveInfinity);

    private static int Check(QuantileReadings readings, IReadOnlyList<double>? response, string name)
    {
        ArgumentNullException.ThrowIfNull(readings, name);
        ArgumentNullException.ThrowIfNull(readings.Lower, name);
        ArgumentNullException.ThrowIfNull(readings.Upper, name);
        int n = readings.Lower.Count;
        if (readings.Upper.Count != n) throw new ArgumentException("One lower and one upper quantile per case are required.", name);
        for (int i = 0; i < n; i++)
            if (!double.IsFinite(readings.Lower[i]) || !double.IsFinite(readings.Upper[i]))
                throw new ArgumentException($"Case {i}'s quantiles must be finite.", name);
        CheckRest(readings.StratumOf, response, n, name);
        return n;
    }

    private static int Check(ResidualReadings readings, IReadOnlyList<double>? response, string name)
    {
        ArgumentNullException.ThrowIfNull(readings, name);
        ArgumentNullException.ThrowIfNull(readings.Center, name);
        int n = readings.Center.Count;
        if (readings.Scale is { } scale && scale.Count != n) throw new ArgumentException("One scale per case is required.", name);
        for (int i = 0; i < n; i++)
        {
            if (!double.IsFinite(readings.Center[i])) throw new ArgumentException($"Case {i}'s centre must be finite.", name);
            if (readings.Scale is { } s && !(s[i] > 0.0 && double.IsFinite(s[i])))
                throw new ArgumentException($"Case {i}'s scale must be positive and finite.", name);
        }
        CheckRest(readings.StratumOf, response, n, name);
        return n;
    }

    private static void CheckRest(int[]? stratumOf, IReadOnlyList<double>? response, int n, string name)
    {
        if (stratumOf is not null)
        {
            if (stratumOf.Length != n) throw new ArgumentException("One stratum per case is required.", name);
            foreach (int h in stratumOf)
                if (h < 0) throw new ArgumentOutOfRangeException(name, "Stratum numbers start at 0.");
        }
        if (response is null) return;
        if (response.Count != n) throw new ArgumentException($"One response per calibration case is required: {response.Count} for {n} cases.", nameof(response));
        for (int i = 0; i < n; i++)
            if (!double.IsFinite(response[i])) throw new ArgumentException($"Calibration case {i}'s response must be finite.", nameof(response));
    }
}
