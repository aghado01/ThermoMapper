using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>Row i contains the fit excluding training unit i: its prediction at that unit and at every test input.</summary>
public sealed record LeaveOneOutReadings(IReadOnlyList<double> Response, IReadOnlyList<double> HeldOut, double[,] Predictions);

/// <summary>A jackknife+ interval, including its attainable finite-sample lower bound (<c>jackknifePlus_coverage</c>).
/// Exchangeable units suffice with symmetric fits, including their randomness, or compact order-reading fits
/// that see no original positions, held-out mask or row-keyed randomness (<c>jackknifePlus_orderReading</c>).</summary>
public sealed record JackknifePlusInterval(int Case, int Units, int Level, double Lower, double Upper)
{
    /// <summary>The weighted rule's α where the interval is the weighted jackknife+'s (slice 14's follow-up 2,
    /// <c>jackknifePlus_likelihoodRatio</c>, <c>jackknifePlus_weighted</c>); null for the counts' interval, equal weights' included.</summary>
    public double? WeightedAlpha { get; init; }

    /// <summary>Whether the interval is the whole line: k = 0 for the counts' interval; both ends unbounded for the weighted one.</summary>
    public bool WholeSpace => WeightedAlpha is null ? Level == 0 : double.IsNegativeInfinity(Lower) && double.IsPositiveInfinity(Upper);
    public bool Empty => Lower > Upper;
    public double Width => Empty ? 0.0 : Upper - Lower;
    public bool Contains(double y) => y >= Lower && y <= Upper;

    /// <summary>max(0, 1 − (2k − 1)/m) for the counts' interval; max(0, 1 − 2α) for the weighted one.</summary>
    public double Coverage => WeightedAlpha is double a ? Math.Max(0.0, 1.0 - 2.0 * a) : CrossConformal.Coverage(Units, Level);
}

public sealed record JackknifePlusIntervals(IReadOnlyList<JackknifePlusInterval> Intervals, double Alpha);

/// <summary>Reduce leave-out readings to conformal sets. The n-fold guarantees
/// (<c>crossFit_coverage</c>, <c>jackknifePlus_coverage</c>) require exchangeability of all n + 1 units.
/// Fits can be symmetric under relabeling, including their randomness, or read compact ordered training rows with the
/// held-out place closed up and no original indices, held-out mask or row-keyed randomness. Ordered sweeps sharing a fixed seed
/// can qualify by this second route. Position-aware fits do not receive that guarantee.
/// <see cref="CvPlus"/> exposes the K-fold bounds under their own assumptions.
/// <para>Each form takes shift weights (<see cref="ShiftWeights"/>, slice 14's follow-up 2): one per training unit, in their order, or per
/// test case, and one per candidate. The candidate is refused when the units it is stranger than carry at least 1 − α of the weight, its
/// own included (<c>strangeWeight</c>); whatever the comparisons, the refused units carry at most 2α of it (<c>sum_weight_strange_le</c>,
/// Lei's lemma). What the weights state decides the guarantee:
/// <list type="bullet">
/// <item><b>Likelihood ratios</b> travelling with the units, the candidate's its own (JAWS; Prinster, Liu and Saria 2022): w(x) for covariate
/// shift, w(y) for response shift, each per candidate. The set misses at most 2α of the weight over the shifted law, with no further term
/// (<c>crossFit_likelihoodRatio</c>, <c>jackknifePlus_likelihoodRatio</c>, <c>tournament_likelihoodRatio</c>).</item>
/// <item><b>Weights fixed by position,</b> the candidate's the largest (Barber, Candès, Ramdas and Tibshirani 2023, Theorem 5 and Appendix D):
/// at most 2α plus Σ_K q_K·TV(P, P∘swap(K, test)), which the layer cannot read; on exchangeable data 2α
/// (<c>crossFit_weighted_swap_tv</c>, <c>crossFit_weighted_exchangeable</c>, <c>jackknifePlus_weighted</c>, <c>tournament_weighted</c>).</item>
/// </list>
/// A case whose weights are all one value takes the counts' set, bit for bit; the others report <c>WeightedAlpha</c> and coverage 1 − 2α.
/// </para></summary>
public static partial class CrossConformal
{
    internal static double Coverage(int units, int level)
        => level == 0 ? 1.0 : Math.Max(0.0, 1.0 - (2.0 * level - 1) / units);

    /// <summary>Closed ends from the k-th smallest μ₋ᵢ(x) − |yᵢ − μ₋ᵢ(xᵢ)| and k-th largest μ₋ᵢ(x) + |yᵢ − μ₋ᵢ(xᵢ)|,
    /// k = ⌊α(n + 1)⌋ (<c>JackRefuses</c>). At k = 0 return the whole line. Rows must come from actual leave-one-out fits.
    /// With <paramref name="weights"/>, one per test case among the candidates, a case whose weights are not all one value takes the
    /// weighted interval (<c>JackRefusesW</c>): y is refused where the upper ends below it, or the lower ends above it, carry at least
    /// 1 − α of the weight, the test's own included, that is Barber et al.'s weighted quantiles with the test's mass at ±∞.</summary>
    public static JackknifePlusIntervals JackknifePlus(LeaveOneOutReadings readings, double alpha = 0.1, ShiftWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(readings);
        ArgumentNullException.ThrowIfNull(readings.Response);
        ArgumentNullException.ThrowIfNull(readings.HeldOut);
        ArgumentNullException.ThrowIfNull(readings.Predictions);
        SplitConformal.CheckAlpha(alpha);
        int n = readings.Response.Count, m = checked(n + 1), k = (int)Math.Floor(alpha * m);
        if (readings.HeldOut.Count != n || readings.Predictions.GetLength(0) != n)
            throw new ArgumentException("One held-out response and one prediction row per leave-one-out fit are required.", nameof(readings));
        double[] residuals = AbsoluteResiduals(readings.Response, readings.HeldOut);
        foreach (double v in readings.Predictions)
            if (!double.IsFinite(v)) throw new ArgumentException("Predictions must be finite.", nameof(readings));
        int cases = readings.Predictions.GetLength(1);
        if (weights is not null) CheckWeights(weights, n, cases, cases, nameof(weights));
        var intervals = new JackknifePlusInterval[cases];
        var lower = new double[n];
        var upper = new double[n];
        for (int t = 0; t < cases; t++)
        {
            bool weighed = weights is not null && !Equal(weights.UnitsFor(t), weights.Candidates, t, 1);
            if (k == 0 && !weighed)
            {
                intervals[t] = new(t, m, k, double.NegativeInfinity, double.PositiveInfinity);
                continue;
            }
            for (int i = 0; i < n; i++)
            {
                lower[i] = readings.Predictions[i, t] - residuals[i];
                upper[i] = readings.Predictions[i, t] + residuals[i];
            }
            if (weighed)
            {
                GridInterval w = WeightedEnds(lower, upper, weights!.UnitsFor(t), weights.Candidates[t], alpha);
                intervals[t] = new(t, m, k, w.Lower, w.Upper) { WeightedAlpha = alpha };
                continue;
            }
            GridInterval ends = OrderedEnds(lower, upper, k);
            intervals[t] = new(t, m, k, ends.Lower, ends.Upper);
        }
        return new(intervals, alpha);
    }

    // The weights serve n training units, `cases` cases and `columns` candidates.
    internal static void CheckWeights(ShiftWeights weights, int units, int cases, int columns, string name)
    {
        if (weights.UnitCount != units)
            throw new ArgumentException($"One weight per training unit is required: {weights.UnitCount} weights for {units} units.", name);
        if (weights.UnitsByCase is { } byCase && byCase.Count != cases)
            throw new ArgumentException($"One set of unit weights per case is required: {byCase.Count} sets for {cases} cases.", name);
        if (weights.Candidates.Count != columns)
            throw new ArgumentException($"One weight per candidate is required: {weights.Candidates.Count} weights for {columns} candidates.", name);
    }

    // Whether a case's unit weights and its candidates' from `first` are all one value, so its set is the counts'.
    internal static bool Equal(IReadOnlyList<double> units, IReadOnlyList<double> candidates, int first, int count)
    {
        double common = units.Count > 0 ? units[0] : candidates[first];
        for (int i = 0; i < units.Count; i++) if (units[i] != common) return false;
        for (int j = 0; j < count; j++) if (candidates[first + j] != common) return false;
        return true;
    }

    /// <summary>
    /// The weighted interval's closed ends (<c>JackRefusesW</c>): the upper end the least upper end whose weight at or below reaches
    /// (1 − α)·W, W the units' weight plus the test's, and the lower end the greatest lower end whose weight at or above reaches it; +∞ or −∞
    /// where none does, the test's own weight being more than α of W. A total weight of 0 refuses nothing: the whole line.
    /// </summary>
    internal static GridInterval WeightedEnds(double[] lower, double[] upper, IReadOnlyList<double> weights, double testWeight, double alpha)
    {
        int n = lower.Length;
        double total = testWeight;
        for (int i = 0; i < n; i++) total += weights[i];
        if (!(total > 0.0)) return new(double.NegativeInfinity, double.PositiveInfinity);
        double need = (1.0 - alpha) * total;
        int[] up = new int[n], down = new int[n];
        for (int i = 0; i < n; i++) up[i] = down[i] = i;
        Array.Sort((double[])upper.Clone(), up);
        Array.Sort(Array.ConvertAll(lower, v => -v), down);
        double hi = double.PositiveInfinity, lo = double.NegativeInfinity, running = 0.0;
        foreach (int i in up)
        {
            running += weights[i];
            if (running >= need) { hi = upper[i]; break; }
        }
        running = 0.0;
        foreach (int i in down)
        {
            running += weights[i];
            if (running >= need) { lo = lower[i]; break; }
        }
        return new(lo, hi);
    }

    // A case's weighted p-values: the weight of the units the candidate is not stranger than, `notStranger` per column, plus its own, over the
    // total; 1 where the total is 0.
    internal static double[] WeightedPValues(double[] notStranger, double units, IReadOnlyList<double> candidates, int first)
    {
        var p = new double[notStranger.Length];
        for (int g = 0; g < p.Length; g++)
        {
            double own = candidates[first + g], total = units + own;
            p[g] = total > 0.0 ? (notStranger[g] + own) / total : 1.0;
        }
        return p;
    }

    private static double[] AbsoluteResiduals(IReadOnlyList<double> response, IReadOnlyList<double> heldOut)
    {
        if (heldOut.Count != response.Count)
            throw new ArgumentException("One held-out prediction per response is required.", nameof(heldOut));
        var residuals = new double[response.Count];
        for (int i = 0; i < residuals.Length; i++)
        {
            if (!double.IsFinite(response[i]) || !double.IsFinite(heldOut[i]))
                throw new ArgumentException("Responses and predictions must be finite.", "readings");
            residuals[i] = Math.Abs(response[i] - heldOut[i]);
        }
        return residuals;
    }

    private static GridInterval OrderedEnds(double[] lower, double[] upper, int level)
    {
        Array.Sort(lower);
        Array.Sort(upper);
        return new(lower[level - 1], upper[upper.Length - level]);
    }

    /// <summary>The n-fold pairwise rule (<c>AcceptsPairwise</c>): candidate score [i,g] and held-out score [i] use the same
    /// fit excluding unit i. Higher scores conform better. A strict candidate &lt; held-out comparison counts as strange, and the
    /// candidate is accepted when the count is below n + 1 − k. Ties count for neither side. With <paramref name="weights"/> (one per
    /// candidate column) a case whose weights are not all one value takes the weighted rule (<c>strangeWeight</c>): π = (weight of the units
    /// the candidate is not stranger than + its own)/W, accepted where π &gt; α.</summary>
    public static ConformalPrediction Pairwise(IReadOnlyList<double> heldOut, double[,] candidateScores,
                                               IReadOnlyList<double[]> grids, double alpha = 0.1, ShiftWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(heldOut);
        ArgumentNullException.ThrowIfNull(candidateScores);
        ArgumentNullException.ThrowIfNull(grids);
        SplitConformal.CheckAlpha(alpha);
        if (grids.Count == 0) throw new ArgumentException("At least one test case is required.", nameof(grids));
        if (candidateScores.GetLength(0) != heldOut.Count)
            throw new ArgumentException("One candidate-score row per leave-one-out fit is required.", nameof(candidateScores));
        foreach (double v in heldOut)
            if (double.IsNaN(v)) throw new ArgumentException("A held-out score is NaN.", nameof(heldOut));
        foreach (double v in candidateScores)
            if (double.IsNaN(v)) throw new ArgumentException("A candidate score is NaN.", nameof(candidateScores));
        int columns = 0;
        var copies = new double[grids.Count][];
        for (int t = 0; t < grids.Count; t++)
        {
            double[] grid = grids[t] ?? throw new ArgumentException("A grid is null.", nameof(grids));
            if (grid.Length == 0) throw new ArgumentException("A grid is empty.", nameof(grids));
            for (int g = 0; g < grid.Length; g++)
                if (!double.IsFinite(grid[g]) || (g > 0 && !(grid[g] > grid[g - 1])))
                    throw new ArgumentException("Grid values must be finite and strictly ascending.", nameof(grids));
            copies[t] = (double[])grid.Clone();
            columns = checked(columns + grid.Length);
        }
        if (columns != candidateScores.GetLength(1))
            throw new ArgumentException("The grids must cover all candidate-score columns.", nameof(candidateScores));
        if (weights is not null) CheckWeights(weights, heldOut.Count, grids.Count, columns, nameof(weights));
        int m = checked(heldOut.Count + 1), k = (int)Math.Floor(alpha * m);
        var sets = new ConformalSet[grids.Count];
        int column = 0;
        for (int t = 0; t < sets.Length; t++)
        {
            double[] grid = copies[t];
            int first = column;
            column += grid.Length;
            IReadOnlyList<double>? unitWeights = weights is not null && !Equal(weights.UnitsFor(t), weights.Candidates, first, grid.Length) ? weights.UnitsFor(t) : null;
            if (k == 0 && unitWeights is null)
            {
                sets[t] = new(t, 0, m, k, grid, null, null, null, null, null, Form: ConformalForm.CrossFitted);
                continue;
            }
            var below = new int[grid.Length];
            var tied = new int[grid.Length];
            var p = new double[grid.Length];
            double[]? notStranger = unitWeights is null ? null : new double[grid.Length];
            double units = 0.0;
            for (int i = 0; i < heldOut.Count; i++)
            {
                if (unitWeights is not null) units += unitWeights[i];
                for (int g = 0; g < grid.Length; g++)
                {
                    double candidate = candidateScores[i, first + g];
                    if (heldOut[i] < candidate) below[g]++;
                    else if (heldOut[i] == candidate) tied[g]++;
                    if (notStranger is not null && heldOut[i] <= candidate) notStranger[g] += unitWeights![i];
                }
            }
            if (notStranger is not null)
            {
                sets[t] = new(t, 0, m, k, grid, WeightedPValues(notStranger, units, weights!.Candidates, first), below, tied, null, null,
                              Form: ConformalForm.CrossFitted) { WeightedAlpha = alpha };
                continue;
            }
            for (int g = 0; g < grid.Length; g++) p[g] = (below[g] + tied[g] + 1) / (double)m;
            sets[t] = new(t, 0, m, k, grid, p, below, tied, null, null, Form: ConformalForm.CrossFitted);
        }
        return new(sets, alpha, ConformalForm.CrossFitted);
    }

    /// <summary>Pairwise sets scored by each leave-one-out fit's predictive density. Each held-out seam holds one response on
    /// one singleton grid. Each candidate seam holds the same test grids. Draws align within each fit and can differ across fits.</summary>
    public static ConformalPrediction Density(IReadOnlyList<CandidateLogDensity> heldOut,
                                              IReadOnlyList<CandidateLogDensity> candidates, double alpha = 0.1, ShiftWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(heldOut);
        ArgumentNullException.ThrowIfNull(candidates);
        SplitConformal.CheckAlpha(alpha);
        if (heldOut.Count == 0 || candidates.Count != heldOut.Count)
            throw new ArgumentException("A nonempty set of paired leave-one-out seams is required.", nameof(heldOut));
        return Density(heldOut.Count, i => (heldOut[i], candidates[i]), alpha, weights);
    }

    /// <summary>Reduce one leave-one-out fit at a time. The callback's matrices can be released after each call, retaining
    /// only the fit-by-candidate scores. Shift weights as for <see cref="Pairwise"/>.</summary>
    public static ConformalPrediction Density(int fits, Func<int, (CandidateLogDensity HeldOut, CandidateLogDensity Candidates)> read,
                                              double alpha = 0.1, ShiftWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(read);
        ArgumentOutOfRangeException.ThrowIfLessThan(fits, 1);
        SplitConformal.CheckAlpha(alpha);
        var own = new double[fits];
        double[,]? scores = null;
        IReadOnlyList<double[]>? grids = null;
        double[]? buffer = null;
        for (int i = 0; i < fits; i++)
        {
            (CandidateLogDensity at, CandidateLogDensity test) = read(i);
            SplitConformal.CheckCalibration(at);
            if (at.Cases != 1) throw new ArgumentException("Each fit must score its one held-out unit.", nameof(read));
            new SplitDraws(at).Check(test);
            if (at.StratumOf[0] != 0 || Array.Exists(test.StratumOf, h => h != 0))
                throw new ArgumentException("Supply one exchangeable group in stratum zero.", nameof(read));
            grids ??= test.Grids;
            scores ??= new double[fits, test.Values.GetLength(1)];
            buffer ??= new double[test.Values.GetLength(1)];
            if (test.Cases != grids.Count) throw new ArgumentException("Every fit must score the same test cases.", nameof(read));
            for (int t = 0; t < grids.Count; t++)
                if (!ConformalBayes.SameSequence(grids[t], test.Grids[t]))
                    throw new ArgumentException("Every fit must score the same test grids.", nameof(read));
            (int[] rows, double[]? drawWeights) = PredictiveChecks.Rows(at.Weights, at.Draws);
            Span<double> value = own.AsSpan(i, 1);
            SplitConformal.LogPredictive(at.Values, rows, drawWeights, 0, 1, value);
            SplitConformal.LogPredictive(test.Values, rows, drawWeights, 0, buffer.Length, buffer);
            for (int g = 0; g < buffer.Length; g++) scores[i, g] = buffer[g];
        }
        return Pairwise(own, scores!, grids!, alpha, weights);
    }
}
