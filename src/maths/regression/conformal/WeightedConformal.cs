using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>
/// Weights for the weighted conformal rule (the conformal layer's slice 14): <see cref="Units"/> at each training or calibration unit,
/// in the units' order, and <see cref="Candidates"/> at each candidate, per column of a candidate seam's grids for the density forms
/// and per test case for the closed-form intervals. A candidate is refused when the weight of the units scoring at most its score, its own
/// weight included, is at most α of the stratum's total with its own (<c>weightLe</c>). What the weights state decides the guarantee:
/// <list type="bullet">
/// <item><b>Covariate shift</b> (Tibshirani, Barber, Candès and Ramdas 2019, Corollary 1): w(x) = dP̃_X/dP_X, the test law's density of
/// the inputs against the training law's, at each unit's and each case's inputs, constant over a case's grid (<see cref="ByCase"/>). With
/// the true ratio the sets cover at least 1 − α under the shifted law (<c>coverage_weighted</c>).</item>
/// <item><b>Response shift</b> (Podkopaev and Ramdas 2021, Theorem 2): w(y), the response's density ratio, the inputs' law given the
/// response unchanged: each unit's at its response and each candidate's at its value, so it moves along a grid. The same guarantee. For a
/// response with finitely many values, ranking within each value's stratum covers under any shift of their proportions with no weights
/// (<c>coverage_stratum</c>).</item>
/// <item><b>An estimated ratio</b> p̂ in place of the true one: the rule refuses at most α plus the total variation between the data set's
/// normalized true and estimated weights, so coverage loses at most its expectation (<c>sum_weight_refused_estimated_le</c>, the layer's
/// statement; the papers' is asymptotic).</item>
/// <item><b>Beyond exchangeability</b> (Barber, Candès, Ramdas and Tibshirani 2023, Theorem 2): weights fixed by position, nonnegative,
/// the candidate's the largest (their test weight 1), for data that drift. A refusal costs at most α + Σ_K q_K·TV(P, P∘swap(K, test)), the
/// distance the swap of unit K with the test unit moves the data's law (<c>refused_le_swap_tv</c>), which the layer cannot read.</item>
/// </list>
/// Every weight is finite and nonnegative. Where a stratum's unit weights and a case's candidate weights all equal one value the rule is
/// CB-1's count, and the layer returns CB-1's set, bit for bit (<c>weighted_eq_cb1</c>).
/// </summary>
public sealed record ShiftWeights
{
    public ShiftWeights(IReadOnlyList<double> units, IReadOnlyList<double> candidates)
    {
        ArgumentNullException.ThrowIfNull(units);
        ArgumentNullException.ThrowIfNull(candidates);
        Units = Checked(units, nameof(units));
        Candidates = Checked(candidates, nameof(candidates));
    }

    /// <summary>Each training or calibration unit's weight, in the units' order.</summary>
    public IReadOnlyList<double> Units { get; }

    /// <summary>Each candidate's weight: per grid column for the density forms, per test case for the closed-form intervals.</summary>
    public IReadOnlyList<double> Candidates { get; }

    /// <summary>Weights given per test case, each spread over its grid's columns: covariate shift's, w(x) constant along a case's grid.</summary>
    public static ShiftWeights ByCase(IReadOnlyList<double> units, IReadOnlyList<double> cases, IReadOnlyList<double[]> grids)
    {
        ArgumentNullException.ThrowIfNull(cases);
        ArgumentNullException.ThrowIfNull(grids);
        if (cases.Count != grids.Count) throw new ArgumentException($"One weight per case is required: {cases.Count} weights for {grids.Count} cases.", nameof(cases));
        var columns = new List<double>();
        for (int t = 0; t < grids.Count; t++)
        {
            double[] grid = grids[t] ?? throw new ArgumentException($"Case {t}'s grid is null.", nameof(grids));
            for (int g = 0; g < grid.Length; g++) columns.Add(cases[t]);
        }
        return new ShiftWeights(units, columns);
    }

    private static double[] Checked(IReadOnlyList<double> weights, string name)
    {
        var copy = new double[weights.Count];
        for (int i = 0; i < copy.Length; i++)
        {
            double w = weights[i];
            if (!(w >= 0.0 && double.IsFinite(w))) throw new ArgumentOutOfRangeException(name, $"Weight {i} is {w}: each weight is finite and nonnegative.");
            copy[i] = w;
        }
        return copy;
    }
}

/// <summary>
/// The weighted rule's calibration (slice 14, <c>weightLe</c>): each stratum's calibration scores ascending, each with its weight, and the
/// running weight of the scores at or below each. A candidate of weight w_c scoring s takes π = (Σ_{sᵢ ≤ s} wᵢ + w_c)/(W_h + w_c), ties
/// counted for it, and is accepted when π &gt; α; its least accepted score is the first sorted score whose running weight with w_c passes α
/// of the total, −∞ where w_c alone does. Where the stratum's weights and w_c are equal it is CB-1's count, read as counts, so it is CB-1's
/// set bit for bit (<c>weighted_eq_cb1</c>). Immutable.
/// </summary>
internal sealed class WeightedCalibration
{
    private readonly double[][] _sorted, _running;   // each stratum's scores ascending, and the weight at or below each
    private readonly double[] _total;
    private readonly double[] _common;   // each stratum's one weight where all are equal, NaN otherwise
    private readonly SplitCalibration _counts;

    internal WeightedCalibration(IReadOnlyList<double> scores, IReadOnlyList<double> weights, int[]? stratumOf, double alpha)
    {
        if (weights.Count != scores.Count)
            throw new ArgumentException($"One weight per calibration unit is required: {weights.Count} weights for {scores.Count} units.", nameof(weights));
        _counts = new SplitCalibration(scores, stratumOf, alpha);
        Alpha = alpha;
        int strata = _counts.Strata;
        var order = new List<int>[strata];
        for (int h = 0; h < strata; h++) order[h] = [];
        for (int i = 0; i < scores.Count; i++) order[stratumOf?[i] ?? 0].Add(i);
        _sorted = new double[strata][];
        _running = new double[strata][];
        _total = new double[strata];
        _common = new double[strata];
        for (int h = 0; h < strata; h++)
        {
            int[] units = [.. order[h]];
            double[] keys = [.. Array.ConvertAll(units, i => scores[i])];
            Array.Sort(keys, units);
            _sorted[h] = keys;
            var running = new double[units.Length];
            double sum = 0.0, common = units.Length > 0 ? weights[units[0]] : double.NaN;
            for (int j = 0; j < units.Length; j++)
            {
                double w = weights[units[j]];
                sum += w;
                running[j] = sum;
                if (w != common) common = double.NaN;
            }
            // Ties take the weight of every unit at their score.
            for (int j = units.Length - 2; j >= 0; j--)
                if (keys[j] == keys[j + 1]) running[j] = running[j + 1];
            _running[h] = running;
            _total[h] = sum;
            _common[h] = common;
        }
    }

    /// <summary>The miscoverage α.</summary>
    public double Alpha { get; }

    /// <summary>The unweighted calibration of the same scores, CB-1's, which equal weights read.</summary>
    public SplitCalibration Counts => _counts;

    /// <summary>Whether the stratum's weights and the candidate's <paramref name="weight"/> are all one value, so the rule is CB-1's count. A
    /// stratum no calibration unit holds is CB-1's at m = 1.</summary>
    public bool Equal(int stratum, double weight)
        => stratum >= _sorted.Length || _sorted[stratum].Length == 0 || _common[stratum] == weight;

    /// <summary>π = (Σ_{sᵢ ≤ s} wᵢ + w_c)/(W_h + w_c), ties counted for the candidate; 1 where the total weight is 0.</summary>
    public double PValue(int stratum, double score, double weight)
    {
        if (double.IsNaN(score)) throw new ArgumentException("The score is NaN.", nameof(score));
        if (stratum >= _sorted.Length) return 1.0;
        double total = _total[stratum] + weight;
        if (!(total > 0.0)) return 1.0;
        int at = AtOrBelow(_sorted[stratum], score);
        double below = at == 0 ? 0.0 : _running[stratum][at - 1];
        return (below + weight) / total;
    }

    /// <summary>The least score accepted at candidate weight <paramref name="weight"/>: the first sorted score whose weight at or below,
    /// with the candidate's, passes α of the total; −∞ where the candidate's weight alone does.</summary>
    public double Threshold(int stratum, double weight)
    {
        if (stratum >= _sorted.Length) return double.NegativeInfinity;
        double total = _total[stratum] + weight;
        if (!(total > 0.0) || weight / total > Alpha) return double.NegativeInfinity;
        double[] running = _running[stratum];
        for (int j = 0; j < running.Length; j++)
            if ((running[j] + weight) / total > Alpha) return _sorted[stratum][j];
        return double.PositiveInfinity;
    }

    // The sorted scores at or below the score.
    private static int AtOrBelow(double[] sorted, double score)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid] <= score) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
}
