using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>The residual scores <see cref="ResidualShrinkage"/> offers, each a nonconformity score of a residual against the other n.</summary>
public enum ResidualScore
{
    /// <summary>Distance to origin, |r|: trusts the predictor unbiased, efficient when it is, wide under a mean shift.</summary>
    DistanceToOrigin,

    /// <summary>Distance to average, |r − r̄| over the other residuals: the Bayes score of a flat prior on the residuals' mean (their
    /// proposition), robust to a mean shift at the price of estimating it.</summary>
    DistanceToAverage,

    /// <summary>RoBAS-EB, |r − â·r̄| with â = υ̂²/(υ̂² + σ̂²/n), σ̂² the others' mean squared deviation and υ̂² = max(r̄² − σ̂²/n, 0) (their
    /// empirical Bayes estimates): shrinks to the origin when the residuals sit near zero and to their average under a mean shift.</summary>
    Shrinkage,
}

/// <summary>
/// Full conformal prediction on a fixed predictor's residuals (the conformal layer's slice 17; Ashouritaklimi, Cortinovis and Caron,
/// RoBAS, arXiv 2607.04236): a predictor f trained apart from the calibration and test units, residuals rᵢ = yᵢ − f(xᵢ), and for a
/// candidate y at a test case the candidate residual y − f(x). Each of the n + 1 residuals is scored against the other n
/// (<see cref="ResidualScore"/>), and CB-1's rule ranks the candidate: accepted where (B + T + 1)/(n + 1) &gt; ⌊α(n + 1)⌋/(n + 1), B and T
/// the residuals scoring above and level with it. Every score is a function of a residual and the multiset of the others, so the
/// exchangeability of the n + 1 units keeps coverage at least 1 − α whatever the predictor (CB-1). RoBAS-EB's sets are intervals for
/// n ≥ 4 (their interval theorem); read here on each case's grid.
/// </summary>
public static class ResidualShrinkage
{
    /// <summary>
    /// The sets for test cases at <paramref name="centres"/> f(x) on <paramref name="grids"/> of candidate y, from the calibration residuals
    /// <paramref name="residuals"/> y − f(x) under <paramref name="score"/>.
    /// </summary>
    public static ConformalPrediction Prediction(IReadOnlyList<double> residuals, IReadOnlyList<double> centres, IReadOnlyList<double[]> grids,
                                                 ResidualScore score = ResidualScore.Shrinkage, double alpha = 0.1)
    {
        ArgumentNullException.ThrowIfNull(residuals);
        ArgumentNullException.ThrowIfNull(centres);
        ArgumentNullException.ThrowIfNull(grids);
        SplitConformal.CheckAlpha(alpha);
        if (!Enum.IsDefined(score)) throw new ArgumentOutOfRangeException(nameof(score));
        if (centres.Count != grids.Count) throw new ArgumentException("One centre per grid is required.", nameof(centres));
        int n = residuals.Count;
        var r = new double[n + 1];
        for (int i = 0; i < n; i++)
        {
            r[i] = residuals[i];
            if (!double.IsFinite(r[i])) throw new ArgumentException($"Residual {i} is {r[i]}: residuals are finite.", nameof(residuals));
        }
        int m = n + 1, level = (int)Math.Floor(alpha * m);
        var sets = new ConformalSet[grids.Count];
        var scores = new double[m];
        for (int t = 0; t < grids.Count; t++)
        {
            double[] grid = grids[t] ?? throw new ArgumentException($"Grid {t} is null.", nameof(grids));
            if (!double.IsFinite(centres[t])) throw new ArgumentException($"Centre {t} is {centres[t]}.", nameof(centres));
            for (int g = 0; g < grid.Length; g++)
                if (!double.IsFinite(grid[g]) || (g > 0 && !(grid[g] > grid[g - 1])))
                    throw new ArgumentException($"Grid {t}'s values must be finite and strictly ascending.", nameof(grids));
            if (level == 0)
            {
                sets[t] = new ConformalSet(t, 0, m, level, (double[])grid.Clone(), null, null, null, null, null);
                continue;
            }
            var below = new int[grid.Length];
            var tied = new int[grid.Length];
            var p = new double[grid.Length];
            var own = new double[grid.Length];
            for (int g = 0; g < grid.Length; g++)
            {
                r[n] = grid[g] - centres[t];
                Score(r, score, scores);
                double candidate = scores[n];
                own[g] = candidate;
                for (int i = 0; i < n; i++)
                {
                    // Nonconformity: a unit scoring above the candidate is less conforming, the candidate's conformity above it.
                    if (scores[i] > candidate) below[g]++;
                    else if (scores[i] == candidate) tied[g]++;
                }
                p[g] = (below[g] + tied[g] + 1) / (double)m;
            }
            sets[t] = new ConformalSet(t, 0, m, level, (double[])grid.Clone(), p, below, tied, own, null);
        }
        return new ConformalPrediction(sets, alpha);
    }

    /// <summary>Each residual's score against the other n of <paramref name="r"/> under <paramref name="score"/>, into
    /// <paramref name="into"/>.</summary>
    internal static void Score(double[] r, ResidualScore score, double[] into)
    {
        int m = r.Length, n = m - 1;
        // The others' mean and spread from the residuals centred at their overall mean, c = r − μ with Σc = 0: the others' mean is
        // μ − c_j/n and their mean squared deviation (Σc² − c_j²)/n − (c_j/n)², which keeps its digits under a large mean shift.
        double mu = 0.0;
        foreach (double v in r) mu += v;
        mu /= m;
        double centred = 0.0;
        foreach (double v in r) centred += (v - mu) * (v - mu);
        for (int j = 0; j < m; j++)
        {
            if (score == ResidualScore.DistanceToOrigin || n == 0) { into[j] = Math.Abs(r[j]); continue; }
            double c = r[j] - mu, mean = mu - c / n;
            if (score == ResidualScore.DistanceToAverage) { into[j] = Math.Abs(r[j] - mean); continue; }
            double sigma2 = Math.Max(0.0, (centred - c * c) / n - (c / n) * (c / n));
            double upsilon2 = Math.Max(mean * mean - sigma2 / n, 0.0), denominator = upsilon2 + sigma2 / n;
            double shrink = denominator > 0.0 ? upsilon2 / denominator : 0.0;
            into[j] = Math.Abs(r[j] - shrink * mean);
        }
    }
}
