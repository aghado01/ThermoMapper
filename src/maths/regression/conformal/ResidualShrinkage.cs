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

/// <summary>A residual score's grid-free set at one test case (<see cref="ResidualShrinkage.Intervals"/>): [<paramref name="Lower"/>,
/// <paramref name="Upper"/>] holds the exact set, each end within the tolerance outside it; the whole line at k = 0.</summary>
/// <param name="Case">The test case.</param>
/// <param name="Units">m = n + 1.</param>
/// <param name="Level">k = ⌊α·m⌋.</param>
/// <param name="Lower">The lower end.</param>
/// <param name="Upper">The upper end.</param>
public sealed record ResidualInterval(int Case, int Units, int Level, double Lower, double Upper)
{
    /// <summary>Whether the set is the whole line: m too small for the level.</summary>
    public bool WholeSpace => Level == 0;

    /// <summary>CB-1's coverage, at least 1 − k/m whatever the predictor. A full-conformal set reports no training-conditional law.</summary>
    public double Coverage => 1.0 - (double)Level / Units;

    /// <summary>The interval's length: +∞ for the whole line.</summary>
    public double Width => Upper - Lower;

    /// <summary>Whether <paramref name="y"/> is in the set.</summary>
    public bool Contains(double y) => y >= Lower && y <= Upper;
}

/// <summary>The grid-free sets, one per test case in the centres' order, at miscoverage <paramref name="Alpha"/> under
/// <paramref name="Score"/>.</summary>
public sealed record ResidualIntervals(IReadOnlyList<ResidualInterval> Intervals, double Alpha, ResidualScore Score);

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

    /// <summary>
    /// The sets' exact ends, grid-free (their bracketed search, Algorithm 1, §3.4), for test cases at <paramref name="centres"/> f(x) from the
    /// calibration residuals <paramref name="residuals"/>, each end within <paramref name="tolerance"/> in the response's units and rounded
    /// outward, so each interval holds the exact set. Each score's set is an interval: distance to origin's and distance to average's at every
    /// n, each unit's acceptance an interval about one shared point; RoBAS-EB's for n ≥ 4 (their interval theorem), below which a level
    /// above 0 is refused. Their search starts from an approximate maximiser of the p-value; here it starts from the candidate residual
    /// whose own score is 0, â·r̄ over the calibration residuals (0 for distance to origin, r̄ for distance to average), accepted since
    /// the candidate's ties count for it. From there each side doubles a step until a candidate is refused, which happens since the
    /// candidate's score grows as |z| and the others' as |z|/n, then bisects.
    /// </summary>
    public static ResidualIntervals Intervals(IReadOnlyList<double> residuals, IReadOnlyList<double> centres, ResidualScore score = ResidualScore.Shrinkage,
                                              double alpha = 0.1, double tolerance = 1e-8)
    {
        ArgumentNullException.ThrowIfNull(residuals);
        ArgumentNullException.ThrowIfNull(centres);
        SplitConformal.CheckAlpha(alpha);
        ConformalBisection.CheckTolerance(tolerance);
        if (!Enum.IsDefined(score)) throw new ArgumentOutOfRangeException(nameof(score));
        int n = residuals.Count, m = n + 1, level = (int)Math.Floor(alpha * m);
        var r = new double[m];
        for (int i = 0; i < n; i++)
        {
            r[i] = residuals[i];
            if (!double.IsFinite(r[i])) throw new ArgumentException($"Residual {i} is {r[i]}: residuals are finite.", nameof(residuals));
        }
        for (int t = 0; t < centres.Count; t++)
            if (!double.IsFinite(centres[t])) throw new ArgumentException($"Centre {t} is {centres[t]}.", nameof(centres));
        if (score == ResidualScore.Shrinkage && n < 4 && level > 0)
            throw new ArgumentException($"RoBAS-EB's sets are intervals for n ≥ 4 (their interval theorem); at n = {n} read the grid form (Prediction).", nameof(residuals));
        var intervals = new ResidualInterval[centres.Count];
        if (level == 0)
        {
            for (int t = 0; t < centres.Count; t++) intervals[t] = new ResidualInterval(t, m, 0, double.NegativeInfinity, double.PositiveInfinity);
            return new ResidualIntervals(intervals, alpha, score);
        }

        // The candidate's own score is 0 at â·r̄, its others the calibration residuals alone.
        double mean = 0.0;
        for (int i = 0; i < n; i++) mean += r[i];
        mean /= n;
        double spread = 0.0, reach = 0.0;
        for (int i = 0; i < n; i++) spread += (r[i] - mean) * (r[i] - mean);
        double sigma2 = spread / n, upsilon2 = Math.Max(mean * mean - sigma2 / n, 0.0), denominator = upsilon2 + sigma2 / n;
        double centre = score switch
        {
            ResidualScore.DistanceToOrigin => 0.0,
            ResidualScore.DistanceToAverage => mean,
            _ => (denominator > 0.0 ? upsilon2 / denominator : 0.0) * mean,
        };
        for (int i = 0; i < n; i++) reach = Math.Max(reach, Math.Abs(r[i] - centre));
        var scores = new double[m];
        bool Accepts(double z)
        {
            r[n] = z;
            Score(r, score, scores);
            int atOrAbove = 0;
            for (int i = 0; i < n; i++) if (scores[i] >= scores[n]) atOrAbove++;
            return atOrAbove + 1 > level;
        }
        if (!Accepts(centre))
            throw new InvalidOperationException("The candidate whose own score is 0 is refused: the calibration scores tie at 0 past the level.");
        double Reach(int sign)
        {
            double step = Math.Max(reach, tolerance), inner = centre;
            for (int doubling = 0; doubling < 2100; doubling++, step *= 2.0)
            {
                double z = centre + sign * step;
                if (!double.IsFinite(z)) break;
                if (!Accepts(z))
                {
                    GridInterval crossing = sign > 0 ? ConformalBisection.Narrow(inner, z, true, Accepts, tolerance)
                                                     : ConformalBisection.Narrow(z, inner, false, Accepts, tolerance);
                    return sign > 0 ? crossing.Upper : crossing.Lower;
                }
                inner = z;
            }
            throw new InvalidOperationException("No refused candidate within the doubles' range.");
        }
        double upper = Reach(+1), lower = Reach(-1);
        for (int t = 0; t < centres.Count; t++) intervals[t] = new ResidualInterval(t, m, level, centres[t] + lower, centres[t] + upper);
        return new ResidualIntervals(intervals, alpha, score);
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
