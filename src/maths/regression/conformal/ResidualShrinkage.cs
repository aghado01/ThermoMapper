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

    /// <summary>
    /// RoBAS-Full, the horseshoe working model's Bayes-assisted score (their (8) and Appendix's derivation): R | θ ~ N(θ, σ²),
    /// θ | τ² ~ N(0, σ²τ²/n), τ ~ C⁺(0, 1). By their marginal-likelihood lemma a unit's score is the marginal likelihood of the other n,
    /// ln p(r₋ⱼ) = −n s²₋ⱼ/(2σ̂²) + ln ₁F₁(1; 3/2; −n r̄²₋ⱼ/(2σ̂²)) up to a constant every unit shares, s²₋ⱼ and r̄₋ⱼ the others' spread
    /// and mean and σ̂² = Σ(r − r̄)²/n over all n + 1 residuals (their plug-in). Near zero it shrinks as the horseshoe does; under a large
    /// mean shift it ranks as distance to average (their asymptotic robustness theorem). Its sets are intervals from n 4 when the
    /// calibration residuals are not all equal (the Lean's <c>robasFull_ordConnected</c>, with no hypothesis left: the four bounds on
    /// Dawson's function it rests on are proved), so the grid-free form serves it there (<see cref="ResidualShrinkage.Intervals"/>). Below
    /// n 4 they are not always: a seeded search found sets in two pieces at n 2 and 3, one about the prior's centre 0 and one about the
    /// data (residuals 0.510, 0.653, 0.546 at k 1: [−0.198, 0.062] and [0.333, 0.668]), and there the grid form serves it
    /// (<see cref="ResidualShrinkage.Prediction"/>).
    /// </summary>
    Horseshoe,
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
    /// n, each unit's acceptance an interval about one shared point; RoBAS-EB's for n ≥ 4 (their interval theorem); RoBAS-Full's for n ≥ 4
    /// with the calibration residuals not all equal (the Lean's <c>robasFull_ordConnected</c>, with no hypothesis left: the four bounds on
    /// Dawson's function it rests on are proved). Below those, a level above 0 is refused. Their search starts from an approximate
    /// maximiser of the p-value; here it starts from a candidate every unit's score reaches, so the p-value is 1: the residual whose own
    /// score is 0, â·r̄ over the calibration residuals (0 for distance to origin, r̄ for distance to average), the candidate's ties
    /// counting for it; for RoBAS-Full the common point (<c>robas_common_point</c>). From there each side doubles a step until a candidate
    /// is refused, which happens since the candidate's score outgrows the others', then bisects.
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
        if (score == ResidualScore.Horseshoe && level > 0)
        {
            if (n < 4)
                throw new ArgumentException($"RoBAS-Full's sets are not always intervals below n 4 (two pieces at n 2 and 3); at n = {n} read the grid form (Prediction).", nameof(residuals));
            if (Array.TrueForAll(r[..n], v => v == r[0]))
                throw new ArgumentException("RoBAS-Full's interval theorem takes calibration residuals not all equal; read the grid form (Prediction).", nameof(residuals));
        }
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
            ResidualScore.Horseshoe => CommonPoint(r[..n]),
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
            if (score == ResidualScore.Horseshoe)
            {
                // The marginal likelihood of the others at the plug-in σ̂², shared by every unit; all residuals equal, every unit ties.
                double plug = centred / n;
                into[j] = plug > 0.0 ? -n * sigma2 / (2.0 * plug) + LogKummer(n * mean * mean / (2.0 * plug)) : 0.0;
                continue;
            }
            double upsilon2 = Math.Max(mean * mean - sigma2 / n, 0.0), denominator = upsilon2 + sigma2 / n;
            double shrink = denominator > 0.0 ? upsilon2 / denominator : 0.0;
            into[j] = Math.Abs(r[j] - shrink * mean);
        }
    }

    /// <summary>
    /// RoBAS-Full's common point (the Lean's <c>robas_common_point</c>): z₀ with n·z₀ = S·g(S²/(2T(z₀))), S the calibration residuals' sum,
    /// T(z) = Σ(r − r̄)² over them and z, and g = K₁/K₀ with K_j(x) = ∫₀¹ v^{2j}e^{xv²} dv. Every unit's score reaches the candidate's
    /// there, so its p-value is 1. Up to x = 1 g is the ratio of the series K_j(x) = Σ xᵏ/(k!(2k + 2j + 1)); past it, by the ODE
    /// 2xK₁ + K₀ = eˣ, g(x) = (1/₁F₁(1; 3/2; −x) − 1)/(2x), where the difference is above 0.9. g is nondecreasing from 1/3 to 1, so z₀ lies
    /// between S/(3n) and S/n, where n·z − S·g changes sign; bisection finds it. The residuals are not all equal.
    /// </summary>
    internal static double CommonPoint(double[] calibration)
    {
        int n = calibration.Length;
        double s = 0.0;
        foreach (double v in calibration) s += v;
        if (s == 0.0) return 0.0;
        var all = new double[n + 1];
        Array.Copy(calibration, all, n);
        double H(double z)
        {
            all[n] = z;
            double mean = 0.0, spread = 0.0;
            foreach (double v in all) mean += v;
            mean /= n + 1;
            foreach (double v in all) spread += (v - mean) * (v - mean);
            return n * z - s * G(s * s / (2.0 * spread));
        }
        static double G(double x)
        {
            if (x > 1.0) return (Math.Exp(-LogKummer(x)) - 1.0) / (2.0 * x);
            double term = 1.0, k0 = 1.0, k1 = 1.0 / 3.0;
            for (int k = 1; k < 100 && term > 1e-18; k++)
            {
                term *= x / k;
                k0 += term / (2 * k + 1);
                k1 += term / (2 * k + 3);
            }
            return k1 / k0;
        }
        double lo = s / (3.0 * n), hi = s / n, low = H(lo);
        for (int step = 0; step < 200 && hi != lo; step++)
        {
            double mid = 0.5 * (lo + hi);
            if (mid == lo || mid == hi) break;
            if (Math.Sign(H(mid)) == Math.Sign(low)) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /// <summary>
    /// ln ₁F₁(1; 3/2; −x) for x ≥ 0, which is D(√x)/√x, D Dawson's integral (Kummer's transformation). Up to √x = 6 it sums
    /// e⁻ˣ Σ xᵏ/(k!(2k + 1)), every term positive; past it, the asymptotic (1/z) Σ (2k − 1)!!/zᵏ with z = 2x, cut at its least term, which
    /// at the seam is about e⁻³⁶. Both keep about 14 digits.
    /// </summary>
    internal static double LogKummer(double x)
    {
        if (!(x >= 0.0)) throw new ArgumentOutOfRangeException(nameof(x), "x is nonnegative.");
        if (x == 0.0) return 0.0;
        if (double.IsPositiveInfinity(x)) return double.NegativeInfinity;
        if (x <= 36.0)
        {
            double term = 1.0, sum = 1.0;
            for (int k = 1; k < 1000; k++)
            {
                term *= x / k;
                double add = term / (2 * k + 1);
                sum += add;
                if (add <= sum * 1e-17) break;
            }
            return -x + Math.Log(sum);
        }
        double z = 2.0 * x, a = 1.0, series = 1.0;
        for (int k = 1; k < 1000; k++)
        {
            double next = a * (2 * k - 1) / z;
            if (next >= a) break;
            a = next;
            series += a;
            if (a <= series * 1e-17) break;
        }
        return Math.Log(series) - Math.Log(z);
    }
}
