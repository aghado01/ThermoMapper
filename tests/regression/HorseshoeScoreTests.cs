using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// RoBAS-Full, the horseshoe residual score (E16, slice 17's follow-up; Ashouritaklimi, Cortinovis and Caron, (8)): ln ₁F₁(1; 3/2; −x)
/// against quadrature of its integral; the closed form against the working model's marginal likelihood integrated over θ and τ directly;
/// coverage under no shift and a shift, narrower than distance to average unshifted and than distance to origin shifted; under a large
/// shift its sets those of distance to average, as their robustness theorem says; the grid-free form refuses it.
/// </summary>
public sealed class HorseshoeScoreTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // Composite Simpson on [a, b] with `intervals` (even) intervals.
    private static double Simpson(Func<double, double> f, double a, double b, int intervals)
    {
        double h = (b - a) / intervals, sum = f(a) + f(b);
        for (int i = 1; i < intervals; i++) sum += (i % 2 == 1 ? 4.0 : 2.0) * f(a + i * h);
        return sum * h / 3.0;
    }

    /// <summary>ln ₁F₁(1; 3/2; −x) by quadrature of ∫₀¹ e^{−x(1 − v²)} dv: with v = 1 − s², ∫₀¹ 2s·e^{−x s²(2 − s²)} ds, which spreads the
    /// mass piled at v = 1; for x past 100, (1/(2x)) ∫₀⁶⁰ e^{−w}(1 − w/x)^{−1/2} dw.</summary>
    internal static double LogKummerByQuadrature(double x)
        => x <= 100.0
            ? Math.Log(Simpson(s => 2.0 * s * Math.Exp(-x * s * s * (2.0 - s * s)), 0.0, 1.0, 20000))
            : Math.Log(Simpson(w => Math.Exp(-w) / Math.Sqrt(1.0 - w / x), 0.0, 60.0, 20000)) - Math.Log(2.0 * x);

    [Fact]
    public void TheKummerFunction_MatchesItsIntegral()
    {
        double worst = 0.0;
        foreach (double x in new[] { 1e-12, 1e-6, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 35.9, 36.0, 36.1, 50.0, 99.0, 101.0, 1e3, 1e6, 1e12 })
        {
            double got = ResidualShrinkage.LogKummer(x), expected = LogKummerByQuadrature(x);
            worst = Math.Max(worst, Math.Abs(got - expected));
            Assert.True(Math.Abs(got - expected) <= 1e-11, $"x {x}: {got:R} against {expected:R}");
        }
        output.WriteLine($"ln ₁F₁(1; 3/2; −x) against quadrature, largest difference {worst:E2}");
        // Dawson's integral at 1, D(1) = 0.53807950691...: ₁F₁(1; 3/2; −1) = D(1).
        Assert.Equal(0.53807950691, Math.Exp(ResidualShrinkage.LogKummer(1.0)), 10);
        Assert.Equal(0.0, ResidualShrinkage.LogKummer(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => ResidualShrinkage.LogKummer(-1.0));
    }

    [Fact]
    public void TheClosedForm_IsTheWorkingModelsMarginalLikelihood()
    {
        // ln p(r₁…rₙ) under R | θ ~ N(θ, σ²), θ | τ ~ N(0, σ²τ²/n), τ ~ C⁺(0, 1), integrated directly: over φ = atan τ, where the half-Cauchy is
        // uniform with density 2/π, and over θ on its conditional's window, each integrand read from the model's own densities.
        var rng = new Xoshiro256PlusPlus(2026_1149);
        const double sigma2 = 1.3;
        foreach ((int n, double shift) in new[] { (5, 0.0), (12, 2.0), (30, -8.0) })
        {
            double[] r = [.. Enumerable.Range(0, n).Select(_ => shift + 1.1 * Gaussian(rng))];
            double mean = r.Average(), s2 = r.Sum(v => (v - mean) * (v - mean)) / n;
            double closed = -0.5 * n * Math.Log(2.0 * Math.PI * sigma2) - n * s2 / (2.0 * sigma2) + Math.Log(2.0 / Math.PI)
                            + ResidualShrinkage.LogKummer(n * mean * mean / (2.0 * sigma2));
            // The log likelihood's top, at θ = r̄, factored out so the integrand stays in range.
            double top = -0.5 * n * Math.Log(2.0 * Math.PI * sigma2) - n * s2 / (2.0 * sigma2);
            double Inner(double phi)
            {
                double tau = Math.Tan(phi), prior = sigma2 * tau * tau / n;
                double precision = n / sigma2 + 1.0 / prior, centre = (n * mean / sigma2) / precision, sd = Math.Sqrt(1.0 / precision);
                return Simpson(theta =>
                {
                    double log = -0.5 * n * Math.Log(2.0 * Math.PI * sigma2) - r.Sum(v => (v - theta) * (v - theta)) / (2.0 * sigma2)
                                 - 0.5 * Math.Log(2.0 * Math.PI * prior) - theta * theta / (2.0 * prior);
                    return Math.Exp(log - top);
                }, centre - 12.0 * sd, centre + 12.0 * sd, 400);
            }
            double direct = Math.Log((2.0 / Math.PI) * Simpson(Inner, 1e-9, 0.5 * Math.PI - 1e-9, 2000)) + top;
            output.WriteLine($"n {n}, shift {shift}: the closed form {closed:F10}, integrated directly {direct:F10}");
            Assert.True(Math.Abs(closed - direct) <= 1e-6, $"n {n}: {closed} against {direct}");
        }
    }

    [Fact]
    public void TheHorseshoeSets_Cover_AndTrackTheBetterReference()
    {
        // As slice 17's residual-score test: f(x) = 0, residuals N(shift, 1), n 15, α 0.1, 2,000 replications, a grid of step 0.05.
        const int n = 15, replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1150);
        double[] grid = [.. Enumerable.Range(0, 401).Select(g => -8.0 + 0.05 * g)];
        ResidualScore[] scores = [ResidualScore.DistanceToOrigin, ResidualScore.DistanceToAverage, ResidualScore.Shrinkage, ResidualScore.Horseshoe];
        foreach (double shift in new[] { 0.0, 3.0 })
        {
            var covered = new int[scores.Length];
            var width = new double[scores.Length];
            for (int r = 0; r < replications; r++)
            {
                double[] residuals = [.. Enumerable.Range(0, n).Select(_ => shift + Gaussian(rng))];
                double yt = shift + Gaussian(rng);
                int nearest = (int)Math.Round((yt + 8.0) / 0.05);
                for (int f = 0; f < scores.Length; f++)
                {
                    ConformalSet set = ResidualShrinkage.Prediction(residuals, [0.0], [grid], scores[f], 0.1).Sets[0];
                    if (nearest >= 0 && nearest < grid.Length && set.Accepts(nearest)) covered[f]++;
                    width[f] += set.Width / replications;
                }
            }
            double se = Math.Sqrt(0.09 / replications);
            output.WriteLine($"shift {shift}: " + string.Join("; ", scores.Select((s, f) => $"{s} coverage {covered[f] / (double)replications:F4}, width {width[f]:F3}"))
                             + $"; se {se:F4}");
            Assert.True(covered[3] / (double)replications >= 0.9 - 3.5 * se, $"RoBAS-Full at shift {shift}");
            if (shift == 0.0) Assert.True(width[3] < width[1], "unbiased: RoBAS-Full should beat distance to average");
            else Assert.True(width[3] < width[0], "shifted: RoBAS-Full should beat distance to origin");
        }
    }

    [Fact]
    public void UnderALargeShift_ItRanksAsDistanceToAverage()
    {
        // Residuals N(40, 1), n 15: their robustness theorem has the score tend to a monotone function of |r − r̄|, so the sets agree with
        // distance to average's on a fine grid; 200 replications.
        var rng = new Xoshiro256PlusPlus(2026_1151);
        double[] grid = [.. Enumerable.Range(0, 1001).Select(g => 35.0 + 0.01 * g)];
        int same = 0;
        double gap = 0.0;
        for (int r = 0; r < 200; r++)
        {
            double[] residuals = [.. Enumerable.Range(0, 15).Select(_ => 40.0 + Gaussian(rng))];
            ConformalSet horseshoe = ResidualShrinkage.Prediction(residuals, [0.0], [grid], ResidualScore.Horseshoe, 0.1).Sets[0];
            ConformalSet average = ResidualShrinkage.Prediction(residuals, [0.0], [grid], ResidualScore.DistanceToAverage, 0.1).Sets[0];
            int disagree = Enumerable.Range(0, grid.Length).Count(g => horseshoe.Accepts(g) != average.Accepts(g));
            if (disagree == 0) same++;
            gap = Math.Max(gap, Math.Abs(horseshoe.Width - average.Width));
        }
        output.WriteLine($"shift 40: identical sets in {same} of 200; largest width difference {gap:F3}");
        Assert.True(gap <= 0.05, $"width difference {gap}");

        // The grid-free form serves it from n 4 (robasFull_interval), and refuses it below, where its sets are not always intervals, and
        // on residuals all equal, outside the theorem.
        Assert.False(ResidualShrinkage.Intervals([.. Enumerable.Range(0, 20).Select(i => 0.1 * i)], [0.0], ResidualScore.Horseshoe).Intervals[0].WholeSpace);
        Assert.Contains("below n 4", Assert.Throws<ArgumentException>(
            () => ResidualShrinkage.Intervals([0.510312, 0.652786, 0.545799], [0.0], ResidualScore.Horseshoe, 0.3)).Message);
        Assert.Contains("not all equal", Assert.Throws<ArgumentException>(
            () => ResidualShrinkage.Intervals([.. Enumerable.Repeat(0.4, 12)], [0.0], ResidualScore.Horseshoe)).Message);
    }

    [Fact]
    public void TheCommonPoint_HasPValueOne()
    {
        // robas_common_point: at z₀ every unit's score reaches the candidate's, so the grid form's p-value is 1 there, on random residuals,
        // shifted, spread and tight, at n 4 to 40.
        var rng = new Xoshiro256PlusPlus(2026_1159);
        foreach (int n in new[] { 4, 5, 8, 15, 40 })
            foreach (double shift in new[] { 0.0, 0.7, -3.0, 25.0 })
            {
                double[] residuals = [.. Enumerable.Range(0, n).Select(_ => shift + (0.05 + rng.NextDouble()) * Gaussian(rng))];
                double z0 = ResidualShrinkage.CommonPoint(residuals);
                double s = residuals.Sum();
                Assert.InRange(z0, Math.Min(s / (3 * n), s / n) - 1e-12, Math.Max(s / (3 * n), s / n) + 1e-12);
                ConformalSet set = ResidualShrinkage.Prediction(residuals, [0.0], [[z0]], ResidualScore.Horseshoe, 0.5).Sets[0];
                Assert.Equal(1.0, set.PValues![0]);
            }
    }

    [Fact]
    public void ItsSets_AreNotAlwaysIntervals()
    {
        // The seeded search's counterexamples, read on the grid form: accepted about the prior's centre 0, refused between, accepted about
        // the data. n 3 at k 1 (α 0.3): pieces [−0.198, 0.062] and [0.333, 0.668]. n 2 at k 2 (α 0.9): [−8.03, 6.89] and [18.02, 19.32].
        ConformalSet three = ResidualShrinkage.Prediction([0.510312, 0.652786, 0.545799], [0.0], [[0.0, 0.2, 0.5]], ResidualScore.Horseshoe, 0.3).Sets[0];
        Assert.Equal(1, three.Level);
        Assert.Equal([true, false, true], [three.Accepts(0), three.Accepts(1), three.Accepts(2)]);
        ConformalSet two = ResidualShrinkage.Prediction([18.0184, 19.3278], [0.0], [[0.0, 12.0, 18.5]], ResidualScore.Horseshoe, 0.9).Sets[0];
        Assert.Equal(2, two.Level);
        Assert.Equal([true, false, true], [two.Accepts(0), two.Accepts(1), two.Accepts(2)]);
        // RoBAS-EB on the n 3 residuals at the same level: one piece through the same three values.
        ConformalSet eb = ResidualShrinkage.Prediction([0.510312, 0.652786, 0.545799], [0.0], [[0.0, 0.2, 0.5]], ResidualScore.Shrinkage, 0.3).Sets[0];
        Assert.False(eb.Accepts(0) && !eb.Accepts(1) && eb.Accepts(2));
    }
}
