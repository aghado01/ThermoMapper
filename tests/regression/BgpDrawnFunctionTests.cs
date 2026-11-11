using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Bgp;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BGP's retained-parameters content (E16, the conformal layer's slice 11): g = m + f drawn jointly at the training and the cases' inputs for
/// each retained draw. A candidate at a training input with that observation's value is the training column, bit for bit, and the full step
/// takes the pair. The draws meet the GP's conditional moments, against a reference by the bordered kriging system, an algebra the engine
/// does not use. With the qualification switch (<c>BASS_QUALIFICATION</c> = 1): held-out coverage of the full step from truth outside the
/// model, the split form beside it at equal data, and the cost row.
/// </summary>
public sealed class BgpDrawnFunctionTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[,] X, double[] Y) GpData(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n, 1];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = (i + 0.5) / n;
            y[i] = 3.0 + Math.Sin(2.0 * Math.PI * x[i, 0]) + 0.2 * Gaussian(rng);
        }
        return (x, y);
    }

    private static GpRegression Gp(double[,] x, double[] y, bool unknownScale) => unknownScale
        ? GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel(), GpInputScaling.None)
        : new GpRegression(x, y, 0.04, new SquaredExponentialKernel(), GpInputScaling.None);

    private static double[][] Singletons(IEnumerable<double> y) => [.. y.Select(v => new[] { v })];

    private static double[,] Column(IReadOnlyList<double> x)
    {
        var c = new double[x.Count, 1];
        for (int i = 0; i < x.Count; i++) c[i, 0] = x[i];
        return c;
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void AtTrainingInputs_TheCandidateIsTheTrainingColumn_BitForBit_AndTheFullStepTakesThePair(bool unknownScale)
    {
        (double[,] x, double[] y) = GpData(40, seed: 8);
        GpRegression gp = Gp(x, y, unknownScale);
        var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
        BgpResult r = sampler.Run(x, draws: 150, burn: 100, seed: 4, chains: 2);
        // The cases: every training input, two new inputs, and the first new input again.
        double[] xs = [.. Enumerable.Range(0, 40).Select(i => x[i, 0]), 0.2137, 0.6403, 0.2137];
        BgpDrawnFunction drawn = sampler.DrawFunction(r, Column(xs), seed: 11);
        Assert.Equal(r.Draws, drawn.Draws);
        Assert.Equal(42, drawn.Points);

        PointwiseLogDensity training = drawn.PointwiseLogDensity(cdf: true);
        Assert.False(training.LeavesUnitsOut);
        CandidateLogDensity atTraining = drawn.CandidateLogDensity(Singletons(y), cdf: true, cases: [.. Enumerable.Range(0, 40)]);
        for (int d = 0; d < drawn.Draws; d++)
        {
            for (int i = 0; i < 40; i++)
            {
                Assert.Equal(training.Values[d, i], atTraining.Values[d, i]);
                Assert.Equal(training.LogCdf![d, i], atTraining.LogCdf![d, i]);
                Assert.Equal(drawn.AtObservation(d, i), drawn.AtCase(d, i));
            }
            Assert.Equal(drawn.AtCase(d, 40), drawn.AtCase(d, 42));
            Assert.Equal(unknownScale ? r.BandwidthDraws.SelectMany(c => c).ElementAt(d).NoiseVariance : 0.04, drawn.NoiseVariance(d));
        }

        // Each term is the Normal of its draw, written out.
        int dr = 17, unit = 23;
        double v = drawn.NoiseVariance(dr), e = y[unit] - drawn.AtObservation(dr, unit);
        Assert.Equal(-0.5 * Math.Log(2.0 * Math.PI * v) - 0.5 * e * e / v, training.Values[dr, unit], 12);

        // The full step and its predictive distribution take the pair; the collapsed content stays refused.
        double[] grid = [.. Enumerable.Range(0, 81).Select(g => 1.0 + 0.05 * g)];
        CandidateLogDensity test = drawn.CandidateLogDensity([grid, grid], cdf: true, cases: [40, 41]);
        Assert.Equal(ConformalForms.Both, test.Forms);
        ConformalPrediction sets = ConformalBayes.FullStep(training, test);
        Assert.All(sets.Sets, s => Assert.Equal(41, s.Units));
        Assert.All(sets.Sets, s => Assert.True(s.Width > 0.0 && s.Width < 4.0));
        Assert.Equal(2, ConformalBayes.Distribution(training, test).Distributions.Count);
        Assert.Contains("leave-unit-out", Assert.Throws<ArgumentException>(
            () => ConformalBayes.FullStep(sampler.PointwiseLogDensity(r), test)).Message);

        // One seed repeats the draws exactly; another draws afresh.
        BgpDrawnFunction again = sampler.DrawFunction(r, Column(xs), seed: 11), other = sampler.DrawFunction(r, Column(xs), seed: 12);
        Assert.Equal(drawn.AtCase(9, 41), again.AtCase(9, 41));
        Assert.NotEqual(drawn.AtCase(9, 41), other.AtCase(9, 41));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void TheDraws_MeetTheGpsConditionalMoments(bool unknownScale)
    {
        // One state held fixed through 6,000 draws: g at three training inputs and three new ones against the universal-kriging mean and
        // covariance from the bordered system [[R, 1], [1ᵀ, 0]], solved by elimination here.
        (double[,] x, double[] y) = GpData(30, seed: 3);
        GpRegression gp = Gp(x, y, unknownScale);
        var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
        BgpResult r = sampler.Run(x, draws: 20, burn: 20, seed: 2, chains: 1);
        const double t = 0.02, eta = 0.05, s2 = 0.7;
        GpFit fit = unknownScale ? gp.Fit(t, eta) : gp.Fit(t);
        double nugget = unknownScale ? eta : 0.04, scale = unknownScale ? s2 : 1.0;
        const int draws = 6000;
        BgpDraw[] chain = [.. Enumerable.Range(0, draws).Select(_ => new BgpDraw(t, fit, [], scale))];
        BgpResult held = r with { BandwidthDraws = [chain], Draws = draws };
        double[] points = [x[4, 0], x[15, 0], x[16, 0], 0.4012, 0.9990, 1.0400];
        BgpDrawnFunction drawn = sampler.DrawFunction(held, Column(points[3..]), seed: 5);

        (double[] mean, double[,] cov) = Kriging(x, y, t, nugget, points);
        double G(int d, int a) => a < 3 ? drawn.AtObservation(d, a switch { 0 => 4, 1 => 15, _ => 16 }) : drawn.AtCase(d, a - 3);
        var sampleMean = new double[6];
        for (int a = 0; a < 6; a++) sampleMean[a] = Enumerable.Range(0, draws).Average(d => G(d, a));
        for (int a = 0; a < 6; a++)
        {
            double variance = scale * cov[a, a];
            double sampleVar = Enumerable.Range(0, draws).Sum(d => (G(d, a) - sampleMean[a]) * (G(d, a) - sampleMean[a])) / (draws - 1);
            output.WriteLine($"point {points[a]:F4}: mean {sampleMean[a]:F5} against {mean[a]:F5}, variance {sampleVar:E4} against {variance:E4}");
            Assert.True(Math.Abs(sampleMean[a] - mean[a]) <= 4.5 * Math.Sqrt(variance / draws), $"mean at {points[a]}");
            Assert.True(Math.Abs(sampleVar / variance - 1.0) <= 4.5 * Math.Sqrt(2.0 / draws), $"variance at {points[a]}");
        }
        foreach ((int a, int b) in new[] { (1, 2), (0, 3), (3, 4), (4, 5) })
        {
            double reference = cov[a, b] / Math.Sqrt(cov[a, a] * cov[b, b]);
            double sample = Enumerable.Range(0, draws).Sum(d => (G(d, a) - sampleMean[a]) * (G(d, b) - sampleMean[b]))
                            / Math.Sqrt(Enumerable.Range(0, draws).Sum(d => (G(d, a) - sampleMean[a]) * (G(d, a) - sampleMean[a]))
                                        * Enumerable.Range(0, draws).Sum(d => (G(d, b) - sampleMean[b]) * (G(d, b) - sampleMean[b])));
            output.WriteLine($"correlation of {points[a]:F4} and {points[b]:F4}: {sample:F4} against {reference:F4}");
            Assert.True(Math.Abs(sample - reference) <= 4.5 * (1.0 - reference * reference) / Math.Sqrt(draws) + 1e-3, $"correlation {a}, {b}");
        }
    }

    // The universal-kriging mean and covariance of g = m + f at the points, from the bordered system M = [[R, 1], [1ᵀ, 0]],
    // R = Σ_t + nugget·I: with M[λ; ν] = [k(X, z); 1], the mean is λᵀY and the covariance k(z, z′) − [k(X, z); 1]ᵀM⁻¹[k(X, z′); 1].
    private static (double[] Mean, double[,] Covariance) Kriging(double[,] x, double[] y, double t, double nugget, double[] points)
    {
        int n = y.Length, p = points.Length;
        double K(double a, double b) => Math.Exp(-(a - b) * (a - b) / (2.0 * t));
        var m = new double[n + 1, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) m[i, j] = K(x[i, 0], x[j, 0]) + (i == j ? nugget : 0.0);
            m[i, n] = m[n, i] = 1.0;
        }
        var rhs = new double[p][];
        for (int a = 0; a < p; a++)
        {
            rhs[a] = new double[n + 1];
            for (int i = 0; i < n; i++) rhs[a][i] = K(x[i, 0], points[a]);
            rhs[a][n] = 1.0;
        }
        double[][] solved = [.. rhs.Select(b => Eliminate(m, b))];
        var mean = new double[p];
        var cov = new double[p, p];
        for (int a = 0; a < p; a++)
        {
            for (int i = 0; i < n; i++) mean[a] += solved[a][i] * y[i];
            for (int b = 0; b < p; b++) cov[a, b] = K(points[a], points[b]) - rhs[a].Zip(solved[b], (u, w) => u * w).Sum();
        }
        return (mean, cov);
    }

    // Gaussian elimination with partial pivoting on a copy.
    private static double[] Eliminate(double[,] a0, double[] b0)
    {
        int n = b0.Length;
        var a = (double[,])a0.Clone();
        var b = (double[])b0.Clone();
        for (int c = 0; c < n; c++)
        {
            int pivot = c;
            for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[pivot, c])) pivot = r;
            if (pivot != c)
            {
                for (int k = 0; k < n; k++) (a[c, k], a[pivot, k]) = (a[pivot, k], a[c, k]);
                (b[c], b[pivot]) = (b[pivot], b[c]);
            }
            for (int r = c + 1; r < n; r++)
            {
                double f = a[r, c] / a[c, c];
                if (f == 0.0) continue;
                for (int k = c; k < n; k++) a[r, k] -= f * a[c, k];
                b[r] -= f * b[c];
            }
        }
        var xs = new double[n];
        for (int r = n - 1; r >= 0; r--)
        {
            double s = b[r];
            for (int k = r + 1; k < n; k++) s -= a[r, k] * xs[k];
            xs[r] = s / a[r, r];
        }
        return xs;
    }

    [Fact]
    public void ThePivotedFactor_ReproducesTheNearSingularCovariance()
    {
        // 30 inputs at spacing 1/30 under a bandwidth whose length is about four spacings, and three cases: the posterior covariance's
        // spectrum falls below rounding, where an unpivoted factor with a floor diverges. The pivoted one stops short and reproduces it.
        (double[,] x, double[] y) = GpData(30, seed: 3);
        GpRegression gp = Gp(x, y, unknownScale: false);
        double[] points = [.. Enumerable.Range(0, 30).Select(i => x[i, 0]), 0.4012, 0.9990, 1.0400];
        (_, double[,] c) = gp.JointPosteriorInFrame(gp.Fit(0.02), Column(points));
        int p = points.Length;
        var l = new double[p * p];
        int rank = BgpSampler.PivotedFactor(c, l);
        Assert.All(l, v => Assert.True(double.IsFinite(v)));
        double largest = Enumerable.Range(0, p).Max(a => c[a, a]), worst = 0.0;
        for (int a = 0; a < p; a++)
            for (int b = 0; b <= a; b++)
            {
                double s = 0.0;
                for (int k = 0; k < rank; k++) s += l[a * p + k] * l[b * p + k];
                worst = Math.Max(worst, Math.Abs(s - c[a, b]));
            }
        output.WriteLine($"rank {rank} of {p}; largest gap {worst:E2} against the largest variance {largest:E3}");
        Assert.True(worst <= 1e-10 * largest, $"gap {worst:E2}");
        Assert.True(rank < p);
    }

    [Fact]
    public void WhatTheContentDoesNotServe_IsRefused()
    {
        (double[,] x, double[] y) = GpData(20, seed: 6);
        GpRegression gp = Gp(x, y, unknownScale: true);
        var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
        BgpResult r = sampler.Run(x, draws: 30, burn: 30, seed: 1, chains: 1);
        Assert.Throws<ArgumentException>(() => sampler.DrawFunction(r, new double[2, 2]));
        BgpDrawnFunction drawn = sampler.DrawFunction(r, Column([0.3, 0.7]));
        Assert.Throws<ArgumentOutOfRangeException>(() => drawn.CandidateLogDensity([[1.0]], cases: [2]));
        Assert.Throws<ArgumentException>(() => drawn.CandidateLogDensity([[1.0]]));   // one grid for two cases
        var pairs = new UnitPartition([.. Enumerable.Range(0, 20).Select(i => i / 2)]);
        Assert.Equal(10, drawn.PointwiseLogDensity(pairs).Units.Units);
        Assert.Throws<ArgumentException>(() => drawn.PointwiseLogDensity(pairs, cdf: true));
        // No cases: the training density alone, the content WAIC reads.
        Assert.Equal(0, sampler.DrawFunction(r).Cases);
    }

    // ---- Held-out coverage and the cost row, behind the qualification switch ----

    // Slice 5's truth: f(x) = sin(2πx) + x/2, with noise 0.2ζ plus, one time in seven, 0.8·Exp(1).
    private static double Truth(double x, Xoshiro256PlusPlus rng)
        => Math.Sin(2.0 * Math.PI * x) + 0.5 * x + 0.2 * Gaussian(rng) + (rng.NextDouble() < 1.0 / 7.0 ? -0.8 * Math.Log(1.0 - rng.NextDouble()) : 0.0);

    private static double[] Uniforms(int n, Xoshiro256PlusPlus rng) => [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];

    private static double[][] Grids(double[] truths, double[] grid, int widths) => [.. truths.Select(v => new[] { v }), .. Enumerable.Repeat(grid, widths)];

    [Fact]
    public void HeldOutCoverage_TheFullStepOnTheDrawnFunction_BesideTheSplitForm()
    {
        // 20 fits of 99 points from slice 5's truth, outside the model: the full step on the drawn function (m 100, k 10), and the split form
        // on the same 99 points, 50 to fit and 49 to calibrate (m 50, k 5). 50 held-out cases per fit at their true value, 20 on a grid.
        if (!BassQualification.Enabled) return;
        const int fits = 20, n = 99, tests = 50, widths = 20;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        var full = new double[fits];
        var split = new double[fits];
        double fullWidth = 0.0, splitWidth = 0.0, leastEss = double.PositiveInfinity;
        int fullOpen = 0, splitOpen = 0, fullPieces = 0;
        for (int fit = 0; fit < fits; fit++)
        {
            var rng = new Xoshiro256PlusPlus(20261111 + 97 * fit);
            double[] x = Uniforms(n, rng), y = [.. x.Select(v => Truth(v, rng))];
            double[] xt = Uniforms(tests + widths, rng), yt = [.. xt[..tests].Select(v => Truth(v, rng))];

            GpRegression gp = GpRegression.WithUnknownScale(Column(x), y, new SquaredExponentialKernel());
            var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
            BgpResult r = sampler.Run(Column(x), draws: 600, burn: 200, seed: fit + 1, chains: 2);
            BgpDrawnFunction drawn = sampler.DrawFunction(r, Column(xt), seed: fit + 101);
            ConformalPrediction p = ConformalBayes.FullStep(drawn.PointwiseLogDensity(), drawn.CandidateLogDensity(Grids(yt, grid, widths)));

            GpRegression half = GpRegression.WithUnknownScale(Column(x[..50]), y[..50], new SquaredExponentialKernel());
            var halfSampler = new BgpSampler(half, new EmpiricalBayesBandwidthPrior(half));
            BgpResult hr = halfSampler.Run(Column(x[..50]), draws: 600, burn: 200, seed: fit + 1, chains: 2);
            ConformalPrediction q = SplitConformal.Density(halfSampler.CandidateLogDensity(hr, Column(x[50..]), Singletons(y[50..])),
                                                           halfSampler.CandidateLogDensity(hr, Column(xt), Grids(yt, grid, widths)));
            int fc = 0, sc = 0;
            for (int t = 0; t < tests; t++)
            {
                Assert.Equal(100, p.Sets[t].Units);
                Assert.Equal(50, q.Sets[t].Units);
                if (p.Sets[t].Accepts(0)) fc++;
                if (q.Sets[t].Accepts(0)) sc++;
                leastEss = Math.Min(leastEss, p.Sets[t].EffectiveSampleSizes![0]);
            }
            for (int t = 0; t < widths; t++)
            {
                ConformalSet a = p.Sets[tests + t], b = q.Sets[tests + t];
                fullWidth += a.Width;
                splitWidth += b.Width;
                if (a.OpenBelow || a.OpenAbove) fullOpen++;
                if (b.OpenBelow || b.OpenAbove) splitOpen++;
                if (a.Intervals.Count > 1) fullPieces++;
            }
            full[fit] = fc / (double)tests;
            split[fit] = sc / (double)tests;
        }
        (double Mean, double Se) Of(double[] c)
        {
            double mean = c.Average();
            return (mean, Math.Sqrt(c.Sum(v => (v - mean) * (v - mean)) / (fits - 1) / fits));
        }
        (double fm, double fse) = Of(full);
        (double sm, double sse) = Of(split);
        output.WriteLine($"the full step on the drawn function: coverage {fm:F4} (se {fse:F4}, nominal 0.9000), mean width {fullWidth / (fits * widths):F3}, "
                         + $"{fullPieces} of {fits * widths} sets in two or more pieces, {fullOpen} open at the grid's ends; least ESS at the truth {leastEss:F1}");
        output.WriteLine($"the split form, 50 to fit and 49 to calibrate: coverage {sm:F4} (se {sse:F4}, nominal 0.9000), mean width {splitWidth / (fits * widths):F3}, "
                         + $"{splitOpen} open at the grid's ends");
        int count = fits * tests;
        Assert.True(Math.Abs(fm - 0.9) <= 3.5 * Math.Max(fse, Math.Sqrt(0.09 / count)), $"the full step's coverage {fm:F4}");
        Assert.True(Math.Abs(sm - 0.9) <= 3.5 * Math.Max(sse, Math.Sqrt(0.09 / count)), $"the split form's coverage {sm:F4}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(50)]
    public void CostRow(int cases)
    {
        // n 400, S 1,000, a grid of 200 per case: the drawn function, its training and candidate densities, and the full step, beside the
        // split form's seam (one refit per accepted state) and step. Indicative: the machine is shared.
        if (!BassQualification.Enabled) return;
        const int n = 400;
        var rng = new Xoshiro256PlusPlus(20261111);
        double[] x = Uniforms(n, rng), y = [.. x.Select(v => Truth(v, rng))];
        double[] xt = Uniforms(cases, rng);
        GpRegression gp = GpRegression.WithUnknownScale(Column(x), y, new SquaredExponentialKernel());
        var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
        BgpResult r = sampler.Run(Column(x), draws: 500, burn: 200, seed: 1, chains: 2);
        int states = 0;
        GpFit? previous = null;
        foreach (BgpDraw d in r.BandwidthDraws.SelectMany(c => c))
        {
            if (previous is null || d.Fit.Bandwidth != previous.Bandwidth || d.Fit.NoiseRatio != previous.NoiseRatio) states++;
            previous = d.Fit;
        }
        double[] grid = [.. Enumerable.Range(0, 200).Select(g => -3.0 + 0.03 * g)];
        double[][] grids = [.. Enumerable.Repeat(grid, cases)];
        double Ms(long from) => (BassCost.Tick() - from) / BassCost.TicksPerSecond * 1e3;

        long t0 = BassCost.Tick();
        BgpDrawnFunction drawn = sampler.DrawFunction(r, Column(xt), seed: 3);
        double draw = Ms(t0);
        t0 = BassCost.Tick();
        PointwiseLogDensity training = drawn.PointwiseLogDensity();
        double trainingMs = Ms(t0);
        t0 = BassCost.Tick();
        CandidateLogDensity candidates = drawn.CandidateLogDensity(grids);
        double candidateMs = Ms(t0);
        t0 = BassCost.Tick();
        ConformalBayes.FullStep(training, candidates);
        double step = Ms(t0);

        t0 = BassCost.Tick();
        CandidateLogDensity seam = sampler.CandidateLogDensity(r, Column(xt), grids);
        double seamMs = Ms(t0);
        CandidateLogDensity calibration = sampler.CandidateLogDensity(r, Column(x[..99]), Singletons(y[..99]));
        t0 = BassCost.Tick();
        SplitConformal.Density(calibration, seam);
        double splitStep = Ms(t0);
        output.WriteLine($"n {n}, S {r.Draws} ({states} distinct states), {cases} cases ({drawn.Points} points), a grid of 200 each "
                         + $"({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): the drawn function {draw:F0} ms "
                         + $"({draw / states:F1} ms per state), its training density {trainingMs:F1} ms, its candidates {candidateMs:F1} ms, "
                         + $"the full step {step:F0} ms; the split form's seam {seamMs:F0} ms and its step {splitStep:F1} ms");
    }
}
