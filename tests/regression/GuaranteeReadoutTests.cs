using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The guarantee read out (E16, the conformal layer's slice 18): the law of L⁺ is exact, the Beta law at 0/1 losses and ℓ + β(B − ℓ) at one
/// loss (<c>gaffke_one</c>), and agrees with Dirichlet draws; Gaffke's bound b*_β is calibrated over replicated small calibration sets,
/// P(J ≤ b*_β) ≥ β, discrete losses included (Vlassis and Thomas's Theorem 1); conformal risk control holds its expected loss within α
/// (<c>crc_orderings</c>) and the high-probability rule its risk within α with chance β (<c>fail_le_fixed</c>); the cell scores are proper
/// (<c>logLoss_proper</c>, <c>gneitingLoss_gap</c>).
/// </summary>
public sealed class GuaranteeReadoutTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    private static double Choose(int n, int k)
    {
        double c = 1.0;
        for (int i = 1; i <= k; i++) c = c * (n - k + i) / i;
        return c;
    }

    [Fact]
    public void TheLaw_IsExact()
    {
        var rng = new Xoshiro256PlusPlus(1901);
        // One loss: b*_β = ℓ + β(B − ℓ).
        for (int r = 0; r < 50; r++)
        {
            double bound = 0.5 + 2.0 * rng.NextDouble(), loss = bound * rng.NextDouble(), beta = 0.05 + 0.9 * rng.NextDouble();
            Assert.Equal(loss + beta * (bound - loss), LossBound.Gaffke([loss], bound, beta), 12);
        }
        // 0/1 losses with B = 1: Beta(k + 1, n − k), whose CDF is P(Bin(n, x) ≥ k + 1).
        for (int r = 0; r < 60; r++)
        {
            int n = 1 + rng.NextInt(20), k = rng.NextInt(n + 1);
            double[] losses = [.. Enumerable.Range(0, n).Select(i => i < k ? 1.0 : 0.0)];
            var law = new LossBound(losses, 1.0);
            double x = rng.NextDouble(), expected = 0.0;
            for (int j = k + 1; j <= n; j++) expected += Choose(n, j) * Math.Pow(x, j) * Math.Pow(1 - x, n - j);
            Assert.Equal(expected, law.Cdf(x), 12);
        }
        // General losses, tied on a lattice or not, against 100,000 Dirichlet draws at five points.
        foreach (double[] losses in new[] { new[] { 0.1, 0.4, 0.4, 0.7, 0.2 }, new[] { 0.03, 0.55, 0.91 }, new[] { 0.5, 0.5, 0.5, 0.0, 0.25, 0.75 } })
        {
            var law = new LossBound(losses, 1.0);
            double[] knots = [.. losses, 1.0];
            const int draws = 100_000;
            var values = new double[draws];
            for (int d = 0; d < draws; d++)
            {
                double total = 0.0, sum = 0.0;
                foreach (double t in knots)
                {
                    double e = -Math.Log(1.0 - rng.NextDouble());
                    total += e;
                    sum += e * t;
                }
                values[d] = sum / total;
            }
            foreach (double x in new[] { 0.2, 0.35, 0.5, 0.65, 0.8 })
            {
                double empirical = values.Count(v => v <= x) / (double)draws, exact = law.Cdf(x), se = Math.Sqrt(exact * (1 - exact) / draws) + 1e-9;
                Assert.True(Math.Abs(empirical - exact) < 4.5 * se, $"losses [{string.Join(", ", losses)}], x {x}: exact {exact}, drawn {empirical}");
            }
            Assert.Equal((losses.Sum() + 1.0) / (losses.Length + 1), law.Mean, 12);
            Assert.Equal(1.0, law.Cdf(1.0));
            Assert.Equal(0.0, law.Cdf(-0.1));
        }
    }

    [Fact]
    public void TheBound_IsCalibrated_OverReplicatedCalibrationSets()
    {
        // P(J ≤ b*_β) ≥ β over 2,000 calibration sets each, for a Bernoulli(0.3) loss, a three-point loss on {0, ½, 1} with chances
        // (0.5, 0.3, 0.2), and U², U uniform, all with B = 1; n 3 and 8.
        const int replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1130);
        (string Name, Func<double> Draw, double Risk)[] laws =
        [
            ("Bernoulli(0.3)", () => rng.NextDouble() < 0.3 ? 1.0 : 0.0, 0.3),
            ("three points", () => { double u = rng.NextDouble(); return u < 0.5 ? 0.0 : u < 0.8 ? 0.5 : 1.0; }, 0.35),
            ("U²", () => { double u = rng.NextDouble(); return u * u; }, 1.0 / 3.0),
        ];
        foreach ((string name, Func<double> draw, double risk) in laws)
            foreach (int n in new[] { 3, 8 })
                foreach (double beta in name.StartsWith("Bernoulli") ? new[] { 0.75, 0.9 } : new[] { 0.9 })
                {
                    int held = 0;
                    for (int r = 0; r < replications; r++)
                    {
                        double[] losses = [.. Enumerable.Range(0, n).Select(_ => draw())];
                        if (risk <= LossBound.Gaffke(losses, 1.0, beta)) held++;
                    }
                    double rate = held / (double)replications, se = Math.Sqrt(beta * (1 - beta) / replications);
                    output.WriteLine($"{name}, n {n}, β {beta}: P(J ≤ b*) {rate:F4}; se {se:F4}");
                    Assert.True(rate >= beta - 3.5 * se, $"{name}, n {n}, β {beta}: {rate}");
                }
    }

    [Fact]
    public void TheRiskRules_ControlTheirLosses()
    {
        // Scores uniform on (0, 1), the set {s ≤ λ} on λ = j/100, miscoverage ℓᵢ(λ) = 1{sᵢ > λ}, B = 1, α 0.1, n 30, 1,500 calibration sets.
        // Conformal risk control's expected test loss is at most α; the high-probability rule's risk 1 − λ̂ exceeds α with chance at most 1 − β.
        const int n = 30, thresholds = 101, replications = 1500;
        const double alpha = 0.1, beta = 0.9;
        var rng = new Xoshiro256PlusPlus(2026_1131);
        double mean = 0.0;
        int failed = 0;
        for (int r = 0; r < replications; r++)
        {
            double[] s = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            var losses = new double[n, thresholds];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < thresholds; j++) losses[i, j] = s[i] > j / 100.0 ? 1.0 : 0.0;
            RiskThreshold crc = RiskControl.Mean(losses, 1.0, alpha);
            mean += rng.NextDouble() > crc.Index / 100.0 ? 1.0 : 0.0;
            RiskThreshold hpd = RiskControl.HighProbability(losses, 1.0, alpha, beta);
            Assert.True(hpd.Met);
            if (1.0 - hpd.Index / 100.0 > alpha) failed++;
        }
        double loss = mean / replications, fail = failed / (double)replications;
        double seLoss = Math.Sqrt(alpha * (1 - alpha) / replications), seFail = Math.Sqrt(beta * (1 - beta) / replications);
        output.WriteLine($"conformal risk control's test loss {loss:F4} (α {alpha}, se {seLoss:F4}); the high-probability rule fails {fail:F4} (1 − β {1 - beta:F2}, se {seFail:F4})");
        Assert.True(loss <= alpha + 3.5 * seLoss);
        Assert.True(fail <= 1 - beta + 3.5 * seFail);
    }

    [Fact]
    public void TheCellScores_AreProper()
    {
        // Exactly, on 500 pairs of laws on 2 to 8 cells: E_q[log(p)] − E_q[log(q)] = KL(q ‖ p) ≥ 0, and the quadratic gap Σ (p − q)².
        var rng = new Xoshiro256PlusPlus(1902);
        for (int r = 0; r < 500; r++)
        {
            int k = 2 + rng.NextInt(7);
            double[] p = Normalize([.. Enumerable.Range(0, k).Select(_ => 0.05 + rng.NextDouble())]);
            double[] q = Normalize([.. Enumerable.Range(0, k).Select(_ => 0.05 + rng.NextDouble())]);
            double Expect(Func<double[], int, double> score, double[] law) => Enumerable.Range(0, k).Sum(c => q[c] * score(law, c));
            double logGap = Expect(CellScores.Log, p) - Expect(CellScores.Log, q), kl = Enumerable.Range(0, k).Sum(c => q[c] * Math.Log(q[c] / p[c]));
            Assert.Equal(kl, logGap, 10);
            Assert.True(logGap > 0.0);
            Assert.Equal(Enumerable.Range(0, k).Sum(c => (p[c] - q[c]) * (p[c] - q[c])), Expect(CellScores.Quadratic, p) - Expect(CellScores.Quadratic, q), 12);
        }

        // On the split predictive distribution: the right fit scores better than one that ignores x by the quadratic score, over 400 test
        // cases. A fit wrong only in scale recalibrates to the same distribution, the PITs keeping their ranks, so the wrong fit here misses
        // the mean. The log score is +∞ wherever an observation falls in a cell no calibration PIT rises across: the cell law is a step
        // function of the ranks, and the log-loss's propriety needs the forecast positive where the truth is (logLoss_proper's hypothesis).
        double[] grid = [.. Enumerable.Range(0, 9).Select(g => -4.0 + 1.0 * g)];
        double[] x = [.. Enumerable.Range(0, 199).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] xt = [.. Enumerable.Range(0, 400).Select(_ => Gaussian(rng))], yt = [.. xt.Select(v => 0.5 * v + Gaussian(rng))];
        (double Log, double Quadratic, int Infinite) Scores(double slope)
        {
            double[] beta = [.. Enumerable.Repeat(slope, 20)];
            const double sd = 1.0;
            CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(beta, sd, null, x, [.. y.Select(v => new[] { v })]);
            CandidateLogDensity test = ConformalPredictiveDistributionTests.Seam(beta, sd, null, xt, [.. xt.Select(_ => grid)]);
            ConformalDistributions d = SplitConformal.Distribution(calibration, test);
            double log = 0.0, quadratic = 0.0;
            int infinite = 0;
            for (int t = 0; t < xt.Length; t++)
            {
                Assert.Equal(1.0, d.Distributions[t].CellMasses()!.Sum(), 12);
                double l = d.Distributions[t].LogScore(yt[t]);
                if (double.IsPositiveInfinity(l)) infinite++; else log += l;
                quadratic += d.Distributions[t].QuadraticScore(yt[t]);
            }
            return (log / (xt.Length - infinite), quadratic / xt.Length, infinite);
        }
        (double rightLog, double rightQuadratic, int rightInfinite) = Scores(0.5);
        (double wrongLog, double wrongQuadratic, int wrongInfinite) = Scores(0.0);
        output.WriteLine($"mean cell scores, the right fit against one ignoring x: quadratic {rightQuadratic:F5} and {wrongQuadratic:F5}; log over the finite "
                         + $"cases {rightLog:F4} and {wrongLog:F4}, +∞ in {rightInfinite} and {wrongInfinite} of {xt.Length}");
        Assert.True(rightQuadratic < wrongQuadratic);
    }

    private static double[] Normalize(double[] v)
    {
        double total = v.Sum();
        return [.. v.Select(e => e / total)];
    }

    [Fact]
    public void WhatTheReadoutDoesNotServe_IsRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new LossBound([0.2, 1.3], 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LossBound([0.2, double.NaN], 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LossBound([0.2], double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => LossBound.Gaffke([0.2], 1.0, 1.0));
        Assert.Throws<ArgumentException>(() => RiskControl.Mean(new double[,] { { 0.0, 1.0 } }, 1.0, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => RiskControl.Mean(new double[,] { { 2.0, 1.0 } }, 1.0, 0.1));
        Assert.Throws<ArgumentException>(() => RiskControl.Mean(new double[2, 0], 1.0, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => RiskControl.HighProbability(new double[,] { { 1.0, 0.0 } }, 1.0, 0.1, 0.0));
        // No threshold meets α: the last, flagged.
        RiskThreshold none = RiskControl.Mean(new double[,] { { 1.0, 1.0 }, { 1.0, 1.0 } }, 1.0, 0.1);
        Assert.False(none.Met);
        Assert.Equal(1, none.Index);
    }
}
