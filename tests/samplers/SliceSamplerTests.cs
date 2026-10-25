using System;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The slice sampler. On a log-gamma target, one update from exact draws leaves the target in place, and a chain has its
/// moments whether the width is far too small or far too large; on a well-separated two-component mixture the chain
/// splits its time between the modes in proportion to their weights; and its arguments are checked.
/// </summary>
public sealed class SliceSamplerTests
{
    // u = ln of a Gamma(3, 1) variable: f(u) = 3u − e^u, mean ψ(3), variance ψ′(3).
    private static double LogGammaTarget(double u) => 3.0 * u - Math.Exp(u);
    private const double Digamma3 = 0.9227843350984671, Trigamma3 = 0.3949340668482264;

    [Fact]
    public void OneUpdate_FromExactDraws_StaysOnTheTarget()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1);
        const int n = 400_000;
        double before = 0, after = 0, afterSq = 0;
        for (int i = 0; i < n; i++)
        {
            // An exact Gamma(3, 1) draw as a sum of three exponentials, then its log.
            double g = -Math.Log(1.0 - rng.NextDouble()) - Math.Log(1.0 - rng.NextDouble()) - Math.Log(1.0 - rng.NextDouble());
            double u = Math.Log(g), f = LogGammaTarget(u);
            before += u;
            SliceSampler.Step(ref u, ref f, 1.0, 32, LogGammaTarget, rng);
            Assert.Equal(LogGammaTarget(u), f);
            after += u;
            afterSq += u * u;
        }
        double mean = after / n;
        Assert.InRange(before / n, Digamma3 - 0.004, Digamma3 + 0.004);
        Assert.InRange(mean, Digamma3 - 0.004, Digamma3 + 0.004);
        Assert.InRange(afterSq / n - mean * mean, Trigamma3 * 0.985, Trigamma3 * 1.015);
    }

    [Theory]
    [InlineData(0.01)]
    [InlineData(1.0)]
    [InlineData(100.0)]
    public void Chain_HasTheTargetsMoments_AtAnyWidth(double width)
    {
        var rng = new Xoshiro256PlusPlus(seed: 2);
        double u = 0.0, f = LogGammaTarget(u), sum = 0, sumSq = 0;
        long evaluations = 0;
        const int n = 200_000;
        for (int i = 0; i < n; i++)
        {
            evaluations += SliceSampler.Step(ref u, ref f, width, 1000, LogGammaTarget, rng);
            sum += u;
            sumSq += u * u;
        }
        double mean = sum / n;
        Assert.InRange(mean, Digamma3 - 0.01, Digamma3 + 0.01);
        Assert.InRange(sumSq / n - mean * mean, Trigamma3 * 0.97, Trigamma3 * 1.03);
        // Stepping out costs about (slice length)/w evaluations, shrinkage about log₂(w/slice length): the slice here is a
        // few units long.
        Assert.True(evaluations / (double)n < 15.0 + 3.0 / width, $"{evaluations / (double)n:F1} evaluations per update at width {width}");
    }

    [Fact]
    public void Chain_SplitsItsTimeBetweenModes_ByTheirWeights()
    {
        // 0.3·N(−4, 1) + 0.7·N(4, 1): the width spans the gap, so the chain crosses it.
        static double Mixture(double x) => Math.Log(0.3 * Math.Exp(-0.5 * (x + 4) * (x + 4)) + 0.7 * Math.Exp(-0.5 * (x - 4) * (x - 4)));
        var rng = new Xoshiro256PlusPlus(seed: 3);
        double x = 0.0, f = Mixture(x);
        int right = 0;
        const int n = 200_000;
        for (int i = 0; i < n; i++)
        {
            SliceSampler.Step(ref x, ref f, 10.0, 32, Mixture, rng);
            if (x > 0) right++;
        }
        Assert.InRange(right / (double)n, 0.68, 0.72);
    }

    [Fact]
    public void InvalidArguments_AreRefused()
    {
        var rng = new Xoshiro256PlusPlus(1);
        double x = 0.0, f = 0.0;
        Assert.Throws<ArgumentOutOfRangeException>(() => SliceSampler.Step(ref x, ref f, 0.0, 10, v => 0.0, rng));
        Assert.Throws<ArgumentOutOfRangeException>(() => SliceSampler.Step(ref x, ref f, 1.0, 0, v => 0.0, rng));
        double g = double.NegativeInfinity;
        Assert.Throws<ArgumentException>(() => SliceSampler.Step(ref x, ref g, 1.0, 10, v => 0.0, rng));
    }
}
