using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// GBP mixing (E1's slice 4). The slice step on ln v leaves v's conditional in place: exact draws of that conditional, by rejection from the
/// prior under the likelihood over its greatest value, each stepped once, against the conditional's CDF by a fixed-step trapezoid in u = ln v,
/// over the horseshoe, the Strawderman–Berger and a general GBP, residuals in the conditional's flat part and past it, drifts, powers and the
/// prior at b = 0; a long chain of steps carries the same mean of ln v. The improper conditional at a zero residual is refused. The marginal
/// is the mixture's integral and opts its readers in; the exact laws' steps are their draws, number for number.
/// </summary>
public sealed class GbpMixingTests(ITestOutputHelper output)
{
    private const double KolmogorovCut = 2.6;

    public static IEnumerable<object[]> Cases() =>
    [
        // p, a, b, r, s², m, power
        [1.0, 0.5, 0.5, 0.5, 1.0, 0.0, 1.0],
        [1.0, 0.5, 0.5, 3.0, 1.0, 0.0, 1.0],
        [1.0, 0.5, 0.5, 1e-3, 1.0, 0.0, 1.0],   // deep in the flat part, the wall near ln r²
        [1.0, 0.5, 0.5, -2.0, 0.5, 0.7, 1.0],
        [1.0, 0.5, 0.5, 1.0, 1.0, 0.0, 0.5],
        [1.0, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0],    // the prior
        [1.0, 0.5, 0.5, -2.0, 0.5, 0.7, 0.5],   // a drift on a rung (E1's slice 5)
        [1.0, 0.5, 0.5, 10.0, 1.0, 0.5, 0.5],   // far along the drift's ridge on a rung, r ≈ m·v near v = 20
        [0.5, 1.0, 0.5, 1.0, 2.0, -0.3, 1.0],
        [2.0, 1.5, 3.0, 0.8, 0.3, 0.2, 1.0],
    ];

    // ln of the conditional in u, unnormalized: the GBP's density in u and the likelihood at power b.
    private static Func<double, double> LogConditional(GeneralizedBetaPrime g, double r, double s2, double m, double power)
        => u => g.LogDensityInLog(u) + power * (-0.5 * u - Math.Pow(r - m * Math.Exp(u), 2.0) / (2.0 * s2 * Math.Exp(u)));

    // The conditional's CDF in u on a fixed grid, by the trapezoid, read back by linear interpolation.
    private static Func<double, double> Cdf(Func<double, double> h, double lo = -80.0, double hi = 150.0, double step = 1e-3)
    {
        int n = (int)Math.Round((hi - lo) / step);
        var values = new double[n + 1];
        double top = double.NegativeInfinity;
        for (int i = 0; i <= n; i++) top = Math.Max(top, values[i] = h(lo + i * step));
        var cumulative = new double[n + 1];
        for (int i = 1; i <= n; i++) cumulative[i] = cumulative[i - 1] + 0.5 * (Math.Exp(values[i - 1] - top) + Math.Exp(values[i] - top));
        double total = cumulative[n];
        return v =>
        {
            double at = (Math.Log(v) - lo) / step;
            if (at <= 0.0) return 0.0;
            if (at >= n) return 1.0;
            int k = (int)at;
            return (cumulative[k] + (at - k) * (cumulative[k + 1] - cumulative[k])) / total;
        };
    }

    private static double Kolmogorov(IEnumerable<double> xs, Func<double, double> cdf)
    {
        double[] sorted = [.. xs.Order()];
        double gap = 0.0;
        for (int i = 0; i < sorted.Length; i++)
        {
            double f = cdf(sorted[i]);
            gap = Math.Max(gap, Math.Max(Math.Abs(f - (double)i / sorted.Length), Math.Abs(f - (i + 1.0) / sorted.Length)));
        }
        return Math.Sqrt(sorted.Length) * gap;
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void OneStep_FromExactDraws_StaysOnTheConditional(double p, double a, double b, double r, double s2, double m, double power)
    {
        var mixing = new GbpMixing(new GbpParameters(p, a, b));
        GeneralizedBetaPrime g = GeneralizedBetaPrime.Of(p, a, b);
        var rng = new Xoshiro256PlusPlus(20261124);
        // The likelihood's greatest value over v, at m²v² + s²v − r² = 0.
        double LogLikelihood(double v) => power * (-0.5 * Math.Log(v) - Math.Pow(r - m * v, 2.0) / (2.0 * s2 * v));
        double best = m == 0.0 ? r * r / s2 : (-s2 + Math.Sqrt(s2 * s2 + 4.0 * m * m * r * r)) / (2.0 * m * m);
        double top = LogLikelihood(best);
        const int draws = 20_000;
        var exact = new double[draws];
        for (int k = 0; k < draws;)
        {
            double v = g.Sample(rng);
            if (Math.Log(1.0 - rng.NextDouble()) < LogLikelihood(v) - top) exact[k++] = v;
        }
        Func<double, double> cdf = Cdf(LogConditional(g, r, s2, m, power));
        double reference = Kolmogorov(exact, cdf);
        double[] stepped = [.. exact.Select(v => mixing.StepConditional(v, r, s2, m, power, rng))];
        double step = Kolmogorov(stepped, cdf);
        double moved = exact.Zip(stepped).Count(z => z.First != z.Second) / (double)draws;
        output.WriteLine($"GBP({p}, {a}, {b}), r {r}, s² {s2}, m {m}, b {power}: √n·D {reference:F2} exact, {step:F2} stepped; moved {moved:P1}");
        Assert.True(reference <= KolmogorovCut, $"the reference: √n·D {reference:F2}");
        Assert.True(step <= KolmogorovCut, $"one step: √n·D {step:F2}");
        Assert.True(moved > 0.99, $"moved {moved:P1}");
    }

    [Fact]
    public void AChainOfSteps_CarriesTheConditionalsMean()
    {
        // The horseshoe at a small residual, where the conditional is wide in u: a chain from v = 1, its mean of ln v by batch means.
        var mixing = new GbpMixing(GbpParameters.Horseshoe);
        GeneralizedBetaPrime g = GeneralizedBetaPrime.Of(1.0, 0.5, 0.5);
        const double r = 0.01, s2 = 1.0;
        Func<double, double> h = LogConditional(g, r, s2, 0.0, 1.0);
        double lo = -80.0, hi = 150.0, step = 1e-3, top = double.NegativeInfinity, mass = 0.0, first = 0.0;
        for (double u = lo; u <= hi; u += step) top = Math.Max(top, h(u));
        for (double u = lo; u <= hi; u += step)
        {
            double w = Math.Exp(h(u) - top);
            mass += w;
            first += w * u;
        }
        double expected = first / mass;
        var rng = new Xoshiro256PlusPlus(20261125);
        const int steps = 400_000, batches = 40;
        double v = 1.0;
        var logs = new double[steps];
        for (int k = 0; k < steps; k++) logs[k] = Math.Log(v = mixing.StepConditional(v, r, s2, 0.0, 1.0, rng));
        double[] means = [.. Enumerable.Range(0, batches).Select(k => logs.Skip(k * steps / batches).Take(steps / batches).Average())];
        double mean = means.Average(), se = Math.Sqrt(means.Sum(x => (x - mean) * (x - mean)) / (batches - 1) / batches);
        output.WriteLine($"E[ln v | r = {r}] {expected:F4} against the chain's {mean:F4} ± {se:F4}");
        Assert.True(Math.Abs(mean - expected) <= 4.0 * se, $"z {(mean - expected) / se:F2}");
    }

    [Fact]
    public void AZeroResidual_IsRefusedWhereTheConditionalIsImproper()
    {
        var rng = new Xoshiro256PlusPlus(1);
        Assert.Throws<ArgumentOutOfRangeException>(() => new GbpMixing(GbpParameters.Horseshoe).StepConditional(1.0, 0.0, 1.0, 0.0, 1.0, rng));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GbpMixing(GbpParameters.StrawdermanBerger).StepConditional(1.0, 0.0, 1.0, 0.4, 1.0, rng));
        // pa = 1 > ½: proper at 0; and the prior at b = 0 always.
        double v = new GbpMixing(new GbpParameters(1.0, 1.0, 1.0)).StepConditional(1.0, 0.0, 1.0, 0.0, 1.0, rng);
        Assert.True(v > 0.0 && double.IsFinite(v));
        Assert.True(new GbpMixing(GbpParameters.Horseshoe).StepConditional(1.0, 0.0, 1.0, 0.0, 0.0, rng) > 0.0);
        Assert.Throws<ArgumentOutOfRangeException>(() => new GbpMixing(GbpParameters.Horseshoe).StepConditional(0.0, 1.0, 1.0, 0.0, 1.0, rng));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GbpMixing(new GbpParameters(1.0, 0.0, 1.0)));
    }

    [Fact]
    public void TheMarginal_IsTheMixturesDensity_AndOptsItsReadersIn()
    {
        var mixing = new GbpMixing(new GbpParameters(1.0, 1.0, 3.0));
        foreach ((double r, double s2, double m) in new[] { (0.3, 1.0, 0.0), (-2.0, 0.5, 0.4), (5.0, 2.0, -1.0) })
            Assert.Equal(mixing.Noise(s2, m).LogPdf(r), mixing.MarginalLogDensity(r, s2, m));
        Assert.True(((IMixingLaw)mixing).MarginalByQuadrature);
        Assert.False(((IMixingLaw)new StudentTMixing(4)).MarginalByQuadrature);
        Assert.False(((IMixingLaw)ExponentialMixing.Instance).MarginalByQuadrature);
        Assert.False(((IMixingLaw)GigMixing.NormalWald(1.5)).MarginalByQuadrature);
        // GBP(1, a, b): E[v] = a/(b − 1), E[v²] = a(a + 1)/((b − 1)(b − 2)).
        Assert.Equal(0.5, mixing.Mean, 12);
        Assert.Equal(1.0 - 0.25, mixing.Variance, 12);
        var horseshoe = new GbpMixing(GbpParameters.Horseshoe);
        Assert.Equal(double.PositiveInfinity, horseshoe.Mean);
        Assert.Equal(double.PositiveInfinity, horseshoe.Variance);
    }

    [Fact]
    public void TheExactLaws_StepByTheirDraw_NumberForNumber()
    {
        IMixingLaw[] laws = [new StudentTMixing(4), ExponentialMixing.Instance, GigMixing.NormalWald(1.5), new GigMixing(new GigParameters(0.7, 1.2, 0.8))];
        foreach (IMixingLaw law in laws)
            foreach (double m in new[] { 0.0, 0.4 })
            {
                var a = new Xoshiro256PlusPlus(7);
                var b = new Xoshiro256PlusPlus(7);
                double stepped = law.StepConditional(123.0, 0.8, 1.3, m, 1.0, a);
                double drawn = law switch
                {
                    StudentTMixing t => t.SampleConditional(0.8, 1.3, m, 1.0, b),
                    ExponentialMixing e => e.SampleConditional(0.8, 1.3, m, 1.0, b),
                    GigMixing g => g.SampleConditional(0.8, 1.3, m, 1.0, b),
                    _ => throw new InvalidOperationException(),
                };
                Assert.Equal(drawn, stepped);
            }
    }
}
