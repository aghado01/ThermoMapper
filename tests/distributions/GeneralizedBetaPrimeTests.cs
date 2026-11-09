using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// The generalized beta prime (E1's slice 4) and the normal variance-mean mixtures over it and over the GIG. The GBP: its law against the
/// half-Cauchy's at the horseshoe, where √v is half-Cauchy; its tails and levels against Simpson's rule, and its mass and moments against the
/// trapezoid, over its density in u = ln v on chosen laws; over a seeded sweep across far ranges its tails' complement, its quantiles back
/// through its tails and its density against the beta's; its draws, and its draws on an interval, by Kolmogorov–Smirnov. The mixtures: the
/// density by the integral over v against the generalized hyperbolic's closed form and the horseshoe's exponential integral; density and
/// tails against Simpson's rule at a fixed step in u, refined about the drift's narrow peak, which shares nothing adaptive with the code; the
/// tail table against the integral out to the 10⁻⁹ quantiles and across the pole at 0; the draws.
/// </summary>
public sealed class GeneralizedBetaPrimeTests(ITestOutputHelper output)
{
    private const double KolmogorovCut = 2.6;   // P(√n·D > 2.6) ≈ 3 in a million per case

    // ln ∫ exp(h(u)) du over [lo, hi] by the trapezoid at a fixed step, the sum kept in logs as it goes.
    private static double Trapezoid(Func<double, double> h, double lo, double hi, double step)
    {
        int n = Math.Max(1, (int)Math.Ceiling((hi - lo) / step));
        double width = (hi - lo) / n, top = double.NegativeInfinity, sum = 0.0;
        for (int i = 0; i <= n; i++)
        {
            double value = h(lo + i * width) + (i == 0 || i == n ? -Math.Log(2.0) : 0.0);
            if (double.IsNegativeInfinity(value)) continue;
            if (value > top) { sum = sum * Math.Exp(top - value) + 1.0; top = value; }
            else sum += Math.Exp(value - top);
        }
        return top + Math.Log(sum * width);
    }

    // ln ∫ exp(h(u)) du over [lo, hi] by Simpson's rule at a fixed step: for an integral cut where the integrand is not small, whose cut end the
    // trapezoid gets only to the step's square.
    private static double Simpson(Func<double, double> h, double lo, double hi, double step)
    {
        int n = 2 * Math.Max(1, (int)Math.Ceiling((hi - lo) / (2.0 * step)));
        double width = (hi - lo) / n, top = double.NegativeInfinity, sum = 0.0;
        for (int i = 0; i <= n; i++)
        {
            double value = h(lo + i * width) + Math.Log(i == 0 || i == n ? 1.0 : i % 2 == 1 ? 4.0 : 2.0);
            if (double.IsNegativeInfinity(value)) continue;
            if (value > top) { sum = sum * Math.Exp(top - value) + 1.0; top = value; }
            else sum += Math.Exp(value - top);
        }
        return top + Math.Log(sum * width / 3.0);
    }

    private static double Kolmogorov(double[] xs, Func<double, double> cdf)
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

    // ---- the generalized beta prime ------------------------------------------------------------------------------------------------

    [Fact]
    public void TheHorseshoesMixing_IsTheHalfCauchyOnRootV()
    {
        var gbp = GeneralizedBetaPrime.Of(1.0, 0.5, 0.5);
        var cauchy = HalfCauchy.WithScale(1.0);
        double worst = 0.0;
        for (int k = -24; k <= 24; k++)
        {
            double v = Math.Pow(10.0, 0.5 * k), root = Math.Sqrt(v);
            worst = Math.Max(worst, Math.Abs(gbp.LogCdf(v) - cauchy.LogCdf(root)) / Math.Max(1.0, Math.Abs(cauchy.LogCdf(root))));
            worst = Math.Max(worst, Math.Abs(gbp.LogSurvival(v) - cauchy.LogSurvival(root)) / Math.Max(1.0, Math.Abs(cauchy.LogSurvival(root))));
            worst = Math.Max(worst, Math.Abs(gbp.LogPdf(v) - (cauchy.LogPdf(root) - Math.Log(2.0 * root))) / Math.Max(1.0, Math.Abs(gbp.LogPdf(v))));
        }
        foreach (double p in new[] { 1e-12, 1e-6, 0.01, 0.3, 0.5, 0.7, 0.99, 1.0 - 1e-6, 1.0 - 1e-12 })
            worst = Math.Max(worst, Math.Abs(gbp.Quantile(p) - Math.Pow(cauchy.Quantile(p), 2.0)) / Math.Pow(cauchy.Quantile(p), 2.0));
        output.WriteLine($"GBP(1, ½, ½) against C⁺(0, 1) on √v: worst relative gap {worst:G3}");
        Assert.True(worst <= 1e-11, $"worst {worst:G3}");
        Assert.Equal(0.0, gbp.ModeInLog, 15);
        Assert.Equal(double.PositiveInfinity, gbp.Mean);
        Assert.Equal(double.PositiveInfinity, gbp.Variance);
    }

    public static IEnumerable<object[]> Chosen() =>
    [
        [1.0, 0.5, 0.5], [0.5, 1.0, 0.5], [2.0, 1.5, 3.0], [0.3, 4.0, 0.7], [5.0, 0.2, 2.0], [1.0, 30.0, 40.0], [1.0, 1.0, 3.0],
    ];

    [Theory]
    [MemberData(nameof(Chosen))]
    public void TailsAndMoments_AreTheDensitysIntegral(double p, double a, double b)
    {
        var gbp = GeneralizedBetaPrime.Of(p, a, b);
        // The density in u falls at slope pa below its mode and pb above: 70 nats out on each side, at a step well inside its curvature.
        double mode = gbp.ModeInLog, lo = mode - 70.0 / (p * a) - 5.0, hi = mode + 70.0 / (p * b) + 5.0;
        // Simpson's error goes as (step·rate)⁴ on a tail falling at rate pa or pb in u.
        double step = Math.Min(2e-3, 0.01 / (p * Math.Max(a, b)));
        double worst = 0.0;
        // Each tail's mass q, below and above: the quantile, its tails against the integral, and its mass against q.
        foreach (double q in new[] { 1e-9, 1e-4, 0.01, 0.2, 0.5 })
            foreach (bool fromBelow in new[] { true, false })
            {
                double v = fromBelow ? gbp.InverseLogCdf(Math.Log(q)) : gbp.InverseLogSurvival(Math.Log(q)), u = Math.Log(v);
                double below = Simpson(gbp.LogDensityInLog, lo, u, step), above = Simpson(gbp.LogDensityInLog, u, hi, step);
                worst = Math.Max(worst, Math.Abs(gbp.LogCdf(v) - below) / Math.Max(1.0, Math.Abs(below)));
                worst = Math.Max(worst, Math.Abs(gbp.LogSurvival(v) - above) / Math.Max(1.0, Math.Abs(above)));
                worst = Math.Max(worst, Math.Abs(Math.Exp(fromBelow ? below : above) - q) / q);
            }
        Assert.True(worst <= 1e-9, $"GBP({p}, {a}, {b}): tails against the integral, worst relative gap {worst:G3}");
        double total = Trapezoid(gbp.LogDensityInLog, lo, hi, step);
        Assert.True(Math.Abs(total) <= 1e-10, $"the density's mass, ln {total:G3}");
        foreach (double k in new[] { 1.0, 2.0 })
        {
            double moment = k / p < b ? Trapezoid(u => gbp.LogDensityInLog(u) + k * u, lo, hi + 70.0 / (p * b - k + 1e-9), step) : double.PositiveInfinity;
            double read = gbp.LogMoment(k);
            if (double.IsPositiveInfinity(moment)) Assert.Equal(double.PositiveInfinity, read);
            else Assert.True(Math.Abs(read - moment) <= 1e-8 * Math.Max(1.0, Math.Abs(moment)), $"ln E[v^{k}] {read:R} against {moment:R}");
        }
        output.WriteLine($"GBP({p}, {a}, {b}): tails and levels within {worst:G3}; mean {gbp.Mean:G6}, variance {gbp.Variance:G6}");
    }

    [Fact]
    public void Sweep_KeepsItsTailsQuantilesAndDensity()
    {
        var rng = new Xoshiro256PlusPlus(20261121);
        var failures = new List<string>();
        for (int c = 0; c < 400; c++)
        {
            double p = DistributionSweepTests.LogUniform(rng, 0.02, 50.0), a = DistributionSweepTests.LogUniform(rng, 1e-3, 1e3),
                   b = DistributionSweepTests.LogUniform(rng, 1e-3, 1e3);
            var gbp = GeneralizedBetaPrime.Of(p, a, b);
            var beta = Beta.Shapes(a, b);
            string at = $"GBP({p:R}, {a:R}, {b:R})";
            foreach (double logLevel in new[] { -600.0, -40.0, -5.0, -1.0, -0.05, -1e-6 })
            {
                double below = gbp.InverseLogCdf(logLevel), above = gbp.InverseLogSurvival(logLevel);
                // The tails meet the levels where the quantile is a finite, positive double; at the ends of the double range they need not.
                if (below > 0.0 && double.IsFinite(below) && below > 1e-300 && below < 1e300)
                {
                    double back = gbp.LogCdf(below);
                    if (Math.Abs(back - logLevel) > 1e-8 * Math.Max(1.0, Math.Abs(logLevel))) failures.Add($"{at}: ln F(F⁻¹({logLevel})) = {back:R}");
                }
                if (above > 0.0 && double.IsFinite(above) && above > 1e-300 && above < 1e300)
                {
                    double back = gbp.LogSurvival(above);
                    if (Math.Abs(back - logLevel) > 1e-8 * Math.Max(1.0, Math.Abs(logLevel))) failures.Add($"{at}: ln S(S⁻¹({logLevel})) = {back:R}");
                }
            }
            double x = Math.Exp(gbp.ModeInLog + 0.7 / p);
            if (double.IsFinite(x) && x > 0.0)
            {
                // The complement and the beta's density through y = x^p/(1 + x^p), dy/dx = p·x^(p−1)/(1 + x^p)².
                double f = gbp.LogCdf(x), s = gbp.LogSurvival(x), sum = Math.Max(f, s) + SpecialFunctions.Log1p(Math.Exp(-Math.Abs(f - s)));
                if (Math.Abs(sum) > 1e-12) failures.Add($"{at}: ln(F + S) = {sum:R} at {x:R}");
                double t = p * Math.Log(x), y = 1.0 / (1.0 + Math.Exp(-t));
                double viaBeta = beta.LogPdf(y) + Math.Log(p) + (p - 1.0) * Math.Log(x) - 2.0 * SpecialFunctions.Log1p(Math.Exp(t));
                if (y > 1e-12 && y < 1.0 - 1e-12 && Math.Abs(gbp.LogPdf(x) - viaBeta) > 1e-9 * Math.Max(1.0, Math.Abs(viaBeta)))
                    failures.Add($"{at}: ln f({x:R}) = {gbp.LogPdf(x):R} against the beta's {viaBeta:R}");
            }
        }
        foreach (string f in failures.Take(20)) output.WriteLine(f);
        Assert.True(failures.Count == 0, $"{failures.Count} failures; first: {failures.FirstOrDefault()}");
    }

    [Fact]
    public void Draws_FollowTheLaw_WholeAndOnAnInterval()
    {
        var rng = new Xoshiro256PlusPlus(20261122);
        var failures = new List<string>();
        var cases = new List<(double P, double A, double B)> { (1.0, 0.5, 0.5), (0.5, 1.0, 0.5), (2.0, 1.5, 3.0), (0.3, 4.0, 0.7), (5.0, 0.2, 2.0) };
        for (int c = 0; c < 60; c++)
            cases.Add((DistributionSweepTests.LogUniform(rng, 0.05, 20.0), DistributionSweepTests.LogUniform(rng, 0.05, 50.0),
                       DistributionSweepTests.LogUniform(rng, 0.05, 50.0)));
        foreach ((double p, double a, double b) in cases)
        {
            var gbp = GeneralizedBetaPrime.Of(p, a, b);
            // A law reaching past the double range at either end draws 0 or ∞ there, which no CDF of a double can follow.
            if (!(gbp.Quantile(1e-5) > 1e-300 && gbp.Quantile(1.0 - 1e-5) < 1e300)) continue;
            double[] xs = [.. Enumerable.Range(0, 2_000).Select(_ => gbp.Sample(rng))];
            double whole = Kolmogorov(xs, gbp.Cdf);
            if (whole > KolmogorovCut) failures.Add($"√n·D = {whole:F2}: GBP({p:R}, {a:R}, {b:R})");
            // On the interval between two quantiles, one of them far out.
            double lower = gbp.Quantile(0.3), upper = gbp.Quantile(1.0 - 1e-6);
            double fl = gbp.Cdf(lower), mass = gbp.IntervalProbability(lower, upper);
            double[] on = [.. Enumerable.Range(0, 2_000).Select(_ => gbp.SampleOn(rng, lower, upper))];
            if (on.Any(x => x < lower || x > upper)) failures.Add($"a draw off [{lower:R}, {upper:R}]: GBP({p:R}, {a:R}, {b:R})");
            double truncated = Kolmogorov(on, x => (gbp.Cdf(x) - fl) / mass);
            if (truncated > KolmogorovCut) failures.Add($"√n·D = {truncated:F2} on an interval: GBP({p:R}, {a:R}, {b:R})");
        }
        foreach (string f in failures) output.WriteLine(f);
        Assert.Empty(failures);
    }

    // ---- the mixtures ---------------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData(-0.5, 4.0, 1.0, 0.0, 0.8)]
    [InlineData(-0.5, 0.09, 1.0, 1.2, 0.5)]   // a strong drift and a heavy tail on its side: its table took 20 s before the integral about its peak
    [InlineData(-2.5, 0.0, 5.0, -0.7, 1.0)]
    [InlineData(0.7, 1.2, 0.8, 0.3, 0.7)]
    [InlineData(2.0, 3.0, 0.0, -0.4, 1.0)]
    public void TheIntegratedDensity_MeetsTheGeneralizedHyperbolicsClosedForm(double p, double a, double b, double drift, double scale)
    {
        var gh = GeneralizedHyperbolic.Of(p, a, b, drift, scale);
        double worst = 0.0;
        foreach (double level in new[] { 1e-6, 0.001, 0.1, 0.4, 0.6, 0.9, 0.999, 1.0 - 1e-6 })
        {
            double x = gh.Quantile(level), closed = gh.LogPdf(x);
            worst = Math.Max(worst, Math.Abs(gh.IntegratedLogPdf(x) - closed) / Math.Max(1.0, Math.Abs(closed)));
        }
        output.WriteLine($"GH({p}, {a}, {b}), m {drift}, s {scale}: the density's integral against the closed form, worst {worst:G3}");
        Assert.True(worst <= 1e-9, $"worst {worst:G3}");
    }

    // e^x·E₁(x), x > 0: the series for x ≤ 1, Lentz's continued fraction above (the exponential integral of order 1).
    private static double ScaledE1(double x)
    {
        if (x <= 1.0)
        {
            double sum = 0.0, term = 1.0;
            for (int k = 1; k < 60; k++)
            {
                term *= -x / k;
                sum += term / k;
            }
            return Math.Exp(x) * (-0.57721566490153286061 - Math.Log(x) - sum);
        }
        double tiny = 1e-300, f = tiny, c = f, d = 0.0;
        // e^x·E₁(x) = 1/(x + 1 − 1²/(x + 3 − 2²/(x + 5 − …))).
        for (int k = 0; k < 10_000; k++)
        {
            double an = k == 0 ? 1.0 : -(double)k * k, bn = x + 2.0 * k + 1.0;
            d = bn + an * d;
            d = Math.Abs(d) < tiny ? tiny : d;
            c = bn + an / c;
            c = Math.Abs(c) < tiny ? tiny : c;
            d = 1.0 / d;
            double delta = c * d;
            f *= delta;
            if (Math.Abs(delta - 1.0) < 1e-16) break;
        }
        return f;
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(0.3)]
    [InlineData(7.0)]
    public void TheHorseshoesDensity_IsTheExponentialIntegrals(double scale)
    {
        var horseshoe = NormalGeneralizedBetaPrime.Of(1.0, 0.5, 0.5, 0.0, scale);
        double worst = 0.0;
        for (int k = -16; k <= 16; k++)
            foreach (double sign in new[] { -1.0, 1.0 })
            {
                double x = sign * scale * Math.Pow(10.0, 0.5 * k), c = 0.5 * (x / scale) * (x / scale);
                double exact = -0.5 * Math.Log(2.0 * Math.Pow(Math.PI, 3.0)) + Math.Log(ScaledE1(c)) - Math.Log(scale);
                worst = Math.Max(worst, Math.Abs(horseshoe.LogPdf(x) - exact) / Math.Max(1.0, Math.Abs(exact)));
            }
        output.WriteLine($"horseshoe, s {scale}: the density's integral against (2π³)^(−½)e^c·E₁(c)/s, worst relative gap {worst:G3}");
        Assert.True(worst <= 1e-9, $"worst {worst:G3}");
        Assert.Equal(double.PositiveInfinity, horseshoe.LogPdf(0.0));
    }

    public static IEnumerable<object[]> Mixtures() =>
    [
        [1.0, 0.5, 0.5, 0.0, 1.0],     // the horseshoe, its pole at the median
        [1.0, 0.5, 0.5, 0.6, 0.8],     // skewed: the walk crosses the pole
        [0.5, 1.0, 0.5, 0.0, 1.3],     // the Strawderman–Berger
        [0.5, 1.0, 0.5, -0.4, 1.0],
        [2.0, 1.5, 3.0, 0.3, 0.5],     // no pole, finite variance
        [1.0, 0.3, 2.0, 0.5, 1.0],     // a pole, with a drift
    ];

    // ln ∫ exp(h(u)) du over [−150, 250] by Simpson's rule at step 10⁻³, refined to a fiftieth of s/√(xm) within 40 of it about
    // u* = ln(x/m), where z vanishes and, far out on the drift's side, the integrand narrows to that width.
    private static double Reference(Func<double, double> h, double x, double drift, double scale)
    {
        const double lo = -150.0, hi = 250.0, step = 1e-3;
        if (!(x * drift > 0.0)) return Simpson(h, lo, hi, step);
        double centre = Math.Log(x / drift), width = scale / Math.Sqrt(x * drift);
        if (width >= 0.05) return Simpson(h, lo, hi, step);
        double a = centre - 40.0 * width, b = centre + 40.0 * width;
        double[] parts = [Simpson(h, lo, a, step), Simpson(h, a, b, width / 50.0), Simpson(h, b, hi, step)];
        double top = parts.Max();
        return top + Math.Log(parts.Sum(v => Math.Exp(v - top)));
    }

    [Fact]
    public void PastTheTablesEnds_TheSideNearOne_KeepsTheFarTailsDigits()
    {
        // Past an end of the tail table the side near 1 is read as ln(1 − e^(far tail)): never above 0, and 1 − F the far tail's own integral.
        // Before, it was its own integral, which lost 1 − F's digits (+2.9·10⁻¹⁵ for the generalized hyperbolic at 32 scales).
        var laws = new (string Name, NormalVarianceMeanMixture Law)[]
        {
            ("generalized hyperbolic", GeneralizedHyperbolic.Of(0.7, 1.2, 0.8, 0.5, 1.0)),
            ("Normal-Wald, γ = 1.5", GeneralizedHyperbolic.Of(-0.5, 2.25, 1.0, 0.4, 0.8)),
            ("skewed t, ν = 5", GeneralizedHyperbolic.Of(-2.5, 0.0, 5.0, 0.3, 1.0)),
            ("horseshoe, skewed", NormalGeneralizedBetaPrime.Of(1.0, 0.5, 0.5, 0.3, 1.0)),
        };
        foreach ((string name, NormalVarianceMeanMixture law) in laws)
        {
            double worst = 0.0, integratedAbove = double.NegativeInfinity;
            foreach (double side in new[] { -1.0, 1.0 })
                foreach (double scales in new[] { 8.0, 16.0, 32.0, 64.0, 1e3 })
                {
                    double x = side * scales * law.Scale, f = law.LogCdf(x), s = law.LogSurvival(x);
                    Assert.True(f <= 0.0 && s <= 0.0, $"{name} at {x}: ln F {f:R}, ln S {s:R}");
                    // The side near 1 against the far tail's own mass: near = ln(1 − e^(far)), measured against e^(far) where that is normal.
                    double near = side > 0.0 ? f : s, far = side > 0.0 ? s : f;
                    if (far > -700.0) worst = Math.Max(worst, Math.Abs(near - SpecialFunctions.Log1p(-Math.Exp(far))) / Math.Exp(far));
                    integratedAbove = Math.Max(integratedAbove, side > 0.0 ? law.IntegratedLogCdf(x) : law.IntegratedLogSurvival(x));
                }
            output.WriteLine($"{name}: the near side against ln(1 − e^(far tail)), worst gap over e^(far) {worst:E2}; the near side's own integral, highest {integratedAbove:E2}");
            Assert.True(worst < 1e-9, $"{name}: {worst:E2}");
        }
    }

    [Fact]
    public void NearZero_UnderADrift_TheTableHoldsItsIntegral()
    {
        // Where v's density near 0 goes as v^(α−1), α finite, the mixture's density goes as |x|^(2α−1) at 0: a pole at α ≤ ½, a cusp or kink
        // above it, where the normal score is not smooth. A drift moves the median off 0, and a cubic cell across 0 missed the integral by
        // 10⁻³ of ln F (GBP(1, 1, 1.5) at m/s = −0.2, the conformal chip's find). A seeded sweep over drifts and shapes, the GBP's and the
        // variance gamma's, reads both log tails at 0, beside it and out to s against the integral.
        var rng = new Xoshiro256PlusPlus(20261108);
        double[] offsets = [0.0, 1e-9, -1e-9, 1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3, 1e-2, -1e-2, 0.03, -0.03, 0.1, -0.1, 0.3, -0.3, 0.6, -0.6, 1.0, -1.0];
        double worst = 0.0;
        string worstAt = "";
        for (int trial = 0; trial < 40; trial++)
        {
            double scale = Math.Exp(-2.0 + 4.0 * rng.NextDouble()), drift = scale * (-1.0 + 2.0 * rng.NextDouble());
            NormalVarianceMeanMixture law;
            string name;
            if (trial % 4 == 3)
            {
                // The variance gamma, v ~ Gamma(p, a/2): α = p.
                double shape = 0.55 + 2.5 * rng.NextDouble(), a = 0.5 + 2.0 * rng.NextDouble();
                law = GeneralizedHyperbolic.Of(shape, a, 0.0, drift, scale);
                name = $"VG({shape:F3}, {a:F3})";
            }
            else
            {
                // The GBP, α = p·a from just above the pole's ½ to 3.
                double p = 0.5 + 1.5 * rng.NextDouble(), alpha = 0.52 + 2.5 * rng.NextDouble(), b = 0.4 + 2.5 * rng.NextDouble();
                law = NormalGeneralizedBetaPrime.Of(p, alpha / p, b, drift, scale);
                name = $"GBP({p:F3}, {alpha / p:F3}, {b:F3})";
            }
            double lawWorst = 0.0;
            foreach (double offset in offsets.Append(law.Quantile(0.5) / scale))
            {
                double x = offset * scale, below = law.IntegratedLogCdf(x), above = law.IntegratedLogSurvival(x);
                double gap = Math.Max(Math.Abs(law.LogCdf(x) - below) / Math.Max(1.0, Math.Abs(below)),
                                      Math.Abs(law.LogSurvival(x) - above) / Math.Max(1.0, Math.Abs(above)));
                lawWorst = Math.Max(lawWorst, gap);
                if (gap > worst) (worst, worstAt) = (gap, $"{name}, m/s {drift / scale:F3}, x/s {offset:G3}");
            }
            output.WriteLine($"{name}, m/s {drift / scale,6:F3}, s {scale:F3}: worst gap {lawWorst:E2}");
        }
        output.WriteLine($"The worst relative gap in the log tails, at 0 and out to s: {worst:E2} ({worstAt})");
        Assert.True(worst <= 1e-5, $"{worst:E2} at {worstAt}");
    }

    [Theory]
    [MemberData(nameof(Mixtures))]
    public void DensityAndTails_MeetAFixedStepRule(double p, double a, double b, double drift, double scale)
    {
        var noise = NormalGeneralizedBetaPrime.Of(p, a, b, drift, scale);
        GeneralizedBetaPrime g = noise.Mixing;
        double logRoot = 0.5 * Math.Log(2.0 * Math.PI) + Math.Log(scale);
        double Z(double x, double u) => (x * Math.Exp(-0.5 * u) - drift * Math.Exp(0.5 * u)) / scale;
        double density = 0.0, tails = 0.0, table = 0.0, levels = 0.0;
        foreach (double level in new[] { 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999 })
        {
            double x = noise.Quantile(level);
            if (x == 0.0) continue;
            double f = Reference(u => -0.5 * Math.Pow(Z(x, u), 2.0) - 0.5 * u - logRoot + g.LogDensityInLog(u), x, drift, scale);
            double below = Reference(u => Normal.StandardLogCdf(Z(x, u)) + g.LogDensityInLog(u), x, drift, scale);
            double above = Reference(u => Normal.StandardLogCdf(-Z(x, u)) + g.LogDensityInLog(u), x, drift, scale);
            density = Math.Max(density, Math.Abs(noise.LogPdf(x) - f) / Math.Max(1.0, Math.Abs(f)));
            tails = Math.Max(tails, Math.Abs(noise.IntegratedLogCdf(x) - below) / Math.Max(1.0, Math.Abs(below)));
            tails = Math.Max(tails, Math.Abs(noise.IntegratedLogSurvival(x) - above) / Math.Max(1.0, Math.Abs(above)));
            levels = Math.Max(levels, Math.Abs((level <= 0.5 ? Math.Exp(below) : -SpecialFunctions.Expm1(above)) - level));
        }
        double low = noise.Quantile(1e-9), high = noise.Quantile(1.0 - 1e-9);
        var grid = Enumerable.Range(0, 81).Select(k => low + (high - low) * k / 80.0).Concat([0.0, 1e-12, -1e-12, 1e-6, -1e-6, 1e-3, -1e-3]);
        double worstAt = double.NaN;
        foreach (double x in grid)
        {
            double below = noise.IntegratedLogCdf(x), above = noise.IntegratedLogSurvival(x);
            double gap = Math.Max(Math.Abs(noise.LogCdf(x) - below) / Math.Max(1.0, Math.Abs(below)),
                                  Math.Abs(noise.LogSurvival(x) - above) / Math.Max(1.0, Math.Abs(above)));
            if (gap > table) (table, worstAt) = (gap, x);
        }
        output.WriteLine($"GBP({p}, {a}, {b}), m {drift}, s {scale}: density {density:G3}, tails {tails:G3}, table {table:G3} at {worstAt:G4}, " +
                         $"levels {levels:G3}");
        Assert.True(density <= 1e-8, $"the density against the trapezoid, worst {density:G3}");
        Assert.True(tails <= 1e-8, $"the tails against the trapezoid, worst {tails:G3}");
        Assert.True(table <= 5e-4, $"the table against the integral, worst relative gap in the log tails {table:G3}");
        Assert.True(levels <= 1e-4, $"the table's quantiles against the integral's mass, worst {levels:G3}");
    }

    [Theory]
    [MemberData(nameof(Mixtures))]
    public void TheMixturesDraws_FollowItsLaw(double p, double a, double b, double drift, double scale)
    {
        var noise = NormalGeneralizedBetaPrime.Of(p, a, b, drift, scale);
        var rng = new Xoshiro256PlusPlus(20261123);
        // The table's build, on the first read of a tail: a predictive builds one per draw.
        var watch = System.Diagnostics.Stopwatch.StartNew();
        noise.LogCdf(0.0);
        double build = watch.Elapsed.TotalMilliseconds;
        double[] xs = [.. Enumerable.Range(0, 4_000).Select(_ => noise.Sample(rng))];
        double statistic = Kolmogorov(xs, noise.Cdf);
        double lower = noise.Quantile(0.2), upper = noise.Quantile(0.95), fl = noise.Cdf(lower), mass = noise.IntervalProbability(lower, upper);
        double[] on = [.. Enumerable.Range(0, 2_000).Select(_ => noise.SampleOn(rng, lower, upper))];
        double truncated = Kolmogorov(on, x => (noise.Cdf(x) - fl) / mass);
        output.WriteLine($"GBP({p}, {a}, {b}), m {drift}, s {scale}: √n·D {statistic:F2} whole, {truncated:F2} on [{lower:G4}, {upper:G4}]; " +
                         $"the table built in {build:F0} ms (wall clock)");
        Assert.True(statistic <= KolmogorovCut && truncated <= KolmogorovCut);
        Assert.All(on, x => Assert.InRange(x, lower, upper));
    }
}
