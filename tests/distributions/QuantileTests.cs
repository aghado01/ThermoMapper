using System;
using System.Linq;
using Maths.Distributions;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The quantiles (<see cref="IContinuousDistribution.Quantile"/> and the log forms), against references that share no code
/// with them: closed forms where a shape admits one, far into both tails, and round trips through the log tails they invert,
/// as <see cref="NormalQuantileTests"/> checks Φ⁻¹. R's <c>qgamma</c>, <c>qbeta</c> and <c>qt</c> are in the oracle project.
/// </summary>
public sealed class QuantileTests
{
    private const double Epsilon = 2.220446049250313e-16;
    private const double Ln2 = 0.69314718055994530942;

    // The smaller tail's log probability, from the median to far below the double range.
    private static readonly double[] Tails = [-Ln2, Math.Log(0.3), Math.Log(0.1), Math.Log(1e-3), Math.Log(1e-10), Math.Log(1e-100), -700.0, -5000.0];

    private static double Log1mExp(double u) => u > -Ln2 ? Math.Log(-Expm1(u)) : Log1p(-Math.Exp(u));

    private static double Log1p(double x)
    {
        double u = 1.0 + x;
        return u == 1.0 ? x : Math.Log(u) * x / (u - 1.0);
    }

    private static double Expm1(double x)
    {
        double u = Math.Exp(x);
        if (u == 1.0) return x;
        if (u - 1.0 == -1.0) return -1.0;
        return (u - 1.0) * x / Math.Log(u);
    }

    private static void Near(double expected, double actual, double relative, string what)
        => Assert.True(actual == expected || Math.Abs(actual - expected) <= relative * Math.Abs(expected), $"{what}: {actual:R} against {expected:R}");

    // A log tail's error is its prefactor's, about ε·S with S the size of its terms, on top of x's own rounding.
    private static void RoundTrip(double target, double actual, double prefactor, string what)
    {
        double tolerance = 2e-13 * Math.Max(1.0, Math.Abs(target)) + 8.0 * Epsilon * prefactor;
        Assert.True(Math.Abs(actual - target) <= tolerance, $"{what}: {actual:R} against {target:R}, tolerance {tolerance:E1}");
    }

    [Theory]
    [InlineData(1e-3)]
    [InlineData(0.01)]
    [InlineData(0.1)]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(2.5)]
    [InlineData(10.0)]
    [InlineData(100.0)]
    [InlineData(1000.0)]
    [InlineData(1e5)]
    public void Gamma_InvertsItsLogTails(double shape)
    {
        foreach (double rate in new[] { 1.0, 7.5 })
        {
            Gamma d = Gamma.ShapeRate(shape, rate);
            double logGamma = SpecialFunctions.LogGamma(shape);
            foreach (double tail in Tails)
            {
                double x = d.InverseLogCdf(tail);
                if (x > 1e-300)
                    RoundTrip(tail, d.LogCdf(x), Math.Abs(shape * Math.Log(rate * x)) + rate * x + Math.Abs(logGamma), $"Gamma({shape}, {rate}) lower {tail}");
                else
                    Assert.True((tail + logGamma + Math.Log(shape)) / shape < -600.0, $"Gamma({shape}) lower {tail} underflowed to {x}");

                double above = d.InverseLogSurvival(tail);
                RoundTrip(tail, d.LogSurvival(above), Math.Abs(shape * Math.Log(rate * above)) + rate * above + Math.Abs(logGamma), $"Gamma({shape}, {rate}) upper {tail}");
            }
        }
    }

    [Theory]
    [InlineData(1e-3)]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(4.0)]
    [InlineData(50.0)]
    [InlineData(1000.0)]
    public void InverseGamma_InvertsItsLogTails(double shape)
    {
        foreach (double scale in new[] { 1.0, 0.03 })
        {
            InverseGamma d = InverseGamma.ShapeScale(shape, scale);
            double logGamma = SpecialFunctions.LogGamma(shape);
            foreach (double tail in Tails)
            {
                double x = d.InverseLogCdf(tail), g = scale / x;
                RoundTrip(tail, d.LogCdf(x), Math.Abs(shape * Math.Log(g)) + g + Math.Abs(logGamma), $"IG({shape}, {scale}) lower {tail}");

                double above = d.InverseLogSurvival(tail);
                if (double.IsPositiveInfinity(above))
                    Assert.True((tail + logGamma + Math.Log(shape)) / shape < -600.0, $"IG({shape}) upper {tail} overflowed");
                else
                    RoundTrip(tail, d.LogSurvival(above), Math.Abs(shape * Math.Log(scale / above)) + scale / above + Math.Abs(logGamma), $"IG({shape}, {scale}) upper {tail}");
            }
        }
    }

    [Theory]
    [InlineData(0.01, 0.01)]
    [InlineData(0.1, 5.0)]
    [InlineData(0.5, 0.5)]
    [InlineData(2.0, 5.0)]
    [InlineData(5.0, 0.1)]
    [InlineData(20.0, 30.0)]
    [InlineData(100.0, 1000.0)]
    [InlineData(1000.0, 1000.0)]
    public void Beta_InvertsItsLogTails_AndAgreesWithItsReflection(double a, double b)
    {
        Beta d = Beta.Shapes(a, b), reflected = Beta.Shapes(b, a);
        double logBeta = SpecialFunctions.LogBeta(a, b);
        double Prefactor(double x) => Math.Abs(logBeta) + Math.Abs(a * Math.Log(x)) + Math.Abs(b * Log1p(-x));
        foreach (double tail in Tails)
        {
            // X and 1 − X read the same curve. A U-shaped law's flat centre is the worst conditioned, at about 10⁻¹³.
            double x = d.InverseLogCdf(tail), y = reflected.InverseLogSurvival(tail);
            Assert.True(Math.Abs(x - (1.0 - y)) <= 1e-12, $"Beta({a}, {b}) lower {tail}: {x:R} against 1 − {y:R}");
            if (x > 1e-300 && 1.0 - x > 1e-3)
                RoundTrip(tail, d.LogCdf(x), Prefactor(x), $"Beta({a}, {b}) lower {tail}");

            // The upper tail keeps its digits in the reflection, whose x is the complement.
            double complement = reflected.InverseLogCdf(tail);
            if (complement > 1e-300 && 1.0 - complement > 1e-3)
                RoundTrip(tail, reflected.LogCdf(complement), Prefactor(1.0 - complement), $"Beta({a}, {b}) upper {tail}");
            Assert.True(Math.Abs(1.0 - complement - d.InverseLogSurvival(tail)) <= 1e-12, $"Beta({a}, {b}) upper {tail}: {d.InverseLogSurvival(tail):R} against 1 − {complement:R}");
        }
    }

    [Theory]
    [InlineData(0.1)]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(2.0)]
    [InlineData(3.0)]
    [InlineData(10.0)]
    [InlineData(30.0)]
    [InlineData(1000.0)]
    public void StudentT_InvertsItsLogTails_AndIsSymmetric(double nu)
    {
        StudentT d = StudentT.LocationScale(nu, 2.0, 3.0), standard = StudentT.Standard(nu);
        double prefactor = Math.Abs(SpecialFunctions.LogBeta(0.5 * nu, 0.5));
        foreach (double tail in Tails)
        {
            double t = standard.InverseLogSurvival(tail);
            Assert.Equal(-t, standard.InverseLogCdf(tail));
            Near(2.0 + 3.0 * t, d.InverseLogSurvival(tail), 1e-15, $"t({nu}) location-scale {tail}");
            if (double.IsPositiveInfinity(t))
            {
                Assert.True(tail < -300.0 * nu, $"t({nu}) overflowed at {tail}");
                continue;
            }
            RoundTrip(tail, standard.LogSurvival(t), 0.5 * nu * Math.Abs(tail) + prefactor, $"t({nu}) upper {tail}");
            RoundTrip(tail, standard.LogCdf(-t), 0.5 * nu * Math.Abs(tail) + prefactor, $"t({nu}) lower {tail}");
        }
    }

    [Fact]
    public void ClosedForms_HoldFarIntoBothTails()
    {
        foreach (double tail in new[] { -Ln2, -1.0, Math.Log(1e-10), -100.0, -700.0, -5000.0, -1e6, -1e14, -1e16, -1e100, -1e300 })
        {
            // Gamma(1, β) is the exponential: S = e^(−βx).
            Gamma exponential = Gamma.ShapeRate(1.0, 2.5);
            Near(-tail / 2.5, exponential.InverseLogSurvival(tail), 1e-13, $"exponential S {tail}");
            Near(-Log1mExp(tail) / 2.5, exponential.InverseLogCdf(tail), 1e-13, $"exponential F {tail}");

            // IG(1, β) is its reciprocal: F = e^(−β/x).
            InverseGamma reciprocal = InverseGamma.ShapeScale(1.0, 2.5);
            Near(-2.5 / tail, reciprocal.InverseLogCdf(tail), 1e-13, $"IG(1) F {tail}");
            Near(-2.5 / Log1mExp(tail), reciprocal.InverseLogSurvival(tail), 1e-13, $"IG(1) S {tail}");

            // Beta(a, 1): F = xᵃ. Beta(1, b): S = (1 − x)ᵇ.
            Near(Math.Exp(tail / 3.0), Beta.Shapes(3.0, 1.0).InverseLogCdf(tail), 1e-13, $"Beta(3, 1) F {tail}");
            Near(Math.Exp(Log1mExp(tail) / 3.0), Beta.Shapes(3.0, 1.0).InverseLogSurvival(tail), 1e-13, $"Beta(3, 1) S {tail}");
            Near(-Expm1(tail / 3.0), Beta.Shapes(1.0, 3.0).InverseLogSurvival(tail), 1e-13, $"Beta(1, 3) S {tail}");
        }

        foreach (double tail in new[] { -1.0, Math.Log(1e-10), -100.0, -700.0 })
        {
            // A solve in ln S carries about ε·|ln S| of rounding into ln t.
            double reach = 1e-13 + 4.0 * Epsilon * Math.Abs(tail);
            // The Cauchy: t = −cot(πp) in the lower tail.
            Near(1.0 / Math.Tan(Math.PI * Math.Exp(tail)), StudentT.Standard(1.0).InverseLogSurvival(tail), reach, $"Cauchy {tail}");
            // ν = 2: t = (2p − 1)/√(2p(1 − p)), so far out t = √(1/(2p))·(1 − 2p)/√(1 − p).
            double p = Math.Exp(tail);
            Near((1.0 - 2.0 * p) / Math.Sqrt(2.0 * p * (1.0 - p)), StudentT.Standard(2.0).InverseLogSurvival(tail), reach, $"t(2) {tail}");
            // Beta(½, ½), the arcsine: x = sin²(πp/2).
            double sine = Math.Sin(0.5 * Math.PI * p);
            Near(sine * sine, Beta.Shapes(0.5, 0.5).InverseLogCdf(tail), 1e-13, $"arcsine {tail}");
        }

        foreach (double tail in new[] { -1.0, Math.Log(1e-10), -100.0, -700.0 })
        {
            // Gamma(½, rate ½) is χ²₁, the square of a normal: S(x) = 2Φ(−√x).
            double z = Normal.StandardQuantile(0.5 * Math.Exp(tail));
            Near(z * z, Gamma.ShapeRate(0.5, 0.5).InverseLogSurvival(tail), 1e-13, $"chi-square S {tail}");
            // Gamma(3): S = e^(−g)(1 + g + g²/2), with no incomplete gamma in it.
            double g = Gamma.ShapeRate(3.0, 1.0).InverseLogSurvival(tail);
            Near(tail, -g + Math.Log(1.0 + g + 0.5 * g * g), 1e-13, $"Gamma(3) S {tail}");
        }

        // Past the double range, the tail is still its power law.
        Assert.Equal(0.0, Gamma.ShapeRate(1e-5, 1.0).Quantile(0.5));
        Near(Math.Exp(-5000.0 / 8.0), Gamma.ShapeRate(8.0, 1.0).InverseLogCdf(-5000.0 - SpecialFunctions.LogGamma(9.0)), 1e-11, "Gamma(8) power law");
        Assert.Equal(double.PositiveInfinity, StudentT.Standard(1.0).InverseLogSurvival(-1000.0));
        Assert.Equal(1.0, Beta.Shapes(1.0, 1e-3).Quantile(0.5));
        Near(Math.Exp(-1000.0 * Ln2), Beta.Shapes(1e-3, 1.0).Quantile(0.5), 1e-12, "Beta(1e-3, 1)");
    }

    [Fact]
    public void NearTheCentre_ASmallTKeepsRelativeDigits()
    {
        foreach (double offset in new[] { 1e-3, 1e-8, 1e-12 })
        {
            double p = 0.5 + offset, exact = p - 0.5;
            Near(Math.Tan(Math.PI * exact), StudentT.Standard(1.0).Quantile(p), 1e-14, $"Cauchy centre {offset}");
            Near(2.0 * exact / Math.Sqrt(2.0 * p * (1.0 - p)), StudentT.Standard(2.0).Quantile(p), 1e-14, $"t(2) centre {offset}");
            Near(-Math.Tan(Math.PI * exact), StudentT.Standard(1.0).Quantile(1.0 - p), 1e-14, $"Cauchy centre below {offset}");
        }
        Assert.Equal(4.0, StudentT.LocationScale(5.0, 4.0, 2.0).Quantile(0.5));
    }

    [Fact]
    public void Quantile_IsMonotone_AndReadsTheSameCurveAsTheLogForms()
    {
        IContinuousDistribution[] family =
        [
            Gamma.ShapeRate(0.3, 2.0), Gamma.ShapeRate(40.0, 1.0), InverseGamma.ShapeScale(3.0, 2.0), Beta.Shapes(0.4, 3.0), Beta.Shapes(30.0, 40.0),
            StudentT.LocationScale(3.0, -1.0, 2.0), Normal.MeanSd(1.0, 2.0), Weibull.ShapeScale(0.5, 2.0), HalfCauchy.WithScale(2.0), Uniform.On(-1.0, 3.0),
        ];
        foreach (IContinuousDistribution d in family)
        {
            double previous = double.NegativeInfinity;
            foreach (double p in new[] { 1e-200, 1e-30, 1e-6, 0.01, 0.2, 0.5, 0.8, 0.99, 1.0 - 1e-6, 1.0 - 1e-12 })
            {
                double x = d.Quantile(p);
                Assert.True(x >= previous, $"{d.GetType().Name} at {p}: {x:R} after {previous:R}");
                previous = x;
                Near(x, d.InverseLogCdf(Math.Log(p)), 1e-9, $"{d.GetType().Name} log CDF at {p}");
                Near(x, d.InverseLogSurvival(Log1p(-p)), 1e-9, $"{d.GetType().Name} log survival at {p}");
                Assert.True(Math.Abs(d.Cdf(x) - p) <= 1e-9 * p + 1e-13, $"{d.GetType().Name} CDF at {p}: {d.Cdf(x):R}");
            }
        }
    }

    [Fact]
    public void WeibullAndHalfCauchy_LogFormsMatchTheirClosedForms()
    {
        Weibull w = Weibull.ShapeScale(0.7, 3.0);
        HalfCauchy h = HalfCauchy.WithScale(1e-5);
        foreach (double tail in new[] { -Ln2, -1.0, -50.0, -700.0, -5000.0 })
        {
            Near(tail, w.LogSurvival(w.InverseLogSurvival(tail)), 1e-13, $"Weibull S {tail}");
            double below = w.InverseLogCdf(tail);
            if (below > 0.0) RoundTrip(tail, w.LogCdf(below), 0.0, $"Weibull F {tail}");
            else Assert.True(tail / 0.7 < -700.0, $"Weibull F {tail} underflowed");
            double lowest = h.InverseLogCdf(tail);
            if (lowest > 0.0) RoundTrip(tail, h.LogCdf(lowest), 0.0, $"half-Cauchy F {tail}");
            else Assert.True(tail < -690.0, $"half-Cauchy F {tail} underflowed");
            double x = h.InverseLogSurvival(tail);
            if (tail >= -700.0) RoundTrip(tail, h.LogSurvival(x), 0.0, $"half-Cauchy S {tail}");
            else Assert.True(double.IsInfinity(x) || x > 1e300 * 1e-5, $"half-Cauchy S {tail}: {x:R}");
        }
    }

    [Fact]
    public void Truncated_ReadsItsParentsQuantiles_InTheFarTail()
    {
        // N(0, 1) on [10, ∞): 1 − Φ(10) is below ε, so p·mass only lives in logs.
        var tail = Truncated.Of(Normal.Standard, 10.0, double.PositiveInfinity);
        foreach (double p in new[] { 1e-6, 0.01, 0.5, 0.9, 1.0 - 1e-9 })
        {
            double z = tail.Quantile(p);
            Assert.InRange(z, 10.0, 15.0);
            Assert.Equal(p, tail.Cdf(z), 1e-11 * Math.Max(1.0, 1.0 / (1.0 - p)));
        }
        Assert.Equal(10.0, tail.Quantile(0.0));
        Assert.Equal(double.PositiveInfinity, tail.Quantile(1.0));
        // Median of the far tail: Q(z) = Q(10)/2, so z ≈ 10 + ln 2/10.
        Assert.Equal(10.0 + Ln2 / 10.0, tail.Quantile(0.5), 2e-3);

        foreach (var truncated in new IContinuousDistribution[]
                 {
                     Truncated.Of(Gamma.ShapeRate(2.0, 1.0), 1000.0, 1010.0), Truncated.Of(Gamma.ShapeRate(1000.0, 1.0), 0.0, 1.0),
                     Truncated.Of(InverseGamma.ShapeScale(50.0, 1e-3), 1.0, double.PositiveInfinity), Truncated.Of(StudentT.Standard(3.0), 1e6, double.PositiveInfinity),
                     Truncated.Of(Beta.Shapes(2.0, 5.0), 0.999, 1.0), Truncated.Of(Normal.Standard, -1.0, 2.0), Truncated.Of(Normal.Standard, -40.0, -39.0),
                     Truncated.Of(Normal.Standard, 10.0, 10.0001), Truncated.Of(StudentT.LocationScale(3.0, 2.0, 0.5), 5.0, 5.01),
                 })
            foreach (double p in new[] { 1e-3, 0.1, 0.5, 0.9, 0.999 })
            {
                double x = truncated.Quantile(p);
                Assert.True(Math.Abs(truncated.Cdf(x) - p) <= 1e-9, $"{truncated.GetType().Name}: F({x:R}) = {truncated.Cdf(x):R} against {p}");
            }
    }

    [Fact]
    public void Mixture_QuantileLiesBetweenItsComponents_AndInvertsItsCdf()
    {
        Normal a = Normal.MeanSd(-3.0, 1.0), b = Normal.MeanSd(4.0, 0.5);
        Mixture mix = Mixture.Of([a, b], [1.0, 3.0]);
        foreach (double p in new[] { 1e-12, 0.01, 0.2, 0.25, 0.5, 0.9, 1.0 - 1e-9 })
        {
            double x = mix.Quantile(p);
            Assert.True(Math.Abs(mix.Cdf(x) - p) <= 1e-10, $"F({x:R}) = {mix.Cdf(x):R} against {p}");
            Assert.InRange(x, Math.Min(a.Quantile(p), b.Quantile(p)), Math.Max(a.Quantile(p), b.Quantile(p)));
        }
        // A gap between two truncated components: the bracket holds where the density vanishes.
        Mixture gapped = Mixture.Of([Truncated.Of(Normal.Standard, -1.0, 0.0), Truncated.Of(Normal.Standard, 5.0, 6.0)], [1.0, 1.0]);
        Assert.InRange(gapped.Quantile(0.5), 0.0, 5.0);
        Assert.Equal(0.25, gapped.Cdf(gapped.Quantile(0.25)), 1e-9);
        Assert.Equal(0.75, gapped.Cdf(gapped.Quantile(0.75)), 1e-9);
        // Truncating the mixture reads through both.
        var truncated = Truncated.Of(mix, 3.0, double.PositiveInfinity);
        Assert.Equal(0.3, truncated.Cdf(truncated.Quantile(0.3)), 1e-10);
    }

    [Fact]
    public void Ends_AndRefusals()
    {
        foreach (IContinuousDistribution d in new IContinuousDistribution[]
                 {
                     Gamma.ShapeRate(2.0, 1.0), InverseGamma.ShapeScale(2.0, 1.0), Beta.Shapes(2.0, 3.0), StudentT.Standard(3.0), Normal.Standard,
                     Weibull.ShapeScale(2.0, 1.0), HalfCauchy.WithScale(1.0), Uniform.On(0.0, 1.0), Truncated.Of(Normal.Standard, 0.0, 1.0), Mixture.Of([Normal.Standard], [1.0]),
                 })
        {
            string name = d.GetType().Name;
            foreach (double p in new[] { -0.1, 1.1, double.NaN })
                Assert.Throws<ArgumentOutOfRangeException>(() => d.Quantile(p));
            foreach (double logP in new[] { 0.1, double.NaN })
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => d.InverseLogCdf(logP));
                Assert.Throws<ArgumentOutOfRangeException>(() => d.InverseLogSurvival(logP));
            }
            Assert.True(d.Quantile(0.0) < d.Quantile(1.0), name);
            Assert.Equal(d.Quantile(0.0), d.InverseLogCdf(double.NegativeInfinity));
            Assert.Equal(d.Quantile(1.0), d.InverseLogCdf(0.0));
            Assert.Equal(d.Quantile(1.0), d.InverseLogSurvival(double.NegativeInfinity));
            Assert.Equal(d.Quantile(0.0), d.InverseLogSurvival(0.0));
        }
        Assert.Equal(0.0, Gamma.ShapeRate(2.0, 1.0).Quantile(0.0));
        Assert.Equal(double.PositiveInfinity, Gamma.ShapeRate(2.0, 1.0).Quantile(1.0));
        Assert.Equal(0.0, InverseGamma.ShapeScale(2.0, 1.0).Quantile(0.0));
        Assert.Equal(double.PositiveInfinity, InverseGamma.ShapeScale(2.0, 1.0).Quantile(1.0));
        Assert.Equal(0.0, Beta.Shapes(2.0, 3.0).Quantile(0.0));
        Assert.Equal(1.0, Beta.Shapes(2.0, 3.0).Quantile(1.0));
        Assert.Equal(double.NegativeInfinity, StudentT.Standard(3.0).Quantile(0.0));
        Assert.Equal(double.PositiveInfinity, StudentT.Standard(3.0).Quantile(1.0));
    }
}
