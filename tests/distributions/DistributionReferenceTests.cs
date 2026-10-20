using System;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// Densities, CDFs and interval probabilities against references computed independently with mpmath at 50–60
/// significant digits. Tail cases are checked in log space where the probability underflows, and the interval
/// cases include both same-tail differences and narrow intervals straddling the centre, the two places a naive
/// F(b) − F(a) loses its digits.
/// </summary>
public sealed class DistributionReferenceTests
{
    // Relative agreement for probabilities; absolute agreement, scaled by max(1, |reference|), for logs.
    private static void Close(double expected, double actual, double rel = 1e-13)
    {
        if (double.IsNegativeInfinity(expected)) { Assert.Equal(expected, actual); return; }
        if (expected == 0.0) { Assert.Equal(0.0, actual); return; }
        Assert.True(Math.Abs(actual - expected) <= rel * Math.Abs(expected), $"expected {expected:R}, got {actual:R}");
    }

    private static void CloseLog(double expected, double actual, double tol = 1e-12)
    {
        if (double.IsNegativeInfinity(expected)) { Assert.Equal(expected, actual); return; }
        Assert.True(Math.Abs(actual - expected) <= tol * Math.Max(1.0, Math.Abs(expected)), $"expected {expected:R}, got {actual:R}");
    }

    // The log prefactor of an incomplete gamma or beta is a sum of terms of size S, so any implementation through
    // lnΓ reaches at best about ε·S relative; the tolerance is a small multiple of that floor.
    private static double PrefactorTolerance(double magnitude) => Math.Max(1e-13, 8.0 * 2.220446049250313e-16 * magnitude);

    private static double GammaPrefactor(double a, double x)
        => Math.Abs(SpecialFunctions.LogGamma(a)) + Math.Abs(a * Math.Log(x)) + x;

    private static double BetaPrefactor(double x, double a, double b)
        => Math.Abs(SpecialFunctions.LogGamma(a)) + Math.Abs(SpecialFunctions.LogGamma(b)) + Math.Abs(SpecialFunctions.LogGamma(a + b))
         + Math.Abs(a * Math.Log(x)) + Math.Abs(b * Math.Log(1.0 - x));

    [Theory]
    [InlineData(1e-8, 1.8420680738180209e+1)]
    [InlineData(0.1, 2.252712651734206)]
    [InlineData(0.5, 5.7236494292470009e-1)]
    [InlineData(1.0, 0.0)]
    [InlineData(1.5, -1.2078223763524522e-1)]
    [InlineData(2.0, 0.0)]
    [InlineData(2.5, 2.8468287047291916e-1)]
    [InlineData(10.0, 1.280182748008147e+1)]
    [InlineData(100.0, 3.591342053695754e+2)]
    [InlineData(1e5, 1.0512877089736569e+6)]
    public void LogGamma_MatchesReference(double x, double expected)
        => Assert.True(Math.Abs(SpecialFunctions.LogGamma(x) - expected) <= 2e-14 * Math.Max(1.0, Math.Abs(expected)),
                       $"lnΓ({x}) = {SpecialFunctions.LogGamma(x):R}, expected {expected:R}");

    [Theory]
    [InlineData(0.5, 0.01, 1.1246291601828489e-1, 8.8753708398171511e-1, -2.1851317470723739, -1.193049737373956e-1)]
    [InlineData(0.5, 2.0, 9.5449973610364159e-1, 4.5500263896358414e-2, -4.6567912292390164e-2, -3.0900371531220866)]
    [InlineData(0.5, 800.0, 1.0, 0.0, 0.0, -8.0391529483319384e+2)]
    [InlineData(1.0, 1.0, 6.3212055882855768e-1, 3.6787944117144232e-1, -4.5867514538708189e-1, -1.0)]
    [InlineData(3.5, 2.0, 2.2022259152428408e-1, 7.7977740847571592e-1, -1.5131164644779843, -2.4874677377441855e-1)]
    [InlineData(3.5, 10.0, 9.9443031692705443e-1, 5.5696830729455713e-3, -5.5852515924360062e-3, -5.1904171255943225)]
    [InlineData(50.0, 40.0, 7.0335066659394954e-2, 9.2966493334060505e-1, -2.6544847900314759, -7.2931044484225142e-2)]
    [InlineData(50.0, 60.0, 9.1559331890630817e-1, 8.440668109369183e-2, -8.8182987887198306e-2, -2.4721087206353381)]
    [InlineData(1000.0, 950.0, 5.5054686230738034e-2, 9.4494531376926197e-1, -2.8994282926292987, -5.6628222190107071e-2)]
    [InlineData(1000.0, 1100.0, 9.9894067674607002e-1, 1.0593232539299773e-3, -1.0598847333686436e-3, -6.8501250144362653)]
    [InlineData(0.01, 0.5, 9.9437324380603282e-1, 5.626756193967184e-3, -5.6426460200850016e-3, -5.1802221673589263)]
    [InlineData(2.5, 1e-8, 3.0090110907617638e-21, 1.0, -4.7252675469370845e+1, -3.0090110907617638e-21)]
    public void RegularizedGamma_MatchesReference(double a, double x, double p, double q, double logP, double logQ)
    {
        double rel = PrefactorTolerance(GammaPrefactor(a, x));
        Close(p, SpecialFunctions.RegularizedGammaP(a, x), rel);
        Close(q, SpecialFunctions.RegularizedGammaQ(a, x), rel);
        CloseLog(logP, SpecialFunctions.LogRegularizedGammaP(a, x), rel);
        CloseLog(logQ, SpecialFunctions.LogRegularizedGammaQ(a, x), rel);
    }

    [Theory]
    [InlineData(0.3, 2.0, 5.0, -5.4502894510748891e-1, -8.6708398781950683e-1)]
    [InlineData(0.999, 0.5, 0.5, -2.0340514397558013e-2, -3.9052936169690375)]
    [InlineData(1e-10, 2.5, 0.5, -5.8644818689527257e+1, -3.3953054527483614e-26)]
    [InlineData(0.5, 100.0, 100.0, -6.9314718055994531e-1, -6.9314718055994531e-1)]
    [InlineData(0.2, 1000.0, 3000.0, -3.2450660887386052e+1, -8.0696935495279592e-15)]
    [InlineData(0.97, 3.0, 40.0, -9.8639759782203305e-59, -1.3356363115608764e+2)]
    public void RegularizedBeta_MatchesReference(double x, double a, double b, double logI, double logComplement)
    {
        double rel = PrefactorTolerance(BetaPrefactor(x, a, b));
        CloseLog(logI, SpecialFunctions.LogRegularizedBeta(x, a, b), rel);
        CloseLog(logComplement, SpecialFunctions.LogRegularizedBetaComplement(x, a, b), rel);
        Close(Math.Exp(logI), SpecialFunctions.RegularizedBeta(x, a, b), rel);
    }

    [Theory]
    [InlineData(-40.0, -8.0460844201375379e+2)]
    [InlineData(-38.0, -7.2655721601882013e+2)]
    [InlineData(-10.0, -5.3231285150512471e+1)]
    [InlineData(-5.0, -1.5064998393988726e+1)]
    [InlineData(-1.0, -1.8410216450092635)]
    [InlineData(-1e-10, -6.9314718063973377e-1)]
    [InlineData(0.0, -6.9314718055994531e-1)]
    [InlineData(0.5, -3.6894641528865639e-1)]
    [InlineData(3.0, -1.3508099647481938e-3)]
    [InlineData(8.3, -5.2055697448902853e-17)]
    public void Normal_LogCdf_MatchesReference(double z, double logPhi)
    {
        Normal n = Normal.Standard;
        CloseLog(logPhi, n.LogCdf(z), 1e-13);
        CloseLog(logPhi, n.LogSurvival(-z), 1e-13);
        if (logPhi > -700.0)   // below that the linear value is subnormal or zero, and only the log keeps digits
        {
            Close(Math.Exp(logPhi), n.Cdf(z));
            Close(Math.Exp(logPhi), Normal.MeanSd(2.0, 3.0).Cdf(2.0 + 3.0 * z), 1e-12);
        }
    }

    [Theory]
    [InlineData(-1.0, 2.0, -2.0016629432446258e-1)]
    [InlineData(8.0, 9.0, -3.5013618593437148e+1)]
    [InlineData(-1e-10, 1e-10, -2.3251642282585184e+1)]
    [InlineData(35.0, 36.0, -616.97510126192251)]
    [InlineData(-36.0, -35.0, -616.97510126192251)]
    [InlineData(3.0, double.PositiveInfinity, -6.6077262215103495)]
    [InlineData(double.NegativeInfinity, double.PositiveInfinity, 0.0)]
    [InlineData(double.NegativeInfinity, -37.0, -6.8903058557689059e+2)]
    public void Normal_IntervalProbability_MatchesReference(double a, double b, double logP)
    {
        Normal n = Normal.Standard;
        CloseLog(logP, n.LogIntervalProbability(a, b), 1e-12);
        Close(Math.Exp(logP), n.IntervalProbability(a, b), 1e-11);
    }

    [Theory]
    [InlineData(1.0, -1e6, -1.4960240443814008e+1, -3.1830993684428714e-7)]
    [InlineData(1.0, -50.0, -5.0568862015054332, -6.385694317785949e-3)]
    [InlineData(3.5, -2.0, -2.7635217709491866, -6.5145917605951117e-2)]
    [InlineData(3.5, -1e-8, -6.9314718799539221e-1, -6.9314717312449846e-1)]
    [InlineData(3.5, 0.0, -6.9314718055994531e-1, -6.9314718055994531e-1)]
    [InlineData(3.5, 1.0, -2.1154683736224546e-1, -1.6572182959995702)]
    [InlineData(30.0, 10.0, -2.2876257041409728e-11, -2.4500921553658308e+1)]
    [InlineData(1.0, 1e200, 0.0, -4.6166174848465854e+2)]
    [InlineData(2.5, -1e30, -1.7302330351910826e+2, 0.0)]
    public void StudentT_LogCdf_MatchesReference(double nu, double t, double logF, double logS)
    {
        StudentT d = StudentT.Standard(nu);
        CloseLog(logF, d.LogCdf(t), 1e-12);
        CloseLog(logS, d.LogSurvival(t), 1e-12);
        Close(Math.Exp(logF), d.Cdf(t), 1e-12);
        Close(Math.Exp(logF), StudentT.LocationScale(nu, -1.0, 0.25).Cdf(-1.0 + 0.25 * t), 1e-11);
    }

    [Theory]
    [InlineData(3.5, -1.0, 2.0, -2.9267838401464737e-1)]
    [InlineData(3.5, -1e-9, 1e-9, -2.1019592247645968e+1)]
    [InlineData(1.0, 100.0, 200.0, -6.4431055835572827)]
    [InlineData(5.0, 40.0, 41.0, -1.8356188744903019e+1)]
    [InlineData(5.0, -41.0, -40.0, -1.8356188744903019e+1)]
    public void StudentT_IntervalProbability_MatchesReference(double nu, double a, double b, double logP)
    {
        StudentT d = StudentT.Standard(nu);
        CloseLog(logP, d.LogIntervalProbability(a, b), 1e-12);
    }

    [Theory]
    [InlineData(2.5, 3.0, 0.4, 2.0852587940567533e-1, 7.9147412059432467e-1)]
    [InlineData(2.5, 3.0, 4.0, 9.9978288705654728e-1, 2.1711294345272332e-4)]
    [InlineData(0.3, 0.5, 1e-6, 1.4344013128284508e-2, 9.8565598687171549e-1)]
    public void Gamma_RateAndScale_AreTheSameDistribution(double shape, double rate, double x, double cdf, double survival)
    {
        Gamma byRate = Gamma.ShapeRate(shape, rate);
        Gamma byScale = Gamma.ShapeScale(shape, 1.0 / rate);
        Close(cdf, byRate.Cdf(x));
        Close(survival, byRate.Survival(x));
        Close(cdf, byScale.Cdf(x), 1e-12);
        Close(byRate.LogPdf(x), byScale.LogPdf(x), 1e-12);
        Close(shape / rate, byRate.Mean);
        Close(shape / (rate * rate), byRate.Variance);
    }

    [Theory]
    [InlineData(3.0, 2.0, 0.5, 2.3810330555354434e-1, 1.5888308335967186e-1)]
    [InlineData(0.5, 10.0, 100.0, 6.5472084601857703e-1, -6.4288276754098143)]
    [InlineData(7.0, 0.1, 0.003, 7.6990978844429113e-9, -9.5575362737895352)]
    public void InverseGamma_MatchesReference(double shape, double scale, double x, double cdf, double logPdf)
    {
        InverseGamma d = InverseGamma.ShapeScale(shape, scale);
        Close(cdf, d.Cdf(x), 1e-12);
        CloseLog(logPdf, d.LogPdf(x), 1e-13);
        Close(1.0 - cdf, d.Survival(x), 1e-12);
    }

    [Fact]
    public void InverseGamma_DrawsAreReciprocalGammaDraws_AndMatchTheMoments()
    {
        // IG(α, β) is 1/Gamma(α, rate β): same stream, reciprocal values.
        var a = new Xoshiro256PlusPlus(seed: 41);
        var b = new Xoshiro256PlusPlus(seed: 41);
        InverseGamma ig = InverseGamma.ShapeScale(4.5, 2.0);
        Gamma g = Gamma.ShapeRate(4.5, 2.0);
        for (int i = 0; i < 100; i++) Assert.Equal(1.0 / g.Sample(b), ig.Sample(a));

        var rng = new Xoshiro256PlusPlus(seed: 42);
        const int n = 200_000;
        double sum = 0.0, sumSq = 0.0;
        for (int i = 0; i < n; i++) { double v = ig.Sample(rng); sum += v; sumSq += v * v; }
        double mean = sum / n, variance = sumSq / n - mean * mean;
        Assert.Equal(2.0 / 3.5, ig.Mean, 15);
        Assert.Equal(4.0 / (3.5 * 3.5 * 2.5), ig.Variance, 15);
        Assert.True(Math.Abs(mean - ig.Mean) < 5.0 * Math.Sqrt(ig.Variance / n), $"mean {mean:R} vs {ig.Mean:R}");
        Assert.True(Math.Abs(variance - ig.Variance) < 0.05 * ig.Variance, $"variance {variance:R} vs {ig.Variance:R}");
    }

    [Fact]
    public void InverseGamma_MomentsAreInfiniteBelowTheirShapes()
    {
        Assert.Equal(double.PositiveInfinity, InverseGamma.ShapeScale(1.0, 2.0).Mean);
        Assert.Equal(double.PositiveInfinity, InverseGamma.ShapeScale(2.0, 2.0).Variance);
    }

    [Fact]
    public void StudentT_Draws_MatchTheCdf()
    {
        // Kolmogorov–Smirnov distance between 50,000 draws and the CDF; 1.95/√n is the 0.1% critical value.
        StudentT d = StudentT.LocationScale(3.5, 1.0, 2.0);
        var rng = new Xoshiro256PlusPlus(seed: 43);
        const int n = 50_000;
        var draws = new double[n];
        for (int i = 0; i < n; i++) draws[i] = d.Sample(rng);
        Array.Sort(draws);
        double ks = 0.0;
        for (int i = 0; i < n; i++)
        {
            double f = d.Cdf(draws[i]);
            ks = Math.Max(ks, Math.Max(f - (double)i / n, (i + 1.0) / n - f));
        }
        Assert.True(ks < 1.95 / Math.Sqrt(n), $"KS distance {ks:F5}");
    }

    [Fact]
    public void Truncated_Normal_IsTheConditionalLaw()
    {
        Normal n = Normal.MeanSd(0.5, 2.0);
        Truncated<Normal> t = Truncated.Of(n, -1.0, 3.0);

        Close(n.IntervalProbability(-1.0, 3.0), t.Mass, 1e-14);
        Assert.Equal(0.0, t.Cdf(-1.0));
        Assert.Equal(1.0, t.Cdf(3.0));
        Close(n.IntervalProbability(-1.0, 1.2) / t.Mass, t.Cdf(1.2), 1e-13);
        Close(n.IntervalProbability(1.2, 3.0) / t.Mass, t.Survival(1.2), 1e-13);
        Close(n.LogPdf(1.2) - Math.Log(t.Mass), t.LogPdf(1.2), 1e-14);
        Assert.Equal(double.NegativeInfinity, t.LogPdf(3.5));
        CloseLog(Math.Log(n.IntervalProbability(0.0, 3.0) / t.Mass), t.LogIntervalProbability(0.0, 10.0), 1e-13);
    }

    [Fact]
    public void Truncated_FarTail_KeepsItsDigits()
    {
        // [10, ∞) under N(0,1): mass 7.6e-24, so the density there is φ(x)/Q(10), not φ(x)/(1 − Φ(10)) = φ/0.
        Truncated<Normal> t = Truncated.Of(Normal.Standard, 10.0, double.PositiveInfinity);
        CloseLog(-5.3231285150512471e+1, t.LogMass, 1e-13);
        CloseLog(Normal.Standard.LogPdf(10.5) + 5.3231285150512471e+1, t.LogPdf(10.5), 1e-13);
        Assert.True(t.Cdf(10.5) > 0.99 && t.Cdf(10.5) < 1.0);
    }

    [Fact]
    public void Truncated_WithNoMass_IsRejected()
        => Assert.Throws<ArgumentException>(() => Truncated.Of(Gamma.ShapeRate(2.0, 1.0), -2.0, -1.0));

    [Fact]
    public void Factories_RejectInvalidParameters()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Gamma.ShapeRate(0.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Gamma.ShapeScale(1.0, -1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => InverseGamma.ShapeScale(1.0, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => StudentT.Standard(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => StudentT.LocationScale(3.0, 0.0, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Normal.MeanSd(0.0, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Beta.Shapes(1.0, 0.0));
        Assert.Throws<ArgumentException>(() => Truncated.Of(Normal.Standard, 2.0, 1.0));
    }
}
