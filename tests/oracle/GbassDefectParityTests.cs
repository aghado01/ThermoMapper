using System;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// E1's defects of GBASS 2.0.1, run against the installed package (<c>gbass_defects_oracle.R</c>; E1's build brief, slice 0). Each is
/// confirmed or struck here, with its reference computed on this side where it can be. These are corrections, not parity variants: the
/// engine takes the corrected reference.
/// <list type="bullet">
/// <item><b>Confirmed:</b> <c>nw_est_mom</c>'s kurtosis; <c>nw_gamma_prior</c>'s stray −1; the w step at β ≠ 0 dropping w's prior a and
/// its floor; the npart edge between birth and mutation; the GBP v step scaled from mean(ln v); γ's draw untruncated.</item>
/// <item><b>Struck for the code:</b> the paper's two slips, r = ŷ − y and Table 1's t row, which GBASS's code has right.</item>
/// <item><b>Found here:</b> <c>mu_gig</c> and <c>var_gig</c> return negative values where the moments are infinite, so <c>tbass</c> reports
/// a negative or undefined noise variance at ν ≤ 2; and <c>nw_gamma_prior</c>'s default targets, which no normal on γ meets, corrected
/// or not.</item>
/// </list>
/// </summary>
public sealed class GbassDefectParityTests(ITestOutputHelper output)
{
    private static readonly Lazy<JsonElement?> Fixture = new(() =>
    {
        // Sixty cases on [0, 1]², a smooth surface and t(2) noise, for hbass's GBP step.
        var rng = new Xoshiro256PlusPlus(20261105);
        var t2 = StudentT.Standard(2.0);
        var rows = Enumerable.Range(0, 60).Select(_ =>
        {
            double x1 = rng.NextDouble(), x2 = rng.NextDouble();
            return new[] { x1, x2, Math.Sin(2 * Math.PI * x1) + x2 + 0.2 * t2.Sample(rng) };
        }).ToArray();
        return GbassOracle.Run("gbass_defects_oracle.R", rows, "3000", "2000", "5");
    });

    [Fact]
    public void NwEstMom_Kurtosis_IsCorrected_AgainstTheMixture()
    {
        if (Fixture.Value is not { } r) return;
        foreach (JsonElement c in r.GetProperty("nw_kurtosis").EnumerateArray())
        {
            double alpha = c.GetProperty("alpha").GetDouble(), beta = c.GetProperty("beta").GetDouble(), delta = c.GetProperty("delta").GetDouble();
            double gamma = Math.Sqrt(alpha * alpha - beta * beta), m = delta / gamma, lambda = delta * delta;
            // X − EX = β(V − m) + √V·Z with V ~ IG(mean m, shape λ), by quadrature over ln V.
            double LogIg(double u) { double v = Math.Exp(u); return 0.5 * Math.Log(lambda / (2 * Math.PI)) - 1.5 * u - lambda * (v - m) * (v - m) / (2 * m * m * v) + u; }
            double E(Func<double, double> g) => Math.Exp(GbassOracle.LogIntegral(u => LogIg(u) + Math.Log(g(Math.Exp(u))), Math.Log(m)));
            double c2 = E(v => beta * beta * (v - m) * (v - m) + v);
            double c4 = E(v => Math.Pow(beta * (v - m), 4) + 6 * beta * beta * (v - m) * (v - m) * v + 3 * v * v);
            double mixture = c4 / (c2 * c2), corrected = 3 + 3 * (1 + 4 * beta * beta / (alpha * alpha)) / (delta * gamma);
            double code = c.GetProperty("code").GetDouble();
            output.WriteLine($"NIG({alpha}, {beta}, {delta}): mixture {mixture:F6}, corrected {corrected:F6}, nw_est_mom {code:F6}");
            Assert.Equal(1.0, c.GetProperty("mixture").GetDouble() / mixture, 9);
            Assert.Equal(1.0, corrected / mixture, 9);
            Assert.True(Math.Abs(code - mixture) > 0.01, "nw_est_mom's kurtosis is the mixture's");
        }
    }

    [Fact]
    public void NwGammaPrior_SecondProbability_HasAStrayMinusOne()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement p = r.GetProperty("nw_gamma_prior");
        // P(ξ < q₁) = P(|γ| > q₁⁻² − 1) and P(ξ > q₂) = P(|γ| < q₂⁻² − 1), γ ~ N(100, 100²).
        static double Inside(double t, double m, double s) => Normal.MeanSd(m, s).Cdf(t) - Normal.MeanSd(m, s).Cdf(-t);
        double below = 1 - Inside(Math.Pow(0.1, -2) - 1, 100, 100), above = Inside(Math.Pow(0.9, -2) - 1, 100, 100);
        Assert.Equal(below, p.GetProperty("definition").GetProperty("below_q1").GetDouble(), 12);
        Assert.Equal(above, p.GetProperty("definition").GetProperty("above_q2").GetDouble(), 12);
        Assert.Equal(below, p.GetProperty("code").GetProperty("lhs1").GetDouble(), 12);
        double lhs2 = p.GetProperty("code").GetProperty("lhs2").GetDouble();
        output.WriteLine($"P(xi > 0.9) at (100, 100): {above:F6}, as written {lhs2:F6}");
        Assert.True(Math.Abs(lhs2 - above) > 1e-3);
        Assert.True(p.GetProperty("corrected_changed").GetBoolean());
        // Its default targets (0.5, 0.05) are out of reach of any normal on γ: neither fit meets them.
        foreach (string fit in new[] { "as_written", "corrected" })
        {
            JsonElement f = p.GetProperty(fit);
            double m = f.GetProperty("m_gamma").GetDouble(), s = f.GetProperty("s_gamma").GetDouble();
            double b = 1 - Inside(99, m, s), a = Inside(Math.Pow(0.9, -2) - 1, m, s);
            output.WriteLine($"{fit}: (m, s) = ({m:F3}, {s:F3}), P(xi < 0.1) {b:F4}, P(xi > 0.9) {a:F5}, targets 0.5 and 0.05");
            Assert.Equal(b, f.GetProperty("below_q1").GetDouble(), 9);
            Assert.True(Math.Abs(b - 0.5) > 0.1 || Math.Abs(a - 0.05) > 0.01);
        }
    }

    [Fact]
    public void WStep_AtNonzeroSkew_DropsThePriorsAAndTheFloor()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement w = r.GetProperty("w_step");
        JsonElement[] arms = [.. w.GetProperty("arms").EnumerateArray()];
        output.WriteLine($"floor {w.GetProperty("floor").GetDouble():G6}; beta = 0 kept the old w {arms[0].GetProperty("kept_old").GetDouble():F3}, below the floor {arms[0].GetProperty("below_floor").GetDouble():F3}; " +
                         $"beta = 0.4 below the floor {arms[1].GetProperty("below_floor").GetDouble():F3}");
        Assert.Equal(0.0, arms[0].GetProperty("below_floor").GetDouble());
        Assert.InRange(arms[0].GetProperty("kept_old").GetDouble(), 0.4, 0.6);   // the floor sits at the conditional's median
        Assert.InRange(arms[1].GetProperty("below_floor").GetDouble(), 0.3, 1.0);
        JsonElement[] mhn = [.. w.GetProperty("mhn").EnumerateArray()];
        foreach (string k in new[] { "alpha", "beta", "gamma" })
            Assert.Equal(mhn[0].GetProperty(k).GetDouble(), mhn[1].GetProperty(k).GetDouble());   // a = 0 and a = 5 alike
    }

    [Fact]
    public void Npart_BirthAdmitsWhatMutationRefuses()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement n = r.GetProperty("npart");
        Assert.Equal(n.GetProperty("npart").GetDouble(), n.GetProperty("positive_rows").GetDouble());
        Assert.True(n.GetProperty("birth_admits").GetBoolean());
        Assert.False(n.GetProperty("mutation_admits").GetBoolean());
    }

    [Fact]
    public void GbpStep_IsScaledFromTheMeanOfLogV()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement g = r.GetProperty("gbp_step");
        double step = g.GetProperty("prop_sigma").GetDouble(), meanSd = g.GetProperty("sd_mean_log_v").GetDouble();
        double[] each = GbassOracle.Numbers(g.GetProperty("sd_log_v"));
        double median = each.Order().ElementAt(each.Length / 2), accept = g.GetProperty("acc_rate_v").GetDouble();
        output.WriteLine($"hbass N {g.GetProperty("N").GetInt32()}, {g.GetProperty("seconds").GetDouble():F1} s: step {step:F3}, 2.38·sd(mean ln v) {2.38 * meanSd:F3}, " +
                         $"median sd(ln v_i) {median:F3} (2.38× is {2.38 * median:F3}), acceptance {accept:F3}");
        Assert.InRange(step / (2.38 * meanSd), 0.5, 2.0);   // set from mean(ln v), which the adaptation reads
        Assert.True(step < 0.5 * 2.38 * median, "the step is on each observation's own scale");
        Assert.True(accept > 0.6);
    }

    [Fact]
    public void GammaDraw_PutsMassBelowZero()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement g = r.GetProperty("gamma_sign");
        double s = g.GetProperty("s_gamma").GetDouble(), m = g.GetProperty("m_gamma").GetDouble(), n = g.GetProperty("N").GetDouble(), sumV = g.GetProperty("sum_v").GetDouble();
        double mean = (s * s * n + m) / (s * s * sumV + 1), sd = s / Math.Sqrt(s * s * sumV + 1), below = Normal.MeanSd(mean, sd).Cdf(0.0);
        output.WriteLine($"gamma | v at N {n}, sum v {sumV}: N({mean:F4}, {sd:F4}), mass below 0 {below:F4}");
        Assert.Equal(mean, g.GetProperty("mean").GetDouble(), 12);
        Assert.Equal(below, g.GetProperty("mass_below_zero").GetDouble(), 12);
        Assert.True(below > 0.05);
    }

    [Fact]
    public void TextSlips_AreAbsentFromTheCode()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement t = r.GetProperty("text_slips");
        Assert.Equal("y - yhat", t.GetProperty("residual").GetString());
        JsonElement vp = t.GetProperty("tbass_v_prior");
        Assert.Equal(5.0, vp.GetProperty("b").GetDouble());
        Assert.Equal(-2.5, vp.GetProperty("p").GetDouble());
        double[] x = GbassOracle.Numbers(t.GetProperty("x")), code = GbassOracle.Numbers(t.GetProperty("code_marginal")), table = GbassOracle.Numbers(t.GetProperty("table1_marginal"));
        var t5 = StudentT.Standard(5.0);
        for (int i = 0; i < x.Length; i++)
        {
            Assert.Equal(1.0, code[i] / Math.Exp(t5.LogPdf(x[i])), 9);
            // Table 1's b = ν/2 is a t_ν scaled by 1/√2.
            Assert.Equal(1.0, table[i] / (Math.Sqrt(2) * Math.Exp(t5.LogPdf(Math.Sqrt(2) * x[i]))), 9);
        }
    }

    [Fact]
    public void GigMoments_AreNegativeWhereInfinite_SoTbassReportsNoVariance()
    {
        if (Fixture.Value is not { } r) return;
        foreach (JsonElement c in r.GetProperty("gig_moments").EnumerateArray())
        {
            double df = c.GetProperty("df").GetDouble(), mu = GbassOracle.Number(c.GetProperty("mu_gig")), s2 = GbassOracle.Number(c.GetProperty("s2"));
            double varGig = GbassOracle.Number(c.GetProperty("var_gig"));
            output.WriteLine($"tbass df {df}: mu_gig {mu}, var_gig {varGig}, s2 {s2}, bias {GbassOracle.Number(c.GetProperty("bias"))}");
            if (!c.GetProperty("mean_exists").GetBoolean()) Assert.False(s2 > 0, $"s2 at df {df}");
            else Assert.Equal(1.0, s2, 12);   // w = 1 is the variance under c = (ν − 2)/ν
            if (!c.GetProperty("var_exists").GetBoolean()) Assert.False(varGig > 0, $"var_gig at df {df}");
        }
    }
}
