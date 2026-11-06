using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// GBASS's primitives against their laws (<c>gbass_primitives_oracle.R</c>; E1's build brief, slice 0), each law integrated on both
/// sides:
/// <list type="bullet">
/// <item><b>The GIG.</b> <see cref="GigParameters"/>' kernel integrated here gives R's moments by quadrature and by Bessel K, so the two
/// parameterizations, (p, a, b) and GIGrvg's (λ, ψ, χ), are one. <c>GIGrvg::rgig</c> and GBASS's <c>rgig2</c> carry the law's
/// probabilities at their deciles, the limits a = 0 and b = 0 included. <c>mu_gig</c> and <c>var_gig</c> are the moments where those
/// exist, and where they do not (a = 0, p ≥ −1 for the mean, p ≥ −2 for the variance) they return finite or negative values: a
/// defect, recorded beside the law.</item>
/// <item><b>The modified half-normal.</b> <c>rmhn</c> carries the law of x^(α−1)exp(−βx² + γx) on each of its four branches, the series
/// branch at 0 &lt; α ≤ 1 included, which E1's scratch left unchecked.</item>
/// <item><b>Bessel K.</b> Base R's <c>besselK(x, ν, expon.scaled = TRUE)</c> against K's integral representation over far ranges: the
/// reference a log Bessel K of real order meets in slice 3.</item>
/// </list>
/// </summary>
public sealed class GbassPrimitiveParityTests(ITestOutputHelper output)
{
    private const int Draws = 50_000;   // GBASS's rgig2 draws one value a call, so the GIG cases set the oracle's time

    // kind, then three parameters: 1 GIG(p, a, b), 2 MHN(α, β, γ), 3 K_ν(x) as (x, ν, 0).
    private static readonly double[][] Cases =
    [
        [1, -2.5, 0, 5], [1, -0.75, 0, 1.5], [1, -1.5, 0, 3], [1, 1, 2, 0], [1, 0.5, 2.3, 0.04], [1, 0.5, 2.3, 0],
        [1, -8, 0.5, 30], [1, 2.5, 0.001, 1000], [1, -1, 0.64, 1.3], [1, 20, 50, 0.01], [1, -20, 0.01, 50],
        [2, 520, 300, -15], [2, 5.5, 2, 3], [2, 40, 1, 30], [2, 0.6, 1.3, 2.1], [2, 1, 0.7, 4], [2, 0.3, 2, 0.5], [2, 0.8, 0.5, -2],
        [2, 305, 12.5, 40], [2, 3, 1e-3, 0.5],
        [3, 1.5, 0.5, 0], [3, 1e-6, 3, 0], [3, 800, 40, 0], [3, 20, 250, 0], [3, 0.5, 0, 0], [3, 3, 1, 0], [3, 1e-3, 1, 0],
        [3, 50, 0.5, 0], [3, 2, 80, 0],
    ];

    private static readonly Lazy<JsonElement?> Fixture = new(() => GbassOracle.Run("gbass_primitives_oracle.R", Cases, Draws.ToString(), "20261105"));

    [Fact]
    public void GigLaw_IsOneParameterization_AndTheDrawsCarryIt()
    {
        if (Fixture.Value is not { } r) return;
        int checkedCases = 0;
        foreach (JsonElement c in r.GetProperty("gig").EnumerateArray())
        {
            var g = new GigParameters(c.GetProperty("p").GetDouble(), c.GetProperty("a").GetDouble(), c.GetProperty("b").GetDouble());
            JsonElement m = c.GetProperty("moments");
            // The law in u = ln x, kernel times the Jacobian x.
            double H(double u) => g.LogKernel(Math.Exp(u)) + u;
            double guess = g.A > 0 ? Math.Log(Math.Max((g.P + Math.Sqrt(g.P * g.P + g.A * g.B)) / g.A, 1e-300)) : Math.Log(g.B / (-2.0 * g.P));
            double logZ = GbassOracle.LogIntegral(H, guess);
            Assert.Equal(m.GetProperty("log_mass").GetDouble(), logZ, 9);
            bool meanExists = m.GetProperty("mean_exists").GetBoolean(), varExists = m.GetProperty("var_exists").GetBoolean();
            Assert.Equal(g.A > 0 || g.P < -1, meanExists);
            Assert.Equal(g.A > 0 || g.P < -2, varExists);
            double muGig = GbassOracle.Number(c.GetProperty("mu_gig")), varGig = GbassOracle.Number(c.GetProperty("var_gig"));
            if (meanExists)
            {
                double mean = Math.Exp(GbassOracle.LogIntegral(u => H(u) + u, guess) - logZ);
                Assert.Equal(1.0, m.GetProperty("mean").GetDouble() / mean, 9);
                Assert.Equal(1.0, muGig / mean, 9);
                if (g.A > 0 && g.B > 0) Assert.Equal(1.0, m.GetProperty("mean_bessel").GetDouble() / mean, 9);
                if (varExists)
                {
                    double second = Math.Exp(GbassOracle.LogIntegral(u => H(u) + 2 * u, guess) - logZ);
                    double variance = second - mean * mean;
                    Assert.Equal(1.0, m.GetProperty("var").GetDouble() / variance, 7);
                    Assert.Equal(1.0, varGig / variance, 7);
                }
            }
            // Where a moment is infinite GBASS's formula still returns a number, finite or negative: the defect.
            if (!meanExists) Assert.True(muGig <= 0, $"mu_gig {muGig} at {g}");
            if (!varExists) Assert.True(varGig <= 0, $"var_gig {varGig} at {g}");
            foreach (string sampler in new[] { "rgig", "rgig2" })
            {
                if (c.GetProperty(sampler).ValueKind == JsonValueKind.Null) continue;
                JsonElement s = c.GetProperty(sampler);
                double[] deciles = GbassOracle.Numbers(s.GetProperty("deciles"));
                double[] cdf = [.. deciles.Select(q => Math.Exp(GbassOracle.LogIntegral(H, guess, upper: Math.Log(q)) - logZ))];
                Assert.Equal(GbassOracle.Numbers(s.GetProperty("cdf")), cdf, (x, y) => Math.Abs(x - y) < 1e-8);
                double worst = GbassOracle.DecileScores(cdf, s.GetProperty("n").GetInt32()).Max(Math.Abs);
                output.WriteLine($"GIG({g.P}, {g.A}, {g.B}) {sampler}: worst decile z {worst:F2}; mean {(meanExists ? muGig : double.NaN):G6}, mu_gig {muGig:G6}, var_gig {varGig:G6}");
                Assert.True(worst < 5.0, $"{sampler} at GIG({g.P}, {g.A}, {g.B}): decile z {worst:F2}");
            }
            checkedCases++;
        }
        Assert.Equal(Cases.Count(c => c[0] == 1), checkedCases);
    }

    [Fact]
    public void ModifiedHalfNormalDraws_CarryTheLaw_OnEveryBranch()
    {
        if (Fixture.Value is not { } r) return;
        var branches = new HashSet<string>();
        foreach (JsonElement c in r.GetProperty("mhn").EnumerateArray())
        {
            double alpha = c.GetProperty("alpha").GetDouble(), beta = c.GetProperty("beta").GetDouble(), gamma = c.GetProperty("gamma").GetDouble();
            double H(double u) { double x = Math.Exp(u); return alpha * u - beta * x * x + gamma * x; }
            double guess = Math.Log(Math.Max((gamma + Math.Sqrt(gamma * gamma + 8 * beta * Math.Max(alpha, 1e-3))) / (4 * beta), 1e-300));
            double logZ = GbassOracle.LogIntegral(H, guess);
            Assert.Equal(c.GetProperty("log_mass").GetDouble(), logZ, 9);
            double mean = Math.Exp(GbassOracle.LogIntegral(u => H(u) + u, guess) - logZ);
            double variance = Math.Exp(GbassOracle.LogIntegral(u => H(u) + 2 * u, guess) - logZ) - mean * mean;
            Assert.Equal(1.0, c.GetProperty("mean").GetDouble() / mean, 9);
            Assert.Equal(1.0, c.GetProperty("var").GetDouble() / variance, 6);
            JsonElement s = c.GetProperty("rmhn");
            double[] deciles = GbassOracle.Numbers(s.GetProperty("deciles"));
            double[] cdf = [.. deciles.Select(q => Math.Exp(GbassOracle.LogIntegral(H, guess, upper: Math.Log(q)) - logZ))];
            double worst = GbassOracle.DecileScores(cdf, s.GetProperty("n").GetInt32()).Max(Math.Abs);
            string branch = c.GetProperty("branch").GetString()!;
            branches.Add(branch);
            output.WriteLine($"MHN({alpha}, {beta}, {gamma}) {branch}: mean {mean:G8}, var {variance:G8}, worst decile z {worst:F2}");
            Assert.True(worst < 5.0, $"rmhn at ({alpha}, {beta}, {gamma}), {branch}: decile z {worst:F2}");
        }
        Assert.Equal(["negative", "normal", "series", "sqrt_gamma"], branches.Order());
    }

    [Fact]
    public void BesselK_MatchesItsIntegral_OverFarRanges()
    {
        if (Fixture.Value is not { } r) return;
        int n = 0;
        foreach (JsonElement c in r.GetProperty("bessel").EnumerateArray())
        {
            double x = c.GetProperty("x").GetDouble(), nu = c.GetProperty("nu").GetDouble();
            // ln K_ν(x) = ln ∫₀^∞ exp(−x cosh t) cosh(νt) dt, ln cosh written to hold its digits.
            double H(double t) => -x * Math.Cosh(t) + (nu * t > 20 ? nu * t - Math.Log(2) : Math.Log(Math.Cosh(nu * t)));
            double guess = nu * nu <= x ? 0.0 : Math.Asinh(nu / x);
            double ours = GbassOracle.LogIntegral(H, guess, lower: 0.0);
            double rs = c.GetProperty("log_k_r").GetDouble(), integral = c.GetProperty("log_k_integral").GetDouble();
            output.WriteLine($"ln K_{nu}({x}): besselK {rs:R}, R's integral {integral:R}, ours {ours:R}");
            Assert.False(c.GetProperty("warned").GetBoolean());
            Assert.True(Math.Abs(ours - integral) <= 1e-11 * Math.Max(1.0, Math.Abs(ours)), $"integral at ({x}, {nu})");
            Assert.True(Math.Abs(ours - rs) <= 1e-11 * Math.Max(1.0, Math.Abs(ours)), $"besselK at ({x}, {nu})");
            n++;
        }
        Assert.Equal(Cases.Count(c => c[0] == 3), n);
    }
}
