using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// khaos's samplers against the exact structure law (<c>khaos_exactlaw_oracle.R</c>; E19's build brief, slice 0): on four inputs at
/// degree 3 and order 2, 𝒜 holds 30 multi-indices and the structures of at most two terms number 466, few enough to weigh each, with λ
/// integrated and, under the g-prior, g₀² integrated by quadrature. The corrected samplers (<c>khaos_corrected</c>), the readouts'
/// corrected reference, must reproduce that law within their own Monte Carlo error (batch means): P(M) and each term's inclusion, and
/// under the g-prior g₀²'s CDF at 300, 1000 and 10⁴. (Not E ln g₀²: with no term held, g₀² meets only the intercept under its vague
/// IG(0.001, 1000) prior, so ln g₀² has a long right tail and its mean's batch-means error is unreliable.) The samplers as written miss
/// the law by many standard errors, through the defects
/// <see cref="KhaosDefectParityTests"/> confirms; those two runs are long, so they wait for <c>BASS_QUALIFICATION</c> = 1 and otherwise
/// return at once.
/// </summary>
public sealed class KhaosExactLawTests(ITestOutputHelper output)
{
    private const int Iterations = 8_000;
    private static bool Qualification => Environment.GetEnvironmentVariable("BASS_QUALIFICATION") is "1";

    /// <summary>Forty cases on [0, 1]⁴, y = 0.5ψ₁(x₁) + 0.35ψ₁(x₂)ψ₁(x₃) + N(0, 1): weak enough that the posterior spreads over
    /// structures.</summary>
    internal static double[][] Data()
    {
        var rng = new Xoshiro256PlusPlus(20261109);
        return [.. Enumerable.Range(0, 40).Select(_ =>
        {
            double[] x = [.. Enumerable.Range(0, 4).Select(_ => rng.NextDouble())];
            double y = 0.5 * KhaosOracle.Psi(1, x[0]) + 0.35 * KhaosOracle.Psi(1, x[1]) * KhaosOracle.Psi(1, x[2]) + Normal.SampleStandard(rng);
            return x.Append(y).ToArray();
        })];
    }

    private static readonly Dictionary<string, Lazy<JsonElement?>> Runs = new()
    {
        ["ridge corrected"] = new(() => GbassOracle.Run("khaos_exactlaw_oracle.R", Data(), "ridge", "corrected", Iterations.ToString(), "3")),
        ["gprior corrected"] = new(() => GbassOracle.Run("khaos_exactlaw_oracle.R", Data(), "gprior", "corrected", Iterations.ToString(), "3")),
        ["ridge written"] = new(() => GbassOracle.Run("khaos_exactlaw_oracle.R", Data(), "ridge", "written", Iterations.ToString(), "3")),
        ["gprior written"] = new(() => GbassOracle.Run("khaos_exactlaw_oracle.R", Data(), "gprior", "written", Iterations.ToString(), "3")),
    };

    /// <summary>The worst z over P(M = m), each term's inclusion and, under the g-prior, g₀²'s CDF at three points: the chain's reading
    /// less the exact, over its Monte Carlo standard error by batch means, floored at the independent draws' binomial error.</summary>
    private (double WorstZ, double WorstM) Compare(string run, JsonElement r)
    {
        JsonElement exact = r.GetProperty("exact"), chain = r.GetProperty("chain");
        Assert.Equal(466, exact.GetProperty("structures").GetInt32());
        int draws = chain.GetProperty("draws").GetInt32();
        double[] me = GbassOracle.Numbers(exact.GetProperty("M")), mc = GbassOracle.Numbers(chain.GetProperty("M")), ms = GbassOracle.Numbers(chain.GetProperty("M_se"));
        double[] ie = GbassOracle.Numbers(exact.GetProperty("inclusion")), ic = GbassOracle.Numbers(chain.GetProperty("inclusion"));
        double[] isd = GbassOracle.Numbers(chain.GetProperty("inclusion_se"));
        Assert.Equal(30, ie.Length);
        Assert.Equal(1.0, me.Sum(), 12);
        double Z(double chainValue, double exactValue, double se, bool probability) =>
            Math.Abs(chainValue - exactValue) / (probability ? Math.Max(se, Math.Sqrt(exactValue * (1 - exactValue) / draws)) : se);
        var z = new List<double>();
        z.AddRange(me.Select((e, m) => Z(mc[m], e, ms[m], true)));
        z.AddRange(ie.Select((e, t) => Z(ic[t], e, isd[t], true)));
        bool gprior = r.GetProperty("prior").GetString() == "gprior";
        string cdf = "";
        if (gprior)
        {
            // The exact law's g₀² integrand, read through the pencil, against the dense multivariate t.
            Assert.True(exact.GetProperty("pencil_check").GetDouble() < 1e-9);
            double[] ce = GbassOracle.Numbers(exact.GetProperty("g0sq_cdf")), cc = GbassOracle.Numbers(chain.GetProperty("g0sq_cdf"));
            double[] cs = GbassOracle.Numbers(chain.GetProperty("g0sq_cdf_se")), cuts = GbassOracle.Numbers(exact.GetProperty("g0sq_cuts"));
            z.AddRange(ce.Select((e, k) => Z(cc[k], e, cs[k], true)));
            cdf = $"; P(g₀² < {string.Join(", ", cuts)}) {string.Join(" / ", cc.Select(v => v.ToString("F3")))} against {string.Join(" / ", ce.Select(v => v.ToString("F3")))}";
        }
        double worstM = me.Zip(mc, (a, b) => Math.Abs(a - b)).Max(), worstI = ie.Zip(ic, (a, b) => Math.Abs(a - b)).Max();
        output.WriteLine($"{run}, {draws} kept draws ({chain.GetProperty("seconds").GetDouble():F1} s): P(M = 0, 1, 2) " +
                         $"{string.Join(" / ", mc.Select(v => v.ToString("F3")))} against the exact {string.Join(" / ", me.Select(v => v.ToString("F3")))}; " +
                         $"inclusion off by up to {worstI:F3}" + cdf +
                         $"; worst z {z.Max():F2} over {z.Count} readings; draws holding a term outside 𝒜 {chain.GetProperty("outside").GetDouble():P1}; " +
                         $"ESS of M {GbassOracle.Number(chain.GetProperty("ess_M")):F0}" +
                         (gprior ? $", of ln g₀² {GbassOracle.Number(chain.GetProperty("ess_log_g0sq")):F0}, g₀² moving on {chain.GetProperty("g0sq_moves").GetDouble():P1} of iterations" : ""));
        return (z.Max(), worstM);
    }

    [Theory]
    [InlineData("ridge corrected")]
    [InlineData("gprior corrected")]
    public void CorrectedSampler_ReproducesTheExactStructureLaw(string run)
    {
        if (Runs[run].Value is not { } r) return;
        (double worstZ, _) = Compare(run, r);
        Assert.Equal(0.0, r.GetProperty("chain").GetProperty("outside").GetDouble());
        Assert.True(worstZ < 4.5, $"a reading {worstZ:F2} standard errors from the exact law");
    }

    [Theory]
    [InlineData("ridge written")]
    [InlineData("gprior written")]
    public void WrittenSampler_MissesTheExactStructureLaw(string run)
    {
        if (!Qualification || Runs[run].Value is not { } r) return;
        (double worstZ, double worstM) = Compare(run, r);
        Assert.True(worstZ > 6 && worstM > 0.05, $"P(M) off by only {worstM:F3}, worst z {worstZ:F2}");
    }
}
