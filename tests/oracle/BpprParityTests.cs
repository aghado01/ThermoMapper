using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// E36's ridge block against <c>BayesPPR</c> 0.2.0 on fixed draws (<c>r/oracles/bppr_oracle.R</c>; the BPPR brief, decision 8): two ridges
/// and a proposed third, K = 4 columns each, on four inputs. Read against the corrected ratios, as the brief settles:
/// <list type="bullet">
/// <item>the standardization, each ridge's knots (t₀ and the type-7 quantiles above it) and its basis columns, <c>get_mns_basis</c>'s;</item>
/// <item>the birth's collapsed likelihood ratio, −n/2·(ln sse′ − ln sse) − K/2·ln(g + 1), Zellner–Siow with the intercept inside the prior,
/// under the engine's zero centre, Jeffreys σ² and τ = 1/g held;</item>
/// <item>the birth's active-count and variable-set term with <c>adapt_act_feat</c>, against the engine's π(a)π(𝒥|a)/(q(a)q(𝒥|a)) at
/// ω₀ = υ₀ = 1: <c>BayesPPR</c> already takes 𝒥 uniformly at a = 1, the correction of the paper's eq. (mhb);</item>
/// <item>the count ratio: <c>BayesPPR</c>'s Poisson(10) against the engine's count law at h₁ = 10⁷, h₂ = 10⁶, whose λ is 10 to 10⁻⁶ (BASS's
/// count prior integrates λ; this is the closest it states).</item>
/// </list>
/// Opt-in: it runs where the R toolchain with the pinned renv is present.
/// </summary>
public sealed class BpprParityTests(ITestOutputHelper output)
{
    private const int K = 4;
    private const double G = 60.0;   // BayesPPR's var_coefs, the g-prior's g; the engine's τ = 1/g

    [Fact]
    public void TheRidgeBasisAndTheBirthsRatio_MeetBayesPPR_OnFixedDraws()
    {
        if (!ROracle.IsAvailable) return;
        const int n = 60, p = 4;
        var rng = new Xoshiro256PlusPlus(3661);
        double[][] x = [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray())];
        double[] y = [.. Enumerable.Range(0, n).Select(i => Math.Exp(0.6 * x[0][i] + 0.8 * x[1][i]) + 0.05 * Normal.SampleStandard(rng))];
        var data = new BassData([.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j}", x[j]))], new ScalarResponse(y));
        var config = new BassConfig
        {
            Structure = new BassStructureSpec
            {
                MaxInteraction = 3, MaxBasis = 10, OrderWeight = 1.0, VariableWeight = 1.0, Directional = new DirectionalBasis { Columns = K },
            },
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.Jeffreys,
            Coefficients = BassCoefficientPrior.G(TauPrior.Held(1.0 / G)),
            Count = new BassCountPrior(1e7, 1e6),
            Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 1, Samples = 1 },
        };
        BassModel model = BassModel.Create(config, data);
        DirectionalTermSpace space = model.Prior.Directional!;
        BassTerm Unit(int[] feat, double[] dir, double at)
        {
            var d = new ProjectionDirection(feat, dir);
            KnotRange r = space.KnotRangeOf(d)!.Value;
            return new BassTerm(directional: new DirectionalHinge(d, r.Lower + at * r.Width));
        }
        BassTerm[] current = [Unit([0, 1], [0.6, 0.8], 0.3), Unit([2], [-1.0], 0.6)];
        BassTerm born = Unit([1, 3], [-0.8, 0.6], 0.45);

        object Ridge(BassTerm t) => new { feat = t.Directional!.Direction.Variables, dir = t.Directional.Direction.Weights, knot0 = t.Directional.Knot };
        string fixture = Path.Combine(Path.GetTempPath(), $"bppr_{Guid.NewGuid():N}.json");
        File.WriteAllText(fixture, JsonSerializer.Serialize(new
        {
            x, y, K, var_coefs = G, n_act_max = 3, n_ridge_mean = 10,
            ridges = current.Select(Ridge).ToArray(), born = Ridge(born),
        }));
        JsonElement r = ROracle.Run("oracles/bppr_oracle.R", fixture);
        File.Delete(fixture);
        double[] Doubles(JsonElement e) => [.. e.EnumerateArray().Select(v => v.GetDouble())];

        // The standardization: z = a + b·u at the training cases against bppr()'s (x − mean)/sd.
        double worstZ = 0.0;
        double[][] xSt = [.. r.GetProperty("x_st").EnumerateArray().Select(Doubles)];
        for (int j = 0; j < p; j++)
        {
            (double a, double b) = model.Inputs.Standardizing(j)!.Value;
            ReadOnlySpan<double> u = model.Inputs.Cases.Continuous(j);
            for (int i = 0; i < n; i++) worstZ = Math.Max(worstZ, Math.Abs(a + b * u[i] - xSt[j][i]));
        }

        // Each ridge's knots and basis columns.
        double worstKnot = 0.0, worstColumn = 0.0;
        void Ridges(BassTerm unit, JsonElement theirs)
        {
            double[] knots = space.RidgeKnots(unit.Directional!, K)!, rKnots = Doubles(theirs.GetProperty("knots"));
            for (int k = 0; k < knots.Length; k++) worstKnot = Math.Max(worstKnot, Math.Abs(knots[k] - rKnots[k]) / Math.Max(1.0, Math.Abs(rKnots[k])));
            double[][] rBasis = [.. theirs.GetProperty("basis").EnumerateArray().Select(Doubles)];
            var column = new double[n];
            for (int l = 1; l <= K; l++)
            {
                model.Basis.CaseColumn(new BassTerm(ridge: new RidgeColumn(unit.Directional!.Direction, knots, l)), model.Inputs.Cases, column);
                for (int i = 0; i < n; i++) worstColumn = Math.Max(worstColumn, Math.Abs(column[i] - rBasis[l - 1][i]) / Math.Max(1.0, Math.Abs(rBasis[l - 1][i])));
            }
        }
        JsonElement[] rRidges = [.. r.GetProperty("ridges").EnumerateArray()];
        for (int k = 0; k < current.Length; k++) Ridges(current[k], rRidges[k]);
        Ridges(born, r.GetProperty("born"));

        // The birth's collapsed likelihood ratio.
        double LogMarginal(BassTerm[] units)
        {
            var chain = new BassChain(model, units, Math.Log(1.0 / G), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
            return model.Target.Evaluate(chain.Current).LogMarginal;
        }
        double engineLikelihood = LogMarginal([.. current, born]) - LogMarginal(current), rLikelihood = r.GetProperty("lik_ratio").GetDouble();

        // The birth's active-count and variable-set term, and the count ratio.
        var counts = new BassTermCounts(model.Prior);
        foreach (BassTerm t in current) counts.Add(t);
        double engineActFeat = model.Prior.LogPrior(born) - model.Kernel.LogBirth(born, counts), rActFeat = r.GetProperty("log_act_feat").GetDouble();
        double engineCount = model.Prior.LogCountRatio(current.Length), rCount = r.GetProperty("log_count").GetDouble();

        output.WriteLine($"standardization: worst |z − bppr's| {worstZ:E2}");
        output.WriteLine($"knots: worst {worstKnot:E2}; basis columns: worst {worstColumn:E2} (relative to max(1, |·|))");
        output.WriteLine($"likelihood ratio: engine {engineLikelihood:R}, BayesPPR {rLikelihood:R}, |Δ| {Math.Abs(engineLikelihood - rLikelihood):E2}");
        output.WriteLine($"active count and variables: engine {engineActFeat:R}, BayesPPR {rActFeat:R}, |Δ| {Math.Abs(engineActFeat - rActFeat):E2}");
        output.WriteLine($"count ratio: engine at h₁ = 1e7, h₂ = 1e6 {engineCount:R}, BayesPPR's Poisson(10) {rCount:R}, |Δ| {Math.Abs(engineCount - rCount):E2}");
        Assert.True(worstZ <= 1e-13, $"standardization {worstZ:E2}");
        Assert.True(worstKnot <= 1e-13, $"knots {worstKnot:E2}");
        Assert.True(worstColumn <= 1e-12, $"columns {worstColumn:E2}");
        Assert.True(Math.Abs(engineLikelihood - rLikelihood) <= 1e-9, $"likelihood ratio {Math.Abs(engineLikelihood - rLikelihood):E2}");
        Assert.True(Math.Abs(engineActFeat - rActFeat) <= 1e-12, $"active count and variables {Math.Abs(engineActFeat - rActFeat):E2}");
        Assert.True(Math.Abs(engineCount - rCount) <= 2e-6, $"count ratio {Math.Abs(engineCount - rCount):E2}");
    }
}
