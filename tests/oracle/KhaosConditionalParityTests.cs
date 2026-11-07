using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// KHAOS's ratios and conditionals at frozen states (<c>khaos_conditionals_oracle.R</c>; E19's build brief, slice 0): the corrected
/// references the engine meets from slice 1 on, computed here, beside khaos 2.1.1's expressions as written.
/// <list type="bullet">
/// <item><b>The collapsed marginal,</b> ln p(y | structure) with β and σ² integrated, under the modified g-prior and the raw ridge: the
/// engine's <see cref="Maths.Regression.GPrior.PriorPencil"/> with prior precision D_g⁻¹ΨᵀΨD_g⁻¹ at τ = 1/g₀² (or I at 1/τ²) equals
/// the dense multivariate t over n × n. The ridge's structure ratio as written is that marginal's difference.</item>
/// <item><b>g₀²'s conditional</b> with β and σ² integrated, through the pencil at O(p) per g₀², and given β and σ², the conjugate
/// IG(a_g + (M + 1)/2, b_g + βᵀD_g⁻¹ΨᵀΨD_g⁻¹β/(2σ²)): each is the joint less a constant.</item>
/// <item><b>σ²'s conditional</b> given β under each prior: IG(a + n/2 + (M + 1)/2, b + RSS/2 + βᵀS₀⁻¹β/2), the joint less a constant.</item>
/// <item><b>The proposals.</b> The birth's forward probability of an inclusion vector, Σ_{q₀} J(q₀)·P(χ | q₀)/P(valid | q₀), is the
/// law of khaos's own draw; the degree mutation's ratio under uniform compositions and under random_partition's law; the variable
/// swap's ratio is its own kernel's.</item>
/// </list>
/// Where khaos's expressions as written differ from these, <see cref="KhaosDefectParityTests"/> holds the verdict.
/// </summary>
public sealed class KhaosConditionalParityTests(ITestOutputHelper output)
{
    [Fact]
    public void CollapsedMarginal_ThroughThePencil_IsTheDenseReference()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        (double[][] basis, double[] y) = KhaosOracle.State();
        double[] g = KhaosOracle.StateG();
        Assert.Equal(g, GbassOracle.Numbers(r.GetProperty("g")), (u, v) => Math.Abs(u - v) < 1e-15);
        int checkedStructures = 0;
        foreach (JsonElement s in r.GetProperty("marginal").GetProperty("structures").EnumerateArray())
        {
            int k = s.GetProperty("k").GetInt32();
            double[][] cols = KhaosOracle.Columns(basis, k + 1);
            double gprior = KhaosOracle.PencilMarginal(cols, y, KhaosOracle.GPriorPrecision(cols, g[..(k + 1)]), 1.0 / KhaosOracle.G0Sq, KhaosOracle.ASigma, KhaosOracle.BSigma);
            double ridge = KhaosOracle.PencilMarginal(cols, y, KhaosOracle.Identity(k + 1), 1.0 / KhaosOracle.Tau2, KhaosOracle.ASigma, KhaosOracle.BSigma);
            JsonElement gr = s.GetProperty("gprior"), rr = s.GetProperty("ridge");
            output.WriteLine($"structure {k}: g-prior ln p {gprior:F10} (dense {gr.GetProperty("dense").GetDouble():F10}), ridge {ridge:F10} (dense {rr.GetProperty("dense").GetDouble():F10})");
            Assert.Equal(gr.GetProperty("dense").GetDouble(), gprior, 9);
            Assert.Equal(rr.GetProperty("dense").GetDouble(), ridge, 9);
            Assert.Equal(gr.GetProperty("dense").GetDouble(), gr.GetProperty("closed").GetDouble(), 9);
            Assert.Equal(rr.GetProperty("dense").GetDouble(), rr.GetProperty("closed").GetDouble(), 9);
            checkedStructures++;
        }
        Assert.Equal(KhaosOracle.StateTerms.Length + 1, checkedStructures);
    }

    [Fact]
    public void RidgeStructureRatio_AsWritten_IsTheMarginalsDifference()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        JsonElement m = r.GetProperty("marginal");
        double[] dense = [.. m.GetProperty("structures").EnumerateArray().Select(s => s.GetProperty("ridge").GetProperty("dense").GetDouble())];
        double[] code = [.. m.GetProperty("births").EnumerateArray().Select(b => b.GetProperty("ridge_code").GetDouble())];
        for (int k = 1; k < dense.Length; k++)
        {
            output.WriteLine($"ridge birth {k - 1} → {k}: as written {code[k - 1]:F10}, the marginals' difference {dense[k] - dense[k - 1]:F10}");
            Assert.Equal(dense[k] - dense[k - 1], code[k - 1], 9);
        }
        Assert.Equal(-code[^1], m.GetProperty("death").GetProperty("ridge_code").GetDouble(), 10);
        double[] gcode = [.. m.GetProperty("births").EnumerateArray().Select(b => b.GetProperty("gprior_code").GetDouble())];
        Assert.Equal(-gcode[^1], m.GetProperty("death").GetProperty("gprior_code").GetDouble(), 10);
    }

    [Fact]
    public void G0SqConditional_DerivedForms_AreTheJointLessAConstant()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        (double[][] basis, double[] y) = KhaosOracle.State();
        double[] g = KhaosOracle.StateG();
        JsonElement s = r.GetProperty("g0sq");
        JsonElement[] grid = [.. s.GetProperty("grid").EnumerateArray()];
        double[] theta = [.. grid.Select(e => e.GetProperty("theta").GetDouble())];
        Assert.Equal(KhaosOracle.Theta, theta);
        // With β and σ² integrated: the pencil's marginal at τ = 1/g₀², one decomposition's worth per structure, plus ln IG(g₀²).
        double[,] precision = KhaosOracle.GPriorPrecision(basis, g);
        double[] integrated = [.. theta.Select(t => KhaosOracle.PencilMarginal(basis, y, precision, 1.0 / t, KhaosOracle.ASigma, KhaosOracle.BSigma)
                                                    + KhaosOracle.LogInverseGammaKernel(t, KhaosOracle.AG, KhaosOracle.BG))];
        double spreadIntegrated = KhaosOracle.Spread(grid.Select((e, i) => e.GetProperty("joint_integrated").GetDouble() - integrated[i]));
        // Given β and σ²: the conjugate inverse gamma.
        int columns = g.Length;
        double energy = 0;
        foreach (double[] row in basis)
        {
            double u = 0;
            for (int j = 0; j < columns; j++) u += row[j] * KhaosOracle.StateBeta[j] / g[j];
            energy += u * u;
        }
        double shape = KhaosOracle.AG + columns / 2.0, rate = KhaosOracle.BG + energy / (2 * KhaosOracle.S2);
        Assert.Equal(shape, s.GetProperty("given_shape").GetDouble(), 12);
        Assert.Equal(1.0, s.GetProperty("given_rate").GetDouble() / rate, 12);
        double spreadGiven = KhaosOracle.Spread(grid.Select((e, i) => e.GetProperty("joint_given").GetDouble() - KhaosOracle.LogInverseGammaKernel(theta[i], shape, rate)));
        output.WriteLine($"g₀²'s conditional over {theta.Length} points from {theta[0]} to {theta[^1]}: integrated through the pencil, the joint less it spreads {spreadIntegrated:G3} nats; " +
                         $"given β and σ², IG({shape}, {rate:G6}), {spreadGiven:G3}");
        Assert.True(spreadIntegrated < 1e-9 && spreadGiven < 1e-9);
    }

    [Fact]
    public void Sigma2Conditional_GivenBeta_IsTheJointLessAConstant()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        (double[][] basis, double[] y) = KhaosOracle.State();
        double[] g = KhaosOracle.StateG(), beta = KhaosOracle.StateBeta;
        int n = y.Length, columns = beta.Length;
        double rss = 0, energyG = 0;
        for (int i = 0; i < n; i++)
        {
            double fit = 0, u = 0;
            for (int j = 0; j < columns; j++) { fit += basis[i][j] * beta[j]; u += basis[i][j] * beta[j] / g[j]; }
            rss += (y[i] - fit) * (y[i] - fit);
            energyG += u * u;
        }
        JsonElement s = r.GetProperty("sigma2");
        Assert.Equal(1.0, s.GetProperty("rss").GetDouble() / rss, 12);
        double[] grid = GbassOracle.Numbers(s.GetProperty("grid"));
        var stated = new Dictionary<string, double>
        {
            ["gprior"] = KhaosOracle.BSigma + rss / 2 + energyG / (2 * KhaosOracle.G0Sq),
            ["ridge"] = KhaosOracle.BSigma + rss / 2 + beta.Sum(b => b * b) / (2 * KhaosOracle.Tau2),
        };
        foreach ((string prior, double rate) in stated)
        {
            JsonElement p = s.GetProperty(prior);
            double shape = KhaosOracle.ASigma + n / 2.0 + columns / 2.0;
            Assert.Equal(shape, p.GetProperty("stated_shape").GetDouble(), 12);
            Assert.Equal(1.0, p.GetProperty("stated_rate").GetDouble() / rate, 12);
            double[] joint = GbassOracle.Numbers(p.GetProperty("joint"));
            double spread = KhaosOracle.Spread(grid.Select((x, i) => joint[i] - KhaosOracle.LogInverseGammaKernel(x, shape, rate)));
            output.WriteLine($"σ² given β under the {prior}: IG({shape}, {rate:G6}); the joint less it spreads {spread:G3} nats");
            Assert.True(spread < 1e-9);
        }
    }

    [Fact]
    public void BirthForwardProbability_IsTheLawOfKhaosOwnDraw()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        int checkedCases = 0;
        foreach (JsonElement c in r.GetProperty("birth").EnumerateArray())
        {
            int order = c.GetProperty("order").GetInt32(), draws = c.GetProperty("draws").GetInt32();
            double[] J = GbassOracle.Numbers(c.GetProperty("J"));
            double[][] W = [.. c.GetProperty("weights").EnumerateArray().Select(GbassOracle.Numbers)];
            int p = W[0].Length;
            double[] valid = [.. W.Select(w => KhaosOracle.AllInclusions(p).Where(x => x.Sum() >= 1 && x.Sum() <= order).Sum(x => KhaosOracle.CoinProbability(x, w)))];
            Assert.Equal(valid, GbassOracle.Numbers(c.GetProperty("valid_exact")), (u, v) => Math.Abs(u - v) < 1e-13);
            double total = 0, worstZ = 0;
            foreach (JsonElement e in c.GetProperty("chis").EnumerateArray())
            {
                int[] chi = [.. GbassOracle.Numbers(e.GetProperty("chi")).Select(v => (int)v)];
                double exact = Enumerable.Range(0, order).Sum(q0 => J[q0] * KhaosOracle.CoinProbability(chi, W[q0]) / valid[q0]);
                Assert.Equal(exact, e.GetProperty("exact").GetDouble(), 13);
                double freq = e.GetProperty("freq").GetDouble();
                worstZ = Math.Max(worstZ, Math.Abs(freq - exact) / Math.Sqrt(exact * (1 - exact) / draws));
                total += exact;
            }
            output.WriteLine($"birth at η = {string.Join(",", GbassOracle.Numbers(c.GetProperty("eta")))}, order {order}: P(valid | q₀) = {string.Join(", ", valid.Select(v => v.ToString("F4")))}; " +
                             $"khaos's own draws against the exact law, worst z {worstZ:F2} over {draws}");
            Assert.Equal(1.0, total, 12);
            Assert.True(worstZ < 5.0);
            checkedCases++;
        }
        Assert.Equal(KhaosOracle.ProposalCases.Count(c => c[0] == 1), checkedCases);
    }

    [Fact]
    public void MutationRatios_UnderUniformCompositions_AndTheSwapsOwnKernel()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        int degreeCases = 0, swaps = 0;
        foreach (JsonElement c in r.GetProperty("mutate_degree").EnumerateArray())
        {
            int q = c.GetProperty("q").GetInt32(), dc = c.GetProperty("d_curr").GetInt32(), dn = c.GetProperty("d_cand").GetInt32();
            double[] dp = GbassOracle.Numbers(c.GetProperty("d_probs"));
            double At(int d) => dp[d - q];
            double uniform = Math.Log(At(dc)) - Math.Log(KhaosOracle.Choose(dc - 1, q - 1)) - Math.Log(At(dn)) + Math.Log(KhaosOracle.Choose(dn - 1, q - 1));
            Assert.Equal(uniform, c.GetProperty("uniform").GetDouble(), 12);
            double[] actual = [.. from a in KhaosOracle.Compositions(dc, q) from b in KhaosOracle.Compositions(dn, q)
                                  select Math.Log(At(dc) * KhaosOracle.SortedCutsProbability(a)) - Math.Log(At(dn) * KhaosOracle.SortedCutsProbability(b))];
            Assert.Equal(actual.Min(), c.GetProperty("actual_least").GetDouble(), 12);
            Assert.Equal(actual.Max(), c.GetProperty("actual_greatest").GetDouble(), 12);
            output.WriteLine($"degree mutation, order {q}, {dc} → {dn}, s_d {c.GetProperty("s_d").GetDouble()}: uniform compositions {uniform:F4}; " +
                             $"random_partition's own law {actual.Min():F4} to {actual.Max():F4}; as written {c.GetProperty("code").GetDouble():F4}");
            degreeCases++;
        }
        foreach (JsonElement c in r.GetProperty("mutate_vars").EnumerateArray())
            foreach (JsonElement s in c.GetProperty("swaps").EnumerateArray())
            {
                Assert.Equal(s.GetProperty("exact").GetDouble(), s.GetProperty("code").GetDouble(), 12);
                swaps++;
            }
        Assert.Equal(KhaosOracle.ProposalCases.Count(c => c[0] == 2), degreeCases);
        Assert.True(swaps > 0);
    }
}
