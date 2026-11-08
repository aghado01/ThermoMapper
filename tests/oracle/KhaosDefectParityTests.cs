using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// The KHAOS brief's defects of khaos 2.1.1, each confirmed or struck against a reference that shares nothing with its code
/// (<c>khaos_defects_oracle.R</c>, and <c>khaos_conditionals_oracle.R</c> for the frozen-state rows; E19's build brief, slice 0). These
/// are corrections, not parity variants: the engine takes the corrected reference.
/// <list type="bullet">
/// <item><b>Confirmed:</b> the g-prior's structure ratio without the prior's normalizer, also in both log-evidence cores (which nothing in
/// 2.1.1 calls); g₀²'s conditional without its data term and one power of g₀² short; the birth's forward probability; random_partition
/// and the C(d, q) count; the variable swap onto an active variable; the death's normal approximation.</item>
/// <item><b>Struck for the code:</b> SC-1's printed inversions (every factor of the code's birth, death and degree mutation points the
/// Metropolis–Hastings way) and the paper's slips (the code reads P₁(x) = x, coins over p, degrees up to d_max).</item>
/// <item><b>Found here:</b> the death's prior ratio off by K/(K + 1) against the birth's; the move probabilities at M = 0; σ²'s draw
/// given β without β's prior; the g-prior's current state left at the g₀² of its last acceptance; the degree draw failing where the
/// order is the largest degree; the degree weights indexed from q.</item>
/// </list>
/// </summary>
public sealed class KhaosDefectParityTests(ITestOutputHelper output)
{
    /// <summary>The seeded chain's data: 100 cases on [0, 1]⁵, Ishigami's function on the first three inputs and noise 0.5·N(0, 1).</summary>
    private static readonly Lazy<JsonElement?> Fixture = new(() =>
        GbassOracle.Run("khaos_defects_oracle.R", KhaosReadoutParityTests.Data(100, 20261105), "1500", "5"));

    private static JsonElement? Conditionals => KhaosOracle.Conditionals.Value;

    [Fact]
    public void GPriorStructureRatio_LeavesOutThePriorsNormalizer()
    {
        if (Fixture.Value is not { } d || Conditionals is not { } c) return;
        // At the frozen state: the ratio as written is the marginals' difference less ½Δln|S₀⁻¹|, to rounding.
        JsonElement m = c.GetProperty("marginal");
        JsonElement[] st = [.. m.GetProperty("structures").EnumerateArray()];
        double[] code = [.. m.GetProperty("births").EnumerateArray().Select(b => b.GetProperty("gprior_code").GetDouble())];
        for (int k = 1; k < st.Length; k++)
        {
            double dDense = st[k].GetProperty("gprior").GetProperty("dense").GetDouble() - st[k - 1].GetProperty("gprior").GetProperty("dense").GetDouble();
            double dHalf = st[k].GetProperty("gprior").GetProperty("half_ldet_prior_precision").GetDouble() - st[k - 1].GetProperty("gprior").GetProperty("half_ldet_prior_precision").GetDouble();
            output.WriteLine($"g-prior birth {k - 1} → {k}: as written {code[k - 1]:F6}, the marginals' difference {dDense:F6}, ½Δln|S₀⁻¹| {dHalf:F6}");
            Assert.Equal(dDense - dHalf, code[k - 1], 10);
            Assert.True(Math.Abs(dHalf) > 0.1);
        }
        // On an orthogonal design the code's penalty per term is ½ln n + ½ln(1 + 1/(g₀²g²)), about ½ln n whatever g₀² and g are; the
        // stated model's is ½ln(1 + g₀²g²).
        foreach (JsonElement e in d.GetProperty("normalizer").EnumerateArray())
        {
            double g0sq = e.GetProperty("g0sq").GetDouble(), g = e.GetProperty("g_m").GetDouble(), x = g0sq * g * g;
            double halfLogN = 0.5 * Math.Log(64);
            Assert.Equal(halfLogN, e.GetProperty("half_log_n").GetDouble(), 14);
            Assert.Equal(halfLogN + 0.5 * Math.Log(1 + 1 / x), e.GetProperty("code").GetDouble(), 10);
            Assert.Equal(0.5 * Math.Log(1 + x), e.GetProperty("reference").GetDouble(), 9);
            output.WriteLine($"orthogonal design, g₀² {g0sq}, g {g:F4}: the code's penalty {e.GetProperty("code").GetDouble():F4}, the model's {e.GetProperty("reference").GetDouble():F4} (½ln n {halfLogN:F4})");
        }
        // The log-evidence cores: the fast core omits the same ½Δln|S₀⁻¹| on that design; the full core also reads the least-squares
        // residual for the posterior's quadratic form. Nothing in 2.1.1 calls either.
        JsonElement[] ev = [.. d.GetProperty("evidence").GetProperty("structures").EnumerateArray()];
        for (int k = 1; k < ev.Length; k++)
        {
            double Dk(string key) => ev[k].GetProperty(key).GetDouble() - ev[k - 1].GetProperty(key).GetDouble();
            double half = Dk("half_ldet_prior_precision");
            Assert.Equal(Dk("dense") - half, Dk("fast"), 10);
            output.WriteLine($"log evidence {k - 1} → {k}: dense {Dk("dense"):F4}, fast core {Dk("fast"):F4}, full core {Dk("full"):F4}, ½Δln|S₀⁻¹| {half:F4}");
            Assert.True(Math.Abs(Dk("full") - (Dk("dense") - half)) > 0.5);
        }
        Assert.Equal(0, d.GetProperty("evidence").GetProperty("callers").GetArrayLength());
    }

    [Fact]
    public void G0SqConditional_HasNoDataTerm_AndOnePowerShort()
    {
        if (Conditionals is not { } c) return;
        JsonElement[] grid = [.. c.GetProperty("g0sq").GetProperty("grid").EnumerateArray()];
        double V(JsonElement e, string k) => e.GetProperty(k).GetDouble();
        double withoutData = KhaosOracle.Spread(grid.Select(e => V(e, "joint_integrated") - V(e, "code")));
        double withData = KhaosOracle.Spread(grid.Select(e => V(e, "joint_integrated") - V(e, "code") - V(e, "data_term")));
        double withPower = KhaosOracle.Spread(grid.Select(e => V(e, "joint_integrated") - V(e, "code") - V(e, "data_term") + Math.Log(V(e, "theta"))));
        double given = KhaosOracle.Spread(grid.Select(e => V(e, "joint_given") - V(e, "code")));
        output.WriteLine($"dgsq_full against the joint over g₀² from {V(grid[0], "theta")} to {V(grid[^1], "theta")}: spreads {withoutData:F3} nats as written, " +
                         $"{withData:F3} with the data term, {withPower:G3} with the data term and a further 1/g₀²; against the joint given β and σ², {given:F1}");
        Assert.True(withoutData > 1 && withData > 1 && given > 1);
        Assert.True(withPower < 1e-9);
    }

    [Fact]
    public void BirthForwardProbability_DiffersFromBothWrittenForms()
    {
        if (Conditionals is not { } c) return;
        JsonElement[] cases = [.. c.GetProperty("birth").EnumerateArray()];
        foreach (JsonElement bc in cases)
        {
            int order = bc.GetProperty("order").GetInt32(), p = bc.GetProperty("eta").GetArrayLength();
            double worstBirth = 0, worstDeath = 0, least = 0;
            bool anyRedraw = false;
            foreach (JsonElement e in bc.GetProperty("chis").EnumerateArray())
            {
                double exact = Math.Log(e.GetProperty("exact").GetDouble()), noRedraw = e.GetProperty("no_redraw").GetDouble();
                double death = e.GetProperty("death_reverse").GetDouble();
                JsonElement birth = e.GetProperty("code_birth");
                worstBirth = Math.Max(worstBirth, Math.Abs(noRedraw - exact));
                worstDeath = Math.Max(worstDeath, Math.Abs(death - exact));
                least = Math.Min(least, birth.GetProperty("least").GetDouble() - noRedraw);
                if (birth.GetProperty("least").GetDouble() < noRedraw - 1e-12) anyRedraw = true;
                // When every invalid draw is q = 0 (order = p), the birth's term never moves off the q₀-marginal: an all-zero draw adds ln 1.
                if (order >= p) Assert.Equal(noRedraw, birth.GetProperty("greatest").GetDouble(), 12);
            }
            output.WriteLine($"birth at order {order} of p {p}: |ln forward| off the exact by up to {worstBirth:F3} as the birth writes it without a redraw, " +
                             $"{worstDeath:F3} as the death writes it; redraws past q_max lower the birth's term by up to {Math.Max(0.0, -least):F2} nats");
            Assert.True(worstBirth > 0.01 && worstDeath > 0.01);
            Assert.Equal(order < p, anyRedraw);
        }
    }

    [Fact]
    public void RandomPartition_IsNotUniform_AndTheRatiosCountNeither()
    {
        if (Fixture.Value is not { } d || Conditionals is not { } c) return;
        foreach (JsonElement e in d.GetProperty("partition").EnumerateArray())
        {
            int dd = e.GetProperty("d").GetInt32(), q = e.GetProperty("q").GetInt32();
            double[] law = GbassOracle.Numbers(e.GetProperty("law"));
            double uniform = e.GetProperty("uniform").GetDouble(), count = e.GetProperty("code_count").GetDouble();
            Assert.Equal(1.0 / KhaosOracle.Choose(dd - 1, q - 1), uniform, 14);
            Assert.Equal(1.0 / KhaosOracle.Choose(dd, q), count, 14);
            double tv = 0.5 * law.Sum(v => Math.Abs(v - uniform));
            output.WriteLine($"random_partition({dd}, {q}): law {law.Min():G4} to {law.Max():G4}, uniform {uniform:G4}, total variation {tv:G3}; the ratios' 1/C(d, q) = {count:G4}");
            // One cut or none draws uniformly; two or more weight compositions with tied cuts down.
            if (q <= 2) Assert.True(tv < 1e-14); else Assert.True(tv > 0.05);
            Assert.All(law, v => Assert.NotEqual(count, v, 6));
        }
        // The degree mutation's ratio as written reads C(d, q): off the uniform compositions' ratio, and off the proposal's own where its
        // law is uniform (one cut).
        foreach (JsonElement m in c.GetProperty("mutate_degree").EnumerateArray())
        {
            double code = m.GetProperty("code").GetDouble(), uniform = m.GetProperty("uniform").GetDouble();
            Assert.True(Math.Abs(code - uniform) > 0.1);
            if (m.GetProperty("q").GetInt32() <= 2) Assert.Equal(uniform, m.GetProperty("actual_least").GetDouble(), 12);
        }
    }

    [Fact]
    public void VariableSwap_ProposesAVariableAlreadyInTheTerm()
    {
        if (Fixture.Value is not { } d || Conditionals is not { } c) return;
        foreach (JsonElement m in c.GetProperty("mutate_vars").EnumerateArray())
        {
            int[] vars = [.. GbassOracle.Numbers(m.GetProperty("vars")).Select(v => (int)v)];
            foreach (var position in m.GetProperty("swaps").EnumerateArray().GroupBy(s => s.GetProperty("position").GetInt32()))
            {
                double share = position.Where(s => s.GetProperty("duplicate").GetBoolean()).Sum(s => s.GetProperty("forward").GetDouble());
                output.WriteLine($"swapping out x{vars[position.Key - 1]} of {{{string.Join(", ", vars.Select(v => "x" + v))}}}: a repeated variable with probability {share:F3}");
                Assert.True(share > 0.1);
            }
        }
        JsonElement chain = d.GetProperty("swap_chain");
        int terms = chain.GetProperty("terms").GetInt32(), repeated = chain.GetProperty("repeated").GetInt32();
        output.WriteLine($"adaptive_khaos_gprior, {chain.GetProperty("nmcmc").GetInt32()} iterations ({chain.GetProperty("seconds").GetDouble():F1} s): {repeated} of {terms} kept terms repeat a variable " +
                         $"({(double)repeated / terms:P1}); {chain.GetProperty("draws_with_repeat").GetInt32()} of {chain.GetProperty("draws").GetInt32()} draws hold one; " +
                         $"swaps accepted {chain.GetProperty("swaps_accepted").GetInt32()} of {chain.GetProperty("swaps_proposed").GetInt32()}");
        Assert.True(repeated > terms / 10);
    }

    [Fact]
    public void DeathsNormalApproximation_AgainstThePoissonBinomial()
    {
        if (Fixture.Value is not { } d) return;
        double worst = 0;
        foreach (JsonElement e in d.GetProperty("normal").EnumerateArray())
        {
            int order = e.GetProperty("order").GetInt32();
            foreach (JsonElement q in e.GetProperty("per_q0").EnumerateArray())
            {
                double[] w = GbassOracle.Numbers(q.GetProperty("weights"));
                double[] pmf = KhaosOracle.PoissonBinomial(w);
                double exact = pmf.Skip(order + 1).Sum(), normal = q.GetProperty("normal_over").GetDouble();
                Assert.Equal(exact, q.GetProperty("exact_over").GetDouble(), 13);
                Assert.Equal(pmf[0], q.GetProperty("zero").GetDouble(), 13);
                worst = Math.Max(worst, Math.Abs(normal / exact - 1));
                output.WriteLine($"η = {string.Join(",", GbassOracle.Numbers(e.GetProperty("eta")))}, q_max {order}, q₀ {q.GetProperty("q0").GetInt32()}: " +
                                 $"P(q > q_max) {exact:G4} exact, {normal:G4} by the death's normal");
            }
        }
        Assert.True(worst > 0.05);
    }

    [Fact]
    public void Sigma2Floor_IsRecorded()
    {
        if (Fixture.Value is not { } d) return;
        JsonElement f = d.GetProperty("floor");
        Assert.Equal(0, f.GetProperty("gprior_default").GetDouble());
        Assert.Equal(0, f.GetProperty("ridge_default").GetDouble());
        Assert.StartsWith("max(s2_lower", f.GetProperty("gprior_clamp").GetString());
        Assert.StartsWith("max(s2_lower", f.GetProperty("ridge_clamp").GetString());
        output.WriteLine($"σ²'s floor: s2_lower defaults to 0 in both samplers, applied as {f.GetProperty("gprior_clamp").GetString()}; not ported");
    }

    [Fact]
    public void PapersSlips_AreStruck_InTheCodesOwnReading()
    {
        if (Fixture.Value is not { } d) return;
        JsonElement t = d.GetProperty("typos");
        // P₁(x) = x, not 1.
        Assert.Equal(GbassOracle.Numbers(t.GetProperty("p1").GetProperty("x")), GbassOracle.Numbers(t.GetProperty("p1").GetProperty("code")));
        // The coins run over the p inputs, not M.
        Assert.Equal("stats::rbinom(p, 1, wts)", t.GetProperty("coin_draw").GetString());
        // The degree is drawn from {q, …, d_max}, not {q, …, p}.
        JsonElement ds = t.GetProperty("degree_support");
        Assert.Equal(Enumerable.Range(ds.GetProperty("q").GetInt32(), ds.GetProperty("degree").GetInt32() - ds.GetProperty("q").GetInt32() + 1).Select(v => (double)v),
                     GbassOracle.Numbers(ds.GetProperty("support")));
        // The partition count is C(d, q): neither the paper's binom(d − 1, d − 1) nor the uniform binom(d − 1, q − 1).
        Assert.Contains("lchoose(dtot.curr, nint.curr)", t.GetProperty("mutate_count").GetString());
        // Found: the degree weights index from q, (d − q + 1)^(−s_d), where the paper writes d^(−s_d); both directions read the same.
        JsonElement dw = t.GetProperty("degree_weights");
        double[] code = GbassOracle.Numbers(dw.GetProperty("code"));
        double[] fromQ = [.. Enumerable.Range(1, code.Length).Select(i => 1.0 / i)];
        Assert.Equal(fromQ.Select(v => v / fromQ.Sum()), code, (u, v) => Math.Abs(u - v) < 1e-14);
        output.WriteLine($"degree weights at q 2, d_max 6, s_d 1: as written {string.Join(", ", code.Select(v => v.ToString("F3")))}; the paper's d^(−1) " +
                         $"{string.Join(", ", GbassOracle.Numbers(dw.GetProperty("paper")).Select(v => v.ToString("F3")))}");
    }

    /// <summary>SC-1's table: each factor's exponent in the Metropolis–Hastings ratio and in the paper's print (+1 numerator, −1
    /// denominator). A weight is the degree draw's probability d_probs = d^(−s_d)/c_d; the count is the partition's.</summary>
    private static readonly (string Move, string Factor, int MetropolisHastings, int Printed)[] Sc1 =
    [
        ("birth", "d_weight", -1, 1), ("birth", "count", 1, 1), ("birth", "terms", -1, -1), ("birth", "coins", -1, 1),
        ("death", "d_weight", 1, -1), ("death", "count", -1, -1), ("death", "terms", 1, -1), ("death", "coins", 1, 1),
        ("mutate_degree", "d_weight_curr", 1, -1), ("mutate_degree", "d_weight_cand", -1, 1), ("mutate_degree", "count_curr", -1, 1), ("mutate_degree", "count_cand", 1, -1),
    ];

    [Fact]
    public void Sc1PrintedInversions_AreStruck_ForTheCode()
    {
        if (Fixture.Value is not { } d) return;
        JsonElement dirs = d.GetProperty("directions");
        int inverted = 0;
        foreach ((string move, string factor, int mh, int printed) in Sc1)
        {
            double code = dirs.GetProperty(move).GetProperty(factor).GetDouble();
            Assert.Equal(mh, code, 10);
            if (printed != mh) inverted++;
        }
        output.WriteLine($"Every factor of the code's birth, death and degree mutation takes its Metropolis–Hastings exponent; {inverted} of {Sc1.Length} are inverted in the print");
        Assert.Equal(8, inverted);
    }

    [Fact]
    public void DeathPriorRatio_IsNotTheBirthsReciprocal()
    {
        if (Conditionals is not { } c) return;
        foreach (JsonElement p in c.GetProperty("prior").EnumerateArray())
            foreach (JsonElement e in p.GetProperty("ratios").EnumerateArray())
            {
                int k = e.GetProperty("K").GetInt32();
                output.WriteLine($"K = {k}: birth {e.GetProperty("birth").GetDouble():F4} + death {e.GetProperty("death").GetDouble():F4} = {e.GetProperty("sum").GetDouble():F4}, ln(K/(K + 1)) = {Math.Log(k / (k + 1.0)):F4}");
                Assert.Equal(Math.Log(k / (k + 1.0)), e.GetProperty("sum").GetDouble(), 12);
            }
    }

    [Fact]
    public void MoveProbabilities_AtNoTerms_AreNotThoseTheRatiosRead()
    {
        if (Conditionals is not { } c) return;
        foreach (JsonElement p in c.GetProperty("prior").EnumerateArray())
        {
            double[] at0 = GbassOracle.Numbers(p.GetProperty("moves_at_0")), at1 = GbassOracle.Numbers(p.GetProperty("moves_at_1"));
            Assert.Equal(new[] { 1.0, 0.0, 0.0 }, at0);
            Assert.InRange(at1[1], 1.0 / 3 - 0.02, 1.0 / 3 + 0.02);
            // The exact move parts: a birth from M = 0 is forced, so its forward is 1 and its reverse death 1/3.
            double birthExact = Math.Log(1.0 / 3) - Math.Log(1.0), deathExact = -birthExact;
            double birthCode = p.GetProperty("birth_from_0_code").GetDouble(), deathCode = p.GetProperty("death_to_0_code").GetDouble();
            output.WriteLine($"move parts at M = 0: birth as written {birthCode:F4}, exact {birthExact:F4}; death to M = 0 as written {deathCode:F4}, exact {deathExact:F4}");
            Assert.Equal(0.0, birthCode, 14);
            Assert.Equal(0.0, deathCode, 14);
        }
    }

    [Fact]
    public void GPriorCurrentState_IsStale_AfterG0SqMoves()
    {
        if (Fixture.Value is not { } d) return;
        JsonElement s = d.GetProperty("swap_chain").GetProperty("stale");
        double mean = s.GetProperty("mean_abs").GetDouble(), greatest = s.GetProperty("greatest").GetDouble();
        double overOne = s.GetProperty("over_one").GetDouble();
        output.WriteLine($"the g-prior's current state as its ratios read it, against its value at the iteration's g₀², over {s.GetProperty("iterations").GetInt32()} iterations: " +
                         $"mean |gap| {mean:F2} nats, greatest {greatest:F2}, over 1 nat on {overOne:P1} of them");
        Assert.True(greatest > 1 && overOne > 0.1);
    }

    [Fact]
    public void DegreeDraw_FailsWhereTheOrderIsTheLargestDegree()
    {
        if (Fixture.Value is not { } d) return;
        JsonElement s = d.GetProperty("degree_draw");
        string written = s.GetProperty("written").GetString()!;
        output.WriteLine($"the degree draw at q = d_max = 2: as written \"{written}\"; corrected {s.GetProperty("corrected").GetDouble()}");
        Assert.Contains("incorrect number of probabilities", written);
        Assert.Equal(2.0, s.GetProperty("corrected").GetDouble());
    }

    [Fact]
    public void Sigma2Draw_GivenBeta_LeavesOutBetasPrior()
    {
        if (Conditionals is not { } c) return;
        JsonElement s = c.GetProperty("sigma2");
        double[] grid = GbassOracle.Numbers(s.GetProperty("grid"));
        foreach (string prior in new[] { "gprior", "ridge" })
        {
            JsonElement p = s.GetProperty(prior);
            double[] joint = GbassOracle.Numbers(p.GetProperty("joint"));
            double spread = KhaosOracle.Spread(grid.Select((x, i) => joint[i] - KhaosOracle.LogInverseGammaKernel(x, p.GetProperty("code_shape").GetDouble(), p.GetProperty("code_rate").GetDouble())));
            output.WriteLine($"σ²'s draw as written under the {prior}: IG({p.GetProperty("code_shape").GetDouble()}, {p.GetProperty("code_rate").GetDouble():G5}) against the conditional " +
                             $"IG({p.GetProperty("stated_shape").GetDouble()}, {p.GetProperty("stated_rate").GetDouble():G5}); the joint less it spreads {spread:F2} nats");
            Assert.Equal(p.GetProperty("stated_shape").GetDouble() - (KhaosOracle.StateBeta.Length / 2.0), p.GetProperty("code_shape").GetDouble(), 12);
            Assert.True(spread > 1);
        }
    }
}
