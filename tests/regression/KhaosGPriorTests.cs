using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.BassKernelEnumeration;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's modified g-prior (the KHAOS brief's slice 2): <see cref="ModifiedGForm"/> through the design's prior-precision builder
/// (<see cref="BassPriorPrecision"/>), with g₀² = 1/τ under <see cref="TauPrior.Gamma"/> and moved by τ's slice step.
/// <list type="bullet">
/// <item><b>The prior on τ:</b> the Gamma on τ itself, its density, draw, refusals and serialized form.</item>
/// <item><b>The factor</b> g_m = (1 + q(d + q − 2))^(−ζ/2).</item>
/// <item><b>The target against the model's statement:</b> at frozen structures and a grid of g₀², the engine's log marginal against
/// a dense multivariate t, y | g₀² ~ t with scale I + B S₀ Bᵀ, S₀ = g₀² D_g(BᵀB)⁻¹D_g, built from the columns alone; so g₀²'s
/// conditional, which τ's slice step reads, is the joint's.</item>
/// <item><b>ζ = 0 is the g-form</b> with the intercept inside under the zero centre: its pencil against <see cref="CollapsedGPrior"/>
/// at every enumerated structure and every τ.</item>
/// <item><b>The exact law with g₀² integrated:</b> on slice 1's enumerable space, each structure's weight by Simpson's rule over
/// u = ln τ from the target, and the production chain with g₀² drawn against it: P(M), each inclusion and P(g₀² &lt; c).</item>
/// <item><b>Geweke with g₀² drawn,</b> behind <see cref="BassQualification"/>.</item>
/// </list>
/// </summary>
public sealed class KhaosGPriorTests(ITestOutputHelper output)
{
    private const int P = 4;

    // khaos's adaptive_khaos_gprior defaults: λ ~ Gamma(4, rate 4/n), g₀² ~ IG(0.001, 1,000), ζ = 1, σ² Jeffreys, on slice 1's
    // exact-law space (d_max 3, q_max 2, M ≤ 2).
    internal static BassConfig Config(double zeta = 1.0, TauPrior? scale = null, int maxBasis = 2, int cases = 40, Sigma2Prior? sigma2 = null) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = maxBasis, Legendre = new LegendreBasis { MaxDegree = 3 } },
        Count = new BassCountPrior(4, 4.0 / cases),
        Coefficients = new BassCoefficientPrior(new ModifiedGForm(zeta), scale ?? TauPrior.Gamma(0.001, 1000)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = sigma2 ?? Sigma2Prior.Jeffreys,
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    // ---- τ's prior ----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TauGamma_IsOnTauItself()
    {
        TauPrior prior = TauPrior.Gamma(2.5, 7.0);
        ResolvedTauPrior a = prior.For(1.0), b = prior.For(500.0);
        Assert.Equal((2.5, 7.0), (a.Gamma!.Value.Shape, a.Gamma.Value.Rate));
        Assert.Equal((2.5, 7.0), (b.Gamma!.Value.Shape, b.Gamma.Value.Rate));   // whatever τ₀
        foreach (double u in new[] { -6.0, -1.0, 0.0, 0.7, 2.0 })
        {
            double expected = 2.5 * Math.Log(7.0) - SpecialFunctions.LogGamma(2.5) + 2.5 * u - 7.0 * Math.Exp(u);
            Assert.Equal(expected, a.LogDensity(u), 12);
        }
        var rng = new Xoshiro256PlusPlus(5);
        double mean = Enumerable.Range(0, 200_000).Average(_ => Math.Exp(a.SampleLogTau(rng)));
        Assert.True(Math.Abs(mean - 2.5 / 7.0) < 0.005, $"the draws' mean {mean:F4} against {2.5 / 7.0:F4}");

        Assert.Throws<ArgumentOutOfRangeException>(() => TauPrior.Gamma(0.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => TauPrior.Gamma(1.0, double.PositiveInfinity));
        Assert.Throws<ArgumentException>(() => new TauPrior(TauPriorKind.UnitInformationGamma, 0.5, rate: 2.0));
        Assert.Equal("Gamma(2.5, rate 7)", prior.ToString());

        // The rate is written only where it is declared, so every other kind's form stands.
        Assert.DoesNotContain("Rate", JsonSerializer.Serialize(TauPrior.UnitInformationGamma()));
        Assert.DoesNotContain("Rate", JsonSerializer.Serialize(TauPrior.Held(0.5)));
        string json = JsonSerializer.Serialize(prior);
        Assert.Contains("\"Rate\":7", json);
        Assert.Equal(prior, JsonSerializer.Deserialize<TauPrior>(json));
    }

    // ---- the factor ---------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Factor_ReadsOrderAndDegree()
    {
        var form = new ModifiedGForm(1.0);
        Assert.Equal(1.0, form.Factor(new BassTerm(legendre: [new(0, 1)])));                          // q 1, d 1: 1 + 1·0
        Assert.Equal(Math.Pow(8.0, -0.5), form.Factor(new BassTerm(legendre: [new(2, 8)])), 15);       // q 1, d 8: 1 + 7
        Assert.Equal(Math.Pow(11.0, -0.5), form.Factor(new BassTerm(legendre: [new(0, 2), new(3, 3)])), 15);   // q 2, d 5: 1 + 2·5
        Assert.Equal(Math.Pow(11.0, -1.25), new ModifiedGForm(2.5).Factor(new BassTerm(legendre: [new(0, 2), new(3, 3)])), 15);
        Assert.Equal(1.0, new ModifiedGForm(0.0).Factor(new BassTerm(legendre: [new(0, 2), new(3, 3)])));
        Assert.Throws<ArgumentException>(() => form.Factor(new BassTerm([new HingeFactor(0, 1, 3)])));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ModifiedGForm(-0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ModifiedGForm(double.NaN));
    }

    [Fact]
    public void Config_RefusesAHingeBasis_AndKeepsTheForm()
    {
        BassData data = KhaosChainTests.Data();
        Assert.Throws<ArgumentException>(() => BassModel.Create(Config() with { Structure = new BassStructureSpec { MaxInteraction = 2 } }, data));
        BassModel model = BassModel.Create(Config(2.0), data);
        Assert.Equal(2.0, ((ModifiedGForm)model.Config.Coefficients.Form).Zeta);
        Assert.Contains(model.Resolved.AutoPicks, p => p.Setting == "TauPrior" && p.Value == "Gamma(0.001, rate 1000)");
        string json = model.Config.ToJson();
        Assert.Contains("\"ModifiedG\"", json);
        BassConfig back = BassConfig.FromJson(json);
        Assert.Equal(2.0, ((ModifiedGForm)back.Coefficients.Form).Zeta);
        Assert.Equal(TauPrior.Gamma(0.001, 1000), back.Coefficients.Scale);
        Assert.Equal(json, back.ToJson());
    }

    // ---- the target against the model's statement ---------------------------------------------------------------------------------

    // ln ∫ N(y; 0, σ²C) dσ²/σ², C = I + B S₀ Bᵀ, up to a constant in n alone: −½ ln|C| − (n/2) ln yᵀC⁻¹y.
    private static double DenseLogMarginal(double[][] columns, double[] g, double g0sq, double[] y)
    {
        int n = y.Length, k = columns.Length;
        var btb = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++) btb[i, j] = Enumerable.Range(0, n).Sum(r => columns[i][r] * columns[j][r]);
        double[,] inverse = Inverse(btb);
        var s0 = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++) s0[i, j] = g0sq * g[i] * inverse[i, j] * g[j];
        var c = new double[n, n];
        for (int r = 0; r < n; r++)
            for (int t = 0; t < n; t++)
            {
                double sum = r == t ? 1.0 : 0.0;
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < k; j++) sum += columns[i][r] * s0[i, j] * columns[j][t];
                c[r, t] = sum;
            }
        double[,] l = Cholesky(c);
        double halfLogDet = Enumerable.Range(0, n).Sum(i => Math.Log(l[i, i]));
        double[] z = Forward(l, y);
        return -halfLogDet - 0.5 * n * Math.Log(z.Sum(v => v * v));
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            double d = a[j, j];
            for (int k = 0; k < j; k++) d -= l[j, k] * l[j, k];
            l[j, j] = Math.Sqrt(d);
            for (int i = j + 1; i < n; i++)
            {
                double s = a[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = s / l[j, j];
            }
        }
        return l;
    }

    private static double[] Forward(double[,] l, double[] b)
    {
        int n = b.Length;
        var z = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = b[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * z[k];
            z[i] = s / l[i, i];
        }
        return z;
    }

    private static double[,] Inverse(double[,] a)
    {
        int n = a.GetLength(0);
        double[,] l = Cholesky(a);
        var inverse = new double[n, n];
        for (int c = 0; c < n; c++)
        {
            double[] z = Forward(l, [.. Enumerable.Range(0, n).Select(i => i == c ? 1.0 : 0.0)]);
            for (int i = n - 1; i >= 0; i--)
            {
                double s = z[i];
                for (int k = i + 1; k < n; k++) s -= l[k, i] * inverse[k, c];
                inverse[i, c] = s / l[i, i];
            }
        }
        return inverse;
    }

    // The structure's columns as evaluated, the intercept's first, and their factors (1 for the intercept).
    private static (double[][] Columns, double[] G) Columns(BassModel model, IReadOnlyList<BassTerm> terms, ModifiedGForm form)
    {
        int n = model.Inputs.Cases.Cases;
        var columns = new List<double[]> { Enumerable.Repeat(1.0, n).ToArray() };
        var g = new List<double> { 1.0 };
        foreach (BassTerm t in terms)
        {
            columns.Add([.. Enumerable.Range(0, n).Select(i => t.Legendre.Aggregate(1.0, (acc, f) => acc * ShiftedLegendre.Value(f.Degree, model.Inputs.Cases.Continuous(f.Variable)[i])))]);
            g.Add(form.Factor(t));
        }
        return ([.. columns], [.. g]);
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.5)]
    [InlineData(0.0)]
    public void Target_IsTheDenseMarginal_AndG0sqsConditionalIsTheJoints(double zeta)
    {
        BassData data = KhaosChainTests.Data();
        double[] y = [.. ((ScalarResponse)data.Response).Values];
        BassModel model = BassModel.Create(Config(zeta), data);
        var form = (ModifiedGForm)model.Config.Coefficients.Form;
        List<BassTerm> terms = [.. model.Prior.Legendre!.Terms()];
        var structures = new List<BassTerm[]> { Array.Empty<BassTerm>(), new[] { terms[0] }, new[] { terms[7] }, new[] { terms[2], terms[11] }, new[] { terms[5], terms[20] }, new[] { terms[13], terms[29] } };
        double[] g0sqs = [0.5, 5.0, 50.0, 500.0, 5000.0, 1e6];
        double? offset = null;
        double worstMarginal = 0.0, worstConditional = 0.0;
        foreach (BassTerm[] s in structures)
        {
            (double[][] columns, double[] g) = Columns(model, s, form);
            BassState state = BassMixing.Loaded(model, s, 0.0, out _) ?? throw new InvalidOperationException("outside F");
            double? conditionalOffset = null;
            foreach (double g0sq in g0sqs)
            {
                double u = -Math.Log(g0sq);
                BassEvaluation e = model.Target.Evaluate(state.WithLogTau(u));
                double dense = DenseLogMarginal(columns, g, g0sq, y);
                double gap = e.LogMarginal - dense;
                offset ??= gap;
                worstMarginal = Math.Max(worstMarginal, Math.Abs(gap - offset.Value));
                // g₀²'s conditional in u = ln τ, as τ's slice step reads it: the joint's, the Gamma on τ with its Jacobian.
                double joint = model.Target.LogTauDensity(u) + dense;
                double conditional = e.LogPrior + e.LogLikelihood(1.0) - joint;
                conditionalOffset ??= conditional;
                worstConditional = Math.Max(worstConditional, Math.Abs(conditional - conditionalOffset.Value));
            }
        }
        output.WriteLine($"ζ {zeta}: {structures.Count} structures × {g0sqs.Length} g₀² from 0.5 to 10⁶: the log marginal against the dense t, " +
                         $"worst gap after one constant {worstMarginal:E2}; g₀²'s conditional against the joint, per structure, {worstConditional:E2}");
        Assert.True(worstMarginal < 1e-9, $"the marginal differs from the dense t by {worstMarginal:E2}");
        Assert.True(worstConditional < 1e-9, $"g₀²'s conditional differs from the joint's by {worstConditional:E2}");
    }

    // A vague prior on τ draws subnormal τ about half the time (Gamma(0.001, rate 1,000) puts P(τ < 10⁻³⁰⁰) near ½), where ln(1 + bλ/τ)
    // overflowed and the marginal read −∞ though the density is finite, so a chain started from such a draw could not move. There
    // ln|bA + τD| − ln|τD| tends to Σ ln bλⱼ − k ln τ, so the log marginal rises by ½ per unit of u = ln τ per direction.
    [Fact]
    public void Target_IsFinite_AtSubnormalTau()
    {
        BassModel model = BassModel.Create(Config(), KhaosChainTests.Data());
        List<BassTerm> terms = [.. model.Prior.Legendre!.Terms()];
        foreach (BassTerm[] s in new[] { Array.Empty<BassTerm>(), new[] { terms[0], terms[7] } })
        {
            BassState state = BassMixing.Loaded(model, s, 0.0, out _)!;
            int k = s.Length + 1;   // the intercept inside
            double previous = double.NaN;
            // Below about u = −708 τ is subnormal, so its spacing, 2⁻¹⁰⁷⁴, enters ln τ: 2·10⁻¹⁰ at −721.74, 1% at −740. The asymptote is
            // read where that spacing is negligible, and finiteness all the way down.
            foreach (double u in new[] { -700.0, -712.0, -721.74, -730.0, -740.0 })
            {
                BassEvaluation e = model.Target.Evaluate(state.WithLogTau(u));
                Assert.True(double.IsFinite(e.LogMarginal) && double.IsFinite(e.LogPrior), $"M {s.Length}, u {u}: marginal {e.LogMarginal}, prior {e.LogPrior}");
                if (u < -722.0) continue;
                previous = double.IsNaN(previous) ? e.LogMarginal - 0.5 * k * u : previous;
                Assert.Equal(previous, e.LogMarginal - 0.5 * k * u, 6);
            }
        }
        // Across the band where Goldberg's ln(1 + x) overflowed in its product, bλ/τ between 2.5·10³⁰⁵ and the doubles' end, a wall of −∞
        // that held a chain below it: the log marginal is smooth through it.
        BassState empty = BassMixing.Loaded(model, [], 0.0, out _)!;
        double[] us = [.. Enumerable.Range(0, 81).Select(i => -712.0 + 0.125 * i)];
        double[] m = [.. us.Select(u => model.Target.Evaluate(empty.WithLogTau(u)).LogMarginal)];
        for (int i = 1; i < us.Length; i++)
            Assert.True(double.IsFinite(m[i]) && Math.Abs(m[i] - m[i - 1] - 0.5 * 0.125) < 1e-6, $"u {us[i]}: {m[i]} after {m[i - 1]}");

        // A chain from a subnormal τ climbs to the posterior's τ.
        var chain = new BassChain(model, [], -721.74, new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        for (int t = 0; t < 300; t++) chain.Step();
        Assert.True(chain.Current.LogTau > -100.0, $"the chain's u is {chain.Current.LogTau:F1} after 300 sweeps from −721.74");
    }

    // ---- ζ = 0 is the g-form ------------------------------------------------------------------------------------------------------

    [Fact]
    public void ZetaZero_IsTheGForm_AtEveryStructureAndTau()
    {
        BassData data = KhaosChainTests.Data();
        BassModel modified = BassModel.Create(Config(0.0, TauPrior.Gamma(2.0, 0.1)), data);
        BassModel gForm = BassModel.Create(Config(0.0, TauPrior.Gamma(2.0, 0.1)) with { Coefficients = BassCoefficientPrior.G(TauPrior.Gamma(2.0, 0.1)) }, data);
        Space space = KhaosChainTests.Enumerate(BassModel.Create(Config(0.0, TauPrior.Held(0.05)), data));
        double worst = 0.0;
        int compared = 0;
        foreach (State s in space.States)
        {
            BassState a = BassMixing.Loaded(modified, s.Terms, 0.0, out _)!, b = BassMixing.Loaded(gForm, s.Terms, 0.0, out _)!;
            foreach (double u in new[] { -9.0, -5.0, -2.5, 0.0, 1.5 })
            {
                BassEvaluation ea = modified.Target.Evaluate(a.WithLogTau(u)), eb = gForm.Target.Evaluate(b.WithLogTau(u));
                double la = ea.LogPrior + ea.LogLikelihood(1.0), lb = eb.LogPrior + eb.LogLikelihood(1.0);
                worst = Math.Max(worst, Math.Abs(la - lb) / Math.Max(1.0, Math.Abs(lb)));
                compared++;
            }
        }
        output.WriteLine($"ζ = 0 through the pencil against the g-form's CollapsedGPrior: {space.States.Count} structures, {compared} (structure, τ) pairs, worst relative gap {worst:E2}");
        Assert.True(worst < 1e-11, $"ζ = 0 differs from the g-form by {worst:E2}");
    }

    // ---- the exact law with g₀² integrated ---------------------------------------------------------------------------------------

    /// <summary>g₀²'s cuts at which the exact law reads its posterior CDF, the R oracle's.</summary>
    internal static readonly double[] G0sqCuts = [300.0, 1000.0, 10000.0];

    /// <summary>
    /// The law over <paramref name="space"/>'s structures with g₀² integrated: each structure's log weight, ln ∫ e^{ℓ(S, u)} du over
    /// u = ln τ from the target's log density (the Gamma on τ with its Jacobian, the structure prior and the marginal), by Simpson's
    /// rule on 2,000 panels between g₀² = 10⁻⁶, the cuts and 10¹⁶, as the R oracle integrates; and each structure's P(g₀² &lt; c).
    /// </summary>
    internal static (double[] Pi, double[][] Below) IntegratedLaw(BassModel model, Space space)
    {
        double[] bounds = [.. new[] { 1e16 }.Concat(G0sqCuts.Reverse()).Append(1e-6).Select(g => -Math.Log(g))];   // ascending in u
        var log = new double[space.States.Count];
        var below = new double[space.States.Count][];
        for (int k = 0; k < log.Length; k++)
        {
            BassState state = BassMixing.Loaded(model, space.States[k].Terms, 0.0, out _)!;
            double F(double u)
            {
                BassEvaluation e = model.Target.Evaluate(state.WithLogTau(u));
                return e.LogPrior + e.LogLikelihood(1.0);
            }
            double top = Enumerable.Range(0, 2001).Max(i => F(bounds[0] + (bounds[^1] - bounds[0]) * i / 2000.0));
            double[] parts = [.. Enumerable.Range(0, bounds.Length - 1).Select(j => Simpson(F, bounds[j], bounds[j + 1], top))];
            double total = parts.Sum();
            log[k] = top + Math.Log(total);
            // g₀² < c is u > −ln c: the parts above each cut's bound.
            below[k] = [.. Enumerable.Range(0, G0sqCuts.Length).Select(c => parts.Skip(G0sqCuts.Length - c).Sum() / total)];
        }
        double max = log.Max();
        double[] pi = [.. log.Select(v => Math.Exp(v - max))];
        double sum = pi.Sum();
        return ([.. pi.Select(v => v / sum)], below);
    }

    private static double Simpson(Func<double, double> h, double a, double b, double top)
    {
        double step = (b - a) / 2000.0, sum = 0.0;
        for (int i = 0; i <= 2000; i++)
        {
            double f = Math.Exp(h(a + i * step) - top);
            sum += i == 0 || i == 2000 ? f : i % 2 == 1 ? 4.0 * f : 2.0 * f;
        }
        return step / 3.0 * sum;
    }

    [Theory]
    [InlineData(20261120)]
    public void Chain_WithG0sqDrawn_MeetsTheIntegratedLaw(int seed)
    {
        BassData data = KhaosChainTests.Data();
        BassModel model = BassModel.Create(Config(), data);
        Space space = KhaosChainTests.Enumerate(BassModel.Create(Config(scale: TauPrior.Held(0.05)), data));
        (double[] pi, double[][] below) = IntegratedLaw(model, space);
        var chain = new BassChain(model, [], Math.Log(1.0 / 40.0), new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
        const int burn = 5_000, steps = 400_000, batches = 40;
        for (int t = 0; t < burn; t++) chain.Step();
        var visits = new int[steps];
        var g0sq = new double[steps];
        for (int t = 0; t < steps; t++)
        {
            BassState s = chain.Step();
            visits[t] = space.Find(s.Terms);
            Assert.True(visits[t] >= 0, "the chain holds a structure outside the enumeration");
            g0sq[t] = Math.Exp(-s.LogTau);
        }
        var readings = new List<(string Name, double Exact, double[] Series)>();
        for (int m = 0; m <= 2; m++)
        {
            int mm = m;
            readings.Add(($"P(M = {m})", Enumerable.Range(0, pi.Length).Where(i => space.States[i].Terms.Length == mm).Sum(i => pi[i]),
                          [.. visits.Select(i => space.States[i].Terms.Length == mm ? 1.0 : 0.0)]));
        }
        for (int v = 0; v < P; v++)
        {
            int j = v;
            bool[] holds = [.. space.States.Select(st => st.Terms.Any(t => t.Legendre.Any(f => f.Variable == j)))];
            readings.Add(($"P(x{j} in)", Enumerable.Range(0, pi.Length).Where(i => holds[i]).Sum(i => pi[i]), [.. visits.Select(i => holds[i] ? 1.0 : 0.0)]));
        }
        for (int c = 0; c < G0sqCuts.Length; c++)
        {
            int cc = c;
            readings.Add(($"P(g0² < {G0sqCuts[c]})", Enumerable.Range(0, pi.Length).Sum(i => pi[i] * below[i][cc]), [.. g0sq.Select(g => g < G0sqCuts[cc] ? 1.0 : 0.0)]));
        }
        double worst = 0.0;
        foreach ((string name, double exact, double[] series) in readings)
        {
            double mean = series.Average(), se = KhaosChainTests.BatchError(series, batches);
            double z = se > 0.0 ? (mean - exact) / se : (mean == exact ? 0.0 : double.PositiveInfinity);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{name,-22} chain {mean:F4}  exact {exact:F4}  z {z,6:F2}");
        }
        output.WriteLine($"ζ 1, g₀² ~ IG(0.001, 1,000) drawn, {steps} steps: worst |z| {worst:F2}; moves {string.Join(", ", chain.Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}");
        Assert.True(worst < 4.5, $"the chain's worst |z| {worst:F2} against the integrated law");
    }

    // ---- Geweke with g₀² drawn ----------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;

    private static BassConfig GewekeConfig(double zeta) => Config(zeta, TauPrior.Gamma(2.0, 2.0), maxBasis: 3, cases: GewekeN, sigma2: Sigma2Prior.InverseGamma(3, 2)) with
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 3, Legendre = new LegendreBasis { MaxDegree = 3, DegreeMutation = 0.5 } },
    };

    private static BassData GewekeData(double[] y) => new(KhaosBasisTests.UnitData(GewekeN, 4, 77).Predictors, new ScalarResponse(y));

    private static readonly string[] GewekeNames = ["M", "ln sigma2", "ln g0sq", "atan f0", "atan f6", "x0 in", "x3 in", "degree sum", "z0^2"];

    private static double[] Functions(BassState state, double sigma2, double[] f, double[] y)
        => [state.Count, Math.Log(sigma2), -state.LogTau, Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
            state.Terms.Any(t => t.Legendre.Any(x => x.Variable == 0)) ? 1 : 0, state.Terms.Any(t => t.Legendre.Any(x => x.Variable == 3)) ? 1 : 0,
            state.Terms.Sum(t => t.LegendreDegree), (y[0] - f[0]) * (y[0] - f[0]) / sigma2];

    [Theory]
    [InlineData(1.0, 20261121)]
    [InlineData(1.0, 20261122)]
    [InlineData(0.0, 20261123)]
    public void Geweke_WithG0sqDrawn(double zeta, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 200_000, steps = 300_000, batches = 50;
        BassConfig config = GewekeConfig(zeta);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData([.. Enumerable.Range(0, GewekeN).Select(i => Math.Sin(i))]));
        var sampler = new BassPriorSampler(prior);

        (BassState S, double Sigma2, double[] F, double[] Y) Draw()
        {
            BassState state = sampler.DrawState(rng).State;
            var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
            double[] y = [.. Enumerable.Range(0, GewekeN).Select(i => fit.Fitted[i] + Math.Sqrt(fit.Sigma2) * Normal.SampleStandard(rng))];
            return (chain.Current, fit.Sigma2, fit.Fitted, y);
        }
        var marginal = new double[draws][];
        for (int k = 0; k < draws; k++)
        {
            (BassState s, double sigma2, double[] f, double[] y) = Draw();
            marginal[k] = Functions(s, sigma2, f, y);
        }

        (BassState S, double Sigma2, double[] F, double[] Y) theta = Draw();
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++)
        {
            double[] y = [.. Enumerable.Range(0, GewekeN).Select(i => theta.F[i] + Math.Sqrt(theta.Sigma2) * Normal.SampleStandard(rng))];
            BassModel model = BassModel.Create(config, GewekeData(y));
            var chain = new BassChain(model, theta.S.Terms, theta.S.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            chain.Step();
            BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, model.Target);
            successive[t] = Functions(chain.Current, fit.Sigma2, fit.Fitted, y);
            theta = (chain.Current, fit.Sigma2, fit.Fitted, y);
        }

        double worst = 0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(KhaosChainTests.BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"ζ {zeta}, g₀² ~ IG(2, 2) drawn, seed {seed}: worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }
}
