using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.BassKernelEnumeration;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's basis through the engine's sweep (the KHAOS brief's slice 1 gates), under the normalized ridge and KHAOS's raw ridge:
/// <list type="bullet">
/// <item><b>The law, exactly.</b> Every structure of a small Legendre term space enumerated (p = 4, d_max = 3, q_max = 2, M ≤ 2:
/// 30 terms, 466 structures), π from the production target, checked against a dense marginal computed here from the columns with
/// nothing of the design's; the kernel of each move enumerated from its statement (the birth's inclusion law with validity by
/// brute force, the death's uniform pick, each mutation's draw), with the production code's ratios, in detailed balance and
/// stationary for each move and the sweep.</item>
/// <item><b>The chain against the law.</b> The production chain's frequencies of M, of each variable's inclusion and of the
/// likeliest structures against π, by batch means.</item>
/// <item><b>Geweke's joint-distribution test</b> on twelve cases with proper priors, behind <see cref="BassQualification"/>.</item>
/// <item><b>The config,</b> its refusals and its serialized form.</item>
/// </list>
/// </summary>
public sealed class KhaosChainTests(ITestOutputHelper output)
{
    private const int P = 4, MaxDegree = 3, MaxOrder = 2, Cases = 40;
    private const double KhaosTau2 = 10.0;   // khaos's τ², slice 0's exact-law fixture

    // Slice 0's exact-law fixture, its inputs on [0, 1] exactly: y = 0.5ψ₁(x₀) + 0.35ψ₁(x₁)ψ₁(x₂) + ε, ε ~ N(0, 1).
    internal static BassData Data(int n = Cases, int seed = 20261109)
        => KhaosBasisTests.UnitData(n, P, seed, x => 0.5 * ShiftedLegendre.Value(1, x[0]) + 0.35 * ShiftedLegendre.Value(1, x[1]) * ShiftedLegendre.Value(1, x[2]), 1.0);

    internal static BassConfig Config(bool raw, int maxBasis = 2, Sigma2Prior? sigma2 = null, int seed = 1) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = MaxOrder, MaxBasis = maxBasis, Legendre = new LegendreBasis { MaxDegree = MaxDegree } },
        Count = new BassCountPrior(4, 40.0 / Cases),   // khaos's h₁ = 4, h₂ = 40/n
        Coefficients = raw ? BassCoefficientPrior.RawRidge(TauPrior.Held(1.0 / KhaosTau2)) : BassCoefficientPrior.Ridge(TauPrior.Held(0.5)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = sigma2 ?? Sigma2Prior.Jeffreys,       // khaos's g₁ = g₂ = 0
        Sampling = new BassSampling { Chains = 1, MasterSeed = seed, BurnIn = 0, Samples = 1 },
    };

    /// <summary>F of a Legendre fit at its held τ: every structure of at most M_max distinct terms the design and target admit.</summary>
    internal static Space Enumerate(BassModel model)
    {
        double logTau = Math.Log(model.Target.HeldTau!.Value);
        List<BassTerm> terms = [.. model.Prior.Legendre!.Terms()];
        var ids = new Dictionary<BassTerm, int>();
        for (int k = 0; k < terms.Count; k++) ids.Add(terms[k], k);
        var temperature = new Temperature(1.0, TemperingPath.Likelihood);
        var states = new List<State>();
        void Visit(int from, List<BassTerm> chosen)
        {
            if (BassMixing.Loaded(model, chosen, logTau, out _) is not { } loaded) return;
            BassTerm[] held = [.. loaded.Terms];
            BassEvaluation evaluation = model.Target.Evaluate(loaded);
            double logPrior = model.Prior.LogCollectionWeight(held.Length) + held.Sum(model.Prior.LogPrior);
            states.Add(new State(held, [.. held.Select(t => ids[t]).Order()], loaded, BassMixing.Counts(model, loaded),
                                 logPrior + evaluation.LogLikelihood(1.0), temperature.LogDensity(evaluation)));
            if (chosen.Count == model.MaxBasis) return;
            for (int i = from; i < terms.Count; i++)
            {
                chosen.Add(terms[i]);
                Visit(i + 1, chosen);
                chosen.RemoveAt(chosen.Count - 1);
            }
        }
        Visit(0, []);
        return new Space(model, 1.0, terms, states);
    }

    internal const string DegreeGroup = "degree mutation", SwapGroup = "variable swap";

    // ln D(d | q) as the settings state it: (d − q + 1)^(−s_d) normalized over {q, …, d_max}.
    private static double LogDegree(LegendreBasis settings, int q, int d)
        => -settings.DegreePenalty * Math.Log(d - q + 1)
           - Math.Log(Enumerable.Range(q, settings.MaxDegree - q + 1).Sum(e => Math.Pow(e - q + 1, -settings.DegreePenalty)));

    // Σ over valid χ′ of P(χ′ | η): the validity by enumeration, no Poisson-binomial.
    private static double Validity(double[] eta, int top)
    {
        double total = 0.0;
        for (int mask = 1; mask < 1 << eta.Length; mask++)
        {
            int count = 0;
            double product = 1.0;
            for (int j = 0; j < eta.Length; j++)
            {
                bool on = ((mask >> j) & 1) != 0;
                count += on ? 1 : 0;
                product *= on ? eta[j] : 1.0 - eta[j];
            }
            if (count <= top) total += product;
        }
        return total;
    }

    /// <summary>KHAOS's birth, death and change as the statement gives them: the birth's inclusion law Σ_{q₀} J(q₀)P(χ | q₀)/P(valid | q₀)
    /// with each validity enumerated, its degree and composition; the death's uniform pick; the change's mutations, each draw's
    /// probability. The production code gives the η (<c>make_weights</c>, checked against <c>khaos</c> in the oracle tests) and the
    /// ratios.</summary>
    internal static IEnumerable<Enumerated> Khaos(Space space, string key, State s)
    {
        var proposal = (LegendreTermProposal)space.Model.Kernel;
        LegendreBasis settings = proposal.Settings;
        LegendreTermSpace termSpace = proposal.Space;
        int p = termSpace.Variables, top = termSpace.MaxOrder, m = s.Terms.Length;
        int[] a = s.Counts.AppearanceCounts(BassPart.Legendre);
        switch (key)
        {
            case "birth":
                double[] logJ = [.. Enumerable.Range(1, top).Select(q0 => -settings.OrderPenalty * Math.Log(q0))];
                double h = logJ.Sum(Math.Exp);
                var etas = new double[top][];
                var valid = new double[top];
                for (int q0 = 1; q0 <= top; q0++)
                {
                    etas[q0 - 1] = new double[p];
                    proposal.InclusionProbabilities(a, q0, etas[q0 - 1]);
                    valid[q0 - 1] = Validity(etas[q0 - 1], top);
                }
                foreach (BassTerm t in space.Terms)
                {
                    bool[] chi = new bool[p];
                    foreach (LegendreFactor f in t.Legendre) chi[f.Variable] = true;
                    double inclusion = 0.0;
                    for (int q0 = 1; q0 <= top; q0++)
                    {
                        double product = 1.0;
                        for (int j = 0; j < p; j++) product *= chi[j] ? etas[q0 - 1][j] : 1.0 - etas[q0 - 1][j];
                        inclusion += Math.Exp(logJ[q0 - 1]) / h * product / valid[q0 - 1];
                    }
                    int q = t.Legendre.Count, d = t.LegendreDegree;
                    double logq = Math.Log(inclusion) + LogDegree(settings, q, d) + LegendreTermSpace.LogComposition(d, q);
                    yield return new Enumerated(BirthDeath, [.. s.Terms, t], logq, proposal.LogBirthRatio(t, s.Counts));
                }
                break;
            case "death":
                for (int k = 0; k < m; k++)
                    yield return new Enumerated(BirthDeath, [.. s.Terms.Where((_, i) => i != k)], -Math.Log(m), proposal.LogDeathRatio(s.Terms[k], s.Counts));
                break;
            case "change":
                double mutation = p <= 3 ? 1.0 : settings.DegreeMutation;
                for (int k = 0; k < m; k++)
                {
                    BassTerm t = s.Terms[k];
                    int q = t.Legendre.Count, d = t.LegendreDegree;
                    BassTerm[] Replaced(BassTerm next) => [.. s.Terms.Select((u, i) => i == k ? next : u)];
                    if (mutation > 0.0)
                        for (int next = q; next <= termSpace.MaxDegree; next++)
                            foreach (int[] parts in Compositions(next, q))
                            {
                                var term = new BassTerm(legendre: t.Legendre.Select((f, i) => new LegendreFactor(f.Variable, parts[i])));
                                double logq = -Math.Log(m) + Math.Log(mutation) + LogDegree(settings, q, next) + LegendreTermSpace.LogComposition(next, q);
                                yield return new Enumerated(DegreeGroup, Replaced(term), logq, proposal.LogDegreeRatio(q, d, next));
                            }
                    if (mutation >= 1.0 || q == p) continue;   // no swap, or a null one
                    bool[] held = new bool[p];
                    foreach (LegendreFactor f in t.Legendre) held[f.Variable] = true;
                    double total = Enumerable.Range(0, p).Where(u => !held[u]).Sum(u => settings.SwapPseudoCount + a[u]);
                    for (int pick = 0; pick < q; pick++)
                        for (int u = 0; u < p; u++)
                        {
                            if (held[u]) continue;
                            int from = t.Legendre[pick].Variable;
                            var term = new BassTerm(legendre: t.Legendre.Select((f, i) => i == pick ? new LegendreFactor(u, f.Degree) : f));
                            double logq = -Math.Log(m) + Math.Log(1.0 - mutation) - Math.Log(q) + Math.Log((settings.SwapPseudoCount + a[u]) / total);
                            yield return new Enumerated(SwapGroup, Replaced(term), logq, proposal.LogSwapRatio(t, from, u, a));
                        }
                }
                break;
        }
    }

    private static IEnumerable<int[]> Compositions(int d, int q)
    {
        if (q == 1)
        {
            yield return [d];
            yield break;
        }
        for (int a = 1; a <= d - q + 1; a++)
            foreach (int[] rest in Compositions(d - a, q - 1)) yield return [a, .. rest];
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Kernel_IsInDetailedBalance_WithTheEnumeratedLaw(bool raw)
    {
        BassModel model = BassModel.Create(Config(raw), Data());
        Assert.IsType<LegendreTermProposal>(model.Kernel);
        Space space = Enumerate(model);
        Assert.Equal(30, space.Terms.Count);
        Assert.Equal(1 + 30 + 30 * 29 / 2, space.States.Count);   // no structure of distinct Legendre terms is refused here
        Dictionary<string, Dictionary<(int From, int To), double>> flows = Kernel(space, (key, s) => Khaos(space, key, s));
        (Dictionary<string, int> transitions, double worst) = Check(space, flows);
        foreach (string group in new[] { BirthDeath, DegreeGroup, SwapGroup }) Assert.True(transitions.GetValueOrDefault(group) > 0, $"no {group} transitions");
        output.WriteLine($"{(raw ? "raw" : "normalized")} ridge: {space.States.Count} structures; transitions {string.Join(", ", transitions.Select(kv => $"{kv.Key} {kv.Value}"))}; worst log gap {worst:E2}");

        // The chain the enumeration builds, its spectral gap a reading.
        var chain = ReversibleChain.Of([.. space.States.Select(st => st.LogPi)], flows.Values);
        output.WriteLine($"  the sweep's chain over {chain.Size} states and {chain.Transitions} transitions");
    }

    // ln p(y | S) up to a constant common to every S, from the columns alone: (β₀, β) ~ N(0, σ²V) on the design B, σ² Jeffreys or
    // IG(a, b), so y ~ N(0, σ²C) with C = I + BVBᵀ, and the marginal is −½ln|C| − (a + n/2)ln(b + ½yᵀC⁻¹y).
    private static double DenseLogMarginal(double[][] columns, double[] y, double variance, double a, double b)
    {
        int n = y.Length;
        var c = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = i == j ? 1.0 : 0.0;
                foreach (double[] col in columns) sum += variance * col[i] * col[j];
                c[i, j] = c[j, i] = sum;
            }
        // Cholesky, then the half log determinant and the quadratic form.
        var l = new double[n, n];
        double halfLogDet = 0.0;
        for (int j = 0; j < n; j++)
        {
            double d = c[j, j];
            for (int k = 0; k < j; k++) d -= l[j, k] * l[j, k];
            l[j, j] = Math.Sqrt(d);
            halfLogDet += Math.Log(l[j, j]);
            for (int i = j + 1; i < n; i++)
            {
                double s = c[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = s / l[j, j];
            }
        }
        var z = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = y[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * z[k];
            z[i] = s / l[i, i];
        }
        double quad = z.Sum(v => v * v);
        return -halfLogDet - (a + 0.5 * n) * Math.Log(b + 0.5 * quad);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Pi_IsTheDenseMarginalsLaw(bool raw)
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(raw), data);
        Space space = Enumerate(model);
        double[] y = [.. ((ScalarResponse)data.Response).Values];
        int n = y.Length;
        double tau = model.Target.HeldTau!.Value;
        var dense = new double[space.States.Count];
        for (int k = 0; k < dense.Length; k++)
        {
            BassTerm[] terms = space.States[k].Terms;
            var columns = new List<double[]> { Enumerable.Repeat(1.0, n).ToArray() };   // the intercept, inside either ridge under the zero centre
            foreach (BassTerm t in terms)
            {
                var col = new double[n];
                for (int i = 0; i < n; i++)
                    col[i] = t.Legendre.Aggregate(1.0, (acc, f) => acc * ShiftedLegendre.Value(f.Degree, model.Inputs.Cases.Continuous(f.Variable)[i]));
                if (!raw)
                {
                    // The normalized ridge's: centred, then scaled to unit root mean square over the n cases.
                    double mean = col.Average();
                    for (int i = 0; i < n; i++) col[i] -= mean;
                    double rms = Math.Sqrt(col.Sum(v => v * v) / n);
                    for (int i = 0; i < n; i++) col[i] /= rms;
                }
                columns.Add(col);
            }
            dense[k] = model.Prior.LogCollectionWeight(terms.Length) + terms.Sum(model.Prior.LogPrior) + DenseLogMarginal([.. columns], y, 1.0 / tau, 0.0, 0.0);
        }
        double[] engine = space.Posterior();
        double top = dense.Max();
        double[] expected = [.. dense.Select(v => Math.Exp(v - top))];
        double total = expected.Sum();
        double worst = 0.0;
        for (int k = 0; k < engine.Length; k++) worst = Math.Max(worst, Math.Abs(engine[k] - expected[k] / total));
        output.WriteLine($"{(raw ? "raw" : "normalized")} ridge: π against the dense marginal's law, worst |Δ| {worst:E2}; P(M = 0, 1, 2) {string.Join(" / ", Enumerable.Range(0, 3).Select(mm => space.States.Where((st, i) => st.Terms.Length == mm).Sum(st => engine[space.States.IndexOf(st)]).ToString("F4")))}");
        Assert.True(worst < 1e-11, $"worst |Δπ| {worst:E2}");
    }

    // The fixture's inputs drawn into [0.05, 0.95] with [0, 1] declared: the engine's π is the dense marginal's on ψ at the raw values,
    // where the training range would stretch them back onto [0, 1] and read another law.
    [Fact]
    public void Pi_ReadsTheDeclaredDomain()
    {
        BassData unit = Data();
        BassData Shrunk(InputDomain? domain)
            => new([.. unit.Predictors.Cast<ContinuousPredictor>().Select(c => new ContinuousPredictor(c.Name, [.. c.Values.ToArray().Select(x => 0.05 + 0.9 * x)], domain))],
                   unit.Response);
        BassData data = Shrunk(new InputDomain(0.0, 1.0));
        BassModel model = BassModel.Create(Config(raw: true), data);
        Space space = Enumerate(model);
        double[] y = [.. ((ScalarResponse)data.Response).Values];
        int n = y.Length;
        double tau = model.Target.HeldTau!.Value;
        double[][] x = [.. data.Predictors.Cast<ContinuousPredictor>().Select(c => c.Values.ToArray())];
        var dense = new double[space.States.Count];
        for (int k = 0; k < dense.Length; k++)
        {
            BassTerm[] terms = space.States[k].Terms;
            var columns = new List<double[]> { Enumerable.Repeat(1.0, n).ToArray() };
            foreach (BassTerm t in terms)
                columns.Add([.. Enumerable.Range(0, n).Select(i => t.Legendre.Aggregate(1.0, (acc, f) => acc * ShiftedLegendre.Value(f.Degree, x[f.Variable][i])))]);
            dense[k] = model.Prior.LogCollectionWeight(terms.Length) + terms.Sum(model.Prior.LogPrior) + DenseLogMarginal([.. columns], y, 1.0 / tau, 0.0, 0.0);
        }
        double[] engine = space.Posterior();
        double top = dense.Max();
        double[] expected = [.. dense.Select(v => Math.Exp(v - top))];
        double total = expected.Sum();
        double worst = 0.0;
        for (int k = 0; k < engine.Length; k++) worst = Math.Max(worst, Math.Abs(engine[k] - expected[k] / total));

        // Undeclared, the same values are stretched onto [0, 1], and the law moves.
        Space stretched = Enumerate(BassModel.Create(Config(raw: true), Shrunk(null)));
        double[] other = stretched.Posterior();
        double moved = 0.0;
        for (int k = 0; k < engine.Length; k++)
            if (stretched.Find(space.States[k].Terms) is int at and >= 0) moved = Math.Max(moved, Math.Abs(engine[k] - other[at]));
        output.WriteLine($"inputs in [0.05, 0.95], [0, 1] declared: π against the dense marginal on the raw values, worst |Δ| {worst:E2}; " +
                         $"the training range's law differs by up to {moved:F4}");
        Assert.True(worst < 1e-11, $"worst |Δπ| {worst:E2}");
        Assert.True(moved > 1e-3, $"the domain moved the law by only {moved:E2}");
    }

    // Batch means' standard error of `values`' mean over `batches` batches.
    internal static double BatchError(IReadOnlyList<double> values, int batches)
    {
        int size = values.Count / batches;
        double[] means = [.. Enumerable.Range(0, batches).Select(k => Enumerable.Range(k * size, size).Average(i => values[i]))];
        double grand = means.Average();
        return Math.Sqrt(means.Sum(v => (v - grand) * (v - grand)) / (batches - 1) / batches);
    }

    [Theory]
    [InlineData(true, 20261110)]
    [InlineData(false, 20261111)]
    public void Chain_MeetsTheExactLaw(bool raw, int seed)
    {
        BassModel model = BassModel.Create(Config(raw), Data());
        Space space = Enumerate(model);
        double[] pi = space.Posterior();
        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
        const int burn = 5_000, steps = 400_000, batches = 40;
        for (int t = 0; t < burn; t++) chain.Step();
        var visits = new int[steps];
        for (int t = 0; t < steps; t++)
        {
            BassState s = chain.Step();
            int i = space.Find(s.Terms);
            Assert.True(i >= 0, $"the chain holds a structure outside the enumeration: {string.Join(" + ", s.Terms)}");
            visits[t] = i;
        }
        // Readings: P(M = m), each variable's inclusion, and the six likeliest structures.
        var readings = new List<(string Name, Func<State, bool> Holds)>();
        for (int mm = 0; mm <= 2; mm++) { int m = mm; readings.Add(($"P(M = {m})", st => st.Terms.Length == m)); }
        for (int v = 0; v < P; v++) { int j = v; readings.Add(($"P(x{j} in)", st => st.Terms.Any(t => t.Legendre.Any(f => f.Variable == j)))); }
        foreach (int top in Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(6))
        {
            int k = top;
            readings.Add(($"P({(space.States[k].Terms.Length == 0 ? "∅" : string.Join(" + ", space.States[k].Terms))})", st => ReferenceEquals(st, space.States[k])));
        }
        double worst = 0.0;
        foreach ((string name, Func<State, bool> holds) in readings)
        {
            bool[] indicator = [.. space.States.Select(holds)];
            double exact = Enumerable.Range(0, pi.Length).Where(i => indicator[i]).Sum(i => pi[i]);
            double[] series = [.. visits.Select(i => indicator[i] ? 1.0 : 0.0)];
            double mean = series.Average(), se = BatchError(series, batches);
            double z = se > 0.0 ? (mean - exact) / se : (mean == exact ? 0.0 : double.PositiveInfinity);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{name,-34} chain {mean:F4}  exact {exact:F4}  z {z,6:F2}");
        }
        output.WriteLine($"{(raw ? "raw" : "normalized")} ridge, {steps} steps: worst |z| {worst:F2}; moves {string.Join(", ", chain.Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}");
        Assert.True(worst < 4.5, $"the chain's worst |z| {worst:F2} against the exact law");
    }

    // ---- Geweke --------------------------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;

    // Proper priors that read no response: the zero centre, an inverse gamma on σ², τ held; p = 3 mutates degrees only, so a
    // fourth predictor brings the swap in.
    private static BassConfig GewekeConfig(bool raw) => Config(raw, maxBasis: 3, sigma2: Sigma2Prior.InverseGamma(3, 2)) with
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 3, Legendre = new LegendreBasis { MaxDegree = 3, DegreeMutation = 0.5 } },
    };

    private static BassData GewekeData4(double[] y)
    {
        BassData d = KhaosBasisTests.UnitData(GewekeN, 4, 77);
        return new BassData(d.Predictors, new ScalarResponse(y));
    }

    private static readonly string[] GewekeNames = ["M", "ln sigma2", "atan f0", "atan f6", "x0 in", "x3 in", "degree sum", "z0^2"];

    private static double[] Functions(BassState state, double sigma2, double[] f, double[] y)
        => [state.Count, Math.Log(sigma2), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
            state.Terms.Any(t => t.Legendre.Any(x => x.Variable == 0)) ? 1 : 0, state.Terms.Any(t => t.Legendre.Any(x => x.Variable == 3)) ? 1 : 0,
            state.Terms.Sum(t => t.LegendreDegree), (y[0] - f[0]) * (y[0] - f[0]) / sigma2];

    [Theory]
    [InlineData(true, 20261112)]
    [InlineData(true, 20261114)]
    [InlineData(false, 20261113)]
    public void Geweke_TheSweepLeavesTheJointInvariant(bool raw, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 200_000, steps = 300_000, batches = 50;
        BassConfig config = GewekeConfig(raw);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData4([.. Enumerable.Range(0, GewekeN).Select(i => Math.Sin(i))]));
        var sampler = new BassPriorSampler(prior);
        double logTau = Math.Log(prior.Target.HeldTau!.Value);

        (BassState S, double Sigma2, double[] F, double[] Y) Draw()
        {
            BassState state = sampler.DrawState(rng).State;
            var chain = new BassChain(prior, state.Terms, logTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
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
            BassModel model = BassModel.Create(config, GewekeData4(y));
            var chain = new BassChain(model, theta.S.Terms, logTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
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
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"{(raw ? "raw" : "normalized")} ridge, seed {seed}: worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- the config ----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Config_RefusesWhatALegendreBasisDoesNotRead()
    {
        BassConfig good = Config(raw: true);
        good.Validate();
        BassStructureSpec s = good.Structure;
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Degree = 2 } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Support = new SupportRule(SupportKind.ColumnMass) } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Legendre = new LegendreBasis { InclusionPasses = 0 } } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Sampling = good.Sampling with { LocalSteps = new LocalSteps() } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Sampling = good.Sampling with { FactorMoves = new FactorMoves() } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held()) }).Validate());

        // Data the basis does not serve, and a d_max below q_max.
        BassData categorical = new([new ContinuousPredictor("x", [0, 0.5, 1, 0.2]), new CategoricalPredictor("z", [0, 1, 0, 1], ["a", "b"])],
                                   new ScalarResponse([1, 2, 3, 2.5]));
        Assert.Throws<ArgumentException>(() => BassModel.Create(good, categorical));
        Assert.Throws<ArgumentException>(() => BassModel.Create(good with { Structure = s with { Legendre = new LegendreBasis { MaxDegree = 1 } } }, Data()));

        BassModel model = BassModel.Create(good, Data());
        Assert.Throws<InvalidOperationException>(() => model.Proposal);
        Assert.True(model.Prior.Has(BassPart.Legendre));
        Assert.Equal(30.0, model.Prior.Legendre!.Size);

        // A Legendre column sits near unit root mean square, so the raw ridge's held τ halves its shrinkage where the normalized ridge's does.
        ResolvedTauPrior held = TauPrior.Held(0.002).For(100.0);
        Assert.Equal(0.002 / 400.0, BassResolved.HalvingTau(new RawRidgeForm(), held, 100.0, 4.0));
    }

    [Fact]
    public void Config_SerializesTheBasisAndTheRawRidge_AndAHingeConfigAsItWas()
    {
        Assert.DoesNotContain("Legendre", new BassConfig().ToJson());
        Assert.DoesNotContain("RawRidge", new BassConfig().ToJson());
        BassConfig config = Config(raw: true) with { Structure = Config(raw: true).Structure with { Legendre = new LegendreBasis { MaxDegree = 9, DegreePenalty = 0.4 } } };
        string json = config.ToJson();
        Assert.Contains("\"RawRidge\"", json);
        BassConfig back = BassConfig.FromJson(json);
        Assert.Equal(config.Structure.Legendre, back.Structure.Legendre);
        Assert.IsType<RawRidgeForm>(back.Coefficients.Form);
        Assert.Equal(json, back.ToJson());
    }
}
