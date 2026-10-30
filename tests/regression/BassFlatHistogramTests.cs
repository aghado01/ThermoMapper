using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's flat-histogram sampler (N17 Phase 2; the Phase 2 brief): the configuration and its refusals, and BASS's chain
/// on the multicanonical target π_w ∝ prior(structure, τ)·w(ℓ₁) against π_w by enumeration on fixture A's F, with τ held
/// and with τ sampled, by quadrature over u.
/// </summary>
public sealed class BassFlatHistogramTests
{
    private readonly ITestOutputHelper _out;
    public BassFlatHistogramTests(ITestOutputHelper output) => _out = output;

    // A small run's spec and its pilot.
    internal static BassFlatHistogramSpec Spec(FlatHistogramSpec? run = null, BassFlatHistogramPilot? pilot = null, long? maxIterations = null)
        => new()
        {
            Run = run ?? new FlatHistogramSpec { LearningIterations = 200, ProductionDraws = 100, Start = FlatHistogramStart.Wham },
            Pilot = pilot ?? new BassFlatHistogramPilot { BurnIn = 63, TracedScans = 40, Chains = 2 },
            MaxIterations = maxIterations,
        };

    internal static BassConfig WithFlatHistogram(BassConfig config, BassFlatHistogramSpec spec)
        => config with { Sampling = config.Sampling with { FlatHistogram = spec } };

    // ---- Step 1: the configuration ----

    [Fact]
    public void Config_RoundTripsThroughJson_AndMovesTheFingerprint()
    {
        BassData data = BassSweepTests.Data();
        BassConfig plain = BassSweepTests.Config();
        BassFlatHistogramSpec explicitSpec = Spec(new FlatHistogramSpec
        {
            LearningIterations = 1_000, ProductionDraws = 500, Chains = 3, Window = new FlatHistogramWindow(-40.0, -10.0), LowSlope = 0.0, Nodes = 12,
            GainScale = 50.0, GainStart = 7,
        }) with { Pilot = null, MaxIterations = 300 };
        foreach (BassFlatHistogramSpec spec in new[] { Spec(), explicitSpec })
        {
            BassConfig config = WithFlatHistogram(plain, spec);
            config.Validate();
            string json = config.ToJson();
            _out.WriteLine(json);
            Assert.Equal(config, BassConfig.FromJson(json));
            Assert.Equal(json, BassConfig.FromJson(json).ToJson());
            Assert.DoesNotContain("NeedsPilot", json);
            Assert.NotEqual(Bass.Fingerprint(plain, data), Bass.Fingerprint(config, data));
        }
        Assert.Contains("\"FlatHistogram\":null", new BassConfig().ToJson());
        Assert.Contains("\"Pilot\":{\"Tempering\":{\"Path\":\"Likelihood\",\"Form\":\"Linear\"", WithFlatHistogram(plain, Spec()).ToJson());
        // A field the spec does not have is refused, and a required one missing.
        string good = WithFlatHistogram(plain, Spec()).ToJson();
        Assert.ThrowsAny<Exception>(() => BassConfig.FromJson(good.Replace("\"TracedScans\":40", "\"TracedScans\":40,\"Scans\":4")));
        Assert.ThrowsAny<Exception>(() => BassConfig.FromJson(good.Replace(",\"TracedScans\":40", "")));
    }

    [Fact]
    public void Config_RefusesWhatTheSamplerCannotRun()
    {
        BassConfig plain = BassSweepTests.Config();
        void Refused(BassConfig c, string field) => Assert.Contains(field, Assert.Throws<ArgumentException>(c.Validate).Message);
        BassConfig With(BassFlatHistogramSpec spec, BassConfig? config = null) => WithFlatHistogram(config ?? plain, spec);

        // Another sampler beside it.
        Refused(With(Spec()) with { Sampling = With(Spec()).Sampling with { Smc = new SmcSpec() } }, "Sampling.FlatHistogram");
        Refused(With(Spec()) with { Sampling = With(Spec()).Sampling with { Tempering = new TemperingSpec() } }, "Sampling.Tempering");
        Refused(With(Spec()) with { Sampling = With(Spec()).Sampling with { RHatTarget = 1.01 } }, "Sampling.RHatTarget");

        // The pilot: the integrated form, the powered path, a missing pilot, one nothing reads, and its own fields.
        TemperingSpec linear = new BassFlatHistogramPilot { BurnIn = 1, TracedScans = 2 }.Tempering;
        BassFlatHistogramPilot Pilot(TemperingSpec t) => new() { BurnIn = 63, TracedScans = 40, Tempering = t };
        Refused(With(Spec(pilot: Pilot(linear with { Form = TemperedForm.Integrated }))), "Sampling.FlatHistogram.Pilot.Tempering.Form");
        Refused(With(Spec(pilot: Pilot(new TemperingSpec { Path = TemperingPath.Powered }))), "Sampling.FlatHistogram.Pilot.Tempering.Path");
        Refused(With(Spec(pilot: Pilot(linear with { Levels = 2 }))), "Sampling.FlatHistogram.Pilot.Tempering.Levels");
        Refused(With(Spec() with { Pilot = null }), "Sampling.FlatHistogram.Pilot");
        FlatHistogramSpec declared = new()
        {
            LearningIterations = 10, ProductionDraws = 10, Window = new FlatHistogramWindow(-40.0, -10.0), LowSlope = 0.2, Nodes = 8,
        };
        Refused(With(Spec(declared)), "Sampling.FlatHistogram.Pilot");
        With(Spec(declared) with { Pilot = null }).Validate();
        Refused(With(Spec(pilot: new BassFlatHistogramPilot { BurnIn = 63, TracedScans = 1 })), "Sampling.FlatHistogram.Pilot.TracedScans");
        Refused(With(Spec(pilot: new BassFlatHistogramPilot { BurnIn = -1, TracedScans = 40 })), "Sampling.FlatHistogram.Pilot.BurnIn");
        Refused(With(Spec(pilot: new BassFlatHistogramPilot { BurnIn = 63, TracedScans = 40, Chains = 0 })), "Sampling.FlatHistogram.Pilot.Chains");

        // The engine's own fields, named under the spec; the call budget.
        Refused(With(Spec(new FlatHistogramSpec { LearningIterations = 10, ProductionDraws = 10, Chains = 0, Start = FlatHistogramStart.Wham })),
                "Sampling.FlatHistogram.Run.Chains");
        Refused(With(Spec(new FlatHistogramSpec { LearningIterations = 10, ProductionDraws = 0, Start = FlatHistogramStart.Wham })),
                "Sampling.FlatHistogram.Run.ProductionDraws");
        Refused(With(Spec(new FlatHistogramSpec { LearningIterations = 10, ProductionDraws = 10, LowSlope = 1.5, Start = FlatHistogramStart.Wham })),
                "Sampling.FlatHistogram.Run.LowSlope");
        Refused(With(Spec(maxIterations: 0)), "Sampling.FlatHistogram.MaxIterations");

        // An improper τ prior wherever b_lo = 0: the pilot's reference rung, or a declared b_lo of 0.
        BassConfig improper = plain with { Coefficients = BassCoefficientPrior.G(TauPrior.HyperG(2.0)) };
        Refused(With(Spec(), improper), "Coefficients.Scale");
        Refused(With(Spec(declared with { LowSlope = 0.0 }) with { Pilot = null }, improper), "Coefficients.Scale");
        // From b_lo > 0 it runs: a pilot without the reference rung, or a declared b_lo above 0.
        With(Spec(pilot: Pilot(linear with { ReferenceRung = false })), improper).Validate();
        With(Spec(declared) with { Pilot = null }, improper).Validate();

        // Jeffreys σ² and the flat intercept keep the prior over (structure, τ) proper, so b_lo = 0 runs (with no evidence).
        With(Spec(), plain with { Sigma2 = Sigma2Prior.Jeffreys }).Validate();
        With(Spec(), plain with { Intercept = GPriorCentre.Flat }).Validate();
    }

    // ---- Step 1: the chain on π_w ----

    // Fixture A (BassSweepTests' F) with τ sampled, or held at unit information.
    internal static BassModel FixtureA(bool heldTau)
    {
        BassConfig config = BassSweepTests.Config();
        if (heldTau) config = config with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) };
        return BassModel.Create(config, BassSweepTests.Data());
    }

    // Weights over the middle of fixture A's ℓ₁, at SAMC's limit on its exact density of states, so π_w spreads over F.
    internal static MulticanonicalWeights SpreadWeights(BassFlatHistogramSpace space, int nodes = 16, double lowSlope = 0.2)
    {
        (double[] levels, double[] masses) = space.DensityOfStates();
        double low = BassFlatHistogramSpace.Quantile(levels, masses, 0.02), high = levels.Max() - 0.5;
        return SamcLearning.FixedPoint(MulticanonicalWeights.Flat(low, high, lowSlope, nodes), levels, masses, out _);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Chain_SamplesPiW_OnFixtureA(bool heldTau)
    {
        BassModel model = FixtureA(heldTau);
        var space = BassFlatHistogramSpace.Of(model);
        MulticanonicalWeights weights = SpreadWeights(space);
        double[] law = space.Law(weights.LogWeight);

        int cap = model.MaxBasis;
        var countLaw = new double[cap + 1];
        var inclusion = new Dictionary<BassTerm, double>();
        double meanU = 0.0, meanLevel = 0.0;
        for (int s = 0; s < space.States.Count; s++)
        {
            BassFlatHistogramSpace.Structure st = space.States[s];
            countLaw[st.Terms.Length] += law[s];
            foreach (BassTerm t in st.Terms) inclusion[t] = inclusion.GetValueOrDefault(t) + law[s];
            meanU += law[s] * space.MeanU(st, weights.LogWeight);
            meanLevel += law[s] * space.MeanLevel(st, weights.LogWeight);
        }
        BassTerm[] tracked = [.. inclusion.Where(x => x.Value >= 0.005).Select(x => x.Key)];
        int[] top = [.. Enumerable.Range(0, law.Length).OrderByDescending(s => law[s]).Take(10)];
        var topIndex = new Dictionary<string, int>();
        for (int k = 0; k < top.Length; k++) topIndex.Add(space.States[top[k]].Key, k);
        var termIndex = new Dictionary<BassTerm, int>();
        for (int k = 0; k < tracked.Length; k++) termIndex.Add(tracked[k], k);

        const int chains = 4, burn = 1_000, sweeps = 30_000, batches = 25;
        int items = cap + 1 + tracked.Length + top.Length + 2;
        var batch = new Batches(items, chains, batches, sweeps);
        var target = new MulticanonicalTarget<BassState>(model.Target, weights);
        var client = new BassMulticanonicalModel(model, "fixture A", learningIterations: 0, new BassSweepCount(new BassRunOptions(), 0));
        int[] seeds = SeedTree.Derive(7101, 2 * chains);
        Parallel.For(0, chains, c =>
        {
            var chain = (BassMulticanonicalChain)client.StartChain(target, seeds[c], seeds[chains + c]);
            for (int s = 0; s < burn + sweeps; s++)
            {
                BassState state = chain.Step();
                MulticanonicalEvaluation e = chain.CurrentEvaluation;
                if (s % 97 == 0)
                {
                    // The chain's density is the inner prior plus ln w(ℓ₁), and its evaluation is its state's.
                    Assert.Equal(e.LogPrior + weights.LogWeight(e.UntemperedLogLikelihood), chain.Chain.Moves.CurrentLogDensity, 12);
                    Assert.Equal(model.Target.Evaluate(state).LogMarginal, e.UntemperedLogLikelihood);
                }
                if (s < burn) continue;
                int at = s - burn;
                batch.Add(state.Count, c, at, 1.0);
                foreach (BassTerm t in state.Terms)
                    if (termIndex.TryGetValue(t, out int k)) batch.Add(cap + 1 + k, c, at, 1.0);
                if (topIndex.TryGetValue(BassMixing.Key(state.Terms), out int j)) batch.Add(cap + 1 + tracked.Length + j, c, at, 1.0);
                batch.Add(items - 2, c, at, state.LogTau);
                batch.Add(items - 1, c, at, e.UntemperedLogLikelihood);
            }
        });

        _out.WriteLine($"## π_w on fixture A, τ {(heldTau ? "held" : "sampled")}: {space.States.Count:N0} structures in F; window " +
                       $"[{weights.Low:F2}, {weights.High:F2}], K = {weights.Nodes}, b_lo = {weights.LowSlope}; {chains} × {sweeps:N0} sweeps");
        _out.WriteLine("quantity | π_w | chain ± batch SE | z | allowance");
        _out.WriteLine("--- | --- | --- | --- | ---");
        void Check(string name, double expected, int item, double slack)
        {
            (double mean, double se) = batch.Read(item);
            double allowance = 5 * se + slack;
            _out.WriteLine($"{name} | {expected:F4} | {mean:F4} ± {se:F4} | {(se > 0 ? $"{(mean - expected) / se:F2}" : "—")} | {allowance:F4}");
            Assert.True(Math.Abs(mean - expected) <= allowance, $"{name}: the chain reads {mean:F4} ± {se:F4}, π_w {expected:F4}");
        }
        for (int m = 0; m <= cap; m++) Check($"P(M = {m})", countLaw[m], m, 0.004);
        for (int k = 0; k < tracked.Length; k++) Check($"incl. {BassMixing.Key([tracked[k]])}", inclusion[tracked[k]], cap + 1 + k, 0.005);
        for (int k = 0; k < top.Length; k++) Check($"S = {space.States[top[k]].Key}", law[top[k]], cap + 1 + tracked.Length + k, 0.005);
        if (!heldTau) Check("E[u]", meanU, items - 2, 0.001);
        Check("E[ℓ₁]", meanLevel, items - 1, 0.01);
    }


    [Fact]
    public void Reweigh_ReadsTheCachedLevel_WithNoNewMarginal()
    {
        BassModel model = FixtureA(heldTau: false);
        var counting = new Counting(model.Target);
        MulticanonicalWeights before = MulticanonicalWeights.Flat(-60.0, -10.0, 0.3, 11);
        MulticanonicalWeights after = before.WithLogWeights([.. Enumerable.Range(0, 11).Select(k => 0.4 * k - 0.05 * k * k)]);
        var target = new MulticanonicalTarget<BassState>(counting, before);
        var moved = new MulticanonicalTarget<BassState>(counting, after);
        var client = new BassMulticanonicalModel(model, "fixture A", learningIterations: 0, new BassSweepCount(new BassRunOptions(), 0));
        var started = (BassMulticanonicalChain)client.StartChain(target, 11, 12);
        for (int s = 0; s < 200; s++) started.Step();

        // The same chain, its target's weights moved: Reweigh re-reads ln w from the cached ℓ₁.
        var chain = new BassMulticanonicalChain(client, started.Chain, moved, started.Composition, started.Terms, started.Draws, started.Start, null);
        MulticanonicalEvaluation old = chain.CurrentEvaluation;
        long evaluations = counting.Evaluations;
        chain.Reweigh();
        MulticanonicalEvaluation now = chain.CurrentEvaluation;
        Assert.Equal(evaluations, counting.Evaluations);
        Assert.Same(old.Inner, now.Inner);
        Assert.Same(after, now.Weights);
        Assert.Equal(old.UntemperedLogLikelihood, now.UntemperedLogLikelihood);
        Assert.Equal(after.LogWeight(old.UntemperedLogLikelihood), now.LogWeight);
        Assert.Equal(now.LogPrior + now.LogWeight, chain.Chain.Moves.CurrentLogDensity);
        Assert.NotEqual(old.LogWeight, now.LogWeight);
        // Under the same weights, the evaluation is kept.
        chain.Reweigh();
        Assert.Same(now, chain.CurrentEvaluation);
        Assert.Equal(evaluations, counting.Evaluations);
    }

    // The collapsed target, its evaluations counted.
    private sealed class Counting(BassTarget inner) : ISamplingTarget<BassState>
    {
        private long _evaluations;
        public long Evaluations => System.Threading.Interlocked.Read(ref _evaluations);

        public ITargetEvaluation Evaluate(BassState state)
        {
            System.Threading.Interlocked.Increment(ref _evaluations);
            return inner.Evaluate(state);
        }
    }

    /// <summary>Batch means over chains: each item's per-batch sums, read as the mean and its standard error.</summary>
    internal sealed class Batches(int items, int chains, int batches, int length)
    {
        private readonly double[,,] _sums = new double[items, chains, batches];

        public void Add(int item, int chain, int step, double value) => _sums[item, chain, step / (length / batches)] += value;

        public (double Mean, double Error) Read(int item)
        {
            double size = length / batches;
            double[] means = [.. Enumerable.Range(0, chains).SelectMany(c => Enumerable.Range(0, batches).Select(b => _sums[item, c, b] / size))];
            double mean = means.Average();
            return (mean, Math.Sqrt(means.Sum(m => (m - mean) * (m - mean)) / (means.Length - 1) / means.Length));
        }
    }
}

/// <summary>
/// F enumerated on a small fixture (every structure the design and the target admit, extended from admitted ones only, F
/// being downward closed), each with ln w_S, its prior weight on F up to one constant, and ℓ₁ at each node of a trapezoid
/// grid in u where p(u) holds all but e⁻⁶⁰ of its peak; one node at the held τ. From it, the exact density of states, the
/// laws under any weight of ℓ₁, and ln Z(b).
/// </summary>
internal sealed class BassFlatHistogramSpace
{
    /// <summary>The u grid's nodes.</summary>
    public const int Nodes = 801;

    private BassFlatHistogramSpace(BassModel model, List<Structure> states, double[] u, double[] logQuadrature)
    {
        Model = model;
        States = states;
        U = u;
        LogQuadrature = logQuadrature;
    }

    /// <summary>A structure in F: its terms as loaded, its key, ln w_S and ℓ₁ at each node.</summary>
    public sealed record Structure(BassTerm[] Terms, string Key, double LogW, double[] Levels);

    public BassModel Model { get; }
    public List<Structure> States { get; }

    /// <summary>The u nodes; the held τ's log alone when τ is held.</summary>
    public double[] U { get; }

    /// <summary>Each node's ln of p(u)·du with the trapezoid's end halves; 0 at the held τ.</summary>
    public double[] LogQuadrature { get; }

    public static BassFlatHistogramSpace Of(BassModel model)
    {
        BassTarget target = model.Target;
        double[] u, logQuadrature;
        if (target.SamplesTau)
        {
            double[] coarse = [.. Enumerable.Range(0, 4001).Select(k => -50.0 + 0.025 * k)];
            double[] lp = [.. coarse.Select(target.LogTauDensity)];
            double peak = lp.Max();
            int lo = Array.FindIndex(lp, v => v >= peak - 60), hi = Array.FindLastIndex(lp, v => v >= peak - 60);
            double a = coarse[Math.Max(0, lo - 1)], z = coarse[Math.Min(coarse.Length - 1, hi + 1)];
            u = [.. Enumerable.Range(0, Nodes).Select(k => a + (z - a) * k / (Nodes - 1))];
            double logDu = Math.Log((z - a) / (Nodes - 1));
            logQuadrature = [.. u.Select((v, k) => target.LogTauDensity(v) + logDu + (k == 0 || k == Nodes - 1 ? -Math.Log(2.0) : 0.0))];
        }
        else
        {
            u = [Math.Log(target.HeldTau!.Value)];
            logQuadrature = [0.0];
        }

        List<BassTerm> terms = BassKernelEnumeration.AllTerms(model);
        var states = new List<Structure>();
        void Visit(int from, List<BassTerm> chosen)
        {
            if (BassMixing.Loaded(model, chosen, u[0], out _) is not { } loaded) return;   // no superset is in F
            BassTerm[] held = [.. loaded.Terms];
            var levels = new double[u.Length];
            for (int k = 0; k < u.Length; k++)
                levels[k] = target.SamplesTau ? target.Evaluate(loaded.Fit, u[k]).LogMarginal : target.Evaluate(loaded.Fit).LogMarginal;
            states.Add(new Structure(held, BassMixing.Key(held), target.LogStructurePrior(loaded.Fit), levels));
            if (chosen.Count == model.MaxBasis) return;
            for (int i = from; i < terms.Count; i++)
            {
                chosen.Add(terms[i]);
                Visit(i + 1, chosen);
                chosen.RemoveAt(chosen.Count - 1);
            }
        }
        Visit(0, []);
        return new BassFlatHistogramSpace(model, states, u, logQuadrature);
    }

    /// <summary>The prior's density of states over F × the u nodes: each (structure, node)'s ℓ₁ and its log mass
    /// ln w_S + ln q_k, unnormalized.</summary>
    public (double[] Levels, double[] LogMasses) DensityOfStates()
    {
        int n = U.Length;
        var levels = new double[States.Count * n];
        var masses = new double[levels.Length];
        for (int s = 0; s < States.Count; s++)
            for (int k = 0; k < n; k++)
            {
                levels[s * n + k] = States[s].Levels[k];
                masses[s * n + k] = States[s].LogW + LogQuadrature[k];
            }
        return (levels, masses);
    }

    /// <summary>The law over F of prior·e^{f(ℓ₁)}, normalized: f = ln w for π_w, b·ℓ for the tempered posterior.</summary>
    public double[] Law(Func<double, double> logWeight)
    {
        double[] log = [.. States.Select(s => s.LogW + LogIntegral(s, logWeight))];
        double total = LogSumExp(log);
        return [.. log.Select(v => Math.Exp(v - total))];
    }

    /// <summary>ln Σ_k q_k e^{f(ℓ₁(S, u_k))}.</summary>
    public double LogIntegral(Structure s, Func<double, double> logWeight)
    {
        var v = new double[U.Length];
        for (int k = 0; k < v.Length; k++) v[k] = LogQuadrature[k] + logWeight(s.Levels[k]);
        return LogSumExp(v);
    }

    /// <summary>E[u | S] under p(u)·e^{f(ℓ₁)}.</summary>
    public double MeanU(Structure s, Func<double, double> logWeight) => Mean(s, logWeight, k => U[k]);

    /// <summary>E[ℓ₁ | S] under p(u)·e^{f(ℓ₁)}.</summary>
    public double MeanLevel(Structure s, Func<double, double> logWeight) => Mean(s, logWeight, k => s.Levels[k]);

    /// <summary>ln Z(b) = ln E_prior[e^{b·ℓ₁}], the prior on F normalized over the grid.</summary>
    public double LogEvidence(double b)
        => LogSumExp(States.Select(s => s.LogW + LogIntegral(s, l => b * l))) - LogSumExp(States.Select(s => s.LogW + LogIntegral(s, _ => 0.0)));

    private double Mean(Structure s, Func<double, double> logWeight, Func<int, double> value)
    {
        double top = double.NegativeInfinity;
        for (int k = 0; k < U.Length; k++) top = Math.Max(top, LogQuadrature[k] + logWeight(s.Levels[k]));
        double num = 0.0, den = 0.0;
        for (int k = 0; k < U.Length; k++)
        {
            double w = Math.Exp(LogQuadrature[k] + logWeight(s.Levels[k]) - top);
            num += w * value(k);
            den += w;
        }
        return num / den;
    }

    /// <summary>The <paramref name="p"/>-quantile of a discrete measure.</summary>
    public static double Quantile(double[] levels, double[] logMasses, double p)
    {
        int[] order = [.. Enumerable.Range(0, levels.Length).OrderBy(i => levels[i])];
        double top = logMasses.Max(), total = logMasses.Sum(m => Math.Exp(m - top)), run = 0.0;
        foreach (int i in order)
        {
            run += Math.Exp(logMasses[i] - top) / total;
            if (run >= p) return levels[i];
        }
        return levels[order[^1]];
    }

    public static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        if (double.IsNegativeInfinity(max)) return max;
        return max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }
}
