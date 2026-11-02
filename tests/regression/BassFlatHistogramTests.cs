using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's flat-histogram sampler (N17 Phase 2; the Phase 2 brief).
/// <list type="bullet">
/// <item><b>Step 1.</b> The configuration and its refusals, and BASS's chain on the multicanonical target
/// π_w ∝ prior(structure, τ)·w(ℓ₁) against π_w by enumeration on fixture A's F, with τ held and with τ sampled, by
/// quadrature over u.</item>
/// <item><b>Step 2.</b> The fit on fixture A, τ held and sampled, under the g-form and the ridge: its learned weights against
/// SAMC's limit on the exact density of states, its reweighted posterior against the exact one, and its ln Z(b) against
/// the enumeration's; ln Z(1) on the oracle's fixtures B and C; every view end to end on a scalar and a functional fit;
/// resume across learning and production; the improper options; a chain that never tunnels.</item>
/// <item><b>Step 3.</b> The pilot: linear-form, reference-rung, self-sizing ladders read by
/// <see cref="TemperingPilot.FromLadders"/> and pooled; eight pilots' WHAM starts against SAMC's limit on the exact density
/// of states; learning from the warm start against learning from flat; a pilot whose ĝ has a gap, interpolated or
/// refused.</item>
/// <item><b>The warm start's chip.</b> Why SAMC's fixed point refuses a pilot's ĝ (flat occupancy outside the hull of its
/// atoms' hats); the starts it accepts pinned before the remedy; 48 pilots over their own shapes and the gap pilots, the
/// refused ones tied at ten draws, and a pilot too thin to tie refused.</item>
/// <item><b>The loose ends' chip.</b> Today's fits pinned; a window whose top the learning raises, against the
/// enumeration; a margin above the pilot's draws through the adapter; the readout's tied nodes, and a refusal naming the
/// segment where flat occupancy leaves the hull.</item>
/// </list>
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

    // Fixture A (BassSweepTests' F) with τ sampled, or held at unit information, under the g-form or the normalized ridge.
    internal static BassModel FixtureA(bool heldTau, bool ridge = false)
    {
        TauPrior scale = heldTau ? TauPrior.Held(ridge ? 1.0 : null) : TauPrior.UnitInformationGamma();
        BassConfig config = BassSweepTests.Config() with { Coefficients = ridge ? BassCoefficientPrior.Ridge(scale) : BassCoefficientPrior.G(scale) };
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

    // ---- Step 2: the fit, its evidence and its views ----

    private const int FitChains = 8;

    // A flat-start run of `chains` chains on `shape`'s declared window, b_lo and K.
    internal static FlatHistogramSpec DeclaredRun(MulticanonicalWeights shape, long learning, int production, int chains = FitChains)
        => new()
        {
            Chains = chains, LearningIterations = learning, ProductionDraws = production, Window = new FlatHistogramWindow(shape.Low, shape.High),
            LowSlope = shape.LowSlope, Nodes = shape.Nodes,
        };

    // The window over F's exact density of states: from the prior's 1% quantile of ℓ₁ to the posterior's 99%.
    internal static MulticanonicalWeights ExactShape(BassFlatHistogramSpace space, int nodes, double lowSlope = 0.0)
    {
        (double[] levels, double[] masses) = space.DensityOfStates();
        double[] posterior = [.. masses.Select((m, i) => m + levels[i])];
        return MulticanonicalWeights.Flat(BassFlatHistogramSpace.Quantile(levels, masses, 0.01), BassFlatHistogramSpace.Quantile(levels, posterior, 0.99),
                                          lowSlope, nodes);
    }

    // Each chain's self-normalized estimate of g over its own draws under the fit's weights, the pooled estimate, and its
    // error from the chains' spread, sd/√C.
    private static (double Mean, double Error) Reweighted(BassFit fit, Func<BassDraw, double> g)
    {
        double[] w = fit.DrawWeights!;
        var chains = new double[fit.Draws.Count];
        double num = 0.0, den = 0.0;
        int i = 0;
        for (int c = 0; c < chains.Length; c++)
        {
            double cn = 0.0, cd = 0.0;
            foreach (BassDraw d in fit.Draws[c])
            {
                double weight = w[i++];
                cn += weight * g(d);
                cd += weight;
            }
            chains[c] = cn / cd;
            num += cn;
            den += cd;
        }
        double mean = chains.Average();
        return (num / den, Math.Sqrt(chains.Sum(x => (x - mean) * (x - mean)) / (chains.Length - 1) / chains.Length));
    }

    // Centred RMS between two sets of node values over the nodes not flagged empty.
    private static double CentredRms(IReadOnlyList<double> a, IReadOnlyList<double> b, bool[] empty)
    {
        int[] kept = [.. Enumerable.Range(0, a.Count).Where(k => !empty[k])];
        double shift = kept.Average(k => a[k] - b[k]);
        return Math.Sqrt(kept.Average(k => (a[k] - b[k] - shift) * (a[k] - b[k] - shift)));
    }

    [Theory]
    [InlineData(true, false, 0.0)]
    [InlineData(false, false, 0.0)]
    [InlineData(true, true, 0.0)]
    [InlineData(false, true, 0.0)]
    [InlineData(false, false, 0.3)]
    public void Fit_OnFixtureA_MatchesTheEnumeration(bool heldTau, bool ridge, double lowSlope)
    {
        BassModel model = FixtureA(heldTau, ridge);
        BassData data = model.Data;
        var space = BassFlatHistogramSpace.Of(model);
        MulticanonicalWeights shape = ExactShape(space, 12, lowSlope);
        (double[] levels, double[] masses) = space.DensityOfStates();
        MulticanonicalWeights exact = SamcLearning.FixedPoint(shape, levels, masses, out bool[] empty);
        BassConfig config = model.Config with { Sampling = model.Config.Sampling with { MasterSeed = 7301 } };
        BassFit fit = Bass.Fit(WithFlatHistogram(config, new BassFlatHistogramSpec { Run = DeclaredRun(shape, 30_000, 30_000) }), data);
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        string label = $"τ {(heldTau ? "held" : "sampled")}, {(ridge ? "ridge" : "g-form")}, b_lo = {lowSlope}";
        _out.WriteLine($"## Fixture A, {label}: {space.States.Count:N0} structures; window [{fh.Window.Low:F2}, {fh.Window.High:F2}], K = {fh.Nodes}, " +
                       $"t₀ = {fh.GainScale}, gain from {fh.GainStart}; {FitChains} chains × 30,000 learning + 30,000 production");
        Assert.True(fh.Complete);
        Assert.Equal((shape.Low, shape.High, shape.Nodes, shape.LowSlope), (fh.Weights.Low, fh.Weights.High, fh.Weights.Nodes, fh.Weights.LowSlope));
        Assert.Equal(10.0 * shape.Nodes, fh.GainScale);

        // The learned weights against SAMC's limit on the exact density of states.
        double rms = CentredRms(fh.Weights.LogWeights, exact.LogWeights, empty);
        _out.WriteLine($"learned ln w against θ*: centred RMS {rms:F4} over {empty.Count(e => !e)} nodes (tolerance 0.15); " +
                       $"production flatness {fh.Production!.MaxFlatnessDeviation:F3}, round trips {fh.Production.TotalRoundTrips}, b = 1 Kish ESS " +
                       $"{fh.EffectiveSampleSize:F0} of {FitChains * 30_000}");
        Assert.True(rms <= 0.15, $"learned ln w sits {rms:F4} from θ*");
        Assert.Empty(fh.Production.ChainsWithoutRoundTrip);

        // The reweighted structure posterior against the exact one, within 5 SE of the chains' spread.
        double[] law = space.Law(l => l);
        int cap = model.MaxBasis;
        var countLaw = new double[cap + 1];
        var inclusion = new Dictionary<BassTerm, double>();
        double meanTau = 0.0;
        for (int s = 0; s < law.Length; s++)
        {
            countLaw[space.States[s].Terms.Length] += law[s];
            foreach (BassTerm t in space.States[s].Terms) inclusion[t] = inclusion.GetValueOrDefault(t) + law[s];
            meanTau += law[s] * space.MeanTau(space.States[s], l => l);
        }
        _out.WriteLine("quantity | exact | reweighted ± chains' SE | z | allowance");
        _out.WriteLine("--- | --- | --- | --- | ---");
        void Check(string name, double expected, Func<BassDraw, double> g, double slack)
        {
            (double mean, double se) = Reweighted(fit, g);
            double allowance = 5 * se + slack;
            _out.WriteLine($"{name} | {expected:F4} | {mean:F4} ± {se:F4} | {(se > 0 ? $"{(mean - expected) / se:F2}" : "—")} | {allowance:F4}");
            Assert.True(Math.Abs(mean - expected) <= allowance, $"{name}: reweighted {mean:F4} ± {se:F4}, exact {expected:F4}");
        }
        BassTerm[] Terms(BassDraw d) => [.. d.Structure.TermIds.Select(id => fit.Terms[id])];
        for (int m = 0; m <= cap; m++) Check($"P(M = {m})", countLaw[m], d => d.Structure.TermIds.Length == m ? 1 : 0, 0.002);
        foreach ((BassTerm t, double p) in inclusion.Where(x => x.Value >= 0.005))
            Check($"incl. {BassMixing.Key([t])}", p, d => Terms(d).Contains(t) ? 1 : 0, 0.002);
        foreach (int s in Enumerable.Range(0, law.Length).OrderByDescending(s => law[s]).Take(10))
            Check($"S = {space.States[s].Key}", law[s], d => BassMixing.Key(Terms(d)) == space.States[s].Key ? 1 : 0, 0.002);
        if (heldTau) Assert.All(fit.Draws.SelectMany(c => c), d => Assert.Equal(model.Target.HeldTau!.Value, d.Tau));
        else Check("E[τ]", meanTau, d => d.Tau, 0.005 * meanTau);

        // ln Z(b) − ln Z(b_lo) on the grid against the enumeration, within 5 errors + 0.01; the evidence only from b_lo = 0.
        BassFlatHistogramEvidence curve = fh.Evidence!;
        Assert.Equal(lowSlope, curve.From);
        Assert.Equal(lowSlope == 0.0, curve.IsEvidence);
        double[] grid = curve.Curve.InverseTemperatures;
        Assert.Equal((BassFlatHistogramEvidence.GridPoints, lowSlope, 1.0), (grid.Length, grid[0], grid[^1]));
        _out.WriteLine("b | exact ln Z(b) − ln Z(b_lo) | curve ± error | z");
        _out.WriteLine("--- | --- | --- | ---");
        foreach (int k in new[] { 10, 25, 50, 75, 100 })
        {
            double b = grid[k];
            double expected = space.LogEvidence(b) - space.LogEvidence(lowSlope);
            (double value, double error) = curve.Curve.LogRatio(0, k);
            _out.WriteLine($"{b:F3} | {expected:F4} | {value:F4} ± {error:F4} | {(value - expected) / error:F2}");
            Assert.True(Math.Abs(value - expected) <= 5 * error + 0.01, $"ln Z({b}) − ln Z(b_lo): {value:F4} ± {error:F4} against {expected:F4}");
        }
        Assert.Equal(curve.Curve.LogRatio(0, grid.Length - 1), (curve.LogRatio, curve.MonteCarloError));
        if (lowSlope == 0.0)
        {
            Assert.Equal(0.0, curve.Curve.LogEvidence[0]);
            Assert.Equal(curve.LogRatio, fit.Evidence!.LogEvidence);
            Assert.Equal(curve.MonteCarloError, fit.Evidence.MonteCarloError);
            Assert.Null(fit.Evidence.SteppingStone);
            Assert.Same(curve, fit.Evidence.FlatHistogram);
        }
        else Assert.Null(fit.Evidence);
    }

    // ---- The window's top: raised by the learning, or a margin above the pilot's draws (the loose ends' item 1) ----

    [Fact]
    public void Fit_OnFixtureA_FromALowTop_RaisedByTheLearning_MatchesTheEnumeration()
    {
        // The exact shape's spacing (12 nodes from the prior's 1% quantile of ℓ₁ to the posterior's 99%), its top three
        // spacings lower, τ held: raised by the visits over the first 5,000 of 30,000 iterations, the window ends a whole
        // number of spacings higher, inside F's highest ℓ₁; the learned weights meet SAMC's limit over it, and ln Z(b)
        // the enumeration's, within the tolerances of the declared fits above.
        BassModel model = FixtureA(heldTau: true);
        var space = BassFlatHistogramSpace.Of(model);
        MulticanonicalWeights full = ExactShape(space, 12);
        var low = MulticanonicalWeights.Flat(full.Low, full.High - 3 * full.Spacing, 0.0, 9);
        (double[] levels, double[] masses) = space.DensityOfStates();
        BassConfig config = model.Config with { Sampling = model.Config.Sampling with { MasterSeed = 7311 } };
        FlatHistogramSpec run = DeclaredRun(low, 30_000, 30_000) with { ExtendTopUntil = 5_000 };
        BassFit fit = Bass.Fit(WithFlatHistogram(config, new BassFlatHistogramSpec { Run = run }), model.Data);
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        MulticanonicalWeights w = fh.Weights;
        int added = w.Nodes - low.Nodes;
        Assert.True(added >= 1, "the window should rise");
        Assert.Equal(low.High + added * low.Spacing, w.High, 9);
        Assert.Equal(low.Spacing, w.Spacing, 9);
        Assert.True(w.High <= levels.Max());
        Assert.Equal(low.Nodes, fh.WeightTrace[0].Weights.Nodes);
        Assert.Same(w, fh.Draws!.Weights);

        MulticanonicalWeights exact = SamcLearning.FixedPoint(MulticanonicalWeights.Flat(w.Low, w.High, 0.0, w.Nodes), levels, masses, out bool[] empty);
        double rms = CentredRms(w.LogWeights, exact.LogWeights, empty);
        _out.WriteLine($"window [{low.Low:F2}, {low.High:F2}] raised to [{w.Low:F2}, {w.High:F2}] (+{added} nodes; F's highest ℓ₁ {levels.Max():F2}); learning's " +
                       $"visits above {(double)fh.Learning.AboveWindow / fh.Learning.Visits:P1}; learned ln w to θ* {rms:F4} (tolerance 0.15); production " +
                       $"round trips {fh.Production!.TotalRoundTrips}, b = 1 Kish ESS {fh.EffectiveSampleSize:F0}");
        Assert.True(rms <= 0.15, $"learned ln w sits {rms:F4} from θ*");
        BassFlatHistogramEvidence curve = fh.Evidence!;
        foreach (int k in new[] { 25, 50, 100 })
        {
            double b = curve.Curve.InverseTemperatures[k], expected = space.LogEvidence(b) - space.LogEvidence(0.0);
            (double value, double error) = curve.Curve.LogRatio(0, k);
            _out.WriteLine($"b = {b:F2}: ln Z {value:F4} ± {error:F4} against {expected:F4}");
            Assert.True(Math.Abs(value - expected) <= 5 * error + 0.01, $"ln Z({b}): {value:F4} ± {error:F4} against {expected:F4}");
        }
    }

    [Fact]
    public void Fit_WithATopMargin_StartsOverThePilotsWindowRaisedUnheld()
    {
        // The gap pilot's fit (K = 20) with a margin above of 2 sd, not held inside the pilot's highest draw: the run's
        // start is the window the pilot gives with that margin, and the fit completes with its evidence. From a flat start:
        // on fixture A with τ held the pilot's highest draw is F's highest ℓ₁, so the margin's nodes hold no mass, SAMC's
        // limit there does not exist, and the WHAM start over that window is refused (read in exploration).
        BassModel model = FixtureA(heldTau: true);
        BassConfig gap = GapConfig(model, 7803, GapPilot);
        BassFlatHistogramSpec spec = gap.Sampling.FlatHistogram!;
        BassConfig config = WithFlatHistogram(gap, spec with { Run = spec.Run with { Start = FlatHistogramStart.Flat, WindowTopMargin = 2.0 } });
        BassFit fit = Bass.Fit(config, model.Data);
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        TemperingPilot pilot = PilotOf(model, GapPilot, SeedTree.Derive(7803, 2)[1]).Pilot;
        FlatHistogramWindow held = Wham.Window(pilot, 0.01, 1.0), raised = Wham.Window(pilot, 0.01, 1.0, 2.0);
        MulticanonicalWeights start = fh.WeightTrace[0].Weights;
        _out.WriteLine($"window held [{held.Low:F2}, {held.High:F2}], with the margin [{raised.Low:F2}, {raised.High:F2}]; start [{start.Low:F2}, {start.High:F2}], " +
                       $"K {start.Nodes}; ln Z(1) {fit.Evidence!.LogEvidence:F3} ± {fit.Evidence.MonteCarloError:F3}");
        Assert.Equal((raised.Low, raised.High, 20), (start.Low, start.High, start.Nodes));
        Assert.All(start.LogWeights, v => Assert.Equal(0.0, v));
        Assert.True(raised.High > held.High);
        Assert.True(fh.Complete);
    }

    // The oracle's fixtures B and C, as BassTargetParityTests and the Friedman reference harness build them.
    private static BassConfig OracleConfig(bool categorical)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = categorical ? 1 : 2, Support = new SupportRule(Minimum: 3), MaxBasis = 2 },
            Intercept = GPriorCentre.Mean,
            Sigma2 = Sigma2Prior.ScaledToResponse(1e-3),
            Coefficients = BassCoefficientPrior.G(),
        };

    private static double Hinge(double x, double knot, int sign) => Math.Max(sign * (x - knot), 0.0);

    private static BassData OracleData(bool categorical)
    {
        var rows = new List<double[]>();
        if (!categorical)
        {
            var rng = new Xoshiro256PlusPlus(seed: 7);
            for (int rep = 0; rep < 2; rep++)
                for (int j = 0; j < 4; j++)
                    for (int i = 0; i < 4; i++)
                    {
                        double x1 = i / 3.0, x2 = j / 3.0;
                        rows.Add([1 + 2 * Hinge(x1, 1 / 3.0, 1) + 3 * Hinge(x1, 1 / 3.0, 1) * Hinge(x2, 2 / 3.0, -1) + 0.25 * BassQualification.Gaussian(rng), x1, x2]);
                    }
        }
        else
        {
            var rng = new Xoshiro256PlusPlus(seed: 21);
            for (int z = 0; z < 3; z++)
                for (int j = 0; j < 3; j++)
                    for (int i = 0; i < 3; i++)
                    {
                        double x1 = i * 0.5, x2 = j * 0.5;
                        rows.Add([1 + 1.5 * x1 + 1.5 * Hinge(x2, 0.5, 1) * (z == 0 ? 1 : 0) + 0.4 * BassQualification.Gaussian(rng), x1, x2, z]);
                    }
        }
        var predictors = new List<BassPredictor>
        {
            new ContinuousPredictor("x1", [.. rows.Select(r => r[1])]),
            new ContinuousPredictor("x2", [.. rows.Select(r => r[2])]),
        };
        if (categorical) predictors.Add(new CategoricalPredictor("z", [.. rows.Select(r => (int)r[3])], ["0", "1", "2"]));
        return new BassData(predictors, new ScalarResponse([.. rows.Select(r => r[0])]));
    }

    [Theory]
    [InlineData("B")]
    [InlineData("C")]
    public void Fit_OnFixturesBAndC_ReadsTheOraclesEvidence(string fixture)
    {
        bool categorical = fixture == "C";
        double oracle = categorical ? -30.903 : -28.589;
        BassData data = OracleData(categorical);
        BassConfig config = OracleConfig(categorical) with { Sampling = new BassSampling { MasterSeed = 7401 } };
        var space = BassFlatHistogramSpace.Of(BassModel.Create(config, data));
        double enumerated = space.LogEvidence(1.0);
        Assert.True(Math.Abs(enumerated - oracle) <= 0.002, $"the enumeration reads {enumerated:F4} against the oracle's {oracle}");
        MulticanonicalWeights shape = ExactShape(space, 12);
        BassFit fit = Bass.Fit(WithFlatHistogram(config, new BassFlatHistogramSpec { Run = DeclaredRun(shape, 30_000, 30_000) }), data);
        BassEvidence e = fit.Evidence!;
        double z = (e.LogEvidence - oracle) / e.MonteCarloError;
        _out.WriteLine($"## Fixture {fixture}: {space.States.Count:N0} structures; ln Z(1) {e.LogEvidence:F4} ± {e.MonteCarloError:F4} against the oracle's " +
                       $"{oracle} (enumeration {enumerated:F4}), z {z:F2}; round trips {fit.Diagnostics.FlatHistogram!.Production!.TotalRoundTrips}");
        Assert.True(Math.Abs(e.LogEvidence - oracle) <= 4 * e.MonteCarloError + 0.01, $"ln Z(1) {e.LogEvidence:F4} ± {e.MonteCarloError:F4} against {oracle}");
        Assert.Equal(config.Coefficients, e.Priors.Coefficients);
    }

    // A direct functional response with one output variable: x₀ and x₁ on a four-value grid (n = 16), t on five locations
    // (N = 80), Y = 1 + 2x₀ + 1.5·t·x₁ + 2(x₀ − ½)₊(t − ½)₊ + noise.
    private static BassData FunctionalData()
    {
        var rng = new Xoshiro256PlusPlus(13);
        const int n = 16, m = 5;
        double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i / 4) / 3.0;
            x1[i] = (i % 4) / 3.0;
            for (int r = 0; r < m; r++)
                y[i, r] = 1 + 2 * x0[i] + 1.5 * t[r] * x1[i] + 2 * Math.Max(x0[i] - 0.5, 0) * Math.Max(t[r] - 0.5, 0) + 0.2 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    // A rough window from b_lo = 0: from the 5% quantile of ℓ₁ over 400 exact prior draws to the top a short untempered
    // chain reaches.
    internal static MulticanonicalWeights RoughShape(BassModel model, int nodes)
    {
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(3);
        double[] prior = [.. Enumerable.Range(0, 400).Select(_ => sampler.DrawState(rng).Evaluation.LogMarginal)];
        BassState start = sampler.DrawState(rng).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(4), new Xoshiro256PlusPlus(5));
        double top = double.NegativeInfinity;
        for (int s = 0; s < 2_000; s++) top = Math.Max(top, model.Target.Evaluate(chain.Step()).LogMarginal);
        return MulticanonicalWeights.Flat(BassMixing.Quantile(prior, 0.05), top, 0.0, nodes);
    }

    // A small fit's config: scalar on BassSweepTests' larger grid, or the functional response above.
    private static (BassConfig Config, BassData Data) EndToEnd(bool functional, int seed)
    {
        BassData data = functional ? FunctionalData() : BassSweepTests.Data(large: true);
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 8 },
            Sampling = new BassSampling { MasterSeed = seed },
        };
        return (config, data);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Fit_EndToEnd_EveryViewReadsTheWeightedDraws(bool functional)
    {
        (BassConfig plain, BassData data) = EndToEnd(functional, 7501);
        BassModel model = BassModel.Create(plain, data);
        BassFlatHistogramSpec spec = new() { Run = DeclaredRun(RoughShape(model, 10), 4_000, 1_000, chains: 4) };
        BassConfig config = WithFlatHistogram(plain, spec);
        BassFit fit = Bass.Fit(config, data);
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        const int chains = 4, draws = 1_000;

        // The draws and their weights at b = 1, chain-major, as the readout weighs them.
        Assert.Equal(new RunIdentity(plain.Sampling.MasterSeed, chains, Bass.Fingerprint(config, data)), fit.Identity);
        Assert.All(fit.Draws, store => Assert.Equal(draws, store.Count));
        double[] w = fit.DrawWeights!;
        Assert.Equal(chains * draws, w.Length);
        Assert.Equal(1.0, w.Sum(), 12);
        Assert.Equal(Reweighting.Weigh(fh.Draws!, 1.0).Weights, w);
        Assert.Equal(chains, fh.ChainEffectiveSampleSizes.Length);
        Assert.InRange(fh.EffectiveSampleSize, 1.0, chains * draws);
        Assert.Equal(chains, fh.ProductionChains.Count);
        Assert.Equal(fh.Production!.InWindow, fh.ProductionChains.Sum(c => c.InWindow));
        Assert.Equal(fh.Production.RoundTrips, fh.ProductionChains.Select(c => c.RoundTrips[0]));
        Assert.Equal(spec.Run.LearningIterations, fh.LearningIterations);
        Assert.Equal(spec.Run.LearningIterations * chains, fh.Learning.Visits);
        Assert.Equal(chains * (long)draws, fh.Production.Visits);
        Assert.Equal(0, fh.WeightTrace[0].Iteration);
        Assert.Equal(spec.Run.LearningIterations, fh.WeightTrace[^1].Iteration);
        Assert.Same(fh.Weights, fh.Draws!.Weights);

        // Each draw's weight is e^{ℓ₁}/w(ℓ₁) of its own state, its ℓ₁ rebuilt from its structure and τ, up to the one
        // normalizer: the weights align with the stores' pooled order.
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        double offset = double.NaN;
        for (int s = 0; s < pooled.Length; s += 37)
        {
            BassTerm[] terms = [.. pooled[s].Structure.TermIds.Select(id => fit.Terms[id])];
            double level = model.Target.Evaluate(BassMixing.Loaded(model, terms, Math.Log(pooled[s].Tau), out _)!).LogMarginal;
            double logNormalizer = level - fh.Weights.LogWeight(level) - Math.Log(w[s]);
            if (double.IsNaN(offset)) offset = logNormalizer;
            else Assert.Equal(offset, logNormalizer, 8);
        }

        // π_w's diagnostics over production: the verdict's three groups, and each chain's account.
        BassDiagnostics d = fit.Diagnostics;
        Assert.Equal(draws, d.SamplesUsed);
        Assert.Equal([BassDiagnostics.PredictiveGroup, BassDiagnostics.LevelGroup, BassDiagnostics.StructuralGroup], d.Verdict.Groups.Select(g => g.Name));
        Assert.Equal(chains, d.Chains.Count);
        Assert.True(double.IsFinite(d.ModeSeparation.LogPosteriorGap));
        Assert.InRange(d.MassAtCap, 0.0, 1.0);
        Assert.Equal(d.AcceptanceRate, (double)d.Moves.Sum(m => m.Accepted) / d.Moves.Sum(m => m.Attempts), 12);
        Assert.Equal(d.Moves.Sum(m => m.Attempts), d.Chains.Sum(a => a.Moves.Sum(m => m.Attempts)));
        Assert.Null(d.Tempering);
        Assert.Null(d.Smc);
        _out.WriteLine($"{(functional ? "functional" : "scalar")}: window [{fh.Window.Low:F1}, {fh.Window.High:F1}], K {fh.Nodes}; round trips " +
                       $"{string.Join(", ", fh.Production.RoundTrips)}; flatness {fh.Production.MaxFlatnessDeviation:F2}; ESS {fh.EffectiveSampleSize:F0}; " +
                       $"verdict {d.Verdict.Outcome}; ln Z(1) {fit.Evidence!.LogEvidence:F3} ± {fit.Evidence.MonteCarloError:F3}");

        // Every view runs on the weighted draws.
        BassRunLengths lengths = fit.Metadata().Lengths;
        Assert.Equal((chains, (int?)null, (int?)null, (int?)1, (SmcSpec?)null), (lengths.Chains, lengths.BurnIn, lengths.Samples, lengths.Thin, lengths.Smc));
        Assert.Equal(spec, lengths.FlatHistogram);
        Assert.All(lengths.DrawsPerChain, n => Assert.Equal(draws, n));
        Assert.Equal(Enumerable.Range(0, pooled.Length).Sum(s => w[s] * pooled[s].Sigma2), fit.Summary().Posterior.Sigma2.Mean, 12);
        Assert.Equal(Enumerable.Range(0, pooled.Length).Sum(s => w[s] * pooled[s].Structure.TermIds.Length), fit.Summary().Posterior.Count.Mean, 12);
        foreach (PointwiseContent content in Enum.GetValues<PointwiseContent>())
            Assert.Equal(w, fit.PointwiseLogDensity(data, content: content).Weights!.ToArray());
        BassPrediction prediction = fit.Predict(data, PredictionCases.Of(data));
        Assert.Equal(pooled.Length, prediction.DrawCount);
        Assert.Equal(w, fit.Sensitivity(data).Decompositions[0].Weights);
        fit.Effects(data, effects: new EffectSpec { GridPoints = 10 });
        BassTrainingFit training = fit.TrainingFit(data);
        if (!functional)
        {
            double[][] fitted = BassQualification.Fitted(BassModel.Create(config, data), fit);
            for (int i = 0; i < training.Cases; i++)
                Assert.Equal(Enumerable.Range(0, pooled.Length).Sum(s => w[s] * fitted[s][i]), training.Mean[i, 0], 9);
            Assert.Null(fit.ResidualDependence);   // a scalar response has no output index
        }
        else Assert.NotNull(fit.ResidualDependence);
        Assert.Throws<InvalidOperationException>(() => fit.Traces());
        Assert.Null(fit.Checkpoint);
    }

    [Fact]
    public void Fit_ResumesBitForBit_AcrossLearningAndProduction()
    {
        // Each iteration two sweeps (Sampling.Thin), so the account's restart at production is counted in sweeps.
        (BassConfig thinned, BassData data) = EndToEnd(functional: false, 7601);
        BassConfig plain = thinned with { Sampling = thinned.Sampling with { Thin = 2 } };
        BassModel model = BassModel.Create(plain, data);
        BassFlatHistogramSpec spec = new() { Run = DeclaredRun(RoughShape(model, 8), 300, 200, chains: 3) };
        BassFit whole = Bass.Fit(WithFlatHistogram(plain, spec), data);

        // A call budget of 170: mid-learning, then across into production, then the rest.
        BassConfig cut = WithFlatHistogram(plain, spec with { MaxIterations = 170 });
        var run = new BassRunOptions { Checkpoint = true };
        BassFit first = Bass.Fit(cut, data, run);
        BassFlatHistogramReadout readout = first.Diagnostics.FlatHistogram!;
        Assert.False(readout.Complete);
        Assert.Equal(170L, readout.LearningIterations);
        Assert.Null(readout.Production);
        Assert.All(first.Draws, store => Assert.Empty(store));
        Assert.Null(first.Evidence);
        Assert.Null(first.ResidualDependence);
        Assert.Equal(ConvergenceOutcome.Undetermined, first.Diagnostics.Verdict.Outcome);
        Assert.Throws<InvalidOperationException>(() => first.Summary());
        Assert.Throws<InvalidOperationException>(() => first.Predict(data, PredictionCases.Of(data)));
        Assert.Throws<InvalidOperationException>(() => first.PointwiseLogDensity(data));
        Assert.Throws<InvalidOperationException>(() => first.Sensitivity(data));
        Assert.Throws<InvalidOperationException>(() => first.Traces());
        first.Metadata();
        Assert.NotNull(first.Checkpoint);

        BassFit second = Bass.Resume(first.Checkpoint!, cut, data, run);
        Assert.False(second.Diagnostics.FlatHistogram!.Complete);
        Assert.Equal(300L, second.Diagnostics.FlatHistogram.LearningIterations);
        Assert.All(second.Draws, store => Assert.Equal(40, store.Count));
        Assert.Throws<InvalidOperationException>(() => second.Summary());
        BassFit third = Bass.Resume(second.Checkpoint!, cut, data, run);
        Assert.True(third.Diagnostics.FlatHistogram!.Complete);
        Assert.Null(third.Checkpoint);

        // Bit for bit with the uninterrupted run: the draws and accounts, the weights, the learning and production records.
        Assert.Equal(BassBitIdentityTests.Hash(whole), BassBitIdentityTests.Hash(third));
        Assert.Equal(whole.DrawWeights, third.DrawWeights);
        BassFlatHistogramReadout a = whole.Diagnostics.FlatHistogram!, b = third.Diagnostics.FlatHistogram!;
        Assert.Equal(a.Weights.LogWeights, b.Weights.LogWeights);
        Assert.Equal(a.WeightTrace.Select(t => (t.Iteration, t.Weights.LogWeights.ToArray())), b.WeightTrace.Select(t => (t.Iteration, t.Weights.LogWeights.ToArray())));
        foreach ((FlatHistogramDiagnostics x, FlatHistogramDiagnostics y) in new[] { (a.Learning, b.Learning), (a.Production!, b.Production!) })
        {
            Assert.Equal(x.Occupancy, y.Occupancy);
            Assert.Equal(x.RoundTrips, y.RoundTrips);
            Assert.Equal(x.LastEnd, y.LastEnd);
            Assert.Equal((x.InWindow, x.BelowWindow, x.AboveWindow), (y.InWindow, y.BelowWindow, y.AboveWindow));
        }
        Assert.Equal(whole.Evidence!.LogEvidence, third.Evidence!.LogEvidence);
        Assert.Equal(whole.Diagnostics.Count.RHat, third.Diagnostics.Count.RHat);
        Assert.Equal(2, third.Metadata().Lengths.Thin);
        // The accounts cover production's 400 sweeps (a null proposal is no attempt), the learning's 600 left out.
        Assert.All(third.Diagnostics.Chains, a => Assert.InRange(a.Moves.Sum(m => m.Attempts), 1, 400));
        for (int c = 0; c < 3; c++)
        {
            Assert.Equal(whole.Diagnostics.Chains[c].MoveKinds.Entries, third.Diagnostics.Chains[c].MoveKinds.Entries);
            Assert.Equal(whole.Diagnostics.Chains[c].Support.Entries, third.Diagnostics.Chains[c].Support.Entries);
        }

        // Another config's checkpoint, and an MCMC one, are refused.
        Assert.Throws<ArgumentException>(() => Bass.Resume(second.Checkpoint!, WithFlatHistogram(plain, spec), data));
        BassConfig mcmc = plain with { Sampling = plain.Sampling with { Chains = 2, BurnIn = 20, Samples = 10 } };
        BassFit chains = Bass.Fit(mcmc, data, run);
        Assert.Throws<ArgumentException>(() => Bass.Resume(chains.Checkpoint!, WithFlatHistogram(mcmc, spec with { MaxIterations = 170 }), data));
    }

    [Fact]
    public void Fit_UnderImproperSigmaOrFlatIntercept_RunsFromZero_WithNoEvidence()
    {
        foreach (BassConfig config in new[]
                 {
                     FixtureA(heldTau: false).Config with { Sigma2 = Sigma2Prior.Jeffreys },
                     FixtureA(heldTau: false).Config with { Intercept = GPriorCentre.Flat },
                 })
        {
            BassModel model = BassModel.Create(config, BassSweepTests.Data());
            MulticanonicalWeights shape = ExactShape(BassFlatHistogramSpace.Of(model), 10);
            BassFit fit = Bass.Fit(WithFlatHistogram(config, new BassFlatHistogramSpec { Run = DeclaredRun(shape, 2_000, 500, chains: 2) }), model.Data);
            BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
            Assert.True(fh.Complete);
            Assert.Equal(0.0, fh.LowSlope);
            Assert.Null(fh.Evidence);
            Assert.Null(fit.Evidence);
            Assert.Equal(1.0, fit.DrawWeights!.Sum(), 12);
            fit.Summary();
        }
    }

    [Fact]
    public void Diagnostics_FlagAChainThatNeverTunnels()
    {
        BassModel model = FixtureA(heldTau: false);
        MulticanonicalWeights shape = ExactShape(BassFlatHistogramSpace.Of(model), 10);
        // Production too short to cross the window: three draws per chain.
        BassFit fit = Bass.Fit(WithFlatHistogram(model.Config, new BassFlatHistogramSpec { Run = DeclaredRun(shape, 2_000, 3, chains: 2) }), model.Data);
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        Assert.Equal([0, 1], fh.Production!.ChainsWithoutRoundTrip);
        Assert.All(fh.ProductionChains, c => Assert.Equal([0L], c.RoundTrips));
        Assert.Equal(0, fh.Production.TotalRoundTrips);
    }

    // ---- Step 3: the warm start on BASS ----

    // Fixture A's pilot: two linear-form, reference-rung, self-sizing ladders, burned 255 scans and traced 2,000.
    internal static readonly BassFlatHistogramPilot FixtureAPilot = new() { Chains = 2, BurnIn = 255, TracedScans = 2_000 };

    internal static (TemperingPilot Pilot, BassTemperingReadout Ladders) PilotOf(BassModel model, BassFlatHistogramPilot spec, int seed)
        => BassFlatHistogramFit.Pilot(model, spec, seed, new BassRunOptions(), new BassSweepCount(new BassRunOptions(), 0));

    // The mean over nodes of z², z the pilots' mean departure from θ* over its standard error from their spread, each
    // pilot's ln w centred over the nodes first (WhamTests' measure); and the per-pilot centred RMS.
    private static (double MeanZ2, double Rms) Agreement(IReadOnlyList<MulticanonicalWeights> starts, MulticanonicalWeights exact, bool[] empty)
    {
        int[] kept = [.. Enumerable.Range(0, exact.Nodes).Where(k => !empty[k])];
        double exactCentre = kept.Average(k => exact.LogWeights[k]);
        double[] centres = [.. starts.Select(w => kept.Average(k => w.LogWeights[k]))];
        int p = starts.Count;
        double z2 = 0.0, rms = 0.0;
        foreach (int k in kept)
        {
            double[] d = [.. Enumerable.Range(0, p).Select(i => starts[i].LogWeights[k] - centres[i] - (exact.LogWeights[k] - exactCentre))];
            double mean = d.Average();
            double se = Math.Sqrt(d.Sum(v => (v - mean) * (v - mean)) / (p - 1) / p);
            z2 += mean * mean / (se * se);
            rms += d.Sum(v => v * v) / p;
        }
        return (z2 / kept.Length, Math.Sqrt(rms / kept.Length));
    }

    [Fact]
    public void Pilot_PoolsItsLaddersRungs_AndItsWhamStartMeetsTheExactLimit_WithinThePilotsSpread()
    {
        BassModel model = FixtureA(heldTau: false);
        var space = BassFlatHistogramSpace.Of(model);
        MulticanonicalWeights shape = ExactShape(space, 12);
        (double[] levels, double[] masses) = space.DensityOfStates();
        MulticanonicalWeights exact = SamcLearning.FixedPoint(shape, levels, masses, out bool[] empty);

        // One pilot: each ladder on the linear form, sized itself, frozen with b = 1 and the reference rung's b = 0, which
        // the ladders share and the pilot pools.
        (TemperingPilot pilot, BassTemperingReadout ladders) = PilotOf(model, FixtureAPilot, 7701);
        Assert.Equal((TemperingPath.Likelihood, TemperedForm.Linear, true), (ladders.Path, ladders.Form, ladders.ReferenceRung));
        Assert.All(ladders.Ladders, rungs => Assert.Equal((1.0, 0.0), (rungs[0], rungs[^1])));
        double[] distinct = [.. ladders.Ladders.SelectMany(r => r).Distinct()];
        Assert.Equal(distinct, pilot.InverseTemperatures);
        Assert.Equal(0.0, pilot.HotEnd);
        Assert.Equal(ladders.Levels.Sum() * (long)FixtureAPilot.TracedScans, pilot.Draws);
        Assert.Equal(2 * FixtureAPilot.TracedScans, pilot.Levels(pilot.ColdRung).Length);
        Assert.All(pilot.Inefficiencies, g => Assert.True(g >= 1.0));
        WhamEstimate estimate = Wham.Estimate(pilot);
        Assert.True(estimate.Residual < 1e-8, $"WHAM's fixed point holds to {estimate.Residual:E1}");
        _out.WriteLine($"## Fixture A's pilot: ladders {string.Join(" | ", ladders.Ladders.Select(r => string.Join(", ", r.Select(b => b.ToString("F3")))))}; " +
                       $"{pilot.Rungs} pooled rungs, {pilot.Draws:N0} draws, g {string.Join(", ", pilot.Inefficiencies.Select(g => g.ToString("F1")))}; " +
                       $"{ladders.RoundTripsPerScan:F3} round trips per scan");

        // Eight independent pilots' WHAM starts against θ* on the exact density of states.
        List<MulticanonicalWeights> starts = [.. SeedTree.Derive(7702, 8).Select(s => Wham.WarmStart(PilotOf(model, FixtureAPilot, s).Pilot, shape))];
        (double z2, double rms) = Agreement(starts, exact, empty);
        _out.WriteLine($"8 pilots' WHAM starts against θ*: per-pilot centred RMS {rms:F3} (tolerance 0.25); mean z² over nodes {z2:F2} (tolerance 3; t₇ gives 1.4)");
        Assert.True(z2 < 3.0, $"the warm starts depart from θ* by more than their spread: mean z² {z2:F2}");
        Assert.True(rms < 0.25, $"a warm start departs from θ* by {rms:F3} RMS");
        // Pinned before the warm start's remedy (2026-10-31): a start the fixed point accepts stays as it was, bit for bit.
        Assert.Equal("95ad874b5dce0d4e", StartHash(starts));
    }

    // FNV-1a over the starts' windows, slopes and ln w, bit for bit.
    internal static string StartHash(IEnumerable<MulticanonicalWeights> starts)
    {
        var hash = new RunFingerprint();
        foreach (MulticanonicalWeights w in starts) hash.Add(w.Low).Add(w.High).Add(w.LowSlope).Add([.. w.LogWeights]);
        return hash.ToString();
    }

    // Fixture A, τ sampled: 48 pilots (FixtureAPilot, master seeds from 7905), each with the window, b_lo and K a run
    // resolves from it, its WHAM ĝ, and SAMC's limit on the exact density of states over that shape.
    internal sealed record ResolvedPilot(int Seed, TemperingPilot Pilot, MulticanonicalWeights Shape, WhamEstimate Estimate, MulticanonicalWeights Exact,
                                         bool[] ExactEmpty);

    internal static readonly Lazy<ResolvedPilot[]> FixtureAPilots = new(() =>
    {
        BassModel model = FixtureA(heldTau: false);
        (double[] levels, double[] masses) = BassFlatHistogramSpace.Of(model).DensityOfStates();
        return [.. SeedTree.Derive(7905, 48).Select(seed =>
        {
            TemperingPilot pilot = PilotOf(model, FixtureAPilot, seed).Pilot;
            FlatHistogramWindow window = Wham.Window(pilot, 0.01, 1.0);
            var shape = MulticanonicalWeights.Flat(window.Low, window.High, pilot.HotEnd, Wham.Nodes(pilot, window));
            MulticanonicalWeights exact = SamcLearning.FixedPoint(shape, levels, masses, out bool[] empty);
            return new ResolvedPilot(seed, pilot, shape, Wham.Estimate(pilot), exact, empty);
        })];
    }, System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);

    [Fact]
    public void WarmStart_OnFortyEightFixtureAPilots_TiesTheSevenTheFixedPointRefuses_AndLeavesTheRestAsTheyWere()
    {
        // Fixture A's pilots over the shapes they resolve themselves, the window's margins reaching the hot rung's sparse
        // tail: SAMC's fixed point refuses 7 of 48 ĝ, each infeasible at the window's lowest segments (read in exploration,
        // 2026-10-31). The accepted starts were pinned before the warm start's remedy, which leaves them bit for bit; the
        // refused seven, their nodes under ten draws tied, are accepted. Read in exploration: the tied starts' pooled
        // centred RMS to θ* 0.29, against 0.84 for the 41 accepted as they stood, whose lowest nodes sit up to 17 nats off.
        // Phase 1's 0.25 was read on denser pilots over the exact window, and neither group meets it here; the bound below,
        // 0.35, was set after that reading.
        var accepted = new List<MulticanonicalWeights>();
        var acceptedRms = new List<double>();
        var refused = new List<int>();
        var tiedRms = new List<double>();
        foreach (ResolvedPilot p in FixtureAPilots.Value)
        {
            MulticanonicalWeights? plain = null;
            try { plain = SamcLearning.FixedPoint(p.Shape, p.Estimate.Levels, p.Estimate.LogMasses, out _); }
            catch (InvalidOperationException)
            {
                refused.Add(p.Seed);
                (int at, _) = SamcLearning.Infeasibility(p.Shape, p.Estimate.Levels, p.Estimate.LogMasses);
                Assert.InRange(at, 0, 2);
            }
            MulticanonicalWeights start = Wham.WarmStart(p.Pilot, p.Shape, Wham.ThinDraws, out int[] tied);
            double rms = CentredRms(start.LogWeights, p.Exact.LogWeights, p.ExactEmpty);
            if (plain is null)
            {
                Assert.NotEmpty(tied);
                tiedRms.Add(rms);
                _out.WriteLine($"  seed {p.Seed}: K {p.Shape.Nodes}, tied [{string.Join(", ", tied)}], centred RMS to θ* {rms:F3}");
                continue;
            }
            Assert.Empty(tied);
            Assert.Equal(StartHash([plain]), StartHash([start]));
            accepted.Add(start);
            acceptedRms.Add(rms);
        }
        double Pooled(List<double> values) => Math.Sqrt(values.Average(v => v * v));
        _out.WriteLine($"## Fixture A, 48 pilots over their own shapes: {refused.Count} refused by the fixed point ({string.Join(", ", refused)}), none by " +
                       $"the warm start; pooled centred RMS to θ* {Pooled(tiedRms):F3} tied, {Pooled(acceptedRms):F3} for the {accepted.Count} accepted as they stood");
        Assert.Equal(7, refused.Count);
        // Pinned before the warm start's remedy (2026-10-31): a start the fixed point accepts stays as it was, bit for bit.
        Assert.Equal("e09cb5b3b470f47e", StartHash(accepted));
        Assert.True(Pooled(tiedRms) < 0.35, $"the tied starts depart from θ* by {Pooled(tiedRms):F3} RMS");
    }

    [Fact]
    public void Pilot_PoolsSharedRungs_ByTheirEffectiveDraws()
    {
        var a = new TemperingPilot([1.0, 0.4, 0.0], [[-1.0, -2.0, -3.0, -4.0], [-5.0, -6.0], [-9.0, -8.0]], [2.0, 1.0, 1.5]);
        var b = new TemperingPilot([1.0, 0.5, 0.0], [[-1.5, -2.5], [-5.5, -6.5, -7.5], [-9.5, -8.5, -7.0]], [1.0, 3.0, 1.0]);
        TemperingPilot pooled = BassFlatHistogramFit.Pooled([a, b]);
        Assert.Equal([1.0, 0.4, 0.0, 0.5], pooled.InverseTemperatures);
        Assert.Equal([-1.0, -2.0, -3.0, -4.0, -1.5, -2.5], pooled.Levels(0).ToArray());
        Assert.Equal([-9.0, -8.0, -9.5, -8.5, -7.0], pooled.Levels(2).ToArray());
        Assert.Equal(6.0 / (4.0 / 2.0 + 2.0 / 1.0), pooled.Inefficiencies[0], 15);
        Assert.Equal(1.0, pooled.Inefficiencies[1]);
        Assert.Equal(5.0 / (2.0 / 1.5 + 3.0 / 1.0), pooled.Inefficiencies[2], 15);
        Assert.Equal(3.0, pooled.Inefficiencies[3]);
        Assert.Equal(a.Draws + b.Draws, pooled.Draws);
    }

    // The learned weights' centred RMS to θ* after `learning` iterations of four chains on fixture A, from flat or from the
    // pilot's WHAM start, the gain counted from the pilot's draws over the chains; the mean over three seeds.
    internal static double LearningError(long learning, TemperingPilot? pilot, int seed)
    {
        (BassModel model, MulticanonicalWeights shape, MulticanonicalWeights exact, bool[] empty) = FixtureALimit.Value;
        return SeedTree.Derive(seed, 3).Average(s =>
        {
            FlatHistogramSpec spec = DeclaredRun(shape, learning, 1, chains: 4) with
            {
                Start = pilot is null ? FlatHistogramStart.Flat : FlatHistogramStart.Wham,
            };
            var client = new BassMulticanonicalModel(model, "fixture A", learning, new BassSweepCount(new BassRunOptions(), 0));
            MulticanonicalResult<BassState> run = MulticanonicalRun.Run(spec, model.Target, client, s, pilot, parallel: new ParallelOptions { MaxDegreeOfParallelism = 1 });
            return CentredRms(run.Weights.LogWeights, exact.LogWeights, empty);
        });
    }

    // Fixture A, τ sampled: the model, the exact shape (K = 12 from b_lo = 0) and SAMC's limit on it.
    private static readonly Lazy<(BassModel Model, MulticanonicalWeights Shape, MulticanonicalWeights Exact, bool[] Empty)> FixtureALimit = new(() =>
    {
        BassModel model = FixtureA(heldTau: false);
        var space = BassFlatHistogramSpace.Of(model);
        MulticanonicalWeights shape = ExactShape(space, 12);
        (double[] levels, double[] masses) = space.DensityOfStates();
        return (model, shape, SamcLearning.FixedPoint(shape, levels, masses, out bool[] empty), empty);
    });

    internal static (double Flat, double WarmHalf, double[] FlatCurve, double[] WarmCurve) WarmAgainstFlat(int pilotSeed, int seed, long length)
    {
        TemperingPilot pilot = PilotOf(FixtureALimit.Value.Model, FixtureAPilot, pilotSeed).Pilot;
        double[] flat = [LearningError(length / 4, null, seed), LearningError(length / 2, null, seed), LearningError(length, null, seed)];
        double[] warm = [LearningError(0, pilot, seed), LearningError(length / 4, pilot, seed), LearningError(length / 2, pilot, seed)];
        return (flat[2], warm[2], flat, warm);
    }

    [Fact]
    public void Pilot_LearningFromTheWarmStart_LeadsTheFlatStart_AtEqualIterations()
    {
        // Over 32 chain seeds (ComparativeSweep) the lead holds on 30 at 4,000, while half a flat run's length holds on 17 at
        // 4,000 against 8,000 and on 7 at 16,000 against 32,000: the warm start's learning, its gain counted from the pilot's
        // draws, falls more slowly than flat learning.
        (_, double warmHalf, double[] f, double[] w) = WarmAgainstFlat(7701, 7703, 8_000);
        _out.WriteLine($"## Fixture A, mean centred RMS to θ* over 3 seeds: flat start 2k {f[0]:F3}, 4k {f[1]:F3}, 8k {f[2]:F3}; " +
                       $"warm start 0 {w[0]:F3}, 2k {w[1]:F3}, 4k {w[2]:F3}");
        Assert.True(warmHalf < f[1], $"4,000 iterations from the warm start ({warmHalf:F3}) should beat 4,000 from flat ({f[1]:F3})");
    }

    // The gap pilot: three fixed rungs, b = 1, 0.9 and the reference rung's 0, two ladders traced 500 scans; on fixture A
    // with τ held its ĝ has about 130 atoms, one per structure drawn.
    internal static readonly BassFlatHistogramPilot GapPilot = new()
    {
        Chains = 2, BurnIn = 50, TracedScans = 500,
        Tempering = new TemperingSpec { Form = TemperedForm.Linear, ReferenceRung = true, Levels = 3, Adaptive = false, HotEnd = 0.9 },
    };

    [Fact]
    public void GapPilots_AreRefused_ExactlyWhereNoLimitFlattensTheirHats()
    {
        // The warm start's chip, step 1: why SAMC's fixed point refuses a pilot's ĝ with no node empty or a gap beside a
        // drawn node. ĝ is a discrete measure, its atoms the pilot's draws; flat hat occupancy is reachable by a tilt of it
        // only inside the hull of the atoms' hat vectors (SamcLearning.Infeasibility). Where a node's hat sees ĝ only
        // through atoms far from the node, beside an empty segment or the window's edge, its flat share cannot be met without
        // overfilling its neighbour: no finite θ* exists, G is unbounded below, and Newton's steps grow until none descends.
        // So the refusal is the limit's, not Newton's, and the fixed point refuses exactly the infeasible ĝ, naming the segment
        // where flat occupancy leaves the hull and by how much (the loose ends' item 4).
        BassModel model = FixtureA(heldTau: true);
        foreach ((int seed, bool refused, int segment) in new[] { (7801, true, 13), (7802, true, 16), (7803, false, -1) })
        {
            int[] roots = SeedTree.Derive(seed, 2);   // the run's, the pilot's, as Bass.Fit derives them
            TemperingPilot pilot = PilotOf(model, GapPilot, roots[1]).Pilot;
            FlatHistogramWindow window = Wham.Window(pilot, 0.01, 1.0);
            var shape = MulticanonicalWeights.Flat(window.Low, window.High, pilot.HotEnd, 20);
            WhamEstimate estimate = Wham.Estimate(pilot);
            string verdict = "accepted";
            try { SamcLearning.FixedPoint(shape, estimate.Levels, estimate.LogMasses, out _); }
            catch (InvalidOperationException e) { verdict = e.Message; }
            (int at, double gap) = SamcLearning.Infeasibility(shape, estimate.Levels, estimate.LogMasses);
            _out.WriteLine($"## Gap pilot, seed {seed}: window [{window.Low:F2}, {window.High:F2}], K 20, {estimate.Levels.Length} atoms; {verdict}; " +
                           (at < 0 ? "flat occupancy inside the hull" : $"flat occupancy outside the hull at segment {at}, by {gap:F4}"));
            foreach (int k in at < 0 ? [12, 16] : new[] { at - 2, at - 1, at, at + 1 }.Where(k => k >= 0 && k < 20))
                _out.WriteLine("  " + FlatHistogramLimit.NodeReport(shape, pilot, estimate, k));
            Assert.Equal(refused, verdict != "accepted");
            Assert.Equal(segment, at);
            if (refused) Assert.Contains($"leaves the hull of the measure's hats at segment {segment} (ℓ {shape.Node(segment):G6} to {shape.Node(segment + 1):G6}), by {gap:G3}", verdict);
        }
    }

    [Theory]
    [InlineData(7803, false)]
    [InlineData(7802, true)]
    public void Fit_FromAPilotWithAGap_InterpolatesItsEmptyNodes_OrTiesItsThinOnes(int seed, bool tied)
    {
        // Three fixed rungs, b = 1, 0.9 and the reference rung's 0, leave the levels between the prior's and the posterior's
        // thinly drawn; τ held, so ℓ₁'s levels are discrete. Where SAMC's fixed point accepts the pilot's ĝ, as on 7803 with
        // nodes 13–15 empty, the start reads the empty nodes off the line through their neighbours, as before the warm
        // start's remedy (pinned). Where it refuses ĝ, as on 7802, whose node 15 sees ĝ through one draw at 0.24 of its hat
        // beside empty nodes 13 and 14 (the chip's step 1), the nodes under ten pilot draws are tied to their neighbours and
        // the fit runs from that start. Read on seeds 7801–7803 in exploration: 7801 and 7802 refused before the remedy.
        BassModel model = FixtureA(heldTau: true);
        BassFit fit = Bass.Fit(GapConfig(model, seed, GapPilot), model.Data);
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        BassFlatHistogramPilotReadout read = fh.Pilot!;
        IReadOnlyList<double> start = fh.WeightTrace[0].Weights.LogWeights;
        _out.WriteLine($"window [{fh.Window.Low:F2}, {fh.Window.High:F2}], K {fh.Nodes}; empty nodes {string.Join(", ", read.EmptyNodes)}; start " +
                       string.Join(", ", start.Select(v => v.ToString("F2"))));
        Assert.Equal([1.0, 0.9, 0.0], read.InverseTemperatures);
        Assert.All(read.Ladders.Levels, n => Assert.Equal(3, n));
        Assert.Equal(2L * 3 * 500, read.Draws);
        Assert.Equal(0.0, fh.LowSlope);
        Assert.Equal(20, fh.Nodes);
        Assert.True(fh.Complete);
        Assert.NotNull(fit.Evidence);

        // The run's start is the warm start on the pilot rebuilt from the seed, over the shape the run resolved.
        TemperingPilot pilot = PilotOf(model, GapPilot, SeedTree.Derive(seed, 2)[1]).Pilot;
        var shape = MulticanonicalWeights.Flat(fh.Window.Low, fh.Window.High, fh.LowSlope, fh.Nodes);
        Assert.Equal(StartHash([Wham.WarmStart(pilot, shape, Wham.ThinDraws, out int[] thin)]), StartHash([fh.WeightTrace[0].Weights]));
        Assert.Equal(Wham.NodeDraws(pilot, shape), read.NodeDraws);
        Assert.Equal(thin, read.TiedNodes);   // the readout says whether the start was tied (the loose ends' item 4)
        if (!tied)
        {
            Assert.Empty(thin);
            int[] gap = read.EmptyNodes;
            Assert.Equal([13, 14, 15], gap);
            int below = gap[0] - 1, above = gap[^1] + 1;
            foreach (int k in gap)
                Assert.Equal(start[below] + (start[above] - start[below]) * (k - below) / (above - below), start[k], 9);
            // Pinned before the warm start's remedy (2026-10-31): a start the fixed point accepts stays as it was, bit for bit.
            Assert.Equal("33f442de917867dc", StartHash([fh.WeightTrace[0].Weights]));
            return;
        }

        // Tied: the nodes under ten draws by the readout's own count, each interior run on the line through the drawn nodes
        // either side, and node 0 below the first drawn node at b_lo = 0, level with it.
        _out.WriteLine($"tied nodes [{string.Join(", ", thin)}]");
        Assert.Equal([.. Enumerable.Range(0, fh.Nodes).Where(k => read.NodeDraws[k] < Wham.ThinDraws)], thin);
        Assert.Equal([0, 7, 10, 11, 13, 14, 15, 16], thin);
        Assert.Equal(start[1], start[0]);

        // A fit cut short and resumed reports the same tied nodes, from its checkpoint.
        BassConfig cut = GapConfig(model, seed, GapPilot);
        cut = WithFlatHistogram(cut, cut.Sampling.FlatHistogram! with { MaxIterations = 1_000 });
        var options = new BassRunOptions { Checkpoint = true };
        BassFit first = Bass.Fit(cut, model.Data, options);
        Assert.Equal(thin, first.Diagnostics.FlatHistogram!.Pilot!.TiedNodes);
        BassFit resumed = Bass.Resume(first.Checkpoint!, cut, model.Data, options);
        Assert.Equal(thin, resumed.Diagnostics.FlatHistogram!.Pilot!.TiedNodes);
        foreach ((int below, int above) in new[] { (6, 8), (9, 12), (12, 17) })
            for (int k = below + 1; k < above; k++)
                Assert.Equal(start[below] + (start[above] - start[below]) * (k - below) / (above - below), start[k], 9);
    }

    // Today's fits from the gap pilots' WHAM starts, accepted (7803) and tied (7802), pinned before the window's top became
    // adjustable (the loose ends' item 1, 2026-10-31): with the options left off, a fit stays as it was, bit for bit.
    [Fact]
    public void TodaysFits_FromTheGapPilots_ArePinned()
    {
        BassModel model = FixtureA(heldTau: true);
        string accepted = FitHash(Bass.Fit(GapConfig(model, 7803, GapPilot), model.Data));
        string tied = FitHash(Bass.Fit(GapConfig(model, 7802, GapPilot), model.Data));
        _out.WriteLine($"gap pilot 7803 (accepted start) {accepted}; 7802 (tied start) {tied}");
        Assert.Equal("f323760aeb7d24ee", accepted);
        Assert.Equal("a98e0802bb4b83c1", tied);
    }

    // FNV-1a over a flat-histogram fit, bit for bit: the weights and their trace, both phases' visits, the draws' levels and
    // their weights at b = 1.
    internal static string FitHash(BassFit fit)
    {
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        var hash = new RunFingerprint();
        foreach (LearningSnapshot s in fh.WeightTrace)
            hash.Add((ulong)s.Iteration).Add(s.Weights.Low).Add(s.Weights.High).Add(s.Weights.LowSlope).Add([.. s.Weights.LogWeights]);
        hash.Add([.. fh.Weights.LogWeights]);
        foreach (FlatHistogramDiagnostics? d in new[] { fh.Learning, fh.Production })
        {
            if (d is null) continue;
            hash.Add(d.Occupancy).Add((ulong)d.InWindow).Add((ulong)d.BelowWindow).Add((ulong)d.AboveWindow);
            foreach (long r in d.RoundTrips) hash.Add((ulong)r);
        }
        for (int c = 0; c < fh.Draws!.Chains; c++) hash.Add(fh.Draws.Levels(c));
        return hash.Add(fit.DrawWeights!).ToString();
    }

    [Fact]
    public void Fit_FromAPilotTooThinToTie_IsRefused_NamingItsWindowAndWhy()
    {
        // The gap pilot traced for 3 scans, 18 draws over 20 nodes: the fixed point refuses its ĝ, and no two nodes hold ten
        // draws to tie the rest to, so the warm start is refused and the fit throws, naming the window, K, b_lo, the empty
        // nodes and why (read on seeds 7801–7806 in exploration: every one refused so).
        BassModel model = FixtureA(heldTau: true);
        BassConfig config = GapConfig(model, 7801, GapPilot with { TracedScans = 3 });
        var refusal = Assert.Throws<InvalidOperationException>(() => Bass.Fit(config, model.Data));
        _out.WriteLine(refusal.Message);
        Assert.Contains("WHAM start", refusal.Message);
        Assert.Contains("K = 20, b_lo = 0", refusal.Message);
        Assert.Contains("fewer than two nodes hold 10 of the pilot's draws", refusal.Message);
        Assert.Contains("flat occupancy leaves the hull of the measure's hats at segment 3 (ℓ -32.5261 to -31.5102)", refusal.Message);
        Assert.IsType<InvalidOperationException>(refusal.InnerException);
    }

    // The gap pilot's fit on fixture A, τ held, K = 20, from the WHAM start.
    private static BassConfig GapConfig(BassModel model, int seed, BassFlatHistogramPilot pilot) => model.Config with
    {
        Sampling = model.Config.Sampling with
        {
            MasterSeed = seed,
            FlatHistogram = new BassFlatHistogramSpec
            {
                Run = new FlatHistogramSpec { LearningIterations = 3_000, ProductionDraws = 500, Start = FlatHistogramStart.Wham, Nodes = 20 },
                Pilot = pilot,
            },
        },
    };

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

    /// <summary>E[τ | S] under p(u)·e^{f(ℓ₁)}.</summary>
    public double MeanTau(Structure s, Func<double, double> logWeight) => Mean(s, logWeight, k => Math.Exp(U[k]));

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

/// <summary>
/// A node's view of a pilot's ĝ (the warm start's chip, step 1); whether SAMC's limit exists on a discrete measure is the
/// engine's <see cref="SamcLearning.Infeasibility"/>.
/// </summary>
internal static class FlatHistogramLimit
{
    /// <summary>Node k of a pilot's ĝ over a shape: the pilot's draws under its hat, and on each side the segment's draws
    /// and atoms, their fractions u and their log mass under the hat.</summary>
    internal static string NodeReport(MulticanonicalWeights shape, TemperingPilot pilot, WhamEstimate estimate, int k)
    {
        double draws = 0.0;
        int below = 0, above = 0;
        for (int i = 0; i < pilot.Rungs; i++)
            foreach (double level in pilot.Levels(i))
            {
                if (!shape.TryLocate(level, out int s, out double u)) continue;
                if (s == k) { draws += 1.0 - u; above++; }
                else if (s == k - 1) { draws += u; below++; }
            }
        string Side(int segment, bool upper, int segmentDraws)
        {
            if (segment < 0 || segment > shape.Nodes - 2) return "none (the window's edge)";
            var u = new List<double>();
            var logHat = new List<double>();
            for (int j = 0; j < estimate.Levels.Length; j++)
                if (shape.TryLocate(estimate.Levels[j], out int s, out double f) && s == segment)
                {
                    u.Add(f);
                    logHat.Add(estimate.LogMasses[j] + Math.Log(upper ? 1.0 - f : f));
                }
            if (u.Count == 0) return $"{Count(segmentDraws, "draw")}, no atom";
            double max = logHat.Max();
            double mass = max + Math.Log(logHat.Sum(v => Math.Exp(v - max)));
            return $"{Count(segmentDraws, "draw")}, {Count(u.Count, "atom")} at u {u.Min():F3}–{u.Max():F3}, ln mass under the hat {mass:F2}";
        }
        static string Count(int n, string noun) => n == 1 ? $"1 {noun}" : $"{n} {noun}s";
        return $"node {k} (ℓ {shape.Node(k):F2}): {draws:F1} pilot draws under its hat; below, {Side(k - 1, false, below)}; above, {Side(k, true, above)}";
    }
}
