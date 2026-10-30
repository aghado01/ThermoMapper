using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Friedman: a path from a pilot reference (the Friedman reference brief). SMC along the geometric path from a reference q
/// to BASS's posterior, ln γ_b = (1 − b)·ln q + b·ln π̃, with q(S, u) = [ε·p₀(S) + (1 − ε)·e(S)]·p(u): D01's unconstrained
/// structure measure beside the pilots' empirical structures, τ from its prior. The simplest discrete form of Surjanovic
/// et al. 2022's variational reference, read on S10's Friedman fixture against SMC from the prior.
/// <list type="bullet">
/// <item><b>The target.</b> γ_b on the likelihood path, ln q the untempered prior and Δ = ln π̃ − ln q the likelihood, so
/// the SMC core tempers it unchanged. <see cref="BassTarget"/> states the structure prior as D02's count weight w on F,
/// unnormalized; π̃ takes it as p₀ = w/C restricted to F, C the count law's normalizer, so the run's normalizer is
/// Z_F·p(y), Z_F being F's share of p₀, and ln p(y) is the normalizer less ln Z_F, read from
/// <see cref="BassPriorSampler"/>'s acceptance with its error carried. A draw of p₀ outside F has γ_b = 0 at every b &gt; 0,
/// so it falls at the first stage.</item>
/// <item><b>The kernel.</b> C4's sweep over each particle's own design (<see cref="BassSmcKernel"/>): D01's palette through a
/// reversible-jump chain, then τ's slice step, with γ_b in the model target's place. The proposals and their ratios are
/// D01's, unchanged.</item>
/// <item><b>The checks.</b> The kernel run as a chain on fixture A at b = 0.3 and 0.7 against γ_b by enumeration (τ by
/// quadrature); q's draws against its density, each component and mixed; SMC from a reference on the oracle's fixtures B
/// and C against their exact ln Z; and ε = 1, SMC from p₀, against C4's from-prior SMC on fixture B.</item>
/// <item><b>The runs.</b> The pilots (M4's four untempered chains, S10's tempered fit's cold rungs, fitting SMC from the
/// prior at 4,000 particles) make e, with or without the SMC pilot; ε 0.05 and 0.5; 1,000 and 4,000 particles; four
/// replicates from master seed 53, on seeds no pilot uses.</item>
/// </list>
/// It reads, it does not judge (<see cref="BassQualification.Enabled"/>). It asserts its own soundness: the kernel's law,
/// q's, the evidence on the enumerable fixtures, the weights normalized and the basins' fractions summing to one.
/// </summary>
public sealed class BassFriedmanReferenceTests
{
    private const int Replicates = 4;
    private const string A1 = "4 × x1·x2 + 2 × x3 + x4 + x5", A2 = "4 × x1·x2 + x2 + 2 × x3 + x4 + x5";
    private static readonly double[] Marks = [0.0, 0.25, 0.5, 0.75];   // 0: the first stage past b = 0
    private static readonly SmcSpec Defaults = new();   // C4's: conditional ESS, sweeps, resampling, the stage budget
    private readonly ITestOutputHelper _out;
    public BassFriedmanReferenceTests(ITestOutputHelper output) => _out = output;

    // ---- The checks ----

    [Theory]
    [InlineData(0.3)]
    [InlineData(0.7)]
    public void Kernel_SamplesTheBridgeOnFixtureA(double b)
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassModel model = BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data());
        var reference = new Reference(model, 0.5, ChainPilot(model, 4101));
        Enumerated space = Enumerated.Of(model);
        double[] gamma = space.Bridge(reference, b);
        double[] plain = space.Bridge(new Reference(model, 1.0, []), b);   // the likelihood path's linear form, q = p₀

        // The enumerated law: M, each term's inclusion, the ten most probable structures, E[u].
        int cap = model.MaxBasis;
        var countLaw = new double[cap + 1];
        var inclusion = new Dictionary<BassTerm, double>();
        double meanU = 0.0;
        for (int s = 0; s < space.States.Count; s++)
        {
            Enumerated.Structure st = space.States[s];
            countLaw[st.Terms.Length] += gamma[s];
            foreach (BassTerm t in st.Terms) inclusion[t] = inclusion.GetValueOrDefault(t) + gamma[s];
            meanU += gamma[s] * space.MeanU(st, b);
        }
        BassTerm[] tracked = [.. inclusion.Where(x => x.Value >= 0.005).Select(x => x.Key)];
        int[] top = [.. Enumerable.Range(0, gamma.Length).OrderByDescending(s => gamma[s]).Take(10)];
        var topIndex = new Dictionary<string, int>();
        for (int k = 0; k < top.Length; k++) topIndex.Add(space.States[top[k]].Key, k);
        var termIndex = new Dictionary<BassTerm, int>();
        for (int k = 0; k < tracked.Length; k++) termIndex.Add(tracked[k], k);

        // The kernel as a chain at b from a draw of q in F, four chains.
        const int chains = 4, burn = 1_000, sweeps = 30_000, batches = 25;
        int items = cap + 1 + tracked.Length + top.Length + 1;
        var batch = new Batches(items, chains, batches, sweeps);
        var bridge = new Bridge(model, reference);
        var target = new ParticleTarget(bridge);
        var kernel = new Kernel(bridge);
        var draws = new Particles(model, reference);
        var stage = new SmcStage<RefParticle>(target, new Temperature(b, TemperingPath.Likelihood));
        int[] seeds = SeedTree.Derive(4102, chains);
        var clock = Stopwatch.StartNew();
        Parallel.For(0, chains, c =>
        {
            var rng = new Xoshiro256PlusPlus(seeds[c]);
            RefParticle p;
            do p = draws.Draw(rng); while (p.Inside is null);
            ITargetEvaluation e = target.Evaluate(p);
            for (int s = 0; s < burn + sweeps; s++)
            {
                SmcMove<RefParticle> move = kernel.Move(p, e, stage, rng);
                (p, e) = (move.State, move.Evaluation);
                if (s < burn) continue;
                int at = s - burn;
                BassState state = p.Inside!.State;
                batch.Add(state.Count, c, at, 1.0);
                foreach (BassTerm t in state.Terms)
                    if (termIndex.TryGetValue(t, out int k)) batch.Add(cap + 1 + k, c, at, 1.0);
                if (topIndex.TryGetValue(BassMixing.Key(state.Terms), out int j)) batch.Add(cap + 1 + tracked.Length + j, c, at, 1.0);
                batch.Add(items - 1, c, at, state.LogTau);
                if (at % 97 == 0) Assert.Equal(stage.LogDensity(bridge.Evaluate(state)), stage.LogDensity(e));   // the move's evaluation is its state's
            }
        });

        _out.WriteLine($"## Kernel on fixture A at b = {b}: {space.States.Count:N0} structures in F; reference ε = 0.5 over " +
                       $"{reference.Structures} pilot structures; {chains} × {sweeps:N0} sweeps in {clock.Elapsed.TotalSeconds:F1} s; " +
                       $"γ_b against the likelihood path's linear form from p₀: total variation {0.5 * gamma.Zip(plain).Sum(x => Math.Abs(x.First - x.Second)):F3}");
        _out.WriteLine("quantity | γ_b | chain ± batch SE | z | allowance");
        _out.WriteLine("--- | --- | --- | --- | ---");
        void Check(string name, double expected, int item, double slack)
        {
            (double mean, double se) = batch.Read(item);
            double allowance = 5 * se + slack;
            _out.WriteLine($"{name} | {expected:F4} | {mean:F4} ± {se:F4} | {(se > 0 ? $"{(mean - expected) / se:F2}" : "—")} | {allowance:F4}");
            Assert.True(Math.Abs(mean - expected) <= allowance, $"{name}: the chain reads {mean:F4} ± {se:F4}, γ_b {expected:F4}");
        }
        for (int m = 0; m <= cap; m++) Check($"P(M = {m})", countLaw[m], m, 0.004);
        for (int k = 0; k < tracked.Length; k++) Check($"incl. {BassMixing.Key([tracked[k]])}", inclusion[tracked[k]], cap + 1 + k, 0.005);
        for (int k = 0; k < top.Length; k++) Check($"S = {space.States[top[k]].Key}", gamma[top[k]], cap + 1 + tracked.Length + k, 0.005);
        Check("E[u]", meanU, items - 1, 0.001);
    }

    [Fact]
    public void Reference_DrawsMatchItsDensityOnFixtureA()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassModel model = BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data());
        List<PilotDraw> pilot = ChainPilot(model, 4201);
        Enumerated space = Enumerated.Of(model);
        double logZF = space.LogMass - Reference.CountNormalizer(model);
        _out.WriteLine($"## q's draws on fixture A: {space.States.Count:N0} structures in F, Z_F = {Math.Exp(logZF):F6} by enumeration");

        // u's quadrature against its own normalization, then its mean and quartiles.
        double logOne = space.LogIntegral(null, 0.0);
        Assert.True(Math.Abs(logOne) <= 1e-6, $"p(u) integrates to e^{logOne:R} on the grid");
        (double meanU, double[] quartiles) = space.TauMoments();

        int[] components = SeedTree.Derive(4202, 3);
        foreach ((string name, double epsilon, int seed) in new[] { ("p₀", 1.0, components[0]), ("e", 0.0, components[1]), ("mixed", 0.5, components[2]) })
        {
            var reference = new Reference(model, epsilon, pilot);
            // q_S on F, which with ε(1 − Z_F) outside sums to one.
            double[] q = [.. space.States.Select(s => Math.Exp(reference.LogStructure(s.Terms, s.LogW - reference.LogC)))];
            double outside = epsilon * (1.0 - Math.Exp(logZF));
            Assert.True(Math.Abs(q.Sum() + outside - 1.0) <= 1e-9, $"{name}: q sums to {q.Sum() + outside:R} over F and its complement");

            const int streams = 8, perStream = 50_000, n = streams * perStream;
            int cap = model.MaxBasis;
            int[] top = [.. Enumerable.Range(0, q.Length).OrderByDescending(s => q[s]).Take(15)];
            var index = new Dictionary<string, int>();
            for (int k = 0; k < top.Length; k++) index.Add(space.States[top[k]].Key, k);
            var counts = new long[streams, cap + 2 + top.Length + 4];   // M in F, outside, top structures, u's quartiles, then u's sum
            var sums = new double[streams, 2];
            int[] seeds = SeedTree.Derive(seed, streams);
            Parallel.For(0, streams, k =>
            {
                var rng = new Xoshiro256PlusPlus(seeds[k]);
                for (int i = 0; i < perStream; i++)
                {
                    (_, BassTerm[]? terms) = reference.DrawStructure(rng);
                    double u = reference.DrawLogTau(rng);
                    sums[k, 0] += u;
                    sums[k, 1] += u * u;
                    for (int j = 0; j < 3; j++)
                        if (u <= quartiles[j]) counts[k, cap + 2 + top.Length + j]++;
                    if (terms is null || BassMixing.Loaded(model, terms, u, out _) is not { } loaded)
                    {
                        counts[k, cap + 1]++;
                        continue;
                    }
                    counts[k, loaded.Count]++;
                    if (index.TryGetValue(BassMixing.Key(loaded.Terms), out int t)) counts[k, cap + 2 + t]++;
                }
            });
            long Total(int column) => Enumerable.Range(0, streams).Sum(k => counts[k, column]);
            _out.WriteLine("");
            _out.WriteLine($"### {name} (ε = {epsilon}), {n:N0} draws");
            _out.WriteLine("quantity | q | draws ± SE | z");
            _out.WriteLine("--- | --- | --- | ---");
            void Check(string label, double expected, double observed, double se)
            {
                _out.WriteLine($"{label} | {expected:F5} | {observed:F5} ± {se:F5} | {(observed - expected) / Math.Max(se, 1e-12):F2}");
                Assert.True(Math.Abs(observed - expected) <= 5 * se + 1e-12, $"{name}, {label}: {observed:F5} ± {se:F5} against {expected:F5}");
            }
            void Share(string label, double p, long hits) => Check(label, p, (double)hits / n, Math.Sqrt(p * (1 - p) / n));
            var countLaw = new double[cap + 1];
            for (int s = 0; s < q.Length; s++) countLaw[space.States[s].Terms.Length] += q[s];
            for (int m = 0; m <= cap; m++) Share($"M = {m} in F", countLaw[m], Total(m));
            Share("outside F", outside, Total(cap + 1));
            for (int t = 0; t < top.Length; t++) Share($"S = {space.States[top[t]].Key}", q[top[t]], Total(cap + 2 + t));
            for (int j = 0; j < 3; j++) Share($"u ≤ its {25 * (j + 1)}% point", 0.25 * (j + 1), Total(cap + 2 + top.Length + j));
            double sum = Enumerable.Range(0, streams).Sum(k => sums[k, 0]), squares = Enumerable.Range(0, streams).Sum(k => sums[k, 1]);
            double mean = sum / n, sd = Math.Sqrt(squares / n - mean * mean);
            Check("E[u]", meanU, mean, sd / Math.Sqrt(n));
        }
    }

    [Theory]
    [InlineData("B", 0.05)]
    [InlineData("B", 0.5)]
    [InlineData("C", 0.05)]
    [InlineData("C", 0.5)]
    public void Evidence_MatchesTheExactLogZ(string fixture, double epsilon)
    {
        if (!BassQualification.Enabled) return;   // heavy
        bool categorical = fixture == "C";
        double oracle = categorical ? -30.903 : -28.589;
        BassConfig config = OracleConfig(categorical);
        BassData data = OracleData(categorical ? FixtureC() : FixtureB(), categorical);
        BassModel model = BassModel.Create(config, data);
        Enumerated space = Enumerated.Of(model);
        (double exact, double exactZF) = Exact(space);
        Assert.True(Math.Abs(exact - oracle) <= 0.0015, $"the enumeration's ln Z {exact:F4} against the oracle's {oracle:F3}");
        ZF zf = EstimateZF(model, SeedTree.Derive(2026, 7)[6]);
        Assert.True(Math.Abs(zf.Log - exactZF) <= 5 * zf.Error, $"ln Z_F {zf.Log:F5} ± {zf.Error:F5} against {exactZF:F5}");

        BassFit pilotFit = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 2, BurnIn = 500, Samples = 500, MasterSeed = 3001 } }, data);
        var reference = new Reference(model, epsilon, [.. FromFit("pilot", pilotFit)]);
        int[] seeds = SeedTree.Derive(SeedTree.Derive(2026, 5)[4], Replicates);
        var clock = Stopwatch.StartNew();
        Replicate[] reps = RunReplicates(model, reference, 2000, seeds);
        double seconds = clock.Elapsed.TotalSeconds;
        Assert.All(reps, r => Assert.True(r.Reached));

        // The exact posterior's mass on e's structures: what a population held on them would read, ln Z + ln P(e).
        double[] posterior = space.Bridge(new Reference(model, 1.0, []), 1.0);
        double onPilots = space.States.Select((s, k) => reference.Empirical.ContainsKey(s.Key) ? posterior[k] : 0.0).Sum();
        _out.WriteLine($"## Fixture {fixture}, ε = {epsilon}: exact ln Z {exact:F4} by enumeration (the oracle's {oracle:F3}); ln Z_F {exactZF:F5} " +
                       $"exact, {zf.Log:F5} ± {zf.Error:F5} from {zf.Accepted:N0} of {zf.Attempts:N0} prior attempts; reference over {reference.Structures} " +
                       $"pilot structures, holding {onPilots:P2} of the exact posterior (ln Z + ln P(e) = {exact + Math.Log(onPilots):F4}); " +
                       $"4 × 2,000 particles in {seconds:F1} s");
        _out.WriteLine("replicate | stages | normalizer ± genealogy | ln p(y) | z against the oracle | families | p₀ draws in F | weight on e's structures | distinct structures");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | ---");
        foreach (Replicate r in reps)
        {
            double logp = r.LogNormalizer - zf.Log, error = Math.Sqrt(r.Error * r.Error + zf.Error * zf.Error);
            _out.WriteLine($"{r.Seed} | {r.Stages} | {r.LogNormalizer:F4} ± {r.Error:F4} | {logp:F4} | {(logp - oracle) / error:F2} | {r.Families} | " +
                           $"{r.PriorInF:N0} of {r.PriorDraws:N0} | {r.OnPilots:P1} | {r.Distinct}");
            Assert.True(Math.Abs(logp - oracle) <= 4 * error + 0.01, $"{logp:F4} ± {error:F4} against {oracle:F3}");
        }
        Pooled pooled = Pool([.. reps.Select(r => r.LogNormalizer - zf.Log)], [.. reps.Select(r => r.Error)], zf.Error);
        _out.WriteLine($"pooled {pooled.Mean:F4} ± {pooled.Error:F4} (within {pooled.Within:F4}, between {pooled.Between:F4}, Z_F {zf.Error:F4}); " +
                       $"with the exact Z_F {pooled.Mean + zf.Log - exactZF:F4}");
        Assert.True(Math.Abs(pooled.Mean - oracle) <= 4 * pooled.Error + 0.01, $"pooled {pooled.Mean:F4} ± {pooled.Error:F4} against {oracle:F3}");
    }

    [Fact]
    public void PriorReference_AgreesWithC4sSmcOnFixtureB()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassConfig config = OracleConfig(categorical: false);
        BassData data = OracleData(FixtureB(), categorical: false);
        BassModel model = BassModel.Create(config, data);
        (double exact, double exactZF) = Exact(Enumerated.Of(model));
        ZF zf = EstimateZF(model, SeedTree.Derive(2026, 7)[6]);
        int[] seeds = SeedTree.Derive(SeedTree.Derive(2026, 5)[4], Replicates);
        Replicate[] mine = RunReplicates(model, new Reference(model, 1.0, []), 2000, seeds);
        Pooled pooled = Pool([.. mine.Select(r => r.LogNormalizer - zf.Log)], [.. mine.Select(r => r.Error)], zf.Error);
        (double meanM, double seM) = MeanCount(mine.Select(r => r.Result.Particles.Select(p => p.Inside?.State.Count ?? 0).ToArray()).ToArray(),
                                               [.. mine.Select(r => r.Result.Weights)]);

        // The share of p₀'s draws in F, over every replicate's initial population, against Z_F.
        long drawn = mine.Sum(r => r.PriorDraws), inF = mine.Sum(r => r.PriorInF);
        double share = (double)inF / drawn, zfExact = Math.Exp(exactZF), shareSe = Math.Sqrt(zfExact * (1 - zfExact) / drawn);
        _out.WriteLine($"## ε = 1 on fixture B: p₀'s draws in F {inF:N0} of {drawn:N0} = {share:F4} ± {shareSe:F4} against Z_F {zfExact:F4}; " +
                       $"exact ln Z {exact:F4}");
        Assert.True(Math.Abs(share - zfExact) <= 5 * shareSe, $"p₀'s share in F {share:F4} against {zfExact:F4}");
        _out.WriteLine($"this harness, q = p₀: ln p(y) {pooled.Mean:F4} ± {pooled.Error:F4} (within {pooled.Within:F4}, between {pooled.Between:F4}); " +
                       $"replicates [{string.Join(", ", mine.Select(r => $"{r.LogNormalizer - zf.Log:F4} ± {r.Error:F4}"))}]; E[M] {meanM:F4} ± {seM:F4}");

        foreach (TemperedForm form in new[] { TemperedForm.Integrated, TemperedForm.Linear })
        {
            BassFit fit = Bass.Fit(config with { Sampling = config.Sampling with { MasterSeed = 2026, Smc = new SmcSpec { Particles = 2000, Form = form } } }, data);
            BassEvidence evidence = fit.Evidence!;
            int n = fit.Draws[0].Count;
            int[][] counts = [.. fit.Draws.Select(store => store.Select(d => d.Structure.TermIds.Length).ToArray())];
            double[][] weights = [.. Enumerable.Range(0, fit.Draws.Count).Select(c => fit.DrawWeights![(c * n)..((c + 1) * n)])];
            (double theirM, double theirSe) = MeanCount(counts, weights);
            double gap = pooled.Mean - evidence.LogEvidence, allowance = 4 * Math.Sqrt(pooled.Error * pooled.Error + evidence.MonteCarloError * evidence.MonteCarloError) + 0.01;
            _out.WriteLine($"C4's from-prior SMC, {form} form: {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4}, E[M] {theirM:F4} ± {theirSe:F4}; " +
                           $"ln p(y)'s gap {gap:F4} (allowance {allowance:F4}), E[M]'s {meanM - theirM:F4}");
            Assert.True(Math.Abs(gap) <= allowance, $"{form}: {pooled.Mean:F4} against C4's {evidence.LogEvidence:F4}");
            Assert.True(Math.Abs(meanM - theirM) <= 5 * Math.Sqrt(seM * seM + theirSe * theirSe) + 0.01, $"{form}: E[M] {meanM:F4} against C4's {theirM:F4}");
        }
        Assert.True(Math.Abs(pooled.Mean - exact) <= 4 * pooled.Error + 0.01, $"{pooled.Mean:F4} ± {pooled.Error:F4} against the exact {exact:F4}");
    }

    // ---- The runs on Friedman ----

    [Fact]
    public void Friedman_PilotsAndTheirReference()
    {
        if (!BassQualification.Enabled) return;   // heavy
        FriedmanPilots pilots = Pilots.Value;
        BassMixing.Fixture fixture = pilots.Fixture;
        pilots.Legend.Print(_out);

        _out.WriteLine("");
        _out.WriteLine($"## The pilots (master seed 53): M4's chains {pilots.M4Seconds:F0} s, S10's tempered fit {pilots.TemperedSeconds:F0} s, " +
                       $"SMC from the prior at 4,000 {pilots.Smc4000Seconds:F0} s and at 1,000 {pilots.Smc1000Seconds:F0} s (timing only)");
        foreach (BassFit fit in new[] { pilots.Smc1000, pilots.Smc4000 })
        {
            BassSmcReadout smc = fit.Diagnostics.Smc!;
            _out.WriteLine($"from-prior SMC, {fit.Draws[0].Count:N0} particles: ln Ẑ [{string.Join(", ", smc.Replicates.Select(x => $"{x.LogNormalizer:F2} ± {x.LogNormalizerError:F2}"))}], " +
                           $"pooled {fit.Evidence!.LogEvidence:F2} ± {fit.Evidence.MonteCarloError:F2}; stages [{string.Join(", ", smc.Replicates.Select(x => x.StageCount))}]; " +
                           $"families [{string.Join(", ", smc.Replicates.Select(x => x.Families))}]");
        }
        string[] sources = ["M4", "S10", "SMC"];
        var sets = new List<(string Name, Dictionary<string, double> Fractions)>();
        foreach (string source in sources)
        {
            PilotDraw[] own = [.. pilots.Draws.Where(d => d.Source == source)];
            sets.Add(($"{source} ({own.Length:N0} draws, {own.Select(d => BassMixing.Key(d.Terms)).Distinct().Count():N0} structures)",
                      Fractions(own.Select(d => (fixture.Basin(d.Terms), d.Mass)))));
        }
        foreach (bool withSmc in new[] { true, false })
        {
            var reference = new Reference(pilots.Model, 0.5, pilots.Draws.Where(d => withSmc || d.Source != "SMC"));
            sets.Add(($"e, {(withSmc ? "all pilots" : "no SMC pilot")} ({reference.Structures:N0} structures)",
                      Fractions(pilots.Draws.Where(d => withSmc || d.Source != "SMC").Select(d => (fixture.Basin(d.Terms), d.Mass)))));
        }
        BasinTable("the pilots' basins (each source's draws pooled; SMC's weighted)", pilots.Legend, [.. sets], logMass: null);

        _out.WriteLine("");
        _out.WriteLine($"Z_F: {Math.Exp(pilots.ZF.Log):F5}, ln Z_F {pilots.ZF.Log:F5} ± {pilots.ZF.Error:F5}, from {pilots.ZF.Accepted:N0} of {pilots.ZF.Attempts:N0} " +
                       $"attempts of BassPriorSampler ({pilots.ZFSeconds:F1} s); ln C {Reference.CountNormalizer(pilots.Model):F5}");

        // The height e lifts q over ε·p₀ on its own structures, h = ln[(1 − ε)e(S)/(ε·p₀(S))]: a move off them pays (1 − b)·h at b.
        _out.WriteLine("");
        _out.WriteLine("e's height over ε·p₀ on its structures, h = ln[(1 − ε)·e(S)/(ε·p₀(S))], weighted by e: its 10%, 50% and 90% points");
        foreach (bool withSmc in new[] { true, false })
            foreach (double epsilon in new[] { 0.05, 0.5 })
            {
                var reference = new Reference(pilots.Model, epsilon, pilots.Draws.Where(d => withSmc || d.Source != "SMC"));
                var heights = new DrawSample();
                foreach ((BassTerm[] terms, double logE) in reference.Support)
                {
                    double logP0 = pilots.Model.Prior.LogCollectionWeight(terms.Length) + terms.Sum(pilots.Model.Prior.LogPrior) - reference.LogC;
                    heights.Add(Math.Log(1 - epsilon) + logE - Math.Log(epsilon) - logP0, Math.Exp(logE));
                }
                _out.WriteLine($"{(withSmc ? "all pilots" : "no SMC pilot")}, ε = {epsilon}: {heights.Quantile(0.1):F1}, {heights.Quantile(0.5):F1}, {heights.Quantile(0.9):F1}");
            }

        // The evidence over e's own structures, ln Σ_{S∈e} p₀(S)·∫ p(u) L(S, u) du − ln Z_F, by quadrature in u: what a
        // population held on them reads, with no move.
        _out.WriteLine("");
        var grid = Enumerated.Grid(pilots.Model.Target);
        foreach (bool withSmc in new[] { true, false })
        {
            var reference = new Reference(pilots.Model, 0.5, pilots.Draws.Where(d => withSmc || d.Source != "SMC"));
            (BassTerm[] Terms, double LogE)[] support = [.. reference.Support];
            var logMass = new double[support.Length];
            var basins = new string[support.Length];
            Parallel.For(0, support.Length, k =>
            {
                BassState loaded = BassMixing.Loaded(pilots.Model, support[k].Terms, 0.0, out _) ?? throw new InvalidOperationException("A pilot's structure lies outside F.");
                logMass[k] = reference.LogP0(loaded.Fit) + Enumerated.LogIntegral(pilots.Model.Target, loaded.Fit, grid);
                basins[k] = fixture.Basin(support[k].Terms);
            });
            double total = LogSumExp(logMass) - pilots.ZF.Log;
            var byBasin = Enumerable.Range(0, support.Length).GroupBy(k => basins[k])
                                    .Select(g => (Basin: g.Key, Count: g.Count(), Log: LogSumExp(g.Select(k => logMass[k])) - pilots.ZF.Log))
                                    .OrderByDescending(x => x.Log).Take(4).ToArray();
            _out.WriteLine($"the evidence over e's own structures, {(withSmc ? "all pilots" : "no SMC pilot")} ({support.Length:N0}): {total:F2}; by basin " +
                           string.Join(", ", byBasin.Select(x => $"{pilots.Legend.Name(x.Basin)} ({x.Count:N0} structures) {x.Log:F2}")) +
                           $"; the best single structure {logMass.Max() - pilots.ZF.Log:F2}");
        }

        // q on Friedman: p₀'s share in F from this harness's own sampler against Z_F, and e's draws against e.
        var all = new Reference(pilots.Model, 0.5, pilots.Draws);
        const int streams = 16, perStream = 2_500, n = streams * perStream;
        var inF = new long[streams];
        int[] seeds = SeedTree.Derive(SeedTree.Derive(53, 11)[10], streams);
        Parallel.For(0, streams, k =>
        {
            var rng = new Xoshiro256PlusPlus(seeds[k]);
            for (int i = 0; i < perStream; i++)
                if (all.DrawPrior(rng) is { } terms && BassMixing.Loaded(pilots.Model, terms, 0.0, out _) is not null) inF[k]++;
        });
        double z = Math.Exp(pilots.ZF.Log), share = (double)inF.Sum() / n;
        double se = Math.Sqrt(z * (1 - z) / n + z * z * pilots.ZF.Error * pilots.ZF.Error);
        _out.WriteLine($"this harness's p₀ in F: {share:F5} ± {se:F5} of {n:N0} draws against Z_F {z:F5} (z {(share - z) / se:F2})");
        Assert.True(Math.Abs(share - z) <= 5 * se, $"p₀'s share in F {share:F5} against {z:F5}");

        const int draws = 400_000;
        var rngE = new Xoshiro256PlusPlus(SeedTree.Derive(53, 12)[11]);
        var hits = new Dictionary<string, long>();
        for (int i = 0; i < draws; i++)
        {
            BassTerm[] terms = all.DrawEmpirical(rngE);
            string key = BassMixing.Key(terms);
            hits[key] = hits.GetValueOrDefault(key) + 1;
        }
        _out.WriteLine("e's ten heaviest structures: e | draws ± SE | z");
        foreach ((string key, double logE) in all.Empirical.OrderByDescending(x => x.Value).Take(10))
        {
            double p = Math.Exp(logE), f = (double)hits.GetValueOrDefault(key) / draws, sd = Math.Sqrt(p * (1 - p) / draws);
            _out.WriteLine($"{p:F5} | {f:F5} ± {sd:F5} | {(f - p) / sd:F2}");
            Assert.True(Math.Abs(f - p) <= 5 * sd, $"e's draw frequency {f:F5} against {p:F5}");
        }
        Assert.Equal(1.0, all.Empirical.Values.Sum(Math.Exp), 9);
    }

    [Theory]
    [InlineData("all pilots", 0.05, 1000)]
    [InlineData("all pilots", 0.5, 1000)]
    [InlineData("no SMC pilot", 0.05, 1000)]
    [InlineData("no SMC pilot", 0.5, 1000)]
    [InlineData("all pilots", 0.05, 4000)]
    [InlineData("all pilots", 0.5, 4000)]
    [InlineData("no SMC pilot", 0.05, 4000)]
    [InlineData("no SMC pilot", 0.5, 4000)]
    public void Friedman_PathFromPilotReference(string pilotSet, double epsilon, int particles)
    {
        if (!BassQualification.Enabled) return;   // heavy
        FriedmanPilots pilots = Pilots.Value;
        BassMixing.Fixture fixture = pilots.Fixture;
        bool withSmc = pilotSet == "all pilots";
        var reference = new Reference(pilots.Model, epsilon, pilots.Draws.Where(d => withSmc || d.Source != "SMC"));
        var clock = Stopwatch.StartNew();
        Replicate[] reps = RunReplicates(pilots.Model, reference, particles, ReplicateSeeds, fixture.Basin);
        double seconds = clock.Elapsed.TotalSeconds;
        double fromPrior = particles == 1000 ? pilots.Smc1000Seconds : pilots.Smc4000Seconds;
        Report($"Friedman, reference over {pilotSet} ({reference.Structures:N0} structures), ε = {epsilon}, {particles:N0} particles",
               pilots, reps, seconds, fromPrior);
    }

    // The test's replicates: four seeds under master seed 53's ninth child, which no pilot's run derives (a fit takes the
    // first 2·chains children; fitting SMC the first R).
    private static int[] ReplicateSeeds
    {
        get
        {
            int[] taken = SeedTree.Derive(53, 8);
            int root = SeedTree.Derive(53, 9)[8];
            int[] seeds = SeedTree.Derive(root, Replicates);
            Assert.DoesNotContain(root, taken);
            Assert.Empty(seeds.Intersect(taken));
            return seeds;
        }
    }

    private void Report(string title, FriedmanPilots pilots, Replicate[] reps, double seconds, double fromPriorSeconds)
    {
        BassMixing.Fixture fixture = pilots.Fixture;
        Legend legend = pilots.Legend;
        int r = reps.Length;
        double[] logp = [.. reps.Select(x => x.LogNormalizer - pilots.ZF.Log)];
        Pooled pooled = Pool(logp, [.. reps.Select(x => x.Error)], pilots.ZF.Error);

        _out.WriteLine($"## {title}: {seconds:F0} s for {r} replicates in parallel, against {fromPriorSeconds:F0} s from the prior at the same count; " +
                       $"reached b = 1: {reps.All(x => x.Reached)}; stages [{string.Join(", ", reps.Select(x => x.Stages))}] of {Defaults.MaxStages:N0}");
        _out.WriteLine("replicate | seed | stages | normalizer ± genealogy | ln p(y) | families at the first stage, 0.25, 0.5, 0.75, 1 | b at ≤ 2 families | " +
                       "b at 1 | resampled | Kish ESS | p₀ draws in F | σ's 95% interval | weight on e's structures | distinct structures | " +
                       "mean log post (M4's) | s");
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", 15)));
        var pooledSigma = new DrawSample();
        for (int c = 0; c < r; c++)
        {
            Replicate x = reps[c];
            IReadOnlyList<SmcStageReport> stages = x.Result.Stages;
            string At(double b)   // at the first stage reaching b; mark 0 is the first stage past b = 0
            {
                if (b == 0.0) return $"{stages[Math.Min(1, stages.Count - 1)].Families}";
                foreach (SmcStageReport st in stages)
                    if (st.Temperature.InverseTemperature >= b - 1e-12) return $"{st.Families}";
                return "—";
            }
            string First(int families)
            {
                foreach (SmcStageReport st in stages)
                    if (st.Families <= families) return $"{st.Temperature.InverseTemperature:F3}";
                return "—";
            }
            var sigma = new DrawSample();
            double[] w = x.Result.Weights;
            for (int i = 0; i < w.Length; i++)
                if (w[i] > 0.0)
                {
                    sigma.Add(Math.Sqrt(x.Sigma2[i]), w[i]);
                    pooledSigma.Add(Math.Sqrt(x.Sigma2[i]), w[i] / r);
                }
            double kish = 1.0 / w.Sum(v => v * v);
            _out.WriteLine($"{c} | {x.Seed} | {x.Stages} | {x.LogNormalizer:F2} ± {x.Error:F2} | {logp[c]:F2} | " +
                           $"{string.Join(", ", Marks.Select(At))}, {stages[^1].Families} | {First(2)} | {First(1)} | {stages.Count(s => s.Resampled)} | {kish:F0} | " +
                           $"{x.PriorInF:N0} of {x.PriorDraws:N0} | [{sigma.Quantile(0.025):F4}, {sigma.Quantile(0.975):F4}] | {x.OnPilots:P1} | " +
                           $"{x.Distinct:N0} | {x.MeanLogPosterior:F1} | {x.Seconds:F0}");
        }
        _out.WriteLine($"ln p(y): mean {pooled.Mean:F2}, sd {pooled.Sd:F2}, range {pooled.Range:F2}; pooled {pooled.Mean:F2} ± {pooled.Error:F2} " +
                       $"(within {pooled.Within:F2}, between {pooled.Between:F2}, Z_F {pilots.ZF.Error:F3}); σ's 95% interval pooled " +
                       $"[{pooledSigma.Quantile(0.025):F4}, {pooledSigma.Quantile(0.975):F4}] against 1");

        // The final populations' basins, their log masses ln p(y) + ln f, and the share by the count of x1·x2 terms.
        var sets = new (string Name, Dictionary<string, double> Fractions)[r];
        for (int c = 0; c < r; c++)
        {
            Replicate x = reps[c];
            var pairs = new List<(string, double)>();
            for (int i = 0; i < x.Result.Weights.Length; i++)
            {
                if (x.Result.Weights[i] == 0.0) continue;
                BassParticle inside = x.Result.Particles[i].Inside ?? throw new InvalidOperationException("A particle outside F holds weight at b = 1.");
                pairs.Add((fixture.Basin(inside.State.Terms), x.Result.Weights[i]));
            }
            sets[c] = ($"replicate {c}", Fractions(pairs));
        }
        _out.WriteLine("");
        _out.WriteLine("replicate | largest (log mass) | second (log mass) | a1 | a2");
        _out.WriteLine("--- | --- | --- | --- | ---");
        for (int c = 0; c < r; c++)
        {
            var ranked = sets[c].Fractions.OrderByDescending(x => x.Value).Take(2).ToArray();
            string Cell(KeyValuePair<string, double> b) => $"{legend.Name(b.Key)} {b.Value:P1} ({logp[c] + Math.Log(b.Value):F2})";
            _out.WriteLine($"{c} | {Cell(ranked[0])} | {(ranked.Length > 1 ? Cell(ranked[1]) : "—")} | " +
                           $"{sets[c].Fractions.GetValueOrDefault(A1):P2} | {sets[c].Fractions.GetValueOrDefault(A2):P2}");
        }
        BasinTable("the replicates' final populations (b = 1), weighted", legend, sets, logp);

        // Along the path: each mark's population, after resampling and before the stage's mutation.
        _out.WriteLine("");
        _out.WriteLine("### Along the path: the weight on e's own structures, the three largest basins, and the share by the count of x1·x2 terms");
        _out.WriteLine("replicate | mark | b | on e's structures | basins | x1·x2 terms");
        _out.WriteLine("--- | --- | --- | --- | --- | ---");
        for (int c = 0; c < r; c++)
            foreach (View v in reps[c].Views)
                _out.WriteLine($"{c} | {(v.Mark == 0 ? "first" : $"{v.Mark:F2}")} | {v.InverseTemperature:F4} | {v.OnPilots:P1} | " +
                               string.Join(", ", v.Fractions.OrderByDescending(x => x.Value).Take(3).Select(x => $"{legend.Name(x.Key)} {x.Value:P1}")) + " | " +
                               string.Join(", ", v.Fractions.GroupBy(x => Pairs(x.Key)).OrderBy(g => g.Key).Select(g => $"{g.Key}: {g.Sum(x => x.Value):P0}")));

        // Origins: the final weight by its lineage's first draw, its component and its basin.
        _out.WriteLine("");
        _out.WriteLine("replicate | final weight from p₀'s draws | the three heaviest origin basins");
        _out.WriteLine("--- | --- | ---");
        for (int c = 0; c < r; c++)
        {
            Replicate x = reps[c];
            double fromPrior = 0.0;
            var origins = new Dictionary<string, double>();
            for (int i = 0; i < x.Result.Weights.Length; i++)
            {
                double w = x.Result.Weights[i];
                if (w == 0.0) continue;
                RefParticle p = x.Result.Particles[i];
                if (p.FromPrior) fromPrior += w;
                string origin = fixture.Basin(p.Origin!);
                origins[origin] = origins.GetValueOrDefault(origin) + w;
            }
            _out.WriteLine($"{c} | {fromPrior:P2} | " + string.Join(", ", origins.OrderByDescending(o => o.Value).Take(3).Select(o => $"{legend.Name(o.Key)} {o.Value:P1}")));
        }

        // The kernel's acceptance by kind.
        _out.WriteLine("");
        BassKindCount[] kinds = [.. reps[0].Account.Kinds.Entries.Where(k => k.Attempts > 0)];
        _out.WriteLine("replicate | sweeps | " + string.Join(" | ", kinds.Select(k => $"{k.Kind}: proposals, accepted, refused")));
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", 1 + kinds.Length)));
        for (int c = 0; c < r; c++)
        {
            BassSmcAccount account = reps[c].Account;
            _out.WriteLine($"{c} | {account.Sweeps:N0} | " + string.Join(" | ", kinds.Select(k => account.Kinds[k.Kind])
                .Select(k => $"{k.Attempts:N0}, {k.AcceptanceRate:P2}, {(k.Attempts > 0 ? (double)k.Refused / k.Attempts : double.NaN):P2}")));
        }
    }

    // Basins at 1% or more in some sample, with each sample's fraction, and, given each sample's ln p(y), its log mass; the
    // rest as "other"; then the share by the count of x1·x2 terms, and a1's and a2's.
    private void BasinTable(string title, Legend legend, (string Name, Dictionary<string, double> Fractions)[] sets, double[]? logMass)
    {
        foreach ((string name, Dictionary<string, double> f) in sets)
            Assert.True(Math.Abs(f.Values.Sum() - 1.0) <= 1e-9, $"{name}'s fractions sum to {f.Values.Sum():R}");
        string[] shown = [.. sets.SelectMany(s => s.Fractions).Where(x => x.Value >= 0.01).Select(x => x.Key).Distinct().OrderBy(legend.Order, StringComparer.Ordinal)];
        _out.WriteLine("");
        _out.WriteLine($"### Basins: {title}");
        _out.WriteLine($"basin | key | {string.Join(" | ", sets.Select(s => s.Name))}{(logMass is null ? "" : " | " + string.Join(" | ", sets.Select(s => $"log mass, {s.Name}")))}");
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", 1 + sets.Length * (logMass is null ? 1 : 2))));
        foreach (string basin in shown)
            _out.WriteLine($"{legend.Name(basin)} | {basin} | {string.Join(" | ", sets.Select(s => $"{s.Fractions.GetValueOrDefault(basin):P1}"))}" +
                           (logMass is null ? "" : " | " + string.Join(" | ", sets.Select((s, k) => s.Fractions.GetValueOrDefault(basin) > 0 ? $"{logMass[k] + Math.Log(s.Fractions[basin]):F2}" : "—"))));
        _out.WriteLine($"other | — | {string.Join(" | ", sets.Select(s => $"{s.Fractions.Where(x => !shown.Contains(x.Key)).Sum(x => x.Value):P1}"))}" +
                       (logMass is null ? "" : string.Concat(Enumerable.Repeat(" | —", sets.Length))));
        if (logMass is not null)
            for (int k = 0; k < sets.Length; k++)
            {
                double total = LogSumExp(sets[k].Fractions.Values.Where(f => f > 0).Select(f => logMass[k] + Math.Log(f)));
                Assert.True(Math.Abs(total - logMass[k]) <= 1e-9 * Math.Max(1.0, Math.Abs(logMass[k])), $"{sets[k].Name}'s basins' log masses add to {total:R}");
            }
        _out.WriteLine("");
        _out.WriteLine("x1·x2 terms | " + string.Join(" | ", sets.Select(s => s.Name)));
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", sets.Length)));
        foreach (int k in sets.SelectMany(s => s.Fractions.Keys).Select(Pairs).Distinct().Order())
            _out.WriteLine($"{k} | " + string.Join(" | ", sets.Select(s => $"{s.Fractions.Where(x => Pairs(x.Key) == k).Sum(x => x.Value):P2}")));
        _out.WriteLine($"a1 | {string.Join(" | ", sets.Select(s => $"{s.Fractions.GetValueOrDefault(A1):P2}"))}");
        _out.WriteLine($"a2 | {string.Join(" | ", sets.Select(s => $"{s.Fractions.GetValueOrDefault(A2):P2}"))}");
    }

    // ---- The reference, the bridge and the kernel ----

    /// <summary>One pilot draw: its source, its structure and its mass (1 per MCMC draw, N·w per SMC particle).</summary>
    private sealed record PilotDraw(string Source, BassTerm[] Terms, double Mass);

    /// <summary>
    /// The reference q(S, u) = [ε·p₀(S) + (1 − ε)·e(S)]·p(u), drawn exactly and scored exactly. p₀ is the sampler's proposal
    /// in <see cref="BassPriorSampler"/>: M from the count law truncated at the cap, then M terms from D01's term measure, so
    /// a set S of M distinct terms has p₀(S) = w(S)/C. e is the pilots' structures, each by its pooled mass.
    /// </summary>
    private sealed class Reference
    {
        private readonly BassModel _model;
        private readonly double[] _countCdf;
        private readonly Dictionary<string, double> _logE = [];
        private readonly BassTerm[][] _structures;
        private readonly double[] _cumulative;
        private readonly ConditionalWeakTable<BassTerm[], StrongBox<double>> _cache = new();   // ln q_S by a state's term array

        public Reference(BassModel model, double epsilon, IEnumerable<PilotDraw> pilots)
        {
            if (!(epsilon >= 0.0 && epsilon <= 1.0)) throw new ArgumentOutOfRangeException(nameof(epsilon));
            if (model.Target.SamplesTau && !model.Target.Scale.IsProper) throw new ArgumentException("q draws τ from its prior.", nameof(model));
            _model = model;
            Epsilon = epsilon;
            int cap = model.MaxBasis;
            var log = new double[cap + 1];
            for (int m = 0; m <= cap; m++) log[m] = model.Prior.LogCollectionWeight(m) - SpecialFunctions.LogGamma(m + 1.0);
            LogC = LogSumExp(log);
            _countCdf = new double[cap + 1];
            double total = 0.0;
            for (int m = 0; m <= cap; m++) _countCdf[m] = total += Math.Exp(log[m] - LogC);
            _countCdf[cap] = 1.0;

            var mass = new Dictionary<string, (BassTerm[] Terms, double Mass)>();
            foreach (PilotDraw d in pilots)
            {
                if (!(d.Mass > 0.0)) continue;
                string key = BassMixing.Key(d.Terms);
                mass[key] = mass.TryGetValue(key, out var held) ? (held.Terms, held.Mass + d.Mass) : (d.Terms, d.Mass);
            }
            if (epsilon < 1.0 && mass.Count == 0) throw new ArgumentException("ε < 1 needs pilots.", nameof(pilots));
            double sum = mass.Values.Sum(x => x.Mass);
            var keys = mass.Keys.Order(StringComparer.Ordinal).ToArray();
            _structures = [.. keys.Select(k => mass[k].Terms)];
            _cumulative = new double[keys.Length];
            double running = 0.0;
            for (int k = 0; k < keys.Length; k++)
            {
                _logE.Add(keys[k], Math.Log(mass[keys[k]].Mass / sum));
                _cumulative[k] = running += mass[keys[k]].Mass / sum;
            }
            if (_cumulative.Length > 0) _cumulative[^1] = 1.0;
        }

        public double Epsilon { get; }

        /// <summary>ln C, C = Σ_M Γ(h₁ + M)/(Γ(h₁) M! (h₂ + 1)^M) over M ≤ M_max: the count weight's normalizer.</summary>
        public double LogC { get; }

        /// <summary>e's distinct structures.</summary>
        public int Structures => _structures.Length;

        /// <summary>ln e by structure key.</summary>
        public IReadOnlyDictionary<string, double> Empirical => _logE;

        /// <summary>e's structures with their ln e.</summary>
        public IEnumerable<(BassTerm[] Terms, double LogE)> Support => _structures.Select(s => (s, _logE[BassMixing.Key(s)]));

        public static double CountNormalizer(BassModel model) => new Reference(model, 1.0, []).LogC;

        /// <summary>ln p₀ of a structure in F: its count weight and terms' measure, less ln C.</summary>
        public double LogP0(BassStructureFit fit) => _model.Target.LogStructurePrior(fit) - LogC;

        /// <summary>ln q_S = ln[ε·p₀(S) + (1 − ε)·e(S)] for the structure <paramref name="terms"/> with ln p₀ =
        /// <paramref name="logP0"/>.</summary>
        public double LogStructure(IReadOnlyList<BassTerm> terms, double logP0)
        {
            double prior = Math.Log(Epsilon) + logP0;
            if (Epsilon >= 1.0 || !_logE.TryGetValue(BassMixing.Key(terms), out double logE)) return prior;
            double empirical = Math.Log(1.0 - Epsilon) + logE;
            double top = Math.Max(prior, empirical);
            return top + Math.Log(Math.Exp(prior - top) + Math.Exp(empirical - top));
        }

        // LogStructure kept by the state's term array, which a τ step shares.
        public double LogStructureOf(BassTerm[] terms, double logP0)
        {
            if (_cache.TryGetValue(terms, out StrongBox<double>? box)) return box.Value;
            double value = LogStructure(terms, logP0);
            _cache.AddOrUpdate(terms, new StrongBox<double>(value));
            return value;
        }

        /// <summary>A structure from q: the component, and its terms, null when p₀'s draw holds the intercept.</summary>
        public (bool FromPrior, BassTerm[]? Terms) DrawStructure(Xoshiro256PlusPlus rng)
            => Epsilon >= 1.0 || rng.NextDouble() < Epsilon ? (true, DrawPrior(rng)) : (false, DrawEmpirical(rng));

        /// <summary>A collection from p₀: M by inversion of the truncated count law, then M terms from the term measure;
        /// null when a term is the intercept.</summary>
        public BassTerm[]? DrawPrior(Xoshiro256PlusPlus rng)
        {
            double v = rng.NextDouble();
            int m = 0;
            while (v >= _countCdf[m]) m++;
            var terms = new BassTerm[m];
            for (int j = 0; j < m; j++)
            {
                if (_model.Prior.DrawTerm(rng) is not { } term) return null;
                terms[j] = term;
            }
            return terms;
        }

        /// <summary>A structure from e, by inversion of its cumulative masses.</summary>
        public BassTerm[] DrawEmpirical(Xoshiro256PlusPlus rng)
        {
            double v = rng.NextDouble();
            int k = Array.BinarySearch(_cumulative, v);
            k = k < 0 ? ~k : k + 1;
            return _structures[Math.Min(k, _structures.Length - 1)];
        }

        /// <summary>u from τ's prior, or ln of the held τ.</summary>
        public double DrawLogTau(Xoshiro256PlusPlus rng)
            => _model.Target.SamplesTau ? _model.Target.Scale.SampleLogTau(rng) : Math.Log(_model.Target.HeldTau!.Value);
    }

    /// <summary>
    /// γ_b's evaluation on the likelihood path: ln q as the untempered prior and b·Δ, Δ = ln π̃ − ln q, as the likelihood.
    /// Off F at b &gt; 0 it is −∞ throughout; a p₀ draw outside F (<see cref="Outside"/>) holds γ₀ only.
    /// </summary>
    private sealed class Bridged(BassEvaluation? bass, double logPrior, double delta) : ITargetEvaluation
    {
        /// <summary>A p₀ draw outside F: weight 1/N at b = 0 and 0 from the first b &gt; 0, so its ln q, set to 0 here, is
        /// never read by a ratio.</summary>
        public static readonly Bridged Outside = new(null, 0.0, double.NegativeInfinity);

        public BassEvaluation? Bass { get; } = bass;
        public double LogPrior { get; } = logPrior;
        public double Delta { get; } = delta;

        public double LogLikelihood(double inverseTemperature)
        {
            if (double.IsNegativeInfinity(LogPrior)) return double.NegativeInfinity;
            return inverseTemperature == 0.0 ? 0.0 : inverseTemperature * Delta;
        }

        public double LogLikelihoodSlope(double inverseTemperature) => Delta;
    }

    /// <summary>γ_b over BASS's states, the kernel's target at b &gt; 0: a candidate off F, refused, is off the support.</summary>
    private sealed class Bridge(BassModel model, Reference reference) : ISamplingTarget<BassState>
    {
        public BassModel Model { get; } = model;
        public Reference Reference { get; } = reference;

        public Bridged Evaluate(BassState state)
        {
            BassTarget target = Model.Target;
            BassEvaluation e = target.Evaluate(state);
            if (e.OffSupport) return new Bridged(e, double.NegativeInfinity, double.NegativeInfinity);
            double logP0 = Reference.LogP0(state.Fit);
            double logQ = Reference.LogStructureOf(state.TermArray, logP0);
            double logTau = target.SamplesTau ? target.LogTauDensity(state.LogTau) : 0.0;
            return new Bridged(e, logQ + logTau, logP0 - logQ + e.LogMarginal);
        }

        ITargetEvaluation ISamplingTarget<BassState>.Evaluate(BassState state) => Evaluate(state);
    }

    /// <summary>A particle: a BASS particle in F, or none for a p₀ draw outside it; with its lineage's first draw (its
    /// component and its collection), which forks carry.</summary>
    private sealed class RefParticle(BassParticle? inside, bool fromPrior, BassTerm[]? origin)
    {
        public BassParticle? Inside { get; } = inside;
        public bool FromPrior { get; } = fromPrior;
        public BassTerm[]? Origin { get; } = origin;
    }

    private sealed class ParticleTarget(Bridge bridge) : ISamplingTarget<RefParticle>
    {
        public Bridge Bridge { get; } = bridge;
        public ITargetEvaluation Evaluate(RefParticle particle) => particle.Inside is { } inside ? Bridge.Evaluate(inside.State) : Bridged.Outside;
    }

    /// <summary>q's particles: an exact draw of q, loaded into its own design when in F; forks as C4's
    /// (<see cref="BassParticles.Fork"/>). Counts p₀'s draws and those in F.</summary>
    private sealed class Particles(BassModel model, Reference reference) : ISmcParticles<RefParticle>
    {
        private readonly BassParticles _bass = new(model);
        private long _prior, _priorInF, _empirical;

        public long PriorDraws => Interlocked.Read(ref _prior);
        public long PriorInF => Interlocked.Read(ref _priorInF);
        public long EmpiricalDraws => Interlocked.Read(ref _empirical);

        public RefParticle Draw(Xoshiro256PlusPlus rng)
        {
            (bool fromPrior, BassTerm[]? terms) = reference.DrawStructure(rng);
            double u = reference.DrawLogTau(rng);
            if (fromPrior) Interlocked.Increment(ref _prior);
            else Interlocked.Increment(ref _empirical);
            if (terms is null || BassMixing.Loaded(model, terms, u, out _) is not { } loaded)
            {
                if (!fromPrior) throw new InvalidOperationException("A pilot's structure lies outside F.");
                return new RefParticle(null, true, terms);
            }
            if (fromPrior) Interlocked.Increment(ref _priorInF);
            var moves = new BassStructureMoves(model, loaded.Terms, u);
            return new RefParticle(new BassParticle(moves, moves.Start), fromPrior, terms);
        }

        public RefParticle Fork(RefParticle particle)
            => particle.Inside is null ? particle : new RefParticle(_bass.Fork(particle.Inside), particle.FromPrior, particle.Origin);
    }

    /// <summary>A stage's population at a mark: its weighted basins, and its weight on e's own structures.</summary>
    private sealed record View(double Mark, double InverseTemperature, Dictionary<string, double> Fractions, double OnPilots);

    /// <summary>
    /// The mutation at γ_b: <see cref="BassSmcKernel"/>'s sweep (D04) with the bridge in the model target's place. The
    /// structure move on D01's schedule through a reversible-jump chain over the particle's own moves, then τ's slice step,
    /// both against γ_b; the moves' tallies join the account after each sweep. At each mark it records the population's
    /// basins, and it checks that no particle outside F holds weight past b = 0.
    /// </summary>
    private sealed class Kernel(Bridge bridge, Func<IReadOnlyList<BassTerm>, string>? basin = null) : ISmcKernel<RefParticle>
    {
        private readonly BassTauBlock? _tau = bridge.Model.Target.SamplesTau ? new BassTauBlock(bridge) : null;
        private readonly Lock _gate = new();

        public BassSmcAccount Account { get; } = new(bridge.Model.Prior);
        public List<View> Views { get; } = [];

        public void Prepare(SmcStage<RefParticle> stage, IReadOnlyList<RefParticle> particles, ReadOnlySpan<double> weights)
        {
            for (int i = 0; i < particles.Count; i++)
                if (particles[i].Inside is null && weights[i] > 0.0)
                    throw new InvalidOperationException($"A p₀ draw outside F holds weight {weights[i]} at b = {stage.Temperature.InverseTemperature}.");
            if (basin is null) return;
            double b = stage.Temperature.InverseTemperature;
            while (Views.Count < Marks.Length && b >= Marks[Views.Count])
            {
                var fractions = new Dictionary<string, double>();
                double onPilots = 0.0;
                for (int i = 0; i < particles.Count; i++)
                    if (weights[i] > 0.0)
                    {
                        IReadOnlyList<BassTerm> terms = particles[i].Inside!.State.Terms;
                        string key = basin(terms);
                        fractions[key] = fractions.GetValueOrDefault(key) + weights[i];
                        if (bridge.Reference.Empirical.ContainsKey(BassMixing.Key(terms))) onPilots += weights[i];
                    }
                Views.Add(new View(Marks[Views.Count], b, fractions, onPilots));
            }
        }

        public SmcMove<RefParticle> Move(RefParticle particle, ITargetEvaluation evaluation, SmcStage<RefParticle> stage, Xoshiro256PlusPlus rng)
        {
            if (stage.Target is not ParticleTarget target || !ReferenceEquals(target.Bridge, bridge))
                throw new ArgumentException("The stage's target is not this kernel's bridge.", nameof(stage));
            BassParticle inside = particle.Inside ?? throw new InvalidOperationException("A p₀ draw outside F is never moved: its weight is 0 past b = 0.");
            BassStructureMoves moves = inside.Moves;
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, bridge, inside.State, rng, stage.Temperature);
            chain.Step();
            bool accepted = chain.Accepted > 0;
            _tau?.Update(chain, rng);
            inside.State = chain.Current;
            moves.DesignAt(inside.State);   // commits an accepted candidate, which tallies its acceptance
            lock (_gate) Account.Add(chain.MoveStats(), moves);
            moves.ResetTallies();
            return new SmcMove<RefParticle>(particle, chain.CurrentEvaluation, accepted);
        }
    }

    /// <summary>One replicate: the run's result, the kernel's account and marks, each particle's σ² composed at b = 1, q's
    /// draws by component, its time, and the final population's weight on e's own structures, its distinct structures and
    /// its weighted mean log posterior (M4's: the log marginal plus the log structure prior).</summary>
    private sealed record Replicate(int Seed, SmcResult<RefParticle> Result, BassSmcAccount Account, List<View> Views, double[] Sigma2,
                                    long PriorDraws, long PriorInF, long EmpiricalDraws, double Seconds, double OnPilots, int Distinct,
                                    double MeanLogPosterior)
    {
        public double LogNormalizer => Result.LogNormalizerRatio;
        public double Error => Result.LogNormalizerError;
        public int Stages => Result.Stages.Count - 1;
        public int Families => Result.Stages[^1].Families;
        public bool Reached => Result.Stages[^1].Temperature.InverseTemperature == 1.0;
    }

    // SMC from q at C4's defaults, the replicates in parallel within the machine, each particle count's own seeds.
    private static Replicate[] RunReplicates(BassModel model, Reference reference, int particles, int[] seeds,
                                             Func<IReadOnlyList<BassTerm>, string>? basin = null)
    {
        int bound = Environment.ProcessorCount, outer = Math.Min(seeds.Length, bound), inner = Math.Max(1, bound / outer);
        var reps = new Replicate[seeds.Length];
        Parallel.For(0, seeds.Length, new ParallelOptions { MaxDegreeOfParallelism = outer },
                     c => reps[c] = RunOne(model, reference, particles, seeds[c], inner, basin));
        return reps;
    }

    private static Replicate RunOne(BassModel model, Reference reference, int particles, int seed, int parallelism,
                                    Func<IReadOnlyList<BassTerm>, string>? basin)
    {
        var clock = Stopwatch.StartNew();
        int[] streams = SeedTree.Derive(seed, 2);   // the run's, the composition's
        var bridge = new Bridge(model, reference);
        var target = new ParticleTarget(bridge);
        var draws = new Particles(model, reference);
        var kernel = new Kernel(bridge, basin);
        var options = new SmcOptions
        {
            Particles = particles, Resampling = Defaults.Resampling, ResampleBelow = Defaults.ResampleBelow, MutationSteps = Defaults.MutationSteps,
            MaxStages = Defaults.MaxStages, MaxDegreeOfParallelism = parallelism,
        };
        SmcRun<RefParticle> run = SmcRun<RefParticle>.Start(new SmcStage<RefParticle>(target, new Temperature(0.0, TemperingPath.Likelihood)),
                                                            draws, kernel, options, streams[0]);
        while (run.Stage.Temperature.InverseTemperature < 1.0 && run.StagesDone < Defaults.MaxStages)
            run.Advance(new SmcStage<RefParticle>(target, new Temperature(SmcSampler.NextInverseTemperature(run, Defaults.TargetConditionalEss),
                                                                          TemperingPath.Likelihood)));
        SmcResult<RefParticle> result = run.Result();
        var sigma2 = new double[particles];
        bool reached = run.Stage.Temperature.InverseTemperature == 1.0;
        int[] composition = SeedTree.Derive(streams[1], particles);
        Parallel.For(0, particles, new ParallelOptions { MaxDegreeOfParallelism = parallelism }, i =>
            sigma2[i] = reached && result.Weights[i] > 0.0 && result.Particles[i].Inside is { } p
                ? p.Moves.Compose(p.State, new Xoshiro256PlusPlus(composition[i]), new BassTermTable()).Sigma2
                : double.NaN);
        double onPilots = 0.0, level = 0.0, mass = 0.0;
        var distinct = new HashSet<string>();
        for (int i = 0; i < particles; i++)
        {
            if (result.Weights[i] == 0.0 || result.Particles[i].Inside is not { } p) continue;
            string key = BassMixing.Key(p.State.Terms);
            distinct.Add(key);
            if (reference.Empirical.ContainsKey(key)) onPilots += result.Weights[i];
            level += result.Weights[i] * BassMixing.LogPosterior(model, p.State);
            mass += result.Weights[i];
        }
        return new Replicate(seed, result, kernel.Account, kernel.Views, sigma2, draws.PriorDraws, draws.PriorInF, draws.EmpiricalDraws,
                             clock.Elapsed.TotalSeconds, onPilots, distinct.Count, level / mass);
    }

    /// <summary>The replicates' evidence: their mean, sd and range, and the pooled error, the larger of the within-run
    /// (genealogy) and between-run errors (<c>LadderEvidence.Pool</c>'s rule), with Z_F's, shared by every replicate, in
    /// quadrature.</summary>
    private readonly record struct Pooled(double Mean, double Sd, double Range, double Within, double Between, double Error);

    private static Pooled Pool(double[] values, double[] errors, double shared)
    {
        int r = values.Length;
        double mean = values.Average();
        double sd = r > 1 ? Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (r - 1)) : double.NaN;
        double within = Math.Sqrt(errors.Sum(e => e * e)) / r, between = r > 1 ? sd / Math.Sqrt(r) : 0.0;
        double run = Math.Max(within, between);
        return new Pooled(mean, sd, values.Max() - values.Min(), within, between, Math.Sqrt(run * run + shared * shared));
    }

    /// <summary>Z_F, F's share of p₀, from <see cref="BassPriorSampler"/>'s acceptance: ln of accepted over attempted, its
    /// error √((1 − Z)/accepted) (the attempts to a fixed count of acceptances).</summary>
    private readonly record struct ZF(double Log, double Error, long Attempts, long Accepted);

    private static ZF EstimateZF(BassModel model, int root, int samplers = 16, int perSampler = 2_500)
    {
        int[] seeds = SeedTree.Derive(root, samplers);
        var attempts = new long[samplers];
        var accepted = new long[samplers];
        Parallel.For(0, samplers, k =>
        {
            var sampler = new BassPriorSampler(model);
            var rng = new Xoshiro256PlusPlus(seeds[k]);
            for (int i = 0; i < perSampler; i++) sampler.DrawState(rng);
            attempts[k] = sampler.Attempts;
            accepted[k] = sampler.Accepted;
        });
        long a = attempts.Sum(), n = accepted.Sum();
        double z = (double)n / a;
        return new ZF(Math.Log(z), Math.Sqrt((1 - z) / n), a, n);
    }

    // The weighted mean of M per replicate, pooled, and its error from the replicates' spread.
    private static (double Mean, double Error) MeanCount(int[][] counts, double[][] weights)
    {
        double[] means = [.. counts.Select((m, c) => m.Select((v, i) => v * weights[c][i]).Sum() / weights[c].Sum())];
        double mean = means.Average();
        return (mean, Math.Sqrt(means.Sum(v => (v - mean) * (v - mean)) / (means.Length - 1) / means.Length));
    }

    // ---- The enumeration ----

    /// <summary>
    /// F enumerated (every structure the design and the target admit, extended from admitted ones only, F being downward
    /// closed), each with ln w and ℓ₁ on a grid of u, where p(u) holds all but e⁻⁶⁰ of its peak: the bridge's law at any b,
    /// the exact evidence and Z_F, by the trapezoid rule in u.
    /// </summary>
    private sealed class Enumerated
    {
        private const int Nodes = 1601;

        private Enumerated(BassModel model, List<Structure> states, double[] u, double[] logTau, double logDu)
        {
            Model = model;
            States = states;
            U = u;
            LogTau = logTau;
            LogDu = logDu;
        }

        public sealed record Structure(BassTerm[] Terms, string Key, double LogW, double[] LogMarginal);

        public BassModel Model { get; }
        public List<Structure> States { get; }
        public double[] U { get; }
        public double[] LogTau { get; }
        public double LogDu { get; }

        /// <summary>ln Σ_F w.</summary>
        public double LogMass => LogSumExp(States.Select(s => s.LogW));

        /// <summary>The u grid: <see cref="Nodes"/> nodes spanning where ln p(u) is within 60 nats of its peak.</summary>
        public static (double[] U, double[] LogTau, double LogDu) Grid(BassTarget target)
        {
            if (!target.SamplesTau) throw new ArgumentException("The enumeration integrates a sampled τ.", nameof(target));
            double[] coarse = [.. Enumerable.Range(0, 4001).Select(k => -50.0 + 0.025 * k)];
            double[] lp = [.. coarse.Select(target.LogTauDensity)];
            double peak = lp.Max();
            int lo = Array.FindIndex(lp, v => v >= peak - 60), hi = Array.FindLastIndex(lp, v => v >= peak - 60);
            double a = coarse[Math.Max(0, lo - 1)], z = coarse[Math.Min(coarse.Length - 1, hi + 1)];
            double[] u = [.. Enumerable.Range(0, Nodes).Select(k => a + (z - a) * k / (Nodes - 1))];
            return (u, [.. u.Select(target.LogTauDensity)], Math.Log((z - a) / (Nodes - 1)));
        }

        /// <summary>ln ∫ p(u) e^{ℓ₁(S, u)} du for a structure's fit, by the trapezoid rule on the grid.</summary>
        public static double LogIntegral(BassTarget target, BassStructureFit fit, (double[] U, double[] LogTau, double LogDu) grid)
        {
            var v = new double[Nodes];
            for (int k = 0; k < Nodes; k++)
                v[k] = grid.LogTau[k] + target.Evaluate(fit, grid.U[k]).LogMarginal + (k == 0 || k == Nodes - 1 ? -Math.Log(2.0) : 0.0);
            return LogSumExp(v) + grid.LogDu;
        }

        public static Enumerated Of(BassModel model)
        {
            BassTarget target = model.Target;
            (double[] u, double[] logTau, double logDu) = Grid(target);

            List<BassTerm> terms = BassKernelEnumeration.AllTerms(model);
            var states = new List<Structure>();
            void Visit(int from, List<BassTerm> chosen)
            {
                if (BassMixing.Loaded(model, chosen, u[0], out _) is not { } loaded) return;   // no superset is in F
                BassTerm[] held = [.. loaded.Terms];
                var marginal = new double[Nodes];
                for (int k = 0; k < Nodes; k++) marginal[k] = target.Evaluate(loaded.Fit, u[k]).LogMarginal;
                states.Add(new Structure(held, BassMixing.Key(held), target.LogStructurePrior(loaded.Fit), marginal));
                if (chosen.Count == model.MaxBasis) return;
                for (int i = from; i < terms.Count; i++)
                {
                    chosen.Add(terms[i]);
                    Visit(i + 1, chosen);
                    chosen.RemoveAt(chosen.Count - 1);
                }
            }
            Visit(0, []);
            return new Enumerated(model, states, u, logTau, logDu);
        }

        /// <summary>ln ∫ p(u) e^{b·ℓ₁(S, u)} du; with no structure, ln ∫ p(u) du.</summary>
        public double LogIntegral(Structure? s, double b)
        {
            var v = new double[Nodes];
            for (int k = 0; k < Nodes; k++)
                v[k] = LogTau[k] + (s is null ? 0.0 : b * s.LogMarginal[k]) + (k == 0 || k == Nodes - 1 ? -Math.Log(2.0) : 0.0);
            return LogSumExp(v) + LogDu;
        }

        /// <summary>E[u | S] under p(u)·e^{b·ℓ₁(S, u)}.</summary>
        public double MeanU(Structure s, double b)
        {
            double top = double.NegativeInfinity;
            for (int k = 0; k < Nodes; k++) top = Math.Max(top, LogTau[k] + b * s.LogMarginal[k]);
            double num = 0.0, den = 0.0;
            for (int k = 0; k < Nodes; k++)
            {
                double w = Math.Exp(LogTau[k] + b * s.LogMarginal[k] - top) * (k == 0 || k == Nodes - 1 ? 0.5 : 1.0);
                num += w * U[k];
                den += w;
            }
            return num / den;
        }

        /// <summary>p(u)'s mean and its quartiles, by the trapezoid rule and linear interpolation of its CDF.</summary>
        public (double Mean, double[] Quartiles) TauMoments()
        {
            var cdf = new double[Nodes];
            double du = Math.Exp(LogDu), mean = 0.0;
            for (int k = 1; k < Nodes; k++)
            {
                double left = Math.Exp(LogTau[k - 1]), right = Math.Exp(LogTau[k]);
                cdf[k] = cdf[k - 1] + 0.5 * du * (left + right);
                mean += 0.5 * du * (left * U[k - 1] + right * U[k]);
            }
            var quartiles = new double[3];
            for (int j = 0; j < 3; j++)
            {
                double p = 0.25 * (j + 1) * cdf[^1];
                int k = Array.FindIndex(cdf, c => c >= p);
                quartiles[j] = U[k - 1] + (U[k] - U[k - 1]) * (p - cdf[k - 1]) / (cdf[k] - cdf[k - 1]);
            }
            return (mean / cdf[^1], quartiles);
        }

        /// <summary>γ_b over F, normalized: q_S^{1−b}·p₀^b·∫ p(u) e^{b·ℓ₁} du.</summary>
        public double[] Bridge(Reference reference, double b)
        {
            double[] log = [.. States.Select(s =>
            {
                double logP0 = s.LogW - reference.LogC;
                return (1 - b) * reference.LogStructure(s.Terms, logP0) + b * logP0 + LogIntegral(s, b);
            })];
            double total = LogSumExp(log);
            return [.. log.Select(v => Math.Exp(v - total))];
        }
    }

    // The exact ln p(y) under the prior on F and ln Z_F, by enumeration.
    private static (double LogZ, double LogZF) Exact(Enumerated space)
    {
        double mass = space.LogMass;
        return (LogSumExp(space.States.Select(s => s.LogW + space.LogIntegral(s, 1.0))) - mass, mass - Reference.CountNormalizer(space.Model));
    }

    // ---- The pilots ----

    // A short pilot on a small fixture: one untempered chain from a prior draw, its sweeps after burn-in.
    private static List<PilotDraw> ChainPilot(BassModel model, int seed, int burn = 300, int samples = 500)
    {
        int[] seeds = SeedTree.Derive(seed, 3);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(seeds[0])).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
        for (int s = 0; s < burn; s++) chain.Step();
        var draws = new List<PilotDraw>();
        for (int s = 0; s < samples; s++) draws.Add(new PilotDraw("pilot", [.. chain.Step().Terms], 1.0));
        return draws;
    }

    // A fit's draws as pilot draws: an MCMC draw's mass 1, an SMC particle's N·R·w, so each replicate's sum is N.
    private static IEnumerable<PilotDraw> FromFit(string source, BassFit fit)
    {
        double[]? w = fit.DrawWeights;
        int i = 0;
        foreach (DrawStore<BassDraw> store in fit.Draws)
            foreach (BassDraw d in store)
            {
                double mass = w is null ? 1.0 : w[i] * w.Length;
                i++;
                yield return new PilotDraw(source, [.. d.Structure.TermIds.Select(id => fit.Terms[id])], mass);
            }
    }

    private sealed record FriedmanPilots(BassMixing.Fixture Fixture, BassModel Model, Legend Legend, List<PilotDraw> Draws, ZF ZF,
                                         BassFit Smc1000, BassFit Smc4000, double M4Seconds, double TemperedSeconds, double Smc4000Seconds,
                                         double Smc1000Seconds, double ZFSeconds);

    private static readonly Lazy<FriedmanPilots> Pilots = new(BuildPilots, LazyThreadSafetyMode.ExecutionAndPublication);

    // S10's Friedman runs at master seed 53: M4's four untempered chains (burned, then 2,000 sweeps each), S10's tempered fit
    // (its cold rungs), and fitting SMC from the prior at 4,000 particles (R2's setting) and at 1,000 (for the time); Z_F.
    private static FriedmanPilots BuildPilots()
    {
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        var draws = new List<PilotDraw>();

        var clock = Stopwatch.StartNew();
        BassChain[] chains = BassMixing.Burned(model);
        var basins = new string[chains.Length][];
        var levels = new double[chains.Length][];
        var states = new BassTerm[chains.Length][][];
        Parallel.For(0, chains.Length, c =>
        {
            basins[c] = new string[BassMixing.Samples];
            levels[c] = new double[BassMixing.Samples];
            states[c] = new BassTerm[BassMixing.Samples][];
            for (int s = 0; s < BassMixing.Samples; s++)
            {
                BassState state = chains[c].Step();
                states[c][s] = [.. state.Terms];
                basins[c][s] = fixture.Basin(state.Terms);
                levels[c][s] = BassMixing.LogPosterior(model, state);
            }
        });
        double m4 = clock.Elapsed.TotalSeconds;
        foreach (BassTerm[][] chain in states)
            foreach (BassTerm[] terms in chain) draws.Add(new PilotDraw("M4", terms, 1.0));

        clock.Restart();
        BassMixing.Fixture tempered = BassMixing.Get(BassMixing.Friedman, new TemperingSpec());
        draws.AddRange(FromFit("S10", Bass.Fit(tempered.Config, tempered.Data)));
        double temperedSeconds = clock.Elapsed.TotalSeconds;

        BassConfig Smc(int particles) => fixture.Config with { Sampling = fixture.Config.Sampling with { Smc = new SmcSpec { Particles = particles } } };
        clock.Restart();
        BassFit smc4000 = Bass.Fit(Smc(4000), fixture.Data);
        double smc4000Seconds = clock.Elapsed.TotalSeconds;
        draws.AddRange(FromFit("SMC", smc4000));
        clock.Restart();
        BassFit smc1000 = Bass.Fit(Smc(1000), fixture.Data);
        double smc1000Seconds = clock.Elapsed.TotalSeconds;

        clock.Restart();
        ZF zf = EstimateZF(model, SeedTree.Derive(53, 10)[9]);
        double zfSeconds = clock.Elapsed.TotalSeconds;
        return new FriedmanPilots(fixture, model, new Legend(basins, levels), draws, zf, smc1000, smc4000, m4, temperedSeconds,
                                  smc4000Seconds, smc1000Seconds, zfSeconds);
    }

    /// <summary>The Friedman measurement's legend: the basins M4's chains hold in 1% or more of some chain's samples,
    /// labelled A, B, … from the best chain's (highest median log posterior) by share; then a1 and a2, the four-x1·x2
    /// basins its 16,000-particle replicates found.</summary>
    private sealed class Legend
    {
        private readonly Dictionary<string, string> _labels = [];
        private readonly string[][] _basins;
        private readonly double[][] _levels;

        public Legend(string[][] basins, double[][] levels)
        {
            _basins = basins;
            _levels = levels;
            Medians = [.. levels.Select(h => BassMixing.Quantile(h, 0.5))];
            Best = Array.IndexOf(Medians, Medians.Max());
            IEnumerable<int> order = [Best, .. Enumerable.Range(0, basins.Length).Where(c => c != Best)];
            foreach (int c in order)
                foreach (var g in basins[c].GroupBy(x => x).OrderByDescending(g => g.Count()).ThenBy(g => g.Key, StringComparer.Ordinal))
                    if (g.Count() >= 0.01 * basins[c].Length && !_labels.ContainsKey(g.Key)) _labels.Add(g.Key, Letter(_labels.Count));
            _labels.TryAdd(A1, "a1");
            _labels.TryAdd(A2, "a2");
        }

        public double[] Medians { get; }
        public int Best { get; }

        public string Name(string basin) => _labels.TryGetValue(basin, out string? label) ? label : basin;

        public string Order(string basin) => _labels.TryGetValue(basin, out string? label) ? $"0{label.Length:D2}{label}" : $"1{basin}";

        private static string Letter(int k) => k < 26 ? $"{(char)('A' + k)}" : $"{(char)('A' + k / 26 - 1)}{(char)('A' + k % 26)}";

        public void Print(ITestOutputHelper output)
        {
            output.WriteLine($"### The basins M4's chains hold (chain {Best} the best, median log posterior {Medians[Best]:F1}; the others' " +
                             $"{string.Join(", ", Enumerable.Range(0, Medians.Length).Where(c => c != Best).Select(c => $"chain {c} {Medians[c]:F1}"))})");
            output.WriteLine($"label | basin | {string.Join(" | ", Enumerable.Range(0, _basins.Length).Select(c => $"share in chain {c}"))} | median log post where held");
            output.WriteLine("--- | ---" + string.Concat(Enumerable.Repeat(" | ---", _basins.Length + 1)));
            foreach ((string basin, string label) in _labels)
            {
                double[] held = [.. _basins.SelectMany((h, c) => h.Select((b, s) => (b, s, c))).Where(x => x.b == basin).Select(x => _levels[x.c][x.s])];
                output.WriteLine($"{label} | {basin} | " + string.Join(" | ", _basins.Select(h => $"{(double)h.Count(x => x == basin) / h.Length:P1}")) +
                                 $" | {(held.Length > 0 ? BassMixing.Quantile(held, 0.5).ToString("F1") : "—")}");
            }
        }
    }

    // ---- The oracle's fixtures B and C, as BassTargetParityTests builds them ----

    private static BassConfig OracleConfig(bool categorical)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = categorical ? 1 : 2, Support = new SupportRule(Minimum: 3), MaxBasis = 2 },
            Intercept = GPriorCentre.Mean,
            Sigma2 = Sigma2Prior.ScaledToResponse(1e-3),
            Coefficients = BassCoefficientPrior.G(),
        };

    private static double Hinge(double x, double knot, int sign) => Math.Max(sign * (x - knot), 0.0);

    private static double[][] FixtureB()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        var rows = new List<double[]>();
        for (int rep = 0; rep < 2; rep++)
            for (int j = 0; j < 4; j++)
                for (int i = 0; i < 4; i++)
                {
                    double x1 = i / 3.0, x2 = j / 3.0;
                    rows.Add([1 + 2 * Hinge(x1, 1 / 3.0, 1) + 3 * Hinge(x1, 1 / 3.0, 1) * Hinge(x2, 2 / 3.0, -1) + 0.25 * BassQualification.Gaussian(rng), x1, x2]);
                }
        return [.. rows];
    }

    private static double[][] FixtureC()
    {
        var rng = new Xoshiro256PlusPlus(seed: 21);
        var rows = new List<double[]>();
        for (int z = 0; z < 3; z++)
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 3; i++)
                {
                    double x1 = i * 0.5, x2 = j * 0.5;
                    rows.Add([1 + 1.5 * x1 + 1.5 * Hinge(x2, 0.5, 1) * (z == 0 ? 1 : 0) + 0.4 * BassQualification.Gaussian(rng), x1, x2, z]);
                }
        return [.. rows];
    }

    private static BassData OracleData(double[][] rows, bool categorical)
    {
        var predictors = new List<BassPredictor>
        {
            new ContinuousPredictor("x1", [.. rows.Select(r => r[1])]),
            new ContinuousPredictor("x2", [.. rows.Select(r => r[2])]),
        };
        if (categorical) predictors.Add(new CategoricalPredictor("z", [.. rows.Select(r => (int)r[3])], ["0", "1", "2"]));
        return new BassData(predictors, new ScalarResponse([.. rows.Select(r => r[0])]));
    }

    // ---- Helpers ----

    /// <summary>Batch means over chains: each item's per-batch sums, read as the mean and its standard error.</summary>
    private sealed class Batches(int items, int chains, int batches, int length)
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

    private static Dictionary<string, double> Fractions(IEnumerable<(string Basin, double Weight)> items)
    {
        var f = new Dictionary<string, double>();
        double total = 0.0;
        foreach ((string basin, double w) in items)
        {
            f[basin] = f.GetValueOrDefault(basin) + w;
            total += w;
        }
        foreach (string k in f.Keys.ToArray()) f[k] /= total;
        return f;
    }

    // The basin's count of x1·x2 terms.
    private static int Pairs(string basin)
        => basin.Split(" + ").Sum(t => t == "x1·x2" ? 1 : t.EndsWith(" × x1·x2", StringComparison.Ordinal) ? int.Parse(t[..t.IndexOf(' ')]) : 0);

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        if (double.IsNegativeInfinity(max)) return max;
        return max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }
}
