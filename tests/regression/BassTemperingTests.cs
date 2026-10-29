using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Evidence;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's tempering (stage 4's S8): <see cref="TemperingSpec"/>'s contract and refusals, <see cref="BassLadder"/>'s
/// replicas, swaps and design hand-over, its rungs' laws against enumeration, the run's tempering readout, and its
/// evidence against the exact ln Z.
/// <para>The exact references enumerate F (<c>BassTermSpace.F</c>) on small fixtures at M_max = 2 and integrate τ by the
/// trapezoid on u = ln τ through the target's own ℓ_b, which S5's parity tests check against R. The
/// oracle's parity tests (<c>BassTemperingParityTests</c>) check the rungs and the evidence against R's own enumeration.</para>
/// </summary>
public sealed class BassTemperingTests
{
    private readonly ITestOutputHelper _out;
    public BassTemperingTests(ITestOutputHelper output) => _out = output;

    internal static BassConfig Config(TemperingSpec? tempering, int maxBasis = 2, int chains = 3, int burn = 300, int samples = 200, int seed = 11)
    {
        BassConfig config = BassRunTests.Config(chains, burn, samples, seed: seed, maxBasis: maxBasis);
        return config with { Sampling = config.Sampling with { Tempering = tempering } };
    }

    /// <summary>One rung's exact law on F: ln Z_b against the prior normalized on F, the count law, E[τ] and E[ln τ].</summary>
    internal sealed record RungLaw(double LogZ, double[] Count, double TauMean, double LogTauMean);

    /// <summary>The exact law of <paramref name="model"/>'s target at each of <paramref name="rungs"/>, over F and τ.</summary>
    internal static RungLaw[] Exact(BassModel model, params double[] rungs)
        => Enumerate(model, rungs, -60.0, (e, b) => e.LogPrior + e.LogLikelihood(b), withLogZ: true);

    // The rungs' laws over F and τ, with log π_b(S, u) from `logDensity`, by the trapezoid on u from `lowest` up to 8
    // (step 0.05; the integrands are smooth and decay at both ends, so the rule converges geometrically). Each state's
    // marginal and τ moments are reduced first, then the states. ln Z, against the prior normalized on F, only `withLogZ`.
    private static RungLaw[] Enumerate(BassModel model, double[] rungs, double lowest, Func<BassEvaluation, double, double> logDensity, bool withLogZ)
    {
        BassTarget target = model.Target;
        List<(BassTerm[] Terms, double LogPrior)> states = BassTermSpace.F(model);
        BassDesign design = model.NewDesign();
        const double step = 0.05;
        double[] nodes = target.SamplesTau
            ? [.. Enumerable.Range(0, (int)Math.Round((8.0 - lowest) / step) + 1).Select(i => lowest + step * i)]
            : [Math.Log(target.FixedTau!.Value)];
        double logStep = target.SamplesTau ? Math.Log(step) : 0.0;
        var marginal = new double[rungs.Length][];
        var tauGiven = new double[rungs.Length][];
        var logTauGiven = new double[rungs.Length][];
        for (int r = 0; r < rungs.Length; r++) (marginal[r], tauGiven[r], logTauGiven[r]) = (new double[states.Count], new double[states.Count], new double[states.Count]);
        var values = new double[nodes.Length];
        for (int s = 0; s < states.Count; s++)
        {
            BassTerm[] terms = states[s].Terms;
            Assert.Null(design.Load(terms));
            BassStructureFit fit = BassStructureFit.Of(design, terms.Sum(model.Prior.LogPrior));
            BassEvaluation[] e = [.. nodes.Select(u => target.SamplesTau ? target.Evaluate(fit, u) : target.Evaluate(fit))];
            for (int r = 0; r < rungs.Length; r++)
            {
                for (int k = 0; k < nodes.Length; k++) values[k] = e[k].OffSupport ? double.NegativeInfinity : logDensity(e[k], rungs[r]) + logStep;
                double m = LogSumExp(values), tau = 0.0, u = 0.0;
                for (int k = 0; k < nodes.Length; k++)
                {
                    double w = Math.Exp(values[k] - m);
                    tau += w * Math.Exp(nodes[k]);
                    u += w * nodes[k];
                }
                (marginal[r][s], tauGiven[r][s], logTauGiven[r][s]) = (m, tau, u);
            }
        }
        double logPriorMass = LogSumExp(states.Select(x => x.LogPrior));
        var laws = new RungLaw[rungs.Length];
        for (int r = 0; r < rungs.Length; r++)
        {
            double z = LogSumExp(marginal[r]);
            var count = new double[model.MaxBasis + 1];
            double tau = 0.0, logTau = 0.0;
            for (int s = 0; s < states.Count; s++)
            {
                double w = Math.Exp(marginal[r][s] - z);
                count[states[s].Terms.Length] += w;
                tau += w * tauGiven[r][s];
                logTau += w * logTauGiven[r][s];
            }
            laws[r] = new RungLaw(withLogZ ? z - logPriorMass : double.NaN, count, tau, logTau);
        }
        return laws;
    }

    internal static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        return max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }

    // ── The contract ──────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Spec_HasD05sDefaults_RoundTrips_AndMovesTheFingerprint()
    {
        var spec = new TemperingSpec();
        Assert.Equal((TemperingPath.Likelihood, TemperedForm.Integrated, (int?)null, 64, (double?)null, true, SwapSchedule.DeterministicEvenOdd, false),
                     (spec.Path, spec.Form, spec.Levels, spec.MaxLevels, spec.HotEnd, spec.Adaptive, spec.Schedule, spec.ReferenceRung));
        Assert.Equal(8, TemperingSpec.InitialLevels);
        Assert.Equal((2, 3), (TemperingSpec.MinimumLevels(false), TemperingSpec.MinimumLevels(true)));
        Assert.Null(new BassSampling().Tempering);
        Assert.Contains("\"Tempering\":{\"Path\":\"Likelihood\",\"Form\":\"Integrated\",\"Levels\":null,\"MaxLevels\":64,\"HotEnd\":null,", Config(spec).ToJson());

        BassConfig config = Config(spec with { Form = TemperedForm.Linear, ReferenceRung = true, Levels = 6, HotEnd = 0.02, Adaptive = false });
        string json = config.ToJson();
        Assert.Contains("\"Tempering\":{\"Path\":\"Likelihood\",\"Form\":\"Linear\",\"Levels\":6,\"MaxLevels\":64,\"HotEnd\":0.02,\"Adaptive\":false,\"Schedule\":\"DeterministicEvenOdd\",\"ReferenceRung\":true}", json);
        Assert.Equal(config, BassConfig.FromJson(json));
        Assert.Equal(Config(spec), BassConfig.FromJson(Config(spec).ToJson()));
        Assert.Equal(new TemperingSpec { ReferenceRung = true }, BassConfig.FromJson("{\"Sampling\":{\"Tempering\":{\"ReferenceRung\":true}}}").Sampling.Tempering);
        Assert.ThrowsAny<System.Text.Json.JsonException>(() => BassConfig.FromJson("{\"Sampling\":{\"Tempering\":{\"Rungs\":4}}}"));

        BassData data = BassSweepTests.Data();
        string baseline = Bass.Fingerprint(Config(spec), data);
        Assert.NotEqual(Bass.Fingerprint(Config(null), data), baseline);
        foreach (TemperingSpec other in new[]
                 {
                     spec with { Path = TemperingPath.Powered }, spec with { Form = TemperedForm.Linear }, spec with { Levels = 9 },
                     spec with { MaxLevels = 32 }, spec with { HotEnd = 2e-3 }, spec with { Levels = 8, Adaptive = false }, spec with { ReferenceRung = true },
                 })
            Assert.NotEqual(baseline, Bass.Fingerprint(Config(other), data));
    }

    public static TheoryData<string, TemperingSpec> Unrunnable => new()
    {
        { "Sampling.Tempering.Path", new TemperingSpec { Path = (TemperingPath)7 } },
        { "Sampling.Tempering.Form", new TemperingSpec { Form = (TemperedForm)7 } },
        { "Sampling.Tempering.Form", new TemperingSpec { Path = TemperingPath.Powered, Form = TemperedForm.Linear } },
        { "Sampling.Tempering.ReferenceRung", new TemperingSpec { Path = TemperingPath.Powered, ReferenceRung = true } },
        { "Sampling.Tempering.Levels", new TemperingSpec { Levels = 1 } },
        { "Sampling.Tempering.Levels", new TemperingSpec { Levels = 2, ReferenceRung = true } },
        { "Sampling.Tempering.HotEnd", new TemperingSpec { HotEnd = 0.0 } },
        { "Sampling.Tempering.HotEnd", new TemperingSpec { HotEnd = 1.0 } },
        { "Sampling.Tempering.HotEnd", new TemperingSpec { HotEnd = -0.1 } },
        { "Sampling.Tempering.HotEnd", new TemperingSpec { HotEnd = double.NaN } },
        { "Sampling.Tempering.Adaptive", new TemperingSpec { Adaptive = false } },
        { "Sampling.Tempering.MaxLevels", new TemperingSpec { MaxLevels = 1 } },
        { "Sampling.Tempering.MaxLevels", new TemperingSpec { MaxLevels = 2, ReferenceRung = true } },
        { "Sampling.Tempering.MaxLevels", new TemperingSpec { Levels = 10, MaxLevels = 9 } },
        { "Sampling.Tempering.Schedule", new TemperingSpec { Schedule = SwapSchedule.StochasticEvenOdd } },
        { "Sampling.Tempering.Schedule", new TemperingSpec { Schedule = SwapSchedule.Sequential } },
    };

    [Theory]
    [MemberData(nameof(Unrunnable))]
    public void Validate_RefusesALadderThatCannotRun(string field, TemperingSpec spec)
    {
        var e = Assert.Throws<ArgumentException>(() => Config(spec).Validate());
        Assert.StartsWith(field + " ", e.Message);
        Assert.Equal("config", e.ParamName);
        Assert.Throws<ArgumentException>(() => Bass.Fit(Config(spec), BassSweepTests.Data()));
        Config(new TemperingSpec { Levels = 2 }).Validate();
        Config(new TemperingSpec { Levels = 3, ReferenceRung = true }).Validate();
        Config(new TemperingSpec { Levels = 4, Adaptive = false }).Validate();
        Config(new TemperingSpec { MaxLevels = 3, ReferenceRung = true }).Validate();
        Config(new TemperingSpec { Levels = 70, MaxLevels = 70 }).Validate();
    }

    [Fact]
    public void Resolve_RefusesRungsWhereTheFormIsUndefined()
    {
        BassData data = BassSweepTests.Data();   // N = 18
        BassConfig Flat(TemperingSpec t, bool jeffreys = false) => Config(t) with
        {
            Intercept = GPriorCentre.Flat, Sigma2 = jeffreys ? Sigma2Prior.Jeffreys : Sigma2Prior.ScaledToResponse(1e-3),
        };
        double limit = (1 - 2e-3) / 18;   // the flat intercept's integrated form needs c + (Nb − 1)/2 > 0

        var below = Assert.Throws<ArgumentException>(() => BassResolved.Resolve(Flat(new TemperingSpec { HotEnd = 0.99 * limit }), data));
        Assert.StartsWith("Sampling.Tempering.HotEnd", below.Message);
        Assert.Contains("the flat intercept", below.Message);
        Assert.Equal("config", below.ParamName);
        BassResolved.Resolve(Flat(new TemperingSpec { HotEnd = 1.01 * limit }), data);
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(Flat(new TemperingSpec { HotEnd = 0.99 / 18 }, jeffreys: true), data));
        BassResolved.Resolve(Flat(new TemperingSpec { HotEnd = 1.01 / 18 }, jeffreys: true), data);
        // The linear form, and the powered path, which reads ℓ₁ alone, reach every b.
        BassResolved.Resolve(Flat(new TemperingSpec { HotEnd = 1e-6, Form = TemperedForm.Linear }), data);
        BassResolved.Resolve(Flat(new TemperingSpec { HotEnd = 1e-6, Path = TemperingPath.Powered }), data);

        // The reference rung needs b = 0, which the integrated form lacks under either improper option.
        foreach (BassConfig improper in new[]
                 {
                     Config(new TemperingSpec { ReferenceRung = true }) with { Sigma2 = Sigma2Prior.Jeffreys },
                     Flat(new TemperingSpec { ReferenceRung = true, HotEnd = 0.1 }),
                 })
        {
            var e = Assert.Throws<ArgumentException>(() => BassResolved.Resolve(improper, data));
            Assert.StartsWith("Sampling.Tempering.ReferenceRung", e.Message);
            Assert.Throws<ArgumentException>(() => Bass.Fit(improper, data));
            BassResolved.Resolve(improper with { Sampling = improper.Sampling with { Tempering = improper.Sampling.Tempering! with { Form = TemperedForm.Linear } } }, data);
        }
        BassResolved.Resolve(Config(new TemperingSpec { ReferenceRung = true }), data);
    }

    [Fact]
    public void Model_TakesTheConfigsForm()
    {
        BassData data = BassSweepTests.Data();
        Assert.Equal(TemperedForm.Integrated, BassModel.Create(Config(null), data).Target.Form);
        Assert.Equal(TemperedForm.Linear, BassModel.Create(Config(new TemperingSpec { Form = TemperedForm.Linear }), data).Target.Form);
        Assert.Equal(TemperedForm.Linear, BassModel.Create(Config(null), data, TemperedForm.Linear).Target.Form);
    }

    // ── The ladder ────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Adopt_TakesTheGiversDesignBitForBit_AndContinuesAsTheGiverWould()
    {
        BassModel model = BassModel.Create(Config(null, maxBasis: 8), BassSweepTests.Data(large: true));
        var sampler = new BassPriorSampler(model);
        BassState a = sampler.DrawState(new Xoshiro256PlusPlus(1)).State, b = sampler.DrawState(new Xoshiro256PlusPlus(2)).State;
        var giver = new BassStructureMoves(model, a.Terms, a.LogTau);
        var chain = new ReversibleJumpChain<BassState>(giver.Moves, model.Target, giver.Start, new Xoshiro256PlusPlus(3));
        for (int t = 0; t < 400 || chain.Current.Count < 2; t++)
            if (t == 20_000) Assert.Fail("The chain never held two terms.");
            else chain.Step();
        BassState x = chain.Current;
        Assert.True(giver.Holds(x));

        var receiver = new BassStructureMoves(model, b.Terms, b.LogTau);
        Assert.False(receiver.Holds(x));
        receiver.Adopt(x, giver.DesignFor(x));
        Assert.True(receiver.Holds(x));
        AssertSameDesign(giver.DesignAt(x), receiver.DesignAt(x));

        // From x on the same stream, the two continue as one.
        var left = new ReversibleJumpChain<BassState>(giver.Moves, model.Target, x, new Xoshiro256PlusPlus(4));
        var right = new ReversibleJumpChain<BassState>(receiver.Moves, model.Target, x, new Xoshiro256PlusPlus(4));
        for (int t = 0; t < 300; t++)
        {
            BassState l = left.Step(), r = right.Step();
            Assert.Equal(l.Terms, r.Terms);
            Assert.Equal(l.Fit, r.Fit);
        }
        AssertSameDesign(giver.DesignAt(left.Current), receiver.DesignAt(right.Current));

        // A design for another structure is refused before anything changes.
        Assert.Throws<ArgumentException>(() => receiver.Adopt(b, giver.DesignFor(left.Current)));
        Assert.True(receiver.Holds(right.Current));
        AssertSameDesign(giver.DesignAt(left.Current), receiver.DesignAt(right.Current));
    }

    internal static void AssertSameDesign(BassDesign a, BassDesign b)
    {
        Assert.Equal(a.Terms, b.Terms);
        Assert.Equal(a.ResidualSquares, b.ResidualSquares);
        Assert.Equal(a.FittedSquares, b.FittedSquares);
        Assert.Equal(a.Coefficients.ToArray(), b.Coefficients.ToArray());
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(false, true)]
    public void Ladder_EveryReplicasDesignFollowsItsState(bool reference, bool fixedG)
    {
        BassConfig config = Config(null, maxBasis: 8);
        if (fixedG) config = config with { Coefficients = new FixedG() };
        BassModel model = BassModel.Create(config, BassSweepTests.Data(large: true));
        var spec = new TemperingSpec { Levels = 5, HotEnd = 0.01, ReferenceRung = reference };
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(5)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 6);
        BassDesign fresh = model.NewDesign();
        long fromReference = 0, handedToHot = 0;
        for (int t = 0; t < 1500; t++)
        {
            ladder.Step();
            for (int i = 0; i < spec.Levels; i++)
            {
                // A replica holds its state's design, or the state is a reference draw, whose design no replica holds.
                BassState state = ladder.Tempering.Replicas[i].Current;
                bool held = ladder.Structures[i].Holds(state);
                if (!held)
                {
                    fromReference++;
                    Assert.True(reference && ladder.Structures.All(s => !s.Holds(state)), $"replica {i} holds a state it was not handed with its design");
                }
                else if (reference && i == spec.Levels - 1) handedToHot++;
                BassDesign design = ladder.Structures[i].DesignAt(state);
                Assert.Null(fresh.Load(state.Terms));
                double scale = fresh.ResidualSquares + fresh.FittedSquares;
                Assert.True(Math.Abs(design.ResidualSquares - fresh.ResidualSquares) <= 1e-10 * scale);
                Assert.True(Math.Abs(design.FittedSquares - fresh.FittedSquares) <= 1e-10 * scale);
                Assert.True(Math.Abs(state.Fit.ResidualSquares - fresh.ResidualSquares) <= 1e-10 * scale);
            }
        }
        _out.WriteLine($"reference {reference}, fixed g {fixedG}: {ladder.Tempering.SwapAccepts} swaps of {ladder.Tempering.SwapAttempts}, " +
                       $"{fromReference} reference draws, {handedToHot} states handed up to the reference rung");
        Assert.True(ladder.Tempering.SwapAccepts > 300);
        if (reference) Assert.True(fromReference > 0 && handedToHot > 0);
        Assert.Equal(1.0, ladder.Tempering.Betas[0]);
    }

    public static TheoryData<string, double[]> RungCases => new()
    {
        { "integrated", [1.0, 0.4, 0.0] },
        { "linear", [1.0, 0.4, 0.0] },
        { "integrated, reference", [1.0, 0.4, 0.1, 0.0] },
        { "linear, flat + jeffreys", [1.0, 0.4, 0.0] },
        { "fixed g", [1.0, 0.4, 0.0] },
        { "powered", [1.0, 0.6, 0.3] },
        { "powered, flat + jeffreys", [1.0, 0.6, 0.3] },   // the improper constants cancel from b·(prior + ℓ₁)'s ratios
        { "flat, near the limit", [1.0, 0.3, 0.06] },   // the integrated form needs b > (1 − 2c)/N = 0.0554 here
    };

    // Each rung's law, from the ladder's own replicas, against enumeration: b = 1 the posterior and b = 0 the prior on F,
    // whether by local moves or by the reference rung's exact draws, and an interior rung its tempered law, which swaps
    // that missed ℓ_b, or moves that tempered the prior, would shift.
    [Theory]
    [MemberData(nameof(RungCases))]
    public void Ladder_EachRungSamplesItsTemperedLaw(string option, double[] betas)
    {
        bool reference = option.Contains("reference"), powered = option.StartsWith("powered");
        var spec = new TemperingSpec
        {
            Form = option.Contains("linear") ? TemperedForm.Linear : TemperedForm.Integrated, Levels = betas.Length,
            HotEnd = betas[^1] > 0.0 ? betas[^1] : betas[^2],
            ReferenceRung = reference, Path = powered ? TemperingPath.Powered : TemperingPath.Likelihood,
        };
        BassConfig config = Config(spec);
        if (option.Contains("flat")) config = config with { Intercept = GPriorCentre.Flat };
        if (option.Contains("jeffreys")) config = config with { Sigma2 = Sigma2Prior.Jeffreys };
        if (option == "fixed g") config = config with { Coefficients = new FixedG() };
        BassModel model = BassModel.Create(config, BassSweepTests.Data());
        const int ladders = 4, burn = 2_000, scans = 40_000, batches = 20;
        int levels = betas.Length;
        var count = new double[levels][,];
        var logTau = new double[levels][];
        for (int r = 0; r < levels; r++)
        {
            count[r] = new double[model.MaxBasis + 1, ladders * batches];
            logTau[r] = new double[ladders * batches];
        }
        var sampler = new BassPriorSampler(model);
        for (int l = 0; l < ladders; l++)
        {
            BassState start = sampler.DrawState(new Xoshiro256PlusPlus(100 + l)).State;
            var ladder = new BassLadder(model, spec, betas, start.Terms, start.LogTau, seed: 200 + l);
            for (int t = 0; t < burn; t++) ladder.Step();
            for (int t = 0; t < scans; t++)
            {
                ladder.Step();
                int batch = l * batches + t / (scans / batches);
                for (int r = 0; r < levels; r++)
                {
                    BassState s = ladder.Tempering.Replicas[r].Current;
                    count[r][s.Count, batch] += 1.0 / (scans / batches);
                    logTau[r][batch] += s.LogTau / (scans / batches);
                }
            }
        }

        RungLaw[] exact = powered ? ExactPowered(model, betas) : Exact(model, betas);
        for (int r = 0; r < levels; r++)
        {
            for (int m = 0; m <= model.MaxBasis; m++)
            {
                (double mean, double error) = Summary([.. Enumerable.Range(0, ladders * batches).Select(k => count[r][m, k])]);
                _out.WriteLine($"{option}, b = {betas[r]:G3}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
                Assert.True(Math.Abs(mean - exact[r].Count[m]) <= 5 * error + 0.004, $"b = {betas[r]}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
            }
            if (!model.Target.SamplesTau) continue;
            // At b = 0 τ's law is its prior, Gamma(a, rate aN): E[ln τ] = ψ(a) − ln(aN), which the enumeration's rule reproduces.
            if (betas[r] == 0.0)
                Assert.Equal(SpecialFunctions.Digamma(0.5) - Math.Log(0.5 * model.Target.Observations), exact[r].LogTauMean, 1e-9);
            (double u, double uError) = Summary(logTau[r]);
            _out.WriteLine($"{option}, b = {betas[r]:G3}, E[ln τ]: {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
            Assert.True(Math.Abs(u - exact[r].LogTauMean) <= 5 * uError + 0.01, $"b = {betas[r]}, E[ln τ] {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
        }
    }

    // The powered path's rung laws, (prior · p₁)^b with τ powered in u = ln τ, by the same enumeration; its left tail in u
    // decays as e^{b(a + k/2)u}, so the rule starts lower.
    private static RungLaw[] ExactPowered(BassModel model, double[] rungs)
        => Enumerate(model, rungs, -150.0, (e, b) => b * Temperature.LogPosterior(e), withLogZ: false);

    private static (double Mean, double Error) Summary(double[] batchMeans)
    {
        double mean = batchMeans.Average();
        double variance = batchMeans.Sum(v => (v - mean) * (v - mean)) / (batchMeans.Length - 1);
        return (mean, Math.Sqrt(variance / batchMeans.Length));
    }

    [Fact]
    public void Ladder_RefusesWhatItCannotDo()
    {
        BassData data = BassSweepTests.Data();
        var spec = new TemperingSpec { Levels = 4, HotEnd = 0.05, ReferenceRung = true };
        BassModel model = BassModel.Create(Config(spec), data);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 2);

        // The model's target must be in the spec's form.
        var form = Assert.Throws<ArgumentException>(() => new BassLadder(model, spec with { Form = TemperedForm.Linear }, start.Terms, start.LogTau, 2));
        Assert.Equal("spec", form.ParamName);
        Assert.Throws<ArgumentNullException>(() => new BassLadder(model, null!, start.Terms, start.LogTau, 2));
        Assert.Throws<ArgumentNullException>(() => new BassLadder(null!, spec, start.Terms, start.LogTau, 2));
        foreach (TemperingSpec bad in new[] { spec with { Levels = 2 }, spec with { HotEnd = 0.0 }, spec with { HotEnd = 1.0 }, spec with { Path = TemperingPath.Powered } })
            Assert.Equal("spec", Assert.Throws<ArgumentException>(() => new BassLadder(model, bad, start.Terms, start.LogTau, 2)).ParamName);
        Assert.Equal("betas", Assert.Throws<ArgumentException>(() => new BassLadder(model, spec, [1.0, 0.5, 0.0], start.Terms, start.LogTau, 2)).ParamName);
        Assert.Throws<ArgumentException>(() => new BassLadder(model, spec, [0.9, 0.5, 0.1, 0.0], start.Terms, start.LogTau, 2));
        Assert.Throws<InvalidOperationException>(() => BassTemperingReadout.Of([ladder]));   // burn-in has not ended

        // Adaptation runs before burn-in ends, and only then; burn-in ends once.
        for (int t = 0; t < 20; t++) ladder.Step();
        Assert.False(ladder.Frozen);
        ladder.Adapt();
        Assert.Throws<InvalidOperationException>(() => ladder.Adapt());   // no scans since the last adaptation
        ladder.Step();
        ladder.EndBurnIn(traceRungs: true);
        Assert.True(ladder.Frozen);
        Assert.Throws<InvalidOperationException>(() => ladder.Adapt());
        Assert.Throws<InvalidOperationException>(() => ladder.EndBurnIn(traceRungs: false));

        // A snapshot restores only under the spec it ran with.
        var snapshot = ladder.Snapshot();
        Assert.Throws<ArgumentException>(() => BassLadder.Restore(model, spec with { Levels = 5 }, snapshot));
        Assert.Throws<ArgumentException>(() => BassLadder.Restore(model, spec with { ReferenceRung = false }, snapshot));
        Assert.Throws<ArgumentException>(() => BassLadder.Restore(model, spec with { Path = TemperingPath.Powered, ReferenceRung = false }, snapshot));
        Assert.Throws<ArgumentNullException>(() => BassLadder.Restore(model, spec, null!));
        Assert.Throws<ArgumentException>(() => BassLadder.Restore(model, spec, snapshot with { Structures = snapshot.Structures[..^1] }));
        Assert.True(BassLadder.Restore(model, spec, snapshot).Frozen);
    }

    // ── The run ───────────────────────────────────────────────────────────────────────────────────────────────────

    public static TheoryData<string> RunOptions => new()
    {
        "reference", "linear", "powered", "fixed ladder", "fixed g", "functional", "categorical", "two rungs", "three rungs, reference",
        "early stop", "self-sizing", "self-sizing, linear", "self-sizing, powered", "self-sizing, fixed g", "self-sizing, functional",
        "self-sizing, categorical", "self-sizing, flat",
    };

    internal static (BassConfig Config, BassData Data) RunCase(string option, int samples = 200)
    {
        var spec = new TemperingSpec { Levels = 5, HotEnd = 0.01 };
        spec = option switch
        {
            "reference" => spec with { ReferenceRung = true },
            "linear" => spec with { Form = TemperedForm.Linear, ReferenceRung = true },
            "powered" => spec with { Path = TemperingPath.Powered },
            "fixed ladder" => spec with { Adaptive = false, ReferenceRung = true },
            "flat + jeffreys" => spec with { Form = TemperedForm.Linear, ReferenceRung = true },
            "thinned" => spec with { ReferenceRung = true },
            "two rungs" => spec with { Levels = 2, HotEnd = 0.2 },
            "three rungs, reference" => spec with { Levels = 3, HotEnd = 0.2, ReferenceRung = true },
            "self-sizing, linear" => new TemperingSpec { Form = TemperedForm.Linear, ReferenceRung = true },
            "self-sizing, powered" => new TemperingSpec { Path = TemperingPath.Powered },
            "self-sizing, flat" => new TemperingSpec(),
            _ when option.StartsWith("self-sizing", StringComparison.Ordinal) => new TemperingSpec { ReferenceRung = true },
            _ => spec,
        };
        BassConfig config = Config(spec, maxBasis: 6, samples: samples);
        config = option switch
        {
            "fixed g" or "self-sizing, fixed g" => config with { Coefficients = new FixedG() },
            "flat + jeffreys" => config with { Intercept = GPriorCentre.Flat, Sigma2 = Sigma2Prior.Jeffreys },
            "self-sizing, flat" => config with { Intercept = GPriorCentre.Flat },
            "thinned" => config with { Sampling = config.Sampling with { Thin = 2 } },
            "early stop" => config with { Sampling = config.Sampling with { RHatTarget = 1.5, BatchSize = 50, MinEssPerChain = 5 } },
            _ => config,
        };
        BassData data = option switch
        {
            "functional" or "self-sizing, functional" => BassRunTests.FunctionalData(),
            "categorical" or "self-sizing, categorical" => BassPriorSamplerTests.MixedData(),
            _ => BassSweepTests.Data(large: true),
        };
        return (config, data);
    }

    [Theory]
    [MemberData(nameof(RunOptions))]
    public void Fit_ReportsTheLaddersReadout(string option)
    {
        (BassConfig config, BassData data) = RunCase(option);
        TemperingSpec spec = config.Sampling.Tempering!;
        BassFit fit = Bass.Fit(config, data);
        BassTemperingReadout readout = fit.Diagnostics.Tempering!;
        Assert.Equal((spec.Path, spec.Schedule, spec.ReferenceRung), (readout.Path, readout.Schedule, readout.ReferenceRung));
        Assert.Equal(spec.Path == TemperingPath.Likelihood ? spec.Form : TemperedForm.Integrated, readout.Form);
        int chains = config.Sampling.Chains;
        Assert.Equal(chains, readout.Ladders.Length);
        BassResolved resolved = fit.Resolved;
        double hot = spec.HotEnd ?? resolved.HotEnd!.Value;
        Assert.Equal(hot, resolved.HotEnd);
        int start = spec.Levels ?? TemperingSpec.InitialLevels;
        double[] declared = spec.ReferenceRung
            ? ParallelTempering<BassState>.ReferenceLadder(start, hot)
            : ParallelTempering<BassState>.GeometricLadder(start, hot);
        Assert.Equal(chains, readout.Levels.Length);
        Assert.Equal((chains, chains), (readout.PairRejection.Length, readout.PairAcceptance.Length));
        double worst = 0.0, ideal = 0.0;
        for (int c = 0; c < chains; c++)
        {
            double[] ladder = readout.Ladders[c];
            Assert.Equal(readout.Levels[c], ladder.Length);
            if (spec.Levels is int levels) Assert.Equal(levels, ladder.Length);
            else Assert.InRange(ladder.Length, TemperingSpec.MinimumLevels(spec.ReferenceRung), spec.MaxLevels);
            Assert.Equal(1.0, ladder[0]);
            Assert.Equal(declared[^1], ladder[^1]);
            for (int i = 1; i < ladder.Length; i++) Assert.True(ladder[i] < ladder[i - 1]);
            if (!spec.Adaptive) Assert.Equal(declared, ladder);
            Assert.Equal((ladder.Length - 1, ladder.Length - 1), (readout.PairRejection[c].Length, readout.PairAcceptance[c].Length));
            for (int i = 0; i < ladder.Length - 1; i++)
            {
                double rejection = readout.PairRejection[c][i], acceptance = readout.PairAcceptance[c][i];
                Assert.InRange(rejection, 0.0, 1.0);
                Assert.InRange(acceptance, 0.0, 1.0);
                // Rao–Blackwellized and realized acceptance read the same rate.
                Assert.True(Math.Abs(acceptance - (1 - rejection)) <= 0.2, $"chain {c}, pair {i}: acceptance {acceptance:F3}, rejection {rejection:F3}");
                worst = Math.Max(worst, rejection);
            }
            ideal += TemperingLadder.IdealRoundTripRate(readout.PairRejection[c], spec.Schedule);
        }
        Assert.Equal(readout.PairRejection.Average(r => r.Sum()), readout.Barrier, 12);
        Assert.Equal(worst, readout.MaxPairRejection);
        Assert.Equal(ideal / chains, readout.IdealRoundTripRate, 15);
        Assert.False(readout.AtMaxLevels);
        Assert.True(readout.RoundTripsPerScan >= 0 && readout.AnnealedRestartsPerScan >= 0);
        if (spec.ReferenceRung) Assert.True(readout.ReferenceAttemptsPerDraw >= 1.0);
        else Assert.True(double.IsNaN(readout.ReferenceAttemptsPerDraw));
        Assert.All(fit.Draws, d => Assert.Equal(fit.Diagnostics.SamplesUsed, d.Count));
        Assert.Equal(option == "early stop", fit.Diagnostics.StoppedEarly);
        _out.WriteLine($"{option}: round trips {readout.RoundTripsPerScan:F4}/scan (ideal {readout.IdealRoundTripRate:F4}), barrier {readout.Barrier:F2}, " +
                       $"levels [{string.Join(", ", readout.Levels)}], acceptance [{string.Join("; ", readout.PairAcceptance.Select(c => string.Join(", ", c.Select(a => a.ToString("F2")))))}], " +
                       $"evidence {fit.Evidence?.LogEvidence.ToString("F3") ?? "null"}");
        Assert.Null(Bass.Fit(config with { Sampling = config.Sampling with { Tempering = null } }, data).Diagnostics.Tempering);
    }

    [Fact]
    public void Fit_CountsALaddersScanAsOneSweep()
    {
        (BassConfig config, BassData data) = RunCase("reference", samples: 40);
        config = config with { Sampling = config.Sampling with { BurnIn = 30, Thin = 2 } };
        var reports = new List<BassProgress>();
        Bass.Fit(config, data, new BassRunOptions { Progress = new Inline(p => { lock (reports) reports.Add(p); }), ProgressInterval = 7 });
        long planned = config.Sampling.Chains * (30 + 40L * 2);
        Assert.Equal(new BassProgress(planned, planned), reports.MaxBy(r => r.Sweeps));
    }

    private sealed class Inline(Action<BassProgress> report) : IProgress<BassProgress>
    {
        public void Report(BassProgress value) => report(value);
    }

    // ── The evidence ──────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    [InlineData("zero centre")]
    [InlineData("fixed g")]
    [InlineData("declared σ²")]
    public void Evidence_MatchesTheExactLogZ(string option)
    {
        var spec = new TemperingSpec { ReferenceRung = true, Form = option == "linear" ? TemperedForm.Linear : TemperedForm.Integrated };
        BassConfig config = Config(spec, chains: 4, burn: 1_000, samples: 4_000);
        config = option switch
        {
            "zero centre" => config with { Intercept = GPriorCentre.Zero },
            "fixed g" => config with { Coefficients = new FixedG() },
            "declared σ²" => config with { Sigma2 = Sigma2Prior.InverseGamma(2.0, 0.5) },
            _ => config,
        };
        BassData data = BassSweepTests.Data();
        double exact = Exact(BassModel.Create(config, data), 1.0)[0].LogZ;
        BassFit fit = Bass.Fit(config, data);
        BassEvidence evidence = fit.Evidence!;
        _out.WriteLine($"{option}: stepping-stone {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4}, exact {exact:F4}; per chain " +
                       string.Join(", ", evidence.SteppingStone.Chains.Select(c => $"{c.LogRatio:F3} ± {c.MonteCarloError:F3}")));
        Assert.Equal(LadderEstimator.SteppingStone, evidence.SteppingStone.Method);
        Assert.Equal(evidence.SteppingStone.LogEvidence, evidence.LogEvidence);
        Assert.Equal(evidence.SteppingStone.MonteCarloError, evidence.MonteCarloError);
        Assert.Equal(4, evidence.SteppingStone.Chains.Count);
        Assert.Equal(evidence.SteppingStone.Chains.Average(c => c.LogRatio), evidence.LogEvidence, 12);
        // The reported error is the larger of the chains' own errors pooled and the standard error of their mean.
        double within = Math.Sqrt(evidence.SteppingStone.Chains.Sum(c => c.MonteCarloError * c.MonteCarloError)) / 4;
        double between = Math.Sqrt(evidence.SteppingStone.Chains.Sum(c => Math.Pow(c.LogRatio - evidence.LogEvidence, 2)) / 12);
        Assert.Equal(within, evidence.SteppingStone.WithinChainError, 12);
        Assert.Equal(between, evidence.SteppingStone.BetweenChainError, 12);
        Assert.Equal(Math.Max(evidence.SteppingStone.WithinChainError, evidence.SteppingStone.BetweenChainError), evidence.MonteCarloError);
        // The default ladder sizes itself, from the hot end it auto-picks: 1/N, and under fixed g the held τ, 1/N too.
        Assert.Equal(1.0 / 18, fit.Resolved.HotEnd);
        _out.WriteLine($"{option}: N = {fit.Resolved.Observations}, settled levels [{string.Join(", ", fit.Diagnostics.Tempering!.Levels)}], cold chain's mean M " +
                       $"{fit.Draws[0].Average(d => d.Structure.TermIds.Length):F3}, max rejection {evidence.MaxPairRejection:F3}, " +
                       $"between/within {evidence.SteppingStone.BetweenToWithin:F2}");
        Assert.True(double.IsNaN(evidence.SteppingStone.DiscretizationError));
        Assert.True(Math.Abs(evidence.LogEvidence - exact) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact:F4}");

        if (spec.Form == TemperedForm.Linear)
        {
            BassEvidenceEstimate ps = evidence.PathSampling!;
            _out.WriteLine($"{option}: path sampling {ps.LogEvidence:F4} ± {ps.MonteCarloError:F4}, diagnostic {ps.DiscretizationError:F4}");
            Assert.Equal(LadderEstimator.PathSampling, ps.Method);
            // The diagnostic reads the trapezoid's error on the linear form, where ℓ_b is linear in b.
            Assert.True(Math.Abs(ps.LogEvidence + ps.DiscretizationError - exact) <= 4 * ps.MonteCarloError + 0.5 * Math.Abs(ps.DiscretizationError) + 0.01,
                        $"{ps.LogEvidence:F4} + {ps.DiscretizationError:F4} against {exact:F4}");
        }
        else Assert.Null(evidence.PathSampling);

        BassResolved r = fit.Resolved;
        BassEvidencePriors p = evidence.Priors;
        Assert.Equal(BassEvidencePriors.Of(config, r), p);
        Assert.Equal((config.Coefficients, r.TauPrior, r.FixedTau, config.Count, r.MaxBasis), (p.Coefficients, p.TauPrior, p.FixedTau, p.Count, p.MaxBasis));
        Assert.Equal((config.Structure.CategoricalSubsets, r.MaxInteraction, r.MaxCategoricalInteraction, r.MaxFunctionalInteraction),
                     (p.CategoricalSubsets, p.MaxInteraction, p.MaxCategoricalInteraction, p.MaxFunctionalInteraction));
        Assert.Equal((config.Structure.Support, r.SupportMinimum, r.SupportMinimumFunctional, config.Structure.AngleFloor, config.Structure.MaxVarianceInflation),
                     (p.Support, p.SupportMinimum, p.SupportMinimumFunctional, p.AngleFloor, p.MaxVarianceInflation));
        Assert.Equal((config.Intercept, r.CentreLevel, config.Sigma2, r.Sigma2!.Value), (p.Centre, p.CentreLevel, p.Sigma2Prior, p.Sigma2));
    }

    // A direct functional response small enough to enumerate: x on {0, ½, 1} three times (n = 9), t on four locations
    // (N = 36), Y = 1 + 2x + 1.5·t·x + 2(x − ½)₊(t − ½)₊ + noise.
    internal static BassData SmallFunctionalData()
    {
        var rng = new Xoshiro256PlusPlus(17);
        double[] x = [.. Enumerable.Range(0, 9).Select(i => (i % 3) / 2.0)];
        double[] t = [0.0, 1 / 3.0, 2 / 3.0, 1.0];
        var y = new double[x.Length, t.Length];
        for (int i = 0; i < x.Length; i++)
            for (int r = 0; r < t.Length; r++)
            {
                double noise = Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
                y[i, r] = 1 + 2 * x[i] + 1.5 * t[r] * x[i] + 2 * Math.Max(x[i] - 0.5, 0) * Math.Max(t[r] - 0.5, 0) + 0.3 * noise;
            }
        return new BassData([new ContinuousPredictor("x", x)], new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    public void Evidence_MatchesTheExactLogZ_OnADirectFunctionalResponse(string form)
    {
        BassData data = SmallFunctionalData();
        var spec = new TemperingSpec { Form = form == "linear" ? TemperedForm.Linear : TemperedForm.Integrated, ReferenceRung = true };
        BassConfig config = Config(spec, chains: 4, burn: 1_000, samples: 4_000);
        BassModel model = BassModel.Create(config, data);
        double exact = Exact(model, 1.0)[0].LogZ;
        BassEvidence evidence = Bass.Fit(config, data).Evidence!;
        _out.WriteLine($"{form}: stepping-stone {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4}, exact {exact:F4} over " +
                       $"{BassTermSpace.F(model).Count} structures");
        Assert.True(Math.Abs(evidence.LogEvidence - exact) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact:F4}");
        Assert.Equal(1, evidence.Priors.MaxFunctionalInteraction);
        Assert.NotNull(evidence.Priors.SupportMinimumFunctional);
    }

    [Theory]
    [InlineData("jeffreys")]
    [InlineData("flat")]
    [InlineData("not asked")]
    [InlineData("no reference rung")]
    [InlineData("powered")]
    [InlineData("untempered")]
    public void Evidence_IsNullWhereUndefinedOrUnasked(string option)
    {
        var spec = new TemperingSpec { Levels = 4, HotEnd = 0.1, ReferenceRung = true, Form = TemperedForm.Linear };
        BassConfig config = Config(spec, samples: 50, burn: 50);
        config = option switch
        {
            "jeffreys" => config with { Sigma2 = Sigma2Prior.Jeffreys },
            "flat" => config with { Intercept = GPriorCentre.Flat },
            "not asked" => config with { Sampling = config.Sampling with { ReportEvidence = false } },
            "no reference rung" => config with { Sampling = config.Sampling with { Tempering = spec with { ReferenceRung = false } } },
            "powered" => config with { Sampling = config.Sampling with { Tempering = spec with { ReferenceRung = false, Path = TemperingPath.Powered, Form = TemperedForm.Integrated } } },
            _ => config with { Sampling = config.Sampling with { Tempering = null } },
        };
        BassFit fit = Bass.Fit(config, BassSweepTests.Data());
        Assert.Null(fit.Evidence);
        Assert.Equal(option == "untempered", fit.Diagnostics.Tempering is null);
        if (option is "jeffreys" or "flat")
            Assert.Throws<ArgumentException>(() => BassEvidencePriors.Of(config with { Sigma2 = Sigma2Prior.Jeffreys }, fit.Resolved with { Sigma2 = null }));
    }
}
