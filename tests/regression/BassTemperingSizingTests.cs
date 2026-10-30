using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Evidence;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's self-sizing ladder (stage 4, chip E): <see cref="TemperingSpec.Levels"/> null sizes each chain's ladder during
/// burn-in by its barrier within <see cref="TemperingSpec.MaxLevels"/>, and <see cref="TemperingSpec.HotEnd"/> null
/// auto-picks the hot end on the data (<see cref="BassResolved.HotEnd"/>).
/// <list type="bullet">
/// <item>A declared count and hot end run bit for bit as S8's did.</item>
/// <item>On the N = 175 functional fixture, where S8's eight rungs saturated, no pair rejects as S8's did, the two forms
/// agree, the chains agree, and a fixed 32-rung ladder at equal work agrees.</item>
/// <item>The hot end's auto-pick under each centre, the flat intercept's pole and fixed g, logged; every refusal.</item>
/// <item>After resizes every replica holds its state's design and no two share a mid-run term array; a cap that binds is
/// flagged and the evidence still reported; short burn-ins size with the rounds they have; the slopes are traced on the
/// linear form only.</item>
/// </list>
/// Resume across a resize is in <c>BassTemperingResumeTests</c>, and fixture B's evidence against the oracle's exact
/// ln Z, self-sized on both forms, in <c>BassTemperingParityTests</c>.
/// </summary>
public sealed class BassTemperingSizingTests
{
    private readonly ITestOutputHelper _out;
    public BassTemperingSizingTests(ITestOutputHelper output) => _out = output;

    // Recorded at 2bdc22a, before the ladder sized itself: a declared count and hot end run exactly as S8's did.
    [Theory]
    [InlineData("reference", "c0d3fc64b28d76b8")]
    [InlineData("linear", "3bf1c3e139267020")]
    [InlineData("powered", "48cd2d9e1a42b59c")]
    [InlineData("fixed ladder", "473c10b9a2c3a814")]
    [InlineData("fixed g", "7a5ad501887e8497")]
    [InlineData("functional", "5bf78e576a989948")]
    [InlineData("categorical", "5306a0a8fe37b061")]
    public void FixedLevels_RunAsS8Did(string option, string pinned)
    {
        (BassConfig config, BassData data) = BassTemperingTests.RunCase(option, samples: 120);
        Assert.Equal(pinned, Pin(Bass.Fit(config, data)));
    }

    // The functional fixture as S8 ran it (M_max 3), with the reference rung, four chains of 4,000 burn-in and 4,000 draws.
    private static BassConfig FunctionalConfig(TemperingSpec spec) => BassTemperingTests.Config(spec, maxBasis: 3, chains: 4, burn: 4_000, samples: 4_000);

    [Fact]
    public void SelfSizing_OnTheFunctionalFixture_NoPairRejectsAsS8s_AndTheFormsAndChainsAgree()
    {
        BassData data = BassRunTests.FunctionalData();
        var fits = new Dictionary<string, BassFit>();
        foreach ((string name, TemperingSpec spec) in new[]
                 {
                     ("integrated", new TemperingSpec { ReferenceRung = true }),
                     ("linear", new TemperingSpec { ReferenceRung = true, Form = TemperedForm.Linear }),
                     ("integrated, 32 fixed rungs", new TemperingSpec { ReferenceRung = true, Levels = 32, HotEnd = 1e-3 }),
                 })
        {
            BassFit fit = Bass.Fit(FunctionalConfig(spec), data);
            fits[name] = fit;
            BassTemperingReadout r = fit.Diagnostics.Tempering!;
            BassEvidence e = fit.Evidence!;
            double meanM = fit.Draws[0].Average(d => d.Structure.TermIds.Length);
            _out.WriteLine($"{name}: N = {fit.Resolved.Observations}, hot end {fit.Resolved.HotEnd:G4}, settled levels [{string.Join(", ", r.Levels)}], " +
                           $"cold chain's mean M {meanM:F3}, barrier {r.Barrier:F2}, max rejection {r.MaxPairRejection:F3}, round trips {r.RoundTripsPerScan:F4}/scan");
            _out.WriteLine($"  evidence {e.LogEvidence:F3} ± {e.MonteCarloError:F3} (within {e.SteppingStone!.WithinChainError:F3}, between " +
                           $"{e.SteppingStone.BetweenChainError:F3}); chains {string.Join(", ", e.SteppingStone.Chains.Select(c => $"{c.LogRatio:F2} ± {c.MonteCarloError:F2}"))}" +
                           (e.PathSampling is { } ps ? $"; path sampling {ps.LogEvidence:F3} ± {ps.MonteCarloError:F3}" : ""));
            // S8's eight rungs rejected 0.96–1.00 here, its chains spreading 122–183 against errors of 2–5.
            Assert.True(r.MaxPairRejection < 0.95, $"{name}: a pair rejected {r.MaxPairRejection:F3}");
            Assert.Equal(r.MaxPairRejection, e.MaxPairRejection);
            Assert.Equal(e.MaxPairRejection > BassEvidence.RejectionThreshold, e.Saturated);
            Assert.True(e.SteppingStone.BetweenToWithin < 2.0, $"{name}: between/within {e.SteppingStone.BetweenToWithin:F2}");
            Assert.False(e.AtMaxLevels);
            if (spec.Levels is null) Assert.All(r.Levels, l => Assert.InRange(l, TemperingSpec.InitialLevels + 1, 32));
        }
        void Agree(string a, string b)
        {
            BassEvidence x = fits[a].Evidence!, y = fits[b].Evidence!;
            double error = Math.Sqrt(x.MonteCarloError * x.MonteCarloError + y.MonteCarloError * y.MonteCarloError);
            Assert.True(Math.Abs(x.LogEvidence - y.LogEvidence) <= 3 * error, $"{a} {x.LogEvidence:F3} against {b} {y.LogEvidence:F3} (± {error:F3})");
        }
        Agree("integrated", "linear");
        Agree("integrated", "integrated, 32 fixed rungs");
        Agree("linear", "integrated, 32 fixed rungs");
    }

    [Fact]
    public void HotEnd_AutoPicksOnTheData_AndIsLogged()
    {
        BassData data = BassSweepTests.Data();   // N = 18
        double Picked(BassConfig config, out BassAutoPick? logged)
        {
            BassResolved resolved = BassResolved.Resolve(config, data);
            logged = resolved.AutoPicks.SingleOrDefault(p => p.Setting == "Tempering.HotEnd");
            return resolved.HotEnd!.Value;
        }
        BassConfig With(TemperingSpec spec) => BassTemperingTests.Config(spec);

        // 1/N, where E[τ] = 1/N halves the shrinkage, under each centre and on each form.
        foreach (BassConfig config in new[]
                 {
                     With(new TemperingSpec()), With(new TemperingSpec { Form = TemperedForm.Linear }), With(new TemperingSpec { ReferenceRung = true }),
                     With(new TemperingSpec()) with { Intercept = GPriorCentre.Zero },
                     With(new TemperingSpec { Form = TemperedForm.Linear }) with { Intercept = GPriorCentre.Flat },
                     With(new TemperingSpec { Path = TemperingPath.Powered }) with { Intercept = GPriorCentre.Flat },
                     With(new TemperingSpec()) with { Sigma2 = Sigma2Prior.Jeffreys },
                 })
        {
            Assert.Equal(1.0 / 18, Picked(config, out BassAutoPick? logged));
            Assert.Equal(("Tempering.HotEnd", (1.0 / 18).ToString("R")), (logged!.Setting, logged.Value));
            Assert.StartsWith("1/N", logged.Reason);
        }
        // 2/N on the integrated form under the flat intercept, past its pole at (1 − 2c)/N.
        BassConfig flat = With(new TemperingSpec()) with { Intercept = GPriorCentre.Flat };
        Assert.Equal(2.0 / 18, Picked(flat, out BassAutoPick? pole));
        Assert.StartsWith("2/N", pole!.Reason);
        Assert.Equal(2.0 / 18, Picked(flat with { Sigma2 = Sigma2Prior.Jeffreys }, out _));
        // Under a held τ the smaller of those and the b at which it halves a column's shrinkage: τ under the g-form,
        // τ/N under the ridge.
        Assert.Equal(1.0 / 18, Picked(With(new TemperingSpec()) with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) }, out _));
        Assert.Equal(1.0 / 18, Picked(With(new TemperingSpec()) with { Coefficients = BassCoefficientPrior.G(TauPrior.Held(1.0)) }, out BassAutoPick? above));
        Assert.EndsWith("below the b at which the held τ halves that shrinkage", above!.Reason);
        Assert.Equal(0.01, Picked(With(new TemperingSpec()) with { Coefficients = BassCoefficientPrior.G(TauPrior.Held(0.01)) }, out BassAutoPick? held));
        Assert.StartsWith("the b at which the held τ halves", held!.Reason);
        Assert.Equal(1.0 / 18, Picked(With(new TemperingSpec()) with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.Held(1.0)) }, out _));
        Assert.Equal(0.09 / 18, Picked(With(new TemperingSpec()) with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.Held(0.09)) }, out _));
        Assert.Equal(1.0 / 18, Picked(With(new TemperingSpec()) with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) }, out _));
        Assert.Equal(0.08, Picked(flat with { Coefficients = BassCoefficientPrior.G(TauPrior.Held(0.08)) }, out _));   // the held τ, below 2/N and past the pole at 0.0554

        // A declared hot end is no auto-pick; a self-sizing ladder's start is logged; an untempered run has none.
        BassResolved declared = BassResolved.Resolve(With(new TemperingSpec { HotEnd = 0.02, Levels = 5 }), data);
        Assert.Equal(0.02, declared.HotEnd);
        Assert.DoesNotContain(declared.AutoPicks, p => p.Setting is "Tempering.HotEnd" or "Tempering.Levels");
        Assert.Contains(BassResolved.Resolve(With(new TemperingSpec()), data).AutoPicks, p => p.Setting == "Tempering.Levels" && p.Value == "8, self-sizing within 64");
        Assert.Contains(BassResolved.Resolve(With(new TemperingSpec { MaxLevels = 5 }), data).AutoPicks, p => p.Setting == "Tempering.Levels" && p.Value == "5, self-sizing within 5");
        Assert.Equal((5, 3, 6), (new TemperingSpec { MaxLevels = 5 }.StartLevels(), new TemperingSpec { MaxLevels = 3, ReferenceRung = true }.StartLevels(),
                                 new TemperingSpec { Levels = 6 }.StartLevels()));
        Assert.Null(BassResolved.Resolve(With(null!), data).HotEnd);

        // The ladder's own pick, on a model resolved without it, is the same.
        BassModel model = BassModel.Create(With(null!) with { Coefficients = BassCoefficientPrior.G(TauPrior.Held(0.01)) }, data);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        Assert.Equal(0.01, new BassLadder(model, new TemperingSpec(), start.Terms, start.LogTau, 2).Tempering.InverseTemperatures[^1]);
        Assert.Equal(0.01, new BassLadder(model, new TemperingSpec { ReferenceRung = true }, start.Terms, start.LogTau, 2).Tempering.InverseTemperatures[^2]);
    }

    [Fact]
    public void HotEnd_AnAutoPickTheFormCannotTake_IsRefused()
    {
        // Tiny data: two cases, so the flat intercept's 2/N is 1, outside (0, 1).
        var tiny = new BassData([new ContinuousPredictor("x", [0.0, 1.0])], new ScalarResponse([0.3, 1.1]));
        BassConfig flat = BassTemperingTests.Config(new TemperingSpec()) with { Intercept = GPriorCentre.Flat };
        var e = Assert.Throws<ArgumentException>(() => BassResolved.Resolve(flat, tiny));
        Assert.StartsWith("Sampling.Tempering.HotEnd auto-picks 1 on N = 2", e.Message);
        Assert.Contains("outside (0, 1); set it.", e.Message);
        Assert.Equal("config", e.ParamName);
        BassResolved.Resolve(flat with { Sampling = flat.Sampling with { Tempering = new TemperingSpec { HotEnd = 0.9 } } }, tiny);
        // A ladder built on a model resolved without tempering makes the same pick, and refuses it in its spec's name.
        BassModel untempered = BassModel.Create(flat with { Sampling = flat.Sampling with { Tempering = null } }, tiny);
        BassState start = new BassPriorSampler(untempered).DrawState(new Xoshiro256PlusPlus(1)).State;
        Assert.Equal("spec", Assert.Throws<ArgumentException>(() => new BassLadder(untempered, new TemperingSpec(), start.Terms, start.LogTau, 2)).ParamName);
        Assert.Equal(1.0, new BassLadder(untempered, new TemperingSpec { HotEnd = 0.5 }, start.Terms, start.LogTau, 2).Tempering.InverseTemperatures[0]);

        // Under fixed g a held τ below the flat intercept's pole takes the hot end there, which the integrated form lacks.
        BassData data = BassSweepTests.Data();
        BassConfig small = BassTemperingTests.Config(new TemperingSpec()) with { Intercept = GPriorCentre.Flat, Coefficients = BassCoefficientPrior.G(TauPrior.Held(1e-3)) };
        var pole = Assert.Throws<ArgumentException>(() => BassResolved.Resolve(small, data));
        Assert.StartsWith("Sampling.Tempering.HotEnd auto-picks 0.001", pole.Message);
        Assert.Contains("set it above", pole.Message);
        BassResolved.Resolve(small with { Sampling = small.Sampling with { Tempering = new TemperingSpec { Form = TemperedForm.Linear } } }, data);
    }

    [Fact]
    public void Ladder_RefusesWhatItCannotSize()
    {
        BassData data = BassSweepTests.Data();
        BassModel model = BassModel.Create(BassTemperingTests.Config(null), data);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var fixedLadder = new BassLadder(model, new TemperingSpec { Levels = 4, HotEnd = 0.05 }, start.Terms, start.LogTau, 2);
        fixedLadder.Step();
        Assert.Contains("fixed", Assert.Throws<InvalidOperationException>(() => fixedLadder.Resize()).Message);

        var spec = new TemperingSpec { ReferenceRung = true, MaxLevels = 12 };
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, 2);
        Assert.Equal(8, ladder.Tempering.Replicas.Count);
        Assert.Contains("No scans", Assert.Throws<InvalidOperationException>(() => ladder.Resize()).Message);
        ladder.Step();
        ladder.Resize();
        Assert.Throws<InvalidOperationException>(() => ladder.Resize());   // no scans since
        ladder.Step();
        ladder.EndBurnIn(traceRungs: false);
        Assert.Contains("frozen", Assert.Throws<InvalidOperationException>(() => ladder.Resize()).Message);

        // The spec's own refusals, in its name; a count outside the budget; a snapshot of a count the spec does not run.
        foreach (TemperingSpec bad in new[] { spec with { Adaptive = false }, spec with { MaxLevels = 2 }, new TemperingSpec { Levels = 6, MaxLevels = 5 } })
            Assert.Equal("spec", Assert.Throws<ArgumentException>(() => new BassLadder(model, bad, start.Terms, start.LogTau, 2)).ParamName);
        Assert.Equal("inverseTemperatures", Assert.Throws<ArgumentException>(() =>
            new BassLadder(model, spec with { MaxLevels = 4 }, [1.0, 0.5, 0.2, 0.1, 0.0], start.Terms, start.LogTau, 2)).ParamName);
        Assert.Equal("inverseTemperatures", Assert.Throws<ArgumentException>(() => new BassLadder(model, spec, [1.0, 0.0], start.Terms, start.LogTau, 2)).ParamName);
        var snapshot = ladder.Snapshot();
        Assert.Throws<ArgumentException>(() => BassLadder.Restore(model, spec with { MaxLevels = ladder.Tempering.Replicas.Count - 1 }, snapshot));
        Assert.Throws<ArgumentException>(() => BassLadder.Restore(model, spec with { Levels = ladder.Tempering.Replicas.Count + 1 }, snapshot));
        Assert.Equal(ladder.Tempering.InverseTemperatures, BassLadder.Restore(model, spec, snapshot).Tempering.InverseTemperatures);
        Assert.Equal(ladder.Tempering.InverseTemperatures, BassLadder.Restore(model, spec with { Levels = ladder.Tempering.Replicas.Count }, snapshot).Tempering.InverseTemperatures);
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void AfterResizes_EveryReplicaHoldsItsStatesDesign_AndNoTwoShareAMidRunArray(bool reference, bool fixedG)
    {
        BassConfig config = BassTemperingTests.Config(null, maxBasis: 8);
        if (fixedG) config = config with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) };
        BassModel model = BassModel.Create(config, BassSweepTests.Data(large: true));
        var spec = new TemperingSpec { ReferenceRung = reference, MaxLevels = 24 };
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(5)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 6);
        BassDesign fresh = model.NewDesign();
        var counts = new List<int>();
        int resized = 0;
        for (int round = 1; round <= 256; round *= 2)
        {
            for (int t = 0; t < round; t++) ladder.Step();
            int before = ladder.Tempering.Replicas.Count;
            ladder.Resize();
            resized += ladder.Tempering.Replicas.Count != before ? 1 : 0;
            counts.Add(ladder.Tempering.Replicas.Count);
            IReadOnlyList<ReversibleJumpChain<BassState>> replicas = ladder.Tempering.Replicas;
            Assert.Equal(replicas.Count, ladder.Structures.Count);
            for (int i = 0; i < replicas.Count; i++)
            {
                BassState state = replicas[i].Current;
                // Each holds its own state's design, or the state is a reference draw, whose design no replica holds.
                Assert.True(ladder.Structures[i].Holds(state) || reference && ladder.Structures.All(s => !s.Holds(state)),
                            $"round {round}: replica {i} holds a state without its design");
                BassDesign design = ladder.Structures[i].DesignAt(state);
                Assert.Null(fresh.Load(state.Terms));
                double scale = fresh.ResidualSquares + fresh.FittedSquares;
                Assert.True(Math.Abs(design.ResidualSquares - fresh.ResidualSquares) <= 1e-10 * scale);
                for (int j = 0; j < i; j++)
                    Assert.True(!replicas[j].Current.SharesStructure(state) || state.SharesStructure(ladder.Structures[0].Start),
                                $"round {round}: replicas {j} and {i} share a mid-run term array");
            }
            // The rungs stay a ladder on the spec's ends.
            double[] inverseTemperatures = ladder.Tempering.InverseTemperatures;
            Assert.Equal(1.0, inverseTemperatures[0]);
            if (reference) Assert.Equal(0.0, inverseTemperatures[^1]);
            else Assert.Equal(inverseTemperatures[^1], model.Resolved.HotEnd ?? 1.0 / model.Resolved.Observations);
        }
        _out.WriteLine($"reference {reference}, fixed g {fixedG}: counts [{string.Join(", ", counts)}], {resized} resizes, " +
                       $"{ladder.Tempering.SwapAccepts} swaps of {ladder.Tempering.SwapAttempts}");
        Assert.True(resized >= 1);
        // The ladder runs on as one: the swaps still hand designs over.
        for (int t = 0; t < 500; t++) ladder.Step();
        for (int i = 0; i < ladder.Tempering.Replicas.Count; i++)
        {
            BassState state = ladder.Tempering.Replicas[i].Current;
            Assert.True(ladder.Structures[i].Holds(state) || reference && ladder.Structures.All(s => !s.Holds(state)));
        }
    }

    [Fact]
    public void Resize_ForksTheReferenceRungsState_WithItsDesign()
    {
        // Three rungs whose reference pair carries most of the barrier: a grown ladder's new rung near b = 0 is nearest the
        // reference rung on the barrier, and starts from its exact draw, whose design its replica holds only once loaded.
        BassModel model = BassModel.Create(BassTemperingTests.Config(null, maxBasis: 8), BassSweepTests.Data(large: true));
        var spec = new TemperingSpec { ReferenceRung = true, MaxLevels = 24 };
        BassDesign fresh = model.NewDesign();
        int forks = 0;
        for (int seed = 0; seed < 10 && forks == 0; seed++)
        {
            BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(30 + seed)).State;
            var ladder = new BassLadder(model, spec, [1.0, 0.5, 0.0], start.Terms, start.LogTau, seed: 40 + seed);
            for (int t = 0; t < 4; t++) ladder.Step();
            BassState hot = ladder.Tempering.Replicas[^1].Current;
            ladder.Resize();
            for (int k = 1; k < ladder.Tempering.Replicas.Count - 1; k++)
            {
                BassState state = ladder.Tempering.Replicas[k].Current;
                if (!(state.LogTau == hot.LogTau && state.Terms.SequenceEqual(hot.Terms))) continue;
                // A fork: the draw's terms and τ on an array of its own, with a design it holds, as a whole load gives it.
                forks++;
                Assert.False(state.SharesStructure(hot));
                Assert.True(ladder.Structures[k].Holds(state));
                BassDesign design = ladder.Structures[k].DesignAt(state);
                Assert.Null(fresh.Load(state.Terms));
                BassTemperingTests.AssertSameDesign(fresh, design);
            }
            for (int t = 0; t < 50; t++) ladder.Step();   // and the ladder runs on
        }
        Assert.True(forks > 0, "no resize forked the reference rung's state");
    }

    [Fact]
    public void ACapThatBinds_IsFlagged_AndTheEvidenceStillReported()
    {
        var spec = new TemperingSpec { ReferenceRung = true, MaxLevels = 4 };
        BassFit fit = Bass.Fit(BassTemperingTests.Config(spec, maxBasis: 3, chains: 3, burn: 600, samples: 400), BassRunTests.FunctionalData());
        BassTemperingReadout r = fit.Diagnostics.Tempering!;
        Assert.All(r.Levels, l => Assert.Equal(4, l));
        Assert.True(r.AtMaxLevels);
        BassEvidence e = fit.Evidence!;
        Assert.True(e.AtMaxLevels);
        Assert.Equal(r.MaxPairRejection, e.MaxPairRejection);
        Assert.True(e.Saturated, $"max rejection {e.MaxPairRejection:F3}");
        Assert.True(double.IsFinite(e.LogEvidence) && e.MonteCarloError > 0);
        _out.WriteLine($"cap 4: max rejection {e.MaxPairRejection:F3}, evidence {e.LogEvidence:F2} ± {e.MonteCarloError:F2}, between/within {e.SteppingStone!.BetweenToWithin:F2}");

        // A fixed ladder at its budget is no binding cap; an unsaturated one is not flagged.
        BassFit fixedFit = Bass.Fit(BassTemperingTests.Config(new TemperingSpec { ReferenceRung = true, Levels = 6, HotEnd = 0.05, MaxLevels = 6 }, samples: 200),
                                    BassSweepTests.Data());
        Assert.False(fixedFit.Diagnostics.Tempering!.AtMaxLevels);
        Assert.False(fixedFit.Evidence!.AtMaxLevels);
        Assert.False(fixedFit.Evidence.Saturated, $"max rejection {fixedFit.Evidence.MaxPairRejection:F3}");
    }

    [Theory]
    [InlineData(0, 0)]
    [InlineData(1, 0)]
    [InlineData(3, 1)]
    [InlineData(7, 2)]
    [InlineData(100, 5)]
    public void ShortBurnIns_SizeWithTheRoundsTheyHave(int burn, int sizingRounds)
    {
        // Rounds of 1, 2, 4, … fit the burn-in, and every round but the last sizes; so 0 and 1 sweeps keep the start's count.
        var spec = new TemperingSpec { ReferenceRung = true };
        BassConfig config = BassTemperingTests.Config(spec, chains: 2, burn: burn, samples: 60);
        BassData data = BassRunTests.FunctionalData();
        BassFit fit = Bass.Fit(config, data);
        int[] levels = fit.Diagnostics.Tempering!.Levels;
        _out.WriteLine($"burn-in {burn}: levels [{string.Join(", ", levels)}]");
        if (burn <= 1) Assert.All(levels, l => Assert.Equal(TemperingSpec.InitialLevels, l));
        Assert.NotNull(fit.Evidence);

        // Chain 0's burn-in replayed by hand, from its own seeds: the run sized and equalized in the same rounds, and the
        // last round kept the count the one before it settled.
        BassModel model = BassModel.Create(config, data);
        int[] seeds = SeedTree.Derive(SeedTree.Derive(config.Sampling.MasterSeed, 2 * config.Sampling.Chains)[0], 4);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(seeds[0])).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seeds[3]);
        int done = 0, sized = 0;
        for (int round = 1; done + round <= burn; round *= 2)
        {
            for (int t = 0; t < round; t++) ladder.Step();
            done += round;
            int settled = ladder.Tempering.Replicas.Count;
            if (done + 2 * round <= burn)
            {
                ladder.Resize();
                sized++;
            }
            else
            {
                ladder.Adapt();
                Assert.Equal(settled, ladder.Tempering.Replicas.Count);
            }
        }
        Assert.Equal(sizingRounds, sized);
        Assert.Equal(fit.Diagnostics.Tempering.Ladders[0], ladder.Tempering.InverseTemperatures);
    }

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    public void EndBurnIn_TracesSlopesOnTheLinearFormOnly(string form)
    {
        var spec = new TemperingSpec { ReferenceRung = true, Form = form == "linear" ? TemperedForm.Linear : TemperedForm.Integrated };
        BassModel model = BassModel.Create(BassTemperingTests.Config(spec), BassSweepTests.Data());
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(3)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, 4);
        for (int t = 0; t < 20; t++) ladder.Step();
        ladder.EndBurnIn(traceRungs: true);
        for (int t = 0; t < 20; t++) ladder.Step();
        Assert.Equal(form == "linear", ladder.Tempering.TracesSlopes);
        Assert.All(ladder.Tempering.RungTraces!, trace => Assert.Equal(form == "linear", trace.Draws(RungFunctional.Slope).All(double.IsFinite)));
        Assert.All(ladder.Tempering.RungTraces!.Skip(1), trace => Assert.All(trace.Draws(RungFunctional.LogRatio), v => Assert.True(double.IsFinite(v))));
    }

    // A hash of everything a fixed ladder's run produces bit for bit: every chain's draws, each chain's frozen rungs and
    // the round-trip rate, and each chain's evidence estimate.
    internal static string Pin(BassFit fit)
    {
        ulong h = 14695981039346656037UL;
        void Add(double v) => Word(BitConverter.DoubleToInt64Bits(v));
        void Word(long w)
        {
            unchecked
            {
                h ^= (ulong)w;
                h *= 1099511628211UL;
            }
        }
        foreach (DrawStore<BassDraw> chain in fit.Draws)
            foreach (BassDraw d in chain)
            {
                foreach (int id in d.Structure.TermIds) Word(id);
                Add(d.Tau);
                Add(d.Sigma2);
                foreach (double a in d.Coefficients) Add(a);
                Add(d.Lambda);
            }
        BassTemperingReadout readout = fit.Diagnostics.Tempering!;
        foreach (double[] ladder in readout.Ladders)
            foreach (double b in ladder) Add(b);
        Add(readout.RoundTripsPerScan);
        if (fit.Evidence is { } e)
            foreach (BassEvidenceEstimate estimate in new[] { e.SteppingStone, e.PathSampling }.OfType<BassEvidenceEstimate>())
                foreach (LadderEstimate c in estimate.Chains)
                {
                    Add(c.LogRatio);
                    Add(c.MonteCarloError);
                }
        return h.ToString("x16");
    }
}
