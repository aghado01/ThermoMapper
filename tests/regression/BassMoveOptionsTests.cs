using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's Phase 1 options as configuration and in the diagnostics (the local-move options brief, Steps 1 and 5).
/// <list type="bullet">
/// <item>The local steps and the fixed-M factor moves on <see cref="BassSampling"/>, the split rule on
/// <see cref="TemperingSpec"/>: each off by default, its record's defaults the decided values; each validates by name,
/// resolves on the data and round-trips through JSON.</item>
/// <item>The moves by kind (<see cref="BassMoveTally"/>): the kinds of each palette move add up to the chain's own count of
/// its attempts and acceptances, chain by chain and untempered or at a ladder's b = 1; refusals by kind match the
/// reasons' tally; a duplicate term and a term short of support are refused as their kind; a null move leaves the state,
/// the design and the tallies consistent; and a resumed run keeps them.</item>
/// </list>
/// </summary>
public sealed class BassMoveOptionsTests
{
    private readonly ITestOutputHelper _out;
    public BassMoveOptionsTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Options_AreOffByDefault_AndTheirRecordsCarryTheDecidedValues()
    {
        var config = new BassConfig();
        Assert.Null(config.Sampling.LocalSteps);
        Assert.Null(config.Sampling.FactorMoves);
        Assert.Null(new TemperingSpec().Split);
        Assert.Equal((0.5, 0.02), (new LocalSteps().Probability, new LocalSteps().WindowFraction));
        Assert.Equal((1.0, 1.0, 1.0), (new FactorMoves().Change, new FactorMoves().PruneOrGrow, new FactorMoves().VariableChange));
        Assert.Equal(0.9, new SplitRule().Threshold);

        BassResolved resolved = BassResolved.Resolve(config, BassSweepTests.Data());
        Assert.Null(resolved.StepWindows);
        Assert.Null(resolved.FactorShares);
        Assert.Null(resolved.SplitThreshold);
        Assert.DoesNotContain(resolved.AutoPicks, p => p.Setting == "LocalSteps windows");
    }

    [Theory]
    [InlineData(1, 1)]
    [InlineData(3, 1)]
    [InlineData(20, 1)]
    [InlineData(50, 1)]
    [InlineData(51, 2)]
    [InlineData(60, 2)]
    [InlineData(150, 3)]   // 0.02·150 is 3.0000000000000004 in binary; decimal keeps it 3
    [InlineData(500, 10)]
    public void Window_IsTheCeilingOfTheFraction_AtLeastOne(int values, int expected)
        => Assert.Equal(expected, new LocalSteps().Window(values));

    [Fact]
    public void Options_RefuseOutOfRangeFields_ByName()
    {
        string Reason(BassSampling sampling) => Assert.Throws<ArgumentException>(new BassConfig { Sampling = sampling }.Validate).Message;
        Assert.StartsWith("Sampling.LocalSteps.Probability", Reason(new BassSampling { LocalSteps = new LocalSteps(0) }));
        Assert.StartsWith("Sampling.LocalSteps.Probability", Reason(new BassSampling { LocalSteps = new LocalSteps(1.5) }));
        Assert.StartsWith("Sampling.LocalSteps.WindowFraction", Reason(new BassSampling { LocalSteps = new LocalSteps(0.5, 0) }));
        Assert.StartsWith("Sampling.LocalSteps.WindowFraction", Reason(new BassSampling { LocalSteps = new LocalSteps(0.5, double.NaN) }));
        Assert.StartsWith("Sampling.FactorMoves.Change", Reason(new BassSampling { FactorMoves = new FactorMoves(Change: -1) }));
        Assert.StartsWith("Sampling.FactorMoves.PruneOrGrow", Reason(new BassSampling { FactorMoves = new FactorMoves(PruneOrGrow: double.PositiveInfinity) }));
        Assert.StartsWith("Sampling.FactorMoves.VariableChange", Reason(new BassSampling { FactorMoves = new FactorMoves(VariableChange: double.NaN) }));
        Assert.StartsWith("Sampling.FactorMoves", Reason(new BassSampling { FactorMoves = new FactorMoves(1, 0, 0) }));
        Assert.StartsWith("Sampling.FactorMoves.Change",
                          Reason(new BassSampling { FactorMoves = new FactorMoves(Change: 0), LocalSteps = new LocalSteps() }));
        Assert.StartsWith("Sampling.Tempering.Split.Threshold", Reason(new BassSampling { Tempering = new TemperingSpec { Split = new SplitRule(1) } }));
        Assert.StartsWith("Sampling.Tempering.Split.Threshold", Reason(new BassSampling { Tempering = new TemperingSpec { Split = new SplitRule(0) } }));
        Assert.StartsWith("Sampling.Tempering.Split", Reason(new BassSampling { Tempering = new TemperingSpec { Levels = 6, Split = new SplitRule() } }));

        // A family without the change keeps the prune or grow and the variable change, and needs no local steps.
        new BassConfig { Sampling = new BassSampling { FactorMoves = new FactorMoves(Change: 0) } }.Validate();
    }

    [Fact]
    public void Options_ResolveOnTheData()
    {
        BassConfig config = BassSweepTests.Config() with
        {
            Sampling = new BassSampling
            {
                LocalSteps = new LocalSteps(0.5, 0.5), FactorMoves = new FactorMoves(2, 1, 1),
                Tempering = new TemperingSpec { Split = new SplitRule(0.8) },
            },
        };
        BassResolved resolved = BassResolved.Resolve(config, BassSweepTests.Data());   // three knots per variable
        Assert.Equal([2, 2], resolved.StepWindows!.Continuous);
        Assert.Empty(resolved.StepWindows.Functional);
        Assert.Equal(new BassFactorShares(0.5, 0.125, 0.125, 0.25), resolved.FactorShares);
        Assert.Equal(0.8, resolved.SplitThreshold);
        Assert.Contains(resolved.AutoPicks, p => p.Setting == "LocalSteps windows" && p.Value == "continuous [2, 2]");

        // The functional fixture's output variable has 20 locations; the windows follow each variable's own knots.
        BassMixing.Fixture functional = BassMixing.Get(BassMixing.Functional);
        BassResolved f = BassResolved.Resolve(functional.Config with { Sampling = functional.Config.Sampling with { LocalSteps = new LocalSteps() } },
                                              functional.Data);
        var model = BassModel.Create(functional.Config, functional.Data);
        Assert.Equal(Enumerable.Range(0, 2).Select(v => Math.Max(1, (int)Math.Ceiling(0.02m * model.Inputs.ContinuousKnots(v).Length))),
                     f.StepWindows!.Continuous);
        Assert.Equal([1], f.StepWindows.Functional);
        Assert.Null(f.FactorShares);
        Assert.Null(f.SplitThreshold);
    }

    [Fact]
    public void Options_RoundTripThroughJson()
    {
        var config = new BassConfig
        {
            Sampling = new BassSampling
            {
                LocalSteps = new LocalSteps(0.25, 0.05), FactorMoves = new FactorMoves(1, 2, 0.5),
                Tempering = new TemperingSpec { Split = new SplitRule(0.85) },
            },
        };
        string json = config.ToJson();
        Assert.Equal(config, BassConfig.FromJson(json));
        Assert.Equal(json, BassConfig.FromJson(json).ToJson());
        Assert.Contains("\"LocalSteps\":{\"Probability\":0.25,\"WindowFraction\":0.05}", json);
        Assert.Contains("\"FactorMoves\":{\"Change\":1,\"PruneOrGrow\":2,\"VariableChange\":0.5}", json);
        Assert.Contains("\"Split\":{\"Threshold\":0.85}", json);

        // Off reads back off; a record left empty takes its decided values.
        Assert.Equal(new BassConfig(), BassConfig.FromJson(new BassConfig().ToJson()));
        BassConfig partial = BassConfig.FromJson("{\"Sampling\":{\"LocalSteps\":{},\"FactorMoves\":{}}}");
        Assert.Equal(new LocalSteps(), partial.Sampling.LocalSteps);
        Assert.Equal(new FactorMoves(), partial.Sampling.FactorMoves);

        // The fingerprint moves with each option.
        BassData data = BassSweepTests.Data();
        string baseline = Bass.Fingerprint(new BassConfig(), data);
        Assert.NotEqual(baseline, Bass.Fingerprint(new BassConfig { Sampling = new BassSampling { LocalSteps = new LocalSteps() } }, data));
        Assert.NotEqual(baseline, Bass.Fingerprint(new BassConfig { Sampling = new BassSampling { FactorMoves = new FactorMoves() } }, data));
    }

    // Every palette move's kinds against the chain's own count of it: the redraw, the local step and the toggle make the
    // change; a birth, a death, a prune, a grow and a variable change are their own.
    private static void AssertKindsAddUp(IReadOnlyList<MoveStat> moves, BassMoveTally kinds)
    {
        BassMoveKind[] For(string key) => key switch
        {
            "birth" => [BassMoveKind.Birth],
            "death" => [BassMoveKind.Death],
            "change" => [BassMoveKind.Redraw, BassMoveKind.LocalStep, BassMoveKind.Toggle],
            "prune" => [BassMoveKind.Prune],
            "grow" => [BassMoveKind.Grow],
            _ => [BassMoveKind.VariableChange],
        };
        foreach (MoveStat m in moves)
        {
            BassKindCount[] counts = [.. For(m.Key).Select(k => kinds[k])];
            Assert.Equal(m.Attempts, counts.Sum(c => c.Attempts));
            Assert.Equal(m.Accepted, counts.Sum(c => c.Accepted));
        }
        Assert.All(kinds.Entries, c => Assert.True(c.Rejected >= 0, $"{c.Kind}: {c.Rejected} rejected"));
        Assert.Equal(moves.Sum(m => m.Attempts), kinds.Entries.Sum(c => c.Attempts));
    }

    [Theory]
    [InlineData("continuous", "none", false)]
    [InlineData("continuous", "all", false)]
    [InlineData("categorical, R's measure", "all", false)]
    [InlineData("functional", "all", false)]
    [InlineData("continuous", "all", true)]
    [InlineData("categorical, uniform", "all", true)]
    public void Kinds_AddUpToTheChainsMoveCounts(string fixture, string options, bool tempered)
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture(fixture, "equal", BassLocalStepTests.OptionsFor(options));
        config = config with
        {
            Coefficients = BassCoefficientPrior.G(),
            Sampling = config.Sampling with
            {
                Chains = 2, MasterSeed = 9, BurnIn = 400, Samples = 600,
                Tempering = tempered ? new TemperingSpec { Split = new SplitRule() } : null,
            },
        };
        BassFit fit = Bass.Fit(config, data);
        BassDiagnostics d = fit.Diagnostics;
        foreach (BassChainAccount chain in d.Chains) AssertKindsAddUp(chain.Moves, chain.MoveKinds);
        AssertKindsAddUp(d.Moves, d.MoveKinds);
        _out.WriteLine($"{fixture}, {options}, tempered {tempered}: " + string.Join("; ", d.MoveKinds.Entries.Where(c => c.Attempts > 0)
            .Select(c => $"{c.Kind} {c.Accepted}/{c.Attempts} (null {c.Null}, refused {c.Refused})")));
        if (options == "none")
        {
            Assert.All(new[] { BassMoveKind.LocalStep, BassMoveKind.Toggle, BassMoveKind.Prune, BassMoveKind.Grow, BassMoveKind.VariableChange },
                       k => Assert.Equal(0, d.MoveKinds[k].Attempts));
            Assert.All(d.MoveKinds.Entries, c => Assert.Equal(0, c.Null));
            return;
        }
        Assert.True(d.MoveKinds[BassMoveKind.LocalStep].Attempts > 0);
        Assert.True(d.MoveKinds[BassMoveKind.Prune].Attempts > 0 && d.MoveKinds[BassMoveKind.Grow].Attempts > 0 && d.MoveKinds[BassMoveKind.VariableChange].Attempts > 0);
        if (fixture.StartsWith("categorical", StringComparison.Ordinal)) Assert.True(d.MoveKinds[BassMoveKind.Toggle].Attempts > 0);
        Assert.True(d.MoveKinds.Entries.Sum(c => c.Null) > 0);
    }

    [Fact]
    public void Refusals_ByKind_MatchTheReasons()
    {
        // Without prunes, whose refusal of a term's only factor has no design reason, every refusal is by a reason.
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("categorical, R's measure", "equal",
                                                                        s => s with { LocalSteps = new LocalSteps(), FactorMoves = new FactorMoves(1, 0, 1) });
        var model = BassModel.Create(config, data);
        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(71), new Xoshiro256PlusPlus(72));
        for (int i = 0; i < 20_000; i++) chain.Step();
        BassStructureMoves moves = chain.Structure;
        Assert.Equal(moves.Refusals.Total, moves.Kinds.Entries.Sum(c => c.Refused));
        Assert.True(moves.Refusals.Count(BassRefusal.Support) > 0 && moves.Refusals.Count(BassRefusal.Numerical) > 0);
        AssertKindsAddUp(chain.Moves.MoveStats(), moves.Kinds);

        // With them, a prune of a term's only factor is refused as a prune, beyond the reasons.
        (config, data) = BassLocalStepTests.Fixture("continuous", "equal", s => s with { FactorMoves = new FactorMoves() });
        model = BassModel.Create(config, data);
        chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(73), new Xoshiro256PlusPlus(74));
        for (int i = 0; i < 20_000; i++) chain.Step();
        moves = chain.Structure;
        Assert.True(moves.Kinds[BassMoveKind.Prune].Refused > 0);
        Assert.True(moves.Kinds.Entries.Sum(c => c.Refused) > moves.Refusals.Total);
        AssertKindsAddUp(chain.Moves.MoveStats(), moves.Kinds);
    }

    [Fact]
    public void DuplicateAndSupport_AreRefusedAsTheirKind_TheDesignUntouched()
    {
        // x0+ at knot 1 beside x1+ at knot 1 on the 2,263-state fixture: changing x0 to x1 can draw the held term, a
        // duplicate; a step of x0+ to knot 2, the largest value, leaves no support.
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal",
                                                                        s => s with { LocalSteps = new LocalSteps(1.0, 0.02), FactorMoves = new FactorMoves() });
        var model = BassModel.Create(config, data);
        BassTerm a = new([new HingeFactor(0, 1, 1)]), b = new([new HingeFactor(1, 1, 1)]);
        IRjMove<BassState> Move(BassStructureMoves moves, string key) => moves.Moves.Single(m => m.Key == key);
        (BassRefusal Reason, BassMoveKind Kind, string Key)[] cases =
            [(BassRefusal.Numerical, BassMoveKind.VariableChange, "variable"), (BassRefusal.Support, BassMoveKind.LocalStep, "change")];
        foreach ((BassRefusal reason, BassMoveKind kind, string key) in cases)
        {
            bool seen = false;
            for (int seed = 0; seed < 400 && !seen; seed++)
            {
                var moves = new BassStructureMoves(model, [a, b], 0);
                BassState start = moves.Start;
                BassTerm[] before = [.. moves.DesignAt(start).Terms];
                Proposal<BassState> p = Move(moves, key).Propose(start, new Xoshiro256PlusPlus(seed))!.Value;
                if (p.Candidate.Fit.Refusal != reason) continue;
                seen = true;
                Assert.Equal(1, moves.Refusals.Count(reason));
                Assert.Equal((1L, 1L, 0L), (moves.Kinds[kind].Attempts, moves.Kinds[kind].Refused, moves.Kinds[kind].Accepted));
                Assert.True(double.IsNegativeInfinity(model.Target.Evaluate(p.Candidate).LogPrior));
                Assert.Equal(before, moves.DesignAt(start).Terms);
                Assert.False(moves.Holds(p.Candidate));
            }
            Assert.True(seen, $"no {kind} was refused for {reason}");
        }
    }

    [Fact]
    public void NullMoves_LeaveTheStateTheDesignAndTheTalliesConsistent()
    {
        // x0·x1 alone on the 2,263-state fixture: its part is full, so a variable change, its own reverse, is a null move.
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s with { FactorMoves = new FactorMoves() });
        var model = BassModel.Create(config, data);
        var full = new BassTerm([new HingeFactor(0, 1, 0), new HingeFactor(1, -1, 2)]);
        var moves = new BassStructureMoves(model, [full], 0);
        BassState start = moves.Start;
        IRjMove<BassState>[] nullMoves = [.. moves.Moves.Where(m => m.Key is "variable")];
        var chain = new ReversibleJumpChain<BassState>(nullMoves, model.Target, start, new Xoshiro256PlusPlus(5));
        for (int i = 0; i < 100; i++) Assert.Same(start, chain.Step());
        Assert.Equal((100L, 0L), (chain.Attempts, chain.Accepted));
        Assert.Equal((100L, 100L), (moves.Kinds[BassMoveKind.VariableChange].Attempts, moves.Kinds[BassMoveKind.VariableChange].Null));
        Assert.Equal(0, moves.Kinds.Entries.Sum(c => c.Accepted + c.Refused));
        Assert.Equal(0, moves.Refusals.Total);
        Assert.True(moves.Holds(start));
        Assert.Equal([full], moves.DesignAt(start).Terms);
        AssertKindsAddUp(chain.MoveStats(), moves.Kinds);
    }

    [Fact]
    public void Resume_KeepsTheKinds()
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("categorical, R's measure", "equal", BassLocalStepTests.OptionsFor("all"));
        config = config with { Coefficients = BassCoefficientPrior.G(), Sampling = config.Sampling with { Chains = 2, MasterSeed = 3, BurnIn = 200, Samples = 400 } };
        BassFit whole = Bass.Fit(config, data);
        BassConfig half = config with { Sampling = config.Sampling with { Samples = 200 } };
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(whole, resumed);
        Assert.Equal(whole.Diagnostics.MoveKinds.Entries, resumed.Diagnostics.MoveKinds.Entries);
        for (int c = 0; c < 2; c++) Assert.Equal(whole.Diagnostics.Chains[c].MoveKinds.Entries, resumed.Diagnostics.Chains[c].MoveKinds.Entries);
    }
}
