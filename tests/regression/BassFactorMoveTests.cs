using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The fixed-M factor moves (the local-move options brief, Step 3): a factor's prune and grow within a term's part, and a
/// factor's variable change, sharing the change's slot.
/// <list type="bullet">
/// <item>The ratios as built, from the proposals both ways, against the target's prior ratio: they cancel where the
/// design draft derives it (K ≥ 1 left in the part, and the variable change) and, by the uniform order law on
/// {0, …, cap}, at a part's K = 1 → 0, over every term and move of each fixture.</item>
/// <item>The production draws against the statement; null moves and refusals at the boundaries; the palette's weights.</item>
/// <item>Detailed balance by enumeration (<see cref="BassKernelEnumeration"/>) on the three fixtures, the family alone and
/// with the local steps, on both schedules, at b = 1 and 0.3, both subset measures; with a negative control.</item>
/// <item>The real chain's visit frequencies at a fixed τ against the exact posterior, the family alone and all options on.</item>
/// </list>
/// </summary>
public sealed class BassFactorMoveTests
{
    private readonly ITestOutputHelper _out;
    public BassFactorMoveTests(ITestOutputHelper output) => _out = output;

    private static BassModel Model(string fixture, string schedule = "equal", LocalSteps? steps = null)
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture(fixture, schedule, s => s with { FactorMoves = new FactorMoves(), LocalSteps = steps });
        return BassModel.Create(config, data);
    }

    [Theory]
    [InlineData("continuous")]
    [InlineData("categorical, R's measure")]
    [InlineData("categorical, uniform")]
    [InlineData("functional")]
    public void Ratios_CancelThePriorRatio_AsDerived(string fixture)
    {
        BassModel model = Model(fixture);
        BassTermProposal proposal = model.Proposal;
        BassPart[] parts = [.. Enum.GetValues<BassPart>().Where(model.Prior.Has)];
        var cases = new Dictionary<string, int>();
        double worst = 0.0;
        void Cancels(string kind, BassTerm from, BassTermMove move)
        {
            double gap = Math.Abs(move.LogProposalRatio + model.Prior.LogPrior(move.Term!) - model.Prior.LogPrior(from));
            worst = Math.Max(worst, gap);
            Assert.True(gap <= 1e-12, $"{kind} of {from} to {move.Term}: the ratios leave {gap:E2}");
            cases[kind] = cases.GetValueOrDefault(kind) + 1;
        }
        foreach (BassTerm t in BassKernelEnumeration.AllTerms(model))
            foreach (BassPart part in parts)
            {
                int k = BassStructurePrior.Order(t, part);
                HashSet<int> used = [.. BassTermCounts.Variables(t, part)];
                int[] unused = [.. Enumerable.Range(0, model.Prior.Variables(part)).Where(v => !used.Contains(v))];
                for (int f = 0; f < k; f++)
                {
                    BassTermMove pruned = proposal.Pruned(t, part, f);
                    if (!pruned.LeavesIntercept) Cancels(k == 1 ? "prune, part emptied (K = 1 → 0)" : "prune, K ≥ 1 left", t, pruned);
                    foreach (int v in unused)
                        foreach ((HingeFactor h, CategoricalFactor? c, _) in BassKernelEnumeration.FactorsOf(model, part, v))
                            Cancels("variable change", t, c is null ? proposal.Exchanged(t, part, f, h) : proposal.Exchanged(t, part, f, c));
                }
                if (k >= model.Prior.Cap(part)) continue;
                foreach (int v in unused)
                    foreach ((HingeFactor h, CategoricalFactor? c, _) in BassKernelEnumeration.FactorsOf(model, part, v))
                        Cancels(k == 0 ? "grow into an empty part (0 → 1)" : "grow, K ≥ 1 held", t, c is null ? proposal.Grown(t, part, h) : proposal.Grown(t, part, c));
            }
        _out.WriteLine($"{fixture}: {string.Join("; ", cases.Select(c => $"{c.Key} {c.Value}"))}; worst {worst:E2}");
        Assert.Contains("variable change", cases.Keys);
        if (fixture == "continuous") Assert.Contains("prune, K ≥ 1 left", cases.Keys);   // the other fixtures' caps are 1
        else   // one part: a term's only part is the term, and its K = 1 → 0 is the intercept
        {
            Assert.Contains("prune, part emptied (K = 1 → 0)", cases.Keys);
            Assert.Contains("grow into an empty part (0 → 1)", cases.Keys);
        }
    }

    [Theory]
    [InlineData("continuous", "equal", 1.0, false)]
    [InlineData("continuous", "cancelling", 1.0, false)]
    [InlineData("continuous", "equal", 0.3, false)]
    [InlineData("continuous", "equal", 1.0, true)]
    [InlineData("continuous", "cancelling", 0.3, true)]
    [InlineData("categorical, R's measure", "equal", 1.0, false)]
    [InlineData("categorical, R's measure", "cancelling", 0.3, true)]
    [InlineData("categorical, uniform", "equal", 1.0, true)]
    [InlineData("categorical, uniform", "cancelling", 0.3, false)]
    [InlineData("functional", "equal", 1.0, false)]
    [InlineData("functional", "cancelling", 0.3, true)]
    public void DetailedBalance_Enumerated(string fixture, string schedule, double inverseTemperature, bool withSteps)
    {
        LocalSteps? steps = withSteps ? new LocalSteps() : null;
        BassModel model = Model(fixture, schedule, steps);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model, inverseTemperature);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.All(space, key, s, steps));
        (Dictionary<string, int> transitions, double worst) = BassKernelEnumeration.Check(space, flows);
        _out.WriteLine($"{fixture}, {schedule}, b = {inverseTemperature}, steps {(withSteps ? "on" : "off")}: |F| = {space.States.Count}; " +
                       $"{string.Join(", ", transitions.Select(t => $"{t.Key} {t.Value}"))}; worst gap {worst:E2}");
        Assert.True(transitions[BassKernelEnumeration.PruneGrow] > 0);
        Assert.True(transitions[BassKernelEnumeration.VariableChange] > 0);
    }

    // The check's power at the case the decisions left to derivation: a prune that empties a part (and its grow back),
    // scored as though the order law gave order 0 twice the mass, breaks detailed balance.
    [Fact]
    public void DetailedBalance_Enumerated_RefusesAWrongPartEmptyingRatio()
    {
        BassModel model = Model("categorical, R's measure");
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        BassKernelEnumeration.Enumerated Wrong(BassKernelEnumeration.State s, BassKernelEnumeration.Enumerated e)
        {
            if (e.Group != BassKernelEnumeration.PruneGrow || e.Candidate is null) return e;
            BassTerm[] added = [.. e.Candidate.Except(s.Terms)], removed = [.. s.Terms.Except(e.Candidate)];
            if (added.Length != 1 || removed.Length != 1) return e;
            bool emptiesAPart = Enum.GetValues<BassPart>().Any(p => BassStructurePrior.Order(removed[0], p) == 1 && BassStructurePrior.Order(added[0], p) == 0);
            return emptiesAPart ? e with { LogRatio = e.LogRatio + Math.Log(2.0) } : e;
        }
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.All(space, key, s, null).Select(e => Wrong(s, e)));
        Assert.ThrowsAny<Xunit.Sdk.XunitException>(() => BassKernelEnumeration.Check(space, flows));
    }

    [Fact]
    public void Proposals_DrawAsStated()
    {
        // On the categorical fixture, x0·z2: each part holds one factor, and z1 is unused.
        BassModel model = Model("categorical, R's measure");
        BassTermProposal proposal = model.Proposal;
        var term = new BassTerm([new HingeFactor(0, 1, 1)], [new CategoricalFactor(1, [0, 3])]);
        var expected = new Dictionary<string, Dictionary<string, double>>
        {
            ["prune"] = new() { [$"{new BassTerm(categorical: term.Categorical)}"] = 0.5, [$"{new BassTerm(term.Continuous)}"] = 0.5 },
            ["grow"] = new() { ["null"] = 1.0 },            // each part is at its cap of one
            ["variable"] = new() { ["null"] = 0.5 },        // no unused continuous variable
        };
        foreach (int[] subset in BassKernelEnumeration.ProperSubsets(2))
            expected["variable"][$"{new BassTerm(term.Continuous, [new CategoricalFactor(0, subset)])}"] = 0.5 * Math.Exp(BassKernelEnumeration.LogSubsetMeasure(CategoricalSubsetMeasure.SizeUniform, 2, 1));
        Draws(proposal.ProposePrune, term, expected["prune"], 50_000, 61);
        Draws(proposal.ProposeGrow, term, expected["grow"], 50_000, 62);
        Draws(proposal.ProposeVariableChange, term, expected["variable"], 50_000, 63);

        // On the 2,263-state fixture: x1 alone grows x0 at every sign and knot, and changes to x0 likewise; x0·x1 is full.
        BassModel continuous = Model("continuous");
        var single = new BassTerm([new HingeFactor(1, -1, 2)]);
        var grows = new Dictionary<string, double>();
        var changes = new Dictionary<string, double>();
        foreach (int s in new[] { -1, 1 })
            for (int k = 0; k < 3; k++)
            {
                grows[$"{new BassTerm([new HingeFactor(0, s, k), single.Continuous[0]])}"] = 1.0 / 6;
                changes[$"{new BassTerm([new HingeFactor(0, s, k)])}"] = 1.0 / 6;
            }
        Draws(continuous.Proposal.ProposeGrow, single, grows, 60_000, 64);
        Draws(continuous.Proposal.ProposeVariableChange, single, changes, 60_000, 65);
        Draws(continuous.Proposal.ProposePrune, single, new() { ["refused"] = 1.0 }, 1_000, 66);
        var full = new BassTerm([new HingeFactor(0, 1, 0), new HingeFactor(1, -1, 2)]);
        Draws(continuous.Proposal.ProposeGrow, full, new() { ["null"] = 1.0 }, 1_000, 67);
        Draws(continuous.Proposal.ProposeVariableChange, full, new() { ["null"] = 1.0 }, 1_000, 68);
        Draws(continuous.Proposal.ProposePrune, full, new() { [$"{new BassTerm([full.Continuous[1]])}"] = 0.5, [$"{new BassTerm([full.Continuous[0]])}"] = 0.5 }, 50_000, 69);
    }

    // Each outcome's frequency over `draws` proposals from `term` against its stated probability, within five standard errors.
    private void Draws(Func<Xoshiro256PlusPlus, BassTerm, BassTermMove> propose, BassTerm term, Dictionary<string, double> expected, int draws, int seed)
    {
        Assert.Equal(1.0, expected.Values.Sum(), 12);
        var rng = new Xoshiro256PlusPlus(seed);
        var seen = new Dictionary<string, int>();
        for (int k = 0; k < draws; k++)
        {
            BassTermMove move = propose(rng, term);
            string key = move.LeavesIntercept ? "refused" : move.IsNull ? "null" : $"{move.Term}";
            seen[key] = seen.GetValueOrDefault(key) + 1;
        }
        Assert.Empty(seen.Keys.Except(expected.Keys));
        foreach ((string key, double p) in expected)
        {
            double observed = seen.GetValueOrDefault(key) / (double)draws, error = Math.Sqrt(p * (1 - p) / draws);
            Assert.True(Math.Abs(observed - p) <= 5 * error + 1e-12, $"{term}, {key}: {observed:F5} against {p:F5} (± {error:F5})");
        }
    }

    [Fact]
    public void Palette_SharesTheChangeSlot_AndOffKeepsD01s()
    {
        foreach (MoveSchedule schedule in new MoveSchedule[] { new EqualOverAvailable(), new PriorCancelling(0.4) })
        {
            BassConfig config = BassSweepTests.Config(schedule, maxBasis: 2);
            var off = new BassStructureMoves(BassModel.Create(config, BassSweepTests.Data()), [], 0);
            var on = new BassStructureMoves(BassModel.Create(config with { Sampling = config.Sampling with { FactorMoves = new FactorMoves(2, 1, 1) } },
                                                             BassSweepTests.Data()), [], 0);
            Assert.Equal(["birth", "death", "change"], off.Moves.Select(m => m.Key));
            Assert.Equal(["birth", "death", "change", "prune", "grow", "variable"], on.Moves.Select(m => m.Key));
            Assert.Equal(["death", "birth", "change", "grow", "prune", "variable"], on.Moves.Select(m => m.ReverseKey));
            BassTerm one = new([new HingeFactor(0, 1, 1)]), two = new([new HingeFactor(1, 1, 1)]);
            foreach (BassTerm[] terms in new BassTerm[][] { [], [one], [one, two] })
            {
                BassState state = new BassStructureMoves(BassModel.Create(config, BassSweepTests.Data()), terms, 0).Start;
                double[] w = [.. off.Moves.Select(m => m.Weight(state))], v = [.. on.Moves.Select(m => m.Weight(state))];
                Assert.Equal(w[0], v[0]);
                Assert.Equal(w[1], v[1]);
                Assert.Equal(w[2], v[2] + v[3] + v[4] + v[5], 15);
                Assert.Equal([w[2] * 0.5, w[2] * 0.125, w[2] * 0.125, w[2] * 0.25], v[2..]);
            }
        }
    }

    [Fact]
    public void AtNoTerms_UnderThePriorCancellingSchedule_TheFamilyIsANullMove()
    {
        BassConfig config = BassSweepTests.Config(new PriorCancelling(0.4)) with
        {
            Sampling = new BassSampling { Moves = new PriorCancelling(0.4), FactorMoves = new FactorMoves() },
        };
        var model = BassModel.Create(config, BassSweepTests.Data());
        var moves = new BassStructureMoves(model, [], Math.Log(0.05));
        foreach (IRjMove<BassState> move in moves.Moves.Skip(3))
        {
            Assert.True(move.Weight(moves.Start) > 0.0);
            Proposal<BassState> p = move.Propose(moves.Start, new Xoshiro256PlusPlus(3))!.Value;
            Assert.Same(moves.Start, p.Candidate);
            Assert.Equal(double.NegativeInfinity, p.LogProposalRatio);
        }
        // The chain attempts it and rejects it, never leaving the empty structure that way.
        var chain = new ReversibleJumpChain<BassState>(moves.Moves.Skip(3).ToArray(), model.Target, moves.Start, new Xoshiro256PlusPlus(4));
        for (int i = 0; i < 200; i++) Assert.Same(moves.Start, chain.Step());
        Assert.Equal((200L, 0L), (chain.Attempts, chain.Accepted));
    }

    [Theory]
    [InlineData("family")]
    [InlineData("all")]
    public void VisitFrequencies_MatchTheExactPosterior(string options)
        => BassLocalStepTests.VisitFrequencies(_out, options, BassLocalStepTests.OptionsFor(options));
}
