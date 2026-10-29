using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Lever 2, local steps inside a change (the local-move options brief, Step 2): with probability ρ a change steps its
/// factor instead of redrawing it, a hinge's knot by δ with its sign kept, a subset by one level's toggle.
/// <list type="bullet">
/// <item>The production draws against the statement: each outcome's frequency from a fixed term against its stated
/// probability, and the null moves at the knots' edges and at an empty or full subset.</item>
/// <item>Detailed balance by enumeration (<see cref="BassKernelEnumeration"/>) at a fixed τ on the 2,263-state fixture,
/// on the categorical fixture under both subset measures and on the functional one, on both schedules, at b = 1 and at
/// b = 0.3 on the default path and form, with D01's window and a wider one.</item>
/// <item>The real chain's visit frequencies at a fixed τ against the exact posterior on the 2,263-state fixture.</item>
/// </list>
/// The constant-likelihood laws are <see cref="BassSweepTests.ConstantLikelihood_SamplesThePriorOnF"/>'s, with the
/// options on.
/// </summary>
public sealed class BassLocalStepTests
{
    private readonly ITestOutputHelper _out;
    public BassLocalStepTests(ITestOutputHelper output) => _out = output;

    /// <summary>The 2,263-state fixture at a fixed τ, unit information.</summary>
    internal static BassConfig Continuous(MoveSchedule? moves = null)
        => BassSweepTests.Config(moves) with { Coefficients = new FixedG() };

    /// <summary>The categorical fixture (<see cref="BassKernelEnumeration.CategoricalData"/>): each part's cap 1, M ≤ 2.</summary>
    internal static BassConfig Categorical(CategoricalSubsetMeasure measure, MoveSchedule? moves = null) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxCategoricalInteraction = 1, MaxBasis = 2, CategoricalSubsets = measure },
        Coefficients = new FixedG(),
        Sampling = new BassSampling { Moves = moves ?? new EqualOverAvailable() },
    };

    /// <summary>The functional fixture (<see cref="BassKernelEnumeration.FunctionalData"/>): each part's cap 1, M ≤ 2.</summary>
    internal static BassConfig Functional(MoveSchedule? moves = null) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxFunctionalInteraction = 1, MaxBasis = 2 },
        Coefficients = new FixedG(),
        Sampling = new BassSampling { Moves = moves ?? new EqualOverAvailable() },
    };

    /// <summary>A fixture by name, with the options in <paramref name="sampling"/>'s place.</summary>
    internal static (BassConfig Config, BassData Data) Fixture(string name, string schedule, Func<BassSampling, BassSampling> options)
    {
        MoveSchedule moves = schedule == "equal" ? new EqualOverAvailable() : new PriorCancelling(0.4);
        (BassConfig config, BassData data) = name switch
        {
            "continuous" => (Continuous(moves), BassSweepTests.Data()),
            "categorical, R's measure" => (Categorical(CategoricalSubsetMeasure.SizeUniform, moves), BassKernelEnumeration.CategoricalData()),
            "categorical, uniform" => (Categorical(CategoricalSubsetMeasure.Uniform, moves), BassKernelEnumeration.CategoricalData()),
            "functional" => (Functional(moves), BassKernelEnumeration.FunctionalData()),
            _ => throw new ArgumentOutOfRangeException(nameof(name)),
        };
        return (config with { Sampling = options(config.Sampling) }, data);
    }

    [Fact]
    public void Enumeration_FindsTheSweepTestsF()
    {
        var model = BassModel.Create(Continuous(), BassSweepTests.Data());
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        List<(BassTerm[] Terms, double LogPrior)> f = BassSweepTests.EnumerateF(model);
        Assert.Equal(2_263, space.States.Count);
        Assert.Equal(f.Count, space.States.Count);
        foreach ((BassTerm[] terms, double logPrior) in f)
        {
            int i = space.Find(terms);
            Assert.True(i >= 0);
            Assert.Equal(logPrior, space.States[i].LogPi - model.Target.Evaluate(space.States[i].Loaded).LogMarginal, 1e-10);
        }
    }

    [Theory]
    [InlineData("continuous", "equal", 1.0, 0.02)]
    [InlineData("continuous", "cancelling", 1.0, 0.02)]
    [InlineData("continuous", "equal", 0.3, 0.02)]
    [InlineData("continuous", "cancelling", 0.3, 0.02)]
    [InlineData("continuous", "equal", 1.0, 0.7)]       // w = 3 over 3 knots: ±3 always steps off
    [InlineData("categorical, R's measure", "equal", 1.0, 0.02)]
    [InlineData("categorical, R's measure", "cancelling", 0.3, 0.02)]
    [InlineData("categorical, uniform", "equal", 1.0, 0.02)]
    [InlineData("categorical, uniform", "cancelling", 0.3, 0.02)]
    [InlineData("functional", "equal", 1.0, 0.02)]
    [InlineData("functional", "cancelling", 0.3, 0.6)]
    public void DetailedBalance_Enumerated(string fixture, string schedule, double inverseTemperature, double window)
    {
        var steps = new LocalSteps(0.5, window);
        (BassConfig config, BassData data) = Fixture(fixture, schedule, s => s with { LocalSteps = steps });
        var model = BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model, inverseTemperature);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, steps));
        (Dictionary<string, int> transitions, double worst) = BassKernelEnumeration.Check(space, flows);
        _out.WriteLine($"{fixture}, {schedule}, b = {inverseTemperature}, ω = {window}: |F| = {space.States.Count}; " +
                       $"{string.Join(", ", transitions.Select(t => $"{t.Key} {t.Value}"))}; worst gap {worst:E2}");
        Assert.True(transitions[BassKernelEnumeration.Step] > 0);
        Assert.True(transitions[BassKernelEnumeration.Redraw] > 0);
    }

    // The check's power: a step whose ratio is off by 0.05, or a toggle that leaves out the subset measure's prior ratio
    // under R's measure (which a uniform measure would not notice), breaks detailed balance.
    [Theory]
    [InlineData("continuous")]
    [InlineData("categorical, R's measure")]
    public void DetailedBalance_Enumerated_RefusesAWrongRatio(string fixture)
    {
        var steps = new LocalSteps(0.5, 0.02);
        (BassConfig config, BassData data) = Fixture(fixture, "equal", s => s with { LocalSteps = steps });
        var model = BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        BassKernelEnumeration.Enumerated Perturbed(BassKernelEnumeration.State s, BassKernelEnumeration.Enumerated e)
        {
            if (e.Group != BassKernelEnumeration.Step || e.Candidate is null) return e;
            if (fixture == "continuous") return e with { LogRatio = e.LogRatio + 0.05 };
            BassTerm[] added = [.. e.Candidate.Except(s.Terms)], removed = [.. s.Terms.Except(e.Candidate)];
            if (added.Length != 1 || removed.Length != 1) return e;   // onto a held term: refused either way
            return e with { LogRatio = e.LogRatio - (model.Prior.LogPrior(added[0]) - model.Prior.LogPrior(removed[0])) };
        }
        IEnumerable<BassKernelEnumeration.Enumerated> Wrong(string key, BassKernelEnumeration.State s)
            => BassKernelEnumeration.D01(space, key, s, steps).Select(e => Perturbed(s, e));
        var flows = BassKernelEnumeration.Kernel(space, Wrong);
        Assert.ThrowsAny<Xunit.Sdk.XunitException>(() => BassKernelEnumeration.Check(space, flows));
    }

    [Fact]
    public void ProposeChange_DrawsAsStated()
    {
        // Five knots per variable (N = 75) and ω = ½: w = 3, so the edges' null moves and the window's full width show.
        var steps = new LocalSteps(0.4, 0.5);
        BassConfig config = BassSweepTests.Config(maxBasis: 10) with { Sampling = new BassSampling { LocalSteps = steps } };
        var model = BassModel.Create(config, BassSweepTests.Data(large: true));
        Assert.Equal(3, model.Proposal.Window(BassPart.Continuous, 0));
        var term = new BassTerm([new HingeFactor(0, 1, 1), new HingeFactor(1, -1, 4)]);

        // The statement: a factor ½ each; the redraw 1 − ρ over 2·5 signs and knots; the step ρ over δ ∈ {±1, ±2, ±3}.
        var expected = new Dictionary<string, double>();
        void Add(string key, double p) => expected[key] = expected.GetValueOrDefault(key) + p;
        for (int pick = 0; pick < 2; pick++)
        {
            HingeFactor h = term.Continuous[pick];
            foreach (int sign in new[] { -1, 1 })
                for (int knot = 0; knot < 5; knot++)
                {
                    BassTerm t = BassKernelEnumeration.WithHinge(term, BassPart.Continuous, pick, new HingeFactor(h.Variable, sign, knot));
                    Add($"Redraw {t}", 0.5 * 0.6 / 10);
                }
            for (int delta = -3; delta <= 3; delta++)
            {
                if (delta == 0) continue;
                int knot = h.Knot + delta;
                Add(knot is < 0 or >= 5 ? "LocalStep null" : $"LocalStep {BassKernelEnumeration.WithHinge(term, BassPart.Continuous, pick, h with { Knot = knot })}",
                    0.5 * 0.4 / 6);
            }
        }
        Assert.Equal(1.0, expected.Values.Sum(), 12);

        const int draws = 400_000;
        var rng = new Xoshiro256PlusPlus(41);
        var seen = new Dictionary<string, int>();
        for (int k = 0; k < draws; k++)
        {
            BassTermMove move = model.Proposal.ProposeChange(rng, term);
            Assert.Equal(0.0, move.IsNull ? 0.0 : move.LogProposalRatio);
            string key = move.IsNull ? $"{move.Kind} null" : $"{move.Kind} {move.Term}";
            seen[key] = seen.GetValueOrDefault(key) + 1;
        }
        Assert.Empty(seen.Keys.Except(expected.Keys));
        foreach ((string key, double p) in expected)
        {
            double observed = seen.GetValueOrDefault(key) / (double)draws, error = Math.Sqrt(p * (1 - p) / draws);
            Assert.True(Math.Abs(observed - p) <= 5 * error, $"{key}: {observed:F5} against {p:F5} (± {error:F5})");
        }
    }

    [Fact]
    public void ProposeChange_TogglesAsStated_UnderBothMeasures()
    {
        foreach (CategoricalSubsetMeasure measure in new[] { CategoricalSubsetMeasure.SizeUniform, CategoricalSubsetMeasure.Uniform })
        {
            var steps = new LocalSteps(0.5, 0.02);
            (BassConfig config, BassData data) = Fixture(measure == CategoricalSubsetMeasure.SizeUniform ? "categorical, R's measure" : "categorical, uniform",
                                                         "equal", s => s with { LocalSteps = steps });
            var model = BassModel.Create(config, data);
            var term = new BassTerm(categorical: [new CategoricalFactor(1, [0, 2, 3])]);   // z2, three of four levels
            var expected = new Dictionary<string, double>();
            void Add(string key, double p) => expected[key] = expected.GetValueOrDefault(key) + p;
            foreach (int[] subset in BassKernelEnumeration.ProperSubsets(4))
                Add($"Redraw {BassKernelEnumeration.WithSubset(term, 0, new CategoricalFactor(1, subset))}",
                    0.5 * Math.Exp(BassKernelEnumeration.LogSubsetMeasure(measure, 4, subset.Length)));
            Add("Toggle null", 0.5 / 4);   // level 1 fills the set
            foreach (int[] toggled in new[] { new[] { 2, 3 }, [0, 3], [0, 2] })
                Add($"Toggle {BassKernelEnumeration.WithSubset(term, 0, new CategoricalFactor(1, toggled))}", 0.5 / 4);

            const int draws = 300_000;
            var rng = new Xoshiro256PlusPlus(43);
            var seen = new Dictionary<string, int>();
            for (int k = 0; k < draws; k++)
            {
                BassTermMove move = model.Proposal.ProposeChange(rng, term);
                if (move.Kind == BassMoveKind.Toggle && !move.IsNull) Assert.Equal(0.0, move.LogProposalRatio);
                if (move.Kind == BassMoveKind.Redraw)
                    Assert.Equal(model.Prior.LogPrior(term) - model.Prior.LogPrior(move.Term!), move.LogProposalRatio, 12);
                string key = move.IsNull ? $"{move.Kind} null" : $"{move.Kind} {move.Term}";
                seen[key] = seen.GetValueOrDefault(key) + 1;
            }
            Assert.Empty(seen.Keys.Except(expected.Keys));
            foreach ((string key, double p) in expected)
            {
                double observed = seen.GetValueOrDefault(key) / (double)draws, error = Math.Sqrt(p * (1 - p) / draws);
                Assert.True(Math.Abs(observed - p) <= 5 * error, $"{measure}, {key}: {observed:F5} against {p:F5} (± {error:F5})");
            }
        }
    }

    [Fact]
    public void Steps_AtTheEdges_AndOffByDefault()
    {
        var model = BassModel.Create(BassSweepTests.Config() with { Sampling = new BassSampling { LocalSteps = new LocalSteps() } }, BassSweepTests.Data());
        Assert.Equal(1, model.Proposal.Window(BassPart.Continuous, 0));   // ⌈0.02·3⌉
        var low = new BassTerm([new HingeFactor(0, 1, 0)]);
        var high = new BassTerm([new HingeFactor(0, -1, 2)]);
        Assert.True(model.Proposal.Stepped(low, 0, -1).IsNull);
        Assert.Equal(new BassTerm([new HingeFactor(0, 1, 1)]), model.Proposal.Stepped(low, 0, 1).Term);
        Assert.True(model.Proposal.Stepped(high, 0, 1).IsNull);
        Assert.Equal(new BassTerm([new HingeFactor(0, -1, 1)]), model.Proposal.Stepped(high, 0, -1).Term);
        Assert.Throws<ArgumentOutOfRangeException>(() => model.Proposal.Stepped(low, 0, 0));
        Assert.Throws<ArgumentException>(() => model.Proposal.Toggled(low, 0, 0));

        var categorical = BassModel.Create(Categorical(CategoricalSubsetMeasure.SizeUniform), BassKernelEnumeration.CategoricalData());
        var single = new BassTerm(categorical: [new CategoricalFactor(0, [1])]);   // z1 of two levels
        Assert.True(categorical.Proposal.Toggled(single, 0, 1).IsNull);   // empty
        Assert.True(categorical.Proposal.Toggled(single, 0, 0).IsNull);   // full
        Assert.Throws<ArgumentException>(() => categorical.Proposal.Stepped(single, 0, 1));

        // Off: no window, and the change draws exactly what D01's change draws.
        var off = BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data());
        Assert.Null(off.Proposal.Steps);
        Assert.Throws<InvalidOperationException>(() => off.Proposal.Window(BassPart.Continuous, 0));
        var a = new Xoshiro256PlusPlus(5);
        var b = new Xoshiro256PlusPlus(5);
        var term = new BassTerm([new HingeFactor(0, 1, 1), new HingeFactor(1, -1, 2)]);
        for (int k = 0; k < 1_000; k++)
        {
            BassTerm drawn = off.Proposal.DrawChange(a, term, out double ratio);
            BassTermMove move = off.Proposal.ProposeChange(b, term);
            Assert.Equal((BassMoveKind.Redraw, drawn, ratio), (move.Kind, move.Term, move.LogProposalRatio));
        }
        Assert.Equal(a.NextUInt64(), b.NextUInt64());
    }

    [Theory]
    [InlineData("steps")]
    public void VisitFrequencies_MatchTheExactPosterior(string options)
        => VisitFrequencies(_out, options, OptionsFor(options));

    /// <summary>The options a visit-frequency or constant-likelihood run turns on: "steps", "family" or "all".</summary>
    internal static Func<BassSampling, BassSampling> OptionsFor(string options) => options switch
    {
        "none" => s => s,
        "steps" => s => s with { LocalSteps = new LocalSteps() },
        "family" => s => s with { FactorMoves = new FactorMoves() },
        "all" => s => s with { LocalSteps = new LocalSteps(), FactorMoves = new FactorMoves() },
        _ => throw new ArgumentOutOfRangeException(nameof(options)),
    };

    /// <summary>
    /// The real chain at a fixed τ, the 2,263-state fixture: four chains' visit frequencies against the exact posterior,
    /// the count law, the order law and the ten likeliest structures, each within five batch-means errors.
    /// </summary>
    internal static void VisitFrequencies(ITestOutputHelper output, string name, Func<BassSampling, BassSampling> options)
    {
        BassConfig config = Continuous();
        config = config with { Sampling = options(config.Sampling) };
        var model = BassModel.Create(config, BassSweepTests.Data());
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        double[] pi = space.Posterior();
        int[] top = [.. Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(10)];

        const int chains = 4, burn = 2_000, sweeps = 150_000;
        var traces = new List<double[]>[top.Length + model.MaxBasis + 1];
        for (int k = 0; k < traces.Length; k++) traces[k] = [];
        for (int c = 0; c < chains; c++)
        {
            var chain = new BassChain(model, [], space.States[0].Loaded.LogTau, new Xoshiro256PlusPlus(500 + c), new Xoshiro256PlusPlus(600 + c));
            for (int i = 0; i < burn; i++) chain.Step();
            var hits = new double[traces.Length][];
            for (int k = 0; k < hits.Length; k++) hits[k] = new double[sweeps];
            for (int i = 0; i < sweeps; i++)
            {
                BassState s = chain.Step();
                int at = space.Find(s.Terms);
                Assert.True(at >= 0, "the chain left F");
                int rank = Array.IndexOf(top, at);
                if (rank >= 0) hits[rank][i] = 1;
                hits[top.Length + s.Count][i] = 1;
            }
            for (int k = 0; k < traces.Length; k++) traces[k].Add(hits[k]);
        }
        void Compare(string what, List<double[]> trace, double exact)
        {
            double observed = trace.Sum(t => t.Sum()) / (chains * sweeps), error = BatchError(trace);
            output.WriteLine($"{name}, {what}: {observed:F4} against {exact:F4} (± {error:F4})");
            Assert.True(Math.Abs(observed - exact) <= 5 * error + 0.002, $"{name}, {what}: {observed:F4} against {exact:F4} (± {error:F4})");
        }
        for (int m = 0; m <= model.MaxBasis; m++)
            Compare($"M = {m}", traces[top.Length + m], Enumerable.Range(0, pi.Length).Where(i => space.States[i].Terms.Length == m).Sum(i => pi[i]));
        for (int r = 0; r < top.Length; r++)
            Compare($"state {r} ({string.Join(" + ", space.States[top[r]].Terms.Select(t => t.ToString()))})", traces[r], pi[top[r]]);
    }

    // A batch-means standard error over every chain's trace.
    private static double BatchError(IReadOnlyList<double[]> traces, int batches = 25)
    {
        var means = new List<double>();
        foreach (double[] t in traces)
        {
            int size = t.Length / batches;
            for (int b = 0; b < batches; b++) means.Add(t.Skip(b * size).Take(size).Average());
        }
        double mean = means.Average();
        return Math.Sqrt(means.Sum(m => (m - mean) * (m - mean)) / (means.Count - 1) / means.Count);
    }
}
