using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// What F-M2's Phase 0 harnesses share (the local-move measurements brief): S10's three fixtures, their data and configs
/// as <see cref="BassQualificationRecoveryTests"/> builds them; each chain started, seeded and burned as the run does it,
/// so a fixed state is S10's chain's state at the end of burn-in; a proposal's log ratio split into its parts beside the
/// ratio the chain forms; and the check that a harness proposing from a fixed state committed nothing. They read, they do
/// not judge: qualification runs (<see cref="BassQualification.Enabled"/>).
/// </summary>
internal static class BassMixing
{
    /// <summary>S10's lengths: sweeps of burn-in, then sweeps sampled.</summary>
    internal const int BurnIn = 5000, Samples = 2000;

    /// <summary>S10's fixtures by name, for a theory's data.</summary>
    internal const string Friedman = "Friedman", Mixed = "mixed categorical", Functional = "direct functional";

    /// <summary>One of S10's fixtures: its data and config, each part's variable names, and on Friedman's surface whether a
    /// term uses an inert input (null where the truth has none).</summary>
    internal sealed record Fixture(string Name, BassData Data, BassConfig Config, Func<BassTerm, bool>? Inert)
    {
        public string[] Continuous => [.. Data.Predictors.OfType<ContinuousPredictor>().Select(p => p.Name)];
        public string[] Categorical => [.. Data.Predictors.OfType<CategoricalPredictor>().Select(p => p.Name)];
        public string[] Outputs => Data.Response is FunctionalResponse f ? [.. f.OutputVariables.Select(v => v.Name)] : [];

        /// <summary>The term's variables by name, in its parts' order: "x1·x2", "x3·z", "x0·t".</summary>
        public string Variables(BassTerm term) => string.Join("·",
            term.Continuous.Select(h => Continuous[h.Variable]).Order()
                .Concat(term.Categorical.Select(c => Categorical[c.Variable]).Order())
                .Concat(term.Functional.Select(h => Outputs[h.Variable]).Order()));

        /// <summary>The term by name, each hinge with its sign and its knot's quantile j/(n_v − 1) among the variable's
        /// distinct values, each subset with its levels' indices: "x1+0.14 · x2−0.65", "x3+0.50 · z{0,1}".</summary>
        public string Describe(BassModel model, BassTerm term) => string.Join(" · ",
            term.Continuous.Select(h => Hinge(Continuous[h.Variable], h, model.Inputs.ContinuousKnots(h.Variable).Length))
                .Concat(term.Categorical.Select(c => $"{Categorical[c.Variable]}{{{string.Join(",", c.Levels)}}}"))
                .Concat(term.Functional.Select(h => Hinge(Outputs[h.Variable], h, model.Inputs.FunctionalKnots(h.Variable).Length))));

        private static string Hinge(string name, HingeFactor h, int knots)
            => $"{name}{(h.Sign > 0 ? "+" : "−")}{(double)h.Knot / Math.Max(1, knots - 1):F2}";

        /// <summary>The term's tag: its orders, and on Friedman whether it is inert.</summary>
        public string Tag(BassTerm term) => Inert is null ? Orders(term) : $"{Orders(term)} {(Inert(term) ? "inert" : "active")}";
    }

    /// <summary>S10's fixture <paramref name="name"/>, tempered by <paramref name="tempering"/> when given.</summary>
    internal static Fixture Get(string name, TemperingSpec? tempering = null) => name switch
    {
        Friedman => new(name, BassQualificationRecoveryTests.Friedman().Data, Recovery(53, tempering), t => t.Continuous.Any(h => h.Variable >= 5)),
        Mixed => new(name, BassQualificationRecoveryTests.Mixed().Data, Recovery(57, tempering), null),
        Functional => new(name, BassQualificationRecoveryTests.Functional().Data, Recovery(61, tempering), null),
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    // BassQualificationRecoveryTests' config, whose lengths are BurnIn and Samples.
    private static BassConfig Recovery(int seed, TemperingSpec? tempering)
    {
        BassConfig config = BassQualificationRecoveryTests.Config(seed, tempering);
        if (config.Sampling.BurnIn != BurnIn || config.Sampling.Samples != Samples)
            throw new InvalidOperationException("S10's recovery lengths moved; BassMixing's BurnIn and Samples must follow.");
        return config;
    }

    /// <summary>Chain <paramref name="chain"/>'s start and its four streams' seeds, as the run takes them: the shell's
    /// chain seed from the master's tree, fanned out into the start's, the moves', τ's and a ladder's.</summary>
    internal static (BassState Start, int[] Seeds) Start(BassModel model, int chain)
    {
        BassSampling sampling = model.Config.Sampling;
        int[] chainSeeds = SeedTree.Derive(sampling.MasterSeed, 2 * sampling.Chains);
        int[] seeds = SeedTree.Derive(chainSeeds[chain], 4);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(seeds[0])).State;
        return (start, seeds);
    }

    /// <summary>Untempered chain <paramref name="chain"/>, started as the run starts it.</summary>
    internal static BassChain Chain(BassModel model, int chain)
    {
        (BassState start, int[] seeds) = Start(model, chain);
        return new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
    }

    /// <summary>Chain <paramref name="chain"/>'s ladder under <paramref name="spec"/>, started as the run starts it; on the
    /// declared rungs <paramref name="inverseTemperatures"/> when given.</summary>
    internal static BassLadder Ladder(BassModel model, int chain, TemperingSpec spec, double[]? inverseTemperatures = null)
    {
        (BassState start, int[] seeds) = Start(model, chain);
        return inverseTemperatures is null
            ? new BassLadder(model, spec, start.Terms, start.LogTau, seeds[3])
            : new BassLadder(model, spec, inverseTemperatures, start.Terms, start.LogTau, seeds[3]);
    }

    /// <summary>Every chain of <paramref name="model"/> burned as the run burns it, in parallel: its sweeps, then its
    /// tallies restarted.</summary>
    internal static BassChain[] Burned(BassModel model, int sweeps = BurnIn)
    {
        var chains = new BassChain[model.Config.Sampling.Chains];
        Parallel.For(0, chains.Length, c =>
        {
            BassChain chain = Chain(model, c);
            for (int s = 0; s < sweeps; s++) chain.Step();
            chain.Structure.ResetTallies();
            chains[c] = chain;
        });
        return chains;
    }

    /// <summary>A ladder's burn-in as the run takes it: an adaptive ladder re-spaced after rounds of 1, 2, 4, … sweeps, a
    /// self-sizing one sized after every round but the last; then frozen, its rungs untraced.</summary>
    internal static void Burn(BassLadder ladder, int steps = BurnIn)
    {
        int done = 0;
        if (ladder.Spec.Adaptive)
            for (int round = 1; done + round <= steps; round *= 2)
            {
                for (int k = 0; k < round; k++) ladder.Step();
                done += round;
                if (ladder.SizesItself && done + 2L * round <= steps) ladder.Resize();
                else ladder.Adapt();
            }
        for (; done < steps; done++) ladder.Step();
        ladder.EndBurnIn(traceRungs: false);
    }

    /// <summary>The term counts of <paramref name="state"/>'s structure, built afresh.</summary>
    internal static BassTermCounts Counts(BassModel model, BassState state)
    {
        var counts = new BassTermCounts(model.Prior);
        foreach (BassTerm t in state.Terms) counts.Add(t);
        return counts;
    }

    /// <summary>A term's orders by part, K continuous, Z categorical, T functional: "K2", "K1Z1", "K1T1".</summary>
    internal static string Orders(BassTerm term) => string.Concat(
        term.Continuous.Count > 0 ? $"K{term.Continuous.Count}" : "",
        term.Categorical.Count > 0 ? $"Z{term.Categorical.Count}" : "",
        term.Functional.Count > 0 ? $"T{term.Functional.Count}" : "");

    /// <summary>A move of <paramref name="moves"/> by its key.</summary>
    internal static IRjMove<BassState> Move(IReadOnlyList<IRjMove<BassState>> moves, string key) => moves.Single(m => m.Key == key);

    /// <summary>
    /// A scored proposal: refused (with whether its term repeats one the structure holds), or its log Metropolis–Hastings
    /// ratio split into Δℓ (the log marginal's change), the count prior's ratio, the term prior's, the within-move ratio
    /// and the selection ratio, beside <see cref="Ratio"/>, the ratio the chain forms.
    /// </summary>
    internal readonly record struct Scored(BassRefusal? Refusal, bool Repeated, double LogMarginalChange, double CountPrior,
                                           double TermPrior, double WithinMove, double Selection, double Ratio)
    {
        public bool Admitted => Refusal is null;

        /// <summary>min(1, e^ratio), 0 for a refusal.</summary>
        public double Acceptance => Admitted ? Math.Min(1.0, Math.Exp(Ratio)) : 0.0;

        /// <summary>The Δℓ at which the ratio reaches 0: what the other parts ask of the likelihood.</summary>
        public double Threshold => -(CountPrior + TermPrior + WithinMove + Selection);
    }

    /// <summary>
    /// Scores <paramref name="proposal"/>, which <paramref name="move"/> made from <paramref name="current"/>: refused by the
    /// design or the target, or its ratio's parts, the term prior's <paramref name="termPrior"/> and the within-move
    /// <paramref name="withinMove"/> taken from the caller, each found independently of the move. The ratio the chain
    /// forms is <see cref="ReversibleJumpChain{TState}.Step"/>'s, from the same target, temperature and weights; the parts
    /// must sum to it to rounding.
    /// </summary>
    internal static Scored Score(BassModel model, IReadOnlyList<IRjMove<BassState>> moves, IRjMove<BassState> move, BassState current,
                                 Proposal<BassState> proposal, double termPrior, double withinMove, bool repeated = false)
    {
        BassTarget target = model.Target;
        BassState candidate = proposal.Candidate;
        if (target.Refusal(candidate.Fit) is { } reason)
            return new Scored(reason, repeated, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NegativeInfinity);
        BassEvaluation now = target.Evaluate(current), next = target.Evaluate(candidate);
        double reverse = Move(moves, move.ReverseKey).Weight(candidate), totalNext = TotalWeight(moves, candidate);
        double selection = Math.Log(reverse / totalNext) - Math.Log(move.Weight(current) / TotalWeight(moves, current));
        double chain = reverse > 0.0 && totalNext > 0.0
            ? Temperature.Posterior.LogDensity(next) - Temperature.Posterior.LogDensity(now) + proposal.LogProposalRatio + proposal.LogJacobian + selection
            : double.NegativeInfinity;
        double dl = next.LogMarginal - now.LogMarginal;
        double count = model.Prior.LogCollectionWeight(candidate.Count) - model.Prior.LogCollectionWeight(current.Count);
        double sum = dl + count + termPrior + withinMove + selection;
        double scale = 1.0 + Math.Abs(dl) + Math.Abs(count) + Math.Abs(termPrior) + Math.Abs(withinMove) + Math.Abs(selection)
                       + Math.Abs(next.LogPrior) + Math.Abs(now.LogPrior) + Math.Abs(next.LogMarginal) + Math.Abs(now.LogMarginal);
        Assert.True(double.IsFinite(chain), $"The chain's ratio is {chain}.");
        Assert.True(Math.Abs(sum - chain) <= 1e-10 * scale, $"The parts sum to {sum:R}, the chain forms {chain:R}.");
        return new Scored(null, false, dl, count, termPrior, withinMove, selection, chain);
    }

    private static double TotalWeight(IReadOnlyList<IRjMove<BassState>> moves, BassState state)
    {
        double total = 0.0;
        foreach (IRjMove<BassState> m in moves)
            if (m.Weight(state) is var w && w > 0.0) total += w;
        return total;
    }

    /// <summary>The type-7 quantile <paramref name="q"/> of <paramref name="values"/>; NaN for none.</summary>
    internal static double Quantile(IEnumerable<double> values, double q)
    {
        double[] sorted = [.. values.Order()];
        return sorted.Length == 0 ? double.NaN : PosteriorInterval.Quantile(sorted, q);
    }

    /// <summary>
    /// A chain's design and counts at its current state, taken before a harness proposes from that state without stepping;
    /// <see cref="Check"/> asserts they came back as they were. The design's terms, statistics, coefficients, residual and
    /// factor bit for bit; its next proposals against a design restored from the copy taken here, bit for bit, and against
    /// a freshly loaded design, to rounding; and the counts, through births drawn by the moves and drawn afresh from counts
    /// rebuilt from the structure, which agree term for term and ratio for ratio only if the moves' counts are those.
    /// </summary>
    internal sealed class Untouched
    {
        private readonly BassModel _model;
        private readonly BassChain _chain;
        private readonly BassState _state;
        private readonly BassDesignSnapshot _before;

        internal Untouched(BassModel model, BassChain chain)
        {
            _model = model;
            _chain = chain;
            _state = chain.Current;
            _before = chain.Structure.DesignFor(_state);
        }

        /// <summary>The design's copy as taken: its statistics, factor, β and residual copied, its columns shared.</summary>
        internal BassDesignSnapshot Before => _before;

        internal void Check()
        {
            Assert.Same(_state, _chain.Current);
            BassDesign design = _chain.Structure.DesignAt(_state);
            Assert.Equal(_before.Terms, design.Terms);
            Assert.Equal(Bits(_before.ResidualSquares), Bits(design.ResidualSquares));
            Assert.Equal(Bits(_before.FittedSquares), Bits(design.FittedSquares));
            Assert.Equal(_before.Beta.Select(Bits), design.Coefficients.ToArray().Select(Bits));
            Assert.Equal(_before.Residual.Select(Bits), design.Residual.ToArray().Select(Bits));
            var factor = (BorderedCholesky)design.Factor;
            BorderedCholesky was = _before.Factor;
            Assert.Equal(was.Dimension, factor.Dimension);
            Assert.Equal(was.UpdatesSinceRebuild, factor.UpdatesSinceRebuild);
            for (int i = 0; i < factor.Dimension; i++)
            {
                Assert.Equal(Bits(was.InverseDiagonal(i)), Bits(factor.InverseDiagonal(i)));
                for (int j = 0; j <= i; j++)
                {
                    Assert.Equal(Bits(was.Entry(i, j)), Bits(factor.Entry(i, j)));
                    Assert.Equal(Bits(was.Factor(i, j)), Bits(factor.Factor(i, j)));
                }
            }

            // The next proposals: a birth the design admits (below the cap), a change of the first term and a death, against
            // the copy restored and a fresh load.
            var finder = _model.NewDesign();
            finder.Restore(_before);
            var rng = new Xoshiro256PlusPlus(seed: 97);
            BassTermCounts counts = Counts(_model, _state);
            BassTerm? birthProbe = null, changeProbe = null;
            for (int k = 0; k < 1000 && (birthProbe is null || changeProbe is null && _state.Count > 0); k++)
            {
                if (_model.Proposal.DrawBirth(rng, counts) is not { } t) continue;
                if (birthProbe is null && finder.ProposeBirth(t).Admitted) birthProbe = t;
                if (changeProbe is null && _state.Count > 0 && finder.ProposeChange(0, t).Admitted) changeProbe = t;
            }
            Assert.True(birthProbe is not null || _state.Count == _model.MaxBasis, "No birth below the cap is admitted.");
            Assert.True(changeProbe is not null || _state.Count == 0, "No change is admitted.");
            var restored = _model.NewDesign();
            restored.Restore(_before);
            var loaded = _model.NewDesign();
            Assert.Null(loaded.Load(_state.Terms));
            if (birthProbe is not null)
            {
                DesignProposal birth = design.ProposeBirth(birthProbe);
                Assert.Equal(restored.ProposeBirth(birthProbe), birth);
                Close(loaded.ProposeBirth(birthProbe), birth);
            }
            if (changeProbe is not null)
            {
                DesignProposal change = design.ProposeChange(0, changeProbe);
                Assert.Equal(restored.ProposeChange(0, changeProbe), change);
                Close(loaded.ProposeChange(0, changeProbe), change);
            }
            if (_state.Count > 0)
            {
                DesignProposal death = design.ProposeDeath(0);
                Assert.Equal(restored.ProposeDeath(0), death);
                Close(loaded.ProposeDeath(0), death);
            }

            // The counts: the moves' births from one stream against fresh counts' draws from its twin.
            IRjMove<BassState> births = Move(_chain.Structure.Moves, "birth");
            var moveRng = new Xoshiro256PlusPlus(seed: 101);
            var freshRng = new Xoshiro256PlusPlus(seed: 101);
            for (int k = 0; k < 50; k++)
            {
                Proposal<BassState>? p = births.Propose(_state, moveRng);
                BassTerm? fresh = _model.Proposal.DrawBirth(freshRng, counts);
                Assert.Equal(fresh is null, p is null);
                if (p is not { } made) continue;
                Assert.Equal(fresh, made.Candidate.Terms[^1]);
                Assert.Equal(Bits(_model.Proposal.LogBirthRatio(fresh!, counts)), Bits(made.LogProposalRatio));
            }
        }

        private static void Close(DesignProposal fresh, DesignProposal held)
        {
            Assert.Equal(fresh.Refusal, held.Refusal);
            if (!held.Admitted) return;
            Assert.Equal(fresh.Columns, held.Columns);
            Assert.True(Math.Abs(fresh.ResidualSquares - held.ResidualSquares) <= 1e-9 * held.ResidualSquares,
                        $"RSS {held.ResidualSquares:R} against a fresh load's {fresh.ResidualSquares:R}");
            Assert.True(Math.Abs(fresh.FittedSquares - held.FittedSquares) <= 1e-9 * (held.FittedSquares + held.ResidualSquares),
                        $"Q {held.FittedSquares:R} against a fresh load's {fresh.FittedSquares:R}");
        }
    }

    internal static long Bits(double value) => BitConverter.DoubleToInt64Bits(value);
}
