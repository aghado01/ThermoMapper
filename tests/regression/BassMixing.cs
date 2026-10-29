using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Distributions;
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
/// ratio the chain forms; and the check that a harness proposing from a fixed state committed nothing. The mode-exit
/// pricing's harnesses add the fixed states' check against their chains' modes, structures scored at fixed τ, the design
/// draft's fixed-M factor moves and removals of several terms. They read, they do not judge: qualification runs
/// (<see cref="BassQualification.Enabled"/>).
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
    /// must sum to it to rounding. At a rung's <paramref name="temperature"/> (the likelihood path; the posterior by
    /// default) Δℓ is the tempered log likelihood's change, ℓ_b, and the prior's parts stay whole.
    /// </summary>
    internal static Scored Score(BassModel model, IReadOnlyList<IRjMove<BassState>> moves, IRjMove<BassState> move, BassState current,
                                 Proposal<BassState> proposal, double termPrior, double withinMove, bool repeated = false,
                                 Temperature? temperature = null)
    {
        BassTarget target = model.Target;
        BassState candidate = proposal.Candidate;
        Temperature at = temperature ?? Temperature.Posterior;
        if (at.Path != TemperingPath.Likelihood) throw new ArgumentException("The parts split the likelihood path's ratio.", nameof(temperature));
        double b = at.InverseTemperature;
        if (target.Refusal(candidate.Fit) is { } reason)
            return new Scored(reason, repeated, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NegativeInfinity);
        BassEvaluation now = target.Evaluate(current), next = target.Evaluate(candidate);
        double reverse = Move(moves, move.ReverseKey).Weight(candidate), totalNext = TotalWeight(moves, candidate);
        double selection = Math.Log(reverse / totalNext) - Math.Log(move.Weight(current) / TotalWeight(moves, current));
        double chain = reverse > 0.0 && totalNext > 0.0
            ? at.LogDensity(next) - at.LogDensity(now) + proposal.LogProposalRatio + proposal.LogJacobian + selection
            : double.NegativeInfinity;
        double dl = next.LogLikelihood(b) - now.LogLikelihood(b);
        double count = model.Prior.LogCollectionWeight(candidate.Count) - model.Prior.LogCollectionWeight(current.Count);
        double sum = dl + count + termPrior + withinMove + selection;
        double scale = 1.0 + Math.Abs(dl) + Math.Abs(count) + Math.Abs(termPrior) + Math.Abs(withinMove) + Math.Abs(selection)
                       + Math.Abs(next.LogPrior) + Math.Abs(now.LogPrior) + Math.Abs(next.LogLikelihood(b)) + Math.Abs(now.LogLikelihood(b));
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
        private readonly BassStructureMoves _structure;
        private readonly Func<BassState> _current;
        private readonly BassState _state;
        private readonly BassDesignSnapshot _before;

        internal Untouched(BassModel model, BassChain chain) : this(model, chain.Structure, () => chain.Current)
        {
        }

        /// <summary>The check for a replica of a ladder: its structure moves, and its chain's current state.</summary>
        internal Untouched(BassModel model, BassStructureMoves structure, Func<BassState> current)
        {
            _model = model;
            _structure = structure;
            _current = current;
            _state = current();
            _before = structure.DesignFor(_state);
        }

        /// <summary>The design's copy as taken: its statistics, factor, β and residual copied, its columns shared.</summary>
        internal BassDesignSnapshot Before => _before;

        internal void Check()
        {
            Assert.Same(_state, _current());
            BassDesign design = _structure.DesignAt(_state);
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
            IRjMove<BassState> births = Move(_structure.Moves, "birth");
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

    // F-M2's mode exits (the mode-exit pricing brief): the fixed states checked against their chains' modes, structures
    // scored at fixed τ through the target, the fixed-M factor moves, and removals of several terms at once.

    /// <summary>A worse mode's least gap below the best chain's median log posterior, in nats: M4's worse chains sat 22–42
    /// below, the others 0–12.5.</summary>
    internal const double WorseGap = 20.0;

    /// <summary>
    /// A chain's fixed state, checked against the chain's own sampling phase (2,000 sweeps past burn-in, M4's): the burned
    /// state's log posterior and its rank among the phase's, the phase's median and range, and the median's gap below the
    /// best chain's. <see cref="Chain"/> holds the state the harnesses read: the burned chain, or, when its state lies
    /// outside the phase's range, a chain started at the phase's most frequent structure (<see cref="Substituted"/>).
    /// </summary>
    internal sealed record FixedState(int Index, BassChain Chain, double LogPosterior, double Rank, double Median, double Min, double Max,
                                      double Gap, bool Substituted)
    {
        public bool Worse => Gap >= WorseGap;

        public string Label => Worse ? $"{Index} (worse)" : $"{Index}";
    }

    /// <summary>
    /// Every chain of <paramref name="model"/> burned as the run burns it, each state checked against its chain's mode: a
    /// twin restored from the chain's snapshot continues it bit for bit through the sampling phase, and the run's retained
    /// draws must hold the twins' M draw for draw, so the states are S10's and the phases M4's.
    /// </summary>
    internal static FixedState[] FixedStates(BassModel model)
    {
        BassChain[] chains = Burned(model);
        int count = chains.Length;
        var sampled = new double[count][];
        var terms = new int[count][];
        var frequent = new (BassTerm[] Terms, double LogTau)[count];
        Parallel.For(0, count, c =>
        {
            BassChain twin = BassChain.Restore(model, chains[c].Snapshot());
            var held = new Dictionary<string, (int Count, BassTerm[] Terms, double LogTau)>();
            sampled[c] = new double[Samples];
            terms[c] = new int[Samples];
            for (int s = 0; s < Samples; s++)
            {
                BassState state = twin.Step();
                sampled[c][s] = LogPosterior(model, state);
                terms[c][s] = state.Count;
                string key = Key(state.Terms);
                held[key] = held.TryGetValue(key, out var h) ? (h.Count + 1, h.Terms, h.LogTau) : (1, [.. state.Terms], state.LogTau);
            }
            (int _, BassTerm[] top, double logTau) = held.Values.MaxBy(h => h.Count);
            frequent[c] = (top, logTau);
        });
        BassFit fit = Bass.Fit(model.Config, model.Data);
        for (int c = 0; c < count; c++) Assert.Equal(terms[c], fit.Draws[c].Select(d => d.Structure.TermIds.Length));

        double[] medians = [.. sampled.Select(s => Quantile(s, 0.5))];
        var states = new FixedState[count];
        for (int c = 0; c < count; c++)
        {
            double lp = LogPosterior(model, chains[c].Current), min = sampled[c].Min(), max = sampled[c].Max();
            bool inside = lp >= min && lp <= max;
            BassChain chain = chains[c];
            if (!inside)
            {
                int[] seeds = Start(model, c).Seeds;
                chain = new BassChain(model, frequent[c].Terms, frequent[c].LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
            }
            double rank = sampled[c].Count(x => x <= lp) / (double)Samples;
            states[c] = new FixedState(c, chain, lp, rank, medians[c], min, max, medians.Max() - medians[c], !inside);
        }
        return states;
    }

    /// <summary>The fixed states' check as a table: each state's log posterior against its chain's sampling phase.</summary>
    internal static IEnumerable<string> FixedTable(BassModel model, FixedState[] states)
    {
        yield return "chain | M | log post | rank in its phase | phase median | phase range | gap below the best median | read at";
        yield return "--- | --- | --- | --- | --- | --- | --- | ---";
        foreach (FixedState s in states)
            yield return $"{s.Label} | {s.Chain.Current.Count} | {s.LogPosterior:F1} | {s.Rank:P1} | {s.Median:F1} | [{s.Min:F1}, {s.Max:F1}] | {s.Gap:F1} | " +
                         (s.Substituted ? $"the phase's most frequent structure ({LogPosterior(model, s.Chain.Current):F1})" : "the burned state");
    }

    /// <summary>The unnormalized log posterior M4 reads: the log marginal plus the log structure prior, τ's density left out.</summary>
    internal static double LogPosterior(BassModel model, BassState state)
        => model.Target.Evaluate(state).LogMarginal + model.Target.LogStructurePrior(state.Fit);

    /// <summary>A structure's exact identity, its terms' factors in a canonical order.</summary>
    internal static string Key(IEnumerable<BassTerm> terms) => string.Join(" + ", terms.Select(Key).Order(StringComparer.Ordinal));

    private static string Key(BassTerm t) => string.Join("·",
        t.Continuous.Select(h => $"x{h.Variable}{(h.Sign > 0 ? "+" : "-")}{h.Knot}")
            .Concat(t.Categorical.Select(c => $"z{c.Variable}{{{string.Join(",", c.Levels)}}}"))
            .Concat(t.Functional.Select(h => $"t{h.Variable}{(h.Sign > 0 ? "+" : "-")}{h.Knot}")));

    /// <summary>The sum of the terms' ln π.</summary>
    internal static double LogTermPrior(BassModel model, IEnumerable<BassTerm> terms)
    {
        double sum = 0.0;
        foreach (BassTerm t in terms) sum += model.Prior.LogPrior(t);
        return sum;
    }

    /// <summary>The structure <paramref name="terms"/> at u = <paramref name="logTau"/>, loaded whole on a fresh design in
    /// the given order; null, with the design's or the target's refusal, when it lies outside F.</summary>
    internal static BassState? Loaded(BassModel model, IReadOnlyList<BassTerm> terms, double logTau, out BassRefusal? refusal)
    {
        BassDesign design = model.NewDesign();
        refusal = design.Load(terms);
        if (refusal is not null) return null;
        BassTerm[] loaded = [.. design.Terms];
        var fit = BassStructureFit.Of(design, LogTermPrior(model, loaded));
        refusal = model.Target.Refusal(fit);
        return refusal is null ? new BassState(loaded, fit, logTau) : null;
    }

    /// <summary>The candidate <paramref name="terms"/> a design proposed (<paramref name="proposal"/>) at u =
    /// <paramref name="logTau"/>; null, with its refusal, when it lies outside F.</summary>
    internal static BassState? Proposed(BassModel model, BassTerm[] terms, DesignProposal proposal, double logTau, out BassRefusal? refusal)
    {
        var fit = BassStructureFit.Of(proposal, proposal.Admitted ? LogTermPrior(model, terms) : double.NaN);
        refusal = model.Target.Refusal(fit);
        return refusal is null ? new BassState(terms, fit, logTau) : null;
    }

    /// <summary>
    /// A structure move scored at fixed τ against the state it leaves, at inverse temperature b on the likelihood path:
    /// Δℓ_b, the tempered log likelihood's change in the target's form (the log marginal's at b = 1); the target's log
    /// structure prior ratio, untempered; the caller's proposal ratio; and <see cref="Ratio"/>, their sum, the move's log
    /// ratio with its selection taken as matched. Refused: the reason, and whether the term repeats one held.
    /// </summary>
    internal readonly record struct Exit(BassRefusal? Refusal, bool Repeated, double LogLikelihoodChange, double LogPriorRatio, double LogProposalRatio)
    {
        public bool Admitted => Refusal is null;

        public double Ratio => Admitted ? LogLikelihoodChange + LogPriorRatio + LogProposalRatio : double.NegativeInfinity;

        /// <summary>min(1, e^ratio), 0 for a refusal.</summary>
        public double Acceptance => Admitted ? Math.Min(1.0, Math.Exp(Ratio)) : 0.0;

        public static Exit Refused(BassRefusal reason, bool repeated) => new(reason, repeated, double.NaN, double.NaN, double.NaN);
    }

    /// <summary>The move from <paramref name="from"/> to <paramref name="to"/>, both at the same τ, scored at
    /// <paramref name="inverseTemperature"/>.</summary>
    internal static Exit ExitTo(BassModel model, BassState from, BassState to, double logProposalRatio, double inverseTemperature = 1.0)
    {
        BassTarget target = model.Target;
        BassEvaluation now = target.Evaluate(from), next = target.Evaluate(to);
        double b = inverseTemperature;
        Assert.Equal(Bits(from.LogTau), Bits(to.LogTau));
        return new Exit(null, false, next.LogLikelihood(b) - now.LogLikelihood(b),
                        target.LogStructurePrior(to.Fit) - target.LogStructurePrior(from.Fit), logProposalRatio);
    }

    /// <summary>Term <paramref name="index"/> of <paramref name="state"/> replaced by <paramref name="term"/> through
    /// <paramref name="design"/>, which holds the state (its proposal left pending, uncommitted), scored at
    /// <paramref name="inverseTemperature"/>; the candidate, null when refused.</summary>
    internal static (Exit Exit, BassState? Candidate) Replace(BassModel model, BassDesign design, BassState state, int index, BassTerm term,
                                                              double logProposalRatio, double inverseTemperature = 1.0)
    {
        bool repeated = state.Terms.Contains(term);
        DesignProposal proposal = design.ProposeChange(index, term);
        // The design moves a changed term to the end, so the candidate's order follows it.
        BassTerm[] terms = [.. state.Terms.Take(index), .. state.Terms.Skip(index + 1), term];
        BassState? candidate = Proposed(model, terms, proposal, state.LogTau, out BassRefusal? refusal);
        return candidate is null
            ? (Exit.Refused(refusal!.Value, repeated), null)
            : (ExitTo(model, state, candidate, logProposalRatio, inverseTemperature), candidate);
    }

    /// <summary>A term's factor: its part, and its index among that part's factors.</summary>
    internal readonly record struct FactorAt(BassPart Part, int Index);

    /// <summary>The term's factors, continuous, then categorical, then functional.</summary>
    internal static IEnumerable<FactorAt> Factors(BassTerm term)
        => Enumerable.Range(0, term.Continuous.Count).Select(i => new FactorAt(BassPart.Continuous, i))
            .Concat(Enumerable.Range(0, term.Categorical.Count).Select(i => new FactorAt(BassPart.Categorical, i)))
            .Concat(Enumerable.Range(0, term.Functional.Count).Select(i => new FactorAt(BassPart.Functional, i)));

    /// <summary>The term's order in <paramref name="part"/>.</summary>
    internal static int Order(BassTerm term, BassPart part) => BassStructurePrior.Order(term, part);

    /// <summary>The factor's variable.</summary>
    internal static int Variable(BassTerm term, FactorAt f) => f.Part switch
    {
        BassPart.Continuous => term.Continuous[f.Index].Variable,
        BassPart.Categorical => term.Categorical[f.Index].Variable,
        _ => term.Functional[f.Index].Variable,
    };

    /// <summary>The term with factor <paramref name="f"/> dropped; null when it was the term's only factor (the intercept).</summary>
    internal static BassTerm? Without(BassTerm term, FactorAt f)
    {
        if (term.Continuous.Count + term.Categorical.Count + term.Functional.Count == 1) return null;
        return new BassTerm(f.Part == BassPart.Continuous ? term.Continuous.Where((_, i) => i != f.Index) : term.Continuous,
                            f.Part == BassPart.Categorical ? term.Categorical.Where((_, i) => i != f.Index) : term.Categorical,
                            f.Part == BassPart.Functional ? term.Functional.Where((_, i) => i != f.Index) : term.Functional);
    }

    /// <summary>The term with a hinge (<paramref name="hinge"/>) or a subset (<paramref name="subset"/>) added to
    /// <paramref name="part"/>.</summary>
    internal static BassTerm With(BassTerm term, BassPart part, HingeFactor hinge, CategoricalFactor? subset) => part switch
    {
        BassPart.Continuous => new BassTerm([.. term.Continuous, hinge], term.Categorical, term.Functional),
        BassPart.Categorical => new BassTerm(term.Continuous, [.. term.Categorical, subset!], term.Functional),
        _ => new BassTerm(term.Continuous, term.Categorical, [.. term.Functional, hinge]),
    };

    /// <summary>The term with factor <paramref name="f"/> replaced by a hinge (<paramref name="hinge"/>) or a subset
    /// (<paramref name="subset"/>) of the same part.</summary>
    internal static BassTerm Replaced(BassTerm term, FactorAt f, HingeFactor hinge, CategoricalFactor? subset)
    {
        switch (f.Part)
        {
            case BassPart.Continuous:
                HingeFactor[] continuous = [.. term.Continuous];
                continuous[f.Index] = hinge;
                return new BassTerm(continuous, term.Categorical, term.Functional);
            case BassPart.Categorical:
                CategoricalFactor[] categorical = [.. term.Categorical];
                categorical[f.Index] = subset!;
                return new BassTerm(term.Continuous, categorical, term.Functional);
            default:
                HingeFactor[] functional = [.. term.Functional];
                functional[f.Index] = hinge;
                return new BassTerm(term.Continuous, term.Categorical, functional);
        }
    }

    /// <summary>The hinge factor <paramref name="f"/> stepped by <paramref name="delta"/> distinct values, its sign kept;
    /// null when the knot would leave its variable's values.</summary>
    internal static BassTerm? Stepped(BassModel model, BassTerm term, FactorAt f, int delta)
    {
        if (f.Part == BassPart.Categorical) throw new ArgumentException("A subset has no knot.", nameof(f));
        HingeFactor h = f.Part == BassPart.Continuous ? term.Continuous[f.Index] : term.Functional[f.Index];
        int knot = h.Knot + delta;
        if (knot < 0 || knot >= Knots(model, f.Part, h.Variable)) return null;
        return Replaced(term, f, new HingeFactor(h.Variable, h.Sign, knot), null);
    }

    /// <summary>The distinct values a hinge on <paramref name="variable"/> of <paramref name="part"/> may take its knot at.</summary>
    internal static int Knots(BassModel model, BassPart part, int variable)
        => part == BassPart.Continuous ? model.Inputs.ContinuousKnots(variable).Length : model.Inputs.FunctionalKnots(variable).Length;

    /// <summary>
    /// ln π of one factor under D01's term measure, found here from the inputs, not through the prior: a hinge's sign and
    /// knot, −ln 2 − ln n_v; a subset of <paramref name="size"/> levels by the declared measure, R's −ln(n − 1) −
    /// ln C(n, size) or the paper's −ln(2ⁿ − 2).
    /// </summary>
    internal static double LogFactorPrior(BassModel model, BassPart part, int variable, int size = 0)
    {
        if (part != BassPart.Categorical) return -Math.Log(2.0) - Math.Log(Knots(model, part, variable));
        int n = model.Inputs.LevelCount(variable);
        return model.Prior.Subsets == CategoricalSubsetMeasure.SizeUniform
            ? -Math.Log(n - 1.0) - LogChoose(n, size)
            : -Math.Log(Math.Pow(2.0, n) - 2.0);
    }

    /// <summary>The factor's ln π (<see cref="LogFactorPrior"/>).</summary>
    internal static double LogFactorPrior(BassModel model, BassTerm term, FactorAt f)
        => LogFactorPrior(model, f.Part, Variable(term, f), f.Part == BassPart.Categorical ? term.Categorical[f.Index].Levels.Count : 0);

    /// <summary>
    /// The design draft's prune of factor <paramref name="f"/> (K → K − 1 within its part), its proposal ratio: the
    /// reverse grow's draw of that variable among the p − K + 1 unused and of its factor, over the prune's pick 1/K,
    /// ln K − ln(p − K + 1) + ln π_f. With the term prior's C(p, K)/C(p, K − 1)/π_f it cancels.
    /// </summary>
    internal static double LogPruneRatio(BassModel model, BassTerm term, FactorAt f)
    {
        int p = model.Prior.Variables(f.Part), k = Order(term, f.Part);
        return Math.Log(k) - Math.Log(p - k + 1.0) + LogFactorPrior(model, term, f);
    }

    /// <summary>The draft's grow of a factor of ln π <paramref name="logFactor"/> into <paramref name="part"/> at order K,
    /// its proposal ratio: the reverse prune's pick 1/(K + 1) over the grow's draw, ln(p − K) − ln(K + 1) − ln π_f.</summary>
    internal static double LogGrowRatio(BassModel model, BassTerm term, BassPart part, double logFactor)
    {
        int p = model.Prior.Variables(part), k = Order(term, part);
        return Math.Log(p - k) - Math.Log(k + 1.0) - logFactor;
    }

    /// <summary>ln C(n, k).</summary>
    internal static double LogChoose(int n, int k)
        => SpecialFunctions.LogGamma(n + 1.0) - SpecialFunctions.LogGamma(k + 1.0) - SpecialFunctions.LogGamma(n - k + 1.0);

    /// <summary>A factor for <paramref name="variable"/> of <paramref name="part"/> from the prior: a hinge's sign and
    /// knot uniformly, or a subset from the declared measure. Returns its ln π.</summary>
    internal static double DrawFactor(BassModel model, Xoshiro256PlusPlus rng, BassPart part, int variable, out HingeFactor hinge, out CategoricalFactor? subset)
    {
        if (part == BassPart.Categorical)
        {
            hinge = default;
            subset = model.Prior[BassPart.Categorical]!.DrawSubset(rng, variable, model.Prior.Subsets);
            return LogFactorPrior(model, part, variable, subset.Levels.Count);
        }
        hinge = new HingeFactor(variable, BassDraws.Sign(rng), BassDraws.Below(rng, Knots(model, part, variable)));
        subset = null;
        return LogFactorPrior(model, part, variable);
    }

    /// <summary>
    /// A removal of the terms <see cref="Terms"/> (indices into the state) at once, scored at fixed τ; its proposal ratio
    /// is the reverse's: the k terms born in sequence by D01's birth from the reduced structure, summed over their k!
    /// orders, over the removal's pick 1/C(M, k), with the move's selection taken as matched. For one term it is the
    /// death's ratio.
    /// </summary>
    internal sealed record Removal(int[] Terms, Exit Exit, BassState? Reduced);

    /// <summary>Every removal of <paramref name="size"/> terms from <paramref name="state"/>, each reduced structure loaded
    /// whole.</summary>
    internal static Removal[] Removals(BassModel model, BassState state, int size)
    {
        int m = state.Count;
        var removals = new List<Removal>();
        foreach (int[] set in Subsets(m, size))
        {
            BassTerm[] rest = [.. state.Terms.Where((_, i) => Array.IndexOf(set, i) < 0)];
            BassState? reduced = Loaded(model, rest, state.LogTau, out BassRefusal? refusal);
            if (reduced is null)
            {
                removals.Add(new Removal(set, Exit.Refused(refusal!.Value, false), null));
                continue;
            }
            BassTerm[] removed = [.. set.Select(i => state.Terms[i])];
            double reverse = LogSequentialBirths(model, rest, removed) + LogChoose(m, size);
            removals.Add(new Removal(set, ExitTo(model, state, reduced, reverse), reduced));
        }
        return [.. removals];
    }

    /// <summary>ln Σ over the orders of <paramref name="born"/> of Π q_b(tⱼ | <paramref name="reduced"/> + the terms born
    /// before it): D01's births of those terms in sequence, in any order.</summary>
    internal static double LogSequentialBirths(BassModel model, IReadOnlyList<BassTerm> reduced, BassTerm[] born)
    {
        double total = double.NegativeInfinity;
        foreach (int[] order in Permutations(born.Length))
        {
            BassTermCounts counts = Counts(model, reduced);
            double lq = 0.0;
            foreach (int j in order)
            {
                lq += model.Proposal.LogBirth(born[j], counts);
                counts.Add(born[j]);
            }
            total = LogAdd(total, lq);
        }
        return total;
    }

    /// <summary>The counts of the structure <paramref name="terms"/>.</summary>
    internal static BassTermCounts Counts(BassModel model, IEnumerable<BassTerm> terms)
    {
        var counts = new BassTermCounts(model.Prior);
        foreach (BassTerm t in terms) counts.Add(t);
        return counts;
    }

    private static double LogAdd(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double hi = Math.Max(a, b);
        return hi + Math.Log(Math.Exp(a - hi) + Math.Exp(b - hi));
    }

    // The k-subsets of {0, …, m − 1}, ascending.
    private static IEnumerable<int[]> Subsets(int m, int k)
    {
        var set = new int[k];
        for (int i = 0; i < k; i++) set[i] = i;
        while (k <= m)
        {
            yield return [.. set];
            int j = k - 1;
            while (j >= 0 && set[j] == m - k + j) j--;
            if (j < 0) yield break;
            set[j]++;
            for (int i = j + 1; i < k; i++) set[i] = set[i - 1] + 1;
        }
    }

    // The orders of {0, …, k − 1}.
    private static IEnumerable<int[]> Permutations(int k)
    {
        if (k == 1) { yield return [0]; yield break; }
        foreach (int[] rest in Permutations(k - 1))
            for (int at = 0; at < k; at++)
                yield return [.. rest[..at], k - 1, .. rest[at..]];
    }
}
