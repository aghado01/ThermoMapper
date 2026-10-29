using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's structure kernels, enumerated on small fixtures (the local-move options brief's exactness checks): every term
/// and every structure of F, π(S | τ, y) from the log structure prior plus the collapsed log marginal, and each move's
/// transition probabilities K(S, S′), at a fixed τ.
/// <list type="bullet">
/// <item><b>The ground truth is the moves' statement.</b> Each proposal is enumerated from the selection rules as the
/// decisions state them (the term, the part, the factor, the variable, the sign and knot, δ or the level, with their
/// probabilities), and its candidate is built here. The production code contributes the acceptance only: the move
/// weights of the schedule, the target's log density, and the move's log ratio through its internal accessors, whose
/// candidate must equal the one built here.</item>
/// <item><b>The checks.</b> For each kernel, π(S)K(S, S′) = π(S′)K(S′, S) for every pair, in logs to
/// <see cref="Tolerance"/>, and πK = π; and πK = π for the whole structure kernel, births and deaths included.</item>
/// </list>
/// The fixtures: <see cref="BassSweepTests.Data"/>'s 2,263 states, continuous only, and two with several parts:
/// <see cref="CategoricalData"/> (a continuous predictor and two categorical ones, of two and four levels) and
/// <see cref="FunctionalData"/> (a continuous predictor and two output variables).
/// </summary>
internal static class BassKernelEnumeration
{
    /// <summary>The detailed-balance and stationarity tolerance, in logs: 10⁻⁹ relative.</summary>
    internal const double Tolerance = 1e-9;

    /// <summary>A kernel's group: a move, or a move and its reverse.</summary>
    internal const string BirthDeath = "birth/death", Redraw = "redraw", Step = "local step", PruneGrow = "prune/grow",
                          VariableChange = "variable change";

    /// <summary>One structure of F: its terms in the design's order, its term ids sorted, the loaded state and its
    /// counts, the enumeration's log π, and the production target's log density at the space's b.</summary>
    internal sealed record State(BassTerm[] Terms, int[] Ids, BassState Loaded, BassTermCounts Counts, double LogPi, double LogTarget);

    /// <summary>The enumerated space: the model at its fixed τ and b, every term, and F.</summary>
    internal sealed class Space
    {
        private readonly Dictionary<long, int> _index = [];

        internal Space(BassModel model, double inverseTemperature, List<BassTerm> terms, List<State> states)
        {
            Model = model;
            InverseTemperature = inverseTemperature;
            Terms = terms;
            for (int k = 0; k < terms.Count; k++) TermIds.Add(terms[k], k);
            States = states;
            for (int i = 0; i < states.Count; i++) _index.Add(Key(states[i].Ids), i);
        }

        public BassModel Model { get; }
        public double InverseTemperature { get; }
        public List<BassTerm> Terms { get; }
        public Dictionary<BassTerm, int> TermIds { get; } = [];
        public List<State> States { get; }

        /// <summary>The state holding <paramref name="terms"/>, in any order; −1 outside F (repeated terms included).</summary>
        public int Find(IEnumerable<BassTerm> terms)
        {
            int[] ids = [.. terms.Select(t => TermIds[t]).Order()];
            for (int k = 1; k < ids.Length; k++)
                if (ids[k] == ids[k - 1]) return -1;
            return _index.TryGetValue(Key(ids), out int i) ? i : -1;
        }

        /// <summary>π normalized, by state.</summary>
        public double[] Posterior()
        {
            double top = States.Max(s => s.LogPi);
            double[] p = [.. States.Select(s => Math.Exp(s.LogPi - top))];
            double total = p.Sum();
            return [.. p.Select(v => v / total)];
        }

        private static long Key(int[] ids)
        {
            long key = ids.Length;
            foreach (int id in ids) key = key * (1L << 20) + id + 1;
            return key;
        }
    }

    // Two categorical predictors of two and four levels beside a continuous one on {0, ½, 1}, each combination twice (n = 24).
    internal static BassData CategoricalData(int seed = 17)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int n = 24;
        var x = new double[n];
        var z1 = new int[n];
        var z2 = new int[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i % 3) / 2.0;
            z1[i] = (i / 3) % 2;
            z2[i] = (i / 6) % 4;
            y[i] = 1 + x[i] + (z2[i] <= 1 ? 1.0 : 0.0) + 0.5 * (z1[i] == 0 ? x[i] : 0.0) + 0.3 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x), new CategoricalPredictor("z1", z1, ["a", "b"]),
                             new CategoricalPredictor("z2", z2, ["p", "q", "r", "s"])], new ScalarResponse(y));
    }

    // A direct functional response: a continuous predictor on {0, ½, 1} twice (n = 6), two output variables over four
    // locations, t1 on three values and t2 on two (N = 24).
    internal static BassData FunctionalData(int seed = 19)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] t1 = [0, 0.5, 1, 0.5], t2 = [0, 0, 1, 1];
        const int n = 6;
        var x = new double[n];
        var y = new double[n, t1.Length];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i % 3) / 2.0;
            for (int r = 0; r < t1.Length; r++) y[i, r] = 1 + x[i] + t1[r] * x[i] + 0.5 * t2[r] + 0.2 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t1", t1), new ContinuousPredictor("t2", t2)]));
    }

    /// <summary>Every term of <paramref name="model"/>'s term space: each part's order up to its cap, its variables, and
    /// each variable's signs and knots or proper subsets; the empty term, the intercept, left out.</summary>
    internal static List<BassTerm> AllTerms(BassModel model)
    {
        BassStructurePrior prior = model.Prior;
        List<HingeFactor[]> HingeSets(BassPart part, Func<int, int> knots)
        {
            var sets = new List<HingeFactor[]> { Array.Empty<HingeFactor>() };
            if (!prior.Has(part)) return sets;
            foreach (int[] variables in Subsets(prior.Variables(part), prior.Cap(part)))
            {
                IEnumerable<HingeFactor[]> built = [[]];
                foreach (int v in variables)
                    built = built.SelectMany(prefix => new[] { -1, 1 }.SelectMany(s => Enumerable.Range(0, knots(v)).Select(k => (HingeFactor[])[.. prefix, new HingeFactor(v, s, k)])));
                sets.AddRange(built);
            }
            return sets;
        }
        var subsetSets = new List<CategoricalFactor[]> { Array.Empty<CategoricalFactor>() };
        if (prior.Has(BassPart.Categorical))
            foreach (int[] variables in Subsets(prior.Variables(BassPart.Categorical), prior.Cap(BassPart.Categorical)))
            {
                IEnumerable<CategoricalFactor[]> built = [[]];
                foreach (int v in variables)
                    built = built.SelectMany(prefix => ProperSubsets(model.Inputs.LevelCount(v)).Select(levels => (CategoricalFactor[])[.. prefix, new CategoricalFactor(v, levels)]));
                subsetSets.AddRange(built);
            }
        List<HingeFactor[]> continuous = HingeSets(BassPart.Continuous, v => model.Inputs.ContinuousKnots(v).Length);
        List<HingeFactor[]> functional = HingeSets(BassPart.Functional, v => model.Inputs.FunctionalKnots(v).Length);
        var terms = new List<BassTerm>();
        foreach (HingeFactor[] c in continuous)
            foreach (CategoricalFactor[] z in subsetSets)
                foreach (HingeFactor[] f in functional)
                    if (c.Length + z.Length + f.Length > 0) terms.Add(new BassTerm(c, z, f));
        return terms;
    }

    // The variable sets of 1 … cap of p variables, ascending.
    private static IEnumerable<int[]> Subsets(int p, int cap)
    {
        for (int k = 1; k <= cap; k++)
            foreach (int[] s in Choose(p, k)) yield return s;
    }

    private static IEnumerable<int[]> Choose(int p, int k, int from = 0)
    {
        if (k == 0)
        {
            yield return [];
            yield break;
        }
        for (int v = from; v <= p - k; v++)
            foreach (int[] rest in Choose(p, k - 1, v + 1)) yield return [v, .. rest];
    }

    /// <summary>Every proper, nonempty subset of <paramref name="levels"/> levels.</summary>
    internal static IEnumerable<int[]> ProperSubsets(int levels)
    {
        for (int mask = 1; mask < (1 << levels) - 1; mask++)
            yield return [.. Enumerable.Range(0, levels).Where(l => (mask >> l & 1) != 0)];
    }

    /// <summary>The subset measure's probability of a subset of <paramref name="size"/> of <paramref name="levels"/>
    /// levels, as D01 states it: R's (the size uniform on 1 … n − 1, then uniform among the subsets of that size), or the
    /// paper's (uniform over the 2ⁿ − 2 proper subsets).</summary>
    internal static double LogSubsetMeasure(CategoricalSubsetMeasure measure, int levels, int size)
        => measure == CategoricalSubsetMeasure.SizeUniform
            ? -Math.Log(levels - 1) - (SpecialFunctions.LogGamma(levels + 1.0) - SpecialFunctions.LogGamma(size + 1.0) - SpecialFunctions.LogGamma(levels - size + 1.0))
            : -Math.Log(Math.Pow(2, levels) - 2);

    /// <summary>
    /// F of <paramref name="model"/>, whose coefficient prior holds τ fixed, at inverse temperature
    /// <paramref name="inverseTemperature"/> on the likelihood path: every structure of at most M_max distinct terms the
    /// design and the target admit, found by extending admitted structures only (F is downward closed, D06), with
    /// log π = ln of the count weight and the terms' measure (the prior's own) plus ℓ_b.
    /// </summary>
    internal static Space Enumerate(BassModel model, double inverseTemperature = 1.0)
    {
        double logTau = Math.Log(model.Target.FixedTau ?? throw new ArgumentException("The enumeration holds τ fixed.", nameof(model)));
        var temperature = new Temperature(inverseTemperature, TemperingPath.Likelihood);
        List<BassTerm> terms = AllTerms(model);
        var ids = new Dictionary<BassTerm, int>();
        for (int k = 0; k < terms.Count; k++) ids.Add(terms[k], k);
        var states = new List<State>();
        void Visit(int from, List<BassTerm> chosen)
        {
            if (BassMixing.Loaded(model, chosen, logTau, out _) is not { } loaded) return;   // no superset is in F
            BassTerm[] held = [.. loaded.Terms];
            BassEvaluation evaluation = model.Target.Evaluate(loaded);
            double logPrior = model.Prior.LogCollectionWeight(held.Length) + held.Sum(model.Prior.LogPrior);
            states.Add(new State(held, [.. held.Select(t => ids[t]).Order()], loaded, BassMixing.Counts(model, loaded),
                                 logPrior + evaluation.LogLikelihood(inverseTemperature), temperature.LogDensity(evaluation)));
            if (chosen.Count == model.MaxBasis) return;
            for (int i = from; i < terms.Count; i++)
            {
                chosen.Add(terms[i]);
                Visit(i + 1, chosen);
                chosen.RemoveAt(chosen.Count - 1);
            }
        }
        Visit(0, []);
        return new Space(model, inverseTemperature, terms, states);
    }

    /// <summary>One enumerated proposal: its kernel's group, its candidate (null for a null move or a refusal, which
    /// stay), its probability given the move, and the production code's within-move log ratio.</summary>
    internal readonly record struct Enumerated(string Group, BassTerm[]? Candidate, double LogProbability, double LogRatio);

    /// <summary>
    /// The transitions K(S, S′), S′ ≠ S, of each kernel group in logs, by (from, to), over the production palette of
    /// <paramref name="space"/>'s model: its moves' weights at each state select the move, and
    /// <paramref name="proposals"/> enumerates each move's proposals by its key.
    /// </summary>
    internal static Dictionary<string, Dictionary<(int From, int To), double>> Kernel(
        Space space, Func<string, State, IEnumerable<Enumerated>> proposals)
    {
        BassModel model = space.Model;
        var palette = new BassStructureMoves(model, [], space.States[0].Loaded.LogTau).Moves;
        Dictionary<string, IRjMove<BassState>> byKey = palette.ToDictionary(m => m.Key);
        var temperature = new Temperature(space.InverseTemperature, TemperingPath.Likelihood);
        double Total(BassState s) => palette.Sum(m => Math.Max(0.0, m.Weight(s)));
        var flows = new Dictionary<string, Dictionary<(int, int), double>>();
        for (int i = 0; i < space.States.Count; i++)
        {
            State from = space.States[i];
            double total = Total(from.Loaded);
            foreach (IRjMove<BassState> move in palette)
            {
                double weight = move.Weight(from.Loaded);
                if (!(weight > 0.0)) continue;
                double logSelect = Math.Log(weight / total);
                foreach (Enumerated e in proposals(move.Key, from))
                {
                    if (e.Candidate is null || double.IsNegativeInfinity(e.LogProbability)) continue;
                    int j = space.Find(e.Candidate);
                    if (j < 0 || j == i) continue;   // refused, or the same state
                    State to = space.States[j];
                    double reverse = byKey[move.ReverseKey].Weight(to.Loaded);
                    Assert.True(reverse > 0.0, $"{move.Key} from {i} to {j}: the reverse has no weight");
                    double logSelection = Math.Log(reverse / Total(to.Loaded)) - logSelect;
                    double r = to.LogTarget - from.LogTarget + e.LogRatio + logSelection;
                    Assert.False(double.IsNaN(r), $"{move.Key} from {i} to {j}: a NaN ratio");
                    double logK = logSelect + e.LogProbability + Math.Min(0.0, r);
                    if (!flows.TryGetValue(e.Group, out Dictionary<(int, int), double>? group)) flows[e.Group] = group = [];
                    group[(i, j)] = group.TryGetValue((i, j), out double held) ? LogAdd(held, logK) : logK;
                }
            }
        }
        return flows;
    }

    /// <summary>
    /// Detailed balance for each group, π(S)K(S, S′) = π(S′)K(S′, S) for every transition with its reverse present, and
    /// πK = π for each group and for their union, all in logs to <see cref="Tolerance"/>. Returns each group's
    /// transitions and the largest log discrepancy.
    /// </summary>
    internal static (Dictionary<string, int> Transitions, double Worst) Check(Space space, Dictionary<string, Dictionary<(int From, int To), double>> flows)
    {
        double worst = 0.0;
        var counts = new Dictionary<string, int>();
        double[] lp = [.. space.States.Select(s => s.LogPi)];
        void Stationary(string name, IEnumerable<KeyValuePair<(int From, int To), double>> transitions)
        {
            var inflow = new double[lp.Length];
            var outflow = new double[lp.Length];
            Array.Fill(inflow, double.NegativeInfinity);
            Array.Fill(outflow, double.NegativeInfinity);
            foreach (((int from, int to), double logK) in transitions)
            {
                inflow[to] = LogAdd(inflow[to], lp[from] + logK);
                outflow[from] = LogAdd(outflow[from], lp[from] + logK);
            }
            for (int j = 0; j < lp.Length; j++)
            {
                if (double.IsNegativeInfinity(inflow[j]) && double.IsNegativeInfinity(outflow[j])) continue;
                double gap = Math.Abs(inflow[j] - outflow[j]);
                worst = Math.Max(worst, gap);
                Assert.True(gap <= Tolerance, $"{name}: πK ≠ π at state {j} (in {inflow[j]:R}, out {outflow[j]:R})");
            }
        }
        foreach ((string name, Dictionary<(int From, int To), double> group) in flows)
        {
            counts[name] = group.Count;
            foreach (((int from, int to), double logK) in group)
            {
                Assert.True(group.TryGetValue((to, from), out double back), $"{name}: {from} → {to} has no reverse");
                double gap = Math.Abs(lp[from] + logK - (lp[to] + back));
                worst = Math.Max(worst, gap);
                Assert.True(gap <= Tolerance, $"{name}: π(S)K(S, S′) ≠ π(S′)K(S′, S) for {from} → {to} by {gap:E2} in logs");
            }
            Stationary(name, group);
        }
        Stationary("the whole sweep", flows.Values.SelectMany(g => g));
        return (counts, worst);
    }

    /// <summary>
    /// D01's birth, death and change as the statement gives them, with the local steps <paramref name="steps"/> inside
    /// the change when on: births from the production proposal's probabilities (D01's, checked against Wallenius and the
    /// enumeration elsewhere), deaths uniform over the terms, and changes enumerated factor by factor.
    /// </summary>
    internal static IEnumerable<Enumerated> D01(Space space, string key, State s, LocalSteps? steps)
    {
        BassModel model = space.Model;
        BassTermProposal proposal = model.Proposal;
        int m = s.Terms.Length;
        switch (key)
        {
            case "birth":
                foreach (BassTerm t in space.Terms)
                {
                    if (s.Terms.Contains(t)) continue;   // a repeated term, refused
                    yield return new Enumerated(BirthDeath, [.. s.Terms, t], proposal.LogBirth(t, s.Counts), proposal.LogBirthRatio(t, s.Counts));
                }
                break;
            case "death":
                foreach (BassTerm t in s.Terms)
                    yield return new Enumerated(BirthDeath, [.. s.Terms.Where(u => !u.Equals(t))], -Math.Log(m), proposal.LogDeathRatio(t, s.Counts));
                break;
            case "change":
                if (m == 0) yield break;
                foreach (Enumerated e in Changes(space, s, steps)) yield return e;
                break;
        }
    }

    // Each term 1/M, each of its factors 1/K; then the redraw with probability 1 − ρ (1 with the steps off), every sign and
    // knot at 1/(2n_v) or every proper subset at the measure's probability, and the local step with probability ρ, every
    // δ ∈ {±1, …, ±w} at 1/(2w) or every level's toggle at 1/n_v.
    private static IEnumerable<Enumerated> Changes(Space space, State s, LocalSteps? steps)
    {
        BassModel model = space.Model;
        BassTermProposal proposal = model.Proposal;
        int m = s.Terms.Length;
        double logRedraw = steps is null ? 0.0 : Math.Log(1.0 - steps.Probability), logStep = steps is null ? double.NegativeInfinity : Math.Log(steps.Probability);
        foreach (BassTerm t in s.Terms)
        {
            BassTerm[] others = [.. s.Terms.Where(u => !u.Equals(t))];
            List<(BassPart Part, int Index)> factors = Factors(t);
            double logPick = -Math.Log(m) - Math.Log(factors.Count);
            for (int pick = 0; pick < factors.Count; pick++)
            {
                (BassPart part, int index) = factors[pick];
                if (part == BassPart.Categorical)
                {
                    CategoricalFactor old = t.Categorical[index];
                    int levels = model.Inputs.LevelCount(old.Variable);
                    if (logRedraw > double.NegativeInfinity)
                        foreach (int[] subset in ProperSubsets(levels))
                        {
                            BassTerm changed = WithSubset(t, index, new CategoricalFactor(old.Variable, subset));
                            if (changed.Equals(t)) continue;   // the same state
                            yield return new Enumerated(Redraw, [.. others, changed],
                                                        logPick + logRedraw + LogSubsetMeasure(model.Prior.Subsets, levels, subset.Length),
                                                        proposal.LogRedrawRatio(t, pick, changed));
                        }
                    if (steps is null) continue;
                    for (int level = 0; level < levels; level++)
                    {
                        int[] toggled = old.Contains(level) ? [.. old.Levels.Where(l => l != level)] : [.. old.Levels, level];
                        BassTermMove move = proposal.Toggled(t, pick, level);
                        Assert.Equal(BassMoveKind.Toggle, move.Kind);
                        if (toggled.Length == 0 || toggled.Length == levels)
                        {
                            Assert.True(move.IsNull, $"toggling level {level} of {old} is a null move");
                            continue;
                        }
                        BassTerm changed = WithSubset(t, index, new CategoricalFactor(old.Variable, toggled));
                        Assert.Equal(changed, move.Term);
                        yield return new Enumerated(Step, [.. others, changed], logPick + logStep - Math.Log(levels), move.LogProposalRatio);
                    }
                    continue;
                }
                HingeFactor h = BassStructurePrior.Hinges(t, part)[index];
                int knots = part == BassPart.Continuous ? model.Inputs.ContinuousKnots(h.Variable).Length : model.Inputs.FunctionalKnots(h.Variable).Length;
                if (logRedraw > double.NegativeInfinity)
                    foreach (int sign in new[] { -1, 1 })
                        for (int knot = 0; knot < knots; knot++)
                        {
                            BassTerm changed = WithHinge(t, part, index, new HingeFactor(h.Variable, sign, knot));
                            if (changed.Equals(t)) continue;
                            yield return new Enumerated(Redraw, [.. others, changed], logPick + logRedraw - Math.Log(2.0 * knots),
                                                        proposal.LogRedrawRatio(t, pick, changed));
                        }
                if (steps is null) continue;
                int w = Math.Max(1, (int)Math.Ceiling((decimal)steps.WindowFraction * knots));
                for (int size = 1; size <= w; size++)
                    foreach (int delta in new[] { -size, size })
                    {
                        BassTermMove move = proposal.Stepped(t, pick, delta);
                        Assert.Equal(BassMoveKind.LocalStep, move.Kind);
                        int knot = h.Knot + delta;
                        if (knot < 0 || knot >= knots)
                        {
                            Assert.True(move.IsNull, $"stepping {h} by {delta} is a null move");
                            continue;
                        }
                        BassTerm changed = WithHinge(t, part, index, h with { Knot = knot });
                        Assert.Equal(changed, move.Term);
                        yield return new Enumerated(Step, [.. others, changed], logPick + logStep - Math.Log(2.0 * w), move.LogProposalRatio);
                    }
            }
        }
    }

    /// <summary>
    /// The fixed-M family as Phase 1's decisions state it, by its keys "prune", "grow" and "variable": a term 1/M, a part
    /// 1/P among the fit's parts; a prune's factor 1/K (an empty part a null move, a term left empty refused); a grow's
    /// variable 1/(p − K) among the unused, at its cap or with none unused a null move; a variable change's factor 1/K and
    /// its new variable 1/(p − K), an empty part or none unused a null move; each new factor from its variable's prior,
    /// a sign and knot at 1/(2n_v) or a subset at the measure's probability. Every term at M = 0 is a null move.
    /// </summary>
    internal static IEnumerable<Enumerated> Family(Space space, string key, State s)
    {
        BassModel model = space.Model;
        BassTermProposal proposal = model.Proposal;
        int m = s.Terms.Length;
        if (m == 0) yield break;
        BassPart[] parts = [.. Enum.GetValues<BassPart>().Where(model.Prior.Has)];
        foreach (BassTerm t in s.Terms)
        {
            BassTerm[] others = [.. s.Terms.Where(u => !u.Equals(t))];
            foreach (BassPart part in parts)
            {
                double logPick = -Math.Log(m) - Math.Log(parts.Length);
                int k = BassStructurePrior.Order(t, part), cap = model.Prior.Cap(part);
                HashSet<int> used = [.. BassTermCounts.Variables(t, part)];
                int[] unused = [.. Enumerable.Range(0, model.Prior.Variables(part)).Where(v => !used.Contains(v))];
                switch (key)
                {
                    case "prune":
                        for (int f = 0; f < k; f++)
                        {
                            BassTermMove move = proposal.Pruned(t, part, f);
                            if (Factors(t).Count == 1)
                            {
                                Assert.True(move.LeavesIntercept, $"pruning {t}'s only factor leaves the intercept");
                                continue;
                            }
                            BassTerm pruned = Removed(t, part, f);
                            Assert.Equal(pruned, move.Term);
                            yield return new Enumerated(PruneGrow, [.. others, pruned], logPick - Math.Log(k), move.LogProposalRatio);
                        }
                        break;
                    case "grow":
                        if (k >= cap || unused.Length == 0) break;
                        foreach (int v in unused)
                            foreach ((HingeFactor hinge, CategoricalFactor? subset, double logFactor) in FactorsOf(model, part, v))
                            {
                                BassTermMove move = subset is null ? proposal.Grown(t, part, hinge) : proposal.Grown(t, part, subset);
                                BassTerm grown = Added(t, part, hinge, subset);
                                Assert.Equal(grown, move.Term);
                                yield return new Enumerated(PruneGrow, [.. others, grown], logPick - Math.Log(unused.Length) + logFactor, move.LogProposalRatio);
                            }
                        break;
                    case "variable":
                        if (k == 0 || unused.Length == 0) break;
                        for (int f = 0; f < k; f++)
                            foreach (int v in unused)
                                foreach ((HingeFactor hinge, CategoricalFactor? subset, double logFactor) in FactorsOf(model, part, v))
                                {
                                    BassTermMove move = subset is null ? proposal.Exchanged(t, part, f, hinge) : proposal.Exchanged(t, part, f, subset);
                                    BassTerm changed = Replaced(t, part, f, hinge, subset);
                                    Assert.Equal(changed, move.Term);
                                    yield return new Enumerated(VariableChange, [.. others, changed],
                                                                logPick - Math.Log(k) - Math.Log(unused.Length) + logFactor, move.LogProposalRatio);
                                }
                        break;
                }
            }
        }
    }

    /// <summary>D01's moves with the local steps <paramref name="steps"/> when on, and the fixed-M family by its keys.</summary>
    internal static IEnumerable<Enumerated> All(Space space, string key, State s, LocalSteps? steps)
        => key is "prune" or "grow" or "variable" ? Family(space, key, s) : D01(space, key, s, steps);

    /// <summary>Each factor variable <paramref name="v"/> of <paramref name="part"/> can take, with its prior's log
    /// probability as D01 states it: every sign and knot at 1/(2n_v), or every proper subset at the measure's.</summary>
    internal static IEnumerable<(HingeFactor Hinge, CategoricalFactor? Subset, double LogProbability)> FactorsOf(BassModel model, BassPart part, int v)
    {
        if (part == BassPart.Categorical)
        {
            int levels = model.Inputs.LevelCount(v);
            foreach (int[] subset in ProperSubsets(levels))
                yield return (default, new CategoricalFactor(v, subset), LogSubsetMeasure(model.Prior.Subsets, levels, subset.Length));
            yield break;
        }
        int knots = part == BassPart.Continuous ? model.Inputs.ContinuousKnots(v).Length : model.Inputs.FunctionalKnots(v).Length;
        foreach (int sign in new[] { -1, 1 })
            for (int knot = 0; knot < knots; knot++) yield return (new HingeFactor(v, sign, knot), null, -Math.Log(2.0 * knots));
    }

    /// <summary><paramref name="t"/> without factor <paramref name="index"/> of <paramref name="part"/>; it must keep a factor.</summary>
    internal static BassTerm Removed(BassTerm t, BassPart part, int index) => part switch
    {
        BassPart.Continuous => new BassTerm(t.Continuous.Where((_, i) => i != index), t.Categorical, t.Functional),
        BassPart.Categorical => new BassTerm(t.Continuous, t.Categorical.Where((_, i) => i != index), t.Functional),
        _ => new BassTerm(t.Continuous, t.Categorical, t.Functional.Where((_, i) => i != index)),
    };

    /// <summary><paramref name="t"/> with a factor added to <paramref name="part"/>: <paramref name="subset"/> for the
    /// categorical part, else <paramref name="hinge"/>.</summary>
    internal static BassTerm Added(BassTerm t, BassPart part, HingeFactor hinge, CategoricalFactor? subset) => part switch
    {
        BassPart.Continuous => new BassTerm([.. t.Continuous, hinge], t.Categorical, t.Functional),
        BassPart.Categorical => new BassTerm(t.Continuous, [.. t.Categorical, subset!], t.Functional),
        _ => new BassTerm(t.Continuous, t.Categorical, [.. t.Functional, hinge]),
    };

    /// <summary><paramref name="t"/> with factor <paramref name="index"/> of <paramref name="part"/> replaced:
    /// <paramref name="subset"/> for the categorical part, else <paramref name="hinge"/>.</summary>
    internal static BassTerm Replaced(BassTerm t, BassPart part, int index, HingeFactor hinge, CategoricalFactor? subset) => part switch
    {
        BassPart.Categorical => WithSubset(t, index, subset!),
        _ => WithHinge(t, part, index, hinge),
    };

    /// <summary>A term's factors as the change counts them: continuous, categorical, then functional.</summary>
    internal static List<(BassPart Part, int Index)> Factors(BassTerm t)
        => [.. Enumerable.Range(0, t.Continuous.Count).Select(i => (BassPart.Continuous, i)),
            .. Enumerable.Range(0, t.Categorical.Count).Select(i => (BassPart.Categorical, i)),
            .. Enumerable.Range(0, t.Functional.Count).Select(i => (BassPart.Functional, i))];

    internal static BassTerm WithHinge(BassTerm t, BassPart part, int index, HingeFactor h)
    {
        HingeFactor[] hinges = [.. BassStructurePrior.Hinges(t, part)];
        hinges[index] = h;
        return part == BassPart.Continuous ? new BassTerm(hinges, t.Categorical, t.Functional) : new BassTerm(t.Continuous, t.Categorical, hinges);
    }

    internal static BassTerm WithSubset(BassTerm t, int index, CategoricalFactor c)
    {
        CategoricalFactor[] subsets = [.. t.Categorical];
        subsets[index] = c;
        return new BassTerm(t.Continuous, subsets, t.Functional);
    }

    private static double LogAdd(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double top = Math.Max(a, b);
        return top + Math.Log(Math.Exp(a - top) + Math.Exp(b - top));
    }
}
