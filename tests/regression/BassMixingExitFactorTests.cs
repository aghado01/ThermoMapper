using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's P2 (the mode-exit pricing brief): the design draft's structural moves at fixed M, scored as test-side proposals
/// from each chain's fixed state (<see cref="BassMixing.FixedStates"/>) through the chain's design, none committed.
/// <list type="bullet">
///   <item><b>Prune:</b> every factor of every term whose part holds K ≥ 2, dropped (K → K − 1). On a fit of several parts,
///   also each part's only factor where the term keeps another part (K = 1 → 0), the reverse of a grow into an empty part,
///   tagged apart.</item>
///   <item><b>Grow:</b> 5,000 per state, the term uniformly, its part uniformly among the fit's, the variable uniformly among
///   the part's unused ones, the sign and knot (or the subset) from the prior; at the part's cap, a refusal.</item>
///   <item><b>Variable change:</b> 5,000 per state, the term and then its factor uniformly (as today's change picks them),
///   v → v′ uniformly among the part's unused variables, the new factor from v′'s prior.</item>
/// </list>
/// Each is scored by its Δℓ and its full log ratio, Δℓ plus the target's structure prior ratio plus the proposal ratio,
/// which the draft derives to cancel (<see cref="BassMixing.LogPruneRatio"/>, <see cref="BassMixing.LogGrowRatio"/>); the
/// choice of the move and of the part is taken as matched. The two-step profile: from each prune or change costing under
/// 20 nats, the best single and pair death of the state it reaches. It reads, it does not judge; it asserts its own
/// soundness: nothing committed (<see cref="BassMixing.Untouched"/>), the cancellation to 1e-12 relative (printed before
/// it is asserted), and the counts adding up.
/// </summary>
public sealed class BassMixingExitFactorTests
{
    private const int Proposals = 5_000;
    private const double Cheap = -20.0, Tolerance = 1e-12;
    private readonly ITestOutputHelper _out;
    public BassMixingExitFactorTests(ITestOutputHelper output) => _out = output;

    // One proposal: its kind ("prune", "prune, part emptied", "grow", "change"), part, activity tag (Friedman), orders before
    // and after, the term it changes and the term it proposes (null when unavailable), its score, the candidate state, and
    // the cancellation's residuals, relative (the term's and the target's), NaN where not scored.
    private sealed record Proposed(string Kind, BassPart Part, string Activity, string Orders, BassTerm Old, BassTerm? New, BassMixing.Exit Exit,
                                   BassState? Candidate, double TermResidual, double TargetResidual)
    {
        public bool Available => New is not null;
    }

    // A cheap re-representation's next step: the state it reaches, and that state's best single and pair death among the
    // other terms (removing the re-represented term reaches the fixed state's own death of the old one), each with the same
    // removal's Δℓ at the fixed state.
    private sealed record TwoStep(Proposed Step, BassMixing.Removal Single, double SingleBefore, BassMixing.Removal? Pair, double PairBefore);

    // A state's two-step profile, beside its own best single and pair death (P1's).
    private sealed record TwoSteps(TwoStep[] Steps, double FixedSingle, double FixedPair);

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void FactorMoves_FromFixedStates(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassMixing.Fixture fixture = BassMixing.Get(name);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassMixing.FixedState[] states = BassMixing.FixedStates(model);
        double setup = clock.Elapsed.TotalSeconds;
        clock.Restart();
        var moves = new Proposed[states.Length][];
        var twoSteps = new TwoSteps[states.Length];
        int[] seeds = SeedTree.Derive(113, 2 * states.Length);
        Parallel.For(0, states.Length, c =>
        {
            BassChain chain = states[c].Chain;
            var untouched = new BassMixing.Untouched(model, chain);
            moves[c] = [.. Prunes(fixture, model, chain), .. Grows(fixture, model, chain, new Xoshiro256PlusPlus(seeds[2 * c])),
                        .. Changes(fixture, model, chain, new Xoshiro256PlusPlus(seeds[2 * c + 1]))];
            untouched.Check();
            twoSteps[c] = Profile(model, chain.Current, moves[c]);
            untouched.Check();
        });
        _out.WriteLine($"## P2, {name}: fixed-M factor moves from each chain's fixed state (states and phases {setup:F1} s, moves {clock.Elapsed.TotalSeconds:F1} s)");
        _out.WriteLine("");
        foreach (string line in BassMixing.FixedTable(model, states)) _out.WriteLine(line);

        // The counts: every prune enumerated once, and Proposals grows and changes per state.
        for (int c = 0; c < states.Length; c++)
        {
            BassState state = states[c].Chain.Current;
            int prunable = state.Terms.Sum(t => BassMixing.Factors(t).Count(f => Prunable(model, t, f)));
            Assert.Equal(prunable, moves[c].Count(m => m.Kind.StartsWith("prune", StringComparison.Ordinal)));
            Assert.Equal(Proposals, moves[c].Count(m => m.Kind == "grow"));
            Assert.Equal(Proposals, moves[c].Count(m => m.Kind == "change"));
        }
        Cancellation(moves);
        PerState(states, moves);
        Proposed[] pooled = [.. moves.SelectMany(m => m)];
        Groups("by kind and part", pooled, m => $"{m.Kind}, {m.Part.ToString().ToLowerInvariant()}");
        if (fixture.Inert is not null) Groups("by the factor's variable, active or inert", pooled, m => $"{m.Kind}, {m.Activity}");
        Groups("by the term's orders before and after", pooled, m => $"{m.Kind}, {m.Orders}");
        for (int c = 0; c < states.Length; c++) Cheapest(fixture, model, states[c], moves[c]);
        Profile(fixture, model, states, twoSteps);

        // The derivation's check, asserted after its residuals are printed.
        Assert.All(pooled.Where(m => !double.IsNaN(m.TermResidual)), m => Assert.True(m.TermResidual <= Tolerance, $"term-prior residual {m.TermResidual:G3}"));
        Assert.All(pooled.Where(m => !double.IsNaN(m.TargetResidual)), m => Assert.True(m.TargetResidual <= Tolerance, $"target residual {m.TargetResidual:G3}"));
    }

    // Whether factor f of the term is pruned: its part holds K ≥ 2, or (several parts) it is its part's only factor and the
    // term keeps another part.
    private static bool Prunable(BassModel model, BassTerm term, BassMixing.FactorAt f)
        => BassMixing.Order(term, f.Part) >= 2 || model.Prior.MinimumOrder(f.Part) == 0 && BassMixing.Without(term, f) is not null;

    private static IEnumerable<Proposed> Prunes(BassMixing.Fixture fixture, BassModel model, BassChain chain)
    {
        BassState state = chain.Current;
        BassDesign design = chain.Structure.DesignAt(state);
        for (int i = 0; i < state.Count; i++)
        {
            BassTerm term = state.Terms[i];
            foreach (BassMixing.FactorAt f in BassMixing.Factors(term).Where(f => Prunable(model, term, f)))
            {
                BassTerm pruned = BassMixing.Without(term, f)!;
                double ratio = BassMixing.LogPruneRatio(model, term, f);
                string kind = BassMixing.Order(term, f.Part) >= 2 ? "prune" : "prune, part emptied";
                yield return Score(fixture, model, design, state, i, kind, f.Part, Activity(fixture, model, f.Part, BassMixing.Variable(term, f)), term, pruned, ratio);
            }
        }
    }

    private static IEnumerable<Proposed> Grows(BassMixing.Fixture fixture, BassModel model, BassChain chain, Xoshiro256PlusPlus rng)
    {
        BassState state = chain.Current;
        BassDesign design = chain.Structure.DesignAt(state);
        BassPart[] parts = [.. Enum.GetValues<BassPart>().Where(model.Prior.Has)];
        for (int k = 0; k < Proposals; k++)
        {
            int i = BassDraws.Below(rng, state.Count);
            BassTerm term = state.Terms[i];
            BassPart part = parts[BassDraws.Below(rng, parts.Length)];
            if (BassMixing.Order(term, part) >= model.Prior.Cap(part))
            {
                yield return new Proposed("grow", part, "at the cap", Orders(term, null), term, term, BassMixing.Exit.Refused(BassRefusal.Cap, false), null, double.NaN, double.NaN);
                continue;
            }
            int[] unused = Unused(model, term, part);
            int v = unused[BassDraws.Below(rng, unused.Length)];
            double logFactor = BassMixing.DrawFactor(model, rng, part, v, out HingeFactor hinge, out CategoricalFactor? subset);
            BassTerm grown = BassMixing.With(term, part, hinge, subset);
            double ratio = BassMixing.LogGrowRatio(model, term, part, logFactor);
            yield return Score(fixture, model, design, state, i, "grow", part, Activity(fixture, model, part, v), term, grown, ratio);
        }
    }

    private static IEnumerable<Proposed> Changes(BassMixing.Fixture fixture, BassModel model, BassChain chain, Xoshiro256PlusPlus rng)
    {
        BassState state = chain.Current;
        BassDesign design = chain.Structure.DesignAt(state);
        for (int k = 0; k < Proposals; k++)
        {
            int i = BassDraws.Below(rng, state.Count);
            BassTerm term = state.Terms[i];
            BassMixing.FactorAt[] factors = [.. BassMixing.Factors(term)];
            BassMixing.FactorAt f = factors[BassDraws.Below(rng, factors.Length)];
            int[] unused = Unused(model, term, f.Part);
            int old = BassMixing.Variable(term, f);
            if (unused.Length == 0)
            {
                yield return new Proposed("change", f.Part, "no unused variable", Orders(term, null), term, null, default, null, double.NaN, double.NaN);
                continue;
            }
            int v = unused[BassDraws.Below(rng, unused.Length)];
            double logFactor = BassMixing.DrawFactor(model, rng, f.Part, v, out HingeFactor hinge, out CategoricalFactor? subset);
            BassTerm changed = BassMixing.Replaced(term, f, hinge, subset);
            double ratio = BassMixing.LogFactorPrior(model, term, f) - logFactor;
            string activity = $"{Activity(fixture, model, f.Part, old)} → {Activity(fixture, model, f.Part, v)}";
            yield return Score(fixture, model, design, state, i, "change", f.Part, activity, term, changed, ratio);
        }
    }

    // A proposal scored through the chain's design; the cancellation's residuals, relative to the log priors it differences.
    private static Proposed Score(BassMixing.Fixture fixture, BassModel model, BassDesign design, BassState state, int index, string kind, BassPart part,
                                  string activity, BassTerm old, BassTerm proposed, double logProposalRatio)
    {
        double before = model.Prior.LogPrior(old), after = model.Prior.LogPrior(proposed);
        double termResidual = Math.Abs(after - before + logProposalRatio) / (1.0 + Math.Abs(before) + Math.Abs(after));
        (BassMixing.Exit exit, BassState? candidate) = BassMixing.Replace(model, design, state, index, proposed, logProposalRatio);
        double targetResidual = double.NaN;
        if (candidate is not null)
        {
            double from = model.Target.LogStructurePrior(state.Fit), to = model.Target.LogStructurePrior(candidate.Fit);
            targetResidual = Math.Abs(exit.LogPriorRatio + exit.LogProposalRatio) / (1.0 + Math.Abs(from) + Math.Abs(to));
        }
        return new Proposed(kind, part, activity, Orders(old, proposed), old, proposed, exit, candidate, termResidual, targetResidual);
    }

    private static int[] Unused(BassModel model, BassTerm term, BassPart part)
    {
        var used = new HashSet<int>(BassMixing.Factors(term).Where(f => f.Part == part).Select(f => BassMixing.Variable(term, f)));
        return [.. Enumerable.Range(0, model.Prior.Variables(part)).Where(v => !used.Contains(v))];
    }

    // On Friedman, whether a continuous variable is inert (x6–x10); elsewhere the part.
    private static string Activity(BassMixing.Fixture fixture, BassModel model, BassPart part, int variable)
        => fixture.Inert is null || part != BassPart.Continuous ? part.ToString().ToLowerInvariant() : variable >= 5 ? "inert" : "active";

    private static string Orders(BassTerm old, BassTerm? proposed) => proposed is null || proposed.Equals(old)
        ? BassMixing.Orders(old)
        : $"{BassMixing.Orders(old)} → {BassMixing.Orders(proposed)}";

    // The two-step profile: each distinct state reached by a prune or change costing under 20 nats, and its best removals of
    // the other terms. The design moves the re-represented term last, so the others keep their terms.
    private static TwoSteps Profile(BassModel model, BassState state, Proposed[] moves)
    {
        var before = new Dictionary<string, double>();
        BassMixing.Removal[] singles = BassMixing.Removals(model, state, 1), pairs = BassMixing.Removals(model, state, 2);
        foreach (BassMixing.Removal r in singles.Concat(pairs))
            before[BassMixing.Key(r.Terms.Select(i => state.Terms[i]))] = r.Exit.LogLikelihoodChange;
        var seen = new HashSet<string>();
        var profile = new List<TwoStep>();
        foreach (Proposed m in moves.Where(m => m.Kind != "grow" && m.Available && m.Exit.Admitted && m.Exit.LogLikelihoodChange > Cheap)
                                    .OrderByDescending(m => m.Exit.LogLikelihoodChange))
        {
            BassState reached = m.Candidate!;
            if (!seen.Add(BassMixing.Key(reached.Terms))) continue;
            int changed = reached.Count - 1;
            Assert.Equal(m.New, reached.Terms[changed]);
            string Key(BassMixing.Removal r) => BassMixing.Key(r.Terms.Select(i => reached.Terms[i]));
            BassMixing.Removal single = BassMixing.Removals(model, reached, 1).Where(r => !r.Terms.Contains(changed)).MaxBy(r => r.Exit.LogLikelihoodChange)!;
            BassMixing.Removal? pair = reached.Count < 3 ? null
                : BassMixing.Removals(model, reached, 2).Where(r => !r.Terms.Contains(changed)).MaxBy(r => r.Exit.LogLikelihoodChange);
            profile.Add(new TwoStep(m, single, before[Key(single)], pair, pair is null ? double.NaN : before[Key(pair)]));
        }
        return new TwoSteps([.. profile], singles.Max(r => r.Exit.LogLikelihoodChange), pairs.Max(r => r.Exit.LogLikelihoodChange));
    }

    private void Cancellation(Proposed[][] moves)
    {
        Proposed[] pooled = [.. moves.SelectMany(m => m)];
        double[] term = [.. pooled.Select(m => m.TermResidual).Where(r => !double.IsNaN(r))];
        double[] target = [.. pooled.Select(m => m.TargetResidual).Where(r => !double.IsNaN(r))];
        _out.WriteLine("");
        _out.WriteLine($"The cancellation (prior ratio + proposal ratio, relative to the log priors it differences): the term prior's over {term.Length:N0} proposals, " +
                       $"largest {term.Max():G3}, over {Tolerance:G0} {term.Count(r => r > Tolerance)}; the target's over {target.Length:N0} admitted, largest {target.Max():G3}, over {Tolerance:G0} {target.Count(r => r > Tolerance)}");
    }

    // Per state: each kind's proposals, unavailable and refused shares, mean α and its largest, and the cheapest Δℓ.
    private void PerState(BassMixing.FixedState[] states, Proposed[][] moves)
    {
        _out.WriteLine("");
        _out.WriteLine("chain | M | kind | proposals | unavailable | refused | mean α | largest α | α > 0.01 | largest Δℓ");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < states.Length; c++)
            foreach (IGrouping<string, Proposed> g in moves[c].GroupBy(m => m.Kind))
            {
                Proposed[] set = [.. g];
                Proposed[] admitted = [.. set.Where(m => m.Available && m.Exit.Admitted)];
                _out.WriteLine($"{states[c].Label} | {states[c].Chain.Current.Count} | {g.Key} | {set.Length} | {Share(set.Count(m => !m.Available), set.Length)} | " +
                               $"{Share(set.Count(m => m.Available && !m.Exit.Admitted), set.Length)} | {G(set.Average(m => m.Available ? m.Exit.Acceptance : 0.0))} | " +
                               $"{(admitted.Length > 0 ? G(admitted.Max(m => m.Exit.Acceptance)) : "—")} | {admitted.Count(m => m.Exit.Acceptance > 0.01)} | " +
                               $"{(admitted.Length > 0 ? F(admitted.Max(m => m.Exit.LogLikelihoodChange)) : "—")}");
            }
    }

    // Pooled over the states, by group: proposals, unavailable, refusals by reason, mean α, and Δℓ's quantiles.
    private void Groups(string title, Proposed[] pooled, Func<Proposed, string> key)
    {
        _out.WriteLine("");
        _out.WriteLine($"Pooled, {title} | proposals | unavailable | refused: support | repeated | angle | cap | mean α | Δℓ 5% | 25% | median | 75% | 95% | max");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        foreach (IGrouping<string, Proposed> g in pooled.GroupBy(key).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            Proposed[] set = [.. g];
            int n = set.Length;
            bool Refused(Proposed m, BassRefusal r) => m.Available && m.Exit.Refusal == r;
            double[] dl = [.. set.Where(m => m.Available && m.Exit.Admitted).Select(m => m.Exit.LogLikelihoodChange)];
            _out.WriteLine($"{g.Key} | {n} | {Share(set.Count(m => !m.Available), n)} | {Share(set.Count(m => Refused(m, BassRefusal.Support)), n)} | " +
                           $"{Share(set.Count(m => Refused(m, BassRefusal.Numerical) && m.Exit.Repeated), n)} | {Share(set.Count(m => Refused(m, BassRefusal.Numerical) && !m.Exit.Repeated), n)} | " +
                           $"{Share(set.Count(m => Refused(m, BassRefusal.Cap)), n)} | {G(set.Average(m => m.Available ? m.Exit.Acceptance : 0.0))} | " +
                           string.Join(" | ", new[] { 0.05, 0.25, 0.5, 0.75, 0.95 }.Select(q => F(BassMixing.Quantile(dl, q)))) + $" | {(dl.Length > 0 ? F(dl.Max()) : "—")}");
        }
    }

    // A state's ten cheapest exits by family (the kind, the old term, the new term's variables), each by its best member,
    // with the family's proposals and their mean α.
    private void Cheapest(BassMixing.Fixture fixture, BassModel model, BassMixing.FixedState state, Proposed[] moves)
    {
        _out.WriteLine("");
        _out.WriteLine($"Chain {state.Label}, M {state.Chain.Current.Count}: the ten cheapest exits, by family | kind | Δℓ | ratio | α | the family's proposals | their mean α");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | ---");
        var families = moves.Where(m => m.Available).GroupBy(m => (m.Kind, m.Old, Variables: fixture.Variables(m.New!)))
                            .Where(g => g.Any(m => m.Exit.Admitted))
                            .Select(g => (Best: g.Where(m => m.Exit.Admitted).MaxBy(m => m.Exit.Ratio)!, Count: g.Count(), Mean: g.Average(m => m.Exit.Acceptance)));
        foreach ((Proposed m, int count, double mean) in families.OrderByDescending(f => f.Best.Exit.Ratio).Take(10))
            _out.WriteLine($"{fixture.Describe(model, m.Old)} → {fixture.Describe(model, m.New!)} | {m.Kind} | {F(m.Exit.LogLikelihoodChange)} | {F(m.Exit.Ratio)} | " +
                           $"{G(m.Exit.Acceptance)} | {count} | {G(mean)}");
    }

    // The two-step profile: per state, the cheap re-representations reached and whether they open a cheap removal; then each
    // state's best two-step paths by their total Δℓ.
    private void Profile(BassMixing.Fixture fixture, BassModel model, BassMixing.FixedState[] states, TwoSteps[] twoSteps)
    {
        static double Total(TwoStep t, BassMixing.Removal? r) => r is null ? double.NaN : t.Step.Exit.LogLikelihoodChange + r.Exit.LogLikelihoodChange;
        _out.WriteLine("");
        _out.WriteLine($"Two-step profile, the other terms' removals: chain | M | distinct states reached under {-Cheap} nats (prunes, changes) | " +
                       "best single there | the same at the fixed state | best total, step + single | the fixed state's best single | best pair there | the same at the fixed state | " +
                       "best total, step + pair | the fixed state's best pair | states whose best single, pair is cheaper than at the fixed state | totals over −10 (single, pair) | over −20");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < states.Length; c++)
        {
            TwoStep[] p = twoSteps[c].Steps;
            string count = $"{p.Length} ({p.Count(t => t.Step.Kind.StartsWith("prune", StringComparison.Ordinal))}, {p.Count(t => t.Step.Kind == "change")})";
            if (p.Length == 0)
            {
                _out.WriteLine($"{states[c].Label} | {states[c].Chain.Current.Count} | 0 | — | — | — | {F(twoSteps[c].FixedSingle)} | — | — | — | {F(twoSteps[c].FixedPair)} | — | — | —");
                continue;
            }
            TwoStep single = p.MaxBy(t => t.Single.Exit.LogLikelihoodChange)!, pair = p.Where(t => t.Pair is not null).MaxBy(t => t.Pair!.Exit.LogLikelihoodChange)!;
            _out.WriteLine($"{states[c].Label} | {states[c].Chain.Current.Count} | {count} | {F(single.Single.Exit.LogLikelihoodChange)} | {F(single.SingleBefore)} | " +
                           $"{F(p.Max(t => Total(t, t.Single)))} | {F(twoSteps[c].FixedSingle)} | {F(pair.Pair!.Exit.LogLikelihoodChange)} | {F(pair.PairBefore)} | " +
                           $"{F(p.Where(t => t.Pair is not null).Max(t => Total(t, t.Pair)))} | {F(twoSteps[c].FixedPair)} | " +
                           $"{p.Count(t => t.Single.Exit.LogLikelihoodChange > t.SingleBefore)}, {p.Count(t => t.Pair is not null && t.Pair.Exit.LogLikelihoodChange > t.PairBefore)} | " +
                           $"{p.Count(t => Total(t, t.Single) > -10)}, {p.Count(t => Total(t, t.Pair) > -10)} | {p.Count(t => Total(t, t.Single) > Cheap)}, {p.Count(t => Total(t, t.Pair) > Cheap)}");
        }
        for (int c = 0; c < states.Length; c++)
        {
            if (twoSteps[c].Steps.Length == 0) continue;
            _out.WriteLine("");
            _out.WriteLine($"Chain {states[c].Label}: best two-step paths | step Δℓ | then the death | its Δℓ | the same at the fixed state | total Δℓ | the death's ratio | " +
                           "the end state's log posterior less the fixed state's");
            _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
            double start = BassMixing.LogPosterior(model, states[c].Chain.Current);
            var paths = twoSteps[c].Steps.SelectMany(t => new[] { (Step: t, Removal: (BassMixing.Removal?)t.Single, Before: t.SingleBefore), (Step: t, Removal: t.Pair, Before: t.PairBefore) })
                                   .Where(x => x.Removal is not null);
            foreach ((TwoStep t, BassMixing.Removal? r, double before) in paths.OrderByDescending(x => Total(x.Step, x.Removal)).Take(5))
            {
                BassState reached = t.Step.Candidate!;
                _out.WriteLine($"{t.Step.Kind}: {fixture.Describe(model, t.Step.Old)} → {fixture.Describe(model, t.Step.New!)} | {F(t.Step.Exit.LogLikelihoodChange)} | " +
                               $"{string.Join(" + ", r!.Terms.Select(i => fixture.Describe(model, reached.Terms[i])))} | {F(r.Exit.LogLikelihoodChange)} | {F(before)} | " +
                               $"{F(Total(t, r))} | {F(r.Exit.Ratio)} | {F(BassMixing.LogPosterior(model, r.Reduced!) - start)}");
            }
        }
    }

    private static string Share(double part, double whole) => whole > 0 ? $"{100 * part / whole:F1}%" : "—";
    private static string F(double v) => double.IsNaN(v) ? "—" : v.ToString("F2");
    private static string G(double v) => double.IsNaN(v) ? "—" : v.ToString("G3");
}
