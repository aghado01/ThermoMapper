using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's M2 (the local-move measurements brief): which births fail, and why. From each chain's state at the end of S10's
/// burn-in, on each of its fixtures, 10,000 births are proposed through the chain's own birth move and scored with the
/// target, none committed. Each is an empty draw, a refusal by reason, or its log ratio split into Δℓ, the count prior's
/// ratio, the term prior's, the within-move ratio and the selection ratio (<see cref="BassMixing.Score"/>). Each proposed
/// term is tagged by its orders, on Friedman's surface active or inert, with its knots' quantiles and signs. The deaths:
/// each term of the state removed in turn, its Δℓ and its full ratio. It reads, it does not judge; it asserts its own
/// soundness: nothing committed (<see cref="BassMixing.Untouched"/>), the parts summing to the chain's ratio, the
/// deaths' accounts matching the death move's own, and the counts adding up and matching the moves' tallies.
/// </summary>
public sealed class BassMixingBirthTests
{
    private const int Proposals = 10_000;
    private readonly ITestOutputHelper _out;
    public BassMixingBirthTests(ITestOutputHelper output) => _out = output;

    // One proposed birth: empty (no term), or its term, tag, hinge knots' quantiles and signs, and its score.
    private sealed record Birth(BassTerm? Term, string Tag, bool? Inert, BassMixing.Scored Score);

    // One death of a state's term: its term, tag and score.
    private sealed record Death(BassTerm Term, string Tag, BassMixing.Scored Score);

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void BirthsAndDeaths_FromBurnedStates(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassMixing.Fixture fixture = BassMixing.Get(name);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassChain[] chains = BassMixing.Burned(model);
        var births = new Birth[chains.Length][];
        var deaths = new Death[chains.Length][];
        int[] seeds = SeedTree.Derive(107, chains.Length);
        Parallel.For(0, chains.Length, c =>
        {
            var untouched = new BassMixing.Untouched(model, chains[c]);
            births[c] = Births(fixture, model, chains[c], new Xoshiro256PlusPlus(seeds[c]));
            untouched.Check();
            deaths[c] = Deaths(fixture, model, chains[c], new Xoshiro256PlusPlus(seeds[c]));
            untouched.Check();
        });
        _out.WriteLine($"## M2, {name}: {Proposals:N0} births from each chain's state after {BassMixing.BurnIn:N0} sweeps ({clock.Elapsed.TotalSeconds:F1} s)");
        Losses(chains, births);
        ByTag(births);
        Quantiles(fixture, births);
        Best(fixture, model, chains, births);
        MainEffects(fixture, model, chains, births);
        Knots(fixture, model, births);
        DeathTable(fixture, model, deaths);
    }

    private static Birth[] Births(BassMixing.Fixture fixture, BassModel model, BassChain chain, Xoshiro256PlusPlus rng)
    {
        BassState state = chain.Current;
        IReadOnlyList<IRjMove<BassState>> moves = chain.Structure.Moves;
        IRjMove<BassState> birth = BassMixing.Move(moves, "birth");
        BassTermCounts counts = BassMixing.Counts(model, state);
        chain.Structure.ResetTallies();
        var births = new Birth[Proposals];
        for (int k = 0; k < Proposals; k++)
        {
            if (birth.Propose(state, rng) is not { } p)
            {
                births[k] = new Birth(null, "empty", null, default);
                continue;
            }
            BassTerm term = p.Candidate.Terms[^1];
            Assert.Equal(state.Count + 1, p.Candidate.Count);
            Assert.True(p.Candidate.Terms.Take(state.Count).SequenceEqual(state.Terms));
            double within = -Math.Log(state.Count + 1.0) - model.Proposal.LogBirth(term, counts);
            Assert.Equal(BassMixing.Bits(within), BassMixing.Bits(p.LogProposalRatio));
            BassMixing.Scored score = BassMixing.Score(model, moves, birth, state, p, model.Prior.LogPrior(term), within,
                                                       repeated: state.Terms.Contains(term));
            births[k] = new Birth(term, fixture.Tag(term), fixture.Inert?.Invoke(term), score);
        }
        // The moves' own tallies count what the harness counted.
        Assert.Equal(births.LongCount(b => b.Term is null), chain.Structure.EmptyBirths);
        foreach (BassRefusal reason in Enum.GetValues<BassRefusal>())
            Assert.Equal(births.LongCount(b => b.Term is not null && b.Score.Refusal == reason), chain.Structure.Refusals.Count(reason));
        return births;
    }

    // Each term of the state removed in turn, through the design, its account checked against the death move's own.
    private static Death[] Deaths(BassMixing.Fixture fixture, BassModel model, BassChain chain, Xoshiro256PlusPlus rng)
    {
        BassState state = chain.Current;
        IReadOnlyList<IRjMove<BassState>> moves = chain.Structure.Moves;
        IRjMove<BassState> death = BassMixing.Move(moves, "death");
        BassTermCounts counts = BassMixing.Counts(model, state);
        int m = state.Count;
        var deaths = new Death[m];
        var proposals = new Proposal<BassState>[m];
        for (int i = 0; i < m; i++)
        {
            BassTerm term = state.Terms[i];
            DesignProposal proposal = chain.Structure.DesignAt(state).ProposeDeath(i);
            BassTerm[] rest = [.. state.Terms.Take(i), .. state.Terms.Skip(i + 1)];
            double logTermPrior = 0.0;
            foreach (BassTerm t in rest) logTermPrior += model.Prior.LogPrior(t);
            var candidate = new BassState(rest, BassStructureFit.Of(proposal, logTermPrior), state.LogTau);
            // ln q_b(term | the structure without it) + ln M, from counts rebuilt without it.
            var without = new BassTermCounts(model.Prior);
            foreach (BassTerm t in rest) without.Add(t);
            double within = model.Proposal.LogBirth(term, without) + Math.Log(m);
            Assert.Equal(BassMixing.Bits(model.Proposal.LogDeathRatio(term, counts)), BassMixing.Bits(within));
            proposals[i] = new Proposal<BassState>(candidate, within, 0.0);
            deaths[i] = new Death(term, fixture.Tag(term), BassMixing.Score(model, moves, death, state, proposals[i], -model.Prior.LogPrior(term), within));
        }
        // The death move's own proposals, each term it picks, match the harness's account of that term bit for bit.
        var seen = new bool[m];
        for (int k = 0; k < 40 * m; k++)
        {
            if (death.Propose(state, rng) is not { } p) continue;
            int i = 0;
            while (i < m - 1 && ReferenceEquals(p.Candidate.Terms[i], state.Terms[i])) i++;
            seen[i] = true;
            Assert.Equal(proposals[i].Candidate.Terms, p.Candidate.Terms);
            Assert.Equal(proposals[i].Candidate.Fit, p.Candidate.Fit);
            Assert.Equal(BassMixing.Bits(proposals[i].LogProposalRatio), BassMixing.Bits(p.LogProposalRatio));
        }
        Assert.All(seen, s => Assert.True(s));
        return deaths;
    }

    // The share of births lost to each cause, per chain and pooled: they add up to one.
    private void Losses(BassChain[] chains, Birth[][] births)
    {
        (string Cause, Func<Birth, double> Weight)[] causes =
        [
            ("empty draw", b => b.Term is null ? 1 : 0),
            ("refused: support", b => Refused(b, BassRefusal.Support) ? 1 : 0),
            ("refused: a repeated term", b => Refused(b, BassRefusal.Numerical) && b.Score.Repeated ? 1 : 0),
            ("refused: angle", b => Refused(b, BassRefusal.Numerical) && !b.Score.Repeated ? 1 : 0),
            ("refused: cap", b => Refused(b, BassRefusal.Cap) ? 1 : 0),
            ("refused: variance inflation", b => Refused(b, BassRefusal.VarianceInflation) ? 1 : 0),
            ("refused: exact fit", b => Refused(b, BassRefusal.ExactFit) ? 1 : 0),
            ("rejected, 1 − α", b => b.Term is not null && b.Score.Admitted ? 1 - b.Score.Acceptance : 0),
            ("accepted, α", b => b.Term is not null ? b.Score.Acceptance : 0),
        ];
        Birth[] pooled = [.. births.SelectMany(b => b)];
        _out.WriteLine("");
        _out.WriteLine($"Births lost by cause (share of draws) | {string.Join(" | ", chains.Select((c, i) => $"chain {i} (M {c.Current.Count})"))} | pooled");
        _out.WriteLine(string.Concat(Enumerable.Repeat("--- | ", chains.Length + 1)) + "---");
        foreach ((string cause, Func<Birth, double> weight) in causes)
            _out.WriteLine($"{cause} | {string.Join(" | ", births.Select(b => Share(b.Sum(weight), b.Length)))} | {Share(pooled.Sum(weight), pooled.Length)}");
        foreach (Birth[] b in births.Append(pooled))
            Assert.True(Math.Abs(causes.Sum(cause => b.Sum(cause.Weight)) - b.Length) <= 1e-9 * b.Length, "The causes do not add up.");
        _out.WriteLine($"α per attempt (non-empty) | {string.Join(" | ", births.Append(pooled).Select(b => Share(b.Where(x => x.Term is not null).Sum(x => x.Score.Acceptance), b.Count(x => x.Term is not null))))}");
    }

    private static bool Refused(Birth b, BassRefusal reason) => b.Term is not null && b.Score.Refusal == reason;

    // Acceptance and Δℓ by tag, pooled over the chains.
    private void ByTag(Birth[][] births)
    {
        Birth[] made = [.. births.SelectMany(b => b).Where(b => b.Term is not null)];
        _out.WriteLine("");
        _out.WriteLine("tag | births | share | refused | mean α | mean α, admitted | Δℓ median | Δℓ 95% | Δℓ max");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | ---");
        foreach (IGrouping<string, Birth> g in made.GroupBy(b => b.Tag).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            Birth[] admitted = [.. g.Where(b => b.Score.Admitted)];
            double[] dl = [.. admitted.Select(b => b.Score.LogMarginalChange)];
            _out.WriteLine($"{g.Key} | {g.Count()} | {Share(g.Count(), made.Length)} | {Share(g.Count() - admitted.Length, g.Count())} | " +
                           $"{G(g.Average(b => b.Score.Acceptance))} | {(admitted.Length > 0 ? G(admitted.Average(b => b.Score.Acceptance)) : "—")} | " +
                           $"{F(BassMixing.Quantile(dl, 0.5))} | {F(BassMixing.Quantile(dl, 0.95))} | {(dl.Length > 0 ? F(dl.Max()) : "—")}");
        }
    }

    // Δℓ's quantiles over the admitted births: active and inert apart where the truth says, and all; with the threshold the
    // other parts set (the Δℓ at which the ratio is 0).
    private void Quantiles(BassMixing.Fixture fixture, Birth[][] births)
    {
        Birth[] admitted = [.. births.SelectMany(b => b).Where(b => b.Term is not null && b.Score.Admitted)];
        var groups = new List<(string, Birth[])> { ("all admitted", admitted) };
        if (fixture.Inert is not null)
        {
            groups.Add(("active", [.. admitted.Where(b => b.Inert == false)]));
            groups.Add(("inert", [.. admitted.Where(b => b.Inert == true)]));
        }
        _out.WriteLine("");
        _out.WriteLine("Δℓ over | births | 5% | 25% | median | 75% | 95% | 99% | max | threshold median | threshold 5%–95%");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        foreach ((string label, Birth[] set) in groups)
        {
            double[] dl = [.. set.Select(b => b.Score.LogMarginalChange)];
            double[] th = [.. set.Select(b => b.Score.Threshold)];
            _out.WriteLine($"{label} | {set.Length} | {string.Join(" | ", new[] { 0.05, 0.25, 0.5, 0.75, 0.95, 0.99 }.Select(q => F(BassMixing.Quantile(dl, q))))} | " +
                           $"{(dl.Length > 0 ? F(dl.Max()) : "—")} | {F(BassMixing.Quantile(th, 0.5))} | {F(BassMixing.Quantile(th, 0.05))}–{F(BassMixing.Quantile(th, 0.95))}");
        }
    }

    // Per chain: the largest Δℓ among its births, against the threshold the other parts set for that birth.
    private void Best(BassMixing.Fixture fixture, BassModel model, BassChain[] chains, Birth[][] births)
    {
        _out.WriteLine("");
        _out.WriteLine("chain | M | admitted | largest Δℓ | its term | its threshold | its ratio | largest ratio | births with α > 0.5 | with α > 0.05");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < chains.Length; c++)
        {
            Birth[] admitted = [.. births[c].Where(b => b.Term is not null && b.Score.Admitted)];
            Birth best = admitted.MaxBy(b => b.Score.LogMarginalChange)!;
            _out.WriteLine($"{c} | {chains[c].Current.Count} | {admitted.Length} | {F(best.Score.LogMarginalChange)} | {fixture.Describe(model, best.Term!)} ({best.Tag}) | " +
                           $"{F(best.Score.Threshold)} | {F(best.Score.Ratio)} | {F(admitted.Max(b => b.Score.Ratio))} | " +
                           $"{admitted.Count(b => b.Score.Acceptance > 0.5)} | {admitted.Count(b => b.Score.Acceptance > 0.05)}");
        }
    }

    // Births of K = 1 in each part, by variable: their share and acceptance, and what weighting the variable by use
    // (lever 3a: z_u ∝ w₂ + appearances of u, from each chain's counts) would make of the part's mean acceptance, the
    // per-variable acceptance held as measured.
    private void MainEffects(BassMixing.Fixture fixture, BassModel model, BassChain[] chains, Birth[][] births)
    {
        double w2 = model.Proposal.VariableWeight;
        foreach (BassPart part in Enum.GetValues<BassPart>())
        {
            if (!model.Prior.Has(part)) continue;
            int p = model.Prior.Variables(part);
            string[] names = part switch { BassPart.Continuous => fixture.Continuous, BassPart.Categorical => fixture.Categorical, _ => fixture.Outputs };
            _out.WriteLine("");
            _out.WriteLine($"{part} K = 1: variable | births | refused | mean α | {string.Join(" | ", chains.Select((_, c) => $"3a weight, chain {c}"))}");
            _out.WriteLine("--- | --- | --- | ---" + string.Concat(Enumerable.Repeat(" | ---", chains.Length)));
            var alpha = new double[chains.Length][];
            var weights = new double[chains.Length][];
            for (int c = 0; c < chains.Length; c++)
            {
                BassTermCounts counts = BassMixing.Counts(model, chains[c].Current);
                weights[c] = [.. Enumerable.Range(0, p).Select(u => w2 + counts.Appearances(part, u))];
                double total = weights[c].Sum();
                for (int u = 0; u < p; u++) weights[c][u] /= total;
                alpha[c] = new double[p];
            }
            for (int u = 0; u < p; u++)
            {
                Birth[] set = [.. births.SelectMany(b => b).Where(b => b.Term is { } t && Single(t, part) == u)];
                for (int c = 0; c < chains.Length; c++)
                {
                    Birth[] own = [.. births[c].Where(b => b.Term is { } t && Single(t, part) == u)];
                    alpha[c][u] = own.Length > 0 ? own.Average(b => b.Score.Acceptance) : 0.0;
                }
                _out.WriteLine($"{names[u]} | {set.Length} | {Share(set.Count(b => !b.Score.Admitted), set.Length)} | {G(set.Length > 0 ? set.Average(b => b.Score.Acceptance) : double.NaN)} | " +
                               string.Join(" | ", weights.Select(w => $"{w[u]:F3}")));
            }
            string uniform = string.Join(", ", alpha.Select(a => G(a.Average())));
            string weighted = string.Join(", ", alpha.Select((a, c) => G(a.Select((x, u) => x * weights[c][u]).Sum())));
            _out.WriteLine($"{part} K = 1, mean α by chain: uniform variable {uniform}; use-weighted (3a) {weighted}");
            if (fixture.Inert is not null && part == BassPart.Continuous)
                _out.WriteLine($"{part} K = 1, the active share by chain: uniform 0.500; use-weighted (3a) " +
                               string.Join(", ", weights.Select(w => w.Take(5).Sum().ToString("F3"))));
        }
    }

    // The variable of the term's one factor in `part` when it has exactly one there; −1 otherwise.
    private static int Single(BassTerm term, BassPart part) => part switch
    {
        BassPart.Continuous => term.Continuous.Count == 1 ? term.Continuous[0].Variable : -1,
        BassPart.Categorical => term.Categorical.Count == 1 ? term.Categorical[0].Variable : -1,
        _ => term.Functional.Count == 1 ? term.Functional[0].Variable : -1,
    };

    // The pure continuous main effects (one factor, continuous), active ones where the truth says, by their knot's quantile
    // within the variable's distinct values, j/(n_v − 1), in fifths, and by sign.
    private void Knots(BassMixing.Fixture fixture, BassModel model, Birth[][] births)
    {
        Birth[] mains = [.. births.SelectMany(b => b).Where(b => b.Term is { } t && t.Continuous.Count == 1 && t.Categorical.Count == 0
                                                                   && t.Functional.Count == 0 && b.Inert != true)];
        _out.WriteLine("");
        _out.WriteLine($"{(fixture.Inert is null ? "" : "active ")}continuous main effects: knot quantile | sign | births | refused | mean α | Δℓ median | Δℓ 95% | Δℓ max");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
        foreach (int sign in new[] { -1, 1 })
            for (int fifth = 0; fifth < 5; fifth++)
            {
                Birth[] set = [.. mains.Where(b => b.Term!.Continuous[0].Sign == sign && Fifth(model, b.Term.Continuous[0]) == fifth)];
                double[] dl = [.. set.Where(b => b.Score.Admitted).Select(b => b.Score.LogMarginalChange)];
                _out.WriteLine($"[{fifth / 5.0:F1}, {(fifth + 1) / 5.0:F1}) | {(sign > 0 ? "+" : "−")} | {set.Length} | {Share(set.Count(b => !b.Score.Admitted), set.Length)} | " +
                               $"{G(set.Length > 0 ? set.Average(b => b.Score.Acceptance) : double.NaN)} | {F(BassMixing.Quantile(dl, 0.5))} | {F(BassMixing.Quantile(dl, 0.95))} | {(dl.Length > 0 ? F(dl.Max()) : "—")}");
            }
    }

    private static int Fifth(BassModel model, HingeFactor h)
    {
        int n = model.Inputs.ContinuousKnots(h.Variable).Length;
        return Math.Min(4, (int)(5.0 * h.Knot / Math.Max(1, n - 1)));
    }

    // Each chain's deaths: the state's terms, the Δℓ of removing each and its full ratio.
    private void DeathTable(BassMixing.Fixture fixture, BassModel model, Death[][] deaths)
    {
        _out.WriteLine("");
        _out.WriteLine("chain | term | tag | Δℓ | ratio | α");
        _out.WriteLine("--- | --- | --- | --- | --- | ---");
        for (int c = 0; c < deaths.Length; c++)
            foreach (Death d in deaths[c].OrderByDescending(d => d.Score.Ratio))
                _out.WriteLine($"{c} | {fixture.Describe(model, d.Term)} | {d.Tag} | {F(d.Score.LogMarginalChange)} | {F(d.Score.Ratio)} | {G(d.Score.Acceptance)}");
        Death[] all = [.. deaths.SelectMany(d => d)];
        _out.WriteLine($"deaths: {all.Length}; Δℓ median {F(BassMixing.Quantile(all.Select(d => d.Score.LogMarginalChange), 0.5))}, " +
                       $"largest {F(all.Max(d => d.Score.LogMarginalChange))}; ratio largest {F(all.Max(d => d.Score.Ratio))}; mean α over terms {G(all.Average(d => d.Score.Acceptance))}");
    }

    private static string Share(double part, double whole) => whole > 0 ? $"{100 * part / whole:F2}%" : "—";
    private static string F(double v) => double.IsNaN(v) ? "—" : v.ToString("F2");
    private static string G(double v) => double.IsNaN(v) ? "—" : v.ToString("G3");
}
