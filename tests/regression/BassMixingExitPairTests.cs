using System;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's P1 (the mode-exit pricing brief): pair deaths. From each chain's fixed state (<see cref="BassMixing.FixedStates"/>)
/// on each of S10's fixtures, every unordered pair of terms is removed at once and scored at fixed τ: its Δℓ, and the pair
/// death's full log ratio, the target's (log marginal and log structure prior) plus ln q_pair(t₁, t₂ | S′) + ln C(M, 2),
/// with S′ the reduced structure and q_pair the two D01 births from it in either order (<see cref="BassMixing.Removals"/>);
/// the move's selection is taken as 1. Beside each pair, its terms' single deaths, scored alike. Where M ≤ 12 and a state's
/// best pair costs under 20 nats, every triple's too. It reads, it does not judge; it asserts its own soundness: nothing
/// committed (<see cref="BassMixing.Untouched"/>), the single deaths matching the chain's design and death ratio, and the
/// removals counting C(M, k).
/// </summary>
public sealed class BassMixingExitPairTests
{
    private const double Near = -10.0, Far = -20.0;
    private readonly ITestOutputHelper _out;
    public BassMixingExitPairTests(ITestOutputHelper output) => _out = output;

    // One state's removals: every single, every pair, and every triple where read.
    private sealed record Read(BassState State, BassMixing.Removal[] Singles, BassMixing.Removal[] Pairs, BassMixing.Removal[]? Triples);

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void PairDeaths_FromFixedStates(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassMixing.Fixture fixture = BassMixing.Get(name);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassMixing.FixedState[] states = BassMixing.FixedStates(model);
        double setup = clock.Elapsed.TotalSeconds;
        clock.Restart();
        var reads = new Read[states.Length];
        Parallel.For(0, states.Length, c =>
        {
            BassChain chain = states[c].Chain;
            var untouched = new BassMixing.Untouched(model, chain);
            reads[c] = Removals(model, chain);
            untouched.Check();
        });
        _out.WriteLine($"## P1, {name}: pair deaths from each chain's fixed state (states and phases {setup:F1} s, removals {clock.Elapsed.TotalSeconds:F1} s)");
        _out.WriteLine("");
        foreach (string line in BassMixing.FixedTable(model, states)) _out.WriteLine(line);
        Summary(states, reads);
        for (int c = 0; c < reads.Length; c++) Best(fixture, model, states[c], reads[c]);
    }

    private static Read Removals(BassModel model, BassChain chain)
    {
        BassState state = chain.Current;
        int m = state.Count;
        BassMixing.Removal[] singles = BassMixing.Removals(model, state, 1), pairs = BassMixing.Removals(model, state, 2);
        Assert.Equal(m, singles.Length);
        Assert.Equal(m * (m - 1) / 2, pairs.Length);
        Assert.All(singles, s => Assert.True(s.Exit.Admitted, "A death never leaves F."));
        Assert.All(pairs, s => Assert.True(s.Exit.Admitted, "A death never leaves F."));

        // The single deaths against the chain's own: its design's account of each, and the death move's within-move ratio.
        BassDesign design = chain.Structure.DesignAt(state);
        BassTermCounts counts = BassMixing.Counts(model, state);
        foreach (BassMixing.Removal s in singles)
        {
            int i = s.Terms[0];
            DesignProposal death = design.ProposeDeath(i);
            BassStructureFit fit = s.Reduced!.Fit;
            Assert.True(Math.Abs(death.ResidualSquares - fit.ResidualSquares) <= 1e-9 * fit.ResidualSquares,
                        $"RSS {fit.ResidualSquares:R} against the design's {death.ResidualSquares:R}");
            Assert.True(Math.Abs(death.FittedSquares - fit.FittedSquares) <= 1e-9 * (fit.FittedSquares + fit.ResidualSquares),
                        $"Q {fit.FittedSquares:R} against the design's {death.FittedSquares:R}");
            double ratio = model.Proposal.LogDeathRatio(state.Terms[i], counts);
            Assert.True(Math.Abs(ratio - s.Exit.LogProposalRatio) <= 1e-12 * (1.0 + Math.Abs(ratio)),
                        $"The death's ratio {ratio:R} against the removal's {s.Exit.LogProposalRatio:R}");
        }
        // Each pair's reverse sums both orders: it is at least either order's births.
        foreach (BassMixing.Removal p in pairs)
        {
            BassTerm a = state.Terms[p.Terms[0]], b = state.Terms[p.Terms[1]];
            BassTerm[] rest = [.. state.Terms.Where((_, i) => i != p.Terms[0] && i != p.Terms[1])];
            BassTermCounts reduced = BassMixing.Counts(model, rest), withA = reduced.Clone(), withB = reduced.Clone();
            withA.Add(a);
            withB.Add(b);
            double ab = model.Proposal.LogBirth(a, reduced) + model.Proposal.LogBirth(b, withA);
            double ba = model.Proposal.LogBirth(b, reduced) + model.Proposal.LogBirth(a, withB);
            double pair = Math.Log(Math.Exp(ab) + Math.Exp(ba)) + BassMixing.LogChoose(state.Count, 2);
            Assert.True(Math.Abs(pair - p.Exit.LogProposalRatio) <= 1e-10 * (1.0 + Math.Abs(pair)),
                        $"q_pair {pair:R} against the removal's {p.Exit.LogProposalRatio:R}");
        }

        BassMixing.Removal[]? triples = null;
        if (m <= 12 && pairs.Length > 0 && pairs.Max(p => p.Exit.LogLikelihoodChange) > Far)
        {
            triples = BassMixing.Removals(model, state, 3);
            Assert.Equal(m * (m - 1) * (m - 2) / 6, triples.Length);
        }
        return new Read(state, singles, pairs, triples);
    }

    // Per chain: the best single and pair, the pairs near the state, the co-adapted pairs, and the triples where read.
    private void Summary(BassMixing.FixedState[] states, Read[] reads)
    {
        _out.WriteLine("");
        _out.WriteLine($"chain | M | best single Δℓ | best pair Δℓ | its ratio | its α | pairs Δℓ > {Near} | > {Far} | of C(M, 2) | " +
                       "pairs costing less than their cheaper single | largest relief | largest Δℓ_pair − Δℓ₁ − Δℓ₂ | " +
                       $"best triple Δℓ | triples > {Near} | > {Far} | of C(M, 3)");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < reads.Length; c++)
        {
            Read r = reads[c];
            BassMixing.Removal best = r.Pairs.MaxBy(p => p.Exit.LogLikelihoodChange)!;
            string triples = r.Triples is { } t
                ? $"{F(t.Max(x => x.Exit.LogLikelihoodChange))} | {t.Count(x => x.Exit.LogLikelihoodChange > Near)} | {t.Count(x => x.Exit.LogLikelihoodChange > Far)} | {t.Length}"
                : "— | — | — | —";
            _out.WriteLine($"{states[c].Label} | {r.State.Count} | {F(r.Singles.Max(s => s.Exit.LogLikelihoodChange))} | {F(best.Exit.LogLikelihoodChange)} | " +
                           $"{F(best.Exit.Ratio)} | {G(best.Exit.Acceptance)} | {r.Pairs.Count(p => p.Exit.LogLikelihoodChange > Near)} | " +
                           $"{r.Pairs.Count(p => p.Exit.LogLikelihoodChange > Far)} | {r.Pairs.Length} | {r.Pairs.Count(p => Relief(r, p) > 0.0)} | " +
                           $"{F(r.Pairs.Max(p => Relief(r, p)))} | {F(r.Pairs.Max(p => Interaction(r, p)))} | {triples}");
        }
    }

    // A state's six cheapest pairs, its three most redundant beside them (the largest Δℓ_pair − Δℓ₁ − Δℓ₂), and its three
    // cheapest triples where read.
    private void Best(BassMixing.Fixture fixture, BassModel model, BassMixing.FixedState fixedState, Read r)
    {
        string Name(int i) => fixture.Describe(model, r.State.Terms[i]);
        double Single(int i) => r.Singles[i].Exit.LogLikelihoodChange;
        _out.WriteLine("");
        _out.WriteLine($"Chain {fixedState.Label}, M {r.State.Count}: pairs | Δℓ | ratio | α | single Δℓ, each | relief (pair less the cheaper single) | Δℓ_pair − Δℓ₁ − Δℓ₂");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | ---");
        BassMixing.Removal[] cheapest = [.. r.Pairs.OrderByDescending(p => p.Exit.LogLikelihoodChange).Take(6)];
        BassMixing.Removal[] redundant = [.. r.Pairs.Except(cheapest).OrderByDescending(p => Interaction(r, p)).Take(3)];
        foreach ((BassMixing.Removal p, string mark) in cheapest.Select(p => (p, "")).Concat(redundant.Select(p => (p, " (redundant)"))))
        {
            int a = p.Terms[0], b = p.Terms[1];
            _out.WriteLine($"{Name(a)} + {Name(b)}{mark} | {F(p.Exit.LogLikelihoodChange)} | {F(p.Exit.Ratio)} | {G(p.Exit.Acceptance)} | " +
                           $"{F(Single(a))}, {F(Single(b))} | {F(Relief(r, p))} | {F(p.Exit.LogLikelihoodChange - Single(a) - Single(b))}");
        }
        if (r.Triples is not { } triples) return;
        foreach (BassMixing.Removal t in triples.OrderByDescending(t => t.Exit.LogLikelihoodChange).Take(3))
            _out.WriteLine($"triple {string.Join(" + ", t.Terms.Select(Name))} | {F(t.Exit.LogLikelihoodChange)} | {F(t.Exit.Ratio)} | {G(t.Exit.Acceptance)} | " +
                           $"{string.Join(", ", t.Terms.Select(i => F(Single(i))))} | {F(t.Exit.LogLikelihoodChange - t.Terms.Max(Single))} | —");
    }

    // How much less the pair's removal costs than its cheaper single's: positive when removing both costs less than either.
    private static double Relief(Read r, BassMixing.Removal p)
        => p.Exit.LogLikelihoodChange - Math.Max(r.Singles[p.Terms[0]].Exit.LogLikelihoodChange, r.Singles[p.Terms[1]].Exit.LogLikelihoodChange);

    // How much less the pair's removal costs than its singles' together: positive when the two carry shared structure.
    private static double Interaction(Read r, BassMixing.Removal p)
        => p.Exit.LogLikelihoodChange - r.Singles[p.Terms[0]].Exit.LogLikelihoodChange - r.Singles[p.Terms[1]].Exit.LogLikelihoodChange;

    private static string F(double v) => double.IsNaN(v) ? "—" : v.ToString("F2");
    private static string G(double v) => double.IsNaN(v) ? "—" : v.ToString("G3");
}
