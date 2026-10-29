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
/// F-M2's M3 (the local-move measurements brief): changes by knot displacement. From each chain's state at the end of
/// S10's burn-in, on each of its fixtures, 10,000 changes are proposed through the chain's own change move and scored as
/// M2 scores births (<see cref="BassMixing.Score"/>), none committed. A hinge change is grouped by its part, by whether
/// its sign flipped and by its knot's displacement |j′ − j| in the variable's distinct-value index; a categorical change
/// by its subset's size change. A redraw equal to the current factor is the same state, which the chain accepts as a
/// null move, and is read apart. It reads, it does not judge; it asserts its own soundness: nothing committed, the parts
/// summing to the chain's ratio, one factor changed per proposal, and the counts matching the move's tallies.
/// </summary>
public sealed class BassMixingChangeTests
{
    private const int Proposals = 10_000;

    private static readonly (string Label, int Low, int High)[] Buckets =
        [("0", 0, 0), ("1", 1, 1), ("2–3", 2, 3), ("4–7", 4, 7), ("8–15", 8, 15), ("16–31", 16, 31), ("≥ 32", 32, int.MaxValue)];

    private readonly ITestOutputHelper _out;
    public BassMixingChangeTests(ITestOutputHelper output) => _out = output;

    // One proposed change: its part ("x" continuous, "z" categorical, "t" functional, or "same" for a redraw equal to the
    // current factor), a hinge's displacement and sign flip, a subset's size change, and its score.
    private sealed record Change(string Part, int Displacement, bool Flip, int SizeChange, BassMixing.Scored Score);

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void Changes_ByKnotDisplacement(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassMixing.Fixture fixture = BassMixing.Get(name);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassChain[] chains = BassMixing.Burned(model);
        var changes = new Change[chains.Length][];
        int[] seeds = SeedTree.Derive(109, chains.Length);
        Parallel.For(0, chains.Length, c =>
        {
            var untouched = new BassMixing.Untouched(model, chains[c]);
            changes[c] = Changes(model, chains[c], new Xoshiro256PlusPlus(seeds[c]));
            untouched.Check();
        });
        _out.WriteLine($"## M3, {name}: {Proposals:N0} changes from each chain's state after {BassMixing.BurnIn:N0} sweeps ({clock.Elapsed.TotalSeconds:F1} s)");
        Change[] pooled = [.. changes.SelectMany(c => c)];
        Hinges(pooled);
        Subsets(pooled);
        PerChain(chains, changes);
    }

    private static Change[] Changes(BassModel model, BassChain chain, Xoshiro256PlusPlus rng)
    {
        BassState state = chain.Current;
        int m = state.Count;
        Assert.True(m > 0, "The state holds no term to change.");
        IReadOnlyList<IRjMove<BassState>> moves = chain.Structure.Moves;
        IRjMove<BassState> change = BassMixing.Move(moves, "change");
        chain.Structure.ResetTallies();
        var changes = new Change[Proposals];
        for (int k = 0; k < Proposals; k++)
        {
            Proposal<BassState> p = change.Propose(state, rng) ?? throw new InvalidOperationException("A change from M > 0 is always available.");
            if (ReferenceEquals(p.Candidate, state))
            {
                Assert.Equal(0.0, p.LogProposalRatio);
                changes[k] = new Change("same", 0, false, 0, new BassMixing.Scored(null, false, 0, 0, 0, 0, 0, 0));
                continue;
            }
            // The design moves the changed term last: the candidate is the state less term i, then the new term.
            BassTerm added = p.Candidate.Terms[^1];
            int i = 0;
            while (i < m - 1 && ReferenceEquals(p.Candidate.Terms[i], state.Terms[i])) i++;
            BassTerm old = state.Terms[i];
            BassTerm[] expected = [.. state.Terms.Take(i), .. state.Terms.Skip(i + 1), added];
            Assert.Equal(expected, p.Candidate.Terms);
            (string part, int factor) = Changed(old, added);
            double within = 0.0;
            int displacement = 0, sizeChange = 0;
            bool flip = false;
            if (part == "z")
            {
                CategoricalFactor was = old.Categorical[factor], now = added.Categorical[factor];
                double[] log = model.Prior[BassPart.Categorical]!.LogSubset![was.Variable];
                within = log[was.Levels.Count] - log[now.Levels.Count];
                sizeChange = now.Levels.Count - was.Levels.Count;
            }
            else
            {
                HingeFactor was = part == "x" ? old.Continuous[factor] : old.Functional[factor];
                HingeFactor now = part == "x" ? added.Continuous[factor] : added.Functional[factor];
                Assert.Equal(was.Variable, now.Variable);
                displacement = Math.Abs(now.Knot - was.Knot);
                flip = now.Sign != was.Sign;
            }
            Assert.Equal(BassMixing.Bits(within), BassMixing.Bits(p.LogProposalRatio));
            BassMixing.Scored score = BassMixing.Score(model, moves, change, state, p, model.Prior.LogPrior(added) - model.Prior.LogPrior(old),
                                                       within, repeated: state.Terms.Contains(added));
            changes[k] = new Change(part, displacement, flip, sizeChange, score);
        }
        foreach (BassRefusal reason in Enum.GetValues<BassRefusal>())
            Assert.Equal(changes.LongCount(x => x.Score.Refusal == reason), chain.Structure.Refusals.Count(reason));
        return changes;
    }

    // The one factor a change redrew: its part and its index there.
    private static (string Part, int Factor) Changed(BassTerm old, BassTerm added)
    {
        var differ = new List<(string, int)>();
        Assert.Equal(old.Continuous.Count, added.Continuous.Count);
        Assert.Equal(old.Categorical.Count, added.Categorical.Count);
        Assert.Equal(old.Functional.Count, added.Functional.Count);
        for (int f = 0; f < old.Continuous.Count; f++)
            if (old.Continuous[f] != added.Continuous[f]) differ.Add(("x", f));
        for (int f = 0; f < old.Categorical.Count; f++)
            if (!old.Categorical[f].Equals(added.Categorical[f])) differ.Add(("z", f));
        for (int f = 0; f < old.Functional.Count; f++)
            if (old.Functional[f] != added.Functional[f]) differ.Add(("t", f));
        Assert.Single(differ);
        return differ[0];
    }

    // Hinge changes by part, sign flip and displacement; the redraws equal to the current factor apart.
    private void Hinges(Change[] pooled)
    {
        _out.WriteLine("");
        _out.WriteLine("part | sign | displacement | changes | share | refused | mean α | its se | mean α, admitted | Δℓ median");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        foreach (string part in new[] { "x", "t" }.Where(part => pooled.Any(x => x.Part == part)))
            foreach (bool flip in new[] { false, true })
                foreach ((string label, int low, int high) in Buckets)
                {
                    Change[] set = [.. pooled.Where(x => x.Part == part && x.Flip == flip && x.Displacement >= low && x.Displacement <= high)];
                    Row($"{part} | {(flip ? "flipped" : "kept")} | {label}", set, pooled.Length);
                }
        Row("same factor redrawn | — | —", [.. pooled.Where(x => x.Part == "same")], pooled.Length);
        Row("all changes | — | —", pooled, pooled.Length);
    }

    // Categorical changes by the subset's size change.
    private void Subsets(Change[] pooled)
    {
        Change[] subsets = [.. pooled.Where(x => x.Part == "z")];
        if (subsets.Length == 0) return;
        _out.WriteLine("");
        _out.WriteLine("subset size change | changes | share | refused | mean α | its se | mean α, admitted | Δℓ median");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
        foreach (IGrouping<int, Change> g in subsets.GroupBy(x => x.SizeChange).OrderBy(g => g.Key))
            Row($"{g.Key:+0;−0;0}", [.. g], pooled.Length);
    }

    // Per chain: the state's M, the mean α of all its changes, and of the kept-sign hinge changes near and far.
    private void PerChain(BassChain[] chains, Change[][] changes)
    {
        _out.WriteLine("");
        _out.WriteLine("chain | M | mean α, all | same factor | kept sign, displacement 1–3 (n) | 4–15 (n) | ≥ 32 (n) | flipped (n)");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < chains.Length; c++)
        {
            Change[] hinges = [.. changes[c].Where(x => x.Part is "x" or "t")];
            string Mean(Func<Change, bool> keep)
            {
                Change[] set = [.. hinges.Where(keep)];
                return set.Length > 0 ? $"{set.Average(x => x.Score.Acceptance):G3} ({set.Length})" : "— (0)";
            }
            _out.WriteLine($"{c} | {chains[c].Current.Count} | {changes[c].Average(x => x.Score.Acceptance):G3} | " +
                           $"{(double)changes[c].Count(x => x.Part == "same") / changes[c].Length:P2} | {Mean(x => !x.Flip && x.Displacement is >= 1 and <= 3)} | " +
                           $"{Mean(x => !x.Flip && x.Displacement is >= 4 and <= 15)} | {Mean(x => !x.Flip && x.Displacement >= 32)} | {Mean(x => x.Flip)}");
        }
    }

    private void Row(string label, Change[] set, int total)
    {
        Change[] admitted = [.. set.Where(x => x.Score.Admitted)];
        double mean = set.Length > 0 ? set.Average(x => x.Score.Acceptance) : double.NaN;
        double se = set.Length > 1 ? Math.Sqrt(set.Sum(x => (x.Score.Acceptance - mean) * (x.Score.Acceptance - mean)) / (set.Length - 1) / set.Length) : double.NaN;
        _out.WriteLine($"{label} | {set.Length} | {100.0 * set.Length / total:F2}% | {(set.Length > 0 ? $"{100.0 * (set.Length - admitted.Length) / set.Length:F2}%" : "—")} | " +
                       $"{G(mean)} | {G(se)} | {(admitted.Length > 0 ? G(admitted.Average(x => x.Score.Acceptance)) : "—")} | " +
                       $"{(admitted.Length > 0 ? BassMixing.Quantile(admitted.Select(x => x.Score.LogMarginalChange), 0.5).ToString("F2") : "—")}");
    }

    private static string G(double v) => double.IsNaN(v) ? "—" : v.ToString("G3");
}
