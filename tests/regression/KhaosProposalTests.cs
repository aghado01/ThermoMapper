using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's proposal on F-X4's contract (the KHAOS brief's slice 1): the inclusion's exact law against brute force over every
/// inclusion vector, validity included; the birth's probabilities summing to one over the term space and its draws following
/// them; the death's reverse; the two mutations' draws against their statement and their ratios reciprocal; the swap never
/// onto a held variable; and <c>make_weights</c>' properties.
/// </summary>
public sealed class KhaosProposalTests
{
    internal static (BassStructurePrior Prior, LegendreTermProposal Proposal) Build(int p, int d, int q, LegendreBasis? settings = null)
    {
        var prior = new BassStructurePrior(new LegendreTermSpace(p, d, q), new BassCountPrior(4, 0.5));
        return (prior, new LegendreTermProposal(prior, (settings ?? new LegendreBasis()) with { MaxDegree = d }));
    }

    private static BassTermCounts Counts(BassStructurePrior prior, params BassTerm[] terms)
    {
        var counts = new BassTermCounts(prior);
        foreach (BassTerm t in terms) counts.Add(t);
        return counts;
    }

    private static BassTerm T(params (int Variable, int Degree)[] factors) => new(legendre: factors.Select(f => new LegendreFactor(f.Variable, f.Degree)));

    // Σ over valid χ′ of Π η^χ′(1 − η)^(1 − χ′), by enumeration: the validity with no Poisson-binomial.
    private static double BruteValidity(double[] eta, int top)
    {
        double total = 0.0;
        for (int mask = 1; mask < 1 << eta.Length; mask++)
        {
            int count = 0;
            double product = 1.0;
            for (int j = 0; j < eta.Length; j++)
            {
                bool on = ((mask >> j) & 1) != 0;
                count += on ? 1 : 0;
                product *= on ? eta[j] : 1.0 - eta[j];
            }
            if (count <= top) total += product;
        }
        return total;
    }

    [Fact]
    public void Inclusion_IsTheExactMixtureOverQ0_AgainstBruteForce()
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(5, 6, 3);
        int[] appearances = [0, 2, 1, 0, 3];
        double[] j = [.. Enumerable.Range(1, 3).Select(q0 => 1.0 / q0)];
        double h = j.Sum();
        var etas = new double[3][];
        for (int q0 = 1; q0 <= 3; q0++)
        {
            etas[q0 - 1] = new double[5];
            proposal.InclusionProbabilities(appearances, q0, etas[q0 - 1]);
        }
        double total = 0.0, worst = 0.0;
        for (int mask = 1; mask < 32; mask++)
        {
            bool[] chi = [.. Enumerable.Range(0, 5).Select(v => ((mask >> v) & 1) != 0)];
            if (chi.Count(c => c) > 3) continue;
            double exact = 0.0;
            for (int q0 = 1; q0 <= 3; q0++)
            {
                double[] eta = etas[q0 - 1];
                double product = 1.0;
                for (int v = 0; v < 5; v++) product *= chi[v] ? eta[v] : 1.0 - eta[v];
                exact += j[q0 - 1] / h * product / BruteValidity(eta, 3);
            }
            double engine = Math.Exp(proposal.LogInclusion(appearances, chi));
            worst = Math.Max(worst, Math.Abs(engine / exact - 1.0));
            total += engine;
        }
        Assert.True(worst < 1e-12, $"worst relative error {worst:E2}");
        Assert.Equal(1.0, total, 1e-13);
    }

    [Theory]
    [InlineData(4, 3, 2)]
    [InlineData(5, 5, 3)]
    [InlineData(3, 4, 3)]
    public void Birth_SumsToOneOverTheSpace(int p, int d, int q)
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(p, d, q, new LegendreBasis { DegreePenalty = 0.7, OrderPenalty = 1.3 });
        foreach (BassTermCounts counts in new[] { Counts(prior), Counts(prior, T((0, 1)), T((0, 2), (p - 1, 1))) })
        {
            double total = prior.Legendre!.Terms().Sum(t => Math.Exp(proposal.LogBirth(t, counts)));
            Assert.Equal(1.0, total, 1e-12);
        }
    }

    [Fact]
    public void BirthDraws_FollowTheStatedLaw()
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(4, 3, 2, new LegendreBasis { DegreePenalty = 0.5 });
        BassTermCounts counts = Counts(prior, T((1, 2)), T((1, 1), (3, 1)));
        List<BassTerm> terms = [.. prior.Legendre!.Terms()];
        double[] expected = [.. terms.Select(t => Math.Exp(proposal.LogBirth(t, counts)))];
        var tally = terms.ToDictionary(t => t, _ => 0);
        var rng = new Xoshiro256PlusPlus(20261107);
        const int draws = 400_000;
        for (int k = 0; k < draws; k++) tally[proposal.DrawBirth(rng, counts)!]++;
        double chi = 0.0;
        for (int i = 0; i < terms.Count; i++)
        {
            double e = expected[i] * draws;
            chi += (tally[terms[i]] - e) * (tally[terms[i]] - e) / e;
        }
        int df = terms.Count - 1;
        Assert.True(chi < df + 6 * Math.Sqrt(2.0 * df), $"χ² {chi:F1} on {df} degrees of freedom");
    }

    [Fact]
    public void Death_IsTheBirthUnderTheReducedCounts()
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(5, 6, 3);
        BassTerm dying = T((0, 2), (4, 1));
        BassTerm[] others = [T((0, 1)), T((2, 3), (4, 2))];
        BassTermCounts with = Counts(prior, [.. others, dying]), without = Counts(prior, others);
        Assert.Equal(proposal.LogBirth(dying, without) + Math.Log(3.0), proposal.LogDeathRatio(dying, with), 1e-13);
        Assert.Equal(-Math.Log(3.0) - proposal.LogBirth(dying, without), proposal.LogBirthRatio(dying, without), 1e-13);
        Assert.Throws<ArgumentException>(() => proposal.LogDeathRatio(T((1, 1)), without));   // no term uses variable 1
    }

    [Fact]
    public void DegreeMutation_DrawsItsStatement_AndItsRatiosAreReciprocal()
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(4, 6, 3, new LegendreBasis { DegreePenalty = 0.8 });
        BassTerm term = T((0, 1), (2, 2));
        var tally = new Dictionary<BassTerm, int>();
        var rng = new Xoshiro256PlusPlus(3);
        const int draws = 200_000;
        for (int k = 0; k < draws; k++)
        {
            BassTermMove move = proposal.ProposeDegree(rng, term);
            Assert.Equal(BassMoveKind.DegreeChange, move.Kind);
            BassTerm next = move.Term!;
            Assert.Equal([0, 2], next.Legendre.Select(f => f.Variable));
            tally[next] = tally.GetValueOrDefault(next) + 1;
        }
        // The statement: d′ ∝ (d′ − q + 1)^(−s_d) on {2, …, 6}, then a uniform composition.
        double[] weights = [.. Enumerable.Range(2, 5).Select(dd => Math.Pow(dd - 1, -0.8))];
        double sum = weights.Sum(), chi = 0.0;
        int cells = 0;
        for (int dd = 2; dd <= 6; dd++)
            for (int a = 1; a < dd; a++)
            {
                double e = draws * weights[dd - 2] / sum / (dd - 1);
                int seen = tally.GetValueOrDefault(T((0, a), (2, dd - a)));
                chi += (seen - e) * (seen - e) / e;
                cells++;
            }
        Assert.Equal(cells, tally.Count);
        Assert.True(chi < cells - 1 + 6 * Math.Sqrt(2.0 * (cells - 1)), $"χ² {chi:F1} on {cells - 1} degrees of freedom");
        for (int from = 2; from <= 6; from++)
            for (int to = 2; to <= 6; to++)
                Assert.Equal(-proposal.LogDegreeRatio(2, to, from), proposal.LogDegreeRatio(2, from, to), 1e-14);
        Assert.Equal(Math.Log(weights[1] / sum / 2) - Math.Log(weights[3] / sum / 4), proposal.LogDegreeRatio(2, 3, 5), 1e-13);
    }

    [Fact]
    public void Swap_DrawsAmongTheVariablesLeftOut_AndItsRatiosAreReciprocal()
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(6, 6, 3, new LegendreBasis { SwapPseudoCount = 0.7 });
        BassTerm term = T((1, 2), (4, 1));
        BassTermCounts counts = Counts(prior, term, T((0, 1)), T((0, 2), (5, 1)), T((4, 3)));
        int[] a = counts.AppearanceCounts(BassPart.Legendre);
        var tally = new Dictionary<(int From, int To), int>();
        var rng = new Xoshiro256PlusPlus(8);
        const int draws = 200_000;
        for (int k = 0; k < draws; k++)
        {
            BassTermMove move = proposal.ProposeSwap(rng, term, a);
            Assert.Equal(BassMoveKind.VariableChange, move.Kind);
            int[] before = [.. term.Legendre.Select(f => f.Variable)], after = [.. move.Term!.Legendre.Select(f => f.Variable)];
            int from = before.Except(after).Single(), to = after.Except(before).Single();
            Assert.Equal(term.Legendre.Single(f => f.Variable == from).Degree, move.Term.Legendre.Single(f => f.Variable == to).Degree);
            tally[(from, to)] = tally.GetValueOrDefault((from, to)) + 1;
            Assert.Equal(proposal.LogSwapRatio(term, from, to, a), move.LogProposalRatio, 1e-15);
        }
        int[] left = [0, 2, 3, 5];
        double total = left.Sum(u => 0.7 + a[u]), chi = 0.0;
        foreach (int from in new[] { 1, 4 })
            foreach (int to in left)
            {
                double e = draws * 0.5 * (0.7 + a[to]) / total;
                int seen = tally.GetValueOrDefault((from, to));
                chi += (seen - e) * (seen - e) / e;
            }
        Assert.Equal(8, tally.Count);
        Assert.True(chi < 7 + 6 * Math.Sqrt(14.0), $"χ² {chi:F1} on 7 degrees of freedom");

        // Reciprocal: the swapped term's swap back, under the counts after the swap.
        foreach ((int from, int to) in tally.Keys)
        {
            BassTerm swapped = T([.. term.Legendre.Select(f => f.Variable == from ? (to, f.Degree) : (f.Variable, f.Degree))]);
            int[] after = (int[])a.Clone();
            after[from]--;
            after[to]++;
            Assert.Equal(-proposal.LogSwapRatio(swapped, to, from, after), proposal.LogSwapRatio(term, from, to, a), 1e-13);
        }
        Assert.Throws<ArgumentException>(() => proposal.LogSwapRatio(term, 1, 4, a));
    }

    [Fact]
    public void Change_MixesTheMutations_AndATermOverEveryVariableHasNoSwap()
    {
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(4, 8, 4, new LegendreBasis { DegreeMutation = 0.3 });
        BassTerm full = T((0, 1), (1, 1), (2, 1), (3, 2));
        BassTermCounts counts = Counts(prior, full, T((2, 2)));
        var rng = new Xoshiro256PlusPlus(12);
        int degree = 0, nulls = 0;
        const int draws = 50_000;
        for (int k = 0; k < draws; k++)
        {
            BassTermMove move = proposal.ProposeChange(rng, full, counts);
            if (move.Kind == BassMoveKind.DegreeChange) degree++;
            else
            {
                Assert.True(move.IsNull);
                nulls++;
            }
        }
        Assert.Equal(0.3, (double)degree / draws, 0.01);
        Assert.Equal(draws, degree + nulls);
        (_, LegendreTermProposal small) = Build(3, 4, 2, new LegendreBasis { DegreeMutation = 0.1 });
        Assert.Equal(1.0, small.DegreeMutation);   // KHAOS's rule: p ≤ 3 mutates the degree only
        Assert.Equal(0.3, proposal.DegreeMutation);
    }

    [Fact]
    public void InclusionProbabilities_SumNearQ0_AndStayInTheUnitInterval()
    {
        (_, LegendreTermProposal proposal) = Build(8, 10, 5);
        var rng = new Xoshiro256PlusPlus(21);
        var eta = new double[8];
        for (int c = 0; c < 200; c++)
        {
            int[] a = [.. Enumerable.Range(0, 8).Select(_ => (int)(rng.NextDouble() * 12))];
            for (int q0 = 1; q0 <= 5; q0++)
            {
                proposal.InclusionProbabilities(a, q0, eta);
                Assert.All(eta, e => Assert.True(e > 0.0 && e < 1.0));
                Assert.True(Math.Abs(eta.Sum() - q0) < 0.05 * q0, $"Σ η = {eta.Sum()} at q₀ = {q0}");
            }
        }
        (_, LegendreTermProposal one) = Build(1, 3, 1);
        var single = new double[1];
        one.InclusionProbabilities([5], 1, single);
        Assert.Equal(1.0, single[0]);
        (_, LegendreTermProposal everyVariable) = Build(3, 5, 3);
        var all = new double[3];
        everyVariable.InclusionProbabilities([1, 0, 4], 3, all);
        Assert.All(all, e => Assert.Equal(1.0, e, 1e-15));
    }

    [Fact]
    public void Proposal_RefusesAHingePrior_AndForeignCounts()
    {
        BassData data = KhaosBasisTests.UnitData(12, 2, 2);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var hinges = new BassStructurePrior(inputs, 2, 0, 0, CategoricalSubsetMeasure.SizeUniform, new BassCountPrior());
        Assert.Throws<ArgumentException>(() => new LegendreTermProposal(hinges, new LegendreBasis()));
        (BassStructurePrior prior, LegendreTermProposal proposal) = Build(3, 4, 2);
        (BassStructurePrior other, _) = Build(3, 4, 2);
        Assert.Throws<ArgumentException>(() => proposal.LogBirth(T((0, 1)), new BassTermCounts(other)));
        Assert.Throws<ArgumentException>(() => new LegendreTermProposal(prior, new LegendreBasis { MaxDegree = 5 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LegendreTermProposal(prior, new LegendreBasis { MaxDegree = 4, InclusionPseudoCount = 0 }));
    }
}
