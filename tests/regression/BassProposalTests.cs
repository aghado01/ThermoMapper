using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's term measure and structure proposals (stage 4, S3; D01, D02). Every term space here is small enough to
/// enumerate, and the enumeration is the test's own: the prior and the birth proposal sum to one over it, the empty
/// draw included where the fit has several parts; draws follow the probabilities stated for them; a death's ratio
/// reverses its birth's; a change's ratio is its factor's prior ratio; the count prior is the negative binomial. One
/// part, continuous and categorical, continuous and functional, and all three; both subset measures.
/// </summary>
public sealed class BassProposalTests
{
    private readonly ITestOutputHelper _out;
    public BassProposalTests(ITestOutputHelper output) => _out = output;

    private sealed record Space(string Name, BassInputs Inputs, BassStructurePrior Prior, BassTermProposal Proposal, List<BassTerm> Terms);

    // Continuous predictors with the given distinct values, categorical ones with the given levels, and an output
    // variable with the given distinct values (null: a scalar response).
    private static BassInputs Inputs(int[] distinct, int[] levels, int? outputs)
    {
        BassData data = Data(distinct, levels, outputs);
        return BassInputs.Prepare(data, BassInputTransform.Fit(data));
    }

    private static BassData Data(int[] distinct, int[] levels, int? outputs)
    {
        const int n = 12;
        var predictors = new List<BassPredictor>();
        for (int v = 0; v < distinct.Length; v++)
        {
            int d = distinct[v], shift = v;
            predictors.Add(new ContinuousPredictor($"x{v}", [.. Enumerable.Range(0, n).Select(i => (double)((i + shift) % d))]));
        }
        for (int v = 0; v < levels.Length; v++)
        {
            int l = levels[v], shift = v;
            predictors.Add(new CategoricalPredictor($"z{v}", [.. Enumerable.Range(0, n).Select(i => (i + shift) % l)],
                                                    [.. Enumerable.Range(0, l).Select(c => $"L{c}")]));
        }
        BassResponse response;
        if (outputs is int m)
        {
            var y = new double[n, m];
            for (int i = 0; i < n; i++) for (int r = 0; r < m; r++) y[i, r] = Math.Sin(i + 0.3 * r);
            response = new FunctionalResponse(y, [new ContinuousPredictor("t", [.. Enumerable.Range(0, m).Select(r => (double)r)])]);
        }
        else response = new ScalarResponse([.. Enumerable.Range(0, n).Select(i => Math.Sin(i))]);
        return new BassData(predictors, response);
    }

    private static Space Make(string name, int[] distinct, int[] levels, int? outputs, int capC, int capL, int capF,
                              CategoricalSubsetMeasure measure = CategoricalSubsetMeasure.SizeUniform)
    {
        BassInputs inputs = Inputs(distinct, levels, outputs);
        var prior = new BassStructurePrior(inputs, capC, capL, capF, measure, new BassCountPrior());
        return new(name, inputs, prior, new BassTermProposal(prior, 5, 5), Enumerate(prior));
    }

    public static TheoryData<string> SpaceNames() => new() { "continuous", "continuous-3", "mixed", "mixed-uniform", "functional", "all" };

    private static Space Named(string name) => name switch
    {
        "continuous" => Make(name, [3, 2, 4], [], null, 2, 0, 0),
        "continuous-3" => Make(name, [3, 2, 2, 2], [], null, 3, 0, 0),
        "mixed" => Make(name, [3, 2, 2], [3, 4], null, 2, 2, 0),
        "mixed-uniform" => Make(name, [3, 2, 2], [3, 4], null, 2, 2, 0, CategoricalSubsetMeasure.Uniform),
        "functional" => Make(name, [3, 2], [], 3, 2, 0, 1),
        "all" => Make(name, [2, 2], [3], 2, 2, 1, 1),
        _ => throw new ArgumentException(name),
    };

    // The term space, from the definitions: per part every order from its smallest to its cap, every variable set,
    // every sign and knot or proper subset; the product of the parts, less the all-empty term.
    private static List<BassTerm> Enumerate(BassStructurePrior prior)
    {
        List<HingeFactor[]> continuous = HingeParts(prior, BassPart.Continuous);
        List<HingeFactor[]> functional = HingeParts(prior, BassPart.Functional);
        var categorical = new List<CategoricalFactor[]>();
        if (!prior.Has(BassPart.Categorical)) categorical.Add([]);
        else
            for (int k = prior.MinimumOrder(BassPart.Categorical); k <= prior.Cap(BassPart.Categorical); k++)
                foreach (int[] set in Combinations(prior.Variables(BassPart.Categorical), k))
                    foreach (CategoricalFactor[] choice in Product(set.Select(v => ProperSubsets(v, prior[BassPart.Categorical]!.Sizes[v])).ToList()))
                        categorical.Add(choice);
        var terms = new List<BassTerm>();
        foreach (HingeFactor[] c in continuous)
            foreach (CategoricalFactor[] l in categorical)
                foreach (HingeFactor[] f in functional)
                    if (c.Length + l.Length + f.Length > 0) terms.Add(new BassTerm(c, l, f));
        Assert.Equal(terms.Count, terms.Distinct().Count());
        return terms;
    }

    private static List<HingeFactor[]> HingeParts(BassStructurePrior prior, BassPart part)
    {
        var parts = new List<HingeFactor[]>();
        if (!prior.Has(part)) { parts.Add([]); return parts; }
        int[] sizes = prior[part]!.Sizes;
        for (int k = prior.MinimumOrder(part); k <= prior.Cap(part); k++)
            foreach (int[] set in Combinations(sizes.Length, k))
                foreach (HingeFactor[] choice in Product(set.Select(v =>
                             (from s in new[] { -1, 1 } from knot in Enumerable.Range(0, sizes[v]) select new HingeFactor(v, s, knot)).ToList()).ToList()))
                    parts.Add(choice);
        return parts;
    }

    private static List<CategoricalFactor> ProperSubsets(int variable, int levels)
        => [.. Enumerable.Range(1, (1 << levels) - 2).Select(mask => new CategoricalFactor(variable, Enumerable.Range(0, levels).Where(l => (mask >> l & 1) != 0)))];

    private static IEnumerable<int[]> Combinations(int n, int k)
    {
        if (k == 0) { yield return []; yield break; }
        for (int first = 0; first < n; first++)
            foreach (int[] rest in Combinations(n, k - 1))
                if (rest.Length == 0 || rest[0] > first) yield return [first, .. rest];
    }

    private static IEnumerable<T[]> Product<T>(List<List<T>> choices)
    {
        if (choices.Count == 0) { yield return []; yield break; }
        foreach (T head in choices[0])
            foreach (T[] tail in Product(choices.Skip(1).ToList()))
                yield return [head, .. tail];
    }

    // The all-empty combination's prior mass where the fit has several parts: each part's order 0 with 1/(cap + 1).
    private static double EmptyPriorMass(BassStructurePrior prior)
    {
        BassPart[] present = [.. Enum.GetValues<BassPart>().Where(prior.Has)];
        return present.Length < 2 ? 0.0 : present.Aggregate(1.0, (m, p) => m / (prior.Cap(p) + 1));
    }

    // A structure of M distinct terms drawn from the space, skewed toward its first terms so some variables dominate.
    private static (List<BassTerm> Terms, BassTermCounts Counts) Structure(Space s, int m, Xoshiro256PlusPlus rng)
    {
        var chosen = new List<BassTerm>();
        var counts = new BassTermCounts(s.Prior);
        while (chosen.Count < m)
        {
            double u = rng.NextDouble();
            BassTerm t = s.Terms[(int)(u * u * s.Terms.Count)];
            if (chosen.Contains(t)) continue;
            chosen.Add(t);
            counts.Add(t);
        }
        return (chosen, counts);
    }

    [Theory]
    [MemberData(nameof(SpaceNames))]
    public void TermPrior_SumsToOne_OverTheTermSpace(string name)
    {
        Space s = Named(name);
        double total = s.Terms.Sum(t => Math.Exp(s.Prior.LogPrior(t)));
        double empty = EmptyPriorMass(s.Prior);
        _out.WriteLine($"{name}: {s.Terms.Count} terms, Σπ = {total:R}, empty {empty:R}");
        Assert.True(Math.Abs(total + empty - 1.0) <= 1e-12, $"{name}: Σπ + empty = {total + empty:R}");
    }

    [Fact]
    public void TermPrior_IsD01sConstant()
    {
        // One part: Francom et al.'s c_K = (1/2)^K C(p, K)⁻¹ Π 1/n_v, times 1/K_max.
        Space c = Named("continuous");
        var term = new BassTerm([new HingeFactor(0, 1, 1), new HingeFactor(2, -1, 3)]);
        Close(Math.Log(1.0 / 2 * 1.0 / 3 * 0.25 * (1.0 / 3) * (1.0 / 4)), c.Prior.LogPrior(term));

        // Several parts: each part's order uniform on {0, …, cap}, then as above per part (the stage-4 oracle's
        // cont_lp + cat_lp). R's subset measure: size uniform on {1, …, n − 1}, then uniform among C(n, size).
        Space mixed = Named("mixed");
        var hinge = new HingeFactor(1, 1, 0);
        var subset = new CategoricalFactor(1, [0, 2]);
        Close(Math.Log(1.0 / 3 * (1.0 / 3) * 0.5 * 0.5 * (1.0 / 3) * 0.5 * (1.0 / 3) * (1.0 / 6)),
              mixed.Prior.LogPrior(new BassTerm([hinge], [subset])));
        Close(Math.Log(1.0 / 3 * (1.0 / 3) * 0.5 * 0.5 * (1.0 / 3)), mixed.Prior.LogPrior(new BassTerm([hinge])));
        Close(Math.Log(1.0 / 3 * (1.0 / 3) * 0.5 * (1.0 / 3) * (1.0 / 6)), mixed.Prior.LogPrior(new BassTerm(categorical: [subset])));

        // The paper's measure: uniform over the 2⁴ − 2 proper subsets.
        Space uniform = Named("mixed-uniform");
        Close(Math.Log(1.0 / 3 * (1.0 / 3) * 0.5 * (1.0 / 14)), uniform.Prior.LogPrior(new BassTerm(categorical: [subset])));
    }

    [Theory]
    [MemberData(nameof(SpaceNames))]
    public void BirthProposal_SumsToOne_ForEveryCountState(string name)
    {
        Space s = Named(name);
        var rng = new Xoshiro256PlusPlus(seed: 31 + name.Length);
        double worst = 0.0;
        for (int m = 0; m <= 7; m++)
        {
            (_, BassTermCounts counts) = Structure(s, m, rng);
            double total = s.Terms.Sum(t => Math.Exp(s.Proposal.LogBirth(t, counts)));
            double empty = Math.Exp(s.Proposal.LogEmptyBirth(counts));
            if (EmptyPriorMass(s.Prior) == 0.0) Assert.Equal(0.0, empty);
            worst = Math.Max(worst, Math.Abs(total + empty - 1.0));
        }
        _out.WriteLine($"{name}: worst |Σq − 1| = {worst:E2}");
        Assert.True(worst <= 1e-12, $"{name}: |Σq − 1| = {worst:E2}");
    }

    [Theory]
    [InlineData("mixed")]
    [InlineData("mixed-uniform")]
    [InlineData("all")]
    [InlineData("continuous-3")]
    public void BirthDraws_FollowTheirProbabilities(string name)
    {
        Space s = Named(name);
        var rng = new Xoshiro256PlusPlus(seed: 7);
        (_, BassTermCounts counts) = Structure(s, 5, rng);
        const int draws = 400_000;
        var seen = new Dictionary<BassTerm, int>();
        int empty = 0;
        for (int d = 0; d < draws; d++)
        {
            if (s.Proposal.DrawBirth(rng, counts) is not { } t) { empty++; continue; }
            seen[t] = seen.GetValueOrDefault(t) + 1;
        }
        Assert.All(seen.Keys, t => Assert.Contains(t, s.Terms));

        // Pearson's statistic over the term space and the empty draw; under the stated law it is χ² with as many
        // degrees of freedom as cells less one.
        double chi2 = 0.0;
        int cells = 0;
        foreach (BassTerm t in s.Terms)
        {
            double expected = draws * Math.Exp(s.Proposal.LogBirth(t, counts)), o = seen.GetValueOrDefault(t);
            chi2 += (o - expected) * (o - expected) / expected;
            cells++;
        }
        double emptyExpected = draws * Math.Exp(s.Proposal.LogEmptyBirth(counts));
        if (emptyExpected > 0) { chi2 += (empty - emptyExpected) * (empty - emptyExpected) / emptyExpected; cells++; }
        else Assert.Equal(0, empty);
        int df = cells - 1;
        double z = (chi2 - df) / Math.Sqrt(2.0 * df);
        _out.WriteLine($"{name}: {cells} cells, χ² = {chi2:F1} on {df}, z = {z:F2}");
        Assert.True(z < 5.0, $"{name}: χ² = {chi2:F1} on {df} (z = {z:F2})");
    }

    [Fact]
    public void BirthProbability_ByHand()
    {
        // One part, p = 3, K_max = 2. The structure: x0, x0·x1, x0·x2.
        Space s = Named("continuous");
        var counts = new BassTermCounts(s.Prior);
        counts.Add(new BassTerm([new HingeFactor(0, 1, 0)]));
        counts.Add(new BassTerm([new HingeFactor(0, 1, 0), new HingeFactor(1, -1, 1)]));
        counts.Add(new BassTerm([new HingeFactor(0, -1, 2), new HingeFactor(2, 1, 0)]));
        Assert.Equal(3, counts.Terms);
        Assert.Equal(2, counts.Order(BassPart.Continuous, 2));
        Assert.Equal(3, counts.Appearances(BassPart.Continuous, 0));

        // I(2) = (5 + 2)/(2·5 + 3); z = 5 + (3, 1, 1); the set {1, 2} by brute force over its two draw orders.
        double i2 = 7.0 / 13, z0 = 8, z1 = 6, z2 = 6, total = z0 + z1 + z2;
        double wallenius = z1 / total * (z2 / (total - z1)) + z2 / total * (z1 / (total - z2));
        var born = new BassTerm([new HingeFactor(1, 1, 1), new HingeFactor(2, -1, 2)]);
        Close(Math.Log(i2 * wallenius * 0.5 * 0.5 * (1.0 / 2) * (1.0 / 4)), s.Proposal.LogBirth(born, counts));

        // Order 1: I(1) = 6/13, the variable uniform over 3, whatever the counts (R's genCandBasis).
        Close(Math.Log(6.0 / 13 * (1.0 / 3) * 0.5 * (1.0 / 4)), s.Proposal.LogBirth(new BassTerm([new HingeFactor(2, 1, 0)]), counts));

        // With one part this is coherent in R too: ln q − ln π is logProbChangeMod's ln I + ln W + ln C(p, K) + ln K_max.
        Close(Math.Log(i2) + Math.Log(wallenius) + Math.Log(3) + Math.Log(2),
              s.Proposal.LogBirth(born, counts) - s.Prior.LogPrior(born));

        // The birth's within-move ratio: the reverse death's 1/(M + 1) over the draw.
        Close(-Math.Log(4) - s.Proposal.LogBirth(born, counts), s.Proposal.LogBirthRatio(born, counts));
    }

    [Theory]
    [MemberData(nameof(SpaceNames))]
    public void DeathRatio_ReversesTheBirthRatio(string name)
    {
        Space s = Named(name);
        var rng = new Xoshiro256PlusPlus(seed: 101);
        double worst = 0.0;
        for (int trial = 0; trial < 200; trial++)
        {
            (List<BassTerm> terms, BassTermCounts counts) = Structure(s, (int)(rng.NextDouble() * 6), rng);
            BassTerm born = s.Terms[(int)(rng.NextDouble() * s.Terms.Count)];
            if (terms.Contains(born)) continue;
            double birth = s.Proposal.LogBirthRatio(born, counts);
            BassTermCounts after = counts.Clone();
            after.Add(born);
            double death = s.Proposal.LogDeathRatio(born, after);
            worst = Math.Max(worst, Math.Abs(birth + death));

            // The death's reverse reads the counts without the term, which are the birth's counts.
            Close(s.Proposal.LogBirth(born, counts), death - Math.Log(after.Terms));
        }
        Assert.True(worst <= 1e-13, $"{name}: worst |birth + death| = {worst:E2}");
    }

    [Theory]
    [InlineData("mixed")]
    [InlineData("mixed-uniform")]
    [InlineData("all")]
    public void ChangeProposal_IsItsFactorsPriorRedrawn(string name)
    {
        Space s = Named(name);
        var rng = new Xoshiro256PlusPlus(seed: 5);
        foreach (BassTerm term in s.Terms.Where(t => t.Categorical.Count > 0 && t.Continuous.Count > 0).Take(3))
        {
            // The exact law: one of the F factors uniformly, that factor's value from its prior.
            Dictionary<BassTerm, double> law = ChangeLaw(s, term);
            var backward = new Dictionary<BassTerm, Dictionary<BassTerm, double>>();
            Assert.True(Math.Abs(law.Values.Sum() - 1.0) <= 1e-13);
            const int draws = 100_000;
            var seen = new Dictionary<BassTerm, int>();
            for (int d = 0; d < draws; d++)
            {
                BassTerm next = s.Proposal.DrawChange(rng, term, out double logRatio);
                seen[next] = seen.GetValueOrDefault(next) + 1;
                Assert.True(law.ContainsKey(next), $"{next} is not one factor away from {term}");
                // ln q(back)/q(forward), against the exact law both ways; and the prior's ratio cancels it.
                if (!backward.TryGetValue(next, out Dictionary<BassTerm, double>? back)) backward[next] = back = ChangeLaw(s, next);
                double expected = next.Equals(term) ? 0.0 : Math.Log(back[term]) - Math.Log(law[next]);
                Assert.True(Math.Abs(logRatio - expected) <= 1e-12, $"{term} → {next}: {logRatio:R} against {expected:R}");
                Assert.True(Math.Abs(logRatio - (s.Prior.LogPrior(term) - s.Prior.LogPrior(next))) <= 1e-12);
            }
            double chi2 = law.Sum(kv => Math.Pow(seen.GetValueOrDefault(kv.Key) - draws * kv.Value, 2) / (draws * kv.Value));
            int df = law.Count - 1;
            double z = (chi2 - df) / Math.Sqrt(2.0 * df);
            _out.WriteLine($"{name}, {term}: {law.Count} candidates, χ² = {chi2:F1} on {df}, z = {z:F2}");
            Assert.True(z < 5.0);
        }
    }

    private static Dictionary<BassTerm, double> ChangeLaw(Space s, BassTerm term)
    {
        var law = new Dictionary<BassTerm, double>();
        int factors = term.Continuous.Count + term.Categorical.Count + term.Functional.Count;
        void Put(BassTerm t, double p) => law[t] = law.GetValueOrDefault(t) + p / factors;
        for (int f = 0; f < term.Continuous.Count; f++)
        {
            int v = term.Continuous[f].Variable, knots = s.Inputs.ContinuousKnots(v).Length;
            foreach (int sign in new[] { -1, 1 })
                for (int k = 0; k < knots; k++)
                {
                    HingeFactor[] h = [.. term.Continuous];
                    h[f] = new HingeFactor(v, sign, k);
                    Put(new BassTerm(h, term.Categorical, term.Functional), 1.0 / (2 * knots));
                }
        }
        for (int f = 0; f < term.Categorical.Count; f++)
        {
            int v = term.Categorical[f].Variable, levels = s.Inputs.LevelCount(v);
            foreach (CategoricalFactor sub in ProperSubsets(v, levels))
            {
                CategoricalFactor[] c = [.. term.Categorical];
                c[f] = sub;
                double p = s.Prior.Subsets == CategoricalSubsetMeasure.Uniform
                    ? 1.0 / ((1 << levels) - 2)
                    : 1.0 / (levels - 1) / Math.Exp(SpecialFunctions.LogGamma(levels + 1.0) - SpecialFunctions.LogGamma(sub.Levels.Count + 1.0)
                                                    - SpecialFunctions.LogGamma(levels - sub.Levels.Count + 1.0));
                Put(new BassTerm(term.Continuous, c, term.Functional), p);
            }
        }
        for (int f = 0; f < term.Functional.Count; f++)
        {
            int v = term.Functional[f].Variable, knots = s.Inputs.FunctionalKnots(v).Length;
            foreach (int sign in new[] { -1, 1 })
                for (int k = 0; k < knots; k++)
                {
                    HingeFactor[] h = [.. term.Functional];
                    h[f] = new HingeFactor(v, sign, k);
                    Put(new BassTerm(term.Continuous, term.Categorical, h), 1.0 / (2 * knots));
                }
        }
        return law;
    }

    [Fact]
    public void CountPrior_IsTheNegativeBinomial()
    {
        Space s = Named("continuous");
        double h1 = s.Prior.Count.Shape, h2 = s.Prior.Count.Rate;
        // p(M) = Γ(h₁ + M)/(Γ(h₁) M!) · (h₂/(h₂ + 1))^h₁ · (1/(h₂ + 1))^M: the collection weight less ln M!.
        double total = 0.0;
        for (int m = 0; m < 400; m++)
        {
            double logP = s.Prior.LogCollectionWeight(m) - SpecialFunctions.LogGamma(m + 1.0) + h1 * Math.Log(h2 / (h2 + 1));
            total += Math.Exp(logP);
            // The difference of two ln Γ near 2000 at M = 400 keeps about 13 digits of each.
            double ratio = s.Prior.LogCollectionWeight(m + 1) - s.Prior.LogCollectionWeight(m) - Math.Log(m + 1.0);
            Assert.True(Math.Abs(ratio - s.Prior.LogCountRatio(m)) <= 1e-13 + 1e-15 * Math.Abs(s.Prior.LogCollectionWeight(m + 1)),
                        $"M = {m}: {s.Prior.LogCountRatio(m):R} against {ratio:R}");
        }
        Assert.True(Math.Abs(total - 1.0) <= 1e-12, $"Σ p(M) = {total:R}");
        // D02's birth count factor at M = 0 under the defaults h₁ = h₂ = 10: 10/11.
        Close(Math.Log(10.0 / 11), s.Prior.LogCountRatio(0));
    }

    [Fact]
    public void Prior_ReadsTheConfigsCapsAndOrders()
    {
        BassData data = Data([3, 2], [3], null);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var config = new BassConfig();
        BassResolved resolved = BassResolved.Resolve(config, data);
        var prior = BassStructurePrior.From(config, resolved, inputs);
        Assert.Equal(2, prior.Cap(BassPart.Continuous));      // 3 declared, 2 continuous predictors
        Assert.Equal(1, prior.Cap(BassPart.Categorical));     // 3 declared, 1 categorical predictor
        Assert.False(prior.Has(BassPart.Functional));
        Assert.Equal(0, prior.MinimumOrder(BassPart.Continuous));
        Assert.Equal(0, prior.MinimumOrder(BassPart.Categorical));
        Assert.Equal(CategoricalSubsetMeasure.SizeUniform, prior.Subsets);
        var proposal = BassTermProposal.From(config, prior);
        Assert.Equal(5, proposal.OrderWeight);
        Assert.Equal(5, proposal.VariableWeight);

        // One part: orders from 1.
        var single = new BassStructurePrior(Inputs([3, 2], [], null), 2, 0, 0, CategoricalSubsetMeasure.SizeUniform, new BassCountPrior());
        Assert.Equal(1, single.MinimumOrder(BassPart.Continuous));
        Assert.Equal(double.NegativeInfinity, new BassTermProposal(single, 5, 5).LogEmptyBirth(new BassTermCounts(single)));
    }

    [Fact]
    public void Counts_AddAndRemove_AndRefusals()
    {
        Space s = Named("mixed");
        var counts = new BassTermCounts(s.Prior);
        var rng = new Xoshiro256PlusPlus(seed: 3);
        (List<BassTerm> terms, BassTermCounts built) = Structure(s, 6, rng);
        foreach (BassTerm t in terms) counts.Add(t);
        foreach (BassTerm t in terms.AsEnumerable().Reverse()) counts.Remove(t);
        Assert.Equal(0, counts.Terms);
        for (int k = 0; k <= 2; k++) Assert.Equal(0, counts.Order(BassPart.Continuous, k));
        for (int v = 0; v < 3; v++) Assert.Equal(0, counts.Appearances(BassPart.Continuous, v));
        Assert.Throws<ArgumentException>(() => counts.Remove(terms[0]));
        Assert.Equal(6, built.Terms);

        // A term held nowhere has no reverse birth.
        BassTerm absent = s.Terms.First(t => !terms.Contains(t) && t.Continuous.Count == 2);
        var lone = new BassTermCounts(s.Prior);
        lone.Add(new BassTerm(categorical: [new CategoricalFactor(0, [1])]));
        Assert.Throws<ArgumentException>(() => s.Proposal.LogDeathRatio(absent, lone));

        // Outside the term space.
        Assert.Throws<ArgumentException>(() => s.Prior.LogPrior(new BassTerm([new HingeFactor(0, 1, 0), new HingeFactor(1, 1, 0), new HingeFactor(2, 1, 0)])));   // above the cap
        Assert.Throws<ArgumentException>(() => s.Prior.LogPrior(new BassTerm([new HingeFactor(3, 1, 0)])));                  // no such variable
        Assert.Throws<ArgumentException>(() => s.Prior.LogPrior(new BassTerm([new HingeFactor(0, 1, 3)])));                  // no such knot
        Assert.Throws<ArgumentException>(() => s.Prior.LogPrior(new BassTerm(categorical: [new CategoricalFactor(0, [0, 1, 2])])));   // not proper
        Assert.Throws<ArgumentException>(() => s.Prior.LogPrior(new BassTerm(categorical: [new CategoricalFactor(1, [4])])));       // no such level
        Assert.Throws<ArgumentException>(() => s.Prior.LogPrior(new BassTerm(functional: [new HingeFactor(0, 1, 0)])));            // a scalar fit
        Assert.Throws<ArgumentException>(() => s.Proposal.LogBirth(terms[0], new BassTermCounts(Named("mixed").Prior)));         // another prior's counts
    }

    [Fact]
    public void Draws_AreExact()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        const int draws = 300_000;
        foreach (int count in new[] { 1, 3, 7 })
        {
            var seen = new int[count];
            for (int d = 0; d < draws; d++) seen[BassDraws.Below(rng, count)]++;
            double chi2 = seen.Sum(o => Math.Pow(o - (double)draws / count, 2) / ((double)draws / count));
            Assert.True(count == 1 || chi2 < count - 1 + 6 * Math.Sqrt(2.0 * (count - 1)), $"count {count}: χ² = {chi2:F1}");
        }
        for (int d = 0; d < 1000; d++)
        {
            int big = BassDraws.Below(rng, int.MaxValue);
            Assert.InRange(big, 0, int.MaxValue - 1);
        }
        int plus = 0;
        for (int d = 0; d < draws; d++) plus += BassDraws.Sign(rng) > 0 ? 1 : 0;
        Assert.InRange((double)plus, draws / 2.0 - 3 * Math.Sqrt(draws), draws / 2.0 + 3 * Math.Sqrt(draws));
    }

    private static void Close(double expected, double actual, double tolerance = 1e-13)
        => Assert.True(Math.Abs(expected - actual) <= tolerance * Math.Max(1.0, Math.Abs(expected)), $"{actual:R} against {expected:R}");
}
