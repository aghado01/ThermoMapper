using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's exact sampler of the prior on F (stage 4's S7; D02, D05). Its term draws against the term measure the prior
/// scores, over every term of three term spaces (two continuous predictors with an order-2 cap and a four-level factor
/// under each subset measure; a direct functional fit). Its collections against the prior on F enumerated exactly: S6's
/// constant-likelihood fixture (<see cref="BassSweepTests.EnumerateF"/>, 2,263 states) for the count law, the order law,
/// the law over every state, τ's prior and the acceptance rate; a mixed fixture with a categorical predictor for the
/// joint order law and the subsets; and the exact fits Jeffreys excludes. Fixed g holds τ, and a draw refuses past its
/// attempt limit.
/// </summary>
public sealed class BassPriorSamplerTests
{
    private readonly ITestOutputHelper _out;
    public BassPriorSamplerTests(ITestOutputHelper output) => _out = output;

    // x₀ and x₁ on {0, ½, 1} and a four-level factor, each combination once (n = 36).
    private static BassData TermData()
    {
        var rng = new Xoshiro256PlusPlus(3);
        var x0 = new List<double>();
        var x1 = new List<double>();
        var z = new List<int>();
        var y = new List<double>();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 4; k++)
                {
                    x0.Add(i / 2.0);
                    x1.Add(j / 2.0);
                    z.Add(k);
                    y.Add(i + 0.5 * j + (k == 2 ? 1 : 0) + rng.NextDouble());
                }
        return new BassData([new ContinuousPredictor("x0", [.. x0]), new ContinuousPredictor("x1", [.. x1]),
                             new CategoricalPredictor("z", [.. z], ["a", "b", "c", "d"])], new ScalarResponse([.. y]));
    }

    // x on five values and a three-level factor, each combination twice (n = 30).
    internal static BassData MixedData(int seed = 7)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[30];
        var z = new int[30];
        var y = new double[30];
        for (int i = 0; i < 30; i++)
        {
            x[i] = (i % 5) / 4.0;
            z[i] = (i / 5) % 3;
            y[i] = 1 + x[i] + (z[i] == 1 ? 0.5 : 0) + 0.3 * rng.NextDouble();
        }
        return new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", z, ["a", "b", "c"])], new ScalarResponse(y));
    }

    // x on five values (n = 10), one output variable t on seven locations.
    private static BassData FunctionalData()
    {
        var rng = new Xoshiro256PlusPlus(11);
        double[] t = [.. Enumerable.Range(0, 7).Select(r => r / 6.0)];
        var x = new double[10];
        var y = new double[10, 7];
        for (int i = 0; i < 10; i++)
        {
            x[i] = (i % 5) / 4.0;
            for (int r = 0; r < 7; r++) y[i, r] = x[i] * t[r] + 0.1 * rng.NextDouble();
        }
        return new BassData([new ContinuousPredictor("x", x)], new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    [Theory]
    [InlineData("mixed", CategoricalSubsetMeasure.SizeUniform)]
    [InlineData("mixed", CategoricalSubsetMeasure.Uniform)]
    [InlineData("functional", CategoricalSubsetMeasure.SizeUniform)]
    public void DrawTerm_FollowsTheTermMeasure(string fixture, CategoricalSubsetMeasure subsets)
    {
        BassModel model = fixture == "mixed"
            ? BassModel.Create(new BassConfig
            {
                Structure = new BassStructureSpec { MaxInteraction = 2, MaxCategoricalInteraction = 1, CategoricalSubsets = subsets },
            }, TermData())
            : BassModel.Create(new BassConfig { Structure = new BassStructureSpec { MaxInteraction = 1, MaxFunctionalInteraction = 1 } },
                               FunctionalData());
        List<BassTerm> space = BassTermSpace.Terms(model);
        var index = space.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        var p = space.Select(t => Math.Exp(model.Prior.LogPrior(t))).ToList();
        double intercept = 1.0;   // every part empty: Πₚ 1/(cap + 1)
        foreach (BassPart part in Enum.GetValues<BassPart>())
            if (model.Prior.Has(part)) intercept /= model.Prior.Cap(part) + 1;
        Assert.Equal(1.0, p.Sum() + intercept, 12);   // the measure is a probability over the space and the intercept

        const int draws = 300_000;
        var observed = new long[space.Count + 1];
        var rng = new Xoshiro256PlusPlus(17);
        for (int d = 0; d < draws; d++)
        {
            BassTerm? term = model.Prior.DrawTerm(rng);
            observed[term is null ? space.Count : index[term]]++;
        }
        (double chi2, int df) = ChiSquare([.. p.Select(v => v * draws), intercept * draws], observed);
        _out.WriteLine($"{fixture}, {subsets}: {space.Count} terms; χ² = {chi2:F1} on {df} degrees of freedom");
        Assert.True(chi2 <= df + 6 * Math.Sqrt(2.0 * df), $"χ² = {chi2:F1} on {df}");
    }

    [Fact]
    public void Draws_MatchThePriorOnF()
    {
        var model = BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data());
        List<(BassTerm[] Terms, double LogPrior)> f = BassSweepTests.EnumerateF(model);
        Assert.Equal(2263, f.Count);
        Assert.Equal(f.Select(s => Key(s.Terms)).OrderBy(k => k), BassTermSpace.F(model).Select(s => Key(s.Terms)).OrderBy(k => k));

        var sampler = new BassPriorSampler(model);
        (Dictionary<string, double> law, double[] countLaw) = Law(f, model.MaxBasis);
        var orderLaw = new Dictionary<string, double>();
        foreach ((BassTerm[] terms, double lp) in f) orderLaw[OrderKey(terms)] = orderLaw.GetValueOrDefault(OrderKey(terms)) + law[Key(terms)];

        const int draws = 150_000;
        var states = new Dictionary<string, long>();
        var counts = new long[countLaw.Length];
        var orders = orderLaw.Keys.ToDictionary(k => k, _ => 0L);
        var logTau = new double[draws];
        var size = new double[draws];
        var rng = new Xoshiro256PlusPlus(23);
        for (int d = 0; d < draws; d++)
        {
            (BassState s, ITargetEvaluation e) = sampler.Draw(rng);
            Assert.False(double.IsNegativeInfinity(e.LogPrior));
            states[Key(s.Terms)] = states.GetValueOrDefault(Key(s.Terms)) + 1;
            counts[s.Count]++;
            orders[OrderKey(s.Terms)]++;
            logTau[d] = s.LogTau;
            size[d] = s.Count;
        }
        Assert.Equal(draws, sampler.Accepted);
        Assert.True(states.Keys.All(law.ContainsKey));

        for (int m = 0; m < countLaw.Length; m++)
        {
            double observed = (double)counts[m] / draws, error = Math.Sqrt(countLaw[m] * (1 - countLaw[m]) / draws);
            _out.WriteLine($"M = {m}: {observed:F4} against {countLaw[m]:F4} (± {error:F4})");
            Assert.True(Math.Abs(observed - countLaw[m]) <= 5 * error, $"M = {m}: {observed:F4} against {countLaw[m]:F4}");
        }
        foreach ((string key, double q) in orderLaw)
        {
            double observed = (double)orders[key] / draws, error = Math.Sqrt(q * (1 - q) / draws);
            Assert.True(Math.Abs(observed - q) <= 5 * error + 1e-9, $"orders ({key}): {observed:F4} against {q:F4}");
        }
        (double chi2, int df) = ChiSquare([.. f.Select(s => law[Key(s.Terms)] * draws)], [.. f.Select(s => states.GetValueOrDefault(Key(s.Terms)))]);
        _out.WriteLine($"χ² over F's {f.Count} states: {chi2:F1} on {df}");
        Assert.True(chi2 <= df + 6 * Math.Sqrt(2.0 * df), $"χ² = {chi2:F1} on {df}");

        // τ from its prior, independently of the structure: E[ln τ] = ψ(½) − ln(N/2), Var[ln τ] = ψ′(½) = π²/2.
        double expected = SpecialFunctions.Digamma(0.5) - Math.Log(0.5 * 18), meanLogTau = logTau.Average();
        Assert.True(Math.Abs(meanLogTau - expected) <= 5 * Math.Sqrt(Math.PI * Math.PI / 2 / draws), $"E[ln τ] = {meanLogTau:F4} against {expected:F4}");
        double meanTau = logTau.Average(Math.Exp), sdTau = Math.Sqrt(2.0) / 18;   // τ ~ Gamma(½, rate 9): mean 1/18, sd √2/18
        Assert.True(Math.Abs(meanTau - 1.0 / 18) <= 5 * sdTau / Math.Sqrt(draws), $"E[τ] = {meanTau:G5} against {1.0 / 18:G5}");
        Assert.True(Math.Abs(Correlation(logTau, size)) <= 5 / Math.Sqrt(draws), "τ is drawn apart from the structure");

        // The acceptance rate: F's share of the truncated proposal, Σ_F w(s) / Σ_{M ≤ M_max} p(M) with p(M) = w(M)/M!.
        double z = f.Max(s => s.LogPrior);
        double inF = f.Sum(s => Math.Exp(s.LogPrior - z));
        double proposal = Enumerable.Range(0, model.MaxBasis + 1)
            .Sum(m => Math.Exp(model.Prior.LogCollectionWeight(m) - SpecialFunctions.LogGamma(m + 1.0) - z));
        double rate = inF / proposal, rateError = Math.Sqrt(rate * (1 - rate) / sampler.Attempts);
        _out.WriteLine($"acceptance {sampler.AcceptanceRate:F4} against {rate:F4} (± {rateError:F4}); refusals " +
                       string.Join(", ", sampler.Refusals.Entries.Where(e => e.Count > 0).Select(e => $"{e.Category} {e.Count}")));
        Assert.True(Math.Abs(sampler.AcceptanceRate - rate) <= 5 * rateError, $"acceptance {sampler.AcceptanceRate:F4} against {rate:F4}");
        Assert.Equal(sampler.Attempts - sampler.Accepted, sampler.Refusals.Total + sampler.InterceptTerms);
        Assert.Equal(0, sampler.InterceptTerms);   // one part: no draw comes out empty
        Assert.Equal(0, sampler.Refusals.Count(BassRefusal.Cap));   // the count law is truncated at the cap
    }

    [Theory]
    [InlineData(CategoricalSubsetMeasure.SizeUniform)]
    [InlineData(CategoricalSubsetMeasure.Uniform)]
    public void Draws_MatchThePriorOnF_WithACategoricalPart(CategoricalSubsetMeasure subsets)
    {
        var model = BassModel.Create(new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 1, MaxCategoricalInteraction = 1, MaxBasis = 2, CategoricalSubsets = subsets },
        }, MixedData());
        List<(BassTerm[] Terms, double LogPrior)> f = BassTermSpace.F(model);
        (Dictionary<string, double> law, double[] countLaw) = Law(f, model.MaxBasis);
        // The joint order law: each term's (J, L), as a multiset over the structure.
        var orderLaw = new Dictionary<string, double>();
        var subsetLaw = new double[3];   // a categorical factor's subset size, weighted by the factors a structure holds
        foreach ((BassTerm[] terms, double _) in f)
        {
            double q = law[Key(terms)];
            orderLaw[JointOrderKey(terms)] = orderLaw.GetValueOrDefault(JointOrderKey(terms)) + q;
            foreach (BassTerm t in terms)
                foreach (CategoricalFactor c in t.Categorical) subsetLaw[c.Levels.Count] += q;
        }

        var sampler = new BassPriorSampler(model);
        const int draws = 150_000;
        var states = new Dictionary<string, long>();
        var counts = new long[countLaw.Length];
        var orders = orderLaw.Keys.ToDictionary(k => k, _ => 0L);
        var sizes = new double[3];
        var rng = new Xoshiro256PlusPlus(29);
        for (int d = 0; d < draws; d++)
        {
            BassState s = sampler.Draw(rng).State;
            states[Key(s.Terms)] = states.GetValueOrDefault(Key(s.Terms)) + 1;
            counts[s.Count]++;
            orders[JointOrderKey(s.Terms)]++;
            foreach (BassTerm t in s.Terms)
                foreach (CategoricalFactor c in t.Categorical) sizes[c.Levels.Count]++;
        }
        _out.WriteLine($"{subsets}: {f.Count} states in F; acceptance {sampler.AcceptanceRate:F4}, {sampler.InterceptTerms} intercept terms");
        Assert.True(sampler.InterceptTerms > 0);   // two parts: the empty draw has mass 1/4 per term
        for (int m = 0; m < countLaw.Length; m++)
        {
            double observed = (double)counts[m] / draws, error = Math.Sqrt(countLaw[m] * (1 - countLaw[m]) / draws);
            Assert.True(Math.Abs(observed - countLaw[m]) <= 5 * error, $"M = {m}: {observed:F4} against {countLaw[m]:F4}");
        }
        foreach ((string key, double q) in orderLaw)
        {
            double observed = (double)orders[key] / draws, error = Math.Sqrt(q * (1 - q) / draws);
            _out.WriteLine($"  orders [{key}]: {observed:F4} against {q:F4}");
            Assert.True(Math.Abs(observed - q) <= 5 * error + 1e-9, $"orders [{key}]: {observed:F4} against {q:F4}");
        }
        for (int s = 1; s < 3; s++)
            Assert.True(Math.Abs(sizes[s] / draws - subsetLaw[s]) <= 0.01, $"subsets of size {s}: {sizes[s] / draws:F4} against {subsetLaw[s]:F4}");
        (double chi2, int df) = ChiSquare([.. f.Select(s => law[Key(s.Terms)] * draws)], [.. f.Select(s => states.GetValueOrDefault(Key(s.Terms)))]);
        _out.WriteLine($"χ² over F's {f.Count} states: {chi2:F1} on {df}");
        Assert.True(chi2 <= df + 6 * Math.Sqrt(2.0 * df), $"χ² = {chi2:F1} on {df}");
    }

    [Fact]
    public void Draws_MatchThePriorOnF_ForADirectFunctionalFit()
    {
        // Two parts, each term's case side and output side counted apart for support (npart, npart.func).
        var model = BassModel.Create(new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 1, MaxFunctionalInteraction = 1, MaxBasis = 2 },
        }, FunctionalData());
        List<(BassTerm[] Terms, double LogPrior)> f = BassTermSpace.F(model);
        (Dictionary<string, double> law, double[] countLaw) = Law(f, model.MaxBasis);
        var sampler = new BassPriorSampler(model);
        const int draws = 150_000;
        var states = new Dictionary<string, long>();
        var counts = new long[countLaw.Length];
        var rng = new Xoshiro256PlusPlus(47);
        for (int d = 0; d < draws; d++)
        {
            BassState s = sampler.Draw(rng).State;
            states[Key(s.Terms)] = states.GetValueOrDefault(Key(s.Terms)) + 1;
            counts[s.Count]++;
        }
        _out.WriteLine($"{f.Count} states in F; acceptance {sampler.AcceptanceRate:F4}; refusals " +
                       string.Join(", ", sampler.Refusals.Entries.Where(e => e.Count > 0).Select(e => $"{e.Category} {e.Count}")));
        Assert.True(sampler.Refusals.Count(BassRefusal.Support) > 0, "hinges facing their end of the output range have no support");
        for (int m = 0; m < countLaw.Length; m++)
        {
            double observed = (double)counts[m] / draws, error = Math.Sqrt(countLaw[m] * (1 - countLaw[m]) / draws);
            Assert.True(Math.Abs(observed - countLaw[m]) <= 5 * error, $"M = {m}: {observed:F4} against {countLaw[m]:F4}");
        }
        Assert.True(states.Keys.All(law.ContainsKey));
        (double chi2, int df) = ChiSquare([.. f.Select(s => law[Key(s.Terms)] * draws)], [.. f.Select(s => states.GetValueOrDefault(Key(s.Terms)))]);
        _out.WriteLine($"χ² over F's {f.Count} states: {chi2:F1} on {df}");
        Assert.True(chi2 <= df + 6 * Math.Sqrt(2.0 * df), $"χ² = {chi2:F1} on {df}");
    }

    [Fact]
    public void Jeffreys_DrawsNoExactFit()
    {
        BassConfig config = BassSweepTests.Config() with { Sigma2 = Sigma2Prior.Jeffreys };
        var model = BassModel.Create(config, BassSweepTests.Data(exact: true));
        List<(BassTerm[] Terms, double LogPrior)> f = BassTermSpace.F(model);
        List<(BassTerm[] Terms, double LogPrior)> admitted = BassSweepTests.EnumerateF(model);
        Assert.True(admitted.Count > f.Count, "the fixture has exact fits for Jeffreys to exclude");
        (Dictionary<string, double> law, double[] countLaw) = Law(f, model.MaxBasis);

        var sampler = new BassPriorSampler(model);
        const int draws = 60_000;
        var counts = new long[countLaw.Length];
        var rng = new Xoshiro256PlusPlus(31);
        for (int d = 0; d < draws; d++)
        {
            (BassState s, BassEvaluation e) = sampler.DrawState(rng);
            Assert.Null(e.Refusal);
            Assert.True(law.ContainsKey(Key(s.Terms)), "an exact fit or a structure outside F");
            counts[s.Count]++;
        }
        _out.WriteLine($"{admitted.Count - f.Count} exact fits excluded; refused {sampler.Refusals.Count(BassRefusal.ExactFit)} in {sampler.Attempts} attempts");
        Assert.True(sampler.Refusals.Count(BassRefusal.ExactFit) > 0);
        for (int m = 0; m < countLaw.Length; m++)
        {
            double observed = (double)counts[m] / draws, error = Math.Sqrt(countLaw[m] * (1 - countLaw[m]) / draws);
            Assert.True(Math.Abs(observed - countLaw[m]) <= 5 * error, $"M = {m}: {observed:F4} against {countLaw[m]:F4}");
        }
    }

    [Fact]
    public void FixedG_HoldsTau()
    {
        var model = BassModel.Create(BassSweepTests.Config() with { Coefficients = new FixedG() }, BassSweepTests.Data());
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(37);
        for (int d = 0; d < 200; d++)
        {
            (BassState s, BassEvaluation e) = sampler.DrawState(rng);
            Assert.Equal(Math.Log(1.0 / 18), s.LogTau);
            Assert.Equal(1.0 / 18, e.Tau);
        }
    }

    [Fact]
    public void Draw_EvaluatesUnderTheModelsTarget()
    {
        var model = BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data());
        var sampler = new BassPriorSampler(model);
        Assert.Equal(TemperingPath.Likelihood, sampler.Path);
        Assert.True(double.IsNaN(sampler.AcceptanceRate));
        var rng = new Xoshiro256PlusPlus(41);
        for (int d = 0; d < 50; d++)
        {
            (BassState s, BassEvaluation e) = sampler.DrawState(rng);
            BassEvaluation fresh = model.Target.Evaluate(s);
            Assert.Equal(fresh.LogPrior, e.LogPrior);
            Assert.Equal(fresh.LogMarginal, e.LogMarginal);
            // The state's fit is the structure's own, as a fresh design loads it.
            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(s.Terms));
            Assert.Equal(design.ResidualSquares, s.Fit.ResidualSquares);
            Assert.Equal(s.Terms.Sum(model.Prior.LogPrior), s.Fit.LogTermPrior, 12);
        }
    }

    [Fact]
    public void Draw_RefusesPastItsAttemptLimit()
    {
        // No term reaches 18 nonzero points on the three-value grid, so F is the intercept alone; a count prior of mean
        // 10⁶ truncated at three terms leaves it about 6·10⁻⁹ of the proposal.
        BassConfig config = BassSweepTests.Config() with
        {
            Count = new BassCountPrior(1000, 1e-3),
            Structure = BassSweepTests.Config().Structure with { Support = new SupportRule(Minimum: 18) },
        };
        var model = BassModel.Create(config, BassSweepTests.Data());
        var sampler = new BassPriorSampler(model, maxAttempts: 500);
        var e = Assert.Throws<InvalidOperationException>(() => sampler.Draw(new Xoshiro256PlusPlus(43)));
        _out.WriteLine(e.Message);
        Assert.Contains("500 draws", e.Message);
        Assert.Contains("Support 500", e.Message);
        Assert.Equal(500, sampler.Attempts);
        Assert.Equal(0, sampler.Accepted);
        // A second refusal reports its own attempts, not the sampler's whole history.
        e = Assert.Throws<InvalidOperationException>(() => sampler.Draw(new Xoshiro256PlusPlus(44)));
        Assert.Contains("Support 500", e.Message);
        Assert.Equal(1000, sampler.Attempts);
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassPriorSampler(model, maxAttempts: 0));
    }

    // Each state's probability on F, by key, and the count law.
    private static (Dictionary<string, double> Law, double[] CountLaw) Law(List<(BassTerm[] Terms, double LogPrior)> f, int maxBasis)
    {
        double z = f.Max(s => s.LogPrior);
        double total = f.Sum(s => Math.Exp(s.LogPrior - z));
        var law = new Dictionary<string, double>();
        var countLaw = new double[maxBasis + 1];
        foreach ((BassTerm[] terms, double lp) in f)
        {
            double q = Math.Exp(lp - z) / total;
            law[Key(terms)] = q;
            countLaw[terms.Length] += q;
        }
        return (law, countLaw);
    }

    internal static string Key(IEnumerable<BassTerm> terms) => string.Join(" | ", terms.Select(t => t.ToString()).OrderBy(s => s, StringComparer.Ordinal));

    private static string OrderKey(IReadOnlyList<BassTerm> terms)
        => $"{terms.Count(t => t.Continuous.Count == 1)},{terms.Count(t => t.Continuous.Count == 2)}";

    private static string JointOrderKey(IEnumerable<BassTerm> terms)
        => string.Join(" ", terms.Select(t => $"({t.Continuous.Count},{t.Categorical.Count})").OrderBy(s => s, StringComparer.Ordinal));

    private static double Correlation(double[] a, double[] b)
    {
        double ma = a.Average(), mb = b.Average(), sab = 0, saa = 0, sbb = 0;
        for (int i = 0; i < a.Length; i++)
        {
            sab += (a[i] - ma) * (b[i] - mb);
            saa += (a[i] - ma) * (a[i] - ma);
            sbb += (b[i] - mb) * (b[i] - mb);
        }
        return sab / Math.Sqrt(saa * sbb);
    }

    // Pearson's χ², the bins with fewer than five expected pooled into one.
    private static (double Chi2, int Df) ChiSquare(IReadOnlyList<double> expected, IReadOnlyList<long> observed)
    {
        double chi2 = 0, pooledExpected = 0, pooledObserved = 0;
        int bins = 0;
        for (int i = 0; i < expected.Count; i++)
        {
            if (expected[i] < 5)
            {
                pooledExpected += expected[i];
                pooledObserved += observed[i];
                continue;
            }
            chi2 += (observed[i] - expected[i]) * (observed[i] - expected[i]) / expected[i];
            bins++;
        }
        if (pooledExpected > 0)
        {
            chi2 += (pooledObserved - pooledExpected) * (pooledObserved - pooledExpected) / pooledExpected;
            bins++;
        }
        return (chi2, bins - 1);
    }
}

/// <summary>
/// A BASS model's term space and its F, enumerated: every term of each part (each order the part allows, each variable
/// subset of that size, each sign and knot or each proper subset of levels; the intercept, every part empty, left out),
/// and every structure of at most M_max distinct terms the design admits and the target does not refuse, with its log
/// prior. F is downward closed (D06, and an exact fit's supersets fit exactly too), so a refused structure's supersets
/// are skipped.
/// </summary>
internal static class BassTermSpace
{
    public static List<BassTerm> Terms(BassModel model)
    {
        var terms = new List<BassTerm>();
        foreach (HingeFactor[]? continuous in Hinges(model, BassPart.Continuous))
            foreach (CategoricalFactor[]? categorical in Subsets(model))
                foreach (HingeFactor[]? functional in Hinges(model, BassPart.Functional))
                    if (continuous is not null || categorical is not null || functional is not null)
                        terms.Add(new BassTerm(continuous, categorical, functional));
        return terms;
    }

    public static List<(BassTerm[] Terms, double LogPrior)> F(BassModel model) => F(model, Terms(model));

    public static List<(BassTerm[] Terms, double LogPrior)> F(BassModel model, IReadOnlyList<BassTerm> terms)
    {
        BassStructurePrior prior = model.Prior;
        BassDesign design = model.NewDesign();
        var states = new List<(BassTerm[], double)>();
        void Visit(int from, List<BassTerm> chosen)
        {
            BassTerm[] structure = [.. chosen];
            double logTerms = structure.Sum(prior.LogPrior);
            if (design.Load(structure) is not null || model.Target.Refusal(BassStructureFit.Of(design, logTerms)) is not null) return;
            states.Add((structure, prior.LogCollectionWeight(structure.Length) + logTerms));
            if (chosen.Count == model.MaxBasis) return;
            for (int i = from; i < terms.Count; i++)
            {
                chosen.Add(terms[i]);
                Visit(i + 1, chosen);
                chosen.RemoveAt(chosen.Count - 1);
            }
        }
        Visit(0, []);
        return states;
    }

    // The part's options: null for the empty part where its order may be 0, then every factor set of each order.
    private static List<HingeFactor[]?> Hinges(BassModel model, BassPart part)
    {
        BassStructurePrior prior = model.Prior;
        if (!prior.Has(part)) return [null];
        var options = new List<HingeFactor[]?>();
        if (prior.MinimumOrder(part) == 0) options.Add(null);
        for (int k = Math.Max(1, prior.MinimumOrder(part)); k <= prior.Cap(part); k++)
            foreach (int[] set in Combinations(prior.Variables(part), k))
                foreach (HingeFactor[] factors in Product(set.Select(v => Factors(model, part, v)).ToList()))
                    options.Add(factors);
        return options;
    }

    private static List<HingeFactor> Factors(BassModel model, BassPart part, int v)
    {
        int knots = part == BassPart.Continuous ? model.Inputs.ContinuousKnots(v).Length : model.Inputs.FunctionalKnots(v).Length;
        return [.. from s in new[] { -1, 1 } from k in Enumerable.Range(0, knots) select new HingeFactor(v, s, k)];
    }

    private static List<CategoricalFactor[]?> Subsets(BassModel model)
    {
        BassStructurePrior prior = model.Prior;
        if (!prior.Has(BassPart.Categorical)) return [null];
        var options = new List<CategoricalFactor[]?>();
        if (prior.MinimumOrder(BassPart.Categorical) == 0) options.Add(null);
        for (int k = Math.Max(1, prior.MinimumOrder(BassPart.Categorical)); k <= prior.Cap(BassPart.Categorical); k++)
            foreach (int[] set in Combinations(prior.Variables(BassPart.Categorical), k))
                foreach (CategoricalFactor[] factors in Product(set.Select(v =>
                {
                    int n = model.Inputs.LevelCount(v);
                    return Enumerable.Range(1, (1 << n) - 2)
                        .Select(mask => new CategoricalFactor(v, Enumerable.Range(0, n).Where(l => (mask >> l & 1) != 0))).ToList();
                }).ToList()))
                    options.Add(factors);
        return options;
    }

    private static IEnumerable<int[]> Combinations(int n, int k)
    {
        if (k == 0) { yield return []; yield break; }
        for (int first = 0; first <= n - k; first++)
            foreach (int[] rest in Combinations(n - first - 1, k - 1))
                yield return [first, .. rest.Select(r => r + first + 1)];
    }

    private static IEnumerable<T[]> Product<T>(List<List<T>> choices)
    {
        IEnumerable<T[]> result = [[]];
        foreach (List<T> choice in choices) result = result.SelectMany(prefix => choice.Select(c => (T[])[.. prefix, c]));
        return result;
    }
}
