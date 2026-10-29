using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Component execution and joint draws (stage 5, O2; D07): each component fitted to its coefficient response under the
/// seed tree's child seed, the same draws whatever the parallelism, per-component configs, and the pairing: by index
/// within chains, thinned evenly to the shortest count when the components' counts differ (fixed lengths and the R̂
/// stop), resampled within each component when weighted, and refused across chain counts.
/// </summary>
public sealed class BassOutputComponentsTests
{
    private readonly ITestOutputHelper _out;
    public BassOutputComponentsTests(ITestOutputHelper output) => _out = output;

    // Curves over t at m locations whose phase and amplitude move with z₀ and z₁ on a 6 × 4 grid (n = 24).
    internal static (double[,] Y, BassPredictor[] Predictors) Problem(int m = 15, int seed = 8)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int n = 24;
        var z0 = new double[n];
        var z1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            z0[i] = (i % 6) / 5.0;
            z1[i] = (i / 6) / 3.0;
            for (int l = 0; l < m; l++)
            {
                double t = m == 1 ? 0.5 : l / (m - 1.0);
                y[i, l] = Math.Sin(2 * Math.PI * (t - 0.3 * z0[i])) * (1 + z1[i]) + 2 * Math.Max(z0[i] - 0.5, 0) * t + 0.05 * Normal.SampleStandard(rng);
            }
        }
        return (y, [new ContinuousPredictor("z0", z0), new ContinuousPredictor("z1", z1)]);
    }

    internal static BassConfig Config(int chains = 2, int samples = 80, int seed = 31)
        => BassRunTests.Config(chains: chains, burn: 150, samples: samples, seed: seed, maxBasis: 5);

    private static readonly Dictionary<TruncationResidual, BassBasisFit> Fits = [];

    internal static BassBasisFit FitOf(TruncationResidual residual = TruncationResidual.LeaveOneCaseOut)
    {
        lock (Fits)
        {
            if (Fits.TryGetValue(residual, out BassBasisFit? cached)) return cached;
            (double[,] y, BassPredictor[] predictors) = Problem();
            BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.95, ScaleByColumnSd = true, Residual = residual });
            BassBasisFit fit = BassComponents.Fit(d, predictors, Config());
            Fits.Add(residual, fit);
            return fit;
        }
    }

    [Fact]
    public void EachComponent_IsFittedToItsCoefficientResponse_UnderTheSeedTree()
    {
        BassBasisFit fit = FitOf();
        int r = fit.Decomposition.Rank;
        _out.WriteLine($"r = {r}, draws per component: {string.Join(", ", fit.Components.Select(c => c.Draws.Sum(x => x.Count)))}");
        Assert.True(r >= 2);
        int[] seeds = SeedTree.Derive(fit.MasterSeed, r);
        for (int k = 0; k < r; k++)
        {
            BassFit component = fit.Components[k];
            Assert.Equal(seeds[k], component.Config.Sampling.MasterSeed);
            Assert.Equal(Bass.Fingerprint(component.Config, fit.ComponentData(k)), component.Identity.Fingerprint);
            Assert.Equal(fit.Decomposition.CoefficientResponse(k), ((ScalarResponse)fit.ComponentData(k).Response).Y.ToArray());
            // The same fit as running the component alone.
            BassFit alone = Bass.Fit(component.Config, fit.ComponentData(k));
            Assert.Equal(alone.Draws.SelectMany(c => c).Select(d => d.Sigma2), component.Draws.SelectMany(c => c).Select(d => d.Sigma2));
        }
    }

    [Fact]
    public void Draws_DoNotDependOnParallelism()
    {
        BassBasisFit fit = FitOf();
        (double[,] _, BassPredictor[] predictors) = Problem();
        foreach (int bound in new[] { 1, 3 })
        {
            BassBasisFit other = BassComponents.Fit(fit.Decomposition, predictors, Config(), new BassRunOptions { MaxDegreeOfParallelism = bound });
            for (int k = 0; k < fit.Decomposition.Rank; k++)
                Assert.Equal(fit.Components[k].Draws.SelectMany(c => c).Select(d => d.Sigma2), other.Components[k].Draws.SelectMany(c => c).Select(d => d.Sigma2));
        }
    }

    [Fact]
    public void EqualCounts_PairByIndexWithinChains()
    {
        BassJointDraws joint = FitOf().JointDraws;
        Assert.False(joint.Resampled);
        Assert.Equal([80, 80], joint.ChainLengths!);
        for (int k = 0; k < joint.Components; k++)
            for (int j = 0; j < joint.Count; j++) Assert.Equal(j, joint.Draw(k, j));
    }

    [Fact]
    public void UnequalCounts_AreThinnedEvenlyToTheShortest()
    {
        (double[,] y, BassPredictor[] predictors) = Problem();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 3 });
        BassConfig[] configs = [Config(samples: 80), Config(samples: 50), Config(samples: 64)];
        BassBasisFit fit = BassComponents.Fit(d, predictors, configs, masterSeed: 4);
        BassJointDraws joint = fit.JointDraws;
        Assert.Equal([50, 50], joint.ChainLengths!);
        int[] lengths = [80, 50, 64];
        for (int k = 0; k < 3; k++)
            for (int c = 0; c < 2; c++)
                for (int j = 0; j < 50; j++)
                {
                    int index = joint.Draw(k, c * 50 + j);
                    Assert.Equal(c * lengths[k] + (2 * j + 1) * lengths[k] / 100, index);
                }
        // The kept draws spread over the whole run: the first and last fifths of the longest chain both contribute.
        int[] kept = [.. Enumerable.Range(0, 50).Select(j => joint.Draw(0, j))];
        Assert.Contains(kept, i => i < 16);
        Assert.Contains(kept, i => i >= 64);
    }

    [Fact]
    public void RHatStoppedComponents_PairOverTheShortest()
    {
        (double[,] y, BassPredictor[] predictors) = Problem();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 3 });
        BassConfig config = Config(chains: 3, samples: 2000) is var c
            ? c with { Sampling = c.Sampling with { RHatTarget = 1.1, BatchSize = 50, MinEssPerChain = 10 } } : null!;
        BassBasisFit fit = BassComponents.Fit(d, predictors, config);
        int[] counts = [.. fit.Components.Select(f => f.Draws[0].Count)];
        _out.WriteLine($"draws per chain under the R̂ stop: {string.Join(", ", counts)}");
        BassJointDraws joint = fit.JointDraws;
        int shortest = counts.Min();
        Assert.Equal([shortest, shortest, shortest], joint.ChainLengths!);
        for (int k = 0; k < 3; k++)
        {
            Assert.All(fit.Components[k].Draws, chain => Assert.Equal(counts[k], chain.Count));
            for (int ch = 0; ch < 3; ch++)
                for (int j = 0; j < shortest; j++)
                {
                    int index = joint.Draw(k, ch * shortest + j);
                    Assert.Equal(ch * counts[k], index - index % counts[k]);   // from the same chain
                    Assert.Equal((int)((2L * j + 1) * counts[k] / (2L * shortest)), index - ch * counts[k]);
                }
        }
    }

    [Fact]
    public void WeightedComponents_AreResampledWithinEach()
    {
        BassBasisFit fit = FitOf();
        BassFit a = fit.Components[0], b = fit.Components[1];
        int pooled = a.Draws.Sum(c => c.Count);
        var weights = new double[pooled];
        weights[5] = 1.0;
        weights[40] = 3.0;
        BassFit weighted = a with { DrawWeights = weights };
        BassJointDraws joint = BassJointDraws.Pair([weighted, b], seed: 9);
        Assert.True(joint.Resampled);
        Assert.Null(joint.ChainLengths);
        Assert.Equal(pooled, joint.Count);
        int fives = Enumerable.Range(0, joint.Count).Count(j => joint.Draw(0, j) == 5);
        Assert.Equal(joint.Count, fives + Enumerable.Range(0, joint.Count).Count(j => joint.Draw(0, j) == 40));
        double expected = pooled / 4.0, sd = Math.Sqrt(pooled * 0.25 * 0.75);
        _out.WriteLine($"draw 5 taken {fives} times, expected {expected:F0} ± {sd:F1}");
        Assert.InRange(fives, expected - 4 * sd, expected + 4 * sd);
        // Shuffled, so the sorted ancestors do not line up with the other component's order.
        Assert.False(Enumerable.Range(1, joint.Count).All(j => j == joint.Count || joint.Draw(0, j - 1) <= joint.Draw(0, j)));
        // The equal-weight component beside it is thinned evenly over its pooled draws: here it keeps them all.
        for (int j = 0; j < joint.Count; j++) Assert.Equal(j, joint.Draw(1, j));
        // The stream is the seed's.
        BassJointDraws again = BassJointDraws.Pair([weighted, b], seed: 9);
        Assert.Equal(Enumerable.Range(0, joint.Count).Select(j => joint.Draw(0, j)), Enumerable.Range(0, again.Count).Select(j => again.Draw(0, j)));
    }

    [Fact]
    public void Refusals_SayWhy()
    {
        BassBasisFit fit = FitOf();
        (double[,] _, BassPredictor[] predictors) = Problem();
        BassFit a = fit.Components[0];
        BassFit oneChain = a with { Draws = [a.Draws[0]] };
        Assert.Contains("chain count", Assert.Throws<ArgumentException>(() => BassJointDraws.Pair([a, oneChain], 1)).Message);
        Assert.Contains("weights do not match", Assert.Throws<ArgumentException>(() => BassJointDraws.Pair([a with { DrawWeights = [1.0] }], 1)).Message);
        Assert.Contains("all zero", Assert.Throws<ArgumentException>(() =>
            BassJointDraws.Pair([a with { DrawWeights = new double[a.Draws.Sum(c => c.Count)] }], 1)).Message);
        Assert.Contains("one config per component", Assert.Throws<ArgumentException>(() =>
            BassComponents.Fit(fit.Decomposition, predictors, [Config()], 1)).Message);
        Assert.Contains("other data", Assert.Throws<ArgumentException>(() =>
            BassComponents.Assemble(fit.Decomposition, predictors, [.. fit.Components.Reverse()], 1)).Message);
        Assert.Contains("one fit per component", Assert.Throws<ArgumentException>(() =>
            BassComponents.Assemble(fit.Decomposition, predictors, [a], 1)).Message);
        Assert.Contains("values for the decomposition", Assert.Throws<ArgumentException>(() =>
            BassComponents.Assemble(fit.Decomposition, [new ContinuousPredictor("z0", [0.0, 1.0])], fit.Components, 1)).Message);
    }
}
