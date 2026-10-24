using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// <see cref="HierarchicalPSpline"/> on the ensemble shell: chains start overdispersed, the shell reports rank-normalized
/// diagnostics on the log variances and a consensus map on the population curve, the pooled posterior agrees in
/// distribution with a long single chain from the original fixed start, a checkpointed run resumes exactly, another
/// problem's checkpoint is refused, and an R̂ target stops sampling early.
/// </summary>
public sealed class HierarchicalEnsembleTests
{
    private static (double[,] Z, double[][] Y) Data()
    {
        var rng = new Xoshiro256PlusPlus(seed: 41);
        int n = 40, reps = 6;
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = (i + 0.5) / n;
        var y = new double[reps][];
        for (int r = 0; r < reps; r++)
        {
            y[r] = new double[n];
            double shift = 0.3 * (rng.NextDouble() - 0.5);
            for (int i = 0; i < n; i++)
            {
                double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
                y[r][i] = Math.Sin(2.0 * Math.PI * x[i]) + shift + 0.2 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            }
        }
        var knots = new double[10];
        for (int k = 0; k < knots.Length; k++) knots[k] = (k + 1.0) / (knots.Length + 1);
        return (new SplineBasis(3).Design(new KnotConfig(knots), x), y);
    }

    [Fact]
    public void Chains_RunOnTheShell_WithDiagnosticsOnTheLogVariancesAndTheCurve()
    {
        (double[,] z, double[][] y) = Data();
        int nu = z.GetLength(1);
        HierarchicalResult h = new HierarchicalPSpline(z, y, new DifferencePenalty(2)).Run(chains: 4, masterSeed: 3, burn: 500, samples: 600);

        Assert.Equal(4, h.PosteriorDraws.Count);
        Assert.All(h.PosteriorDraws, chain => Assert.Equal(600, chain.Count));
        Assert.Equal(2400, h.Draws);
        Assert.Equal(4, h.ChainSeeds.Length);
        Assert.Null(h.Checkpoint);
        foreach (ConvergenceDiagnostic d in new[] { h.NoiseVarianceConvergence, h.PopulationSmoothingVarianceConvergence })
        {
            Assert.Equal(DiagnosticState.Estimated, d.RHat.State);
            Assert.True(d.RHat.Value < 1.05, $"R̂ {d.RHat}");
            Assert.True(d.BulkEss.Value > 100, $"bulk ESS {d.BulkEss}");
        }
        // τ_u² mixes slowly under the centred parameterization (about 2% of draws effective on this fit, against 15–20%
        // for σ² and τ_0²); the shell reports it, and the mixing study that follows is about it.
        Assert.Equal(DiagnosticState.Estimated, h.DeviationVarianceConvergence.RHat.State);
        Assert.True(h.DeviationVarianceConvergence.BulkEss.Value < h.NoiseVarianceConvergence.BulkEss.Value);
        Assert.Equal(nu, h.RHatPopulation.Length);
        Assert.Equal(RHatDiagnostic.Worst(h.RHatPopulation), h.RHatPopulationMax);
        Assert.True(h.RHatPopulationMax.Value < 1.05, $"population R̂ {h.RHatPopulationMax}");
    }

    [Fact]
    public void OverdispersedStarts_SpreadTheFirstDrawsFarWiderThanTheFixedStart()
    {
        // After one sweep, σ² still reflects where its chain started: a start with a tiny τ_u² holds the replicate curves
        // at the random starting population curve, and the residuals follow. From the fixed start every chain begins
        // alike, so the first σ² draws vary only by their own noise.
        (double[,] z, double[][] y) = Data();
        var sampler = new HierarchicalPSpline(z, y, new DifferencePenalty(2));
        double Spread(HierarchicalStart start)
        {
            double[] logSigma2 = sampler.Run(chains: 32, masterSeed: 5, burn: 0, samples: 1, start: start)
                                        .PosteriorDraws.Select(chain => Math.Log(chain[0].NoiseVariance)).ToArray();
            double mean = logSigma2.Average();
            return Math.Sqrt(logSigma2.Sum(v => (v - mean) * (v - mean)) / (logSigma2.Length - 1));
        }
        double overdispersed = Spread(HierarchicalStart.Overdispersed), fixedStart = Spread(HierarchicalStart.Fixed);
        Assert.True(overdispersed > 3.0 * fixedStart, $"sd of log σ² after one sweep: {overdispersed} against {fixedStart}");
    }

    // A posterior mean and its Monte Carlo standard error by batch means: each chain's draws in 20 consecutive batches,
    // the error from the spread of every batch's mean.
    private static (double Mean, double Se) BatchMean(HierarchicalResult h, Func<HierarchicalDraw, double> value)
    {
        var batches = new List<double>();
        foreach (IReadOnlyList<HierarchicalDraw> chain in h.PosteriorDraws)
        {
            int size = chain.Count / 20;
            for (int b = 0; b < 20; b++)
            {
                double s = 0.0;
                for (int i = b * size; i < (b + 1) * size; i++) s += value(chain[i]);
                batches.Add(s / size);
            }
        }
        double mean = batches.Average();
        double sd = Math.Sqrt(batches.Sum(v => (v - mean) * (v - mean)) / (batches.Count - 1));
        return (mean, sd / Math.Sqrt(batches.Count));
    }

    [Fact]
    public void PooledChains_AgreeInDistributionWithALongSingleChainFromTheFixedStart()
    {
        (double[,] z, double[][] y) = Data();
        var sampler = new HierarchicalPSpline(z, y, new DifferencePenalty(2));
        int nu = z.GetLength(1);
        HierarchicalResult pooled = sampler.Run(chains: 4, masterSeed: 11, burn: 1000, samples: 4000);
        HierarchicalResult single = sampler.Run(chains: 1, masterSeed: 12, burn: 1000, samples: 16000, start: HierarchicalStart.Fixed);

        var quantities = new List<(string Name, Func<HierarchicalDraw, double> Value)>
        {
            ("log σ²", d => Math.Log(d.NoiseVariance)),
            ("log τ_u²", d => Math.Log(d.DeviationVariance)),
            ("log τ_0²", d => Math.Log(d.PopulationSmoothingVariance)),
        };
        for (int j = 0; j < nu; j++)
        {
            int jj = j;
            quantities.Add(($"β_0,{j}", d => d.Population[jj]));
        }
        foreach ((string name, Func<HierarchicalDraw, double> value) in quantities)
        {
            (double m1, double s1) = BatchMean(pooled, value);
            (double m2, double s2) = BatchMean(single, value);
            double se = Math.Sqrt(s1 * s1 + s2 * s2);
            Assert.True(Math.Abs(m1 - m2) <= 4.5 * se, $"{name}: {m1} against {m2} (SE {se})");
        }
    }

    public static TheoryData<bool> RetainReplicates => new() { true, false };

    [Theory]
    [MemberData(nameof(RetainReplicates))]
    public void CheckpointedRun_ResumesExactly(bool retainReplicates)
    {
        (double[,] z, double[][] y) = Data();
        HierarchicalPSpline Sampler() => new(z, y, new DifferencePenalty(2));
        HierarchicalResult straight = Sampler().Run(chains: 2, masterSeed: 6, burn: 100, samples: 300, retainReplicates: retainReplicates);

        HierarchicalResult first = Sampler().Run(chains: 2, masterSeed: 6, burn: 100, samples: 120, retainReplicates: retainReplicates,
                                                 checkpoint: true);
        EnsembleCheckpoint checkpoint = first.Checkpoint!;
        Assert.Equal(120, checkpoint.DrawsPerChain);
        for (int attempt = 0; attempt < 2; attempt++)   // resuming leaves the checkpoint usable
        {
            HierarchicalResult resumed = Sampler().Run(chains: 2, masterSeed: 6, burn: 100, samples: 180,
                                                       retainReplicates: retainReplicates, resumeFrom: checkpoint);
            AssertSame(straight, resumed);
        }
    }

    public static TheoryData<HierarchicalScheme> Schemes => new()
    {
        HierarchicalScheme.Centred, HierarchicalScheme.CollapsedPopulation, HierarchicalScheme.Interweaved,
        HierarchicalScheme.CollapsedInterweaved, HierarchicalScheme.ParameterExpanded,
    };

    // The parameter-expanded sweep carries η, ξ and σ_η² besides the curves and variances, so this checks the whole state
    // travels with a checkpoint.
    [Theory]
    [MemberData(nameof(Schemes))]
    public void EveryScheme_ResumesExactly(HierarchicalScheme scheme)
    {
        (double[,] z, double[][] y) = Data();
        HierarchicalPSpline Sampler() => scheme == HierarchicalScheme.ParameterExpanded
            ? new HierarchicalPSpline(z, y, new DifferencePenalty(2), deviationPrior: new HalfTDeviationPrior())
            : new HierarchicalPSpline(z, y, new DifferencePenalty(2));
        HierarchicalResult straight = Sampler().Run(chains: 2, masterSeed: 8, burn: 100, samples: 200, scheme: scheme);
        EnsembleCheckpoint checkpoint = Sampler().Run(chains: 2, masterSeed: 8, burn: 100, samples: 90, scheme: scheme,
                                                      checkpoint: true).Checkpoint!;
        AssertSame(straight, Sampler().Run(chains: 2, masterSeed: 8, burn: 100, samples: 110, scheme: scheme, resumeFrom: checkpoint));
    }

    [Fact]
    public void Scheme_FollowsThePrior()
    {
        (double[,] z, double[][] y) = Data();
        var inverseGamma = new HierarchicalPSpline(z, y, new DifferencePenalty(2));
        var halfT = new HierarchicalPSpline(z, y, new DifferencePenalty(2), deviationPrior: new HalfTDeviationPrior(Degrees: 3.0, Scale: 0.2));

        // The parameter-expanded sweep belongs to the half-t prior, and the half-t prior samples with it alone.
        Assert.Throws<ArgumentException>(() => inverseGamma.Run(chains: 1, burn: 1, samples: 1, scheme: HierarchicalScheme.ParameterExpanded));
        foreach (HierarchicalScheme other in new[] { HierarchicalScheme.Centred, HierarchicalScheme.CollapsedPopulation, HierarchicalScheme.Interweaved,
                                                     HierarchicalScheme.CollapsedInterweaved })
            Assert.Throws<ArgumentException>(() => halfT.Run(chains: 1, burn: 1, samples: 1, scheme: other));

        // With no scheme each prior takes its own: the half-t run matches an explicit parameter-expanded one.
        HierarchicalResult byDefault = halfT.Run(chains: 1, masterSeed: 4, burn: 20, samples: 30);
        HierarchicalResult explicitly = halfT.Run(chains: 1, masterSeed: 4, burn: 20, samples: 30, scheme: HierarchicalScheme.ParameterExpanded);
        Assert.Equal(explicitly.PopulationCoefficients, byDefault.PopulationCoefficients);
        Assert.All(byDefault.PosteriorDraws[0], d => Assert.True(d.DeviationVariance > 0.0));

        Assert.Throws<ArgumentOutOfRangeException>(() => new HierarchicalPSpline(z, y, new DifferencePenalty(2),
                                                                                 deviationPrior: new HalfTDeviationPrior(Degrees: 0.0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => new HierarchicalPSpline(z, y, new DifferencePenalty(2),
                                                                                 deviationPrior: new HalfTDeviationPrior(Scale: -1.0)));
    }

    [Fact]
    public void Checkpoint_OfAnotherProblem_IsRefused()
    {
        (double[,] z, double[][] y) = Data();
        EnsembleCheckpoint checkpoint = new HierarchicalPSpline(z, y, new DifferencePenalty(2))
            .Run(chains: 2, masterSeed: 6, burn: 50, samples: 50, checkpoint: true).Checkpoint!;
        double[][] other = y.Select(r => (double[])r.Clone()).ToArray();
        other[1][3] += 1e-9;
        Assert.Throws<ArgumentException>(() => new HierarchicalPSpline(z, other, new DifferencePenalty(2))
            .Run(chains: 2, masterSeed: 6, burn: 50, samples: 50, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => new HierarchicalPSpline(z, y, new DifferencePenalty(2))
            .Run(chains: 2, masterSeed: 6, burn: 50, samples: 50, retainReplicates: false, resumeFrom: checkpoint));
    }

    [Fact]
    public void RHatTarget_StopsEarly_OnceTheLogVariancesAndTheCurveConverge()
    {
        (double[,] z, double[][] y) = Data();
        HierarchicalResult h = new HierarchicalPSpline(z, y, new DifferencePenalty(2))
            .Run(chains: 4, masterSeed: 2, burn: 500, samples: 20_000, rHatTarget: 1.01, batchSize: 500);

        int used = h.PosteriorDraws[0].Count;
        Assert.True(used < 20_000, $"no early stop ({used} draws per chain)");
        Assert.Equal(0, used % 500);
        foreach (ConvergenceDiagnostic d in new[] { h.NoiseVarianceConvergence, h.DeviationVarianceConvergence, h.PopulationSmoothingVarianceConvergence })
        {
            Assert.True(d.RHat.MeetsTarget(1.01), $"R̂ {d.RHat}");
            Assert.True(d.BulkEss.Value >= 4 * ChainEnsemble.MinEssPerChain && d.TailEss.Value >= 4 * ChainEnsemble.MinEssPerChain);
        }
        Assert.True(h.RHatPopulationMax.MeetsTarget(1.01), $"population R̂ {h.RHatPopulationMax}");
    }

    private static void AssertSame(HierarchicalResult a, HierarchicalResult b)
    {
        Assert.Equal(a.Draws, b.Draws);
        Assert.Equal(a.PopulationCoefficients, b.PopulationCoefficients);
        Assert.Equal(a.RaoBlackwellizedPopulationCoefficients, b.RaoBlackwellizedPopulationCoefficients);
        for (int r = 0; r < a.ReplicateCoefficients.Length; r++)
        {
            Assert.Equal(a.ReplicateCoefficients[r], b.ReplicateCoefficients[r]);
            Assert.Equal(a.RaoBlackwellizedReplicateCoefficients[r], b.RaoBlackwellizedReplicateCoefficients[r]);
        }
        Assert.Equal(a.NoiseSd, b.NoiseSd);
        Assert.Equal(a.DeviationSd, b.DeviationSd);
        Assert.Equal(a.PopulationSmoothingSd, b.PopulationSmoothingSd);
        Assert.Equal(a.NoiseVarianceConvergence, b.NoiseVarianceConvergence);
        Assert.Equal(a.DeviationVarianceConvergence, b.DeviationVarianceConvergence);
        Assert.Equal(a.PopulationSmoothingVarianceConvergence, b.PopulationSmoothingVarianceConvergence);
        Assert.Equal(a.RHatPopulation, b.RHatPopulation);
        Assert.Equal(a.ChainSeeds, b.ChainSeeds);
        for (int c = 0; c < a.PosteriorDraws.Count; c++)
            for (int i = 0; i < a.PosteriorDraws[c].Count; i++)
            {
                HierarchicalDraw u = a.PosteriorDraws[c][i], v = b.PosteriorDraws[c][i];
                Assert.Equal(u.Population, v.Population);
                Assert.Equal(u.PopulationConditionalMean, v.PopulationConditionalMean);
                Assert.Equal(u.NoiseVariance, v.NoiseVariance);
                Assert.Equal(u.DeviationVariance, v.DeviationVariance);
                Assert.Equal(u.PopulationSmoothingVariance, v.PopulationSmoothingVariance);
                Assert.Equal(u.Replicates is null, v.Replicates is null);
                if (u.Replicates is not null)
                    for (int r = 0; r < u.Replicates.Length; r++)
                    {
                        Assert.Equal(u.Replicates[r], v.Replicates![r]);
                        Assert.Equal(u.ReplicateConditionalMeans![r], v.ReplicateConditionalMeans![r]);
                    }
            }
    }
}
