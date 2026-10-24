using System;
using System.Collections.Generic;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The retained hierarchical P-spline draws (<see cref="HierarchicalResult.PosteriorDraws"/>): every draw keeps the
/// population and replicate curves, the three variances and the sweep's conditional means, and the reported summaries
/// are exactly what the draws give. Leaving the replicate draws out changes no summary.
/// </summary>
public sealed class HierarchicalDrawTests
{
    internal static (double[,] Z, double[][] Y) Data()
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
    public void PosteriorDraws_KeepEveryDraw_AndGiveTheSummaries()
    {
        (double[,] z, double[][] y) = Data();
        int nu = z.GetLength(1), reps = y.Length;
        HierarchicalResult h = new HierarchicalPSpline(z, y, new DifferencePenalty(2)).Run(burn: 200, samples: 400, seed: 7);

        Assert.Single(h.PosteriorDraws);
        IReadOnlyList<HierarchicalDraw> chain = h.PosteriorDraws[0];
        Assert.Equal(400, chain.Count);
        var pop = new double[nu];
        var rep = new double[reps][];
        for (int r = 0; r < reps; r++) rep[r] = new double[nu];
        double noise = 0.0, deviation = 0.0, smoothing = 0.0;
        foreach (HierarchicalDraw d in chain)
        {
            Assert.Equal(nu, d.Population.Length);
            Assert.Equal(nu, d.PopulationConditionalMean.Length);
            Assert.Equal(reps, d.Replicates!.Length);
            Assert.Equal(reps, d.ReplicateConditionalMeans!.Length);
            Assert.True(d.NoiseVariance > 0.0 && d.DeviationVariance > 0.0 && d.PopulationSmoothingVariance > 0.0);
            for (int j = 0; j < nu; j++) pop[j] += d.Population[j];
            for (int r = 0; r < reps; r++)
                for (int j = 0; j < nu; j++) rep[r][j] += d.Replicates[r][j];
            noise += Math.Sqrt(d.NoiseVariance);
            deviation += Math.Sqrt(d.DeviationVariance);
            smoothing += Math.Sqrt(d.PopulationSmoothingVariance);
        }
        // Distinct draws are distinct arrays, not reused buffers.
        Assert.NotSame(chain[0].Population, chain[1].Population);
        Assert.NotSame(chain[0].Replicates![0], chain[1].Replicates![0]);

        for (int j = 0; j < nu; j++) Assert.Equal(pop[j] / 400, h.PopulationCoefficients[j], 12);
        for (int r = 0; r < reps; r++)
            for (int j = 0; j < nu; j++) Assert.Equal(rep[r][j] / 400, h.ReplicateCoefficients[r][j], 12);
        Assert.Equal(noise / 400, h.NoiseSd, 12);
        Assert.Equal(deviation / 400, h.DeviationSd, 12);
        Assert.Equal(smoothing / 400, h.PopulationSmoothingSd, 12);
    }

    [Fact]
    public void WithoutReplicateDraws_TheSummariesAreUnchanged()
    {
        (double[,] z, double[][] y) = Data();
        var sampler = new HierarchicalPSpline(z, y, new DifferencePenalty(2));
        HierarchicalResult full = sampler.Run(burn: 200, samples: 400, seed: 7);
        HierarchicalResult lean = sampler.Run(burn: 200, samples: 400, seed: 7, retainReplicates: false);

        Assert.Equal(full.PopulationCoefficients, lean.PopulationCoefficients);
        for (int r = 0; r < y.Length; r++) Assert.Equal(full.ReplicateCoefficients[r], lean.ReplicateCoefficients[r]);
        Assert.Equal(full.NoiseSd, lean.NoiseSd);
        Assert.Equal(full.DeviationSd, lean.DeviationSd);
        Assert.Equal(full.PopulationSmoothingSd, lean.PopulationSmoothingSd);
        foreach (HierarchicalDraw d in lean.PosteriorDraws[0])
        {
            Assert.Null(d.Replicates);
            Assert.Null(d.ReplicateConditionalMeans);
        }
    }
}
