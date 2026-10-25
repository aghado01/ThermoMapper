using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BAPS's pointwise density (OQ17, the retained-parameters content): each draw's term is the Normal likelihood given its
/// coefficients and σ², a curve unit's term is its observations' sum, every scheme's draws are scored the same way, and
/// WAIC and PSIS-LOO read the result.
/// </summary>
public sealed class BapsPointwiseTests
{
    private static (double[,] Z, double[] Y) Data(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x = Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
        double[] y = x.Select(v =>
        {
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            return Math.Sin(2.0 * Math.PI * v) + 0.2 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }).ToArray();
        var knots = new KnotConfig(Enumerable.Range(1, 12).Select(k => k / 13.0).ToArray());
        return (new SplineBasis(3).Design(knots, x), y);
    }

    [Theory]
    [InlineData(BapsScheme.Gibbs)]
    [InlineData(BapsScheme.Collapsed)]
    [InlineData(BapsScheme.LocallyAdaptiveGibbs)]
    public void Terms_AreTheNormalLikelihoodOfEachDraw(BapsScheme scheme)
    {
        (double[,] z, double[] y) = Data(60, seed: 4);
        var baps = new Baps(z, y, new DifferencePenalty(2));
        BapsResult r = baps.Run(scheme, chains: 2, masterSeed: 3, burn: 200, samples: 150);
        PointwiseLogDensity density = baps.PointwiseLogDensity(r);
        Assert.Equal(r.Draws, density.Draws);
        Assert.Equal(60, density.Units.Units);
        Assert.All(density.Latents, l => Assert.Equal(LatentTreatment.Conditioned, l.Treatment));

        int row = 0;
        foreach (var chain in r.PosteriorDraws)
            foreach (BapsDraw d in chain)
            {
                if (row % 37 == 0)
                    for (int i = 0; i < 60; i += 7)
                    {
                        double fit = 0.0;
                        for (int j = 0; j < d.Coefficients.Length; j++) fit += z[i, j] * d.Coefficients[j];
                        double expected = -0.5 * Math.Log(2.0 * Math.PI * d.Sigma2) - 0.5 * (y[i] - fit) * (y[i] - fit) / d.Sigma2;
                        Assert.Equal(expected, density.Values[row, i], 10);
                    }
                row++;
            }

        // Units of six consecutive observations sum their members' terms.
        PointwiseLogDensity blocks = baps.PointwiseLogDensity(r, UnitPartition.ConsecutiveBlocks(10, 6));
        for (int s = 0; s < density.Draws; s += 41)
            for (int u = 0; u < 10; u++)
            {
                double sum = 0.0;
                for (int k = 0; k < 6; k++) sum += density.Values[s, 6 * u + k];
                Assert.Equal(sum, blocks.Values[s, u], 9);
            }
    }

    [Fact]
    public void WaicAndPsisLoo_ReadTheTerms_AndAgreeOnASmoothFit()
    {
        (double[,] z, double[] y) = Data(80, seed: 5);
        var baps = new Baps(z, y, new DifferencePenalty(2));
        BapsResult r = baps.Run(BapsScheme.Gibbs, chains: 4, masterSeed: 2, burn: 500, samples: 1000);
        PointwiseLogDensity density = baps.PointwiseLogDensity(r);
        PredictiveCheck waic = PredictiveChecks.Waic(density), loo = PredictiveChecks.PsisLoo(density);
        Assert.True(double.IsFinite(waic.Elpd) && double.IsFinite(loo.Elpd));
        Assert.True(Math.Abs(waic.Elpd - loo.Elpd) < 1.0, $"WAIC {waic.Elpd:F2} against PSIS-LOO {loo.Elpd:F2}");
        // Retained-parameter terms of a flexible fit leave a few influential points above 0.7; nearly all are below it.
        Assert.True(loo.ParetoK!.Count(k => k > 0.7) <= 3, $"{loo.ParetoK!.Count(k => k > 0.7)} units with k̂ > 0.7");
        Assert.InRange(loo.EffectiveParameters, 3.0, 20.0);   // a smooth fit's effective degrees of freedom, plus σ²
    }

    [Fact]
    public void Partition_MustCoverTheFittedObservations()
    {
        (double[,] z, double[] y) = Data(40, seed: 6);
        var baps = new Baps(z, y, new DifferencePenalty(2));
        BapsResult r = baps.Run(BapsScheme.Gibbs, chains: 1, burn: 20, samples: 20);
        Assert.Throws<ArgumentException>(() => baps.PointwiseLogDensity(r, UnitPartition.EachObservation(39)));
    }
}
