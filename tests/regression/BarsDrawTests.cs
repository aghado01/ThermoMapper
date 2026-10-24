using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The retained draw state (<see cref="BarsResult.Draws"/>): every draw keeps its knots, its coefficient and σ² draw,
/// and its robust weights, and the summaries are views over them. The fit is the mean of the conditional means
/// recomputed from each draw's knots and weights, the peak readout is read off its coefficients, the knot count's mean
/// is the draws' mean count, and σ² recovers the noise variance.
/// </summary>
public sealed class BarsDrawTests
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[] X, double[] Y) Data(int n, double noise, int seed, bool outliers = false)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-40.0 * (x[i] - 0.45) * (x[i] - 0.45)) + noise * Gaussian(rng) + (outliers && i % 13 == 0 ? 1.5 : 0.0);
        }
        return (x, y);
    }

    private static BarsEnsemble Normal(IWeightResampler? resampler = null) =>
        new(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0), resampler: resampler);

    private static double[] Grid(int points)
    {
        var grid = new double[points];
        for (int i = 0; i < points; i++) grid[i] = i / (points - 1.0);
        return grid;
    }

    [Fact]
    public void EveryDraw_KeepsKnotsCoefficientsAndSigma2()
    {
        (double[] x, double[] y) = Data(100, 0.1, seed: 3);
        var basis = new SplineBasis(3);
        BarsResult r = Normal().Run(x, y, Grid(21), chains: 2, masterSeed: 5, burn: 500, samples: 300);

        Assert.Equal(2, r.Draws.Count);
        double knots = 0.0;
        foreach (IReadOnlyList<BarsDraw> chain in r.Draws)
        {
            Assert.Equal(r.SamplesUsed, chain.Count);
            foreach (BarsDraw d in chain)
            {
                Assert.Equal(basis.Design(d.Knots, x).GetLength(1), d.Coefficients.Length);
                Assert.All(d.Coefficients, v => Assert.True(double.IsFinite(v)));
                Assert.True(d.Sigma2 > 0.0 && double.IsFinite(d.Sigma2));
                Assert.Null(d.Weights);
                knots += d.Knots.Count;
            }
        }
        // The knot count's pooled mean is the retained knots' mean count.
        Assert.Equal(r.MeanKnots, knots / (2 * r.SamplesUsed), 12);
    }

    [Fact]
    public void Fit_IsTheMeanOfTheDrawsConditionalMeans_UnderRobustWeights()
    {
        // With robust weights the conditional mean depends on the weights in force at the draw, so recomputing the fit
        // from each draw's knots and weights checks both were retained as the run used them.
        (double[] x, double[] y) = Data(80, 0.15, seed: 19, outliers: true);
        double[] grid = Grid(31);
        var basis = new SplineBasis(3);
        var model = new WeightedNormalModel();
        BarsResult r = Normal(new StudentTWeights(4.0)).Run(x, y, grid, chains: 2, masterSeed: 6, burn: 400, samples: 300);

        var fit = new double[grid.Length];
        int draws = 0;
        var distinctWeights = new HashSet<double[]>(ReferenceEqualityComparer.Instance);
        foreach (BarsDraw d in r.Draws.SelectMany(c => c))
        {
            Assert.NotNull(d.Weights);
            Assert.Equal(x.Length, d.Weights!.Length);
            distinctWeights.Add(d.Weights);
            double[] mean = model.PosteriorMeanCoefficients(basis.Design(d.Knots, x), y, d.Weights);
            double[,] zGrid = basis.Design(d.Knots, grid);
            for (int gi = 0; gi < grid.Length; gi++)
                for (int j = 0; j < mean.Length; j++) fit[gi] += zGrid[gi, j] * mean[j];
            draws++;
        }
        for (int gi = 0; gi < grid.Length; gi++) Assert.Equal(r.Fit[gi], fit[gi] / draws, 10);

        // The weight block runs every 25 steps, so the draws share far fewer arrays than there are draws.
        Assert.InRange(distinctWeights.Count, 2, draws / 10);
    }

    [Fact]
    public void PeakReadout_IsReadOffTheRetainedCoefficients()
    {
        (double[] x, double[] y) = Data(100, 0.1, seed: 3);
        var basis = new SplineBasis(3);
        BarsResult r = Normal().Run(x, y, Grid(21), chains: 2, masterSeed: 5, burn: 500, samples: 200);

        foreach (BarsDraw d in r.Draws.SelectMany(c => c))
        {
            (double location, double height) = SplineExtrema.Argmax(d.Knots, d.Coefficients, basis);
            Assert.Equal(location, d.Peak.ArgmaxLocation);
            Assert.Equal(height, d.Peak.ArgmaxHeight);
            Assert.Equal(SplineExtrema.SignificantPeakSpans(d.Knots, d.Coefficients, basis, 0.1, 0.5).ToArray(), d.Peak.Spans);
        }
    }

    [Fact]
    public void Sigma2_RecoversTheNoiseVariance()
    {
        const double noise = 0.1;
        (double[] x, double[] y) = Data(200, noise, seed: 41);
        BarsResult r = Normal().Run(x, y, Grid(21), chains: 4, masterSeed: 2, burn: 800, samples: 500);

        double[] sigma2 = r.Draws.SelectMany(c => c).Select(d => d.Sigma2).OrderBy(v => v).ToArray();
        double lo = sigma2[(int)(0.025 * sigma2.Length)], hi = sigma2[(int)(0.975 * sigma2.Length)];
        Assert.InRange(noise * noise, lo, hi);
        Assert.InRange(sigma2[sigma2.Length / 2], 0.75 * noise * noise, 1.3 * noise * noise);
    }

    [Fact]
    public void ModelWithoutNoiseScale_RetainsCoefficientsAndANaNSigma2()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        int n = 120;
        var x = new double[n];
        var counts = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double lambda = 3.0 + 5.0 * Math.Exp(-40.0 * (x[i] - 0.5) * (x[i] - 0.5)), p = Math.Exp(-lambda), s = p, u = rng.NextDouble();
            int k = 0;
            while (u > s && k < 200) { k++; p *= lambda / k; s += p; }
            counts[i] = k;
        }
        BarsResult r = new BarsEnsemble(new SplineBasis(3), new PoissonModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0))
            .Run(x, counts, Grid(21), chains: 2, masterSeed: 3, burn: 300, samples: 150);

        foreach (BarsDraw d in r.Draws.SelectMany(c => c))
        {
            Assert.True(double.IsNaN(d.Sigma2));
            Assert.All(d.Coefficients, v => Assert.True(double.IsFinite(v)));
        }
    }

    [Fact]
    public void StepBasis_RetainsCoefficients_WithoutAPeakReadout()
    {
        (double[] x, double[] y) = Data(80, 0.1, seed: 8);
        var basis = new StepBasis();
        BarsResult r = new BarsEnsemble(basis, new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0))
            .Run(x, y, Grid(21), chains: 2, masterSeed: 4, burn: 200, samples: 100);

        foreach (BarsDraw d in r.Draws.SelectMany(c => c))
        {
            Assert.Equal(basis.Design(d.Knots, x).GetLength(1), d.Coefficients.Length);
            Assert.True(d.Sigma2 > 0.0);
            Assert.True(double.IsNaN(d.Peak.ArgmaxLocation));
        }
    }
}
