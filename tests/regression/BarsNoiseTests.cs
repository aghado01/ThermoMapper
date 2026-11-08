using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BARS under a variance-mean noise (OQ20, E1's slice 6):
/// <list type="bullet">
/// <item><b><see cref="StudentTWeights"/> through the shared law,</b> its weights the law's precision draws number for number and its
/// marginal the law's.</item>
/// <item><b>The noise's Normal model:</b> its priors read once from the data, and its marginal given σ² the normal density of the working
/// response with the coefficients integrated, against a dense computation.</item>
/// <item><b>The ensemble:</b> its refusals, the pointwise density under the law, and a quantile fit's curve.</item>
/// </list>
/// The law's joint-distribution test, recovery and cost are behind <c>BASS_QUALIFICATION</c> (<see cref="BarsNoiseQualificationTests"/>);
/// the resume is <see cref="BarsResumeTests"/>' noise setups.
/// </summary>
public sealed class BarsNoiseTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData(0.5)]
    [InlineData(3.0)]
    [InlineData(4.0)]
    [InlineData(30.0)]
    public void StudentTWeights_AreTheLawsPrecisionDraws_NumberForNumber(double nu)
    {
        var data = new Xoshiro256PlusPlus(11);
        double[] residuals = [.. Enumerable.Range(0, 500).Select(_ => 3.0 * Normal.SampleStandard(data))];
        const double sigma2 = 0.7;
        double[] weights = new StudentTWeights(nu).Resample(residuals, sigma2, new Xoshiro256PlusPlus(5));
        // The step as it was written before it went through the law: Gamma((ν + 1)/2, (ν + z²)/2).
        var rng = new Xoshiro256PlusPlus(5);
        for (int i = 0; i < residuals.Length; i++)
        {
            double z2 = residuals[i] * residuals[i] / sigma2;
            Assert.Equal(ConjugateUpdates.NormalPrecision(0.5 * nu, 0.5 * nu, 1.0, z2).Sample(rng), weights[i]);
        }
        var law = new StudentTMixing(nu);
        foreach (double r in residuals.Take(50))
            Assert.Equal(law.MarginalLogDensity(r, sigma2, 0.0), new StudentTWeights(nu).MarginalLogDensity(r, sigma2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new StudentTWeights(double.PositiveInfinity));
    }

    private static (double[] X, double[] Y) Curve(int n, int seed, Func<Xoshiro256PlusPlus, double, double> noise)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + noise(rng, x[i]);
        }
        return (x, y);
    }

    [Fact]
    public void TheNoiseModel_ReadsItsPriorsOnce_AndItsGivenMarginalIsTheWorkingResponsesDensity()
    {
        (double[] x, double[] y) = Curve(30, 3, (rng, _) => 0.3 * Normal.SampleStandard(rng));
        var model = new WeightedNormalModel(sigma2Prior: Sigma2Prior.InverseGamma(2.0, 0.5));
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.3, 0.55, 0.8]), x);
        var wrng = new Xoshiro256PlusPlus(4);
        double[] w = [.. y.Select(_ => 0.2 + 2.0 * wrng.NextDouble())];

        // On the data itself the reading is the model's own, bit for bit; on a shifted working response it keeps the data's level.
        WeightedNormalModel reading = model.ReadingFrom(y);
        Assert.Equal(model.TemperedMarginal(z, y, w).LogMarginalLikelihood(1.0), reading.TemperedMarginal(z, y, w).LogMarginalLikelihood(1.0));
        double[] shifted = [.. y.Select(v => v + 0.4)];
        Assert.NotEqual(model.TemperedMarginal(z, shifted, w).LogMarginalLikelihood(1.0),
                        reading.TemperedMarginal(z, shifted, w).LogMarginalLikelihood(1.0));

        // Given σ², y ~ N(Zμ₀, σ²(W⁻¹ + m·Z G⁻¹ Zᵀ)) with Zμ₀ the median level, the g-prior's covariance mσ²G⁻¹ at τ = 1/m.
        double level = BandedGPrior.Level(y, GPriorCentre.Median);
        foreach (double sigma2 in new[] { 0.05, 0.3, 2.0 })
        {
            double got = reading.Given(sigma2).TemperedMarginal(z, shifted, w).LogMarginalLikelihood(1.0);
            double expected = DenseLogDensity(z, shifted, w, sigma2, level);
            output.WriteLine($"σ² {sigma2}: given marginal {got:R} against the dense density {expected:R}");
            Assert.Equal(expected, got, 9);
        }
        Assert.Throws<InvalidOperationException>(() => model.Given(0.3));
    }

    // ln N(y | level·1, σ²(diag(1/w) + m·Z(ZᵀZ)⁻¹Zᵀ)), by a dense Cholesky.
    private static double DenseLogDensity(double[,] z, double[] y, double[] w, double sigma2, double level)
    {
        int m = z.GetLength(0), p = z.GetLength(1);
        var g = new double[p, p];
        for (int a = 0; a < p; a++)
            for (int b = 0; b < p; b++)
                for (int i = 0; i < m; i++) g[a, b] += z[i, a] * z[i, b];
        double[,] gInverse = Inverse(g);
        var cov = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int k = 0; k < m; k++)
            {
                double h = 0.0;
                for (int a = 0; a < p; a++)
                    for (int b = 0; b < p; b++) h += z[i, a] * gInverse[a, b] * z[k, b];
                cov[i, k] = sigma2 * (m * h + (i == k ? 1.0 / w[i] : 0.0));
            }
        double[,] l = Cholesky(cov);
        var r = new double[m];
        for (int i = 0; i < m; i++) r[i] = y[i] - level;
        var u = new double[m];   // L u = r
        double logDet = 0.0, quad = 0.0;
        for (int i = 0; i < m; i++)
        {
            double s = r[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * u[k];
            u[i] = s / l[i, i];
            logDet += 2.0 * Math.Log(l[i, i]);
            quad += u[i] * u[i];
        }
        return -0.5 * (m * Math.Log(2.0 * Math.PI) + logDet + quad);
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = a[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        return l;
    }

    private static double[,] Inverse(double[,] a)
    {
        int n = a.GetLength(0);
        double[,] l = Cholesky(a);
        var inverse = new double[n, n];
        for (int c = 0; c < n; c++)
        {
            var e = new double[n];
            e[c] = 1.0;
            var u = new double[n];
            for (int i = 0; i < n; i++)
            {
                double s = e[i];
                for (int k = 0; k < i; k++) s -= l[i, k] * u[k];
                u[i] = s / l[i, i];
            }
            for (int i = n - 1; i >= 0; i--)
            {
                double s = u[i];
                for (int k = i + 1; k < n; k++) s -= l[k, i] * inverse[k, c];
                inverse[i, c] = s / l[i, i];
            }
        }
        return inverse;
    }

    [Fact]
    public void TheEnsemble_RefusesWhatTheNoiseCannotRun()
    {
        BarsNoise noise = BarsNoise.Quantile(0.9);
        var basis = new SplineBasis(3);
        var prior = new PoissonPrior(5.0);
        var kernel = new LocalBetaKernel(50.0);
        Assert.Contains("Normal model", Assert.Throws<ArgumentException>(() => new BarsEnsemble(basis, new PoissonModel(), prior, kernel, noise: noise)).Message);
        Assert.Contains("resampler", Assert.Throws<ArgumentException>(() => new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel,
                                                                                            resampler: new StudentTWeights(), noise: noise)).Message);
        Assert.Contains("untempered", Assert.Throws<ArgumentException>(() => new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel,
                                                                                             temperLevels: 3, noise: noise)).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => BarsNoise.Quantile(1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BarsNoise.Skewed(new StudentTMixing(4.0), Normal.MeanSd(0.0, 0.0)));
        Assert.True(BarsNoise.Quantile(0.9).HasSkew);
        Assert.False(BarsNoise.StudentT(4.0).HasSkew);
        Assert.Equal(2.0 / 0.09, BarsNoise.Quantile(0.9).VarianceScale, 12);
        Assert.Equal(-0.8 / 0.09, BarsNoise.Quantile(0.9).Skew, 12);
    }

    [Fact]
    public void AQuantileFit_DrawsItsSkewedNoise_AndReadsTheLawPointwise()
    {
        // Exponential noise, skewed right: its 0.8-quantile is s·(−ln 0.2 − 1) above the sine.
        (double[] x, double[] y) = Curve(200, 7, (rng, _) => 0.3 * (-Math.Log(1.0 - rng.NextDouble()) - 1.0));
        const double q = 0.8;
        var ensemble = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                        noise: BarsNoise.Quantile(q));
        BarsResult result = ensemble.Run(x, y, grid: x, chains: 2, masterSeed: 3, burn: 800, samples: 300);
        double shift = 0.3 * (-Math.Log(1.0 - q) - 1.0), squares = 0.0;
        int below = 0;
        for (int i = 0; i < x.Length; i++)
        {
            double truth = Math.Sin(2.0 * Math.PI * x[i]) + shift;
            squares += (result.Fit[i] - truth) * (result.Fit[i] - truth);
            below += y[i] <= result.Fit[i] ? 1 : 0;
        }
        double rmse = Math.Sqrt(squares / x.Length), share = (double)below / x.Length;
        output.WriteLine($"AL({q}): the curve's RMSE against the true {q}-quantile {rmse:F4}; {share:P1} of the data below it");
        Assert.True(rmse < 0.2, $"RMSE {rmse:F4}");   // a smoke run; the recovery gate is the qualification's
        Assert.InRange(share, q - 0.08, q + 0.08);

        BarsDraw draw = result.Draws[0][0];
        Assert.Equal(BarsNoise.Quantile(q).Skew, draw.Skew);
        PointwiseLogDensity pointwise = ensemble.PointwiseLogDensity(result, x, y);
        double[,] zd = new SplineBasis(3).Design(draw.Knots, x);
        double eta = 0.0;
        for (int j = 0; j < draw.Coefficients.Length; j++) eta += zd[5, j] * draw.Coefficients[j];
        BarsNoise noise = BarsNoise.Quantile(q);
        Assert.Equal(noise.Mixing.MarginalLogDensity(y[5] - eta, draw.Sigma2 * noise.VarianceScale, Math.Sqrt(draw.Sigma2) * draw.Skew),
                     pointwise.Values[0, 5], 12);
        Assert.Contains(pointwise.Latents, l => l.Name == "latent scales vᵢ");
        Assert.Throws<NotSupportedException>(() => ensemble.PointwiseLogDensity(result, x, y, content: PointwiseContent.CollapsedPredictive));

        var gbp = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                   noise: BarsNoise.ScaleMixture(new GbpMixing(GbpParameters.Horseshoe)));
        BarsResult horseshoe = gbp.Run(x, y, grid: x, chains: 1, masterSeed: 3, burn: 100, samples: 20);
        Assert.Throws<NotSupportedException>(() => gbp.PointwiseLogDensity(horseshoe, x, y));
        Assert.Equal(20, gbp.PointwiseLogDensity(horseshoe, x, y, quadrature: true).Values.GetLength(0));
    }
}
