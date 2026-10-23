using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The Normal model's level centres, the median (default) and the mean.
/// <list type="bullet">
/// <item><b>Translation equivariance.</b> Shifting the response by a constant moves only the fit's level. The marginal,
/// the pointwise variance and the knots do not change, and the coefficients shift by the constant, because B-splines sum
/// to one. Under DMGK's zero centre the level's energy entered the scale, and a response far from zero against its noise
/// was under-fitted: at an offset of 10 noise SDs, the knot count fell by more than half.</item>
/// <item><b>Contamination.</b> An extreme response that robust weights discount cannot move the median centre. It moves
/// the mean centre by R/m, and the fit by about R/m².</item>
/// <item><b>Clean asymmetric responses.</b> On a narrow peak over a flat baseline, the median and the mean differ, and the
/// fits differ by exactly (median − mean)/(m + 1). That is far below the noise.</item>
/// </list>
/// </summary>
public sealed class NormalCentreTests
{
    private readonly ITestOutputHelper _out;
    public NormalCentreTests(ITestOutputHelper output) => _out = output;

    private static double[] Even(int n) => Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();

    // The level a centre puts the prior at, computed here independently of the block: the median, the mean, or zero.
    internal static double Level(double[] y, GPriorCentre centre)
    {
        if (centre == GPriorCentre.Zero) return 0.0;
        if (centre == GPriorCentre.Mean) return y.Average();
        double[] s = y.Order().ToArray();
        return s.Length % 2 == 1 ? s[s.Length / 2] : 0.5 * (s[s.Length / 2 - 1] + s[s.Length / 2]);
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void Marginal_Coefficients_AndVariance_AreLocationInvariant()
    {
        double[] x = Even(60);
        var rng = new Xoshiro256PlusPlus(seed: 61);
        double[] y = x.Select(t => Math.Sin(2.0 * Math.PI * t) + 0.2 * Gaussian(rng)).ToArray();
        double[] w = x.Select((_, i) => i % 7 == 3 ? 0.01 : 1.0 + 0.5 * (i % 3)).ToArray();
        double[,] readout = new SplineBasis(3).Design(new KnotConfig([0.3, 0.6]), [0.1, 0.45, 0.9]);
        foreach (GPriorCentre centre in new[] { GPriorCentre.Median, GPriorCentre.Mean })
        foreach (double[] knots in new[] { Array.Empty<double>(), new[] { 0.3, 0.6 }, new[] { 0.1, 0.25, 0.5, 0.8 } })
        {
            var model = new WeightedNormalModel(centre);
            double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
            double[,] r = new SplineBasis(3).Design(new KnotConfig(knots), [0.1, 0.45, 0.9]);
            foreach (double[]? weights in new[] { null, w })
            {
                double marginal = model.LogMarginalLikelihood(z, y, weights);
                double[] beta = model.PosteriorMeanCoefficients(z, y, weights);
                double[] variance = model.ConditionalVariance(z, y, weights, r)!;
                foreach (double offset in new[] { 7.0, -250.0, 1e4 })
                {
                    double[] shifted = y.Select(v => v + offset).ToArray();
                    // The centred response is computed from the shifted one, so about log₁₀|offset| digits go.
                    double tol = 1e-12 * Math.Max(1.0, Math.Abs(offset));
                    Assert.Equal(marginal, model.LogMarginalLikelihood(z, shifted, weights), 1e3 * tol);
                    double[] shiftedBeta = model.PosteriorMeanCoefficients(z, shifted, weights);
                    for (int j = 0; j < beta.Length; j++) Assert.Equal(beta[j] + offset, shiftedBeta[j], 10 * tol);
                    double[] shiftedVariance = model.ConditionalVariance(z, shifted, weights, r)!;
                    for (int t = 0; t < variance.Length; t++) Assert.Equal(variance[t], shiftedVariance[t], 1e-8 * variance[t] + tol);
                }
            }
        }
        Assert.Equal(3, readout.GetLength(0));
    }

    [Fact]
    public void ExtremeResponse_DoesNotMoveTheMedianCentre()
    {
        // The second external review's fixture: an intercept-only design, 39 zeros and one response R whose robust weight
        // is 1e-12. The median centre stays at 0, so the fit is the zeros' own, about 1e-12·R/m. The mean centre sits at R/m,
        // and the fit at R/m², because the prior still reads the outlier its likelihood weight ignores.
        const int m = 40;
        var z = new double[m, 1];
        for (int i = 0; i < m; i++) z[i, 0] = 1.0;
        var w = Enumerable.Repeat(1.0, m).ToArray();
        w[m - 1] = 1e-12;
        foreach (double outlier in new[] { 1e3, 1e4, 1e6 })
        {
            var y = new double[m];
            y[m - 1] = outlier;
            double median = new WeightedNormalModel(GPriorCentre.Median).PosteriorMeanCoefficients(z, y, w)[0];
            double mean = new WeightedNormalModel(GPriorCentre.Mean).PosteriorMeanCoefficients(z, y, w)[0];
            _out.WriteLine($"[contamination] R = {outlier:G3}: median centre {median:G3}, mean centre {mean:G3}");
            Assert.InRange(Math.Abs(median), 0.0, 1e-10 * outlier);
            Assert.Equal(outlier / (m * m), mean, 0.01 * outlier / (m * m));
        }
    }

    [Fact]
    public void RobustBars_IgnoresAGrossGlitch_UnderTheMedianCentre()
    {
        // A sine under noise 0.15 with one reading at 10⁶, fitted with Student-t weights. The weights discount the glitch
        // to about 1e-13. Under the median centre the fit is as good as without the glitch. Under the mean centre the
        // glitch moves the prior's level by 10⁶/m = 10⁴, and the fit by about 10⁶/m² = 100.
        double[] x = Even(100);
        var rng = new Xoshiro256PlusPlus(seed: 12);
        double[] clean = x.Select(t => Math.Sin(2.0 * Math.PI * t) + 0.15 * Gaussian(rng)).ToArray();
        double[] glitched = (double[])clean.Clone();
        glitched[50] += 1e6;
        double[] grid = Enumerable.Range(0, 26).Select(i => i / 25.0).ToArray();
        double Mse(GPriorCentre centre, double[] y)
        {
            BarsResult r = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(centre), new PoissonPrior(4.0), new LocalBetaKernel(50.0),
                                            resampler: new StudentTWeights(4.0))
                .Run(x, y, grid: grid, chains: 2, masterSeed: 5, burn: 500, samples: 1000);
            return grid.Select((t, g) => Math.Pow(r.Fit[g] - Math.Sin(2.0 * Math.PI * t), 2)).Average();
        }
        double reference = Mse(GPriorCentre.Median, clean);
        double median = Mse(GPriorCentre.Median, glitched), mean = Mse(GPriorCentre.Mean, glitched);
        _out.WriteLine($"[glitch] MSE: clean {reference:G3}; glitched, median centre {median:G3}, mean centre {mean:G3}");
        Assert.True(median < 2.0 * reference + 0.005, $"median centre MSE {median:G3} against the clean fit's {reference:G3}");
        Assert.True(mean > 100.0 * median, $"mean centre MSE {mean:G3}: the glitch should move it");
    }

    [Fact]
    public void CleanAsymmetricResponse_CostsTheMedianOnlyItsDistanceFromTheMean()
    {
        // A narrow peak over a flat baseline under noise 0.05 (n = 150): most responses sit near 0, so the median is about
        // 0.16 below the mean. For a basis that sums to one, the conditional mean fit is (m/(m+1))·f_OLS + c/(m+1), so the
        // two centres' fits differ by exactly (median − mean)/(m + 1) everywhere, about 0.001, a fiftieth of the noise.
        double[] x = Even(150);
        var rng = new Xoshiro256PlusPlus(seed: 11);
        double[] y = x.Select(t => Math.Exp(-120.0 * (t - 0.4) * (t - 0.4)) + 0.05 * Gaussian(rng)).ToArray();
        double median = Level(y, GPriorCentre.Median), mean = Level(y, GPriorCentre.Mean);
        Assert.True(mean - median > 0.1, $"the fixture should be asymmetric: median {median:G3}, mean {mean:G3}");
        double[] probe = Enumerable.Range(0, 21).Select(i => i / 20.0).ToArray();
        var basis = new SplineBasis(3);
        foreach (double[] knots in new[] { new[] { 0.4 }, new[] { 0.3, 0.38, 0.42, 0.5 }, new[] { 0.2, 0.35, 0.4, 0.45, 0.6, 0.8 } })
        {
            var config = new KnotConfig(knots);
            double[,] z = basis.Design(config, x);
            double[] bm = new WeightedNormalModel(GPriorCentre.Median).PosteriorMeanCoefficients(z, y, null);
            double[] ba = new WeightedNormalModel(GPriorCentre.Mean).PosteriorMeanCoefficients(z, y, null);
            foreach (double t in probe)
                Assert.Equal((median - mean) / (x.Length + 1), basis.Evaluate(config, bm, t) - basis.Evaluate(config, ba, t), 1e-12);
        }

        // Over the chains, the two centres fit the clean asymmetric data alike.
        double[] grid = Enumerable.Range(0, 51).Select(i => i / 50.0).ToArray();
        (double Knots, double Mse) Run(GPriorCentre centre)
        {
            BarsResult r = new BarsEnsemble(basis, new WeightedNormalModel(centre), new PoissonPrior(5.0), new LocalBetaKernel(50.0))
                .Run(x, y, grid: grid, chains: 4, masterSeed: 3, burn: 1500, samples: 3000);
            return (r.MeanKnots, grid.Select((t, g) => Math.Pow(r.Fit[g] - Math.Exp(-120.0 * (t - 0.4) * (t - 0.4)), 2)).Average());
        }
        (double km, double mm) = Run(GPriorCentre.Median);
        (double ka, double ma) = Run(GPriorCentre.Mean);
        _out.WriteLine($"[asymmetric] median {median:G4} mean {mean:G4}; knots {km:F2} / {ka:F2}; MSE {mm:G3} / {ma:G3}");
        Assert.InRange(km / ka, 0.95, 1.05);
        Assert.InRange(mm / ma, 0.8, 1.25);
    }

    [Fact]
    public void ConstantResponses_AreRefused_AsZeroCentredResponses()
    {
        // A constant response has nothing left once centred on the null fit, so under the improper σ² prior it has no
        // posterior and is refused as ZeroResponse. Before the exact handling, a projected constant left roundoff of
        // about 1e-16, admitted with marginals near 1,300 (the external review's 40-point cubic fixture).
        double[] x = Even(40);
        var model = new WeightedNormalModel();
        double[] w = x.Select((_, i) => i % 5 == 2 ? 0.01 : 1.0 + 0.5 * (i % 3)).ToArray();
        foreach (double[] knots in new[] { Array.Empty<double>(), new[] { 0.5 }, new[] { 0.2, 0.45, 0.7 } })
        {
            double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
            double[,] doubled = Scaled(z, 2.0);   // the same span, but its functions sum to two: the projection path
            foreach (double level in new[] { 0.0, 1.0, 7.0, 0.1, -250.0, 1e4 })
            {
                double[] y = Enumerable.Repeat(level, x.Length).ToArray();
                foreach (double[,] design in new[] { z, doubled })
                    foreach (double[]? weights in new[] { null, w })
                    {
                        Assert.Equal(DesignRejection.ZeroResponse, model.TemperedMarginal(design, y, weights).Rejection);
                        Assert.Equal(double.NegativeInfinity, model.LogMarginalLikelihood(design, y, weights));
                    }
                // The fit is the constant itself.
                double[] beta = model.PosteriorMeanCoefficients(z, y, null);
                Assert.All(beta, b => Assert.Equal(level, b, 1e-12 * Math.Max(1.0, Math.Abs(level))));
            }
            // DMGK's zero centre keeps a nonzero constant's level in the scale, so only the zero response is refused.
            var zero = new WeightedNormalModel(GPriorCentre.Zero);
            Assert.Equal(DesignRejection.None, zero.TemperedMarginal(z, Enumerable.Repeat(7.0, x.Length).ToArray(), null).Rejection);
            Assert.Equal(DesignRejection.ZeroResponse, zero.TemperedMarginal(z, new double[x.Length], null).Rejection);
        }
    }

    private static double[,] Scaled(double[,] z, double s)
    {
        var r = (double[,])z.Clone();
        for (int i = 0; i < r.GetLength(0); i++) for (int j = 0; j < r.GetLength(1); j++) r[i, j] *= s;
        return r;
    }

    [Fact]
    public void Bars_FitsAShiftedResponse_AsTheUnshiftedOne()
    {
        // The offset probe's comparison as a check: the same chains on the response and on the response plus 1,000.
        double[] x = Even(100);
        var rng = new Xoshiro256PlusPlus(seed: 12);
        double[] y = x.Select(t => Math.Sin(2.0 * Math.PI * t) + 0.15 * Gaussian(rng)).ToArray();
        double[] grid = Enumerable.Range(0, 26).Select(i => i / 25.0).ToArray();
        BarsResult Fit(double offset, GPriorCentre centre) =>
            new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(centre), new PoissonPrior(4.0), new LocalBetaKernel(50.0))
                .Run(x, y.Select(v => v + offset).ToArray(), grid: grid, chains: 2, masterSeed: 5, burn: 500, samples: 1000);

        BarsResult plain = Fit(0.0, GPriorCentre.Median), shifted = Fit(1000.0, GPriorCentre.Median);
        BarsResult zeroPlain = Fit(0.0, GPriorCentre.Zero), zeroShifted = Fit(1000.0, GPriorCentre.Zero);
        _out.WriteLine($"[centre] median: knots {plain.MeanKnots:F3} / {shifted.MeanKnots:F3}; zero: {zeroPlain.MeanKnots:F3} / {zeroShifted.MeanKnots:F3}");

        Assert.Equal(plain.MeanKnots, shifted.MeanKnots, 9);
        for (int g = 0; g < grid.Length; g++)
        {
            Assert.Equal(plain.Fit[g] + 1000.0, shifted.Fit[g], 1e-6);
            Assert.Equal(plain.FitSd[g], shifted.FitSd[g], 1e-6);
        }
        // DMGK's zero centre, the setting kept beside it: the shift costs it its knots.
        Assert.True(zeroShifted.MeanKnots < 0.5 * zeroPlain.MeanKnots,
            $"zero centre: {zeroShifted.MeanKnots:F3} knots shifted against {zeroPlain.MeanKnots:F3}");
    }
}
