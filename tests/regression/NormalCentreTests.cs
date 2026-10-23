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
/// The Normal model's null-fit centre: shifting the response by a constant moves only the fit's level. The marginal,
/// the pointwise variance and the knots do not change, and the coefficients shift by the constant, because B-splines
/// sum to one. Under DMGK's zero centre the level's energy entered the scale, and a response far from zero against its
/// noise was under-fitted: at an offset of 10 noise SDs, the knot count fell by more than half.
/// </summary>
public sealed class NormalCentreTests
{
    private readonly ITestOutputHelper _out;
    public NormalCentreTests(ITestOutputHelper output) => _out = output;

    private static double[] Even(int n) => Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();

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
        var model = new WeightedNormalModel();
        foreach (double[] knots in new[] { Array.Empty<double>(), new[] { 0.3, 0.6 }, new[] { 0.1, 0.25, 0.5, 0.8 } })
        {
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

        BarsResult plain = Fit(0.0, GPriorCentre.NullFit), shifted = Fit(1000.0, GPriorCentre.NullFit);
        BarsResult zeroPlain = Fit(0.0, GPriorCentre.Zero), zeroShifted = Fit(1000.0, GPriorCentre.Zero);
        _out.WriteLine($"[centre] null fit: knots {plain.MeanKnots:F3} / {shifted.MeanKnots:F3}; zero: {zeroPlain.MeanKnots:F3} / {zeroShifted.MeanKnots:F3}");

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
