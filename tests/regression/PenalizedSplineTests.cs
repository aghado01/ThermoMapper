using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The penalized P-spline mixed model selects its own smoothing: the REML log-evidence over λ has an interior
/// maximum (a finite penalty beats both the unpenalized overfit and the over-smoothed limit), and the penalized
/// fit at that λ recovers the truth better than either extreme — automatic smoothing, no tuning.
/// </summary>
public sealed class PenalizedSplineTests
{
    private readonly ITestOutputHelper _out;
    public PenalizedSplineTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void RemlEvidence_SelectsInteriorSmoothing_AndRecoversTruth()
    {
        var rng = new Xoshiro256PlusPlus(seed: 73);
        int n = 120;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(2.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.15 * Gaussian(rng);
        }

        // Over-rich fixed P-spline basis: many equally-spaced interior knots.
        const int nKnots = 25;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        double[,] z = basis.Design(new KnotConfig(knots), x);
        int nu = z.GetLength(1);

        var ps = new PenalizedSpline(z, y, new DifferencePenalty(2));

        // Geometric λ grid; find the REML-evidence argmax.
        const int g = 31;
        var lambdas = new double[g];
        var evidence = new double[g];
        int best = 0;
        for (int j = 0; j < g; j++)
        {
            lambdas[j] = Math.Pow(10.0, -3.0 + 8.0 * j / (g - 1));   // 1e-3 … 1e5
            evidence[j] = ps.RemlLogEvidence(lambdas[j]);
            if (evidence[j] > evidence[best]) best = j;
        }
        _out.WriteLine($"[reml] λ*={lambdas[best]:E2} ℓ*={evidence[best]:F3} " +
                       $"ℓ(min)={evidence[0]:F3} ℓ(max)={evidence[g - 1]:F3}");

        // Interior optimum: the best λ strictly beats both ends of the grid.
        Assert.True(best > 0 && best < g - 1, $"REML optimum landed at the grid edge (index {best}).");
        Assert.True(evidence[best] > evidence[0] && evidence[best] > evidence[g - 1],
            "REML evidence should peak in the interior, not at an extreme.");

        // The fit at λ* recovers the truth better than gross over/under-smoothing.
        double Mse(double lambda)
        {
            double[] beta = ps.Coefficients(lambda);
            double s = 0.0;
            for (int i = 0; i < n; i++)
            {
                double fit = 0.0;
                for (int p = 0; p < nu; p++) fit += z[i, p] * beta[p];
                double d = fit - f[i];
                s += d * d;
            }
            return s / n;
        }

        double mseStar = Mse(lambdas[best]);
        double mseUnder = Mse(lambdas[0]);      // λ→0, overfit
        double mseOver = Mse(lambdas[g - 1]);   // λ→∞, over-smoothed
        _out.WriteLine($"[reml] mse λ*={mseStar:F5} under={mseUnder:F5} over={mseOver:F5}");
        Assert.True(mseStar < mseUnder, $"REML fit {mseStar:F5} should beat the overfit {mseUnder:F5}");
        Assert.True(mseStar < mseOver, $"REML fit {mseStar:F5} should beat the over-smoothed {mseOver:F5}");
    }
}
