using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The hierarchical P-spline borrows strength across related curves: with many noisy replicates of a shared
/// population mean, the hierarchical per-replicate fits (shrunk toward the population) beat fitting each curve
/// independently, and the population curve is recovered. The classic shrinkage win.
/// </summary>
public sealed class HierarchicalPSplineTests
{
    private readonly ITestOutputHelper _out;
    public HierarchicalPSplineTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Mse(double[,] z, double[] beta, double[] truth)
    {
        int n = z.GetLength(0), nu = z.GetLength(1);
        double s = 0.0;
        for (int i = 0; i < n; i++)
        {
            double fit = 0.0;
            for (int j = 0; j < nu; j++) fit += z[i, j] * beta[j];
            double d = fit - truth[i];
            s += d * d;
        }
        return s / n;
    }

    [Fact]
    public void Hierarchical_BorrowsStrength_BeatsIndependent()
    {
        var rng = new Xoshiro256PlusPlus(seed: 2025);
        const int reps = 10, n = 50;
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = (i + 0.5) / n;

        var f0 = new double[n];
        for (int i = 0; i < n; i++) f0[i] = Math.Sin(2.0 * Math.PI * x[i]);

        // Each replicate = population + a small smooth deviation, observed with high noise.
        var truth = new double[reps][];
        var y = new double[reps][];
        for (int r = 0; r < reps; r++)
        {
            double amp = 0.12 * (2.0 * rng.NextDouble() - 1.0);
            double phase = 2.0 * Math.PI * rng.NextDouble();
            truth[r] = new double[n];
            y[r] = new double[n];
            for (int i = 0; i < n; i++)
            {
                truth[r][i] = f0[i] + amp * Math.Sin(Math.PI * x[i] + phase);
                y[r][i] = truth[r][i] + 0.30 * Gaussian(rng);   // high per-curve noise
            }
        }

        const int nKnots = 12;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        double[,] z = basis.Design(new KnotConfig(knots), x);
        var penalty = new DifferencePenalty(2);

        // Independent baseline: fit each curve alone, REML-optimal λ.
        double indepRep = 0.0;
        foreach (double[] yr in y)
        {
            var ps = new PenalizedSpline(z, yr, penalty);
            double bestLam = 1.0, bestEv = double.NegativeInfinity;
            for (int j = 0; j < 25; j++)
            {
                double lam = Math.Pow(10.0, -3.0 + 7.0 * j / 24.0);
                double ev = ps.RemlLogEvidence(lam);
                if (ev > bestEv) { bestEv = ev; bestLam = lam; }
            }
            indepRep += Mse(z, ps.Coefficients(bestLam), truth[Array.IndexOf(y, yr)]);
        }
        indepRep /= reps;

        // Hierarchical: shared population mean + shrunken replicate curves.
        HierarchicalResult h = new HierarchicalPSpline(z, y, penalty).Run(burn: 1000, samples: 2500, seed: 7);
        double hierRep = 0.0;
        for (int r = 0; r < reps; r++) hierRep += Mse(z, h.ReplicateCoefficients[r], truth[r]);
        hierRep /= reps;
        double popMse = Mse(z, h.PopulationCoefficients, f0);

        _out.WriteLine($"[hier] per-rep MSE: independent={indepRep:F5} hierarchical={hierRep:F5} " +
                       $"(ratio {hierRep / indepRep:F2}); population MSE={popMse:F5}; " +
                       $"σ̂={h.NoiseSd:F3} τ̂_u={h.DeviationSd:F3} τ̂_0={h.PopulationSmoothingSd:F3}");

        Assert.True(hierRep < indepRep, $"hierarchical per-rep MSE {hierRep:F5} should beat independent {indepRep:F5}");
        Assert.True(popMse < 0.02, $"population fit MSE {popMse:F5} should recover the shared mean");
    }
}
