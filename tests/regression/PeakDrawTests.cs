using System;
using System.Collections.Generic;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The retained per-draw peak readout (<see cref="BarsResult.PeakDraws"/>): every draw of every chain is kept with its
/// exact closed-form argmax and spans, off the grid, with both clipping flags; λ(T) and π(T) are exactly its grid
/// counts over the draw count.
/// </summary>
public sealed class PeakDrawTests
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[] X, double[] Y) Data(Func<double, double> f, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int n = 120;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = f(x[i]) + 0.05 * Gaussian(rng);
        }
        return (x, y);
    }

    private static BarsResult Run(double[] x, double[] y, double[] grid) =>
        new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0))
            .Run(x, y, grid, chains: 2, masterSeed: 4, burn: 800, samples: 400);

    private static double[] Grid(int points)
    {
        var grid = new double[points];
        for (int i = 0; i < points; i++) grid[i] = i / (points - 1.0);
        return grid;
    }

    private static double DistanceToGrid(double[] grid, double t)
    {
        double best = double.PositiveInfinity;
        foreach (double g in grid) best = Math.Min(best, Math.Abs(g - t));
        return best;
    }

    [Fact]
    public void PeakDraws_KeepEveryDraw_OffTheGrid()
    {
        (double[] x, double[] y) = Data(t => Math.Exp(-80.0 * (t - 0.37) * (t - 0.37)), seed: 21);
        double[] grid = Grid(11);
        BarsResult r = Run(x, y, grid);

        Assert.Equal(2, r.PeakDraws.Count);
        double farthestSpan = 0.0, farthestArgmax = 0.0;
        foreach (IReadOnlyList<PeakDraw> chain in r.PeakDraws)
        {
            Assert.Equal(r.SamplesUsed, chain.Count);
            foreach (PeakDraw d in chain)
            {
                farthestArgmax = Math.Max(farthestArgmax, DistanceToGrid(grid, d.ArgmaxLocation));
                for (int p = 0; p < d.Spans.Count; p++)
                {
                    PeakSpan s = d.Spans[p];
                    Assert.InRange(s.Location, s.Left, s.Right);
                    Assert.True(double.IsFinite(s.Height));
                    if (p > 0) Assert.True(s.Location > d.Spans[p - 1].Location);   // location order
                    farthestSpan = Math.Max(farthestSpan, DistanceToGrid(grid, s.Location));
                }
            }
        }
        // On a 0.1-spaced grid, a grid-snapped location would sit on a grid point; the retained ones do not.
        Assert.True(farthestSpan > 1e-3, $"span locations lie on the grid (max distance {farthestSpan})");
        Assert.True(farthestArgmax > 1e-3, $"argmax locations lie on the grid (max distance {farthestArgmax})");
    }

    [Fact]
    public void LambdaAndPi_AreTheDrawsGridCounts()
    {
        (double[] x, double[] y) = Data(t => Math.Exp(-150.0 * (t - 0.3) * (t - 0.3)) + Math.Exp(-150.0 * (t - 0.7) * (t - 0.7)), seed: 11);
        double[] grid = Grid(41);
        BarsResult r = Run(x, y, grid);

        var lambda = new double[grid.Length];
        var pi = new double[grid.Length];
        double draws = 0.0, spans = 0.0;
        foreach (IReadOnlyList<PeakDraw> chain in r.PeakDraws)
            foreach (PeakDraw d in chain)
            {
                draws++;
                spans += d.Spans.Count;
                foreach (PeakSpan s in d.Spans) lambda[Nearest(grid, s.Location)]++;
                for (int gi = 0; gi < grid.Length; gi++)
                    foreach (PeakSpan s in d.Spans)
                        if (grid[gi] >= s.Left && grid[gi] <= s.Right) { pi[gi]++; break; }
            }
        for (int gi = 0; gi < grid.Length; gi++)
        {
            lambda[gi] /= draws;
            pi[gi] /= draws;
        }
        Assert.Equal(r.PeakIntensity, lambda);
        Assert.Equal(r.SpanCoverage, pi);
        Assert.Equal(r.PeakCountMean, spans / draws, 12);
    }

    [Fact]
    public void ClippingFlags_SurviveInTheDraws()
    {
        // A rising ramp peaks at the right edge: the curve never descends to the span level before 1, so the span clips
        // there and not on the left. λ and π cannot show this; the retained spans must.
        (double[] x, double[] y) = Data(t => t, seed: 5);
        BarsResult r = Run(x, y, Grid(21));

        int draws = 0, clippedRight = 0;
        foreach (IReadOnlyList<PeakDraw> chain in r.PeakDraws)
            foreach (PeakDraw d in chain)
            {
                draws++;
                foreach (PeakSpan s in d.Spans)
                    if (s.RightClipped)
                    {
                        clippedRight++;
                        Assert.False(s.LeftClipped);
                        Assert.Equal(s.Location, s.Right);   // the span ends at the boundary apex
                        break;
                    }
            }
        Assert.True(clippedRight > draws / 2, $"{clippedRight} of {draws} draws kept a right-clipped span");
    }

    // The grid point nearest t, ties to the lower index (BarsEnsemble's deposit rule).
    private static int Nearest(double[] grid, double t)
    {
        int best = 0;
        double bd = double.PositiveInfinity;
        for (int i = 0; i < grid.Length; i++) { double d = Math.Abs(grid[i] - t); if (d < bd) { bd = d; best = i; } }
        return best;
    }
}
