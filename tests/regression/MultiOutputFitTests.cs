using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's fits (the multi-output brief, slice 1): <see cref="Bass.Fit"/> on several outputs, each retained draw holding Σ and one column of
/// coefficients per output, its run monitoring ln|Σ| and each monitored case's fitted value at every output.
/// <list type="bullet">
/// <item>At K = 1 a fit is scalar BASS's, draw for draw on one master seed: the structures and τ, σ² against Σ, the coefficients and λ, to
/// rounding.</item>
/// <item>At K = 3 every draw's Σ is positive definite and its coefficients are the structure's, the run's functionals are named by output,
/// and the readers of one response refuse the draws, which E34's readouts (slice 3) read.</item>
/// </list>
/// </summary>
public sealed class MultiOutputFitTests(ITestOutputHelper output)
{
    private static BassConfig Config(int samples) => new()
    {
        Structure = new BassStructureSpec { MaxBasis = 10 },
        Sampling = new BassSampling { Chains = 2, BurnIn = 300, Samples = samples, MasterSeed = 20261111 },
    };

    [Fact]
    public void AtOneOutput_TheFitIsScalarBass_DrawForDraw()
    {
        (List<BassPredictor> predictors, double[,] y) = MultiOutputCoreTests.Data(90, 1, 20261111);
        double[] column = [.. Enumerable.Range(0, 90).Select(i => y[i, 0])];
        BassFit scalar = Bass.Fit(Config(250), new BassData(predictors, new ScalarResponse(column)));
        BassFit several = Bass.Fit(Config(250), new BassData(predictors, new MultiOutputResponse(y, ["y"])));

        double worst = 0.0;
        int compared = 0;
        for (int c = 0; c < scalar.Draws.Count; c++)
        {
            Assert.Equal(scalar.Draws[c].Count, several.Draws[c].Count);
            for (int d = 0; d < scalar.Draws[c].Count; d++)
            {
                BassDraw one = scalar.Draws[c][d], other = several.Draws[c][d];
                Assert.Equal(TermsOf(scalar, one), TermsOf(several, other));
                Assert.Equal(one.Tau, other.Tau);
                Assert.True(double.IsNaN(other.Sigma2));
                Assert.Empty(other.Coefficients);
                worst = Math.Max(worst, Math.Abs(other.Sigma![0, 0] - one.Sigma2) / one.Sigma2);
                for (int j = 0; j < one.Coefficients.Length; j++)
                    worst = Math.Max(worst, Math.Abs(other.OutputCoefficients![j, 0] - one.Coefficients[j]) / Math.Max(1.0, Math.Abs(one.Coefficients[j])));
                worst = Math.Max(worst, Math.Abs(other.Lambda - one.Lambda) / one.Lambda);
                compared++;
            }
        }
        output.WriteLine($"{compared} draws, worst relative gap {worst:E2}");
        Assert.True(worst < 1e-11, $"worst relative gap {worst:E2}");
    }

    // A draw's terms, read through its fit's table, so two fits' tables may number them differently.
    private static BassTerm[] TermsOf(BassFit fit, BassDraw draw) => [.. draw.Structure.TermIds.Select(id => fit.Terms[id])];

    [Fact]
    public void AtThreeOutputs_EveryDrawHoldsSigmaAndEachOutputsCoefficients()
    {
        const int n = 90, k = 3;
        (List<BassPredictor> predictors, double[,] y) = MultiOutputCoreTests.Data(n, k, 20261112);
        var data = new BassData(predictors, new MultiOutputResponse(y, ["a", "b", "c"]));
        BassFit fit = Bass.Fit(Config(200), data);
        foreach (IReadOnlyList<BassDraw> chain in fit.Draws)
            foreach (BassDraw draw in chain)
            {
                double[,] sigma = draw.Sigma!;
                Assert.Equal(k, sigma.GetLength(0));
                Assert.True(PositiveDefinite(sigma));
                Assert.Equal(draw.Structure.TermIds.Length + 1, draw.OutputCoefficients!.GetLength(0));
                Assert.Equal(k, draw.OutputCoefficients.GetLength(1));
            }
        IReadOnlyList<string> names = fit.Diagnostics.Fitted.Select(v => v.Point).Select(p => $"{p.Case},{p.Location}").ToList();
        output.WriteLine($"monitored points: {string.Join("; ", names)}");
        Assert.Equal(3 * k, fit.Diagnostics.Fitted.Count);
        Assert.All(fit.Diagnostics.Fitted, v => Assert.InRange(v.Point.Location!.Value, 0, k - 1));
        Assert.Null(fit.ResidualDependence);
        Assert.Throws<NotSupportedException>(() => fit.PointwiseLogDensity(data));
    }

    private static bool PositiveDefinite(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = a[i, j];
                for (int r = 0; r < j; r++) s -= l[i, r] * l[j, r];
                if (i == j)
                {
                    if (!(s > 0.0)) return false;
                    l[i, i] = Math.Sqrt(s);
                }
                else l[i, j] = s / l[j, j];
            }
        return true;
    }
}
