using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The fit's views over its draws, wired from S9 into S7's <see cref="BassFit"/>: the residual-dependence check the run
/// fills at its end (F11), null for a scalar response, and the pointwise log density on the fit's own data (OQ17), the
/// same matrix as <see cref="BassPointwise.LogDensity"/> on a rebuilt model, with data of another fingerprint refused.
/// </summary>
public sealed class BassFitViewsTests
{
    private static BassConfig Config()
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6 },
            Sampling = new BassSampling { Chains = 2, MasterSeed = 3, BurnIn = 300, Samples = 60, Thin = 2 },
        };

    private static BassData Functional(int seed = 13, double shift = 0.0)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int n = 20, m = 6;
        double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i / 4) / 4.0;
            x1[i] = (i % 4) / 3.0;
            for (int r = 0; r < m; r++)
                y[i, r] = shift + 1 + 2 * x0[i] + Math.Sin(3 * t[r]) * x1[i] + 0.2 * (rng.NextDouble() - 0.5);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    private static BassData Scalar()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double[] x0 = [.. Enumerable.Range(0, 30).Select(i => (i % 5) / 4.0)];
        double[] x1 = [.. Enumerable.Range(0, 30).Select(i => (i / 5 % 3) / 2.0)];
        double[] y = [.. x0.Zip(x1, (a, b) => 1 + 2 * Math.Max(a - 0.5, 0) + b + 0.2 * (rng.NextDouble() - 0.5))];
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    [Fact]
    public void ResidualDependence_IsFilledAtTheEndOfAFunctionalRun()
    {
        BassData data = Functional();
        BassFit fit = Bass.Fit(Config(), data);
        ResidualDependence check = Assert.IsType<ResidualDependence>(fit.ResidualDependence);
        ResidualDependence expected = BassResidualDependence.Of(BassModel.Create(fit.Config, data), fit.Terms, fit.Draws)!;
        Assert.Equal(expected.PerCase, check.PerCase);
        Assert.Equal(expected.Quantiles, check.Quantiles);
        Assert.Equal(20, check.PerCase.Length);

        Assert.Null(Bass.Fit(Config(), Scalar()).ResidualDependence);
    }

    [Theory]
    [InlineData(PointwiseContent.RetainedParameters)]
    [InlineData(PointwiseContent.CollapsedPredictive)]
    public void PointwiseLogDensity_ReadsTheFitsOwnData(PointwiseContent content)
    {
        BassData data = Functional();
        BassFit fit = Bass.Fit(Config(), data);
        PointwiseLogDensity view = fit.PointwiseLogDensity(data, content: content);
        PointwiseLogDensity direct = BassPointwise.LogDensity(BassModel.Create(fit.Config, data), fit.Terms, fit.Draws, content: content);
        Assert.Equal(direct.Values, view.Values);
        Assert.Equal(direct.Units.UnitOf, view.Units.UnitOf);
        Assert.Equal(fit.Draws.Sum(d => d.Count), view.Values.GetLength(0));
        Assert.Equal(20, view.Values.GetLength(1));   // a direct fit's default unit is a case

        var ex = Assert.Throws<ArgumentException>(() => fit.PointwiseLogDensity(Functional(shift: 1.0), content: content));
        Assert.Contains("other data", ex.Message);
    }
}
