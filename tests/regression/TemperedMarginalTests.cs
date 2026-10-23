using System;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The observation model's marginal at an inverse temperature b (step 3.1), the likelihood half of the spline
/// target's evaluation. The Normal model's omits a constant that may depend on b and the sample size but not on the
/// design; the complete marginal from the collapsed g-prior block checks that, across designs of different sizes.
/// Models without a tempered closed form serve b = 1 only.
/// </summary>
public sealed class TemperedMarginalTests
{
    [Theory]
    [InlineData(0.3)]
    [InlineData(0.7)]
    [InlineData(1.0)]
    public void Normal_DiffersFromTheCompleteMarginal_ByADesignFreeConstant(double b)
    {
        (double[] x, double[] y) = Data(60);
        double tau = 1.0 / x.Length;
        foreach (GPriorCentre centre in new[] { GPriorCentre.Median, GPriorCentre.Mean, GPriorCentre.Zero })
        {
            var model = new WeightedNormalModel(centre);
            double? offset = null;
            foreach (double[] knots in new[] { Array.Empty<double>(), new[] { 0.5 }, new[] { 0.2, 0.45, 0.7 } })
            {
                double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
                double complete = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite, centre)!.At(tau, b).LogMarginal(null);
                double d = model.TemperedMarginal(z, y, null).LogMarginalLikelihood(b) - complete;
                offset ??= d;
                Assert.Equal(offset.Value, d, 1e-9);
            }
        }
    }

    [Fact]
    public void Normal_AtOne_IsTheUntemperedMarginal()
    {
        (double[] x, double[] y) = Data(40);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.4]), x);
        var model = new WeightedNormalModel();
        var w = new double[x.Length];
        for (int i = 0; i < w.Length; i++) w[i] = 0.5 + (i % 3) * 0.25;

        Assert.Equal(model.LogMarginalLikelihood(z, y, null), model.TemperedMarginal(z, y, null).LogMarginalLikelihood(1.0));
        Assert.Equal(model.LogMarginalLikelihood(z, y, w), model.TemperedMarginal(z, y, w).LogMarginalLikelihood(1.0));
    }

    [Fact]
    public void Normal_RefusesWhatItCannotTemper()
    {
        (double[] x, double[] y) = Data(40);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.4]), x);
        var model = new WeightedNormalModel();
        var w = new double[x.Length];
        Array.Fill(w, 2.0);

        Assert.Throws<NotSupportedException>(() => model.TemperedMarginal(z, y, w).LogMarginalLikelihood(0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => model.TemperedMarginal(z, y, null).LogMarginalLikelihood(0.0));
    }

    [Fact]
    public void Normal_OverKnotted_IsMinusInfinityAtEveryTemperature()
    {
        (double[] x, double[] y) = Data(6);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.3, 0.6]), x);   // 6 columns, 6 rows
        ITemperedMarginal marginal = new WeightedNormalModel().TemperedMarginal(z, y, null);
        Assert.Equal(double.NegativeInfinity, marginal.LogMarginalLikelihood(1.0));
        Assert.Equal(double.NegativeInfinity, marginal.LogMarginalLikelihood(0.4));
    }

    public static TheoryData<string> LaplaceModels => new() { "poisson", "whittle", "gibbs" };

    [Theory]
    [MemberData(nameof(LaplaceModels))]
    public void LaplaceModel_TempersItsMarginalLinearly(string name)
    {
        (double[] x, double[] y) = Data(50);
        IObservationModel model = Laplace(name);
        double[] response = Response(name, x, y);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.5]), x);

        double l1 = model.LogMarginalLikelihood(z, response, null);
        ITemperedMarginal marginal = model.TemperedMarginal(z, response, null);
        Assert.True(model.CanTemperLikelihood);
        Assert.Equal(l1, marginal.LogMarginalLikelihood(1.0));
        Assert.Equal(0.3 * l1, marginal.LogMarginalLikelihood(0.3));
        Assert.Equal(0.0, marginal.LogMarginalLikelihood(0.0));   // the prior
        Assert.Throws<ArgumentOutOfRangeException>(() => marginal.LogMarginalLikelihood(1.5));
    }

    [Theory]
    [MemberData(nameof(LaplaceModels))]
    public void LaplaceModel_OverKnotted_StaysOffTheSupport_EvenAtZero(string name)
    {
        (double[] x, double[] y) = Data(6);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.3, 0.6]), x);   // 6 columns, 6 rows
        ITemperedMarginal marginal = Laplace(name).TemperedMarginal(z, Response(name, x, y), null);
        Assert.Equal(double.NegativeInfinity, marginal.LogMarginalLikelihood(0.4));
        Assert.Equal(double.NegativeInfinity, marginal.LogMarginalLikelihood(0.0));   // not 0·(−∞) = NaN
    }

    [Fact]
    public void ModelWithoutATemperedMarginal_ServesOnlyTheUntemperedMarginal()
    {
        (double[] x, double[] y) = Data(40);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.5]), x);
        IObservationModel model = new ConstantModel();

        ITemperedMarginal marginal = model.TemperedMarginal(z, y, null);
        Assert.False(model.CanTemperLikelihood);
        Assert.Equal(-2.5, marginal.LogMarginalLikelihood(1.0));
        Assert.Throws<NotSupportedException>(() => marginal.LogMarginalLikelihood(0.5));
    }

    private static IObservationModel Laplace(string name) => name switch
    {
        "poisson" => new PoissonModel(),
        "whittle" => new WhittleModel(),
        _ => new GibbsLossModel(),
    };

    // Counts for Poisson, a positive periodogram for Whittle, the real response for Gibbs-loss.
    private static double[] Response(string name, double[] x, double[] y)
    {
        var r = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
            r[i] = name switch
            {
                "poisson" => Math.Round(Math.Exp(1.0 + y[i])),
                "whittle" => Math.Exp(y[i]) * (1.0 + 0.3 * Math.Cos(23.0 * x[i])),
                _ => y[i],
            };
        return r;
    }

    // A model that implements only the required members, so it takes the interface's defaults.
    private sealed class ConstantModel : IObservationModel
    {
        public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights) => -2.5;
        public int EffectiveDimension(double[,] design, double[]? weights) => design.GetLength(1);
        public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
            => new double[design.GetLength(1)];
    }

    [Fact]
    public void SplineTarget_OnTheLikelihoodPath_KeepsThePriorWhole()
    {
        (double[] x, double[] y) = Data(60);
        var target = new SplineTarget(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(3.0), x, y);
        ITargetEvaluation e = target.Evaluate(new KnotConfig([0.25, 0.6]));
        double b = 0.4;

        Assert.Equal(e.LogPrior + e.LogLikelihood(b), new Temperature(b, TemperingPath.Likelihood).LogDensity(e));
        Assert.Equal(b * (e.LogPrior + e.LogLikelihood(1.0)), new Temperature(b, TemperingPath.Powered).LogDensity(e));
        Assert.Equal(target.LogPosterior(new KnotConfig([0.25, 0.6])), e.LogPrior + e.LogLikelihood(1.0));
    }

    private static bool PositiveDefinite(ISpdFactor factor) => factor.Status.PositiveDefinite;

    private static (double[] X, double[] Y) Data(int n)
    {
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.1 * Math.Cos(17.0 * i);
        }
        return (x, y);
    }
}
