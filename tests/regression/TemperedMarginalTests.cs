using System;
using System.Linq;
using Maths.Distributions;
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
    public void Normal_LinearForm_IsTheUntemperedMarginalTimesB()
    {
        // The linear form is b·ℓ₁, defined at b = 0 under the improper 1/σ² prior; at b = 1 it is the integrated form.
        (double[] x, double[] y) = Data(60);
        var integrated = new WeightedNormalModel();
        var linear = new WeightedNormalModel(form: TemperedForm.Linear);
        Assert.False(integrated.CanTemperToPrior);
        Assert.True(linear.CanTemperToPrior);
        foreach (double[] knots in new[] { Array.Empty<double>(), new[] { 0.5 }, new[] { 0.2, 0.45, 0.7 } })
        {
            double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
            double atOne = integrated.TemperedMarginal(z, y, null).LogMarginalLikelihood(1.0);
            ITemperedMarginal m = linear.TemperedMarginal(z, y, null);
            Assert.Equal(atOne, m.LogMarginalLikelihood(1.0));
            Assert.Equal(0.0, m.LogMarginalLikelihood(0.0));
            Assert.Equal(0.3 * atOne, m.LogMarginalLikelihood(0.3));
            Assert.Equal(DesignRejection.None, m.Rejection);
        }
        double[,] overKnotted = new SplineBasis(3).Design(new KnotConfig([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]), x[..9]);
        ITemperedMarginal refused = linear.TemperedMarginal(overKnotted, y[..9], null);
        Assert.Equal(DesignRejection.TooFewObservations, refused.Rejection);
        Assert.Equal(double.NegativeInfinity, refused.LogMarginalLikelihood(0.0));
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

        var proper = new WeightedNormalModel(sigma2Prior: Sigma2Prior.ScaledToResponse());
        Assert.Throws<NotSupportedException>(() => proper.TemperedMarginal(z, y, w).LogMarginalLikelihood(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => proper.TemperedMarginal(z, y, null).LogMarginalLikelihood(-0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => proper.TemperedMarginal(z, y, null).LogMarginalLikelihood(1.5));
        Assert.True(double.IsFinite(proper.TemperedMarginal(z, y, w).LogMarginalLikelihood(1.0)));
    }

    public static TheoryData<string> ProperPriors => new() { "scaled", "unit", "fixed" };

    [Theory]
    [MemberData(nameof(ProperPriors))]
    public void Normal_ProperPrior_IsTheCompleteMarginal_DownToZero(string name)
    {
        // Under a proper σ² prior the model returns the block's complete marginal in either form, at every b in [0, 1],
        // exactly 0 at b = 0. The scaled prior's scale is the response's mean square about the centre level, computed here
        // from a sort of y, apart from the model's helper.
        (double[] x, double[] y) = Data(60);
        double tau = 1.0 / x.Length;
        foreach (GPriorCentre centre in new[] { GPriorCentre.Median, GPriorCentre.Mean, GPriorCentre.Zero })
        {
            double[] sorted = (double[])y.Clone();
            Array.Sort(sorted);
            double level = centre switch
            {
                GPriorCentre.Median => 0.5 * (sorted[29] + sorted[30]),
                GPriorCentre.Mean => y.Average(),
                _ => 0.0,
            };
            double meanSquare = y.Select(v => (v - level) * (v - level)).Average();
            (Sigma2Prior prior, InverseGamma ig) = name switch
            {
                "scaled" => (Sigma2Prior.ScaledToResponse(), InverseGamma.ShapeScale(1e-3, 1e-3 * meanSquare)),
                "unit" => (Sigma2Prior.ScaledToResponse(0.5), InverseGamma.ShapeScale(0.5, 0.5 * meanSquare)),
                _ => (Sigma2Prior.InverseGamma(2.0, 0.05), InverseGamma.ShapeScale(2.0, 0.05)),
            };
            var integrated = new WeightedNormalModel(centre, sigma2Prior: prior);
            var linear = new WeightedNormalModel(centre, TemperedForm.Linear, prior);
            Assert.True(integrated.CanTemperToPrior);
            foreach (double[] knots in new[] { Array.Empty<double>(), new[] { 0.5 }, new[] { 0.2, 0.45, 0.7 } })
            {
                double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
                CollapsedGPrior block = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite, centre)!;
                ITemperedMarginal m = integrated.TemperedMarginal(z, y, null);
                foreach (double b in new[] { 1.0, 0.7, 0.3, 1e-6 })
                    Assert.Equal(block.At(tau, b).LogMarginal(ig), m.LogMarginalLikelihood(b), 1e-9);
                Assert.Equal(0.0, m.LogMarginalLikelihood(0.0));
                double atOne = m.LogMarginalLikelihood(1.0);
                Assert.Equal(0.3 * atOne, linear.TemperedMarginal(z, y, null).LogMarginalLikelihood(0.3));
            }
        }
    }

    [Fact]
    public void Normal_ScaledPrior_IsScaleFree_WhereARawUnitPriorIsNot()
    {
        // Rescaling y by c shifts every configuration's marginal by the same −(mb/2)·ln c² under the scaled prior, as under
        // the Jeffreys prior, so the posterior over knots is unchanged across twelve decades. IG(10⁻³, 10⁻³) in raw units
        // is not (OQ19): at y × 10⁻³ its scale outweighs the residuals.
        (double[] x, double[] y) = Data(60);
        double[][] knots = [[], [0.5], [0.2, 0.45, 0.7]];
        double[] Differences(WeightedNormalModel model, double scale, double b)
        {
            double[] ys = y.Select(v => scale * v + 3.0 * scale).ToArray();
            double[] lm = knots.Select(k => model.TemperedMarginal(new SplineBasis(3).Design(new KnotConfig(k), x), ys, null)
                                                 .LogMarginalLikelihood(b)).ToArray();
            return [lm[1] - lm[0], lm[2] - lm[0]];
        }

        var scaled = new WeightedNormalModel(sigma2Prior: Sigma2Prior.ScaledToResponse());
        foreach (double b in new[] { 1.0, 0.4 })
        {
            double[] reference = Differences(scaled, 1.0, b);
            foreach (double scale in new[] { 1e-6, 1e-3, 1e3, 1e6 })
            {
                double[] d = Differences(scaled, scale, b);
                for (int j = 0; j < d.Length; j++) Assert.Equal(reference[j], d[j], 1e-7 * Math.Abs(reference[j]) + 1e-9);
            }
        }
        var raw = new WeightedNormalModel(sigma2Prior: Sigma2Prior.InverseGamma(1e-3, 1e-3));
        double[] rawAtOne = Differences(raw, 1.0, 1.0), rawSmall = Differences(raw, 1e-3, 1.0);
        Assert.True(Math.Abs(rawSmall[1] - rawAtOne[1]) > 1.0, $"the raw-unit prior moved the evidence by only {rawSmall[1] - rawAtOne[1]:G3}");
    }

    [Fact]
    public void Normal_VagueScaledPrior_KeepsTheJeffreysPosterior()
    {
        // IG(10⁻³, 10⁻³·s²) adds 10⁻³·s² to S/2 and 10⁻³ to the shape: about a·s²/σ̂² nats off each configuration, nearly
        // alike, so at this signal-to-noise the differences, all that the posterior over knots reads, move by hundredths.
        (double[] x, double[] y) = Data(60);
        var jeffreys = new WeightedNormalModel();
        var vague = new WeightedNormalModel(sigma2Prior: Sigma2Prior.ScaledToResponse());
        double[,] z0 = new SplineBasis(3).Design(new KnotConfig([]), x);
        foreach (double[] knots in new[] { new[] { 0.5 }, new[] { 0.2, 0.45, 0.7 }, new[] { 0.1, 0.3, 0.5, 0.7, 0.9 } })
        {
            double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
            double dj = jeffreys.LogMarginalLikelihood(z, y, null) - jeffreys.LogMarginalLikelihood(z0, y, null);
            double dv = vague.LogMarginalLikelihood(z, y, null) - vague.LogMarginalLikelihood(z0, y, null);
            Assert.True(Math.Abs(dv - dj) < 0.05, $"knots {knots.Length}: {dv:G6} vs {dj:G6}");
        }
    }

    [Fact]
    public void Normal_ConstantResponse_HasAMarginalOnlyUnderAFixedProperPrior()
    {
        // A constant response centres to zero (S = 0): the Jeffreys prior has no posterior for σ, the scaled prior no
        // scale, and a fixed proper prior a finite marginal, an exact fit it can weigh.
        (double[] x, _) = Data(40);
        var y = new double[x.Length];
        Array.Fill(y, 2.5);
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.4]), x);
        Assert.Equal(DesignRejection.ZeroResponse, new WeightedNormalModel().TemperedMarginal(z, y, null).Rejection);
        Assert.Equal(DesignRejection.ZeroResponse,
                     new WeightedNormalModel(sigma2Prior: Sigma2Prior.ScaledToResponse()).TemperedMarginal(z, y, null).Rejection);
        ITemperedMarginal fixedPrior = new WeightedNormalModel(sigma2Prior: Sigma2Prior.InverseGamma(2.0, 0.05)).TemperedMarginal(z, y, null);
        Assert.Equal(DesignRejection.None, fixedPrior.Rejection);
        Assert.True(double.IsFinite(fixedPrior.LogMarginalLikelihood(1.0)));
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
