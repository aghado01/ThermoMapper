using System;
using Maths.LinAlg;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// R08 / Q4, and Q6's numerical gate (step 3.7): every BARS observation model admits a design only when its strict
/// factor is positive definite and its smallest equilibrated pivot is at least 1e-10; anything else scores −∞. The gate
/// does not depend on column scaling. The IRLS models also reject a fit that did not converge.
/// </summary>
public sealed class DesignAdmissionTests
{
    public static TheoryData<string> Models => new() { "normal", "poisson", "whittle", "gibbs" };

    [Theory]
    [MemberData(nameof(Models))]
    public void DesignWithUnsupportedBasisFunctions_IsRejected(string model)
    {
        // Data on [0, 0.4] ∪ [0.6, 1]; knots packed into the gap leave cubic B-splines with no data under them.
        double[] x = GappedInputs();
        var knots = new[] { 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58 };
        double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);

        Assert.Equal(double.NegativeInfinity, Model(model).LogMarginalLikelihood(z, Response(model, x), null));
    }

    [Theory]
    [MemberData(nameof(Models))]
    public void NearlyCollinearDesign_IsRejected(string model)
    {
        // Columns 1, x and x + 1e-6·cos 7x: the last pivot is 2.0e-11 against its diagonal 13.3, so the relative pivot
        // is 1.5e-12, below the admission threshold. (At 1e-5 it is 1.5e-10, just inside it.)
        double[] x = EvenInputs(40);
        var z = new double[x.Length, 3];
        for (int i = 0; i < x.Length; i++)
        {
            z[i, 0] = 1.0;
            z[i, 1] = x[i];
            z[i, 2] = x[i] + 1e-6 * Math.Cos(7.0 * x[i]);
        }

        Assert.Equal(double.NegativeInfinity, Model(model).LogMarginalLikelihood(z, Response(model, x), null));
    }

    [Theory]
    [MemberData(nameof(Models))]
    public void WellPosedDesign_IsAdmitted(string model)
    {
        double[] x = EvenInputs(40);
        double[,] z = new SplineBasis(3).Design(new KnotConfig(new[] { 0.3, 0.6 }), x);

        Assert.True(double.IsFinite(Model(model).LogMarginalLikelihood(z, Response(model, x), null)));
    }

    [Theory]
    [MemberData(nameof(Models))]
    public void ColumnScaling_DoesNotChangeAdmission(string model)
    {
        // One basis column scaled by 1e-8: the pivot ratio fell to about 1e-16 and refused this design; the
        // equilibrated pivots do not move. The Normal model's g-prior marginal is invariant to column scaling as well.
        double[] x = EvenInputs(40);
        double[,] z = new SplineBasis(3).Design(new KnotConfig(new[] { 0.3, 0.6 }), x);
        double[,] scaled = (double[,])z.Clone();
        for (int i = 0; i < x.Length; i++) scaled[i, 2] *= 1e-8;
        double[] y = Response(model, x);

        double plain = Model(model).LogMarginalLikelihood(z, y, null);
        double rescaled = Model(model).LogMarginalLikelihood(scaled, y, null);
        Assert.True(double.IsFinite(rescaled), $"{model}: the rescaled design was refused");
        if (model == "normal") Assert.Equal(plain, rescaled, 6);
    }

    public static TheoryData<string> IrlsModels => new() { "poisson", "whittle", "gibbs" };

    [Theory]
    [MemberData(nameof(IrlsModels))]
    public void IrlsThatDoesNotConverge_IsRejected(string model)
    {
        // A one-iteration cap cannot settle the step from β = 0 (GibbsLoss needs at least a second, reweighted pass).
        double[] x = EvenInputs(40);
        double[,] z = new SplineBasis(3).Design(new KnotConfig(new[] { 0.3, 0.6 }), x);
        IObservationModel capped = model switch
        {
            "poisson" => new PoissonModel(maxIterations: 1),
            "whittle" => new WhittleModel(maxIterations: 1),
            _ => new GibbsLossModel(4.685, maxIterations: 1),
        };

        Assert.Equal(double.NegativeInfinity, capped.LogMarginalLikelihood(z, Response(model, x), null));
    }

    [Fact]
    public void HealthyPivots_DoNotAdmitAnIllConditionedFactor()
    {
        // A = L·Lᵀ with L lower bidiagonal, diagonal 1 and subdiagonal 2 (the external review's counterexample): every
        // relative pivot is at least 0.2, so the pivot check passes, but κ(E) is about 2·10¹⁵ and the estimate refuses it.
        const int n = 25;
        var band = new double[2, n];
        for (int j = 0; j < n; j++)
        {
            band[0, j] = j == 0 ? 1.0 : 5.0;
            if (j < n - 1) band[1, j] = 2.0;
        }
        var factor = new BandCholesky(n, 1, BandFactorization.Ldlt, strict: true);
        factor.DecomposeBanded(band);
        Assert.True(DesignAdmission.PivotsAdmit(factor.Status));
        Assert.False(DesignAdmission.Admits(factor));
    }

    [Fact]
    public void Normal_DesignWhoseFactorBreaksDown_IsRejectedWithoutBuildingTheFit()
    {
        // Four knots crowded at the left edge of 120 even points, as a hot likelihood-tempered rung proposed them. The
        // Gram is not numerically positive definite: the floored factor's pivots ran to −1e155 and its solve was NaN,
        // and the strict factor stops. Admission must judge the factor before any least-squares statistic is formed.
        double[] x = EvenInputs(120);
        var knots = new[] { 0.010348173957892346, 0.012750809965187548, 0.021147317788105583, 0.06881582846542306 };
        double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
        double[] y = Response("normal", x);
        var model = new WeightedNormalModel();

        Assert.Equal(double.NegativeInfinity, model.LogMarginalLikelihood(z, y, null));
        Assert.Equal(double.NegativeInfinity, model.TemperedMarginal(z, y, null).LogMarginalLikelihood(0.5));
        Assert.All(model.PosteriorMeanCoefficients(z, y, null), v => Assert.True(double.IsNaN(v)));
        ConditionalDraw draw = model.SampleConditional(z, y, null, new Maths.Rng.Xoshiro256PlusPlus(1));
        Assert.All(draw.Coefficients, v => Assert.True(double.IsNaN(v)));
        Assert.True(double.IsNaN(draw.Sigma2));
        Assert.All(model.ConditionalVariance(z, y, null, z)!, v => Assert.True(double.IsNaN(v)));
    }

    private static IObservationModel Model(string name) => name switch
    {
        "normal" => new WeightedNormalModel(),
        "poisson" => new PoissonModel(),
        "whittle" => new WhittleModel(),
        _ => new GibbsLossModel(),
    };

    // A response suited to each model: real values, counts, or a positive periodogram.
    private static double[] Response(string model, double[] x)
    {
        var y = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            double s = Math.Sin(2.0 * Math.PI * x[i]) + 0.1 * Math.Cos(17.0 * x[i]);
            y[i] = model switch
            {
                "poisson" => Math.Round(6.0 + 4.0 * s),
                "whittle" => Math.Exp(s) * (1.0 + 0.3 * Math.Cos(23.0 * x[i])),
                _ => s,
            };
        }
        return y;
    }

    private static double[] EvenInputs(int n)
    {
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = (i + 0.5) / n;
        return x;
    }

    private static double[] GappedInputs()
    {
        var x = new double[40];
        for (int i = 0; i < 20; i++) x[i] = 0.4 * (i + 0.5) / 20;
        for (int i = 0; i < 20; i++) x[20 + i] = 0.6 + 0.4 * (i + 0.5) / 20;
        return x;
    }
}
