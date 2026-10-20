using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// R08 / Q4: every BARS observation model admits a design only when its factored matrix is exact (no floored
/// pivot), positive definite and not proven ill-conditioned; anything else scores −∞. The IRLS models also
/// reject a fit that did not converge.
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
        // Columns 1, x and x + 1e-5·cos 7x: the last pivot is ≈1e-9 against ≈40, so the pivot ratio is ≈3e-11,
        // below the admission threshold, while no pivot reaches the absolute floor.
        double[] x = EvenInputs(40);
        var z = new double[x.Length, 3];
        for (int i = 0; i < x.Length; i++)
        {
            z[i, 0] = 1.0;
            z[i, 1] = x[i];
            z[i, 2] = x[i] + 1e-5 * Math.Cos(7.0 * x[i]);
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
