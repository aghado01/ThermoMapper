using System;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The fitted <see cref="GpTransform"/> shared by training and prediction. The response is centered on its training
/// mean and restored in prediction, so a constant offset changes neither the evidence nor the bandwidth posterior.
/// The default isotropic similarity maps the inputs into the unit cube, so the fit does not depend on the units or
/// origin of the coordinates, and a bandwidth t in the mapped frame is s²·t in the original units. Without these,
/// an offset of 20 moved the log-evidence at t = 0.05 from 36.3 to −861.6 and the posterior-mean bandwidth from 0.092
/// to 0.542, and rescaling the inputs to 100x − 7 emptied the prior's support.
/// </summary>
public sealed class BgpTransformTests
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // 50 points on [0,1] × [0,0.4] with a smooth response, plus three test points.
    private static (double[,] X, double[] Y, double[,] XTest) Data()
    {
        const int n = 50;
        var rng = new Xoshiro256PlusPlus(seed: 8);
        var x = new double[n, 2];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = rng.NextDouble();
            x[i, 1] = 0.4 * rng.NextDouble();
            y[i] = Math.Sin(5.0 * x[i, 0]) + Math.Cos(7.0 * x[i, 1]) + 0.05 * Gaussian(rng);
        }
        var xTest = new double[,] { { 0.2, 0.1 }, { 0.5, 0.3 }, { 0.9, 0.05 } };
        return (x, y, xTest);
    }

    private static double[,] Affine(double[,] x, double a, double b)
    {
        var m = (double[,])x.Clone();
        for (int i = 0; i < m.GetLength(0); i++)
            for (int k = 0; k < m.GetLength(1); k++) m[i, k] = a * x[i, k] + b;
        return m;
    }

    private static void Close(double expected, double actual, double rel, string what)
        => Assert.True(Math.Abs(expected - actual) <= rel * Math.Max(1.0, Math.Abs(expected)), $"{what}: {expected} vs {actual}");

    [Fact]
    public void ResponseOffset_IsRemovedInFitting_AndRestoredInPrediction()
    {
        (double[,] x, double[] y, double[,] xTest) = Data();
        var yShift = new double[y.Length];
        for (int i = 0; i < y.Length; i++) yShift[i] = y[i] + 20.0;

        var kernel = new SquaredExponentialKernel();
        var gp = new GpRegression(x, y, 0.0025, kernel);
        var gpShift = new GpRegression(x, yShift, 0.0025, kernel);
        Close(gp.Transform.ResponseOffset + 20.0, gpShift.Transform.ResponseOffset, 1e-12, "response offset");

        GpFit fit = gp.Fit(0.05), fitShift = gpShift.Fit(0.05);
        Close(fit.LogMarginal, fitShift.LogMarginal, 1e-9, "log-evidence");
        double[] pred = gp.PredictMean(fit, xTest), predShift = gpShift.PredictMean(fitShift, xTest);
        for (int i = 0; i < pred.Length; i++) Close(pred[i] + 20.0, predShift[i], 1e-9, $"prediction {i}");

        BgpResult r = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp)).Run(xTest, draws: 300, burn: 100, seed: 2);
        BgpResult rShift = new BgpSampler(gpShift, new EmpiricalBayesBandwidthPrior(gpShift)).Run(xTest, draws: 300, burn: 100, seed: 2);
        Close(r.BandwidthMean, rShift.BandwidthMean, 1e-9, "posterior-mean bandwidth");
        for (int i = 0; i < pred.Length; i++) Close(r.PosteriorMean[i] + 20.0, rShift.PosteriorMean[i], 1e-8, $"posterior mean {i}");
    }

    [Fact]
    public void SimilarityScaling_MakesTheFitIndependentOfUnitsAndOrigin()
    {
        (double[,] x, double[] y, double[,] xTest) = Data();
        double[,] xUnits = Affine(x, 100.0, -7.0), xTestUnits = Affine(xTest, 100.0, -7.0);

        var kernel = new SquaredExponentialKernel();
        var gp = new GpRegression(x, y, 0.0025, kernel);
        var gpUnits = new GpRegression(xUnits, y, 0.0025, kernel);
        Close(100.0 * gp.Transform.InputScale, gpUnits.Transform.InputScale, 1e-12, "input scale");

        var prior = new EmpiricalBayesBandwidthPrior(gp);
        var priorUnits = new EmpiricalBayesBandwidthPrior(gpUnits);
        Close(prior.Tn, priorUnits.Tn, 1e-12, "T_n");
        Close(prior.LowerBound, priorUnits.LowerBound, 1e-12, "support lower bound");

        foreach (double t in new[] { 0.01, 0.05, 0.3 })
        {
            GpFit f = gp.Fit(t), fu = gpUnits.Fit(t);
            Close(f.LogMarginal, fu.LogMarginal, 1e-9, $"log-evidence at t={t}");
            double[] p = gp.PredictMean(f, xTest), pu = gpUnits.PredictMean(fu, xTestUnits);
            for (int i = 0; i < p.Length; i++) Close(p[i], pu[i], 1e-9, $"prediction {i} at t={t}");
        }

        BgpResult r = new BgpSampler(gp, prior).Run(xTest, draws: 300, burn: 100, seed: 4);
        BgpResult ru = new BgpSampler(gpUnits, priorUnits).Run(xTestUnits, draws: 300, burn: 100, seed: 4);
        Close(r.BandwidthMean, ru.BandwidthMean, 1e-9, "posterior-mean bandwidth");
        Close(100.0 * r.InputScale, ru.InputScale, 1e-12, "reported input scale");
    }

    [Fact]
    public void SimilarityBandwidth_ConvertsToOriginalUnits_ByTheSquaredScale()
    {
        // Spans 10 and 1: one isotropic scale s = 10, so a mapped bandwidth t equals s²·t on the raw coordinates.
        var rng = new Xoshiro256PlusPlus(seed: 21);
        const int n = 40;
        var x = new double[n, 2];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = 10.0 * i / (n - 1);
            x[i, 1] = rng.NextDouble();
            y[i] = Math.Sin(x[i, 0]) + 0.1 * Gaussian(rng);
        }
        var kernel = new SquaredExponentialKernel();
        var mapped = new GpRegression(x, y, 0.01, kernel);
        var raw = new GpRegression(x, y, 0.01, kernel, GpInputScaling.None);

        double s = mapped.Transform.InputScale;
        Assert.Equal(10.0, s, 12);
        Assert.Equal(0.0, mapped.Transform.InputOffset[0], 12);
        foreach (double t in new[] { 0.001, 0.01, 0.1 })
            Close(raw.Fit(s * s * t).LogMarginal, mapped.Fit(t).LogMarginal, 1e-9, $"log-evidence at t={t}");
    }

    [Fact]
    public void UnitFrameInputs_MapByTheIdentity_AndNoneKeepsInputsAsGiven()
    {
        const int n = 30;
        var x = new double[n, 2];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = i / (n - 1.0);                 // min 0, max 1: the widest coordinate already spans [0, 1]
            x[i, 1] = 0.4 * ((i * 7) % n) / (n - 1.0);   // min 0
            y[i] = Math.Cos(3.0 * x[i, 0]) + x[i, 1];
        }
        var kernel = new SquaredExponentialKernel();
        var similarity = new GpRegression(x, y, 0.01, kernel);
        var none = new GpRegression(x, y, 0.01, kernel, GpInputScaling.None);

        Assert.Equal(1.0, similarity.Transform.InputScale);
        Assert.Equal(1.0, none.Transform.InputScale);
        Assert.Equal(new[] { 0.0, 0.0 }, none.Transform.InputOffset);
        foreach (double t in new[] { 0.01, 0.1 })
            Assert.Equal(none.Fit(t).LogMarginal, similarity.Fit(t).LogMarginal);
    }

    [Fact]
    public void NonFiniteResponse_IsRejected()
    {
        var x = new double[,] { { 0.0 }, { 0.5 }, { 1.0 } };
        Assert.Throws<ArgumentException>(() => new GpRegression(x, new[] { 0.1, double.NaN, 0.2 }, 0.01, new SquaredExponentialKernel()));
    }
}
