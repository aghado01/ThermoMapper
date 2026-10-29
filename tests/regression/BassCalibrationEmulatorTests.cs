using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Calibration's emulator contract (stage 7, C-X2; D09): each adapter's prediction against its fit's own
/// <see cref="BassFit.Predict"/> or <see cref="BassBasisFit.Predict"/> for the same draw, bit for bit, at random θ
/// and random draws, for a scalar fit with fixed inputs, a categorical θ, a direct functional fit and a basis fit; the
/// level enumeration against single predictions; the noise covariance and the truncation residuals in physical units;
/// and the refusals.
/// </summary>
public sealed class BassCalibrationEmulatorTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationEmulatorTests(ITestOutputHelper output) => _out = output;

    // A basis fit with a continuous x and a categorical z on three levels, curves at 12 locations (n = 18).
    private static readonly Lazy<BassBasisFit> CategoricalBasisFit = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(17);
        const int n = 18, m = 12;
        var x = new double[n];
        var z = new int[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i / 3) / 5.0;
            z[i] = i % 3;
            for (int l = 0; l < m; l++)
            {
                double t = l / (m - 1.0);
                y[i, l] = Math.Sin(2 * Math.PI * (t - 0.3 * x[i])) + (z[i] == 1 ? 0.8 * t : 0) - (z[i] == 2 ? 0.5 : 0) + 0.05 * Normal.SampleStandard(rng);
            }
        }
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.95, ScaleByColumnSd = true });
        return BassComponents.Fit(d, [new ContinuousPredictor("x", x), new CategoricalPredictor("z", z, ["a", "b", "c"])],
                                  BassRunTests.Config(chains: 2, burn: 150, samples: 40, seed: 5, maxBasis: 5));
    });

    private static double[] Uniforms(Xoshiro256PlusPlus rng, int count, double lower, double upper)
        => [.. Enumerable.Range(0, count).Select(_ => lower + (upper - lower) * rng.NextDouble())];

    [Fact]
    public void ScalarFit_WithFixedInputs_PredictsAsTheFit()
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Scalar);
        var theta = new ContinuousParameter("x0", -0.1, 1.2);   // beyond the training range: the hinges extrapolate
        double[] x1 = [0.0, 0.3, 0.55, 0.8, 1.0];
        var emulator = new BassFitEmulator(fit, data, [theta], new PredictionCases([new ContinuousColumn("x1", x1)]));
        Assert.Equal(5, emulator.Outputs);
        Assert.Equal(["x0"], emulator.Extrapolating);
        Assert.Equal(fit.Draws.Sum(c => c.Count), emulator.DrawCount);

        var rng = new Xoshiro256PlusPlus(3);
        var output = new double[5];
        for (int trial = 0; trial < 12; trial++)
        {
            int s = rng.NextInt(emulator.DrawCount);
            double t = theta.ToPhysical(rng.NextDouble());
            emulator.Predict(s, [t], [], output);
            var cases = new PredictionCases([new ContinuousColumn("x0", [.. Enumerable.Repeat(t, 5)]), new ContinuousColumn("x1", x1)]);
            BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Draws = [s], Predictive = false });
            for (int i = 0; i < 5; i++) Assert.Equal(p.Latent[0, i, 0], output[i]);
            Assert.Equal(new IndependentNoise(p.NoiseVariances[0]), emulator.Noise(s));
        }
        Assert.Null(emulator.TruncationResiduals());
    }

    [Fact]
    public void CategoricalTheta_PredictsAsTheFit_AndEnumeratesTheLevels()
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Categorical);
        ContinuousParameter x0 = new("x0", 0.0, 1.0);
        CategoricalParameter z = new("z", ["c", "a", "b"]);   // declared in another order than the training levels
        var emulator = new BassFitEmulator(fit, data, [z, x0]);
        Assert.Equal(1, emulator.Outputs);

        var rng = new Xoshiro256PlusPlus(4);
        var output = new double[1];
        int[][] all = [[0], [1], [2]];
        var enumerated = new double[3, 1];
        int nonzeroCategorical = 0;
        for (int trial = 0; trial < 15; trial++)
        {
            int s = rng.NextInt(emulator.DrawCount);
            double t = rng.NextDouble();
            emulator.PredictLevels(s, [t], all, enumerated);
            for (int level = 0; level < 3; level++)
            {
                emulator.Predict(s, [t], [level], output);
                var cases = new PredictionCases([new ContinuousColumn("x0", [t]), new CategoricalColumn("z", [z.Levels[level]])]);
                BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Draws = [s], Predictive = false });
                Assert.Equal(p.Latent[0, 0, 0], output[0]);
                Assert.Equal(output[0], enumerated[level, 0]);
            }
            if (enumerated[0, 0] != enumerated[1, 0] || enumerated[1, 0] != enumerated[2, 0]) nonzeroCategorical++;
        }
        Assert.True(nonzeroCategorical > 0, "No draw's terms read z, so the enumeration was not exercised.");
    }

    [Fact]
    public void FixedCategoricalInput_PredictsAsTheFit_AndWeightsPassThrough()
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Categorical);
        ContinuousParameter x0 = new("x0", 0.0, 1.0);
        string[] z = ["b", "c", "a", "b"];
        var emulator = new BassFitEmulator(fit, data, [x0], new PredictionCases([new CategoricalColumn("z", z)]));
        Assert.Null(emulator.DrawWeights);
        var rng = new Xoshiro256PlusPlus(12);
        var output = new double[4];
        for (int trial = 0; trial < 10; trial++)
        {
            int s = rng.NextInt(emulator.DrawCount);
            double t = rng.NextDouble();
            emulator.Predict(s, [t], [], output);
            BassPrediction p = fit.Predict(data, new PredictionCases([new ContinuousColumn("x0", [t, t, t, t]), new CategoricalColumn("z", z)]),
                                           new BassPredictOptions { Draws = [s], Predictive = false });
            for (int i = 0; i < 4; i++) Assert.Equal(p.Latent[0, i, 0], output[i]);
        }

        double[] weights = [.. Enumerable.Range(1, emulator.DrawCount).Select(d => (double)d)];
        var weighted = new BassFitEmulator(fit with { DrawWeights = weights }, data, [x0], new PredictionCases([new CategoricalColumn("z", z)]));
        Assert.Equal(weights, weighted.DrawWeights);
    }

    [Fact]
    public void DirectFunctionalFit_PredictsAsTheFit_AtDeclaredLocations()
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Functional);
        ContinuousParameter x0 = new("x0", 0.0, 1.0);
        ContinuousColumn t = new("t", [0.05, 0.4, 0.75, 0.95]);
        var training = new BassFitEmulator(fit, data, [x0], new PredictionCases([new ContinuousColumn("x1", [0.5])]));
        var declared = new BassFitEmulator(fit, data, [x0], new PredictionCases([new ContinuousColumn("x1", [0.2, 0.9])]), [t]);
        Assert.Equal(6, training.Outputs);
        Assert.Equal(8, declared.Outputs);

        var rng = new Xoshiro256PlusPlus(6);
        for (int trial = 0; trial < 10; trial++)
        {
            int s = rng.NextInt(declared.DrawCount);
            double u = rng.NextDouble();
            var output = new double[8];
            declared.Predict(s, [u], [], output);
            var cases = new PredictionCases([new ContinuousColumn("x0", [u, u]), new ContinuousColumn("x1", [0.2, 0.9])]);
            BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Draws = [s], Predictive = false, OutputLocations = [t] });
            for (int i = 0; i < 2; i++)
                for (int l = 0; l < 4; l++) Assert.Equal(p.Latent[0, i, l], output[i * 4 + l]);
            Assert.Equal(new IndependentNoise(p.NoiseVariances[0]), declared.Noise(s));

            var atTraining = new double[6];
            training.Predict(s, [u], [], atTraining);
            BassPrediction q = fit.Predict(data, new PredictionCases([new ContinuousColumn("x0", [u]), new ContinuousColumn("x1", [0.5])]),
                                           new BassPredictOptions { Draws = [s], Predictive = false });
            for (int l = 0; l < 6; l++) Assert.Equal(q.Latent[0, 0, l], atTraining[l]);
        }
    }

    [Fact]
    public void BasisFit_PredictsAsTheFit_WithItsNoiseAndResidualsInPhysicalUnits()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        ContinuousParameter z0 = new("z0", 0.0, 1.0), z1 = new("z1", 0.0, 1.0);
        var emulator = new BassBasisEmulator(fit, [z1, z0]);
        int m = fit.Decomposition.Locations, r = fit.Decomposition.Rank;
        Assert.Equal(m, emulator.Outputs);
        Assert.Equal(fit.JointDraws.Count, emulator.DrawCount);

        var rng = new Xoshiro256PlusPlus(8);
        var output = new double[m];
        for (int trial = 0; trial < 12; trial++)
        {
            int s = rng.NextInt(emulator.DrawCount);
            double a = rng.NextDouble(), b = rng.NextDouble();
            emulator.Predict(s, [a, b], [], output);
            BassBasisPrediction p = fit.Predict(new PredictionCases([new ContinuousColumn("z0", [b]), new ContinuousColumn("z1", [a])]),
                                                new BassPredictOptions { Draws = [s], Predictive = false });
            for (int l = 0; l < m; l++) Assert.Equal(p.Latent[0, 0, l], output[l]);
            var noise = Assert.IsType<ComponentNoise>(emulator.Noise(s));
            for (int k = 0; k < r; k++) Assert.Equal(p.NoiseVariances[0, k], noise.Variances[k]);
            for (int l = 0; l < m; l++)
                for (int k = 0; k < r; k++) Assert.Equal(fit.Decomposition.Map.Scales[l] * fit.Decomposition.Map.Basis(l, k), noise.Factor[l, k]);
        }

        double[,] residuals = emulator.TruncationResiduals()!;
        double[,] standardized = fit.Decomposition.Truncation.Values;
        Assert.Equal(standardized.GetLength(0), residuals.GetLength(0));
        for (int i = 0; i < residuals.GetLength(0); i++)
            for (int l = 0; l < m; l++) Assert.Equal(fit.Decomposition.Map.Scales[l] * standardized[i, l], residuals[i, l]);
        residuals[0, 0] = double.NaN;   // a copy
        Assert.False(double.IsNaN(emulator.TruncationResiduals()![0, 0]));
    }

    [Fact]
    public void BasisFit_WithACategoricalTheta_EnumeratesAsItPredicts()
    {
        BassBasisFit fit = CategoricalBasisFit.Value;
        ContinuousParameter x = new("x", 0.0, 1.0);
        CategoricalParameter z = new("z", ["a", "b", "c"]);
        var emulator = new BassBasisEmulator(fit, [x, z]);
        int m = emulator.Outputs;
        var rng = new Xoshiro256PlusPlus(9);
        int[][] all = [[0], [1], [2]];
        var enumerated = new double[3, m];
        var output = new double[m];
        for (int trial = 0; trial < 8; trial++)
        {
            int s = rng.NextInt(emulator.DrawCount);
            double u = rng.NextDouble();
            emulator.PredictLevels(s, [u], all, enumerated);
            for (int level = 0; level < 3; level++)
            {
                emulator.Predict(s, [u], [level], output);
                BassBasisPrediction p = fit.Predict(new PredictionCases([new ContinuousColumn("x", [u]), new CategoricalColumn("z", [z.Levels[level]])]),
                                                    new BassPredictOptions { Draws = [s], Predictive = false });
                for (int l = 0; l < m; l++)
                {
                    Assert.Equal(p.Latent[0, 0, l], output[l]);
                    Assert.Equal(output[l], enumerated[level, l]);
                }
            }
        }
    }

    [Fact]
    public void Declarations_ThatDoNotMatchTheFit_AreRefused()
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Categorical);
        ContinuousParameter x0 = new("x0", 0.0, 1.0);
        Assert.Contains("neither", Assert.Throws<ArgumentException>(() => new BassFitEmulator(fit, data, [x0])).Message);
        Assert.Contains("no training case", Assert.Throws<ArgumentException>(() =>
            new BassFitEmulator(fit, data, [x0, new CategoricalParameter("z", ["a", "d"])])).Message);
        Assert.Contains("names no predictor", Assert.Throws<ArgumentException>(() =>
            new BassFitEmulator(fit, data, [x0, new CategoricalParameter("z", ["a", "b"]), new ContinuousParameter("w", 0, 1)])).Message);
        Assert.Contains("fixed input", Assert.Throws<ArgumentException>(() =>
            new BassFitEmulator(fit, data, [x0, new CategoricalParameter("z", ["a", "b"])], new PredictionCases([new ContinuousColumn("x0", [0.5])]))).Message);
        Assert.Contains("kind", Assert.Throws<ArgumentException>(() =>
            new BassFitEmulator(fit, data, [x0, new ContinuousParameter("z", 0, 1)])).Message);
        Assert.Contains("twice", Assert.Throws<ArgumentException>(() => new ParameterLayout([x0, x0])).Message);
        Assert.Contains("bounds", Assert.Throws<ArgumentException>(() => new ParameterLayout([new ContinuousParameter("a", 1, 1)])).Message);

        var emulator = new BassFitEmulator(fit, data, [x0, new CategoricalParameter("z", ["a", "b"])]);
        Assert.Throws<ArgumentOutOfRangeException>(() => emulator.Predict(0, [0.5], [2], new double[1]));
        Assert.Throws<ArgumentOutOfRangeException>(() => emulator.Predict(emulator.DrawCount, [0.5], [0], new double[1]));
        Assert.Throws<ArgumentException>(() => emulator.Predict(0, [double.NaN], [0], new double[1]));

        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        Assert.Contains("one case", Assert.Throws<ArgumentException>(() =>
            new BassBasisEmulator(basis, [new ContinuousParameter("z0", 0, 1)], new PredictionCases([new ContinuousColumn("z1", [0.1, 0.2])]))).Message);
        _out.WriteLine("refusals: unmatched predictor, unseen level, unknown parameter, parameter given as fixed, kind mismatch, duplicate, empty bounds");
    }

    [Fact]
    public void Layout_EnumeratesCombinations_FirstParameterSlowest()
    {
        var layout = new ParameterLayout([new CategoricalParameter("a", ["0", "1"]), new ContinuousParameter("x", 0, 1), new CategoricalParameter("b", ["0", "1", "2"])]);
        Assert.Equal(6, layout.Combinations);
        var levels = new int[2];
        var seen = new List<string>();
        for (long c = 0; c < layout.Combinations; c++)
        {
            layout.Levels(c, levels);
            seen.Add($"{levels[0]}{levels[1]}");
        }
        Assert.Equal(["00", "01", "02", "10", "11", "12"], seen);
        Assert.Equal(1, new ParameterLayout([new ContinuousParameter("x", 0, 1)]).Combinations);
    }
}
