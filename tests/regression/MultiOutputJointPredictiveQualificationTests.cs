using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputJointPredictiveQualificationTests(ITestOutputHelper output)
{
    [Fact]
    public void JointBandsOnIndependentData_FromNonlinearSurfacesOutsideTheFiniteBasis()
    {
        if (!BassQualification.Enabled) return;
        const int n = 160, test = 400;
        int[] seeds = SeedTree.Derive(20261223, 4);
        var rng = new Xoshiro256PlusPlus(seeds[0]);
        double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] z = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        var sigma = new double[,] { { 0.04, 0.024 }, { 0.024, 0.09 } };
        var root = ReferenceAlgebra.Cholesky(sigma)!;
        var y = new double[n, 2];
        double[] Surface(double a, double b) => [Math.Sin(4 * a) + 0.4 * b * b, 0.7 * Math.Cos(3 * b) + a * b];
        for (int i = 0; i < n; i++) { var f = Surface(x[i], z[i]); var e = MultiOutputCoverageTests.Gaussian(root, rng); for (int k = 0; k < 2; k++) y[i, k] = f[k] + e[k]; }
        var data = new BassData([new ContinuousPredictor("x", x), new ContinuousPredictor("z", z)], new MultiOutputResponse(y, ["first", "second"]));
        rng = new Xoshiro256PlusPlus(seeds[1]);
        double[] tx = [.. Enumerable.Range(0, test).Select(_ => rng.NextDouble())], tz = [.. Enumerable.Range(0, test).Select(_ => rng.NextDouble())];
        var observed = new double[test, 2];
        for (int i = 0; i < test; i++) { var f = Surface(tx[i], tz[i]); var e = MultiOutputCoverageTests.Gaussian(root, rng); for (int k = 0; k < 2; k++) observed[i, k] = f[k] + e[k]; }
        var cases = new PredictionCases([new ContinuousColumn("x", tx), new ContinuousColumn("z", tz)]);
        foreach (bool hw in new[] { false, true })
        {
            var config = new BassConfig { Structure = new BassStructureSpec { MaxBasis = 12, MaxInteraction = 2 },
                OutputCovariance = hw ? OutputCovariancePrior.HuangWand() : null,
                Sampling = new BassSampling { Chains = 2, BurnIn = 1200, Samples = 1500, Thin = 3, MasterSeed = seeds[2] } };
            var fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
            var prediction = BassMultiOutput.Predict(fit, data, cases, new BassPredictOptions { Seed = seeds[3] });
            var bands = prediction.JointBands(); int inside = 0;
            for (int i = 0; i < test; i++) if (Enumerable.Range(0, 2).All(k => observed[i, k] >= bands.Lower[i, k] && observed[i, k] <= bands.Upper[i, k])) inside++;
            double coverage = inside / (double)test;
            output.WriteLine($"HW={hw}: one fitted dataset, nonlinear held-out joint coverage {coverage:F4}, binomial conditional SE {Math.Sqrt(coverage * (1 - coverage) / test):F4}; not dataset-level coverage qualification");
            Assert.InRange(coverage, 0.87, 0.995);
            var covariance = BassMultiOutput.Covariance(fit, data);
            output.WriteLine($"mean covariance [{covariance.Mean![0, 0]:G5}, {covariance.Mean[1, 0]:G5}, {covariance.Mean[1, 1]:G5}], true [0.04,0.024,0.09]");
        }
    }
}
