using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The split form on BASS fits (E16, slice 3; CB-2, <c>split_coverage</c>), with the qualification switch
/// (<c>BASS_QUALIFICATION</c> = 1). Held-out coverage on data from outside every model (<see cref="BassCandidateTests.Draw"/>) for
/// a Gaussian chain stopped early, a deliberately wrong fit that never sees x₀, and a healthy fit, each under the density score
/// through the candidate seam, CQR's on the fit's predictive quantiles and the standardized residual's; and CQR on a pair of
/// quantile fits, its natural use. Each set's width beside the fit's own uncalibrated 90% interval. One fact per fit, so a
/// qualification run can take one at a time.
/// </summary>
public sealed class BassSplitConformalTests(ITestOutputHelper output)
{
    // m = 200 and k = 20 at α 0.1, nominal 0.900; 30 fits of 100 test cases each.
    private const int Train = 99, CalibrationUnits = 199, Tests = 100, Replications = 30, Widths = 10;

    [Fact]
    public void Coverage_AChainStoppedEarly() => Coverage("a chain stopped early");

    [Fact]
    public void Coverage_AWrongFit() => Coverage("a wrong fit");

    [Fact]
    public void Coverage_AHealthyFit() => Coverage("a healthy fit");

    [Fact]
    public void Coverage_CqrOnAPairOfQuantileFits() => Coverage("a pair of quantile fits");

    private void Coverage(string kind)
    {
        if (!BassQualification.Enabled) return;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        int[] seeds = SeedTree.Derive(2026_1108, Replications);
        bool pair = kind == "a pair of quantile fits";
        string[] scores = pair ? ["CQR"] : ["density", "CQR", "residual"];
        var covered = new int[scores.Length];
        var width = new double[scores.Length];
        double ownCovered = 0, ownWidth = 0, seconds = 0;
        int split = 0, open = 0, empty = 0;
        for (int r = 0; r < Replications; r++)
        {
            long t0Tick = System.Diagnostics.Stopwatch.GetTimestamp();
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(Train, rng);
            (double[] c0, double[] c1, double[] yc) = BassCandidateTests.Draw(CalibrationUnits, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(Tests, rng);
            bool blind = kind == "a wrong fit";
            BassPredictor[] Predictors(double[] a, double[] b)
                => blind ? [new ContinuousPredictor("x1", b)] : [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];
            var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
            PredictionCases calibrationCases = PredictionCases.Of(Predictors(c0, c1)), testCases = PredictionCases.Of(Predictors(t0, t1));
            var lean = new BassPredictOptions { Predictive = false };
            static double[] Column(double[,] a) => [.. Enumerable.Range(0, a.GetLength(0)).Select(i => a[i, 0])];

            ConformalIntervals cqr;
            double[] ownLower, ownUpper;
            if (pair)
            {
                // Quantile BASS at 0.05 and 0.95: each surface is its quantile, read as the posterior mean.
                BassFit lower = Bass.Fit(BassStudentTTests.Config(BassNoise.AsymmetricLaplace(0.05), burnIn: 1000, samples: 500), data);
                BassFit upper = Bass.Fit(BassStudentTTests.Config(BassNoise.AsymmetricLaplace(0.95), burnIn: 1000, samples: 500), data);
                QuantileReadings Surfaces(PredictionCases cases)
                    => new(Column(lower.Predict(data, cases, lean).Mean()), Column(upper.Predict(data, cases, lean).Mean()));
                QuantileReadings atTests = Surfaces(testCases);
                cqr = SplitConformal.Quantiles(Surfaces(calibrationCases), yc, atTests);
                (ownLower, ownUpper) = ([.. atTests.Lower], [.. atTests.Upper]);
            }
            else
            {
                BassConfig config = kind == "a chain stopped early"
                    ? BassStudentTTests.Config(null, burnIn: 0, samples: 20, chains: 1)
                    : BassStudentTTests.Config(null, burnIn: 1000, samples: 500);
                BassFit fit = Bass.Fit(config, data);
                // The density score through the seam: the calibration cases at their responses, the test cases at theirs.
                CandidateLogDensity calibration = fit.CandidateLogDensity(data, calibrationCases, [.. yc.Select(v => new[] { v })]);
                ConformalPrediction density = SplitConformal.Density(calibration, fit.CandidateLogDensity(data, testCases, [.. yt.Select(v => new[] { v })]));
                for (int t = 0; t < Tests; t++)
                    if (density.Sets[t].AcceptsPlain(0)) covered[0]++;
                PredictionCases few = PredictionCases.Of(Predictors(t0[..Widths], t1[..Widths]));
                foreach (ConformalSet set in SplitConformal.Density(calibration, fit.CandidateLogDensity(data, few, [.. Enumerable.Repeat(grid, Widths)])).Sets)
                {
                    width[0] += set.Width;
                    if (set.Intervals.Count > 1) split++;
                    if (set.OpenBelow || set.OpenAbove) open++;
                }
                // CQR on the fit's predictive quantiles; the residual standardized by its interquartile range.
                BassPrediction atCalibration = fit.Predict(data, calibrationCases, lean), atTests = fit.Predict(data, testCases, lean);
                QuantileReadings Quantiles(BassPrediction p) => new(Column(p.Quantile(0.05)), Column(p.Quantile(0.95)));
                ResidualReadings Residuals(BassPrediction p)
                {
                    double[] q1 = Column(p.Quantile(0.25)), q3 = Column(p.Quantile(0.75));
                    return new(Column(p.ObservationMean()), [.. q1.Zip(q3, (a, b) => b - a)]);
                }
                QuantileReadings ownQuantiles = Quantiles(atTests);
                cqr = SplitConformal.Quantiles(Quantiles(atCalibration), yc, ownQuantiles);
                ConformalIntervals residual = SplitConformal.Residuals(Residuals(atCalibration), yc, Residuals(atTests));
                for (int t = 0; t < Tests; t++)
                {
                    if (residual.Intervals[t].Contains(yt[t])) covered[2]++;
                    width[2] += residual.Intervals[t].Width;
                }
                (ownLower, ownUpper) = ([.. ownQuantiles.Lower], [.. ownQuantiles.Upper]);
            }
            int at = pair ? 0 : 1;
            for (int t = 0; t < Tests; t++)
            {
                ConformalInterval interval = cqr.Intervals[t];
                Assert.Equal((200, 20), (interval.Units, interval.Level));
                if (interval.Contains(yt[t])) covered[at]++;
                width[at] += interval.Width;
                if (interval.Empty) empty++;
                if (yt[t] >= ownLower[t] && yt[t] <= ownUpper[t]) ownCovered++;
                ownWidth += Math.Max(0.0, ownUpper[t] - ownLower[t]);
            }
            seconds += (System.Diagnostics.Stopwatch.GetTimestamp() - t0Tick) / (double)System.Diagnostics.Stopwatch.Frequency;
        }
        int count = Replications * Tests;
        double se = SplitConformalTests.CoverageError(200, 20, Replications, Tests);
        output.WriteLine($"{kind}: its own uncalibrated 90% interval covers {ownCovered / count:F4}, mean width {ownWidth / count:F3}; {seconds:F0} s over {Replications} replications");
        for (int c = 0; c < scores.Length; c++)
        {
            double coverage = covered[c] / (double)count;
            double meanWidth = scores[c] == "density" ? width[c] / (Replications * Widths) : width[c] / count;
            output.WriteLine($"  {scores[c]}: coverage {coverage:F4} (nominal 0.9000, se {se:F4}), mean width {meanWidth:F3}");
        }
        if (!pair) output.WriteLine($"  density sets of two or more intervals {split} of {Replications * Widths}, open at the grid's ends {open}");
        output.WriteLine($"  CQR intervals empty {empty} of {count}");
        for (int c = 0; c < scores.Length; c++)
            Assert.True(Math.Abs(covered[c] / (double)count - 0.9) <= 4 * se, $"{kind}, {scores[c]}: {covered[c] / (double)count:F4}");
    }
}
