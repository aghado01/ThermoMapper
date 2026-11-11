using System;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class BassConformalKeptTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData("gaussian")]
    [InlineData("t")]
    [InlineData("wrong")]
    public void KeptBesideCarriedAtSeveralEssFloors(string kind)
    {
        if (!BassQualification.Enabled) return;
        const int n = 99, tests = 30, widths = 4, fits = 5;
        int[] seeds = SeedTree.Derive(2026111103, fits);
        double[] floors = [10, 100, 400];
        var coverage = new double[3]; var inflated = new double[3]; var bounds = new double[3]; var carried = new double[3];
        var flagged = new double[3]; var bands = new double[3]; var width = new double[3]; var inflatedWidth = new double[3];
        double[] grid = Enumerable.Range(0, 161).Select(g => -6 + 0.1 * g).ToArray();
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            var (x0, x1, y) = BassCandidateTests.Draw(n, rng);
            var (t0, t1, yt) = BassCandidateTests.Draw(tests, rng);
            BassPredictor[] Predictors(double[] a, double[] b) => kind == "wrong"
                ? [new ContinuousPredictor("x1", b)] : [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];
            var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
            var fit = Bass.Fit(BassStudentTTests.Config(kind == "t" ? BassNoise.StudentT(3, 1.0 / 3) : null, burnIn: 1000, samples: 500), data);
            var training = fit.PointwiseLogDensity(data);
            var candidates = fit.CandidateLogDensity(data, PredictionCases.Of(Predictors([.. t0, .. t0[..widths]], [.. t1, .. t1[..widths]])),
                [.. yt.Select(v => new[] { v }), .. Enumerable.Repeat(grid, widths)]);
            var old = ConformalBayes.FullStep(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried });
            for (int f = 0; f < floors.Length; f++)
            {
                var options = new ConformalOptions { MarginBounds = MarginBounds.Kept, DegenerateEss = floors[f] };
                var result = ConformalBayes.FullStep(training, candidates, options);
                Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidates, options with { Score = FullStepScore.LeaveOneOut }));
                for (int t = 0; t < tests; t++)
                {
                    var set = result.Sets[t];
                    if (set.AcceptsPlain(0)) coverage[f]++;
                    if (set.AcceptsInflated(0)) inflated[f]++;
                    Assert.Equal(old.Sets[t].PValues, set.PValues);
                }
                for (int t = tests; t < tests + widths; t++)
                {
                    var set = result.Sets[t];
                    bounds[f] += set.CoverageUnderError; carried[f] += old.Sets[t].CoverageUnderError;
                    flagged[f] += set.Margin!.FlaggedCount!.Value; bands[f] += set.Margin.Count!.Value;
                    width[f] += set.Width; inflatedWidth[f] += (set with { Inflated = true }).Width;
                }
            }
        }
        for (int f = 0; f < floors.Length; f++)
            output.WriteLine($"{kind}, floor {floors[f]}, {fits} fits, {fits * tests} truths: coverage {coverage[f] / (fits * tests):F4}, inflated {inflated[f] / (fits * tests):F4}; {fits * widths} fixed grids, width {width[f] / (fits * widths):F3}, inflated {inflatedWidth[f] / (fits * widths):F3}; flagged {flagged[f] / (fits * widths):F2}, band {bands[f] / (fits * widths):F2}, reading {bounds[f] / (fits * widths):F4}, carried {carried[f] / (fits * widths):F4}. Hoff refused: its per-unit ratio bounds and flags are not implemented; these density-score bounds do not apply.");
    }
}
