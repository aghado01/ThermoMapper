using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// CBMA reading each model's evidence from its engine's fit (E16, slice 17's follow-up): exact evidences give the log-weight form's sets
/// and weights bit for bit, with the posterior model probabilities and no fragile case; an evidence with Monte Carlo error widens the
/// probabilities' range and marks the cases whose sets its tilts change, more of them at a larger multiple; evidences read on another
/// response, and malformed ones, are refused; the engines' flags pass through.
/// </summary>
public sealed class ModelEvidenceTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // A model of y − c: each unit's and candidate's density under y − c = βx + N(0, σ²), the grid kept as given.
    private static (PointwiseLogDensity, CandidateLogDensity) Offset(double[] x, double[] y, double[] beta, double sigma, double c, double[] cases, double[] grid)
    {
        (PointwiseLogDensity training, CandidateLogDensity candidates) = ScoresBeyondTheModelTests.Seams(x, [.. y.Select(v => v - c)], beta, sigma, cases,
                                                                                                           [.. grid.Select(v => v - c)]);
        return (training, new CandidateLogDensity(candidates.Values, [.. Enumerable.Repeat(grid, cases.Length)], []));
    }

    // Two noise laws, σ 1 and 1.5, on 8 units so the evidence leaves both plausible (slice 17's fixture), on a grid of step 0.05.
    private static (double[] X, double[] Y, (PointwiseLogDensity, CandidateLogDensity) Narrow, (PointwiseLogDensity, CandidateLogDensity) Wide) Fixture(int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x = [.. Enumerable.Range(0, 8).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.7 * v + 1.2 * Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 161).Select(g => -4.0 + 0.05 * g)];
        double[] cases = [0.5, -1.0, 1.5, 0.0];
        return (x, y, ScoresBeyondTheModelTests.Seams(x, y, ScoresBeyondTheModelTests.Draws(x, y, 1.0, 3000, rng), 1.0, cases, grid),
                ScoresBeyondTheModelTests.Seams(x, y, ScoresBeyondTheModelTests.Draws(x, y, 1.5, 3000, rng), 1.5, cases, grid));
    }

    [Fact]
    public void ExactEvidences_AreTheLogWeightForm()
    {
        var (x, y, narrow, wide) = Fixture(2301);
        string key = ModelEvidence.KeyOf(y);
        double a = ScoresBeyondTheModelTests.LogEvidence(x, y, 1.0), b = ScoresBeyondTheModelTests.LogEvidence(x, y, 1.5);
        // α 0.3: on 8 units k is 2, where α 0.1 would leave k 0 and every set the whole line.
        AveragedPrediction given = ModelAveraging.FullStep([narrow, wide], [a + Math.Log(0.3), b + Math.Log(0.7)], 0.3);
        AveragedPrediction read = ModelAveraging.FullStep([new AveragedModel(narrow.Item1, narrow.Item2, new ModelEvidence(a, 0.0, key), Math.Log(0.3)),
                                                           new AveragedModel(wide.Item1, wide.Item2, new ModelEvidence(b, 0.0, key, ["a flag"]), Math.Log(0.7))], 0.3);
        Assert.Null(given.Evidence);
        for (int t = 0; t < given.Prediction.Sets.Count; t++)
        {
            Assert.False(given.Prediction.Sets[t].WholeSpace);
            Assert.Equal(given.Prediction.Sets[t].PValues, read.Prediction.Sets[t].PValues);
            Assert.Equal(given.ModelWeights[t], read.ModelWeights[t]);
        }
        EvidenceReading evidence = read.Evidence!;
        double exact = 1.0 / (1.0 + Math.Exp(b + Math.Log(0.7) - a - Math.Log(0.3)));
        Assert.Equal(exact, evidence.Posterior[0], 12);
        Assert.Equal(evidence.Posterior, evidence.PosteriorLow);
        Assert.Equal(evidence.Posterior, evidence.PosteriorHigh);
        Assert.All(evidence.Fragile, Assert.False);
        Assert.Equal(["model 1: a flag"], evidence.Flags);
    }

    [Fact]
    public void AnEvidencesError_WidensTheProbabilities_AndMarksTheSetsItMoves()
    {
        // Two models apart in their mean, so their sets differ where their weights move (a noise scale barely moves a set, slice 17):
        // y = βx and y = βx + 1.2, on data from y = 0.7x + 0.6 + N(0, 1), both plausible on 8 units.
        var rng = new Xoshiro256PlusPlus(2302);
        double[] x = [.. Enumerable.Range(0, 8).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.7 * v + 0.6 + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 161).Select(g => -4.0 + 0.05 * g)], cases = [0.5, -1.0, 1.5, 0.0];
        double[] lifted = [.. y.Select(v => v - 1.2)];
        var plain = Offset(x, y, ScoresBeyondTheModelTests.Draws(x, y, 1.0, 3000, rng), 1.0, 0.0, cases, grid);
        var offset = Offset(x, y, ScoresBeyondTheModelTests.Draws(x, lifted, 1.0, 3000, rng), 1.0, 1.2, cases, grid);
        string key = ModelEvidence.KeyOf(y);
        double a = ScoresBeyondTheModelTests.LogEvidence(x, y, 1.0), b = ScoresBeyondTheModelTests.LogEvidence(x, lifted, 1.0);
        int Fragile(double error, double multiple, out EvidenceReading reading)
        {
            reading = ModelAveraging.FullStep([new AveragedModel(plain.Item1, plain.Item2, new ModelEvidence(a, error, key, priors: "β ~ N(0, 1), σ 1")),
                                               new AveragedModel(offset.Item1, offset.Item2, new ModelEvidence(b, error, key))], 0.3, multiple).Evidence!;
            return reading.Fragile.Count(f => f);
        }
        int none = Fragile(1.5, 0.0, out EvidenceReading still), small = Fragile(0.25, 2.0, out EvidenceReading narrowRange), large = Fragile(1.5, 2.0, out EvidenceReading wideRange);
        Assert.Equal(["β ~ N(0, 1), σ 1", ""], wideRange.Priors);
        output.WriteLine($"P(y = βx | Z): {wideRange.Posterior[0]:F3}; over ±0.5 nats [{narrowRange.PosteriorLow[0]:F3}, {narrowRange.PosteriorHigh[0]:F3}], "
                         + $"over ±3 nats [{wideRange.PosteriorLow[0]:F3}, {wideRange.PosteriorHigh[0]:F3}]; fragile cases of 4: {none} at multiple 0, {small} at ±0.5, {large} at ±3");
        Assert.Equal(0, none);
        Assert.Equal(still.Posterior, still.PosteriorHigh);
        Assert.True(wideRange.PosteriorLow[0] < narrowRange.PosteriorLow[0] && narrowRange.PosteriorLow[0] < wideRange.Posterior[0]);
        Assert.True(wideRange.PosteriorHigh[0] > narrowRange.PosteriorHigh[0] && narrowRange.PosteriorHigh[0] > wideRange.Posterior[0]);
        Assert.True(large >= small && large > 0, "a wide tilt should move some set");
    }

    [Fact]
    public void AnotherResponse_AndMalformedEvidence_AreRefused()
    {
        var (x, y, narrow, wide) = Fixture(2303);
        double a = ScoresBeyondTheModelTests.LogEvidence(x, y, 1.0), b = ScoresBeyondTheModelTests.LogEvidence(x, y, 1.5);
        double[] shifted = [.. y.Select(v => v + 1e-12)];
        Assert.Equal(ModelEvidence.KeyOf(y), ModelEvidence.KeyOf([.. y]));
        Assert.NotEqual(ModelEvidence.KeyOf(y), ModelEvidence.KeyOf(shifted));
        // Weights are part of what an evidence reads: weighed or not, and by which weights, gives another key.
        double[] ones = [.. y.Select(_ => 1.0)];
        Assert.NotEqual(ModelEvidence.KeyOf(y), ModelEvidence.KeyOf(y, ones));
        Assert.NotEqual(ModelEvidence.KeyOf(y, ones), ModelEvidence.KeyOf(y, [.. ones.Select((w, i) => i == 0 ? 2.0 : w)]));
        Assert.Equal(ModelEvidence.KeyOf(y, ones), ModelEvidence.KeyOf([.. y], [.. ones]));
        Assert.Contains("another response", Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep(
            [new AveragedModel(narrow.Item1, narrow.Item2, new ModelEvidence(a, 0.0, ModelEvidence.KeyOf(y))),
             new AveragedModel(wide.Item1, wide.Item2, new ModelEvidence(b, 0.0, ModelEvidence.KeyOf(shifted)))])).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => new ModelEvidence(double.NaN, 0.0, "k"));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ModelEvidence(double.PositiveInfinity, 0.0, "k"));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ModelEvidence(1.0, -0.1, "k"));
        Assert.Throws<ArgumentException>(() => new ModelEvidence(1.0, 0.0, ""));
        var model = new AveragedModel(narrow.Item1, narrow.Item2, new ModelEvidence(a, 0.0, "k"));
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([model with { LogPrior = double.PositiveInfinity }]));
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([model with { LogPrior = double.NegativeInfinity }]));
        Assert.Throws<ArgumentOutOfRangeException>(() => ModelAveraging.FullStep([model], errorMultiple: -1.0));
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep(Array.Empty<AveragedModel>()));
    }
}
