using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A basis fit's presentation data (stage 9, O9; P09, P11, P13, P15): the per-component traces; fitted against observed in
/// physical units against <see cref="BassBasisFit.Predict"/> at the training cases, the observed responses the training
/// data the decomposition represents; the metadata round-tripping the decomposition's account, the joint draws and each
/// component's own; and the per-component posterior summaries.
/// </summary>
public sealed class BassPresentationBasisTests
{
    private readonly ITestOutputHelper _out;
    public BassPresentationBasisTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Traces_AreEachComponentsOwn()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        var traces = fit.Traces();
        Assert.Equal(fit.Components.Count, traces.Count);
        for (int k = 0; k < traces.Count; k++)
        {
            BassTraces own = fit.Components[k].Traces();
            Assert.Equal(own.Chains.Count, traces[k].Chains.Count);
            for (int c = 0; c < own.Chains.Count; c++)
            {
                Assert.Equal(own.Chains[c].Count, traces[k].Chains[c].Count);
                Assert.Equal(own.Chains[c].Sigma2, traces[k].Chains[c].Sigma2);
            }
        }
    }

    [Fact]
    public void TrainingFit_IsThePredictionAtTheTrainingCases_InPhysicalUnits()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassTrainingFit t = fit.TrainingFit(seed: 4);
        BassBasisPrediction p = fit.Predict(PredictionCases.Of(predictors), new BassPredictOptions { Seed = 4 });
        int n = y.GetLength(0), m = y.GetLength(1);
        Assert.Equal((n, m), (t.Cases, t.Locations));
        Assert.Equal(IntervalSource.PredictiveDraws, t.Intervals);
        Assert.Null(t.NoiseSd);
        Assert.Null(t.LocationCoordinates);
        Assert.Equal(p.Mean(), t.Mean);
        double worst = 0.0;
        for (int q = 0; q < t.Probabilities.Count; q++)
        {
            double[,] at = p.Quantile(t.Probabilities[q]);
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) Assert.Equal(at[i, l], t.Quantiles[q, i, l]);
        }
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                // The decomposition reconstructs its training data exactly (D07).
                Assert.Equal(y[i, l], t.Observed[i, l], 12);
                Assert.Equal(t.Observed[i, l] - t.Mean[i, l], t.Residuals[i, l]);
                worst = Math.Max(worst, Math.Abs(t.Residuals[i, l]));
            }
        _out.WriteLine($"largest residual {worst:F3}; y spans {y.Cast<double>().Min():F2} to {y.Cast<double>().Max():F2}");
        // Another seed draws other predictive draws, so other quantiles, and the same mean.
        BassTrainingFit other = fit.TrainingFit(seed: 5);
        Assert.Equal(t.Mean, other.Mean);
        Assert.NotEqual(t.Quantiles[0, 0, 0], other.Quantiles[0, 0, 0]);
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.TrainingFit([1.0]));
    }

    [Fact]
    public void Metadata_RoundTripsTheDecompositionTheJointDrawsAndTheComponents()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassBasisMetadata meta = fit.Metadata();
        BassDecomposition d = fit.Decomposition;
        Assert.Equal((24, 15), (meta.Data.Cases, meta.Data.Locations));
        Assert.Equal(new[] { "z0", "z1" }, meta.Data.Predictors.Select(p => p.Name));
        Assert.Empty(meta.Data.OutputVariables);
        Assert.Equal(d.Rank, meta.Components);
        Assert.Same(d.Pca, meta.Pca);
        Assert.Equal(new BassTruncationSummary(d.Truncation.Kind, d.Truncation.Rms, d.Truncation.InSampleRms, d.Truncation.LeaveOneOutRms,
                                               d.Truncation.StabilityRatio), meta.Truncation);
        Assert.Equal(TruncationResidual.LeaveOneCaseOut, meta.Truncation.Kind);
        Assert.Equal(fit.MasterSeed, meta.MasterSeed);
        Assert.Equal(fit.JointDraws.Count, meta.JointDraws);
        Assert.Same(fit.JointDraws.ChainLengths, meta.JointChainLengths);
        Assert.False(meta.Resampled);
        Assert.Equal(fit.Components.Count, meta.ComponentFits.Count);
        for (int k = 0; k < fit.Components.Count; k++)
        {
            Assert.Same(fit.Components[k].Config, meta.ComponentFits[k].Config);
            Assert.Same(fit.Components[k].Identity, meta.ComponentFits[k].Identity);
            Assert.Equal((24, 1), (meta.ComponentFits[k].Data.Cases, meta.ComponentFits[k].Data.Locations));
        }
        _out.WriteLine($"r = {meta.Components}, share reached {meta.Pca!.ShareReached:F3}, stability ratio {meta.Truncation.StabilityRatio:F3}");
    }

    [Fact]
    public void Summary_IsEachComponentsPosteriorSummary()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        double[] probabilities = [0.1, 0.5, 0.9];
        BassBasisSummary summary = fit.Summary(probabilities);
        Assert.Equal(fit.Components.Count, summary.Components.Count);
        for (int k = 0; k < fit.Components.Count; k++)
        {
            BassPosteriorSummary own = fit.Components[k].Summary(probabilities).Posterior;
            foreach ((PosteriorSummary a, PosteriorSummary b) in new[] { (own.Count, summary.Components[k].Count), (own.Sigma2, summary.Components[k].Sigma2),
                                                                         (own.Tau!, summary.Components[k].Tau!), (own.Lambda, summary.Components[k].Lambda) })
            {
                Assert.Equal((a.Mean, a.Sd, a.Min, a.Max, a.Draws), (b.Mean, b.Sd, b.Min, b.Max, b.Draws));
                Assert.Equal(a.Quantiles, b.Quantiles);
                Assert.Equal(probabilities, b.Probabilities);
            }
        }
    }
}
