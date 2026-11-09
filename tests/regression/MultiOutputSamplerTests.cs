using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Samplers;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputSamplerTests(ITestOutputHelper output)
{
    private static BassConfig Config(TemperedForm form = TemperedForm.Integrated)
    {
        var config = MultiOutputTemperingReference.Config();
        return config with { Sampling = new BassSampling {
            Chains = 2, BurnIn = 1500, Samples = 12_000, Thin = 1, MasterSeed = 20261127,
            ReportEvidence = true,
            Tempering = new TemperingSpec { InverseTemperatures = [1, 0.4, 0.08, 0], ReferenceRung = true, Adaptive = false, Form = form }
        } };
    }

    [Theory]
    [InlineData(GPriorCentre.Mean)]
    [InlineData(GPriorCentre.Zero)]
    [InlineData(GPriorCentre.Flat)]
    public void MatrixTTemperatureSlopeAgreesWithFiniteDifferences(GPriorCentre centre)
    {
        var model = BassModel.Create(MultiOutputTemperingReference.Config(centre), MultiOutputTemperingReference.Data());
        var design = model.NewDesign();
        Assert.Null(design.Load([new BassTerm([new HingeFactor(0, 1, 2)])]));
        var e = model.Target.Evaluate(BassStructureFit.Of(design, model.Prior.LogPrior(design.Terms[0])), -2.0);
        foreach (double b in new[] { 0.03, 0.2, 0.6, 0.9 })
        {
            const double h = 1e-5;
            double finite = (e.LogLikelihood(b - 2 * h) - 8 * e.LogLikelihood(b - h) +
                             8 * e.LogLikelihood(b + h) - e.LogLikelihood(b + 2 * h)) / (12 * h);
            double slope = e.LogLikelihoodSlope(b);
            output.WriteLine($"{centre}, b={b}: slope {slope:G12}, difference {slope - finite:G3}");
            Assert.InRange(Math.Abs(slope - finite), 0, 2e-7);
        }
    }

    [Fact]
    public void SamplersValidateTheMatrixTPropernessBoundary()
    {
        var data = MultiOutputTemperingReference.Data();
        var flat = MultiOutputTemperingReference.Config(GPriorCentre.Flat, 1.5);
        BassConfig On(double hot) => flat with { Sampling = new BassSampling {
            Tempering = new TemperingSpec { InverseTemperatures = [1, hot], Adaptive = false }
        } };
        _ = BassModel.Create(On(0.125), data);
        Assert.Throws<ArgumentException>(() => BassModel.Create(On(0.0625), data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(On(0.03125), data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(flat with { Sampling = new BassSampling { Smc = new SmcSpec() } }, data));
        _ = BassModel.Create(Config(), data);
    }

    [Theory]
    [InlineData(TemperedForm.Integrated)]
    [InlineData(TemperedForm.Linear)]
    public void LadderRungsAndEvidenceMeetTheEnumeratedLaw(TemperedForm form)
    {
        var config = Config(form);
        var data = MultiOutputTemperingReference.Data();
        var model = BassModel.Create(config, data);
        double[] rungs = [1, 0.4, 0.08, 0];
        var exact = MultiOutputTemperingReference.Exact(model, 6, rungs);
        var series = BassTemperingTests.RungBatches(model, config.Sampling.Tempering!, rungs, ladders: 2, batches: 20, scans: 20_000);
        void Check(double[] means, double expected)
        {
            double mean = means.Average();
            double se = Math.Sqrt(means.Sum(x => Math.Pow(x - mean, 2)) / (means.Length - 1) / means.Length);
            output.WriteLine($"rung statistic {mean:G7}, exact {expected:G7}, SE {se:G3}");
            Assert.True(Math.Abs(mean - expected) < 5 * se + 0.003);
        }
        for (int r = 0; r < rungs.Length; r++)
        {
            Check(series.LogTau[r], exact[r].LogTauMean);
            for (int m = 0; m <= model.MaxBasis; m++)
                Check([.. Enumerable.Range(0, 40).Select(i => series.Count[r][m, i])], exact[r].Count[m]);
        }
        var fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1, RungWorkers = 1 });
        var evidence = Assert.IsType<BassEvidence>(fit.Evidence);
        output.WriteLine($"{form}: ln Z {evidence.LogEvidence:G9} ± {evidence.MonteCarloError:G3}, exact {exact[0].LogZ:G9}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact[0].LogZ) <= 5 * evidence.MonteCarloError + 0.015);
        Assert.NotNull(evidence.Priors.Outputs);
        Assert.All(fit.Draws.SelectMany(c => c), d => Assert.NotNull(d.Sigma));
    }

    [Theory]
    [InlineData("chain", false)]
    [InlineData("chain", true)]
    [InlineData("ladder", false)]
    [InlineData("ladder", true)]
    [InlineData("smc", false)]
    [InlineData("smc", true)]
    [InlineData("histogram", false)]
    [InlineData("histogram", true)]
    [InlineData("chain", false, "ridge")]
    [InlineData("ladder", true, "raw")]
    [InlineData("smc", false, "modified")]
    [InlineData("histogram", true, "ridge")]
    public void EverySamplerCarriesTheOutputDrawsAndResumesExactly(string mode, bool hierarchy, string form = "g")
    {
        var data = form == "g" ? MultiOutputTemperingReference.Data() : MultiOutputCoefficientTests.Data(2);
        var config = form == "g" ? MultiOutputTemperingReference.Config() : MultiOutputCoefficientTests.Config(form);
        if (hierarchy) config = config with { OutputCovariance = OutputCovariancePrior.HuangWand(2, [1, 0.7]) };
        var sampling = new BassSampling { Chains = 2, BurnIn = 80, Samples = 80, Thin = 1, MasterSeed = 20261203, ReportEvidence = true };
        if (mode == "ladder") sampling = sampling with { Tempering = Config().Sampling.Tempering };
        if (mode == "smc") sampling = sampling with { Smc = new SmcSpec {
            Particles = 120, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.8, TraceStages = true } };
        if (mode == "histogram") sampling = sampling with { FlatHistogram = new BassFlatHistogramSpec { Run = new FlatHistogramSpec {
            Chains = 2, LearningIterations = 120, ProductionDraws = 80, Window = new FlatHistogramWindow(-100, 0), Nodes = 8, LowSlope = 0
        } } };
        config = config with { Sampling = sampling };
        var options = new BassRunOptions { MaxDegreeOfParallelism = 1, RungWorkers = 1, Checkpoint = true };
        var whole = Bass.Fit(config, data, options);
        BassConfig cut = mode switch {
            "smc" => config with { Sampling = sampling with { Smc = sampling.Smc! with { MaxStages = 1 } } },
            "histogram" => config with { Sampling = sampling with { FlatHistogram = sampling.FlatHistogram! with { MaxIterations = 90 } } },
            _ => config with { Sampling = sampling with { Samples = 40 } }
        };
        var resumed = Bass.Fit(cut, data, options);
        if (mode is "smc" or "histogram")
        {
            int rounds = 0;
            while (resumed.Checkpoint is not null)
            {
                resumed = Bass.Resume(resumed.Checkpoint, cut, data, options);
                if (mode == "smc" && resumed.Diagnostics.Smc!.ReachedPosterior) break;
                Assert.True(++rounds < 100);
            }
        }
        else resumed = Bass.Resume(resumed.Checkpoint!, cut, data, options);
        Assert.Equal(whole.DrawWeights, resumed.DrawWeights);
        Assert.Equal(whole.Evidence?.LogEvidence, resumed.Evidence?.LogEvidence);
        Assert.Equal(whole.Draws.Count, resumed.Draws.Count);
        for (int c = 0; c < whole.Draws.Count; c++)
        {
            Assert.Equal(whole.Draws[c].Count, resumed.Draws[c].Count);
            for (int i = 0; i < whole.Draws[c].Count; i++)
            {
                var a = whole.Draws[c][i];
                var b = resumed.Draws[c][i];
                Assert.Equal(a.Tau, b.Tau);
                Assert.Equal(a.Structure.TermIds.Select(id => whole.Terms[id]), b.Structure.TermIds.Select(id => resumed.Terms[id]));
                Assert.Equal(a.Sigma!.Cast<double>(), b.Sigma!.Cast<double>());
                Assert.Equal(a.OutputCoefficients!.Cast<double>(), b.OutputCoefficients!.Cast<double>());
                Assert.Equal(a.OutputRates, b.OutputRates);
                Assert.Equal(hierarchy, a.OutputRates is not null);
            }
        }
        if (mode == "smc")
        {
            Assert.True(whole.Diagnostics.Smc!.ReachedPosterior);
            Assert.Equal(whole.Diagnostics.Smc.Overdispersed, resumed.Diagnostics.Smc!.Overdispersed);
            Assert.Equal(whole.Diagnostics.Smc.Replicates.Select(r => r.JumpedOverACollapse),
                         resumed.Diagnostics.Smc.Replicates.Select(r => r.JumpedOverACollapse));
        }
        output.WriteLine($"{mode}, hierarchical={hierarchy}: {whole.Draws.Sum(c => c.Count)} draws resumed exactly");
    }

    [Theory]
    [InlineData(TemperedForm.Integrated)]
    [InlineData(TemperedForm.Linear)]
    [InlineData(TemperedForm.Integrated, true)]
    [InlineData(TemperedForm.Linear, true)]
    public void SmcEvidenceAndFinalCountLawMatchTheEnumeration(TemperedForm form, bool hierarchy = false)
    {
        var config = (hierarchy ? HuangWandTests.Config() : MultiOutputTemperingReference.Config()) with { Sampling = new BassSampling {
            ReportEvidence = true, MasterSeed = 20261204, Smc = new SmcSpec {
                Particles = 1200, Replicates = 4, TargetConditionalEss = 0.9, MutationSteps = 3, Form = form, TraceStages = true }
        } };
        var data = MultiOutputTemperingReference.Data();
        var model = BassModel.Create(config, data);
        double logZ, count;
        if (hierarchy)
        {
            var reference = HuangWandTests.Integrate(model, 1);
            (logZ, count) = (reference.LogZ, reference.Count[1]);
        }
        else
        {
            var reference = MultiOutputTemperingReference.Exact(model, 6, [1])[0];
            (logZ, count) = (reference.LogZ, reference.Count[1]);
        }
        var fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        var evidence = fit.Evidence!;
        output.WriteLine($"{form}, hierarchical={hierarchy}: SMC ln Z={evidence.LogEvidence:G9} ± {evidence.MonteCarloError:G3}, exact={logZ:G9}");
        Assert.True(Math.Abs(evidence.LogEvidence - logZ) < 5 * evidence.MonteCarloError + 0.02);
        double probability = fit.Draws.SelectMany(c => c).Select((d, i) => d.Structure.TermIds.Length * fit.DrawWeights![i]).Sum();
        output.WriteLine($"SMC P(M=1)={probability:G8}, exact={count:G8}");
        Assert.InRange(Math.Abs(probability - count), 0, 0.04);
    }
}
