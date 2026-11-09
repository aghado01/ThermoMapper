using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputSamplerQualificationTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData("ladder")]
    [InlineData("smc")]
    public void ParallelOutputPencilsAndRatesGiveIdenticalDraws(string mode)
    {
        var config = MultiOutputCoefficientTests.Config("raw") with { OutputCovariance = OutputCovariancePrior.HuangWand(2, [1, 0.7]) };
        var data = MultiOutputCoefficientTests.Data(2);
        var sampling = new BassSampling { Chains = 1, BurnIn = 100, Samples = 100, Thin = 1, MasterSeed = 20261209 };
        sampling = mode == "ladder" ? sampling with { Tempering = new TemperingSpec {
            InverseTemperatures = [1, 0.3, 0.05, 0], ReferenceRung = true, Adaptive = false } }
            : sampling with { Smc = new SmcSpec { Particles = 120, Replicates = 1, MutationSteps = 2, TargetConditionalEss = 0.8 } };
        config = config with { Sampling = sampling };
        var one = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1, RungWorkers = 1 });
        var two = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 2, RungWorkers = 2 });
        Assert.Equal(one.DrawWeights, two.DrawWeights);
        for (int c = 0; c < one.Draws.Count; c++)
            for (int i = 0; i < one.Draws[c].Count; i++)
            {
                var a = one.Draws[c][i];
                var b = two.Draws[c][i];
                Assert.Equal(a.Tau, b.Tau);
                Assert.Equal(a.Structure.TermIds.Select(id => one.Terms[id]), b.Structure.TermIds.Select(id => two.Terms[id]));
                Assert.Equal(a.Sigma!.Cast<double>(), b.Sigma!.Cast<double>());
                Assert.Equal(a.OutputCoefficients!.Cast<double>(), b.OutputCoefficients!.Cast<double>());
                Assert.Equal(a.OutputRates, b.OutputRates);
            }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void HistogramReweightedLawAndEvidenceMeetTheEnumeration(bool hierarchy)
    {
        if (!BassQualification.Enabled) return;
        var data = MultiOutputTemperingReference.Data();
        var config = hierarchy ? HuangWandTests.Config() : MultiOutputTemperingReference.Config();
        var model = BassModel.Create(config, data);
        double logZ, probability;
        if (hierarchy)
        {
            var exact = HuangWandTests.Integrate(model, 1);
            (logZ, probability) = (exact.LogZ, exact.Count[1]);
        }
        else
        {
            var exact = MultiOutputTemperingReference.Exact(model, 6, [1])[0];
            (logZ, probability) = (exact.LogZ, exact.Count[1]);
        }
        var rng = new Xoshiro256PlusPlus(20261205);
        var sampler = new BassPriorSampler(model);
        double[] levels = [.. Enumerable.Range(0, 3000).Select(_ => sampler.DrawState(rng).Evaluation.LogMarginal).Order()];
        var shape = MulticanonicalWeights.Flat(levels[30], levels[^1] + 1, 0, 16);
        config = config with { Sampling = new BassSampling { ReportEvidence = true, MasterSeed = 20261206, Thin = 1,
            FlatHistogram = new BassFlatHistogramSpec { Run = BassFlatHistogramTests.DeclaredRun(shape, 15_000, 20_000, chains: 4) } } };
        var fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        var evidence = fit.Evidence!;
        double measured = fit.Draws.SelectMany(c => c).Select((d, i) => d.Structure.TermIds.Length * fit.DrawWeights![i]).Sum();
        output.WriteLine($"Huang-Wand={hierarchy}: histogram P(M=1)={measured:G8}, exact={probability:G8}; " +
                         $"ln Z={evidence.LogEvidence:G10} ± {evidence.MonteCarloError:G3}, exact={logZ:G10}");
        Assert.InRange(Math.Abs(measured - probability), 0, 0.025);
        Assert.True(Math.Abs(evidence.LogEvidence - logZ) < 5 * evidence.MonteCarloError + 0.03);
    }

    [Fact]
    public void HuangWandRungsAndEvidenceMeetTheIntegratedRateReference()
    {
        if (!BassQualification.Enabled) return;
        var config = HuangWandTests.Config();
        var data = MultiOutputTemperingReference.Data();
        var model = BassModel.Create(config, data);
        double[] rungs = [1, 0.3, 0];
        var cold = HuangWandTests.Integrate(model, 1);
        var warm = HuangWandTests.Integrate(model, 1, 0.3);
        // The exact reference rung is the prior normalized on F.
        var states = BassTermSpace.F(model);
        double priorMass = BassTemperingTests.LogSumExp(states.Select(s => s.LogPrior));
        double[] counts = [cold.Count[1], warm.Count[1],
            states.Where(s => s.Terms.Length == 1).Sum(s => Math.Exp(s.LogPrior - priorMass))];
        double[][] rates = [cold.Rates, warm.Rates, [1.0, 0.49]];
        var spec = new TemperingSpec { InverseTemperatures = rungs, ReferenceRung = true, Adaptive = false };
        var seeds = SeedTree.Derive(20261207, 2);
        var ladder = new BassLadder(model, spec, [], Math.Log(0.1), seeds[0]);
        for (int i = 0; i < 2000; i++) ladder.Step();
        const int steps = 40_000;
        var series = rungs.Select(_ => new[] { new double[steps], new double[steps], new double[steps] }).ToArray();
        for (int i = 0; i < steps; i++)
        {
            ladder.Step();
            for (int r = 0; r < rungs.Length; r++)
            {
                var state = ladder.Tempering.Replicas[r].Current;
                series[r][0][i] = state.Count;
                series[r][1][i] = Math.Exp(state.LogOutputRates[0]);
                series[r][2][i] = Math.Exp(state.LogOutputRates[1]);
            }
        }
        for (int r = 0; r < rungs.Length; r++)
            for (int j = 0; j < 3; j++)
            {
                double exact = j == 0 ? counts[r] : rates[r][j - 1];
                double mean = series[r][j].Average(), se = KhaosChainTests.BatchError(series[r][j], 40);
                output.WriteLine($"b={rungs[r]}, statistic={j}: {mean:G8}, exact={exact:G8}, SE={se:G3}");
                Assert.True(Math.Abs(mean - exact) < 5 * se + 0.002);
            }
        config = config with { Sampling = new BassSampling { Chains = 2, BurnIn = 2000, Samples = 20_000, Thin = 1,
            MasterSeed = seeds[1], ReportEvidence = true, Tempering = spec } };
        var fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1, RungWorkers = 1 });
        output.WriteLine($"Huang-Wand ladder ln Z={fit.Evidence!.LogEvidence:G10} ± {fit.Evidence.MonteCarloError:G3}, exact={cold.LogZ:G10}");
        Assert.True(Math.Abs(fit.Evidence.LogEvidence - cold.LogZ) < 5 * fit.Evidence.MonteCarloError + 0.02);
    }
}
