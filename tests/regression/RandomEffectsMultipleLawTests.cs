using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.RandomEffectsMultipleReference;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

public sealed class RandomEffectsMultipleLawTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData("crossed")]
    [InlineData("nested")]
    public void BothDrawnScales_AndStructures_MeetTheDenseEnumeration(string kind)
    {
        var (config, data) = Fixture(kind);
        BassModel model = BassModel.Create(config, data);
        Law coarse = Enumerate(model, 64), exact = Enumerate(model, 96);
        output.WriteLine($"{kind}: quadrature refinement Δln Z {exact.LogZ - coarse.LogZ:E3}, " +
            $"ΔE[ln s] {string.Join(", ", exact.LogScales.Zip(coarse.LogScales, (a,b) => (a-b).ToString("E3")))}");
        Assert.InRange(Math.Abs(exact.LogZ - coarse.LogZ), 0, 1e-6);
        for (int k = 0; k < 2; k++) Assert.InRange(Math.Abs(exact.LogScales[k] - coarse.LogScales[k]), 0, 0.001);
        BassTerm[][] structures = Structures(model);
        string[] keys = [.. structures.Select(BassPriorSamplerTests.Key)];
        const int batches = 60, size = 1000;
        double[][] observed = [.. Enumerable.Range(0, keys.Length + 2).Select(_ => new double[batches])];
        var chain = new BassChain(model, [], Math.Log(0.7), new Xoshiro256PlusPlus(210301), new Xoshiro256PlusPlus(210302));
        for (int j = 0; j < 5000; j++) chain.Step();
        for (int b = 0; b < batches; b++)
            for (int j = 0; j < size; j++)
            {
                chain.Step();
                int st = Array.IndexOf(keys, BassPriorSamplerTests.Key(chain.Current.Terms));
                Assert.True(st >= 0);
                observed[st][b] += 1.0 / size;
                for (int k = 0; k < 2; k++) observed[keys.Length + k][b] += Math.Log(chain.Current.GroupScales.Span[k]) / size;
            }
        double[] expected = [.. exact.Structures, .. exact.LogScales];
        for (int k = 0; k < expected.Length; k++)
        {
            var (mean, se) = Summary(observed[k]);
            output.WriteLine($"{kind}, {(k < keys.Length ? keys[k] : $"ln s{k-keys.Length}")}: {mean:G6} ± {se:G4}, dense {expected[k]:G6}");
            Assert.InRange(Math.Abs(mean - expected[k]), 0, 5 * se + 0.003);
        }
    }

    [Theory]
    [InlineData("crossed", 210311)]
    [InlineData("crossed", 210312)]
    [InlineData("nested", 210311)]
    [InlineData("nested", 210312)]
    public void Geweke_TwoScalesAndEffects_PreserveTheJointLaw(string kind, int seed)
    {
        if (!BassQualification.Enabled) return;
        const int pairs = 12_000;
        var (config, data) = Fixture(kind, tauDrawn: true);
        BassModel prior = BassModel.Create(config, data);
        var rng = new Xoshiro256PlusPlus(seed);
        var sampler = new BassPriorSampler(prior);
        string[] names = ["ln s0", "ln s1", "ln RMS(u0)", "ln RMS(u1)", "M", "ln τ", "ln σ²"];
        double[][] differences = [.. names.Select(_ => new double[pairs])];
        double[] Values(double[] scales, double[][] effects, int m, double logTau, double sigma2)
            => [.. scales.Select(s => Math.Log(s)), .. effects.Select(u => 0.5 * Math.Log(u.Sum(v => v * v) / u.Length)), m, logTau, Math.Log(sigma2)];
        for (int j = 0; j < pairs; j++)
        {
            BassState state = sampler.DrawPrior(rng).State;
            double[] scales = state.GroupScales.ToArray();
            var start = new BassChain(prior, [], state.LogTau, rng, rng);
            start.Structure.Adopt(state, sampler.LastDesign());
            start.Moves.SetState(state, prior.Target.Evaluate(state));
            BassSweepFit fit = start.Structure.DrawSweepFit(state, rng, prior.Target, power: 0);
            double[][] effects = PriorEffects(prior, scales, fit.Sigma2, rng);
            double[] y = [.. fit.Fitted];
            for (int i = 0; i < y.Length; i++)
            {
                for (int k = 0; k < 2; k++) y[i] += effects[k][prior.RandomEffects!.Layout(k).GroupOf(i)];
                y[i] += Math.Sqrt(fit.Sigma2) * Normal.SampleStandard(rng);
            }
            double[] before = Values(scales, effects, state.Count, state.LogTau, fit.Sigma2);
            BassConfig conditional = config with { RandomIntercepts = [.. config.RandomIntercepts!.Select((s,k) => s with { Scale = scales[k] })] };
            BassModel model = BassModel.Create(conditional, new BassData(data.Predictors, new ScalarResponse(y), data.Groupings));
            var chain = new BassChain(model, state.Terms, state.LogTau, rng, rng);
            for (int t = 0; t < 4; t++) chain.Step();
            BassDraw draw = chain.Compose(rng, new BassTermTable());
            double[] after = Values(draw.Groups!.Select(g => g.Scale).ToArray(), draw.Groups!.Select(g => g.Effects).ToArray(), chain.Current.Count, chain.Current.LogTau, draw.Sigma2);
            for (int k = 0; k < names.Length; k++) differences[k][j] = after[k] - before[k];
        }
        double worst = 0;
        for (int k = 0; k < names.Length; k++)
        {
            var (mean, se) = Summary(differences[k]);
            output.WriteLine($"{kind}, seed {seed}, {names[k]}: paired difference {mean:G5} ± {se:G4}, z={mean/se:F3}");
            worst = Math.Max(worst, Math.Abs(mean / se));
        }
        Assert.InRange(worst, 0, 4.5);
    }

    [Theory]
    [InlineData("crossed")]
    [InlineData("nested")]
    public void Smc_TwoDrawnScales_MeetTheDenseEvidence(string kind)
    {
        if (!BassQualification.Enabled) return;
        var (config, data) = Fixture(kind);
        BassModel model = BassModel.Create(config, data);
        double exact = Enumerate(model, 96).LogZ;
        config = config with { Sampling = config.Sampling with { MasterSeed = 210320, ReportEvidence = true,
            Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 } } };
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
        var evidence = fit.Evidence!;
        var readout = Assert.IsType<BassSmcReadout>(fit.Diagnostics.Smc);
        output.WriteLine($"{kind}: SMC ln Z {evidence.LogEvidence:F6} ± {evidence.MonteCarloError:F6}, full dense quadrature {exact:F6}; flags {readout.Flags.Count}");
        Assert.InRange(Math.Abs(evidence.LogEvidence - exact), 0, 4 * evidence.MonteCarloError + 0.01);
        Assert.Empty(readout.Flags);
    }
}
