using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class HuangWandTests(ITestOutputHelper output)
{
    internal static BassConfig Config() => MultiOutputTemperingReference.Config() with {
        OutputCovariance = OutputCovariancePrior.HuangWand(2, [1.0, 0.7]),
        Coefficients = BassCoefficientPrior.G(TauPrior.Held(0.1))
    };

    [Fact]
    public void ConfigurationAndConditionalPriorReduceToTheScalarHierarchy()
    {
        var config = Config();
        Assert.Equal(config, BassConfig.FromJson(config.ToJson()));
        Assert.False(config.OutputCovariance!.ReadsResponse);
        Assert.True(OutputCovariancePrior.HuangWand().ReadsResponse);
        var data = MultiOutputCoefficientTests.Data(1);
        var one = BassModel.Create(config with { OutputCovariance = OutputCovariancePrior.HuangWand(2, [0.8]) }, data);
        var prior = one.Resolved.Outputs!.Rates!;
        foreach (double logRate in new[] { -25.0, -3, 0, 8 })
        {
            double r = Math.Exp(logRate);
            var state = new BassStructureMoves(one, [], Math.Log(0.1), logOutputRates: [logRate]).Start;
            var evaluation = one.Target.Evaluate(state);
            double expected = one.Target.OutputsMarginal(state.Fit, 0.1).LogMarginal(InverseWishart.Of(2, new double[,] { { 2 * r } }));
            Assert.Equal(expected, evaluation.LogMarginal, 11);
            double gamma = 0.5 * (Math.Log(1.0 / (2 * 0.8 * 0.8)) + logRate) - SpecialFunctions.LogGamma(0.5) - r / (2 * 0.8 * 0.8);
            Assert.Equal(gamma, prior.LogDensity([logRate]), 12);
            var scalar = BassModel.Create(config with { OutputCovariance = null, Sigma2 = Sigma2Prior.InverseGamma(1, r) },
                new BassData(data.Predictors, new ScalarResponse([.. data.Response.Values])));
            var design = scalar.NewDesign();
            Assert.InRange(Math.Abs(scalar.Target.Evaluate(BassStructureFit.Of(design, 0)).LogMarginal - expected), 0, 2e-10);
        }
    }

    [Fact]
    public void PriorDrawsHaveTheHalfTScalesAndUniformCorrelation()
    {
        var model = BassModel.Create(Config(), MultiOutputTemperingReference.Data());
        var rates = model.Resolved.Outputs!.Rates!;
        var rng = new Xoshiro256PlusPlus(20261130);
        var moments = new RunningMoments[4];
        const int count = 100_000;
        for (int i = 0; i < count; i++)
        {
            double[,] sigma = rates.PriorAt(rates.Sample(rng)).Sample(rng);
            double rho = sigma[1, 0] / Math.Sqrt(sigma[0, 0] * sigma[1, 1]);
            moments[0].Add(sigma[0, 0] <= 1 ? 1 : 0);
            moments[1].Add(sigma[1, 1] <= 0.49 ? 1 : 0);
            moments[2].Add(rho * rho);
            moments[3].Add(Math.Pow(rho, 4));
        }
        double[] exact = [1 / Math.Sqrt(3), 1 / Math.Sqrt(3), 1.0 / 3, 0.2];
        for (int i = 0; i < moments.Length; i++)
        {
            double z = (moments[i].Mean - exact[i]) / Math.Sqrt(moments[i].Variance / count);
            output.WriteLine($"prior statistic {i}: {moments[i].Mean:G7}, exact {exact[i]:G7}, z={z:F3}");
            Assert.InRange(Math.Abs(z), 0, 4.5);
        }
    }

    internal static (double LogZ, double[] Count, double[] Rates) Integrate(BassModel model, double width, double b = 1)
    {
        var states = BassTermSpace.F(model);
        var nodes = BassTargetTests.Nodes(-30, 8, width).ToArray();
        int n = nodes.Length, total = n * n;
        var logs = new double[states.Count];
        var rateMeans = new double[states.Count, 2];
        var values = new double[total];
        var design = model.NewDesign();
        for (int s = 0; s < states.Count; s++)
        {
            Assert.Null(design.Load(states[s].Terms));
            var c = MultiOutputCase.GForm(MultiOutputTemperingReference.Response(),
                MultiOutputCoreTests.Columns(model, design), GPriorCentre.Zero, 0.1, b);
            var collapse = MatrixTReference.Collapse(c);
            for (int a = 0; a < n; a++)
                for (int d = 0; d < n; d++)
                {
                    double r0 = Math.Exp(nodes[a].U), r1 = Math.Exp(nodes[d].U);
                    double logHyper = 0.5 * (nodes[a].U + nodes[d].U + Math.Log(0.5) + Math.Log(1.0 / 0.98))
                        - 2 * SpecialFunctions.LogGamma(0.5) - 0.5 * r0 - r1 / 0.98;
                    values[a * n + d] = Math.Log(nodes[a].W * nodes[d].W) + logHyper +
                        MatrixTReference.LogMarginal(collapse, new ReferenceSigmaPrior(3, new double[,] { { 2 * r0, 0 }, { 0, 2 * r1 } }));
                }
            double z = BassTemperingTests.LogSumExp(values);
            logs[s] = states[s].LogPrior + z;
            for (int a = 0; a < n; a++)
                for (int d = 0; d < n; d++)
                {
                    double mass = Math.Exp(values[a * n + d] - z);
                    rateMeans[s, 0] += mass * Math.Exp(nodes[a].U);
                    rateMeans[s, 1] += mass * Math.Exp(nodes[d].U);
                }
        }
        double logZ = BassTemperingTests.LogSumExp(logs);
        var countLaw = new double[model.MaxBasis + 1];
        var rateMean = new double[2];
        for (int s = 0; s < states.Count; s++)
        {
            double mass = Math.Exp(logs[s] - logZ);
            countLaw[states[s].Terms.Length] += mass;
            for (int l = 0; l < 2; l++) rateMean[l] += mass * rateMeans[s, l];
        }
        return (logZ - BassTemperingTests.LogSumExp(states.Select(s => s.LogPrior)), countLaw, rateMean);
    }

    [Fact]
    public void ChainMeetsTheEnumerationWithBothRatesIntegrated()
    {
        var model = BassModel.Create(Config(), MultiOutputTemperingReference.Data());
        var exact = Integrate(model, 2);
        var refined = Integrate(model, 1);
        Assert.InRange(Math.Abs(exact.LogZ - refined.LogZ), 0, 1e-7);
        for (int k = 0; k < 2; k++) Assert.InRange(Math.Abs(exact.Rates[k] - refined.Rates[k]), 0, 1e-7);
        int[] seeds = SeedTree.Derive(20261201, 2);
        var chain = new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(seeds[0]), new Xoshiro256PlusPlus(seeds[1]));
        for (int i = 0; i < 3000; i++) chain.Step();
        var series = new[] { new double[80_000], new double[80_000], new double[80_000] };
        for (int i = 0; i < series[0].Length; i++)
        {
            var state = chain.Step();
            series[0][i] = state.Count;
            series[1][i] = Math.Exp(state.LogOutputRates[0]);
            series[2][i] = Math.Exp(state.LogOutputRates[1]);
        }
        double[] expected = [refined.Count[1], refined.Rates[0], refined.Rates[1]];
        for (int j = 0; j < series.Length; j++)
        {
            double se = KhaosChainTests.BatchError(series[j], 50), mean = series[j].Average();
            output.WriteLine($"integrated-rate statistic {j}: {mean:G8}, exact {expected[j]:G8}, SE={se:G3}");
            Assert.True(Math.Abs(mean - expected[j]) < 5 * se);
        }
        output.WriteLine($"normalized ln Z={refined.LogZ:G12}, grid refinement={refined.LogZ - exact.LogZ:G3}");
    }
}
