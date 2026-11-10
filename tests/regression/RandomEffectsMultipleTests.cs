using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.RandomEffectsMultipleReference;

namespace Maths.Regression.Tests;

public sealed class RandomEffectsMultipleTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData("crossed")]
    [InlineData("nested")]
    [InlineData("alias")]
    public void JointTarget_AndEachScaleInformation_MatchDenseAlgebra(string kind)
    {
        var (config, data) = Fixture(kind);
        BassModel model = BassModel.Create(config, data);
        var blocks = model.RandomEffects!;
        double[] scales = blocks.Scales.ToArray();
        foreach (BassTerm[] terms in Structures(model))
        {
            BassState state = BassMixing.Loaded(model, terms, Math.Log(0.7), out _)!;
            double expected = DenseLogMarginal(model, terms, scales, 0.7);
            Assert.InRange(Math.Abs(model.Target.Evaluate(state).LogLikelihood(1) - expected), 0, 1e-11);
        }
        int n = data.Response.Cases;
        var factor = new DenseCholesky(n);
        factor.Decompose(Covariance(model, scales));
        const double step = 1e-5;
        for (int k = 0; k < scales.Length; k++)
        {
            double[] plus = (double[])scales.Clone(), minus = (double[])scales.Clone();
            plus[k] *= Math.Exp(step); minus[k] *= Math.Exp(-step);
            double[,] vp = Covariance(model, plus), vm = Covariance(model, minus);
            var solved = new double[n, n];
            for (int j = 0; j < n; j++)
            {
                double[] column = [.. Enumerable.Range(0, n).Select(i => (vp[i, j] - vm[i, j]) / (2 * step))];
                double[] x = factor.Solve(column);
                for (int i = 0; i < n; i++) solved[i, j] = x[i];
            }
            double information = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) information += 0.5 * solved[i, j] * solved[j, i];
            output.WriteLine($"{kind}, {blocks.Specifications[k].Grouping}: information {blocks.ScaleInformation(k):G12}, dense {information:G12}");
            Assert.InRange(Math.Abs(blocks.ScaleInformation(k) - information), 0, 1e-8);
        }
        if (kind == "alias") Assert.Equal(3, model.Resolved.AutoPicks.Count(p => p.Setting == "RandomEffectsAlias"));
    }

    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    public void EachScaleStep_PreservesTheOtherScale_AndReadsItsOwnPrior(int group)
    {
        var (config, data) = Fixture("crossed");
        BassModel model = BassModel.Create(config, data);
        var moves = new BassStructureMoves(model, [], Math.Log(0.7));
        var temperature = new Temperature(0.4, TemperingPath.Powered);
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(210305), temperature);
        RandomEffectsBlock before = BassGroupScaleBlock.Reading(chain).Whitening!;
        var step = new BassGroupScaleBlock(model, group);
        double t = Math.Log(1.7), width = step.Scale(before);
        var proposal = step.Propose(chain, before, width, t)!.Value;
        RandomEffectsBlock after = proposal.Candidate.Fit.Working!.Whitening!;
        Assert.Equal(before.Scales.Span[1-group], after.Scales.Span[1-group]);
        Assert.Equal(before.Scales.Span[group] * 1.7, after.Scales.Span[group], 12);
        BassGroupScalePrior prior = before.Specifications[group].ScalePrior!;
        double ratio = prior.LogDensity(after.Scales.Span[group]) - prior.LogDensity(before.Scales.Span[group]) + t;
        var old = model.Target.Evaluate(chain.Current);
        Assert.Equal(ratio, proposal.Evaluation.LogPrior - old.LogPrior, 12);
        double reverse = step.Scale(after);
        double target = 0.4 * (ratio + proposal.Evaluation.LogLikelihood(1) - old.LogLikelihood(1));
        Assert.Equal(target + Math.Log(width/reverse) + t*t*(1/(width*width) - 1/(reverse*reverse))/2, proposal.LogAcceptance, 12);
    }

    [Fact]
    public void HeldGrouping_HasNoScaleStep_AndTheSecondDrawnScaleRefusesHistogram()
    {
        var (config, data) = Fixture("crossed");
        config = config with { RandomIntercepts = [config.RandomIntercepts![0] with { ScalePrior = null }, config.RandomIntercepts[1]] };
        BassModel model = BassModel.Create(config, data);
        Assert.Equal(["group scale reader"], BassSweepSlots.Blocks(model, BassSweepSlot.Tau).Select(b => b.Name));
        config = config with { Sampling = config.Sampling with { FlatHistogram = new BassFlatHistogramSpec
            { Run = new FlatHistogramSpec { LearningIterations = 3000, ProductionDraws = 500, Nodes = 20 } } } };
        Assert.Contains("RandomIntercepts[1].ScalePrior", Assert.Throws<ArgumentException>(config.Validate).Message);
    }

    [Theory]
    [InlineData("crossed", false)]
    [InlineData("nested", false)]
    [InlineData("crossed", true)]
    [InlineData("nested", true)]
    public void NamedDraws_ResumeBitForBit_AndPredictionReadsEveryGrouping(string kind, bool tempered)
    {
        var (config, data) = Fixture(kind);
        config = config with { Sampling = config.Sampling with { Chains = 2, BurnIn = 160, Samples = 80,
            Tempering = tempered ? new TemperingSpec { InverseTemperatures = [1, 0.3, 0], ReferenceRung = true } : null } };
        BassConfig half = config with { Sampling = config.Sampling with { Samples = 40 } };
        BassFit full = Bass.Fit(config, data), first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        var equivalent = new BassData(data.Predictors, new ScalarResponse(data.Response.Values.ToArray()), data.Groupings);
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, equivalent);
        BassRunTests.AssertSameFit(full, resumed);
        BassModel model = BassModel.Create(config, data);
        BassDraw[] draws = [.. full.Draws.SelectMany(c => c)];
        foreach (BassDraw draw in draws)
        {
            Assert.Equal(2, draw.Groups!.Count);
            Assert.Equal(draw.Groups.Select(g => g.Scale), draw.WorkingOn(model).Whitening!.Scales.ToArray());
            foreach (int k in Enumerable.Range(0, 2))
            {
                var layout = model.RandomEffects!.Layout(k);
                int[] parent = config.RandomIntercepts![k].ParentGrouping is { } name
                    ? data.Groupings.Single(g => g.Name == name).Codes.ToArray() : new int[data.Response.Cases];
                foreach (int family in parent.Distinct())
                {
                    double sum = 0;
                    for (int i = 0; i < parent.Length; i++) if (parent[i] == family) sum += draw.Groups[k].Effects[layout.GroupOf(i)];
                    Assert.InRange(Math.Abs(sum), 0, 1e-12);
                }
            }
        }
        Assert.All(Enumerable.Range(0, 2), k => Assert.True(draws.Select(d => d.Groups![k].Scale).Distinct().Count() > 1));
        BassPrediction grouped = full.Predict(data, PredictionCases.Of(data));
        BassPrediction population = full.Predict(data, PredictionCases.Of(data.Predictors), new BassPredictOptions { Predictive = true });
        for (int d = 0; d < grouped.DrawCount; d++)
            for (int i = 0; i < data.Response.Cases; i++)
            {
                double effects = 0;
                for (int k = 0; k < 2; k++) effects += draws[d].Groups![k].Effects[model.RandomEffects!.Layout(k).GroupOf(i)];
                Assert.InRange(Math.Abs(grouped.Latent[d, i, 0] - population.Latent[d, i, 0] - effects), 0, 1e-12);
            }
        // Named retained records can arrive in a different order without changing the reading or fitted values.
        BassDraw reversed = draws[0] with { Groups = draws[0].Groups!.Reverse().ToArray() };
        var reader = new BassDrawFit(model, full.Terms);
        var a = new double[data.Response.Cases]; var b = new double[a.Length];
        reader.Evaluate(draws[0], a); reader.Evaluate(reversed, b);
        Assert.Equal(a, b);
        Assert.Equal(draws[0].WorkingOn(model).Values.ToArray(), reversed.WorkingOn(model).Values.ToArray());
        double normalizedNoise = 0;
        for (int d = 0; d < population.DrawCount; d++)
        {
            double variance = draws[d].Sigma2 * (1 + draws[d].Groups!.Sum(g => g.Scale*g.Scale));
            for (int i = 0; i < data.Response.Cases; i++)
                normalizedNoise += Math.Pow(population.Predictive![d,i,0] - population.Latent[d,i,0], 2) / variance;
        }
        Assert.InRange(normalizedNoise / (population.DrawCount * data.Response.Cases), 0.9, 1.1);
    }
}
