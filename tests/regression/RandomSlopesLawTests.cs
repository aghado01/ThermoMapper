using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.RandomSlopesReference;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

public sealed class RandomSlopesLawTests(ITestOutputHelper output)
{
    private sealed record Law(double LogZ, double[] Structures, double[] LogScales);

    private static Law Enumerate(BassModel model, int nodes)
    {
        BassTerm[][] structures = RandomEffectsMultipleReference.Structures(model);
        double[] prior = [.. structures.Select(t => model.Prior.LogCollectionWeight(t.Length)+t.Sum(model.Prior.LogPrior))];
        Block[] blocks = Blocks(model,[1,1]);
        double[][,] components = [.. blocks.Select(b => Covariance(b))];
        int n = model.Data.Response.Cases;
        foreach (double[,] component in components) for (int i = 0; i < n; i++) component[i,i] -= 1;
        GaussLegendre rule = GaussLegendre.Of(nodes);
        var mass = new double[structures.Length]; var moments = new double[2];
        double top = double.NegativeInfinity, total = 0;
        for (int a = 0; a < nodes; a++)
            for (int b = 0; b < nodes; b++)
            {
                int[] indices = [a,b];
                double[] scales = [.. indices.Select((index,k) => model.RandomEffects!.Specifications[k].ScalePrior!.Scale
                    *Math.Tan(0.25*Math.PI*(1+rule.Nodes[index])))];
                var v = new double[n,n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) v[i,j] = (i == j ? 1 : 0)+scales[0]*scales[0]*components[0][i,j]+scales[1]*scales[1]*components[1][i,j];
                double cell = Math.Log(rule.Weights[a]*rule.Weights[b]/4);
                for (int j = 0; j < structures.Length; j++)
                {
                    double log = prior[j]+cell+LogMarginal(model,structures[j],v,0.7);
                    if (log > top)
                    {
                        double rescale = Math.Exp(top-log); total *= rescale;
                        for (int k = 0; k < mass.Length; k++) mass[k] *= rescale;
                        for (int k = 0; k < moments.Length; k++) moments[k] *= rescale;
                        top = log;
                    }
                    double weight = Math.Exp(log-top); total += weight; mass[j] += weight;
                    for (int k = 0; k < 2; k++) moments[k] += weight*Math.Log(scales[k]);
                }
            }
        return new(top+Math.Log(total)-Math.Log(prior.Select(Math.Exp).Sum()),[.. mass.Select(v => v/total)],[.. moments.Select(v => v/total)]);
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void StructuresAndDrawnScales_MeetTheDenseEnumeration(RandomSlopeScale scaling)
    {
        var (config,data) = Fixture(scaling);
        BassModel model = BassModel.Create(config,data);
        Law coarse = Enumerate(model,64), exact = Enumerate(model,96);
        output.WriteLine($"{scaling}: refinement Δln Z {exact.LogZ-coarse.LogZ:E3}, ΔE[ln s] {string.Join(", ",exact.LogScales.Zip(coarse.LogScales,(a,b) => (a-b).ToString("E3")))}");
        Assert.InRange(Math.Abs(exact.LogZ-coarse.LogZ),0,1e-6);
        for (int k = 0; k < 2; k++) Assert.InRange(Math.Abs(exact.LogScales[k]-coarse.LogScales[k]),0,0.001);
        string[] keys = [.. RandomEffectsMultipleReference.Structures(model).Select(BassPriorSamplerTests.Key)];
        const int batches = 60, size = 800;
        double[][] observed = [.. Enumerable.Range(0,keys.Length+2).Select(_ => new double[batches])];
        var chain = new BassChain(model,[],Math.Log(0.7),new Xoshiro256PlusPlus(210710),new Xoshiro256PlusPlus(210711));
        for (int j = 0; j < 5000; j++) chain.Step();
        for (int b = 0; b < batches; b++)
            for (int j = 0; j < size; j++)
            {
                chain.Step();
                int st = Array.IndexOf(keys,BassPriorSamplerTests.Key(chain.Current.Terms));
                Assert.True(st >= 0); observed[st][b] += 1.0/size;
                for (int k = 0; k < 2; k++) observed[keys.Length+k][b] += Math.Log(chain.Current.GroupScales.Span[k])/size;
            }
        double[] expected = [.. exact.Structures,.. exact.LogScales];
        for (int k = 0; k < expected.Length; k++)
        {
            var (mean,se) = Summary(observed[k]);
            output.WriteLine($"{scaling}, {(k < keys.Length ? keys[k] : $"ln s{k-keys.Length}")}: {mean:G6} ± {se:G4}, dense {expected[k]:G6}");
            Assert.InRange(Math.Abs(mean-expected[k]),0,5*se+0.003);
        }
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit,210720)]
    [InlineData(RandomSlopeScale.Unit,210721)]
    [InlineData(RandomSlopeScale.Standardized,210720)]
    [InlineData(RandomSlopeScale.Standardized,210721)]
    public void Geweke_IndependentInterceptAndSlope_PreserveTheJointLaw(RandomSlopeScale scaling, int seed)
    {
        if (!BassQualification.Enabled) return;
        const int pairs = 10000;
        var (config,data) = Fixture(scaling);
        config = config with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Gamma(2,2)) };
        BassModel prior = BassModel.Create(config,data);
        Block[] reference = Blocks(prior,[1,1]);
        var rng = new Xoshiro256PlusPlus(seed); var sampler = new BassPriorSampler(prior);
        string[] names = ["ln s_a","ln s_b","ln RMS(a)","ln RMS(b)","M","ln τ","ln σ²"];
        double[][] differences = [.. names.Select(_ => new double[pairs])];
        double[] Values(double[] scales,double[][] effects,int m,double logTau,double sigma2)
            => [.. scales.Select(s => Math.Log(s)),.. effects.Select(u => 0.5*Math.Log(u.Sum(v => v*v)/u.Length)),m,logTau,Math.Log(sigma2)];
        for (int j = 0; j < pairs; j++)
        {
            BassState state = sampler.DrawPrior(rng).State;
            double[] scales = state.GroupScales.ToArray();
            var start = new BassChain(prior,[],state.LogTau,rng,rng);
            start.Structure.Adopt(state,sampler.LastDesign()); start.Moves.SetState(state,prior.Target.Evaluate(state));
            BassSweepFit fit = start.Structure.DrawSweepFit(state,rng,prior.Target,power:0);
            double[][] effects = [.. reference.Select((block,k) =>
            {
                int levels = block.Projector.GetLength(0);
                double[] normal = [.. Enumerable.Range(0,levels).Select(_ => Normal.SampleStandard(rng))];
                return Enumerable.Range(0,levels).Select(g => Math.Sqrt(fit.Sigma2)*scales[k]
                    *Enumerable.Range(0,levels).Sum(h => block.Projector[g,h]*normal[h])).ToArray();
            })];
            double[] y = [.. fit.Fitted];
            for (int i = 0; i < y.Length; i++)
            {
                for (int k = 0; k < reference.Length; k++)
                    for (int g = 0; g < effects[k].Length; g++) y[i] += reference[k].Design[i,g]*effects[k][g];
                y[i] += Math.Sqrt(fit.Sigma2)*Normal.SampleStandard(rng);
            }
            double[] before = Values(scales,effects,state.Count,state.LogTau,fit.Sigma2);
            BassConfig conditional = config with { RandomIntercepts = [config.RandomIntercepts![0] with { Scale = scales[0] }],
                RandomSlopes = [config.RandomSlopes![0] with { Scale = scales[1] }] };
            BassModel model = BassModel.Create(conditional,new BassData(data.Predictors,new ScalarResponse(y),data.Groupings));
            var chain = new BassChain(model,state.Terms,state.LogTau,rng,rng);
            for (int t = 0; t < 4; t++) chain.Step();
            BassDraw draw = chain.Compose(rng,new BassTermTable());
            double[] after = Values(draw.Groups!.Select(g => g.Scale).ToArray(),draw.Groups.Select(g => g.Effects).ToArray(),chain.Current.Count,chain.Current.LogTau,draw.Sigma2);
            for (int k = 0; k < names.Length; k++) differences[k][j] = after[k]-before[k];
        }
        for (int k = 0; k < names.Length; k++)
        {
            var (mean,se) = Summary(differences[k]);
            output.WriteLine($"{scaling}, {seed}, {names[k]}: paired difference {mean:G5} ± {se:G4}, z={mean/se:F3}");
            Assert.InRange(Math.Abs(mean/se),0,4.5);
        }
    }
}
