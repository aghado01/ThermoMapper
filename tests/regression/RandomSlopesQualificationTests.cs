using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Tests;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.RandomSlopesReference;

namespace Maths.Regression.Tests;

public sealed class RandomSlopesQualificationTests(ITestOutputHelper output)
{
    private const double Sigma = 0.35;
    private sealed record Generated(BassData Data,double[] Effects,double[] Population,double[] TestMean,
        PredictionCases Cases,PredictionCases Grouped);
    internal sealed record Reading(double Coverage,double Z2,int Levels,bool ScaleCovered,bool SigmaCovered,
        double LevelZ,double SlopeZ,double GroupError,double PopulationError,double PlainPopulationError);

    private static Generated Generate(RandomSlopeScale scaling,bool skewed,int seed)
    {
        const int n = 280, levels = 8, held = 240;
        var rng = new Xoshiro256PlusPlus(seed);
        int[] codes = [.. Enumerable.Range(0,n).Select(i => i < 4*levels ? i%levels : Math.Min(levels-1,(int)(levels*Math.Pow(rng.NextDouble(),1.4))))];
        double[] x = [.. Enumerable.Range(0,n).Select(_ => rng.NextDouble())];
        var grouping = new CategoricalPredictor("g",codes,[.. Enumerable.Range(0,levels).Select(g => $"g{g}")]);
        BassData prototype = new([new ContinuousPredictor("x",x,new InputDomain(0,1))],new ScalarResponse(new double[n]),[grouping]);
        var slope = new BassRandomSlope("g","x",scaling,1);
        double[] z = Input(prototype,slope);
        Block block = Slope(codes,z,1);
        double[] raw = [.. Enumerable.Range(0,levels).Select(_ => skewed ? -0.85*Math.Log(1-rng.NextDouble())-0.85 : rng.NextDouble() < 0.5 ? -0.9 : 0.9)];
        double[] effects = [.. Enumerable.Range(0,levels).Select(g => Enumerable.Range(0,levels).Sum(h => block.Projector[g,h]*raw[h]))];
        double Truth(double value) => 1.1+1.3*value;
        double[] y = [.. x.Select((v,i) => Truth(v)+z[i]*effects[codes[i]]+Sigma*Normal.SampleStandard(rng))];
        double[] test = [.. Enumerable.Range(0,held).Select(_ => rng.NextDouble())];
        int[] testGroups = [.. Enumerable.Range(0,held).Select(_ => (int)(levels*rng.NextDouble()))];
        double mean = scaling == RandomSlopeScale.Unit ? 0 : x.Average();
        double deviation = scaling == RandomSlopeScale.Unit ? 1 : Math.Sqrt(x.Select(v => (v-mean)*(v-mean)).Average());
        var cases = new PredictionCases([new ContinuousColumn("x",test)]);
        var grouped = new PredictionCases([.. cases.Columns,new CategoricalColumn("g",[.. testGroups.Select(g => $"g{g}")])]);
        return new(new BassData(prototype.Predictors,new ScalarResponse(y),prototype.Groupings),effects,[.. test.Select(Truth)],
            [.. test.Select((v,i) => Truth(v)+(v-mean)/deviation*effects[testGroups[i]])],cases,grouped);
    }

    private static (double Low,double High) Interval(double[] values)
    { Array.Sort(values); return (PosteriorInterval.Quantile(values,0.05),PosteriorInterval.Quantile(values,0.95)); }

    private static double Error(BassFit fit,BassData data,PredictionCases cases,double[] truth)
    {
        BassPrediction prediction = fit.Predict(data,cases,new BassPredictOptions { Predictive = false, MeanSurface = true });
        return Math.Sqrt(Enumerable.Range(0,truth.Length).Average(i =>
            (prediction.MeanSurface![i,0]-truth[i])*(prediction.MeanSurface[i,0]-truth[i])));
    }

    internal static Reading Recovery(RandomSlopeScale scaling,bool skewed,int dataSeed,int chainSeed)
    {
        Generated g = Generate(scaling,skewed,dataSeed);
        BassConfig config = new()
        {
            Structure = new BassStructureSpec { MaxBasis = 1, MaxInteraction = 1, Legendre = new LegendreBasis { MaxDegree = 1 } },
            // Hold coefficient shrinkage to isolate the slope convention; Geweke separately draws τ with both effect scales.
            Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(0.7)), Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.InverseGamma(2,0.1), RandomSlopes = [new("g","x",scaling,1,new())],
            Sampling = new BassSampling { Chains = 2, BurnIn = 1000, Samples = 600, MasterSeed = chainSeed },
        };
        BassFit fit = Bass.Fit(config,g.Data,new BassRunOptions { MaxDegreeOfParallelism = 1 });
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        BassModel model = BassModel.Create(config,g.Data);
        int covered = 0; double z2 = 0;
        for (int level = 0; level < g.Effects.Length; level++)
        {
            int index = level;
            double[] values = [.. draws.Select(d => d.Groups![0].Effects[index])];
            double mean = values.Average(), variance = values.Sum(v => (v-mean)*(v-mean))/(values.Length-1);
            z2 += (mean-g.Effects[level])*(mean-g.Effects[level])/variance;
            var interval = Interval(values);
            if (interval.Low <= g.Effects[level] && g.Effects[level] <= interval.High) covered++;
        }
        foreach (BassDraw draw in draws) RandomSlopesTests.CheckConstraints(model,draw.Groups![0].Effects,0);
        var scale = Interval([.. draws.Select(d => d.Groups![0].Scale)]);
        double target = Math.Sqrt(g.Effects.Sum(u => u*u)/model.RandomEffects!.Layout(0).Dimension)/Sigma;
        var sigma = Interval([.. draws.Select(d => Math.Sqrt(d.Sigma2))]);
        double[] x = ((ContinuousPredictor)g.Data.Predictors[0]).Values.ToArray();
        double[] truth = [.. x.Select(v => 1.1+1.3*v)];
        double[] z = Input(g.Data,config.RandomSlopes![0]);
        double meanZ = z.Average(), energy = z.Sum(v => (v-meanZ)*(v-meanZ));
        double targetLevel = truth.Average(), targetSlope = z.Zip(truth,(v,f) => (v-meanZ)*f).Sum()/energy;
        BassPrediction population = fit.Predict(g.Data,PredictionCases.Of(g.Data.Predictors),new BassPredictOptions { Predictive = false });
        var levels = new double[population.DrawCount]; var slopes = new double[population.DrawCount];
        for (int d = 0; d < population.DrawCount; d++)
        {
            levels[d] = Enumerable.Range(0,x.Length).Average(i => population.Latent[d,i,0]);
            slopes[d] = Enumerable.Range(0,x.Length).Sum(i => (z[i]-meanZ)*population.Latent[d,i,0])/energy;
        }
        double StandardizedError(double[] values,double expected)
        {
            double mean = values.Average(), sd = Math.Sqrt(values.Sum(v => (v-mean)*(v-mean))/(values.Length-1));
            return (mean-expected)/sd;
        }
        BassData bare = new(g.Data.Predictors,g.Data.Response);
        BassFit plain = Bass.Fit(config with { RandomSlopes = null },bare,new BassRunOptions { MaxDegreeOfParallelism = 1 });
        return new((double)covered/g.Effects.Length,z2/g.Effects.Length,g.Effects.Length,
            scale.Low <= target && target <= scale.High,sigma.Low <= Sigma && Sigma <= sigma.High,
            StandardizedError(levels,targetLevel),StandardizedError(slopes,targetSlope),Error(fit,g.Data,g.Grouped,g.TestMean),
            Error(fit,g.Data,g.Cases,g.Population),Error(plain,bare,g.Cases,g.Population));
    }

    internal static SweptClaim[] Claims => [.. new[] { RandomSlopeScale.Unit,RandomSlopeScale.Standardized }.SelectMany(scaling =>
        SeedTree.Derive(210730,4).Zip(SeedTree.Derive(210731,4)).Select((seeds,index) => new SweptClaim(
            $"RandomSlopesQualificationTests.NonGaussianSlopes({scaling},{index})",seeds.First,[seeds.Second],(d,c) =>
            {
                Reading r = Recovery(scaling,index%2 == 1,d,c[0]);
                return [new Comparison("population RMSE: block ≤ 1.1 × no block",r.PopulationError,1.1*r.PlainPopulationError)];
            })))];

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void NonGaussianSlopes_AreRecovered_WithUnknownScale(RandomSlopeScale scaling)
    {
        if (!BassQualification.Enabled) return;
        int[] dataSeeds = SeedTree.Derive(210730,4), chainSeeds = SeedTree.Derive(210731,4);
        var readings = new List<Reading>();
        for (int j = 0; j < dataSeeds.Length; j++)
        {
            Reading r = Recovery(scaling,j%2 == 1,dataSeeds[j],chainSeeds[j]); readings.Add(r);
            output.WriteLine($"{scaling}, dataset {j}: coverage {r.Coverage:F3}, z² {r.Z2:F3}, scale covered {r.ScaleCovered}, σ covered {r.SigmaCovered}; level z {r.LevelZ:F3}, slope z {r.SlopeZ:F3}; held-out RMSE {r.GroupError:F4}; population {r.PopulationError:F4}, no block {r.PlainPopulationError:F4}");
            Assert.InRange(Math.Abs(r.LevelZ),0,2); Assert.InRange(Math.Abs(r.SlopeZ),0,2);
            Assert.InRange(r.GroupError,0,Sigma);
            Assert.InRange(r.PopulationError,0,1.1*r.PlainPopulationError);
        }
        Assert.InRange(readings.Average(r => r.Coverage),0.75,1);
        Assert.InRange(readings.Average(r => r.Z2),0.15,2.5);
        Assert.InRange(readings.Count(r => r.ScaleCovered),3,4);
        Assert.InRange(readings.Count(r => r.SigmaCovered),2,4);
    }

    [Fact]
    public void RecoveryComparisons_AcrossSeeds()
    {
        if (SeedSweep.Seeds is int seeds) SeedSweep.Run("RandomSlopes",Claims,seeds,output);
    }
}
