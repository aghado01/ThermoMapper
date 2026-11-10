using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.RandomSlopesReference;

namespace Maths.Regression.Tests;

public sealed class RandomSlopesPredictionTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData(RandomSlopeScale.Unit,4.0,true)]
    [InlineData(RandomSlopeScale.Standardized,4.0,true)]
    [InlineData(RandomSlopeScale.Unit,8.0,true)]
    [InlineData(RandomSlopeScale.Standardized,8.0,true)]
    [InlineData(RandomSlopeScale.Unit,3.7,true)]
    [InlineData(RandomSlopeScale.Standardized,3.7,true)]
    [InlineData(RandomSlopeScale.Unit,4.0,false)]
    [InlineData(RandomSlopeScale.Standardized,4.0,false)]
    [InlineData(RandomSlopeScale.Unit,8.0,false)]
    [InlineData(RandomSlopeScale.Standardized,8.0,false)]
    [InlineData(RandomSlopeScale.Unit,3.7,false)]
    [InlineData(RandomSlopeScale.Standardized,3.7,false)]
    public void RescalingInput_AndItsDomain_PreservesTheFit(RandomSlopeScale scaling, double factor, bool drawn)
    {
        var (config,data) = Fixture(scaling,drawn:drawn);
        config = config with { Sampling = config.Sampling with { BurnIn = 200, Samples = 100 } };
        double[] x = ((ContinuousPredictor)data.Predictors[0]).Values.ToArray();
        var rescaled = new BassData([new ContinuousPredictor("x",[.. x.Select(v => factor*v)],new InputDomain(0,factor))],
            new ScalarResponse(data.Response.Values.ToArray()),data.Groupings);
        var traceA = new RandomSlopesRescaling.Trace(); var traceB = new RandomSlopesRescaling.Trace();
        BassFit a = Bass.Fit(config,data,new BassRunOptions { SweepObserver = traceA.Observe });
        BassFit b = Bass.Fit(config,rescaled,new BassRunOptions { SweepObserver = traceB.Observe });
        bool exact = factor is 4 or 8;
        var path = RandomSlopesRescaling.Compare(traceA,traceB,config.Sampling.BurnIn,config.Sampling.Samples,drawn,exact);
        double worstAbsolute = 0, worstRelative = path.Worst;
        void CompareVector(double[] left,double[] right)
        {
            worstRelative = Math.Max(worstRelative,RandomSlopesRescaling.Relative(left,right));
            foreach (var pair in left.Zip(right)) worstAbsolute = Math.Max(worstAbsolute,Math.Abs(pair.First-pair.Second));
            if (exact) RandomSlopesRescaling.Bits(left,right);
        }
        void ComparePositive(double left,double right)
        {
            worstRelative = Math.Max(worstRelative,RandomSlopesRescaling.Positive(left,right));
            worstAbsolute = Math.Max(worstAbsolute,Math.Abs(left-right));
            if (exact) RandomSlopesRescaling.Bits([left],[right]);
        }
        if (exact) BassRunTests.AssertSameFit(a,b);
        Assert.Equal(a.Terms.Count,b.Terms.Count);
        for (int id = 0; id < a.Terms.Count; id++) Assert.Equal(a.Terms[id],b.Terms[id]);
        for (int d = 0; d < a.Draws[0].Count; d++)
        {
            BassDraw left = a.Draws[0][d], right = b.Draws[0][d];
            Assert.Equal(left.Structure.TermIds,right.Structure.TermIds);
            ComparePositive(left.Sigma2,right.Sigma2); ComparePositive(left.Tau,right.Tau); ComparePositive(left.Lambda,right.Lambda);
            CompareVector(left.Coefficients,right.Coefficients);
            for (int k = 0; k < left.Groups!.Count; k++)
            {
                Assert.Equal(left.Groups[k].Grouping,right.Groups![k].Grouping); Assert.Equal(left.Groups[k].Input,right.Groups[k].Input);
                ComparePositive(left.Groups[k].Scale,right.Groups[k].Scale);
                CompareVector(left.Groups[k].Effects,right.Groups[k].Effects);
            }
        }
        output.WriteLine($"{scaling}, c={factor}, drawn={drawn}: {traceA.Sweeps.Count} sweeps, every decision/RNG state identical; "
            +$"rejected updates burn {path.BurnRejected}, sampling {path.SamplingRejected}; worst retained absolute {worstAbsolute:G17}, full relative {worstRelative:G17}");
        Assert.InRange(worstRelative,0,drawn ? 1e-9 : 1e-12);
    }

    private static void Near(double a, double b) => Assert.InRange(Math.Abs(a-b),0,1e-12);

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void NonexactRescaling_PreservesWhitening_AndFixedStateConditionalDraws(RandomSlopeScale scaling)
    {
        var (config,data) = Fixture(scaling);
        const double factor = 3.7;
        double[] x = ((ContinuousPredictor)data.Predictors[0]).Values.ToArray();
        var scaled = new BassData([new ContinuousPredictor("x",[.. x.Select(v => factor*v)],new InputDomain(0,factor))],
            new ScalarResponse(data.Response.Values.ToArray()),data.Groupings);
        BassModel a = BassModel.Create(config,data), b = BassModel.Create(config,scaled);
        double worst = RandomSlopesRescaling.Relative(a.Working.Values,b.Working.Values);
        foreach (BassTerm[] terms in RandomEffectsMultipleReference.Structures(a))
        {
            var left = new double[x.Length]; var right = new double[x.Length];
            if (terms.Length == 0) { Array.Fill(left,1); Array.Fill(right,1); }
            else { a.Basis.CaseColumn(terms[0],a.Inputs.Cases,left); b.Basis.CaseColumn(terms[0],b.Inputs.Cases,right); }
            a.RandomEffects!.Whiten(left); b.RandomEffects!.Whiten(right);
            worst = Math.Max(worst,RandomSlopesRescaling.Relative(left,right));
        }
        double[] marginalsA = [.. RandomEffectsMultipleReference.Structures(a).Select(terms => a.Target.Evaluate(BassMixing.Loaded(a,terms,Math.Log(0.7),out _)!).LogLikelihood(1))];
        double[] marginalsB = [.. RandomEffectsMultipleReference.Structures(b).Select(terms => b.Target.Evaluate(BassMixing.Loaded(b,terms,Math.Log(0.7),out _)!).LogLikelihood(1))];
        worst = Math.Max(worst,RandomSlopesRescaling.Relative(marginalsA,marginalsB));
        for (int i = 0; i < x.Length; i++)
        {
            var left = new double[x.Length]; left[i] = 1;
            var right = (double[])left.Clone();
            a.RandomEffects!.Whiten(left); b.RandomEffects!.Whiten(right);
            worst = Math.Max(worst,RandomSlopesRescaling.Relative(left,right));
        }
        var ra = new Xoshiro256PlusPlus(210704); var rb = new Xoshiro256PlusPlus(210704);
        for (int d = 0; d < 1000; d++)
        {
            double[][] left = a.RandomEffects!.DrawEffects(a,[],[0.2],1.3,ra);
            double[][] right = b.RandomEffects!.DrawEffects(b,[],[0.2],1.3,rb);
            for (int k = 0; k < left.Length; k++) worst = Math.Max(worst,RandomSlopesRescaling.Relative(left[k],right[k]));
        }
        output.WriteLine($"{scaling}, c=3.7 primitives: worst symmetric relative error {worst:G17}");
        Assert.InRange(worst,0,1e-12);
    }

    [Fact]
    public void RelativeComparison_IsHomogeneous_AndTreatsZerosUsingThePairedRms()
    {
        double[] a = [0,2e-20], b = [1e-33,2e-20];
        double relative = RandomSlopesRescaling.Relative(a,b);
        Assert.InRange(relative,0,1e-12);
        double changedUnits = RandomSlopesRescaling.Relative([.. a.Select(v => v*1e-150)],[.. b.Select(v => v*1e-150)]);
        Assert.InRange(Math.Abs(relative-changedUnits),0,1e-28);
        Assert.Equal(0,RandomSlopesRescaling.Relative([0,0],[0,0]));
        Assert.True(RandomSlopesRescaling.Relative(a,[1e-25,2e-20]) > 1e-9);
    }

    [Fact]
    public void ObservingEverySweep_PreservesTheUnobservedFit()
    {
        var (config,data) = Fixture(RandomSlopeScale.Standardized);
        config = config with { Sampling = config.Sampling with { BurnIn = 200, Samples = 100 } };
        var trace = new RandomSlopesRescaling.Trace();
        BassFit observed = Bass.Fit(config,data,new BassRunOptions { SweepObserver = trace.Observe });
        BassFit plain = Bass.Fit(config,data);
        Assert.Equal(300,trace.Sweeps.Count);
        BassRunTests.AssertSameFit(plain,observed);
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit,false)]
    [InlineData(RandomSlopeScale.Standardized,false)]
    [InlineData(RandomSlopeScale.Unit,true)]
    [InlineData(RandomSlopeScale.Standardized,true)]
    public void NamedSlopeDraws_ResumeBitForBit(RandomSlopeScale scaling, bool tempered)
    {
        var (config,data) = Fixture(scaling);
        config = config with { Sampling = config.Sampling with { Chains = 2, BurnIn = 100, Samples = 40,
            Tempering = tempered ? new TemperingSpec { InverseTemperatures = [1,0.3,0], ReferenceRung = true } : null } };
        BassConfig half = config with { Sampling = config.Sampling with { Samples = 20 } };
        BassFit full = Bass.Fit(config,data), first = Bass.Fit(half,data,new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!,half,data);
        BassRunTests.AssertSameFit(full,resumed);
        BassModel model = BassModel.Create(config,data);
        foreach (BassDraw draw in resumed.Draws.SelectMany(c => c))
        {
            Assert.Null(draw.Groups![0].Input); Assert.Equal("x",draw.Groups[1].Input);
            RandomSlopesTests.CheckConstraints(model,draw.Groups[1].Effects,1);
        }
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void Prediction_Pointwise_AndMeanSurface_ReadTheDeclaredInputAndBlockIdentity(RandomSlopeScale scaling)
    {
        var (config,data) = Fixture(scaling);
        config = config with { Sampling = config.Sampling with { BurnIn = 100, Samples = 60 } };
        BassFit fit = Bass.Fit(config,data);
        BassModel model = BassModel.Create(config,data);
        double[] x = ((ContinuousPredictor)data.Predictors[0]).Values.ToArray();
        CategoricalPredictor grouping = data.Groupings[0];
        string[] labels = [.. grouping.Codes.ToArray().Select(code => grouping.Levels[code])];
        var cases = new PredictionCases([new ContinuousColumn("x",x),new CategoricalColumn("g",labels)]);
        var options = new BassPredictOptions { MeanSurface = true, Predictive = false };
        BassPrediction prediction = fit.Predict(data,cases,options);
        var reader = new BassDrawFit(model,fit.Terms);
        var surface = new BassSurface(model.Basis,fit.Terms,model.Inputs.Cases,null);
        var means = new BassConditionalMeans(model,fit.Terms);
        var expectedSurface = new double[x.Length];
        double[] z = Input(data,config.RandomSlopes![0]);
        for (int d = 0; d < fit.Draws[0].Count; d++)
        {
            BassDraw draw = fit.Draws[0][d];
            var actual = new double[x.Length]; reader.Evaluate(draw,actual);
            for (int i = 0; i < x.Length; i++) Near(actual[i],prediction.Latent[d,i,0]);
            double[] a = means.Coefficients(draw), f = new double[x.Length];
            surface.Evaluate(draw.Structure.TermIds,a,f);
            double[] residual = [.. data.Response.Values.ToArray().Zip(f,(y,value) => y-value)];
            double[][] u = Conditional(Blocks(model,[.. draw.Groups!.Select(g => g.Scale)]),residual,draw.Sigma2).Mean;
            for (int i = 0; i < x.Length; i++) expectedSurface[i] += (f[i]+u[0][labelsIndex(i)]+z[i]*u[1][labelsIndex(i)])/fit.Draws[0].Count;
        }
        int labelsIndex(int i) => grouping.Codes.Span[i]-1;
        for (int i = 0; i < x.Length; i++) Near(expectedSurface[i],prediction.MeanSurface![i,0]);
        var reversed = fit.Draws[0].Select(d => d with { Groups = d.Groups!.Reverse().ToArray() }).ToArray();
        BassPrediction reordered = BassDirectPrediction.Predict(model,fit.Terms,[reversed],null,cases,options);
        Assert.Equal(prediction.Latent,reordered.Latent); Assert.Equal(prediction.MeanSurface,reordered.MeanSurface);
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void PopulationPrediction_UsesSlopeVarianceAtTheFittedTrainingScale(RandomSlopeScale scaling)
    {
        var (config,data) = Fixture(scaling);
        config = config with { Sampling = config.Sampling with { BurnIn = 100, Samples = 40 } };
        BassFit fit = Bass.Fit(config,data);
        double[] x = [0,0.5,1];
        var cases = new PredictionCases([new ContinuousColumn("x",x)]);
        const int seed = 210705;
        BassPrediction prediction = fit.Predict(data,cases,new BassPredictOptions { Predictive = true, Seed = seed });
        double[] unit = ((ContinuousPredictor)data.Predictors[0]).Values.ToArray();
        double mean = scaling == RandomSlopeScale.Unit ? 0 : unit.Average();
        double sd = scaling == RandomSlopeScale.Unit ? 1 : Math.Sqrt(unit.Select(v => (v-mean)*(v-mean)).Average());
        var rng = new Xoshiro256PlusPlus(seed);
        for (int d = 0; d < fit.Draws[0].Count; d++)
        {
            BassDraw draw = fit.Draws[0][d]; double sigma = Math.Sqrt(draw.Sigma2);
            double[] noise = [.. x.Select(_ => sigma*Normal.SampleStandard(rng))];
            foreach (BassGroupDraw group in draw.Groups!)
                for (int i = 0; i < x.Length; i++)
                    noise[i] += sigma*group.Scale*Normal.SampleStandard(rng)*(group.Input is null ? 1 : (x[i]-mean)/sd);
            for (int i = 0; i < x.Length; i++) Near(prediction.Latent[d,i,0]+noise[i],prediction.Predictive![d,i,0]);
        }
        var unknown = new PredictionCases([new ContinuousColumn("x",x),new CategoricalColumn("g",["new","new","new"])]);
        BassPrediction outside = fit.Predict(data,unknown,new BassPredictOptions { Predictive = true, Seed = seed, UnseenLevels = UnseenLevel.OutsideEverySubset });
        Assert.Equal(prediction.Predictive,outside.Predictive);
        Assert.Throws<ArgumentException>(() => fit.Predict(data,unknown));
    }
}
