using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using static Maths.Regression.Tests.RandomSlopesReference;

namespace Maths.Regression.Tests;

public sealed class RandomSlopesTests
{
    [Theory]
    [InlineData(RandomSlopeScale.Unit,false)]
    [InlineData(RandomSlopeScale.Unit,true)]
    [InlineData(RandomSlopeScale.Standardized,false)]
    [InlineData(RandomSlopeScale.Standardized,true)]
    public void Whitening_Marginal_AndConditionalMean_MatchDenseAlgebra(RandomSlopeScale scaling, bool intercept)
    {
        var (config,data) = Fixture(scaling,intercept);
        BassModel model = BassModel.Create(config,data);
        RandomEffectsBlock block = model.RandomEffects!;
        Block[] reference = Blocks(model);
        double[,] v = Covariance(reference);
        int n = data.Response.Cases;
        var factor = new DenseCholesky(n); factor.Decompose(v);
        Assert.InRange(Math.Abs(block.LogDeterminant-factor.LogDet),0,1e-12);
        var white = new double[n][];
        for (int i = 0; i < n; i++) { white[i] = new double[n]; white[i][i] = 1; block.Whiten(white[i]); }
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
            {
                double actual = 0;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) actual += white[a][i]*v[i,j]*white[b][j];
                Assert.InRange(Math.Abs(actual-(a == b ? 1 : 0)),0,1e-12);
            }
        foreach (BassTerm[] terms in RandomEffectsMultipleReference.Structures(model))
        {
            BassState state = BassMixing.Loaded(model,terms,Math.Log(0.7),out _)!;
            Assert.InRange(Math.Abs(model.Target.Evaluate(state).LogLikelihood(1)-LogMarginal(model,terms,reference,0.7)),0,1e-11);
        }
        double[] residual = [.. data.Response.Values.ToArray().Select(y => y-0.2)];
        var expected = Conditional(reference,residual,1.3);
        double[][] means = block.MeanEffects(model,[],[0.2]);
        for (int k = 0; k < means.Length; k++)
            for (int g = 0; g < means[k].Length; g++) Assert.InRange(Math.Abs(means[k][g]-expected.Mean[k][g]),0,1e-12);
        var constant = Enumerable.Repeat(1.0,n).ToArray(); block.Whiten(constant);
        foreach (double value in constant) Assert.InRange(Math.Abs(value-1),0,1e-12);
        if (!intercept)
        {
            double[] z = Input(data,config.RandomSlopes![0]), before = (double[])z.Clone(); block.Whiten(z);
            for (int i = 0; i < n; i++) Assert.InRange(Math.Abs(z[i]-before[i]),0,1e-12);
        }
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void ConditionalDraws_MatchTheDenseJointCovariance_AndBothConstraints(RandomSlopeScale scaling)
    {
        var (config,data) = Fixture(scaling);
        BassModel model = BassModel.Create(config,data);
        var expected = Conditional(Blocks(model),data.Response.Values.ToArray(),1.3);
        int dimension = expected.Covariance.GetLength(0);
        var sum = new double[dimension]; var products = new double[dimension,dimension];
        var rng = new Xoshiro256PlusPlus(210702);
        const int draws = 12000;
        for (int d = 0; d < draws; d++)
        {
            double[][] effects = model.RandomEffects!.DrawEffects(model,[],[0],1.3,rng);
            CheckConstraints(model,effects[1],1);
            double[] centred = [.. effects.SelectMany((u,k) => u.Select((v,g) => v-expected.Mean[k][g]))];
            for (int a = 0; a < dimension; a++)
            {
                sum[a] += centred[a]/draws;
                for (int b = 0; b < dimension; b++) products[a,b] += centred[a]*centred[b]/draws;
            }
        }
        for (int a = 0; a < dimension; a++)
        {
            Assert.InRange(Math.Abs(sum[a]),0,5*Math.Sqrt(expected.Covariance[a,a]/draws));
            for (int b = 0; b < dimension; b++)
                Assert.InRange(Math.Abs(products[a,b]-expected.Covariance[a,b]),0,
                    0.05*Math.Sqrt(expected.Covariance[a,a]*expected.Covariance[b,b])+0.001);
        }
    }

    internal static void CheckConstraints(BassModel model, double[] effects, int block)
    {
        RandomEffectGrouping layout = model.RandomEffects!.Layout(block);
        double first = 0, second = 0;
        for (int i = 0; i < layout.Cases; i++)
        {
            double z = layout.Multiplier(i), effect = effects[layout.GroupOf(i)];
            first += z*effect; second += z*z*effect;
        }
        Assert.InRange(Math.Abs(first),0,1e-12);
        Assert.InRange(Math.Abs(second),0,1e-12);
    }

    [Fact]
    public void Configuration_OwnsOrderedSlopeValues_AndRequiresExplicitScaling()
    {
        BassRandomSlope original = new("g","x",RandomSlopeScale.Standardized,0.6,new(1.2));
        var source = new List<BassRandomSlope> { original };
        var config = new BassConfig { RandomSlopes = source };
        string json = config.ToJson(); int hash = config.GetHashCode();
        source[0] = original with { Scale = 3 }; source.Add(original with { Input = "y" });
        Assert.Equal(json,config.ToJson()); Assert.Equal(hash,config.GetHashCode());
        BassConfig roundtrip = BassConfig.FromJson(json);
        Assert.Equal(config,roundtrip); Assert.Equal(hash,roundtrip.GetHashCode());
        Assert.Equal(config,new BassConfig { RandomSlopes = new[] { original with { } } });
        Assert.False(config.RandomSlopes is IList<BassRandomSlope>);
        Assert.Throws<JsonException>(() => BassConfig.FromJson(json.Replace("\"InputScale\":\"Standardized\",","")));
        Assert.False(typeof(BassRandomSlope).GetConstructors().Single().GetParameters()[2].IsOptional);
    }

    [Fact]
    public void Declarations_RefuseMissingOrNoncontinuousInputs_InvalidScales_AndDuplicatePairs()
    {
        var (config,data) = Fixture(RandomSlopeScale.Unit);
        BassRandomSlope slope = config.RandomSlopes![0];
        Assert.Throws<ArgumentException>(() => (config with { RandomSlopes = [slope,slope] }).Validate());
        Assert.Throws<ArgumentException>(() => (config with { RandomSlopes = [slope with { InputScale = (RandomSlopeScale)99 }] }).Validate());
        Assert.Throws<ArgumentException>(() => (config with { RandomSlopes = [slope with { Scale = 0 }] }).Validate());
        Assert.Throws<ArgumentException>(() => BassModel.Create(config with { RandomSlopes = [slope with { Input = "missing" }] },data));
        var categorical = new BassData([new CategoricalPredictor("x",[0,1],["a","b"])],new ScalarResponse([0,1]),
            [new CategoricalPredictor("g",[0,1],["a","b"])]);
        Assert.Throws<ArgumentException>(() => BassModel.Create(config,categorical));
    }

    [Theory]
    [InlineData(RandomSlopeScale.Unit)]
    [InlineData(RandomSlopeScale.Standardized)]
    public void DependentConstraints_AndAZeroDimensionalSlope_KeepTheirExactRanks(RandomSlopeScale scaling)
    {
        // Identical within-group input moments leave only the population-slope constraint.
        var data = new BassData([new ContinuousPredictor("x",[0,1,0,1,0,1],new InputDomain(0,1))],new ScalarResponse([0,1,2,3,4,5]),
            [new CategoricalPredictor("g",[0,0,1,1,2,2],["a","b","c"])]);
        var config = new BassConfig { RandomSlopes = [new("g","x",scaling,0.7)] };
        BassModel model = BassModel.Create(config,data);
        Assert.Equal(2,model.RandomEffects!.Layout(0).Dimension);
        double[][] draws = model.RandomEffects.DrawEffects(model,[],[0],1,new Xoshiro256PlusPlus(210706));
        CheckConstraints(model,draws[0],0);
        // Two independent constraints on two observed levels fix both effects at zero.
        var small = new BassData([new ContinuousPredictor("x",[0,0.2,0.7,1],new InputDomain(0,1))],new ScalarResponse([0,1,2,3]),
            [new CategoricalPredictor("g",[0,0,1,1],["a","b"])]);
        BassModel zero = BassModel.Create(config,small);
        Assert.Equal(0,zero.RandomEffects!.Rank); Assert.Equal(0,zero.RandomEffects.Layout(0).Dimension);
        Assert.Equal(new double[2],zero.RandomEffects.DrawEffects(zero,[],[0],1,new Xoshiro256PlusPlus(210707))[0]);
    }

    [Fact]
    public void SeveralSlopeInputs_ResolveByName_ThroughMixedPredictorPositions()
    {
        var (config,original) = Fixture(RandomSlopeScale.Standardized);
        int n = original.Response.Cases;
        double[] w = [.. Enumerable.Range(0,n).Select(i => ((7*i)%29+0.5)/30)];
        var data = new BassData([new ContinuousPredictor("w",w,new InputDomain(0,1)),
            new CategoricalPredictor("category",[.. Enumerable.Range(0,n).Select(i => i%2)],["a","b"]),.. original.Predictors],original.Response,original.Groupings);
        config = config with { Structure = config.Structure with { Legendre = null },
            RandomSlopes = [config.RandomSlopes![0],new("g","w",RandomSlopeScale.Unit,1.2)] };
        BassModel model = BassModel.Create(config,data);
        Assert.Equal(3,model.RandomEffects!.Count);
        for (int k = 1; k < 3; k++)
        {
            double[] z = Input(data,config.RandomSlopes![k-1]);
            for (int i = 0; i < n; i++)
            {
                Assert.Equal(z[i],model.RandomEffects.Layout(k).Multiplier(i));
                Assert.Equal(z[i],model.RandomEffects.Layout(k).Multiplier(model.Inputs.Cases,i));
            }
        }
        var factor = new DenseCholesky(n); factor.Decompose(Covariance(Blocks(model)));
        Assert.InRange(Math.Abs(factor.LogDet-model.RandomEffects.LogDeterminant),0,1e-12);
        BassGroupDraw[] groups = model.RandomEffects.DrawGroups(model,[],[0],1,new Xoshiro256PlusPlus(210708));
        Assert.Equal(new string?[] {null,"x","w"},groups.Select(g => g.Input));
        Assert.Equal(groups,model.RandomEffects.Ordered(groups.Reverse().ToArray()));
        Assert.Throws<ArgumentException>(() => model.RandomEffects.Ordered([groups[0],groups[1],groups[1]]));
    }

    [Fact]
    public void A_RoundingSizedFirstMoment_DoesNotAddAnArtificialConstraint()
    {
        var data = new BassData([new ContinuousPredictor("x",[0.1,0.2,0.3,0.3,0.2,0.1,0.2,0.1,0.3],new InputDomain(0,1))],
            new ScalarResponse([0,1,2,3,4,5,6,7,8]),[new CategoricalPredictor("g",[0,0,0,1,1,1,2,2,2],["a","b","c"])]);
        var config = new BassConfig { RandomSlopes = [new("g","x",RandomSlopeScale.Standardized,0.7)] };
        BassModel model = BassModel.Create(config,data);
        Assert.Equal(2,model.RandomEffects!.Layout(0).Dimension);
        CheckConstraints(model,model.RandomEffects.DrawEffects(model,[],[0],1,new Xoshiro256PlusPlus(210709))[0],0);
    }
}
