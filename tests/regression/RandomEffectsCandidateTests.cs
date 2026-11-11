using System;
using System.Linq;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class RandomEffectsCandidateTests(ITestOutputHelper output)
{
    private static (BassFit Fit, BassData Data) Fixture(string law = "gaussian", bool slopes = false)
    {
        var (config, data) = slopes ? RandomSlopesReference.Fixture(RandomSlopeScale.Standardized)
            : RandomEffectsMultipleReference.Fixture("crossed", drawn: false);
        config = config with { Noise = BassCandidateTests.Noise(law), Intercept = GPriorCentre.Mean,
            Sampling = config.Sampling with { BurnIn = 40, Samples = 16 } };
        return (Bass.Fit(config, data), data);
    }

    [Theory]
    [InlineData("gaussian", false)]
    [InlineData("t", false)]
    [InlineData("quantile", false)]
    [InlineData("normal-wald", false)]
    [InlineData("skewed t", false)]
    [InlineData("generalized hyperbolic", false)]
    [InlineData("gaussian", true)]
    [InlineData("t", true)]
    public void KnownCasesMatchTrainingDensityAndCdfBitwise(string law, bool slopes)
    {
        var (fit, data) = Fixture(law, slopes);
        var cases = PredictionCases.Of(data);
        var cells = fit.GroupCells(data, cases);
        int[] identity = Enumerable.Range(0, cells.Training.Length).ToArray();
        var training = fit.PointwiseLogDensity(data, new UnitPartition(identity, cells.Training), quadrature: true, cdf: true);
        var candidate = fit.CandidateLogDensity(data, cases, data.Response.Values.ToArray().Select(y => new[] { y }).ToArray(),
            strata: cells.Candidates, quadrature: true, cdf: true);
        Assert.Equal(training.Latents, candidate.Latents);
        Assert.Equal(training.ChainLengths, candidate.ChainLengths);
        Assert.Null(candidate.NewGroup);
        Assert.Equal(ConformalForms.Both, candidate.Forms);
        for (int s = 0; s < training.Draws; s++)
            for (int t = 0; t < cases.Cases; t++)
            {
                Assert.Equal(BitConverter.DoubleToInt64Bits(training.Values[s,t]), BitConverter.DoubleToInt64Bits(candidate.Values[s,t]));
                Assert.Equal(BitConverter.DoubleToInt64Bits(training.LogCdf![s,t]), BitConverter.DoubleToInt64Bits(candidate.LogCdf![s,t]));
            }
        Assert.Equal(cases.Cases, ConformalBayes.FullStep(training, candidate).Sets.Count);
    }

    [Fact]
    public void CellsFrameArbitraryLabelsAndMissingColumnsSeparately()
    {
        var (fit, data) = Fixture();
        var cases = new PredictionCases([new ContinuousColumn("x", [0.4,0.4,0.4]),
            new CategoricalColumn("school", ["a|new:b", "a", "a|new:b"]),
            new CategoricalColumn("reader", ["c", "b|new:c", "c"])]);
        var cells = fit.GroupCells(data, cases, UnseenLevel.OutsideEverySubset);
        Assert.NotEqual(cells.Candidates![0], cells.Candidates[1]);
        Assert.Equal(cells.Candidates[0], cells.Candidates[2]);
        Assert.All(cells.NewGroup!, Assert.True);
        var missing = new PredictionCases([new ContinuousColumn("x", [0.4,0.4])]);
        var separate = fit.GroupCells(data, missing);
        Assert.NotEqual(separate.Candidates![0], separate.Candidates[1]);
    }

    [Fact]
    public void NestedCellsFrameParentAndChildAndParentPresence()
    {
        var (config,data) = RandomEffectsMultipleReference.Fixture("nested", drawn:false);
        config = config with { RandomIntercepts = [config.RandomIntercepts!.Single(s => s.Grouping == "child")],
            Sampling = config.Sampling with { BurnIn=20, Samples=8 } };
        var fit = Bass.Fit(config,data);
        var named = new PredictionCases([new ContinuousColumn("x",[0.4,0.4,0.4]),
            new CategoricalColumn("parent",["a/b","a","a/b"]), new CategoricalColumn("child",["c","b/c","c"])]);
        var cells = fit.GroupCells(data,named,UnseenLevel.OutsideEverySubset);
        Assert.NotEqual(cells.Candidates![0],cells.Candidates[1]); Assert.Equal(cells.Candidates[0],cells.Candidates[2]);
        // Compare within one call: a literal slash in a child must not impersonate the nested pair.
        var labels = new PredictionCases([new ContinuousColumn("x",[0.4]),new CategoricalColumn("child",["a/b/c"])]);
        var model = BassModel.Create(config,data);
        string missing = labels.GroupLabelsOf(model.RandomEffects!.Layout(0),data.Groupings.Single(g=>g.Name=="child"))![0];
        string present = named.GroupLabelsOf(model.RandomEffects.Layout(0),data.Groupings.Single(g=>g.Name=="child"))![0];
        Assert.NotEqual(missing,present);
    }

    [Theory]
    [InlineData("t")]
    [InlineData("quantile")]
    [InlineData("normal-wald")]
    [InlineData("skewed t")]
    [InlineData("generalized hyperbolic")]
    public void NewNoiseGroupsPublicSeamMatchesIndependentQuadrature(string law)
    {
        var (fit,data)=Fixture(law);
        var cases=new PredictionCases([new ContinuousColumn("x",[0.4])]);
        var density=fit.CandidateLogDensity(data,cases,[[-3.0,0,4.0]],fit.GroupCells(data,cases).Candidates,quadrature:true,cdf:true);
        var model=BassModel.Create(fit.Config,data);
        var (inputs,_)=cases.Prepare(model.Transform,model.Inputs,UnseenLevel.Refuse,data.Groupings);
        var surface=new BassSurface(model.Basis,fit.Terms,inputs,null);
        double worst=0, refinement=0; int s=0;
        foreach(var draw in fit.Draws.SelectMany(d=>d))
        {
            var mixture=(VarianceMeanMixture)model.Observation; var globals=draw.LatentsOn(model);
            var mixing=mixture.MixingAt(globals); var f=new double[1]; surface.Evaluate(draw.Structure.TermIds,draw.Coefficients,f);
            double extra=draw.Sigma2*draw.Groups!.Sum(g=>g.Scale*g.Scale),unit=draw.Sigma2*mixture.VarianceScale,drift=Math.Sqrt(draw.Sigma2)*(globals?.Skew??0);
            for(int g=0;g<3;g++)
            {
                double e=density.Grids[0][g]-f[0];
                double[] a=ReferenceGrid(mixing,e,unit,drift,extra,0.02),b=ReferenceGrid(mixing,e,unit,drift,extra,0.01);
                refinement=Math.Max(refinement,Math.Max(Math.Abs(a[0]-b[0]),Math.Abs(a[1]-b[1])));
                worst=Math.Max(worst,Math.Max(Math.Abs(density.Values[s,g]-b[0]),Math.Abs(density.LogCdf![s,g]-b[1])));
            }
            s++;
        }
        output.WriteLine($"{law}: refined independent grid {refinement:G17}, public density/CDF gap {worst:G17}");
        Assert.InRange(refinement,0,1e-8); Assert.InRange(worst,0,1e-8);
        Assert.Equal(LatentTreatment.Integrated,density.Latents.Single(l=>l.Scope==LatentScope.Group).Treatment);
    }

    private static double[] ReferenceGrid(IMixingLaw mixing,double e,double unit,double drift,double extra,double step)
    {
        double density=0,cdf=0;
        for(double u=-45;u<=45;u+=step)
        {
            double v=Math.Exp(u),variance=unit*v+extra,d=e-drift*v,weight=Math.Exp(mixing.LogPrior(v)+u)*step;
            density+=Math.Exp(-0.5*Math.Log(2*Math.PI*variance)-d*d/(2*variance))*weight;
            cdf+=Math.Exp(Normal.StandardLogCdf(d/Math.Sqrt(variance)))*weight;
        }
        return [Math.Log(density),Math.Log(cdf)];
    }

    [Theory]
    [InlineData("gaussian")]
    [InlineData("t")]
    public void KnownCellSplitCoverageBySimulation(string law)
    {
        if(!BassQualification.Enabled) return;
        const int replicates=24,groups=3,calibration=99,tests=100;
        int[] seeds=SeedTree.Derive(211105,replicates); var rates=new double[replicates,groups];
        for(int r=0;r<replicates;r++)
        {
            var rng=new Xoshiro256PlusPlus(seeds[r]);
            double Error()=>0.4*Normal.SampleStandard(rng)+(rng.NextDouble()<0.15 ? -0.9*Math.Log(1-rng.NextDouble()):0);
            double[] effects=[-1.4,0.1,1.3];
            int[] trainingGroups=Enumerable.Range(0,60).Select(i=>i%groups).ToArray();
            double[] x=Enumerable.Range(0,60).Select(_=>rng.NextDouble()).ToArray();
            var data=new BassData([new ContinuousPredictor("x",x,new InputDomain(0,1))],
                new ScalarResponse(x.Select((v,i)=>0.7+0.9*v+effects[trainingGroups[i]]+Error()).ToArray()),
                [new CategoricalPredictor("g",trainingGroups,["a","b","c"])]);
            var (baseConfig,_)=RandomEffectsMultipleReference.Fixture("crossed",drawn:false);
            var config=baseConfig with { Intercept=GPriorCentre.Mean,Noise=BassCandidateTests.Noise(law),
                RandomIntercepts=[new BassRandomIntercept("g",1.0)],Sampling=baseConfig.Sampling with { MasterSeed=seeds[r],BurnIn=150,Samples=120 } };
            var fit=Bass.Fit(config,data);
            CandidateLogDensity Read(int count)
            {
                int[] codes=Enumerable.Range(0,count*groups).Select(i=>i%groups).ToArray();
                double[] values=codes.Select(_=>rng.NextDouble()).ToArray();
                var cases=new PredictionCases([new ContinuousColumn("x",values),new CategoricalColumn("g",codes.Select(g=>new[]{"a","b","c"}[g]).ToArray())]);
                var cells=fit.GroupCells(data,cases);
                return fit.CandidateLogDensity(data,cases,values.Select((v,i)=>new[]{0.7+0.9*v+effects[codes[i]]+Error()}).ToArray(),cells.Candidates);
            }
            var sets=SplitConformal.Density(Read(calibration),Read(tests)).Sets;
            for(int g=0;g<groups;g++) rates[r,g]=sets.Where((_,i)=>i%groups==g).Count(s=>s.Accepts(0))/(double)tests;
        }
        for(int g=0;g<=groups;g++)
        {
            double[] values=Enumerable.Range(0,replicates).Select(r=>g==groups ? Enumerable.Range(0,groups).Average(k=>rates[r,k]):rates[r,g]).ToArray();
            double mean=values.Average(),se=Math.Sqrt(values.Sum(v=>(v-mean)*(v-mean))/(replicates-1)/replicates);
            output.WriteLine($"{law}, {(g==groups ? "marginal":$"cell {g}")}: coverage {mean:F5}, replicate SE {se:F5}; {replicates} fits, {tests} tests/cell/fit, 99 calibration/cell");
            Assert.True(Math.Abs(mean-0.9)<=4*Math.Max(se,Math.Sqrt(0.09/(replicates*tests))),$"{law} cell {g}: {mean}");
        }
    }

    [Fact]
    public void MixedGroupingReportsMixedAndNewFlagsArePerCase()
    {
        var (fit, data) = Fixture();
        var school = data.Groupings.Single(g => g.Name == "school");
        var reader = data.Groupings.Single(g => g.Name == "reader");
        var cases = new PredictionCases([new ContinuousColumn("x", [0.4,0.4]),
            new CategoricalColumn("school", [school.Levels[school.Codes.Span[0]], "new"]),
            new CategoricalColumn("reader", [reader.Levels[reader.Codes.Span[0]], reader.Levels[reader.Codes.Span[0]]])]);
        var cells = fit.GroupCells(data, cases, UnseenLevel.OutsideEverySubset);
        var density = fit.CandidateLogDensity(data, cases, [[0.0],[0.0]], cells.Candidates, UnseenLevel.OutsideEverySubset);
        Assert.Equal(new[] { false,true }, density.NewGroup);
        Assert.Equal(LatentTreatment.Mixed, density.Latents.Single(l => l.Scope == LatentScope.Group).Treatment);
        Assert.Equal(ConformalForms.Split, density.Forms);
        var oneKnownGrouping = new PredictionCases([new ContinuousColumn("x", [0.4]),
            new CategoricalColumn("school", [school.Levels[school.Codes.Span[0]]])]);
        var mixed = fit.CandidateLogDensity(data, oneKnownGrouping, [[0.0]], fit.GroupCells(data, oneKnownGrouping).Candidates);
        Assert.Equal(LatentTreatment.Mixed, mixed.Latents.Single(l => l.Scope == LatentScope.Group).Treatment);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void NewGaussianGroupsMatchIndependentNormalReference(bool slopes)
    {
        var (fit, data) = Fixture(slopes: slopes);
        var cases = new PredictionCases([new ContinuousColumn("x", [0.2,0.8])]);
        var cells = fit.GroupCells(data, cases);
        double[][] grids = [[-3.0,0,4.0],[-3.0,0,4.0]];
        var candidate = fit.CandidateLogDensity(data, cases, grids, cells.Candidates, cdf:true);
        var model = BassModel.Create(fit.Config, data);
        var (inputs, _) = cases.Prepare(model.Transform, model.Inputs, UnseenLevel.Refuse, data.Groupings);
        var surface = new BassSurface(model.Basis, fit.Terms, inputs, null);
        int s = 0;
        foreach (var draw in fit.Draws.SelectMany(d => d))
        {
            double[] f = new double[2]; surface.Evaluate(draw.Structure.TermIds, draw.Coefficients, f);
            for (int t = 0; t < 2; t++)
            {
                double extra = 0;
                foreach (var group in draw.Groups!)
                {
                    double[] x = ((ContinuousPredictor)data.Predictors[0]).Values.ToArray();
                    double mean = x.Average(), sd = Math.Sqrt(x.Select(v => (v-mean)*(v-mean)).Average());
                    double z = group.Input is null ? 1 : (new[] {0.2,0.8}[t]-mean)/sd;
                    extra += group.Scale*group.Scale*z*z;
                }
                var normal = Normal.MeanVariance(f[t], draw.Sigma2*(1+extra));
                for (int g = 0; g < grids[t].Length; g++)
                {
                    Assert.InRange(Math.Abs(candidate.Values[s,candidate.Offset(t)+g]-normal.LogPdf(grids[t][g])),0,1e-11);
                    Assert.InRange(Math.Abs(candidate.LogCdf![s,candidate.Offset(t)+g]-normal.LogCdf(grids[t][g])),0,1e-11);
                }
            }
            s++;
        }
        Assert.Equal(LatentTreatment.Integrated, candidate.Latents.Single(l => l.Scope == LatentScope.Group).Treatment);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void WidenedMixtureDensityAndCdfMatchRefinedIndependentGrid(bool exponential)
    {
        // Direct trapezoids with explicit IG(2,2) or Exp(1) priors, independent of production's adaptive panels and LogPrior.
        IMixingLaw law = exponential ? ExponentialMixing.Instance : new StudentTMixing(4);
        double worst = 0, refinement = 0;
        foreach (double drift in new[] {-0.7,0,1.2})
            foreach (double extra in new[] {0.0,0.01,2.0,30.0})
                foreach (double e in new[] {-12.0,-2.0,0,3.0,16.0})
                {
                    double[] coarse = Integrate(exponential,e,drift,extra,0.01);
                    double[] fine = Integrate(exponential,e,drift,extra,0.005);
                    for (int j=0;j<2;j++) refinement=Math.Max(refinement,Math.Abs(coarse[j]-fine[j]));
                    worst = Math.Max(worst, Math.Abs(BassWidenedMixture.LogDensity(law,e,0.7,drift,extra)-fine[0]));
                    worst = Math.Max(worst, Math.Abs(BassWidenedMixture.LogCdf(law,e,0.7,drift,extra)-fine[1]));
                }
        output.WriteLine($"exponential={exponential}: grid refinement {refinement:G17}, worst log density/CDF gap {worst:G17}");
        Assert.InRange(refinement,0,2e-8); Assert.InRange(worst,0,2e-8);
    }

    [Theory(Timeout = 10000)]
    [InlineData(0.0)]
    [InlineData(1.0)]
    public async Task WidenedMixtureBoundaryPeakTerminatesAndMatchesStudentTail(double extra) => await Task.Run(() =>
    {
        // For IG(2,2), the log-density integrand peaks at log((e²/2 + 2)/2.5).
        // Put that peak exactly at the initial scan endpoint; keeping only the adjacent interval used to loop forever.
        double e = Math.Sqrt(2.0 * (2.5 * Math.Exp(40.0) - 2.0));
        var reference = StudentT.Standard(4);
        var mixing = new StudentTMixing(4);
        foreach (double sign in new[] { -1.0, 1.0 })
        {
            double residual = sign * e;
            double density = BassWidenedMixture.LogDensity(mixing, residual, 1.0, 0.0, extra);
            double cdf = BassWidenedMixture.LogCdf(mixing, residual, 1.0, 0.0, extra);
            // Adding N(0,1) changes this far-tail density/CDF only by O(1/e²), below double precision here.
            Assert.InRange(Math.Abs(density - reference.LogPdf(residual)), 0.0, 1e-10);
            Assert.InRange(Math.Abs(cdf - reference.LogCdf(residual)), 0.0, 1e-10);
        }
    });

    private static double[] Integrate(bool exponential,double e,double drift,double extra,double step)
    {
        double density=0,cdf=0;
        for (double u=-40;u<=40;u+=step)
        {
            double v=Math.Exp(u), variance=0.7*v+extra, d=e-drift*v;
            double prior=exponential ? -v+u : Math.Log(4)-2*u-2/v;
            density += Math.Exp(-0.5*Math.Log(2*Math.PI*variance)-d*d/(2*variance)+prior)*step;
            cdf += Math.Exp(Normal.StandardLogCdf(d/Math.Sqrt(variance))+prior)*step;
        }
        return [Math.Log(density),Math.Log(cdf)];
    }

    [Fact]
    public void RefusalsNameCellAndQuadratureRequirements()
    {
        var (fit,data)=Fixture(); var cases=PredictionCases.Of(data);
        double[][] grids=data.Response.Values.ToArray().Select(y=>new[]{y}).ToArray();
        Assert.Contains("GroupCells",Assert.Throws<ArgumentException>(()=>fit.CandidateLogDensity(data,cases,grids)).Message);
        Assert.Contains("spans cells",Assert.Throws<ArgumentException>(()=>fit.CandidateLogDensity(data,cases,grids,new int[cases.Cases])).Message);
        var cells=fit.GroupCells(data,cases);
        int[] bad=(int[])cells.Training.Clone(); bad[0]=123; bad[1]=123;
        // A deliberately mixed stratum across all training cells, plus one separate unit.
        Array.Fill(bad,0); bad[0]=1;
        Assert.Contains("spans cells",Assert.Throws<ArgumentException>(()=>fit.PointwiseLogDensity(data,new UnitPartition(Enumerable.Range(0,bad.Length).ToArray(),bad))).Message);
        var (noise,nd)=Fixture("t"); var fresh=new PredictionCases([new ContinuousColumn("x",[0.5])]);
        Assert.Contains("quadrature: true",Assert.Throws<NotSupportedException>(()=>noise.CandidateLogDensity(nd,fresh,[[0.0]],noise.GroupCells(nd,fresh).Candidates)).Message);
    }
}
