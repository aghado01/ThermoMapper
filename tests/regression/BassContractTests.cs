using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's contracts (stage 4, S1): the data refuses its defects and copies what it is given; the input transform maps
/// the training range onto [0, 1] with exact ends and keeps the observed levels; the config round-trips through JSON
/// by value, refuses out-of-range fields by name, and the fingerprint moves with any setting or any data bit; the
/// resolved record matches direct computation and refuses what the intercept convention cannot fit.
/// </summary>
public sealed class BassContractTests
{
    private static ContinuousPredictor X(string name, params double[] v) => new(name, v);
    private static ScalarResponse Y(params double[] v) => new(v);

    private static BassData Scalar(int n, Func<int, double> y)
    {
        double[] x1 = Enumerable.Range(0, n).Select(i => (double)i).ToArray();
        double[] x2 = Enumerable.Range(0, n).Select(i => Math.Sin(i)).ToArray();
        int[] codes = Enumerable.Range(0, n).Select(i => i % 3).ToArray();
        return new BassData(
            [X("x1", x1), new CategoricalPredictor("z", codes, ["a", "b", "c", "unused"]), X("x2", x2)],
            Y(Enumerable.Range(0, n).Select(y).ToArray()));
    }

    private static BassData Functional(int cases, int locations)
    {
        var y = new double[cases, locations];
        for (int i = 0; i < cases; i++)
            for (int r = 0; r < locations; r++) y[i, r] = i + Math.Cos(r);
        double[] t = Enumerable.Range(0, locations).Select(r => r / (double)locations).ToArray();
        double[] x = Enumerable.Range(0, cases).Select(i => i * 0.5).ToArray();
        return new BassData([X("x", x)], new FunctionalResponse(y, [X("t", t)]));
    }

    [Fact]
    public void Data_RefusesItsDefects_WithReasons()
    {
        Assert.Throws<ArgumentException>(() => new BassData([], Y(1, 2)));
        Assert.Contains("twice", Assert.Throws<ArgumentException>(() => new BassData([X("a", 1, 2), X("a", 3, 4)], Y(1, 2))).Message);
        Assert.Contains("2 cases", Assert.Throws<ArgumentException>(() => new BassData([X("a", 1, 2, 3)], Y(1, 2))).Message);
        Assert.Contains("no knot", Assert.Throws<ArgumentException>(() => X("a", 2, 2, 2)).Message);
        Assert.Throws<ArgumentException>(() => X("a", 1, double.NaN));
        Assert.Throws<ArgumentException>(() => X("a"));
        Assert.Throws<ArgumentException>(() => X(" ", 1, 2));
        Assert.Contains("outside", Assert.Throws<ArgumentException>(() => new CategoricalPredictor("z", [0, 2], ["a", "b"])).Message);
        Assert.Contains("twice", Assert.Throws<ArgumentException>(() => new CategoricalPredictor("z", [0, 1], ["a", "a"])).Message);
        Assert.Contains("no information", Assert.Throws<ArgumentException>(() => new CategoricalPredictor("z", [1, 1], ["a", "b"])).Message);
        Assert.Throws<ArgumentException>(() => Y(1, double.PositiveInfinity));
        Assert.Throws<ArgumentException>(() => Y());

        var y = new double[2, 3];
        Assert.Throws<ArgumentException>(() => new FunctionalResponse(y, []));
        Assert.Contains("3 locations", Assert.Throws<ArgumentException>(() => new FunctionalResponse(y, [X("t", 0, 1)])).Message);
        Assert.Contains("twice", Assert.Throws<ArgumentException>(
            () => new BassData([X("t", 1, 2)], new FunctionalResponse(y, [X("t", 0, 1, 2)]))).Message);
        y[1, 2] = double.NaN;
        Assert.Throws<ArgumentException>(() => new FunctionalResponse(y, [X("t", 0, 1, 2)]));
    }

    [Fact]
    public void Data_CopiesWhatItIsGiven()
    {
        double[] x = [1, 2, 3], yv = [4, 5, 6];
        int[] codes = [0, 1, 0];
        var data = new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", codes, ["a", "b"])], new ScalarResponse(yv));
        string before = Bass.Fingerprint(new BassConfig(), data);
        x[0] = 99; yv[0] = 99; codes[0] = 1;
        Assert.Equal(1.0, ((ContinuousPredictor)data.Predictors[0]).Values.Span[0]);
        Assert.Equal(4.0, ((ScalarResponse)data.Response).Y.Span[0]);
        Assert.Equal(0, ((CategoricalPredictor)data.Predictors[1]).Codes.Span[0]);
        Assert.Equal(before, Bass.Fingerprint(new BassConfig(), data));

        var grid = new double[,] { { 1, 2 }, { 3, 4 } };
        var functional = new FunctionalResponse(grid, [X("t", 0, 1)]);
        grid[1, 0] = 99;
        Assert.Equal(3.0, functional[1, 0]);
        Assert.Equal(4, functional.Observations);
    }

    [Fact]
    public void Transform_MapsTheTrainingRange_OntoTheUnitInterval_WithExactEnds()
    {
        BassInputTransform transform = BassInputTransform.Fit(Scalar(10, i => i * i));
        var x1 = (RangeScale)transform.Predictors[0];
        Assert.Equal((0.0, 9.0), (x1.Min, x1.Max));
        Assert.Equal(0.0, x1.Scale(0.0));
        Assert.Equal(1.0, x1.Scale(9.0));
        Assert.Equal(2.0, x1.Scale(18.0));       // linear past the range, where the hinges extrapolate
        Assert.Equal(-1.0 / 9.0, x1.Scale(-1.0));
        Assert.True(x1.Contains(0.0) && x1.Contains(9.0));
        Assert.False(x1.Contains(9.000001) || x1.Contains(-1e-12));

        var x2 = (RangeScale)transform.Predictors[2];
        double[] sines = Enumerable.Range(0, 10).Select(i => Math.Sin(i)).ToArray();
        Assert.Equal(0.0, x2.Scale(sines.Min()));
        Assert.Equal(1.0, x2.Scale(sines.Max()));

        // The training levels are the declared levels some case takes, in the declared order.
        var z = (TrainingLevels)transform.Predictors[1];
        Assert.Equal(["a", "b", "c"], z.Levels);
        Assert.Equal(2, z.IndexOf("c"));
        Assert.Equal(-1, z.IndexOf("unused"));
        Assert.Equal(-1, z.IndexOf("never declared"));
        Assert.Empty(transform.OutputVariables);

        BassInputTransform f = BassInputTransform.Fit(Functional(4, 5));
        Assert.Equal((0.0, 0.8), (f.OutputVariables[0].Min, f.OutputVariables[0].Max));
        Assert.Equal(1.0, f.OutputVariables[0].Scale(0.8));
    }

    [Fact]
    public void Transform_TrainingLevels_SkipUnobservedCodes()
    {
        var factor = new CategoricalPredictor("z", [2, 0, 2], ["a", "gap", "c"]);
        BassInputTransform transform = BassInputTransform.Fit(new BassData([factor], Y(1, 2, 3)));
        var z = (TrainingLevels)transform.Predictors[0];
        Assert.Equal(["a", "c"], z.Levels);
        Assert.Equal(1, z.IndexOfCode(2));
        Assert.Equal(-1, z.IndexOfCode(1));
    }

    [Fact]
    public void Config_HasTheDecidedDefaults()
    {
        var c = new BassConfig();
        Assert.Equal((1, 3, 3, 3, (int?)null), (c.Structure.Degree, c.Structure.MaxInteraction, c.Structure.MaxCategoricalInteraction,
                                                c.Structure.MaxFunctionalInteraction, c.Structure.MaxBasis));
        Assert.Equal(new SupportRule(SupportKind.Count, null, null), c.Structure.Support);
        Assert.Equal((1e-10, (double?)null, CategoricalSubsetMeasure.SizeUniform, 5.0, 5.0),
                     (c.Structure.AngleFloor, c.Structure.MaxVarianceInflation, c.Structure.CategoricalSubsets,
                      c.Structure.OrderWeight, c.Structure.VariableWeight));
        Assert.Equal(new BassCountPrior(10, 10), c.Count);
        Assert.Equal(new SampledG(0.5), c.Coefficients);
        Assert.Equal(Sigma2Prior.ScaledToResponse(1e-3), c.Sigma2);
        Assert.Equal(GPriorCentre.Mean, c.Intercept);
        Assert.Equal((4, 0, 9000, 1000, 1, 0.0, 0, 100.0, true),
                     (c.Sampling.Chains, c.Sampling.MasterSeed, c.Sampling.BurnIn, c.Sampling.Samples, c.Sampling.Thin,
                      c.Sampling.RHatTarget, c.Sampling.BatchSize, c.Sampling.MinEssPerChain, c.Sampling.ReportEvidence));
        Assert.Equal(new EqualOverAvailable(), c.Sampling.Moves);
        c.Validate();
    }

    [Fact]
    public void Config_RoundTripsThroughJson_ByValue()
    {
        var configs = new[]
        {
            new BassConfig(),
            new BassConfig
            {
                Structure = new BassStructureSpec
                {
                    Degree = 2, MaxBasis = 50, Support = new SupportRule(SupportKind.ColumnMass, 4.5, 3),
                    MaxVarianceInflation = 1e6, CategoricalSubsets = CategoricalSubsetMeasure.Uniform, OrderWeight = 2.5,
                },
                Count = new BassCountPrior(3, 0.5),
                Coefficients = new FixedG(0.01),
                Sigma2 = Sigma2Prior.Jeffreys,
                Intercept = GPriorCentre.Flat,
                Sampling = new BassSampling { Chains = 2, MasterSeed = 7, RHatTarget = 1.01, BatchSize = 200, Moves = new PriorCancelling(0.3) },
            },
            new BassConfig { Coefficients = new FixedG(), Sigma2 = Sigma2Prior.InverseGamma(2, 0.5), Intercept = GPriorCentre.Zero },
        };
        foreach (BassConfig config in configs)
        {
            string json = config.ToJson();
            Assert.Equal(config, BassConfig.FromJson(json));
            Assert.Equal(json, BassConfig.FromJson(json).ToJson());
        }
        string defaults = new BassConfig().ToJson();
        Assert.Contains("\"Coefficients\":{\"Kind\":\"SampledG\",\"TauShape\":0.5}", defaults);
        Assert.Contains("\"Intercept\":\"Mean\"", defaults);
        Assert.Contains("\"Moves\":{\"Kind\":\"EqualOverAvailable\"}", defaults);
    }

    [Fact]
    public void Config_RefusesOutOfRangeFields_ByName()
    {
        string Reason(BassConfig c) => Assert.Throws<ArgumentException>(c.Validate).Message;
        Assert.StartsWith("Structure.Degree", Reason(new BassConfig { Structure = new BassStructureSpec { Degree = 0 } }));
        Assert.StartsWith("Structure.MaxBasis", Reason(new BassConfig { Structure = new BassStructureSpec { MaxBasis = 0 } }));
        Assert.StartsWith("Structure.OrderWeight", Reason(new BassConfig { Structure = new BassStructureSpec { OrderWeight = 0 } }));
        Assert.StartsWith("Structure.AngleFloor", Reason(new BassConfig { Structure = new BassStructureSpec { AngleFloor = 1 } }));
        Assert.StartsWith("Structure.Support.Minimum", Reason(new BassConfig { Structure = new BassStructureSpec { Support = new SupportRule(Minimum: -1) } }));
        Assert.StartsWith("Count.Rate", Reason(new BassConfig { Count = new BassCountPrior(10, 0) }));
        Assert.StartsWith("Coefficients.TauShape", Reason(new BassConfig { Coefficients = new SampledG(double.NaN) }));
        Assert.StartsWith("Coefficients.Tau", Reason(new BassConfig { Coefficients = new FixedG(-1) }));
        Assert.StartsWith("Intercept", Reason(new BassConfig { Intercept = GPriorCentre.Median }));
        Assert.StartsWith("Sampling.RHatTarget", Reason(new BassConfig { Sampling = new BassSampling { RHatTarget = 0.9 } }));
        Assert.StartsWith("Sampling.Chains", Reason(new BassConfig { Sampling = new BassSampling { Chains = 1, RHatTarget = 1.1 } }));
        Assert.StartsWith("Sampling.Moves.C", Reason(new BassConfig { Sampling = new BassSampling { Moves = new PriorCancelling(0.6) } }));
        Assert.StartsWith("Sampling", Reason(new BassConfig { Sampling = null! }));

        // A manifest is checked on load as a built config is.
        string bad = new BassConfig().ToJson().Replace("\"Degree\":1", "\"Degree\":0");
        Assert.StartsWith("Structure.Degree", Assert.Throws<ArgumentException>(() => BassConfig.FromJson(bad)).Message);
    }

    [Fact]
    public void Json_RefusesAMisspelledField_AndReadsAVariantInAnyOrder()
    {
        string json = new BassConfig().ToJson();
        string typo = json.Replace("\"Degree\":1", "\"Degre\":2");
        Assert.Contains("Degre", Assert.ThrowsAny<System.Text.Json.JsonException>(() => BassConfig.FromJson(typo)).Message);
        string sigmaTypo = json.Replace("\"Shape\":0.001", "\"Shap\":0.001");
        Assert.ThrowsAny<System.Text.Json.JsonException>(() => BassConfig.FromJson(sigmaTypo));

        string reordered = json.Replace("{\"Kind\":\"SampledG\",\"TauShape\":0.5}", "{\"TauShape\":0.25,\"Kind\":\"SampledG\"}");
        Assert.Equal(new SampledG(0.25), BassConfig.FromJson(reordered).Coefficients);

        // A field left out takes its default.
        BassConfig partial = BassConfig.FromJson("{\"Sampling\":{\"Chains\":2}}");
        Assert.Equal(new BassConfig { Sampling = new BassSampling { Chains = 2 } }, partial);
    }

    [Fact]
    public void Fingerprint_MovesWithAnySetting_AndAnyDataBit()
    {
        BassData data = Scalar(30, i => Math.Exp(-i / 10.0));
        var config = new BassConfig();
        string baseline = Bass.Fingerprint(config, data);
        Assert.Equal(baseline, Bass.Fingerprint(BassConfig.FromJson(config.ToJson()), data));
        Assert.Equal(baseline, Bass.Fingerprint(new BassConfig(), Scalar(30, i => Math.Exp(-i / 10.0))));

        Assert.NotEqual(baseline, Bass.Fingerprint(config with { Sampling = config.Sampling with { MasterSeed = 1 } }, data));
        Assert.NotEqual(baseline, Bass.Fingerprint(config with { Sigma2 = Sigma2Prior.ScaledToResponse(2e-3) }, data));
        Assert.NotEqual(baseline, Bass.Fingerprint(config with { Structure = config.Structure with { AngleFloor = 2e-10 } }, data));
        Assert.NotEqual(baseline, Bass.Fingerprint(config with { Coefficients = new SampledG(0.25) }, data));

        // One ulp in one response value.
        Assert.NotEqual(baseline, Bass.Fingerprint(config, Scalar(30, i => i == 7 ? Math.BitIncrement(Math.Exp(-0.7)) : Math.Exp(-i / 10.0))));
        Assert.NotEqual(Bass.Fingerprint(config, Functional(3, 4)), Bass.Fingerprint(config, Functional(4, 3)));
    }

    [Fact]
    public void Resolved_MatchesDirectComputation()
    {
        BassData data = Scalar(40, i => Math.Sin(i / 5.0) + 3.0);
        BassResolved r = BassResolved.Resolve(new BassConfig(), data);
        double[] y = Enumerable.Range(0, 40).Select(i => Math.Sin(i / 5.0) + 3.0).ToArray();
        double mean = y.Average(), meanSquare = y.Select(v => (v - mean) * (v - mean)).Average();
        Assert.Equal((40, 1, 40, 39), (r.Cases, r.Locations, r.Observations, r.MaxBasis));
        Assert.Equal((2, 1, 0), (r.MaxInteraction, r.MaxCategoricalInteraction, r.MaxFunctionalInteraction));
        Assert.Equal(4.0, r.SupportMinimum);
        Assert.Null(r.SupportMinimumFunctional);
        Assert.Equal(mean, r.CentreLevel, 1e-14);
        Assert.Equal(meanSquare, r.ResponseMeanSquare, 1e-14);
        Assert.Equal(1e-3, r.Sigma2!.Value.Shape);
        Assert.Equal(1e-3 * meanSquare, r.Sigma2!.Value.Scale, 1e-16);
        Assert.Equal((0.5, 20.0), (r.TauPrior!.Value.Shape, r.TauPrior!.Value.Rate));
        Assert.Null(r.FixedTau);
        Assert.Contains(r.AutoPicks, p => p.Setting == "Levels of 'z'" && p.Value == "'unused' left out");
        Assert.Contains(r.AutoPicks, p => p.Setting == "MaxInteraction" && p.Value == "2");
        Assert.Contains(r.AutoPicks, p => p.Setting == "SupportMinimum" && p.Value == "4");

        // A large response caps the basis count at R's 1,000; a declared cap above N − 1 is cut to it, and logged.
        Assert.Equal(1000, BassResolved.Resolve(new BassConfig(), Scalar(5000, i => i % 7)).MaxBasis);
        BassResolved capped = BassResolved.Resolve(new BassConfig { Structure = new BassStructureSpec { MaxBasis = 500 } }, data);
        Assert.Equal(39, capped.MaxBasis);
        Assert.Contains(capped.AutoPicks, p => p.Setting == "MaxBasis" && p.Reason.Contains("exceeds"));
    }

    [Fact]
    public void Resolved_Functional_CountsBothParts()
    {
        BassResolved r = BassResolved.Resolve(new BassConfig(), Functional(30, 250));
        Assert.Equal((30, 250, 7500, 1000), (r.Cases, r.Locations, r.Observations, r.MaxBasis));
        Assert.Equal((1, 0, 1), (r.MaxInteraction, r.MaxCategoricalInteraction, r.MaxFunctionalInteraction));
        Assert.Equal(3.0, r.SupportMinimum);
        Assert.Equal(20.0, r.SupportMinimumFunctional);
        Assert.Equal(0.5 * 7500, r.TauPrior!.Value.Rate);
    }

    [Fact]
    public void Resolved_FollowsTheIntercept_AndTheCoefficientMode()
    {
        BassData data = Scalar(20, i => 0.1 * i);
        BassResolved zero = BassResolved.Resolve(new BassConfig { Intercept = GPriorCentre.Zero, Sigma2 = Sigma2Prior.Jeffreys }, data);
        Assert.Equal(0.0, zero.CentreLevel);
        Assert.Equal(Enumerable.Range(0, 20).Select(i => 0.01 * i * i).Average(), zero.ResponseMeanSquare, 1e-14);
        Assert.Null(zero.Sigma2);

        BassResolved fixedG = BassResolved.Resolve(new BassConfig { Coefficients = new FixedG() }, data);
        Assert.Equal(1.0 / 20, fixedG.FixedTau);
        Assert.Null(fixedG.TauPrior);
        Assert.Equal(0.3, BassResolved.Resolve(new BassConfig { Coefficients = new FixedG(0.3) }, data).FixedTau);
    }

    [Fact]
    public void Resolved_RefusesWhatTheConventionCannotFit()
    {
        BassData constant = Scalar(20, _ => 2.5), zero = Scalar(20, _ => 0.0);
        foreach (GPriorCentre centre in new[] { GPriorCentre.Mean, GPriorCentre.Flat })
            Assert.Contains("constant", Assert.Throws<ArgumentException>(() => BassResolved.Resolve(new BassConfig { Intercept = centre }, constant)).Message);
        var zeroCentre = new BassConfig { Intercept = GPriorCentre.Zero };
        Assert.Equal(2.5, Math.Sqrt(BassResolved.Resolve(zeroCentre, constant).ResponseMeanSquare), 14);
        Assert.Contains("zero", Assert.Throws<ArgumentException>(() => BassResolved.Resolve(zeroCentre, zero)).Message);

        var unreachable = new BassConfig { Structure = new BassStructureSpec { Support = new SupportRule(Minimum: 21) } };
        Assert.Contains("no term could meet it", Assert.Throws<ArgumentException>(() => BassResolved.Resolve(unreachable, Scalar(20, i => i))).Message);
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(new BassConfig { Intercept = GPriorCentre.Median }, Scalar(20, i => i)));
    }
}
