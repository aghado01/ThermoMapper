using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputReadoutTests(ITestOutputHelper output)
{
    internal static (BassFit Fit, BassData Data) Fixed(string form = "modified", GPriorCentre centre = GPriorCentre.Zero, int count = 40, bool hierarchy = false)
    {
        var data = MultiOutputCoefficientTests.Data(2);
        var config = MultiOutputCoefficientTests.Config(form == "g" ? "raw" : form, centre) with {
            Sampling = new BassSampling { Chains = 1, BurnIn = 20, Samples = 20, MasterSeed = 20261216 }
        };
        if (form == "g") config = config with { Coefficients = BassCoefficientPrior.G() };
        if (hierarchy) config = config with { OutputCovariance = OutputCovariancePrior.HuangWand(2, [1, 0.7]) };
        var fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        var model = BassModel.Create(config, data);
        var terms = MultiOutputCoefficientTests.Terms(form);
        var table = new BassTermTable(); var store = new DrawStore<BassDraw>(); var rng = new Xoshiro256PlusPlus(20261217);
        var moves = new BassStructureMoves(model, terms, Math.Log(0.3));
        for (int i = 0; i < count; i++)
        {
            double tau = i == 0 ? 1e-8 : i == count - 1 ? 1e4 : Math.Exp(-4 + 5.0 * i / count);
            var state = moves.Start.WithLogTau(Math.Log(tau));
            if (hierarchy) state = state.WithOutputLogRate(0, -0.4 + 0.8 * i / count).WithOutputLogRate(1, -1.1 + 0.5 * i / count);
            store.Add(moves.Compose(state, rng, table));
        }
        double[] weights = [.. Enumerable.Range(0, count).Select(i => (double)(i % 3 + 1))];
        double sum = weights.Sum(); for (int i = 0; i < weights.Length; i++) weights[i] /= sum;
        return (fit with { Terms = table, Draws = [store], DrawWeights = weights }, data);
    }

    internal static PredictionCases Cases() => new([new ContinuousColumn("x", [0.25, 0.8]), new ContinuousColumn("z", [0.4, 0.65])]);

    [Theory]
    [InlineData("g", GPriorCentre.Zero)]
    [InlineData("g", GPriorCentre.Mean)]
    [InlineData("g", GPriorCentre.Flat)]
    [InlineData("ridge", GPriorCentre.Zero)]
    [InlineData("ridge", GPriorCentre.Mean)]
    [InlineData("ridge", GPriorCentre.Flat)]
    [InlineData("raw", GPriorCentre.Zero)]
    [InlineData("raw", GPriorCentre.Mean)]
    [InlineData("raw", GPriorCentre.Flat)]
    [InlineData("modified", GPriorCentre.Zero)]
    [InlineData("modified", GPriorCentre.Mean)]
    [InlineData("modified", GPriorCentre.Flat)]
    [InlineData("ridge", GPriorCentre.Mean, true)]
    public void ConditionalMatrixTPredictiveMatchesDenseConditioning(string form, GPriorCentre centre, bool hierarchy = false)
    {
        var (fit, data) = Fixed(form, centre, hierarchy: hierarchy);
        var model = BassModel.Create(fit.Config, data); var cases = Cases();
        var result = BassMultiOutput.Conditionals(fit, data, cases);
        var (inputs, _) = cases.Prepare(model.Transform, model.Inputs, UnseenLevel.Refuse);
        var surface = new BassSurface(model.Basis, fit.Terms, inputs, null);
        var y = new double[data.Response.Cases, 2]; var response = (MultiOutputResponse)data.Response;
        for (int i = 0; i < y.GetLength(0); i++) for (int k = 0; k < 2; k++) y[i, k] = response[i, k];
        double worst = 0;
        for (int d = 0; d < fit.Draws[0].Count; d++)
        {
            var draw = fit.Draws[0][d]; var design = model.NewDesign();
            Assert.Null(design.Load(draw.Structure.TermIds.Select(id => fit.Terms[id]).ToArray()));
            var c = form == "g" ? MultiOutputCase.GForm(y, MultiOutputCoreTests.Columns(model, design), centre, draw.Tau, 1)
                                : MultiOutputCoefficientTests.DenseCase(model, design, draw.Tau, 1);
            var prior = hierarchy ? model.Resolved.Outputs!.Rates!.PriorAt(draw.OutputRates!.Select(x => Math.Log(x)).ToArray()) : model.Resolved.Outputs!.Prior;
            var reference = MultiOutputPredictiveReference.Fit(c, new ReferenceSigmaPrior(prior.DegreesOfFreedom, prior.Scale));
            var z = new double[2, c.Columns];
            int offset = centre == GPriorCentre.Flat ? 0 : 1;
            if (offset == 1) for (int i = 0; i < 2; i++) z[i, 0] = 1;
            for (int j = 0; j < design.Count; j++)
            {
                var column = new double[data.Response.Cases]; model.Basis.CaseColumn(design.Terms[j], model.Inputs.Cases, column);
                double mean = column.Average();
                var fresh = surface.Columns(draw.Structure.TermIds[j]).Case;
                for (int i = 0; i < 2; i++) z[i, j + offset] = fresh[i] - mean;
            }
            var (expectedMean, expectedRows) = reference.At(z);
            var actual = result.Components[d];
            foreach (var pair in actual.Location.Cast<double>().Zip(expectedMean.Cast<double>())) worst = Math.Max(worst, Math.Abs(pair.First - pair.Second));
            foreach (var pair in actual.RowScale.Cast<double>().Zip(expectedRows.Cast<double>())) worst = Math.Max(worst, Math.Abs(pair.First - pair.Second));
            var value = new double[,] { { 0.3, -0.1 }, { 0.5, 1.2 } };
            worst = Math.Max(worst, Math.Abs(actual.LogDensity(value) - reference.LogDensity(z, value)));
            Assert.Equal(reference.Nu, actual.Nu);
        }
        output.WriteLine($"{form}, {centre}, HW={hierarchy}: dense predictive maximum absolute error {worst:G5}");
        Assert.InRange(worst, 0, 2e-9);
    }

    [Fact]
    public void JointDrawsBandsProbabilitiesAndCovarianceRetainTheirWeights()
    {
        var (fit, data) = Fixed(count: 4000);
        var result = BassMultiOutput.Predict(fit, data, Cases(), new BassPredictOptions { Seed = 20261218, MeanSurface = true });
        var again = BassMultiOutput.Predict(fit, data, Cases(), new BassPredictOptions { Seed = 20261218 });
        Assert.Equal(result.Predictive!.Cast<double>(), again.Predictive!.Cast<double>());
        for (int d = 0; d < result.Draws.Count; d++) Assert.Equal(fit.DrawWeights![d], result.Weights[d], 14);
        var model = BassModel.Create(fit.Config, data); var (inputs, _) = Cases().Prepare(model.Transform, model.Inputs, UnseenLevel.Refuse);
        var surface = new BassSurface(model.Basis, fit.Terms, inputs, null);
        for (int d = 0; d < 20; d++) for (int k = 0; k < 2; k++)
        {
            var draw = fit.Draws[0][d]; var coefficients = Enumerable.Range(0, draw.OutputCoefficients!.GetLength(0)).Select(i => draw.OutputCoefficients[i, k]).ToArray();
            var expected = new double[2]; surface.Evaluate(draw.Structure.TermIds, coefficients, expected);
            for (int i = 0; i < 2; i++) Assert.Equal(expected[i], result.Latent[d, i, k]);
        }
        var bands = result.JointBands();
        var mass = result.SpecificationProbability(bands.Lower, bands.Upper);
        foreach (double p in mass) Assert.InRange(p, 0.947, 0.953);
        var covariance = BassMultiOutput.Covariance(fit, data);
        Assert.NotNull(covariance.Mean);
        for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++)
        {
            double expected = fit.Draws[0].Select((d, i) => d.Sigma![a, b] * fit.DrawWeights![i]).Sum();
            Assert.Equal(expected, covariance.Mean![a, b], 11);
        }
        var across = result.JointBands(acrossCases: true);
        double counted = 0;
        for (int d = 0; d < result.Draws.Count; d++)
        {
            bool inside = true;
            for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++) inside &= result.Predictive[d, i, k] >= across.Lower[i, k] && result.Predictive[d, i, k] <= across.Upper[i, k];
            if (inside) counted += result.Weights[d];
        }
        Assert.InRange(counted, 0.947, 0.953);
        var selected = BassMultiOutput.Predict(fit, data, Cases(), new BassPredictOptions { Draws = [2, 0, 2], Predictive = false });
        Assert.Equal(new[] { 2, 0, 2 }, selected.Draws); Assert.Null(selected.Predictive);
        Assert.Throws<InvalidOperationException>(() => selected.JointBands());
        Assert.Throws<ArgumentException>(() => BassMultiOutput.Predict(fit, MultiOutputCoefficientTests.Data(3), Cases()));
    }

    [Fact]
    public void ProjectedSobolAndCoActivePairsHaveTheKnownPolynomialLaw()
    {
        var (fit, data) = Fixed();
        double[] a = [1, -0.4], b = [0.3, 1];
        var sobol = BassMultiOutput.Sensitivity(fit, data, a);
        var exact = BassMultiOutput.CoActive(fit, data, new CoActiveSpec { Output = a, OtherOutput = b, Route = ActiveSubspaceRoute.Exact });
        for (int d = 0; d < fit.Draws[0].Count; d++) Assert.Equal(fit.DrawWeights![d], exact.Posterior.Weights[d], 14);
        for (int d = 0; d < fit.Draws[0].Count; d++)
        {
            double[] af = BassMultiOutput.Project(fit.Draws[0][d], a), ag = BassMultiOutput.Project(fit.Draws[0][d], b);
            double variance = af[1] * af[1] + af[2] * af[2] + af[3] * af[3];
            Assert.InRange(Math.Abs(sobol.Decompositions[0].Variance[d] - variance), 0, 1e-8);
            var decomposition = sobol.Decompositions[0];
            for (int subset = 0; subset < decomposition.Subsets.Count; subset++)
            {
                var variables = decomposition.Subsets[subset];
                double expected = variables.Length == 2 ? af[3] * af[3] : variables[0] == 0 ? af[1] * af[1] : af[2] * af[2];
                Assert.InRange(Math.Abs(decomposition.MainVariance[d, subset] - expected), 0, 1e-8);
            }
            var cfg = exact.Posterior.Pairs[d].Cfg;
            Assert.InRange(Math.Abs(cfg[0, 0] - 12 * (af[1] * ag[1] + af[3] * ag[3])), 0, 1e-7);
            Assert.InRange(Math.Abs(cfg[1, 1] - (60 * af[2] * ag[2] + 12 * af[3] * ag[3])), 0, 1e-7);
            Assert.InRange(Math.Abs(cfg[0, 1] - 12 * Math.Sqrt(5) * af[3] * ag[2]), 0, 1e-7);
        }
        var sampled = BassMultiOutput.CoActive(fit, data, new CoActiveSpec { Output = a, OtherOutput = b, Route = ActiveSubspaceRoute.Sampled, Points = 40_000, Seed = 20261219, Draws = [0, 1] });
        for (int d = 0; d < 2; d++) for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
            Assert.InRange(Math.Abs(sampled.Posterior.Pairs[d].Cfg[i, j] - exact.Posterior.Pairs[d].Cfg[i, j]), 0, 5 * sampled.MonteCarloError![d].Cfg[i, j] + 1e-5);
        // Scalar reader fixtures of the same retained surfaces: the structure and coefficients are held, not refitted.
        foreach (double[] functional in new[] { new double[] { 1, 0 }, new double[] { 0, 1 }, a })
        {
            var response = (MultiOutputResponse)data.Response;
            double[] y = [.. Enumerable.Range(0, response.Cases).Select(i => functional[0] * response[i, 0] + functional[1] * response[i, 1])];
            var scalarData = new BassData(data.Predictors, new ScalarResponse(y));
            var scalarConfig = fit.Config with { OutputCovariance = null, Sigma2 = Sigma2Prior.InverseGamma(3, 1) };
            var scalar = Bass.Fit(scalarConfig, scalarData, new BassRunOptions { MaxDegreeOfParallelism = 1 });
            var store = new DrawStore<BassDraw>();
            foreach (var draw in fit.Draws[0]) store.Add(draw with { Coefficients = BassMultiOutput.Project(draw, functional), OutputCoefficients = null, Sigma = null, Sigma2 = 1, OutputRates = null });
            scalar = scalar with { Terms = fit.Terms, Draws = [store], DrawWeights = fit.DrawWeights };
            var expected = scalar.Sensitivity(scalarData).Decompositions[0];
            var projected = BassMultiOutput.Sensitivity(fit, data, functional).Decompositions[0];
            Assert.Equal(expected.MainVariance.Cast<double>(), projected.MainVariance.Cast<double>());
            Assert.Equal(expected.TotalVariance.Cast<double>(), projected.TotalVariance.Cast<double>());
        }
    }

    [Fact]
    public void JointGaussianProbabilityAndBandsAgreeWithIndependentObservations()
    {
        var (fit, data) = Fixed(count: 1);
        const int count = 60_000; const double rho = 0.8;
        var draw = fit.Draws[0][0] with { Sigma = new double[,] { { 1, rho }, { rho, 1 } }, OutputCoefficients = new double[4, 2] };
        var store = new DrawStore<BassDraw>(); for (int i = 0; i < count; i++) store.Add(draw);
        fit = fit with { Draws = [store], DrawWeights = null };
        var cases = new PredictionCases([new ContinuousColumn("x", [0.5]), new ContinuousColumn("z", [0.5])]);
        var prediction = BassMultiOutput.Predict(fit, data, cases, new BassPredictOptions { Seed = 20261220 });
        double probability = prediction.SpecificationProbability(new double[1, 2], new double[,] { { double.PositiveInfinity, double.PositiveInfinity } })[0];
        double expected = 0.25 + Math.Asin(rho) / (2 * Math.PI);
        Assert.InRange(Math.Abs(probability - expected), 0, 5 * Math.Sqrt(expected * (1 - expected) / count));
        var band = prediction.JointBands();
        var rng = new Xoshiro256PlusPlus(20261221); int inside = 0;
        for (int i = 0; i < count; i++)
        {
            double x = Normal.SampleStandard(rng), y = rho * x + Math.Sqrt(1 - rho * rho) * Normal.SampleStandard(rng);
            if (x >= band.Lower[0, 0] && x <= band.Upper[0, 0] && y >= band.Lower[0, 1] && y <= band.Upper[0, 1]) inside++;
        }
        output.WriteLine($"Gaussian positive quadrant {probability:G6}, exact {expected:G6}; independent joint-band coverage {(double)inside / count:G6}");
        Assert.InRange(Math.Abs((double)inside / count - 0.95), 0, 0.006);
        var quantile = prediction.Quantile(0.975);
        Assert.InRange(Math.Abs(quantile[0, 0] - Normal.Standard.Quantile(0.975)), 0, 2e-8);
    }
}
