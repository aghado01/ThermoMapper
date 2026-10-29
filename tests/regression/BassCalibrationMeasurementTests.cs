using System;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Calibration's declared measurement-error model (stage 7, C-X1; D09): Gaussian variance groups bound by output index
/// and validated as a partition; each group's prior as a density on σ², the inverse gamma against its closed form and
/// the half-Cauchy on σ (Gelman 2006) normalized by quadrature in log σ² and scale-equivariant; the half-Cauchy
/// distribution's tails and draws; and the config's JSON round trip and refusals.
/// </summary>
public sealed class BassCalibrationMeasurementTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationMeasurementTests(ITestOutputHelper output) => _out = output;

    private static GaussianMeasurement TwoGroups(VariancePrior first, VariancePrior second)
        => new([new VarianceGroup([0, 2, 4], first), new VarianceGroup([1, 3], second)]);

    // ∫ p(σ²) dσ² by the trapezoid rule in u = ln σ², ∫ p(eᵘ) eᵘ du, on [−60, 60]: both priors' tails are below e⁻²⁵ there.
    private static double Mass(Func<double, double> logDensity)
    {
        const double lower = -60, upper = 60, h = 0.005;
        int n = (int)((upper - lower) / h);
        double sum = 0.0;
        for (int i = 0; i <= n; i++)
        {
            double u = lower + i * h;
            double v = Math.Exp(logDensity(Math.Exp(u)) + u);
            sum += (i == 0 || i == n ? 0.5 : 1.0) * v;
        }
        return sum * h;
    }

    [Fact]
    public void Groups_AreBoundByOutputIndex()
    {
        var model = new MeasurementModel(TwoGroups(new InverseGammaVariance(10, 0.1), new HalfCauchyDeviation(0.2)), 5);
        Assert.Equal(2, model.Groups);
        Assert.Equal([0, 1, 0, 1, 0], Enumerable.Range(0, 5).Select(model.GroupOf));
        Assert.Equal([1, 3], model.OutputsOf(1));
        var perOutput = new double[5];
        model.OutputVariances([0.3, 0.7], perOutput);
        Assert.Equal([0.3, 0.7, 0.3, 0.7, 0.3], perOutput);
        Assert.True(model.IsProper);
        Assert.Equal(model.LogPrior(0, 0.3) + model.LogPrior(1, 0.7), model.LogPrior([0.3, 0.7]));
    }

    [Fact]
    public void Declarations_ThatAreNotAPartitionOrHaveBadPriors_AreRefused()
    {
        VariancePrior ig = new InverseGammaVariance(10, 0.1);
        Assert.Contains("more than one", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0, 1], ig), new VarianceGroup([1, 2], ig)]), 3)).Message);
        Assert.Contains("in no variance group", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0, 2], ig)]), 3)).Message);
        Assert.Contains("there are 3 outputs", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0, 1, 2, 3], ig)]), 3)).Message);
        Assert.Contains("indexed from 0", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([-1, 0, 1, 2], ig)]), 3)).Message);
        Assert.Contains("names no output", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([], ig), new VarianceGroup([0, 1, 2], ig)]), 3)).Message);
        Assert.Contains("at least one", Assert.Throws<ArgumentException>(() => new MeasurementModel(new GaussianMeasurement([]), 3)).Message);
        // R's s2.df = 0 switches to an improper kernel; here an inverse-gamma group needs ν > 0 and is pointed at the half-Cauchy.
        Assert.Contains("half-Cauchy", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0, 1, 2], new InverseGammaVariance(0, 0.1))]), 3)).Message);
        Assert.Contains("scale", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0, 1, 2], new InverseGammaVariance(4, 0))]), 3)).Message);
        Assert.Contains("scale", Assert.Throws<ArgumentException>(() =>
            new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0, 1, 2], new HalfCauchyDeviation(double.NaN))]), 3)).Message);
    }

    [Fact]
    public void InverseGamma_IsRsKernel_Normalized()
    {
        const double nu = 10, s = 0.1;
        var model = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new InverseGammaVariance(nu, s))]), 1);
        double a = nu / 2, b = a * s * s;
        foreach (double sigma2 in new[] { 1e-4, 0.003, 0.01, 0.2, 5.0 })
        {
            double closed = a * Math.Log(b) - SpecialFunctions.LogGamma(a) - (a + 1) * Math.Log(sigma2) - b / sigma2;
            Assert.Equal(closed, model.LogPrior(0, sigma2), 1e-12);
        }
        Assert.Equal(1.0, Mass(x => model.LogPrior(0, x)), 1e-9);
    }

    [Fact]
    public void HalfCauchy_IsNormalized_AndScaleEquivariant()
    {
        foreach (double s in new[] { 1e-3, 0.2, 40.0 })
        {
            var model = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new HalfCauchyDeviation(s))]), 1);
            double mass = Mass(x => model.LogPrior(0, x));
            _out.WriteLine($"half-Cauchy, scale {s}: mass {mass:F12}");
            Assert.Equal(1.0, mass, 1e-9);
        }
        // σ ~ C⁺(0, s) ⇒ cσ ~ C⁺(0, cs): p_{cs}(c²σ²)·c² = p_s(σ²), and the inverse gamma's (ν, ŝ) ↦ (ν, cŝ) likewise.
        foreach (double c in new[] { 0.01, 3.0, 250.0 })
        {
            var hc = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new HalfCauchyDeviation(0.3))]), 1);
            var hcScaled = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new HalfCauchyDeviation(0.3 * c))]), 1);
            var ig = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new InverseGammaVariance(6, 0.3))]), 1);
            var igScaled = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new InverseGammaVariance(6, 0.3 * c))]), 1);
            foreach (double sigma2 in new[] { 1e-3, 0.05, 0.4, 9.0 })
            {
                Assert.Equal(hc.LogPrior(0, sigma2), hcScaled.LogPrior(0, c * c * sigma2) + Math.Log(c * c), 1e-11);
                Assert.Equal(ig.LogPrior(0, sigma2), igScaled.LogPrior(0, c * c * sigma2) + Math.Log(c * c), 1e-10);
            }
        }
    }

    [Fact]
    public void HalfCauchyDistribution_TailsQuantilesAndDraws()
    {
        HalfCauchy d = HalfCauchy.WithScale(2.5);
        Assert.Equal(2.5, d.Median);
        Assert.Equal(0.5, d.Cdf(2.5), 1e-15);
        foreach (double p in new[] { 1e-9, 0.1, 0.5, 0.9, 1 - 1e-9 })
        {
            double x = d.Quantile(p);
            Assert.Equal(p, d.Cdf(x), 1e-12);
            Assert.Equal(1 - p, d.Survival(x), 1e-12 * Math.Max(1, (1 - p) * 1e3));
        }
        // The far upper tail keeps its digits: S(x) ≈ 2s/(πx).
        Assert.Equal(2 * 2.5 / (Math.PI * 1e12), d.Survival(1e12), 1e-24);
        Assert.Equal(d.Cdf(3.0) - d.Cdf(1.0), d.IntervalProbability(1.0, 3.0), 1e-14);
        Assert.Equal(double.NegativeInfinity, d.LogPdf(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => HalfCauchy.WithScale(0));

        var rng = new Xoshiro256PlusPlus(21);
        const int n = 40000;
        double[] draws = [.. Enumerable.Range(0, n).Select(_ => d.Sample(rng))];
        foreach (double p in new[] { 0.1, 0.5, 0.9 })
        {
            double below = draws.Count(x => x <= d.Quantile(p)) / (double)n;
            Assert.InRange(below, p - 4 * Math.Sqrt(p * (1 - p) / n), p + 4 * Math.Sqrt(p * (1 - p) / n));
        }
        // σ²'s draws from the variance prior: P(σ² ≤ ŝ²) = ½ under the half-Cauchy, the inverse gamma's CDF there under it.
        var model = new MeasurementModel(new GaussianMeasurement([new VarianceGroup([0], new HalfCauchyDeviation(0.4)), new VarianceGroup([1], new InverseGammaVariance(8, 0.4))]), 2);
        var v = new double[2];
        int[] under = new int[2];
        for (int i = 0; i < n; i++)
        {
            model.SamplePrior(rng, v);
            if (v[0] <= 0.16) under[0]++;
            if (v[1] <= 0.16) under[1]++;
        }
        double[] expected = [0.5, InverseGamma.ShapeScale(4, 4 * 0.16).Cdf(0.16)];
        for (int g = 0; g < 2; g++)
            Assert.InRange(under[g] / (double)n, expected[g] - 4 * Math.Sqrt(0.25 / n), expected[g] + 4 * Math.Sqrt(0.25 / n));
    }

    [Fact]
    public void Config_RoundTripsThroughJson_AndIsValidated()
    {
        var config = new CalibrationConfig
        {
            Parameters = [new ContinuousParameter("theta", 10, 20), new CategoricalParameter("mode", ["a", "b", "c"])],
            Measurement = TwoGroups(new InverseGammaVariance(10, 0.1), new HalfCauchyDeviation(0.2)),
        };
        string json = config.ToJson();
        CalibrationConfig back = CalibrationConfig.FromJson(json);
        Assert.Equal(json, back.ToJson());
        Assert.IsType<HalfCauchyDeviation>(((GaussianMeasurement)back.Measurement).Groups[1].Prior);
        _out.WriteLine(json);

        Assert.Throws<JsonException>(() => CalibrationConfig.FromJson(json.Replace("\"Parameters\"", "\"Parameter\"")));
        Assert.Throws<JsonException>(() => CalibrationConfig.FromJson("{\"Parameters\":[]}"));
        Assert.Contains("CalibrationConfig.Measurement", Assert.Throws<ArgumentException>(() =>
            CalibrationConfig.FromJson(json.Replace("\"DegreesOfFreedom\":10", "\"DegreesOfFreedom\":0"))).Message);
        Assert.Contains("CalibrationConfig.Parameters", Assert.Throws<ArgumentException>(() =>
            (config with { Parameters = [new ContinuousParameter("theta", 2, 1)] }).Validate()).Message);
    }
}
