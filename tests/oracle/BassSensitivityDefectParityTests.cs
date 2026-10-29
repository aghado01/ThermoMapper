using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Oracle.Tests;

/// <summary>
/// O-O2's corrected expectations (<c>bass_sensitivity_defects_oracle.R</c>), each asserted on the port with R's defective
/// value beside it:
/// <list type="bullet">
/// <item><b>Category weights are read:</b> the moments are sums of the supplied probabilities, where R's are subset sizes
/// over level counts on both paths;</item>
/// <item><b>truncation bounds are clipped:</b> hinge products under U(0.25, 0.75) match quadrature, where R's basis path
/// integrates over the training range;</item>
/// <item><b>effects use the declared measure and degree:</b> the mean's term moments match quadrature under each measure,
/// where R's effects take the U(0, 1) moment, and a degree-2 term's effect curve is the degree-2 hinge, where R's is the
/// degree-1 one;</item>
/// <item><b>basis sensitivity is in output units</b> (O7): on R's mock two-component fit, each location's variance and
/// effects are the physical reference, y_s² times R's decomposition-unit values, and the indices match the unit-free
/// reference.</item>
/// </list>
/// </summary>
public sealed class BassSensitivityDefectParityTests
{
    private static readonly string BassSource = Environment.GetEnvironmentVariable("BASS_R_SOURCE") ?? "D:/aghado01/packages/BASS";

    private static readonly Lazy<JsonElement?> Fixture = new(() =>
    {
        if (!ROracle.IsAvailable || !Directory.Exists(BassSource)) return null;
        // The check script's basis: y_m, y_s, then two basis columns at 12 locations.
        var text = new StringBuilder();
        for (int r = 0; r < Locations; r++)
            text.AppendLine(string.Join(",", BasisRow(r).Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        string csv = Path.Combine(Path.GetTempPath(), $"bass_sensitivity_{Guid.NewGuid():N}.csv");
        File.WriteAllText(csv, text.ToString());
        try
        {
            return ROracle.Run("oracles/bass_sensitivity_defects_oracle.R", csv);
        }
        finally
        {
            File.Delete(csv);
        }
    });

    // x holds every knot the fixture uses, with ends 0 and 1, so declared units are the fit's.
    private static readonly double[] X = [0.0, 0.2, 0.3, 0.4, 0.6, 1.0];

    [Fact]
    public void CategoryWeights_AreRead()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement cw = r.GetProperty("category_weights");
        double[] rSizes = Doubles(cw.GetProperty("r_direct").GetProperty("sub_cnt"));
        Assert.Equal(new[] { 1.0 / 3, 2.0 / 3 }, rSizes);
        string[] levels = ["a", "b", "c"];
        var data = new BassData([new CategoricalPredictor("z", [0, 1, 2, 0, 1, 2], levels)], new ScalarResponse([1.0, 2.0, 3.0, 1.5, 2.5, 3.5]));
        BassTerm d1 = new(categorical: [new CategoricalFactor(0, [0])]), d2 = new(categorical: [new CategoricalFactor(0, [0, 1])]);
        foreach (JsonElement c in cw.GetProperty("cases").EnumerateArray())
        {
            double[] w = Doubles(c.GetProperty("weights"));
            JsonElement e = c.GetProperty("enumerated");
            BassMoments moments = Moments(data, new() { ["z"] = new CategoricalMeasure(levels.Zip(w).ToDictionary(p => p.First, p => p.Second)) });
            Assert.Equal(e.GetProperty("P_D1").GetDouble(), moments.Mean(d1), 15);
            Assert.Equal(e.GetProperty("P_D2").GetDouble(), moments.Mean(d2), 15);
            Assert.Equal(e.GetProperty("P_D1_and_D2").GetDouble(), moments.Cross(d1, d2), 15);
            Assert.NotEqual(rSizes[0], moments.Mean(d1), 6);
        }
        // Two levels, P(a) = 0.9: E[1(z ∈ {a})²] is 0.9, where R gives 0.5.
        JsonElement two = cw.GetProperty("two_level");
        var pair = new BassData([new CategoricalPredictor("z", [0, 1, 0, 1], ["a", "b"])], new ScalarResponse([1.0, 2.0, 1.5, 2.5]));
        BassMoments m2 = Moments(pair, new() { ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.9, ["b"] = 0.1 }) });
        BassTerm a = new(categorical: [new CategoricalFactor(0, [0])]);
        Assert.Equal(two.GetProperty("enumerated").GetDouble(), m2.Cross(a, a), 15);
        Assert.Equal(0.5, two.GetProperty("r_second_moment").GetDouble());
    }

    [Fact]
    public void TruncationBounds_AreClipped()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement tb = r.GetProperty("truncation_bounds");
        double[] bounds = Doubles(tb.GetProperty("measure"));
        BassMoments moments = Moments(Continuous(), new() { ["x"] = new UniformMeasure(bounds[0], bounds[1]) });
        int rWrong = 0;
        foreach (JsonElement c in tb.GetProperty("cases").EnumerateArray())
        {
            double[] t = Doubles(c.GetProperty("t")), s = Doubles(c.GetProperty("s"));
            // On R's scale: each hinge divided by its largest value on [0, 1], as the port's terms are.
            BassTerm h1 = Hinge(t[0], (int)s[0]), h2 = Hinge(t[1], (int)s[1]);
            double reference = c.GetProperty("quadrature").GetDouble();
            Assert.Equal(reference, moments.Cross(h1, h2), 14);
            Assert.Equal(reference, c.GetProperty("r_direct").GetDouble(), 14);
            if (Math.Abs(c.GetProperty("r_basis").GetDouble() - reference) > 1e-6) rWrong++;
        }
        Assert.Equal(5, rWrong);   // R's basis path, 0.45 for the inventory's 0.1105903 among them
    }

    [Fact]
    public void EffectMoments_UseTheDeclaredMeasure()
    {
        if (Fixture.Value is not { } r) return;
        foreach (JsonElement c in r.GetProperty("effects").GetProperty("first_moment").EnumerateArray())
        {
            double[] trunc = Doubles(c.GetProperty("trunc"));
            InputMeasure measure = c.GetProperty("measure").GetString() == "uniform"
                ? new UniformMeasure(trunc[0], trunc[1])
                : new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.1), new NormalComponent(0.6, 0.7, 0.2)]) { Lower = trunc[0], Upper = trunc[1] };
            BassMoments moments = Moments(Continuous(), new() { ["x"] = measure });
            BassTerm term = Hinge(c.GetProperty("t").GetDouble(), c.GetProperty("s").GetInt32());
            double reference = c.GetProperty("quadrature").GetDouble();
            Assert.Equal(reference, moments.Mean(term), 14);
            // R's effects read C1.temp, the U(0, 1) moment, whatever the measure.
            bool unit = trunc[0] == 0.0 && trunc[1] == 1.0 && c.GetProperty("measure").GetString() == "uniform";
            Assert.Equal(unit, Math.Abs(c.GetProperty("r_c1_temp").GetDouble() - reference) < 1e-12);
        }
    }

    [Fact]
    public void EffectCurves_UseTheDeclaredDegree()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement curve = r.GetProperty("effects").GetProperty("degree2_curve");
        Assert.Equal(0.0, curve.GetProperty("r_max_abs_vs_degree1").GetDouble());          // R's degree-2 curve is the degree-1 hinge
        Assert.All(r.GetProperty("effects").GetProperty("makebasis_degree_args").EnumerateArray(), a => Assert.Equal("1", a.GetString()));

        // A degree-2 fit holding one term, (x − 0.3)₊²/0.7², coefficient 1: its main effect under U(0, 1) is that hinge less
        // its mean 0.7³/3/0.7².
        BassData data = Continuous();
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { Degree = 2, MaxBasis = 2 },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 5, BurnIn = 20, Samples = 2 },
        };
        BassFit fit = Bass.Fit(config, data);
        var table = new BassTermTable();
        int id = table.Intern(Hinge(0.3, 1));
        var store = new DrawStore<BassDraw>();
        store.Add(new BassDraw(new BassStructure([id]), 1.0, 1.0, [0.0, 1.0], 1.0));
        EffectFunctions e = BassSensitivity.Effects(fit with { Terms = table, Draws = [store] }, data).At[0];
        Effect main = e.Effects[0];
        double gap1 = 0.0;
        for (int c = 0; c < main.Grids[0].Values.Length; c++)
        {
            double x = main.Grids[0].Values[c];
            double expected = (Math.Pow(Math.Max(x - 0.3, 0.0), 2) - Math.Pow(0.7, 3) / 3) / 0.49;
            Assert.Equal(expected, main.Values[0, c], 1e-14);
            gap1 = Math.Max(gap1, Math.Abs(main.Values[0, c] - (Math.Max(x - 0.3, 0.0) - 0.245) / 0.7));
        }
        Assert.True(gap1 > 0.1, $"the degree-2 effect sits {gap1} from the degree-1 one");
        Assert.Equal(0.7 * 0.7 * 0.7 / 3 / 0.49, e.Mean[0], 1e-14);
    }

    [Fact]
    public void DecompositionUnits_AreThePhysicalOnes()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement du = r.GetProperty("decomposition_units");
        double[] ys = Doubles(du.GetProperty("y_s"));
        double[] physical = Doubles(du.GetProperty("reference_var_total_physical"));
        double[] decomposition = Doubles(du.GetProperty("reference_var_total_decomposition"));
        double[] rTotal = [.. du.GetProperty("r_var_total").EnumerateArray().Select(v => v[0].GetDouble())];
        double[][] reference = Rows(du.GetProperty("reference_effects_decomposition"));
        double[][] rEffects = Rows(du.GetProperty("r_effects"));
        double[][] rIndices = Rows(du.GetProperty("r_indices"));
        Assert.Equal(new[] { "1", "2", "1x2", "other" }, du.GetProperty("names").EnumerateArray().Select(v => v.GetString()));

        // R's mock: on x1, x2 with R's scaled degree-1 hinges and the default uniform measures,
        // w_1 = 0.5 + 1.2·h(x1; +, 0.3) − 0.8·h(x2; −, 0.6) and w_2 = −0.2 + 0.9·h(x1; +, 0.3)·h(x2; +, 0.4) + 0.5·h(x2; +, 0.2).
        BasisSobolAnalysis analysis = UnitsFit().Sensitivity();
        Assert.Equal(new[] { new[] { 0 }, new[] { 1 }, new[] { 0, 1 } }, analysis.Subsets);
        for (int l = 0; l < Locations; l++)
        {
            SobolDecomposition at = analysis.AtLocation(l);
            double v = at.Variance[0], scale = ys[l] * ys[l];
            Assert.Equal(physical[l], v, 1e-13 * physical[l]);
            // R's variances are the decomposition-unit ones, a factor y_s² below.
            Assert.Equal(decomposition[l], rTotal[l], 1e-13 * decomposition[l]);
            Assert.Equal(v, scale * rTotal[l], 1e-13 * v);
            for (int u = 0; u < 3; u++)
            {
                Assert.Equal(scale * reference[l][u], at.MainVariance[0, u], 1e-13 * v);
                Assert.Equal(reference[l][u], rEffects[l][u], 1e-13 * decomposition[l]);
                // The indices are unit-free: the port's and R's both match the reference.
                Assert.Equal(reference[l][u] / decomposition[l], at.MainIndex(0, u), 1e-13);
                Assert.Equal(reference[l][u] / decomposition[l], rIndices[l][u], 1e-13);
            }
            Assert.Equal(v, at.TotalVariance[0, 0] + reference[l][1] * scale, 1e-13 * v);
        }
        Assert.True(du.GetProperty("max_rel_r_vs_physical").GetDouble() > 0.9);
    }

    private const int Locations = 12;

    // Location r's y_m, y_s and basis row: t = r/11, y_m = 10 + t, y_s = 2 + 3t, K = (sin πt, cos πt).
    private static double[] BasisRow(int r)
    {
        double t = r / 11.0;
        return [10 + t, 2 + 3 * t, Math.Sin(Math.PI * t), Math.Cos(Math.PI * t)];
    }

    // The fixture's mock as a basis fit: each component run briefly on its coefficient response, then given R's terms and
    // one draw; x2 holds X reversed, so both inputs have X's knots.
    private static BassBasisFit UnitsFit()
    {
        var k = new double[Locations, 2];
        var ym = new double[Locations];
        var ys = new double[Locations];
        for (int l = 0; l < Locations; l++)
        {
            double[] row = BasisRow(l);
            (ym[l], ys[l], k[l, 0], k[l, 1]) = (row[0], row[1], row[2], row[3]);
        }
        BassPredictor[] predictors = [new ContinuousPredictor("x1", X), new ContinuousPredictor("x2", [.. X.Reverse()])];
        int n = X.Length;
        var w = new double[2, n];
        for (int c = 0; c < 2; c++)
            for (int i = 0; i < n; i++) w[c, i] = Math.Sin(1.0 + i + 2.0 * c);
        BassDecomposition decomposition = BassDecomposition.FromBasis(k, w, ym, ys, residuals: new double[n, Locations]);
        static HingeFactor H(int variable, int sign, double knot) => new(variable, sign, Array.IndexOf(X, knot));
        (BassTerm[] Terms, double[] Beta)[] components =
        [
            ([new([H(0, 1, 0.3)]), new([H(1, -1, 0.6)])], [0.5, 1.2, -0.8]),
            ([new([H(0, 1, 0.3), H(1, 1, 0.4)]), new([H(1, 1, 0.2)])], [-0.2, 0.9, 0.5]),
        ];
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 2 },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 5, BurnIn = 20, Samples = 2 },
        };
        var fits = new BassFit[2];
        for (int c = 0; c < 2; c++)
        {
            BassFit run = Bass.Fit(config, new BassData(predictors, new ScalarResponse(decomposition.CoefficientResponse(c))));
            var table = new BassTermTable();
            int[] ids = [.. components[c].Terms.Select(table.Intern)];
            var store = new DrawStore<BassDraw>();
            store.Add(new BassDraw(new BassStructure(ids), 1.0, 1.0, components[c].Beta, 1.0));
            fits[c] = run with { Terms = table, Draws = [store], DrawWeights = null };
        }
        return BassComponents.Assemble(decomposition, predictors, fits, 7);
    }

    private static double[][] Rows(JsonElement e) => [.. e.EnumerateArray().Select(Doubles)];

    private static BassData Continuous()
        => new([new ContinuousPredictor("x", X)], new ScalarResponse([.. X.Select(v => v * v)]));

    private static BassMoments Moments(BassData data, Dictionary<string, InputMeasure> measures)
    {
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        return BassMoments.On(inputs, 1, ProductMeasure.Resolve(transform, inputs, measures));
    }

    private static BassTerm Hinge(double knot, int sign) => new([new HingeFactor(0, sign, Array.IndexOf(X, knot))]);

    private static double[] Doubles(JsonElement e) => [.. e.EnumerateArray().Select(v => v.GetDouble())];
}
