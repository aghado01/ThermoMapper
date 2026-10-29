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
/// degree-1 one.</item>
/// </list>
/// The fixture's fourth section, physical units, needs basis fits and waits for basis sensitivity (O7).
/// </summary>
public sealed class BassSensitivityDefectParityTests
{
    private static readonly string BassSource = Environment.GetEnvironmentVariable("BASS_R_SOURCE") ?? "D:/aghado01/packages/BASS";

    private static readonly Lazy<JsonElement?> Fixture = new(() =>
    {
        if (!ROracle.IsAvailable || !Directory.Exists(BassSource)) return null;
        // The check script's basis: y_m, y_s, then two basis columns at 12 locations.
        var text = new StringBuilder();
        for (int r = 0; r < 12; r++)
        {
            double t = r / 11.0;
            double[] row = [10 + t, 2 + 3 * t, Math.Sin(Math.PI * t), Math.Cos(Math.PI * t)];
            text.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        }
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
