using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's one-variable hinge moments (O5) against O-O3's quadrature references (<c>bass_moments_oracle.R</c>): all 156
/// cases, first moments with knots inside and outside the support and cross moments with equal and distinct knots and every
/// sign pair, at degrees 1–3, under the oracle's 9 measures, each declared through the sensitivity records and resolved as
/// <see cref="ProductMeasure"/> does: uniform on [0, 1] and on [0.25, 0.75], a normal mixture under both truncation
/// conventions, a narrowly truncated normal, Student ν = 1, 2, 3 and a ν = (1, 3) mixture. R's own Student moments are
/// beside them: with a working 2F1, degrees 2 and 3 of R's cross moment are its degree-1 value, and the port does not
/// follow them. Categorical moments against enumeration, a zero-probability level included.
/// </summary>
public sealed class BassMomentParityTests(ITestOutputHelper output)
{
    private static readonly string BassSource = Environment.GetEnvironmentVariable("BASS_R_SOURCE") ?? "D:/aghado01/packages/BASS";

    private static readonly (string Name, InputMeasure Measure)[] Measures =
    [
        ("uniform_unit", new UniformMeasure()),
        ("uniform_inner", new UniformMeasure(0.25, 0.75)),
        ("normal_mix", new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.1), new NormalComponent(0.6, 0.7, 0.2)])),
        ("normal_mix_conditional", new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.1), new NormalComponent(0.6, 0.7, 0.2)])
        {
            Weighting = MixtureWeighting.ConditionalOnTruncation,
        }),
        ("normal_narrow", new MixtureMeasure([new NormalComponent(1.0, 0.5, 0.3)]) { Lower = 0.4, Upper = 0.6 }),
        ("student_nu1", new MixtureMeasure([new StudentComponent(1.0, 1.0, 0.5, 0.15)])),
        ("student_nu2", new MixtureMeasure([new StudentComponent(1.0, 2.0, 0.5, 0.15)])),
        ("student_nu3", new MixtureMeasure([new StudentComponent(1.0, 3.0, 0.5, 0.15)])),
        ("student_mix", new MixtureMeasure([new StudentComponent(0.5, 1.0, 0.3, 0.1), new StudentComponent(0.5, 3.0, 0.7, 0.2)])),
    ];

    [Fact]
    public void HingeMoments_MatchQuadrature_UnderEveryMeasure()
    {
        if (!ROracle.IsAvailable || !Directory.Exists(BassSource)) return;
        JsonElement r = Run(Cases(), "continuous");
        string[] names = [.. r.GetProperty("measures").EnumerateArray().Select(m => m.GetProperty("name").GetString()!)];
        Assert.Equal(Measures.Select(m => m.Name), names);
        double[][] cases = Rows(r.GetProperty("cases"));
        double?[][] quadrature = NullableRows(r.GetProperty("quadrature")), repaired = NullableRows(r.GetProperty("r_repaired"));
        Assert.Equal(156, cases.Length);

        // One continuous predictor spanning [0, 1] exactly, so the declared units are the fit's.
        var data = new BassData([new ContinuousPredictor("x", [0.0, 0.25, 0.5, 0.75, 1.0])], new ScalarResponse([0.0, 1.0, 0.0, 1.0, 0.0]));
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        int wrongInR = 0;
        for (int j = 0; j < Measures.Length; j++)
        {
            ContinuousMeasure measure = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["x"] = Measures[j].Measure }).Continuous(0)!;
            double worst = 0.0, largestError = 0.0;
            for (int i = 0; i < cases.Length; i++)
            {
                double[] c = cases[i];
                int q = (int)c[4];
                (int, double)[] hinges = c[3] == 0.0 ? [((int)c[1], c[0])] : [((int)c[1], c[0]), ((int)c[3], c[2])];
                (double value, double error) = PiecewisePolynomial.HingeProduct(q, hinges).Expectation(measure.Measure, measure.Lower, measure.Upper);
                double reference = quadrature[i][j]!.Value;
                double gap = Math.Abs(value - reference);
                Assert.True(gap <= 1e-12 * Math.Abs(reference) + 1e-16, $"{names[j]}, case {i} ({string.Join(", ", c)}): {value} against {reference}");
                worst = Math.Max(worst, reference == 0.0 ? gap : gap / Math.Abs(reference));
                largestError = Math.Max(largestError, error);
                if (repaired[i][j] is { } rValue && c[3] != 0.0 && q > 1 && reference != 0.0 && Math.Abs(rValue - reference) > 0.5 * Math.Abs(reference)) wrongInR++;
            }
            output.WriteLine($"{names[j],-24} largest relative gap {worst:E1}, largest quadrature estimate {largestError:E1}");
            // Student components report their Gauss–Legendre estimate; closed forms report none.
            if (names[j].StartsWith("student", StringComparison.Ordinal)) Assert.InRange(largestError, 1e-18, 1e-12);
            else Assert.Equal(0.0, largestError);
        }
        // R's repaired Student forms at degrees 2–3 return the degree-1 cross moment: the port matches quadrature instead.
        output.WriteLine($"R's repaired Student cross moments off by over half at degrees 2-3: {wrongInR} cases");
        Assert.True(wrongInR > 0);
    }

    [Fact]
    public void CategoricalMoments_MatchEnumeration()
    {
        if (!ROracle.IsAvailable) return;
        double[] p = [0.0, 0.2, 0.3, 0.5];
        int[][] subsets = [[0], [0, 1], [1, 2], [2, 3]];
        var rows = new List<double[]> { p };
        rows.AddRange(subsets.Select(s => Enumerable.Range(0, 4).Select(l => s.Contains(l) ? 1.0 : 0.0).ToArray()));
        JsonElement r = Run(rows, "categorical");
        double[] enumerated = [.. r.GetProperty("P").EnumerateArray().Select(v => v.GetDouble())];
        double[][] joint = Rows(r.GetProperty("P_joint"));
        double[] rRatios = [.. r.GetProperty("r_P").EnumerateArray().Select(v => v.GetDouble())];

        string[] levels = ["a", "b", "c", "d"];
        var data = new BassData([new CategoricalPredictor("z", [0, 1, 2, 3, 0, 1], levels)], new ScalarResponse([0.0, 1.0, 2.0, 3.0, 0.5, 1.5]));
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var measure = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure>
        {
            ["z"] = new CategoricalMeasure(levels.Select((l, k) => (l, p[k])).ToDictionary(e => e.l, e => e.Item2)),
        });
        var moments = BassMoments.On(inputs, 1, measure);
        BassTerm[] terms = [.. subsets.Select(s => new BassTerm(categorical: [new CategoricalFactor(0, s)]))];
        for (int i = 0; i < terms.Length; i++)
        {
            Assert.Equal(enumerated[i], moments.Mean(terms[i]), 15);
            for (int k = 0; k < terms.Length; k++) Assert.Equal(joint[i][k], moments.Cross(terms[i], terms[k]), 15);
        }
        // R's subset size over level count ignores the probabilities: 0.25 for the zero-probability level.
        Assert.Equal(0.0, moments.Mean(terms[0]));
        Assert.Equal(0.25, rRatios[0], 15);
    }

    // O-O3's cases, as r/scripts/bass_stage1_checks.R builds them: (t1, s1, t2, s2, q), s2 = 0 for a first moment.
    private static List<double[]> Cases()
    {
        var rows = new List<double[]>();
        foreach (int q in new[] { 1, 2, 3 })
            foreach (int s in new[] { -1, 1 })
                foreach (double t in new[] { -0.1, 0.1, 0.3, 0.5, 0.8, 1.1 }) rows.Add([t, s, 0, 0, q]);
        double[] kn = [0.1, 0.3, 0.5, 0.8];
        foreach (int q in new[] { 1, 2, 3 })
            for (int j = 0; j < 4; j++)
                for (int i = 0; i <= j; i++)
                    foreach ((int s1, int s2) in new[] { (1, 1), (1, -1), (-1, 1), (-1, -1) }) rows.Add([kn[i], s1, kn[j], s2, q]);
        return rows;
    }

    private static JsonElement Run(List<double[]> rows, string mode)
    {
        string csv = Path.Combine(Path.GetTempPath(), $"bass_moments_{Guid.NewGuid():N}.csv");
        var text = new StringBuilder();
        foreach (double[] row in rows) text.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(csv, text.ToString());
        try
        {
            return ROracle.Run("oracles/bass_moments_oracle.R", csv, mode);
        }
        finally
        {
            File.Delete(csv);
        }
    }

    private static double[][] Rows(JsonElement e) => [.. e.EnumerateArray().Select(row => row.EnumerateArray().Select(v => v.GetDouble()).ToArray())];

    private static double?[][] NullableRows(JsonElement e)
        => [.. e.EnumerateArray().Select(row => row.EnumerateArray().Select(v => v.ValueKind == JsonValueKind.Number ? v.GetDouble() : (double?)null).ToArray())];
}
