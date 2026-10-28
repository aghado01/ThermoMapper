using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's admission angles (stage 4, S4; D06; F-O7) against <c>r/oracles/bass_admission_oracle.R</c>, which reads the
/// centred angles by QR of the explicitly centred design (<c>centred</c> mode) and, for direct fits, of the
/// materialized Khatri–Rao design (<c>hadamard</c> mode). The design layer forms its Gram without cancellation,
/// explicitly centred columns for scalar fits and the three-term Hadamard form for direct ones, each summed pairwise,
/// so its angles hold F-O7's accuracy: within 10⁻⁶ relative, on random terms and on F-O7's hard case, an indicator
/// true for all but one of 10⁵ cases, nearly collinear with a second column once both are centred.
/// </summary>
public sealed class BassAdmissionParityTests
{
    private readonly ITestOutputHelper _out;
    public BassAdmissionParityTests(ITestOutputHelper output) => _out = output;

    private static readonly BassAdmission Permissive = new(50, SupportKind.Count, 1, 1, 1e-14, null);

    private static string Csv(double[][] columns)
    {
        string path = Path.Combine(Path.GetTempPath(), $"bass_admission_{Guid.NewGuid():N}.csv");
        var text = new StringBuilder();
        for (int i = 0; i < columns[0].Length; i++)
            text.AppendLine(string.Join(",", columns.Select(c => c[i].ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(path, text.ToString());
        return path;
    }

    private static double[] Loo(JsonElement r) => [.. r.GetProperty("loo").EnumerateArray().Select(v => v.GetDouble())];

    private void Agree(BassDesign design, double[] reference, string label)
    {
        Assert.Equal(design.Count, reference.Length);
        double worst = 0.0;
        for (int j = 0; j < reference.Length; j++)
        {
            double s = design.Angle(j);
            worst = Math.Max(worst, Math.Abs(s - reference[j]) / reference[j]);
            Assert.True(Math.Abs(s - reference[j]) <= 1e-6 * reference[j], $"{label}, column {j}: {s:R} against R's QR {reference[j]:R}");
        }
        _out.WriteLine($"[admission] {label}: smallest angle {reference.Min():E3}, largest relative error {worst:E2}");
    }

    [Fact]
    public void ScalarAngles_MatchTheExplicitlyCentredQr()
    {
        if (!ROracle.IsAvailable) return;
        var rng = new Xoshiro256PlusPlus(seed: 41);
        const int n = 120;
        double[] x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] x2 = [.. Enumerable.Range(0, n).Select(_ => Math.Round(rng.NextDouble() * 20) / 20)];
        int[] z = [.. Enumerable.Range(0, n).Select(i => i % 4)];
        var data = new BassData([new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2), new CategoricalPredictor("z", z, ["a", "b", "c", "d"])],
                                new ScalarResponse([.. x1.Select((v, i) => v * x2[i])]));
        var basis = new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1);
        BassTerm[] terms =
        [
            new([new HingeFactor(0, 1, 30)]), new([new HingeFactor(0, 1, 33)]), new([new HingeFactor(0, -1, 70)]),
            new([new HingeFactor(0, 1, 40), new HingeFactor(1, -1, 12)]), new([new HingeFactor(1, 1, 5)], [new CategoricalFactor(0, [1, 2])]),
            new(categorical: [new CategoricalFactor(0, [0])]),
        ];
        var design = new BassDesign(basis, data, Permissive);
        Assert.Null(design.Load(terms));

        double[][] columns = [Enumerable.Repeat(1.0, n).ToArray(), .. design.Terms.Select(t => Case(basis, t))];
        string csv = Csv(columns);
        try { Agree(design, Loo(ROracle.Run("oracles/bass_admission_oracle.R", csv, "centred")), "random terms"); }
        finally { File.Delete(csv); }
    }

    [Fact]
    public void ANearConstantIndicator_KeepsItsAngle()
    {
        if (!ROracle.IsAvailable) return;
        // F-O7's hard case: an indicator true for all but one case in 10⁵ (ρ = 10⁻⁵), and a hinge that is nearly that one
        // case alone, so the two are nearly collinear once centred.
        const int n = 100_000;
        var rng = new Xoshiro256PlusPlus(seed: 7);
        double[] x = [.. Enumerable.Range(0, n).Select(_ => 0.9 * rng.NextDouble())];
        x[0] = 1.0;
        x[1] = 0.9 + 1e-4;
        x[2] = 0.0;
        x[3] = 0.9;   // the knot, a training value
        int[] z = new int[n];
        z[0] = 1;
        var data = new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", z, ["common", "rare"])],
                                new ScalarResponse([.. x.Select(v => v * v)]));
        var basis = new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1);
        int knot = basis.Inputs.ContinuousKnots(0).IndexOf(0.9);
        Assert.True(knot >= 0, "the knot at 0.9 must be a training value");
        BassTerm indicator = new(categorical: [new CategoricalFactor(0, [0])]);
        BassTerm hinge = new([new HingeFactor(0, 1, knot)]);
        var design = new BassDesign(basis, data, Permissive);
        Assert.Null(design.Load([indicator, hinge]));

        double[][] columns = [Enumerable.Repeat(1.0, n).ToArray(), .. design.Terms.Select(t => Case(basis, t))];
        string csv = Csv(columns);
        try
        {
            double[] loo = Loo(ROracle.Run("oracles/bass_admission_oracle.R", csv, "centred"));
            Assert.True(loo.Min() < 1e-4, $"the fixture must be nearly collinear: {string.Join(", ", loo)}");
            Agree(design, loo, "near-constant indicator");
        }
        finally { File.Delete(csv); }
    }

    [Fact]
    public void DirectAngles_MatchTheMaterializedDesignsQr()
    {
        if (!ROracle.IsAvailable) return;
        var rng = new Xoshiro256PlusPlus(seed: 43);
        const int n = 25, m = 17;
        double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] t = [.. Enumerable.Range(0, m).Select(r => 3 + r / (m - 1.0))];   // far from zero, so ρ is small
        var y = new double[n, m];
        for (int i = 0; i < n; i++) for (int r = 0; r < m; r++) y[i, r] = x[i] * t[r];
        var data = new BassData([new ContinuousPredictor("x", x)], new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
        var basis = new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1);
        BassTerm[] terms =
        [
            new([new HingeFactor(0, 1, 3)]), new(functional: [new HingeFactor(0, 1, 2)]),
            new([new HingeFactor(0, -1, 20)], functional: [new HingeFactor(0, -1, 14)]),
            new([new HingeFactor(0, 1, 10)], functional: [new HingeFactor(0, 1, 5)]), new([new HingeFactor(0, 1, 12)], functional: [new HingeFactor(0, 1, 5)]),
        ];
        var design = new BassDesign(basis, data, Permissive);
        Assert.Null(design.Load(terms));

        double[][] a = [Enumerable.Repeat(1.0, n).ToArray(), .. design.Terms.Select(term => Case(basis, term))];
        double[][] h = [Enumerable.Repeat(1.0, m).ToArray(), .. design.Terms.Select(term => Output(basis, term))];
        string aCsv = Csv(a), hCsv = Csv(h);
        try { Agree(design, Loo(ROracle.Run("oracles/bass_admission_oracle.R", aCsv, "hadamard", hCsv)), "direct terms"); }
        finally { File.Delete(aCsv); File.Delete(hCsv); }
    }

    private static double[] Case(BassBasis basis, BassTerm term)
    {
        var column = new double[basis.Inputs.Cases.Cases];
        basis.CaseColumn(term, basis.Inputs.Cases, column);
        return column;
    }

    private static double[] Output(BassBasis basis, BassTerm term)
    {
        var column = new double[basis.Inputs.Outputs!.Locations];
        basis.OutputColumn(term, basis.Inputs.Outputs, column);
        return column;
    }
}
