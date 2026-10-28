using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's basis against R BASS's own <c>makeBasis</c>, <c>makeBasisCat</c> and <c>scale_range</c> at the pinned
/// commit (<c>r/oracles/bass_basis_oracle.R</c>; stage 4's S2): random terms mixing continuous, categorical and
/// functional factors of every order, at degrees 1, 2 and 3, on predictors with tied values. Parity where R is right
/// (the hub's parity rule). Opt-in: it runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassBasisParityTests
{
    private static readonly string BassSource = Environment.GetEnvironmentVariable("BASS_R_SOURCE") ?? "D:/aghado01/packages/BASS";

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Basis_MatchesRBass(int degree)
    {
        if (!ROracle.IsAvailable || !Directory.Exists(BassSource)) return;
        const int n = 30, m = 9, termCount = 60;
        var rng = new Xoshiro256PlusPlus(seed: 1009 + degree);

        // Ties in every continuous column: values on a coarse grid.
        double[][] x = [.. Enumerable.Range(0, 3).Select(v => Enumerable.Range(0, n).Select(_ => Math.Round(rng.NextDouble() * 12) * (v + 0.7) - 3).ToArray())];
        string[][] levelNames = [["a", "b", "c"], ["p", "q", "r", "s", "u"]];
        int[][] codes = [.. levelNames.Select(l => Enumerable.Range(0, n).Select(i => i < l.Length ? i : (int)(rng.NextDouble() * l.Length)).ToArray())];
        double[][] t = [Enumerable.Range(0, m).Select(r => r * 0.125).ToArray(), Enumerable.Range(0, m).Select(r => Math.Sin(r)).ToArray()];
        var y = new double[n, m];
        for (int i = 0; i < n; i++) for (int r = 0; r < m; r++) y[i, r] = rng.NextDouble();

        var predictors = new List<BassPredictor>();
        predictors.AddRange(x.Select((c, v) => (BassPredictor)new ContinuousPredictor($"x{v}", c)));
        predictors.AddRange(codes.Select((c, v) => (BassPredictor)new CategoricalPredictor($"z{v}", c, levelNames[v])));
        var data = new BassData(predictors, new FunctionalResponse(y, [new ContinuousPredictor("t0", t[0]), new ContinuousPredictor("t1", t[1])]));
        var basis = new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), degree);
        BassInputs inputs = basis.Inputs;

        var terms = new List<BassTerm>();
        var fixtureTerms = new List<object>();
        while (terms.Count < termCount)
        {
            int[] cont = Pick(rng, 3, (int)(rng.NextDouble() * 4));
            int[] cat = Pick(rng, 2, (int)(rng.NextDouble() * 3));
            int[] func = Pick(rng, 2, (int)(rng.NextDouble() * 3));
            if (cont.Length + cat.Length + func.Length == 0) continue;
            var contFactors = cont.Select(v => (v, s: rng.NextDouble() < 0.5 ? -1 : 1, at: (int)(rng.NextDouble() * n))).ToArray();
            var funcFactors = func.Select(v => (v, s: rng.NextDouble() < 0.5 ? -1 : 1, at: (int)(rng.NextDouble() * m))).ToArray();
            var catFactors = cat.Select(v =>
            {
                int levels = levelNames[v].Length, size = 1 + (int)(rng.NextDouble() * (levels - 1));
                return (v, sub: Pick(rng, levels, size));
            }).ToArray();
            terms.Add(new BassTerm(
                contFactors.Select(f => new HingeFactor(f.v, f.s, KnotOf(inputs.ContinuousKnots(f.v), inputs.Cases.Continuous(f.v)[f.at]))),
                catFactors.Select(f => new CategoricalFactor(f.v, f.sub)),
                funcFactors.Select(f => new HingeFactor(f.v, f.s, KnotOf(inputs.FunctionalKnots(f.v), inputs.Outputs!.Variable(f.v)[f.at])))));
            fixtureTerms.Add(new
            {
                cont = contFactors.Select(f => new { var = f.v, sign = f.s, @case = f.at }),
                cat = catFactors.Select(f => new { var = f.v, levels = f.sub.Select(l => levelNames[f.v][l]) }),
                func = funcFactors.Select(f => new { var = f.v, sign = f.s, loc = f.at }),
            });
        }

        string fixture = Path.Combine(Path.GetTempPath(), $"bass_basis_{Guid.NewGuid():N}.json");
        File.WriteAllText(fixture, JsonSerializer.Serialize(new
        {
            degree,
            x,
            cat = codes.Select((c, v) => new { levels = levelNames[v], values = c.Select(code => levelNames[v][code]) }),
            t,
            terms = fixtureTerms,
        }));
        try
        {
            JsonElement r = ROracle.Run("oracles/bass_basis_oracle.R", fixture);
            JsonElement[] scaled = [.. r.GetProperty("xxt_des").EnumerateArray()];
            for (int v = 0; v < x.Length; v++)
                Assert.Equal(Doubles(scaled[v]), inputs.Cases.Continuous(v).ToArray(), (e, a) => Math.Abs(e - a) <= 1e-15);
            JsonElement[] caseColumns = [.. r.GetProperty("case").EnumerateArray()];
            JsonElement[] outputColumns = [.. r.GetProperty("output").EnumerateArray()];
            var a = new double[n];
            var h = new double[m];
            for (int j = 0; j < terms.Count; j++)
            {
                basis.Validate(terms[j]);
                basis.CaseColumn(terms[j], inputs.Cases, a);
                basis.OutputColumn(terms[j], inputs.Outputs!, h);
                double[] ra = Doubles(caseColumns[j]), rh = Doubles(outputColumns[j]);
                for (int i = 0; i < n; i++) Assert.True(Math.Abs(a[i] - ra[i]) <= 1e-14, $"term {terms[j]}, case {i}: {a[i]:R} against R's {ra[i]:R}");
                for (int l = 0; l < m; l++) Assert.True(Math.Abs(h[l] - rh[l]) <= 1e-14, $"term {terms[j]}, location {l}: {h[l]:R} against R's {rh[l]:R}");
            }
        }
        finally
        {
            File.Delete(fixture);
        }
    }

    // count distinct values from 0 … range − 1, in random order.
    private static int[] Pick(Xoshiro256PlusPlus rng, int range, int count)
        => [.. Enumerable.Range(0, range).OrderBy(_ => rng.NextDouble()).Take(count)];

    private static int KnotOf(ReadOnlySpan<double> knots, double value)
    {
        for (int k = 0; k < knots.Length; k++)
            if (knots[k] == value) return k;
        throw new InvalidOperationException($"{value:R} is not a knot.");
    }

    private static double[] Doubles(JsonElement e) => [.. e.EnumerateArray().Select(v => v.GetDouble())];
}
