using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// E8's exact route (<see cref="BassMoments.GradientOuterProduct"/>, slice 1) against coactivity's <c>C_bass</c> on the same draws
/// (<c>active_subspace_oracle.R</c>): an engine fit of degree one on a ridge from outside the model, two inputs on the declared
/// domain [0, 1] so that its transformed units are C_bass's, each draw packed as the mode's rows. Under the uniform law the engine
/// meets <c>C_bass</c> as written; under a truncated normal on x0 it meets <c>C_bass</c> with slice 0's two corrections patched in,
/// and not as written (<see cref="ActiveSubspaceDefectParityTests"/>). Both also meet the mode's own quadrature reference.
/// </summary>
public sealed class ActiveSubspaceParityTests(ITestOutputHelper output)
{
    private static readonly Lazy<(BassData Data, BassFit Fit)> Fixture = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(20261109);
        const int n = 120;
        double[] x0 = new double[n], x1 = new double[n], y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * (0.7 * x0[i] + 0.4 * x1[i])) + 0.05 * Normal.SampleStandard(rng);
        }
        var domain = new InputDomain(0.0, 1.0);
        var data = new BassData([new ContinuousPredictor("x0", x0, domain), new ContinuousPredictor("x1", x1, domain)], new ScalarResponse(y));
        BassConfig config = new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, Degree = 1 },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 7, BurnIn = 1_000, Samples = 6, Thin = 20 },
        };
        return (data, Bass.Fit(config, data));
    });

    [Theory]
    [InlineData("uniform")]
    [InlineData("normal")]
    public void ExactRoute_MeetsCBass_OnTheEnginesDraws(string law)
    {
        if (!ROracle.IsAvailable) return;
        (BassData data, BassFit fit) = Fixture.Value;
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var measures = new Dictionary<string, InputMeasure>();
        if (law == "normal") measures["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.4, 0.25)]);
        var moments = BassMoments.On(inputs, 1, ProductMeasure.Resolve(fit.Transform, inputs, measures));
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];

        var rows = new List<double[]>();
        var exact = new double[draws.Length][,];
        for (int d = 0; d < draws.Length; d++)
        {
            int id = d + 1;
            BassTerm[] terms = [.. draws[d].Structure.TermIds.Select(t => fit.Terms[t])];
            exact[d] = moments.GradientOuterProduct(terms, draws[d].Coefficients.AsSpan(1));
            rows.Add([3, id, 0, 2, 1, 0]);
            for (int m = 0; m < draws[d].Coefficients.Length; m++) rows.Add([1, id, m, draws[d].Coefficients[m], 0, 0]);
            for (int m = 0; m < terms.Length; m++)
                foreach (HingeFactor h in terms[m].Continuous) rows.Add([2, id, m + 1, h.Variable + 1, h.Sign, inputs.ContinuousKnots(h.Variable)[h.Knot]]);
            rows.Add(law == "normal" ? [4, id, 1, 2, 0.4, 0.25] : [4, id, 1, 1, 0.0, 1.0]);
            rows.Add([4, id, 2, 1, 0.0, 1.0]);
        }
        JsonElement? result = GbassOracle.Run("active_subspace_oracle.R", rows);
        if (result is not { } r) return;

        double asWritten = 0.0, corrected = 0.0, reference = 0.0;
        foreach (JsonElement e in r.GetProperty("draws").EnumerateArray())
        {
            int d = e.GetProperty("draw").GetInt32() - 1;
            asWritten = Math.Max(asWritten, Gap(exact[d], Matrix(e.GetProperty("C"))));
            corrected = Math.Max(corrected, Gap(exact[d], Matrix(e.GetProperty("corrected"))));
            reference = Math.Max(reference, Gap(exact[d], Matrix(e.GetProperty("reference"))));
        }
        output.WriteLine($"{law}: {draws.Length} draws, up to {draws.Max(x => x.Structure.TermIds.Length)} terms; the engine against C_bass as written {asWritten:E2}, " +
                         $"corrected {corrected:E2}, the mode's quadrature {reference:E2}");
        Assert.True(corrected <= 1e-12, $"{corrected:E2}");
        Assert.True(reference <= 1e-12, $"{reference:E2}");
        if (law == "uniform") Assert.True(asWritten <= 1e-12, $"{asWritten:E2}");
        else Assert.True(asWritten > 1e-3, $"{asWritten:E2}");
    }

    private static double[,] Matrix(JsonElement e)
    {
        double[][] rows = [.. e.EnumerateArray().Select(GbassOracle.Numbers)];
        var m = new double[rows.Length, rows[0].Length];
        for (int i = 0; i < rows.Length; i++)
            for (int j = 0; j < rows[i].Length; j++) m[i, j] = rows[i][j];
        return m;
    }

    private static double Gap(double[,] engine, double[,] other)
    {
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < engine.GetLength(0); i++)
            for (int j = 0; j < engine.GetLength(1); j++)
            {
                worst = Math.Max(worst, Math.Abs(engine[i, j] - other[i, j]));
                scale = Math.Max(scale, Math.Abs(engine[i, j]));
            }
        return scale == 0.0 ? worst : worst / scale;
    }
}
