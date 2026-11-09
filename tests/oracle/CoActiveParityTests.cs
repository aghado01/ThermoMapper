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
/// E8's slice 3 against coactivity (knrumsey/coactivity@35e4573) through <c>active_subspace_oracle.R</c>:
/// <list type="bullet">
/// <item>two engine fits' co-active pairs against <c>Cfg_bass</c> on the same draw pairs: as written under the uniform law, with slice
/// 0's corrections patched into <c>Cfg_bass_univariate</c> under a truncated normal, and the mode's quadrature under both;</item>
/// <item><c>Cfg_bass</c>'s readings: one pair by default, a list in mod2 alone refused;</item>
/// <item>a direct functional fit's C at each training output location (W = e_ℓe_ℓᵀ) against coactivity's per-location reading,
/// <c>bassfunc2bass_fixed_t</c> then <c>C_bass_univariate</c> (corrected), and the mode's quadrature of the draw folded at t;</item>
/// <item>the functional converters' readings, counted in the pinned bodies or run: <c>func.true</c> unread, the equispaced grid, the
/// default locations between the training ones, the last functional variable only, the version-one stubs, and two stops.</item>
/// </list>
/// Every comparison here read <c>coactivity</c>; <c>concordance</c> was read and not run.
/// </summary>
public sealed class CoActiveParityTests(ITestOutputHelper output)
{
    // Two responses on x0, x1 drawn on the declared domain [0, 1], so the engine's transformed units are coactivity's; degree one.
    private static readonly Lazy<(BassData F, BassFit FitF, BassData G, BassFit FitG)> Pair = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(20261110);
        const int n = 120;
        double[] x0 = new double[n], x1 = new double[n], yf = new double[n], yg = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            yf[i] = Math.Sin(2.0 * Math.PI * (0.7 * x0[i] + 0.4 * x1[i])) + 0.05 * Normal.SampleStandard(rng);
            yg[i] = Math.Cos(2.0 * Math.PI * (0.3 * x0[i] - 0.6 * x1[i])) + x0[i] + 0.05 * Normal.SampleStandard(rng);
        }
        var domain = new InputDomain(0.0, 1.0);
        BassPredictor[] predictors = [new ContinuousPredictor("x0", x0, domain), new ContinuousPredictor("x1", x1, domain)];
        BassData f = new(predictors, new ScalarResponse(yf)), g = new(predictors, new ScalarResponse(yg));
        return (f, Bass.Fit(Config(), f), g, Bass.Fit(Config(), g));
    });

    // A functional response over t at nine equispaced locations on [0, 1], so its transformed units are the output's own.
    private static readonly Lazy<(BassData Data, BassFit Fit)> Functional = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(20261111);
        const int n = 80, m = 9;
        double[] t = [.. Enumerable.Range(0, m).Select(l => l / (m - 1.0))];
        double[] x0 = new double[n], x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            for (int l = 0; l < m; l++)
                y[i, l] = Math.Sin(2.0 * Math.PI * (0.7 * x0[i] + 0.4 * x1[i])) * (1.0 + t[l]) + 0.8 * t[l] * x1[i] + 0.05 * Normal.SampleStandard(rng);
        }
        var domain = new InputDomain(0.0, 1.0);
        var data = new BassData([new ContinuousPredictor("x0", x0, domain), new ContinuousPredictor("x1", x1, domain)],
                                new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
        return (data, Bass.Fit(Config(), data));
    });

    private static BassConfig Config() => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, Degree = 1 },
        Sampling = new BassSampling { Chains = 1, MasterSeed = 7, BurnIn = 1_000, Samples = 6, Thin = 20 },
    };

    // A draw's rows: its header (last, what it is for: 1 a functional draw, 2 a draw read only within a pair), coefficients, design
    // hinges, functional hinges and laws.
    private static void AddDraw(List<double[]> rows, int id, BassFit fit, BassInputs inputs, BassDraw draw, string law, int kind)
    {
        BassTerm[] terms = [.. draw.Structure.TermIds.Select(t => fit.Terms[t])];
        rows.Add([3, id, 0, 2, 1, kind]);
        for (int m = 0; m < draw.Coefficients.Length; m++) rows.Add([1, id, m, draw.Coefficients[m], 0, 0]);
        for (int m = 0; m < terms.Length; m++)
        {
            foreach (HingeFactor h in terms[m].Continuous) rows.Add([2, id, m + 1, h.Variable + 1, h.Sign, inputs.ContinuousKnots(h.Variable)[h.Knot]]);
            foreach (HingeFactor h in terms[m].Functional) rows.Add([6, id, m + 1, 0, h.Sign, inputs.FunctionalKnots(h.Variable)[h.Knot]]);
        }
        rows.Add(law == "normal" ? [4, id, 1, 2, 0.4, 0.25] : [4, id, 1, 1, 0.0, 1.0]);
        rows.Add([4, id, 2, 1, 0.0, 1.0]);
    }

    private static Dictionary<string, InputMeasure> Measures(string law) =>
        law == "normal" ? new() { ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.4, 0.25)]) } : new();

    /// <summary>
    /// Each pair's C_fg on the exact route against <c>Cfg_bass</c> on the same two draws: as written under the uniform law, and with
    /// slice 0's two corrections patched into <c>Cfg_bass_univariate</c> under a truncated normal on x0, where the code as written
    /// misses; and the mode's quadrature of Σ w∇f∇gᵀ under both. Under the uniform law also: <c>Cfg_bass</c> without mcmc.use reads
    /// the last pair only, and a list in mod2 alone is refused.
    /// </summary>
    [Theory]
    [InlineData("uniform")]
    [InlineData("normal")]
    public void ExactRoute_MeetsCfgBass_OnTwoFitsDraws(string law)
    {
        if (!ROracle.IsAvailable) return;
        (BassData df, BassFit f, BassData dg, BassFit g) = Pair.Value;
        BassInputs inputs = BassInputs.Prepare(df, f.Transform);
        BassDraw[] fd = [.. f.Draws.SelectMany(c => c)], gd = [.. g.Draws.SelectMany(c => c)];
        int count = Math.Min(fd.Length, gd.Length);
        int[] selected = [.. Enumerable.Range(0, count)];
        CoActiveAnalysis a = f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Measures = Measures(law), Draws = selected, OtherDraws = selected });
        var rows = new List<double[]>();
        for (int s = 0; s < count; s++)
        {
            AddDraw(rows, s + 1, f, inputs, fd[s], law, 2);
            AddDraw(rows, count + s + 1, g, inputs, gd[s], law, 2);
            rows.Add([5, s + 1, s + 1, count + s + 1, 0, 0]);
        }
        JsonElement? result = GbassOracle.Run("active_subspace_oracle.R", rows);
        if (result is not { } r) return;

        double asWritten = 0.0, corrected = 0.0, reference = 0.0;
        foreach (JsonElement e in r.GetProperty("pairs").EnumerateArray())
        {
            int s = e.GetProperty("pair").GetInt32() - 1;
            double[,] cfg = a.Posterior.Pairs[s].Cfg;
            asWritten = Math.Max(asWritten, Gap(cfg, Matrix(e.GetProperty("Cfg"))));
            corrected = Math.Max(corrected, Gap(cfg, Matrix(e.GetProperty("corrected"))));
            reference = Math.Max(reference, Gap(cfg, Matrix(e.GetProperty("reference"))));
        }
        output.WriteLine($"{law}: {count} pairs; the engine against Cfg_bass as written {asWritten:E2}, corrected {corrected:E2}, the mode's quadrature {reference:E2}");
        Assert.True(corrected <= 1e-12, $"{corrected:E2}");
        Assert.True(reference <= 1e-12, $"{reference:E2}");
        if (law == "uniform") Assert.True(asWritten <= 1e-12, $"{asWritten:E2}");
        else Assert.True(asWritten > 1e-3, $"{asWritten:E2}");

        JsonElement defaults = r.GetProperty("cfg_default_draws");
        double[,] byDefault = Matrix(defaults.GetProperty("default")), first = Matrix(defaults.GetProperty("first")), last = Matrix(defaults.GetProperty("last"));
        Assert.Equal(0.0, Gap(byDefault, last));
        Assert.True(Gap(byDefault, first) > 0.1);
        JsonElement list = r.GetProperty("cfg_list_mod2");
        Assert.False(list.GetProperty("ran").GetBoolean());
        Assert.Contains("mod2 must be", list.GetProperty("error").GetString());
    }

    /// <summary>
    /// A direct functional fit's C at each training output location, the readout under W = e_ℓe_ℓᵀ, against coactivity's per-location
    /// reading (<c>bassfunc2bass_fixed_t</c> at t_ℓ, then <c>C_bass_univariate</c> with slice 0's corrections) and the mode's quadrature
    /// of the draw with its functional hinges folded at t_ℓ, under the uniform law: four draws of two or more terms (the converter's
    /// one-basis-function case is read, not run) at four of the nine locations, the ends among them, within the oracle's two minutes.
    /// </summary>
    [Fact]
    public void PerLocationC_MeetsBassfunc2bass()
    {
        if (!ROracle.IsAvailable) return;
        (BassData data, BassFit fit) = Functional.Value;
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        int[] selected = [.. Enumerable.Range(0, pooled.Length).Where(d => pooled[d].Structure.TermIds.Length >= 2).Take(4)];
        Assert.NotEmpty(selected);
        int m = data.Response.Locations;
        int[] read = [0, 3, 5, m - 1];
        var perLocation = new ActiveSubspaceAnalysis[m];
        for (int l = 0; l < m; l++)
        {
            var e = new double[m];
            e[l] = 1.0;
            perLocation[l] = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Draws = selected, OutputMetric = OutputMetric.Diagonal(e) });
        }
        var rows = new List<double[]>();
        for (int s = 0; s < selected.Length; s++) AddDraw(rows, s + 1, fit, inputs, pooled[selected[s]], "uniform", 1);
        ReadOnlySpan<double> t = ((FunctionalResponse)data.Response).OutputVariables[0].Values.Span;
        foreach (int l in read) rows.Add([7, 0, 0, l + 1, fit.Transform.OutputVariables[0].Scale(t[l]), 0]);
        JsonElement? result = GbassOracle.Run("active_subspace_oracle.R", rows);
        if (result is not { } r) return;

        double corrected = 0.0, reference = 0.0;
        int compared = 0;
        foreach (JsonElement d in r.GetProperty("functional").EnumerateArray())
        {
            int s = d.GetProperty("draw").GetInt32() - 1;
            foreach (JsonElement at in d.GetProperty("locations").EnumerateArray())
            {
                int l = read[at.GetProperty("location").GetInt32() - 1];
                double[,] c = perLocation[l].Posterior.Draws[s].C;
                corrected = Math.Max(corrected, Gap(c, Matrix(at.GetProperty("corrected"))));
                reference = Math.Max(reference, Gap(c, Matrix(at.GetProperty("reference"))));
                compared++;
            }
        }
        output.WriteLine($"{selected.Length} draws × {read.Length} locations ({compared} compared); the engine against bassfunc2bass_fixed_t + C_bass {corrected:E2}, " +
                         $"the mode's quadrature {reference:E2}");
        Assert.Equal(selected.Length * read.Length, compared);
        Assert.True(corrected <= 1e-12, $"{corrected:E2}");
        Assert.True(reference <= 1e-12, $"{reference:E2}");

        // The converters' readings.
        JsonElement readings = r.GetProperty("functional_readings");
        Assert.False(readings.GetProperty("func_true_in_body").GetBoolean());
        Assert.Equal(3, readings.GetProperty("equispaced_grid").GetInt32());
        Assert.Equal(1, readings.GetProperty("get_u_grid").GetInt32());
        Assert.Equal(1, readings.GetProperty("pca_v2_default_func_use").GetInt32());
        Assert.Equal(1, readings.GetProperty("last_functional_variable").GetInt32());
        Assert.True(readings.GetProperty("pca_v1").GetBoolean());
        Assert.True(readings.GetProperty("cfg_pca_v1").GetBoolean());
        // C_bass_pca_v2's default for 15 basis rows: 15 points at multiples of 1/16, of which only 1/2 is on the grid of multiples of 1/14.
        double[] points = GbassOracle.Numbers(readings.GetProperty("default_points_15"));
        Assert.Equal(15, points.Length);
        Assert.Equal(new[] { 0.5 }, points.Where(x => Math.Abs(x * 14 - Math.Round(x * 14)) < 1e-12));
        JsonElement one = r.GetProperty("functional_one_draw"), noDrop = r.GetProperty("functional_no_drop");
        Assert.False(one.GetProperty("ran").GetBoolean());
        Assert.Contains("incorrect number of dimensions", one.GetProperty("error").GetString());
        Assert.False(noDrop.GetProperty("ran").GetBoolean());
        Assert.Contains("subscript out of bounds", noDrop.GetProperty("error").GetString());
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
