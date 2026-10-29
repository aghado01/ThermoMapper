using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Effect functions of basis fits (O8 for basis fits; <see cref="BassBasisEffects"/>): a one-component fit reproduces its
/// component's direct effects scaled by y_s(ℓ)·K_ℓ1, with the mean y_m(ℓ) + y_s(ℓ)·K_ℓ1·E[w]; on a two-component fit each
/// location's main and two-way effects, and the integrated form's, agree with tensor quadrature of the reconstructed latent
/// output; the default set is the union of the components' pairs; and the refusals.
/// </summary>
public sealed class BassSensitivityBasisEffectTests
{
    private readonly ITestOutputHelper _out;
    public BassSensitivityBasisEffectTests(ITestOutputHelper output) => _out = output;

    private static readonly Lazy<BassBasisFit> OneComponent = new(() =>
    {
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 1, ScaleByColumnSd = true });
        return BassComponents.Fit(d, predictors, BassOutputComponentsTests.Config());
    });

    [Fact]
    public void OneComponent_ReproducesItsDirectEffects_ScaledByTheBasis()
    {
        BassBasisFit fit = OneComponent.Value;
        var effects = new EffectSpec { GridPoints = 7 };
        BasisEffectAnalysis analysis = fit.Effects(effects: effects);
        EffectFunctions direct = fit.Components[0].Effects(fit.ComponentData(0), effects: effects).At[0];
        LinearReconstruction map = fit.Decomposition.Map;
        Assert.Equal(1, analysis.Components);
        Assert.Equal(direct.Mean.Length, analysis.Count);   // one component: the joint draws are its own, in order
        foreach (int l in new[] { 0, 6, map.Locations - 1 })
        {
            double g = map.Scales[l] * map.Basis(l, 0);
            EffectFunctions at = analysis.AtLocation(l);
            Assert.Equal(direct.Effects.Count, at.Effects.Count);
            for (int e = 0; e < at.Effects.Count; e++)
            {
                Assert.Equal(direct.Effects[e].Inputs, at.Effects[e].Inputs);
                double[,] a = at.Effects[e].Values, b = direct.Effects[e].Values;
                for (int d = 0; d < analysis.Count; d++)
                    for (int c = 0; c < a.GetLength(1); c++) Assert.Equal(g * b[d, c], a[d, c], 13);
            }
            for (int d = 0; d < analysis.Count; d++) Assert.Equal(map.Offsets[l] + g * direct.Mean[d], at.Mean[d], 12);
        }
        EffectFunctions own = analysis.ComponentEffects(0);
        Assert.Equal(direct.Mean, own.Mean);
        Assert.Equal(direct.Effects[0].Values, own.Effects[0].Values);
    }

    // The effects by tensor quadrature of the reconstructed latent output, under the uniform measure on [0, 1]²: at grid
    // points `grid`, by the midpoint rule on `nodes` points per input, which is exact on each linear piece of a degree-1
    // surface and errs O(h²) only on the cells a knot cuts. Returns location × (main z0 at the grid, main z1, two-way),
    // each draw × cells, and the means, draw × location.
    private static (double[][,] Main0, double[][,] Main1, double[][,] Pair, double[,] Mean) Quadrature(BassBasisFit fit, int[] draws, double[] grid, int nodes)
    {
        double[] mid = [.. Enumerable.Range(0, nodes).Select(j => (j + 0.5) / nodes)];
        double[,,] Latent(IEnumerable<(double Z0, double Z1)> points)
        {
            var list = points.ToArray();
            var cases = new PredictionCases([new ContinuousColumn("z0", [.. list.Select(p => p.Z0)]), new ContinuousColumn("z1", [.. list.Select(p => p.Z1)])]);
            return fit.Predict(cases, new BassPredictOptions { Draws = draws, Predictive = false }).Latent;
        }
        int s = draws.Length, m = fit.Decomposition.Locations, gp = grid.Length;
        // f on the node grid, the grid × nodes lines both ways, and the grid × grid cells.
        double[,,] full = Latent(from a in mid from b in mid select (a, b));
        double[,,] lines0 = Latent(from g in grid from b in mid select (g, b));
        double[,,] lines1 = Latent(from g in grid from a in mid select (a, g));
        double[,,] cells = Latent(from g0 in grid from g1 in grid select (g0, g1));
        var mean = new double[s, m];
        var main0 = new double[m][,];
        var main1 = new double[m][,];
        var pair = new double[m][,];
        for (int l = 0; l < m; l++)
        {
            main0[l] = new double[s, gp];
            main1[l] = new double[s, gp];
            pair[l] = new double[s, gp * gp];
            for (int d = 0; d < s; d++)
            {
                double f0 = 0.0;
                for (int c = 0; c < nodes * nodes; c++) f0 += full[d, c, l];
                f0 /= nodes * nodes;
                mean[d, l] = f0;
                for (int i = 0; i < gp; i++)
                {
                    double a = 0.0, b = 0.0;
                    for (int j = 0; j < nodes; j++)
                    {
                        a += lines0[d, i * nodes + j, l];
                        b += lines1[d, i * nodes + j, l];
                    }
                    main0[l][d, i] = a / nodes - f0;
                    main1[l][d, i] = b / nodes - f0;
                }
                for (int i = 0; i < gp; i++)
                    for (int k = 0; k < gp; k++)
                        pair[l][d, i * gp + k] = cells[d, i * gp + k, l] - main0[l][d, i] - main1[l][d, k] - f0;
            }
        }
        return (main0, main1, pair, mean);
    }

    [Fact]
    public void EachLocationAndTheIntegratedForm_AgreeWithTensorQuadrature()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        Assert.True(fit.Components.Count >= 2);
        const int gp = 5, nodes = 160;
        BasisEffectAnalysis analysis = fit.Effects(effects: new EffectSpec { Subsets = [["z0"], ["z1"], ["z0", "z1"]], GridPoints = gp });
        int[] draws = [0, analysis.Count / 2, analysis.Count - 1];
        double[] grid = [.. Enumerable.Range(0, gp).Select(i => i / (gp - 1.0))];
        var (main0, main1, pair, mean) = Quadrature(fit, draws, grid, nodes);
        int m = analysis.Locations;
        double scale = 0.0, worst = 0.0;
        void Near(double expected, double actual)
        {
            scale = Math.Max(scale, Math.Abs(expected));
            worst = Math.Max(worst, Math.Abs(expected - actual));
        }
        for (int l = 0; l < m; l++)
        {
            EffectFunctions at = analysis.AtLocation(l);
            Assert.Equal(grid, at.Effects[0].Grids[0].Values);
            for (int j = 0; j < draws.Length; j++)
            {
                int d = draws[j];
                Near(mean[j, l], at.Mean[d]);
                for (int c = 0; c < gp; c++)
                {
                    Near(main0[l][j, c], at.Effects[0].Values[d, c]);
                    Near(main1[l][j, c], at.Effects[1].Values[d, c]);
                }
                for (int c = 0; c < gp * gp; c++) Near(pair[l][j, c], at.Effects[2].Values[d, c]);
            }
        }

        // The integrated form: the location-averaged effects, the location's main effect, and each input with the location.
        EffectFunctions integrated = analysis.Integrated();
        Assert.Equal(new[] { "z0", "z1", BasisEffectAnalysis.LocationInput }, integrated.Inputs);
        Assert.Equal(new[] { new[] { 0 }, new[] { 1 }, new[] { 0, 1 }, new[] { 2 }, new[] { 0, 2 }, new[] { 1, 2 } }, integrated.Effects.Select(e => e.Inputs));
        for (int j = 0; j < draws.Length; j++)
        {
            int d = draws[j];
            double f0 = Enumerable.Range(0, m).Average(l => mean[j, l]);
            Near(f0, integrated.Mean[d]);
            for (int l = 0; l < m; l++) Near(mean[j, l] - f0, integrated.Effects[3].Values[d, l]);
            for (int c = 0; c < gp; c++)
            {
                double a0 = Enumerable.Range(0, m).Average(l => main0[l][j, c]), a1 = Enumerable.Range(0, m).Average(l => main1[l][j, c]);
                Near(a0, integrated.Effects[0].Values[d, c]);
                Near(a1, integrated.Effects[1].Values[d, c]);
                for (int l = 0; l < m; l++)
                {
                    Near(main0[l][j, c] - a0, integrated.Effects[4].Values[d, c * m + l]);
                    Near(main1[l][j, c] - a1, integrated.Effects[5].Values[d, c * m + l]);
                }
            }
            for (int c = 0; c < gp * gp; c++) Near(Enumerable.Range(0, m).Average(l => pair[l][j, c]), integrated.Effects[2].Values[d, c]);
        }
        _out.WriteLine($"largest gap against quadrature {worst:E2} on effects up to {scale:F3}");
        Assert.InRange(worst, 0.0, 1e-3 * scale);
    }

    [Fact]
    public void Effects_FollowTheJointDraws_WhenWeightedComponentsAreResampled()
    {
        // Weighted components are resampled before pairing, so joint draw d takes a shuffled pooled draw of each: every
        // location's effect must read component k's effect at JointDraws.Draw(k, d).
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        (_, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassFit[] weighted = [.. fit.Components.Select(c => c with { DrawWeights = [.. Enumerable.Range(0, c.Draws.Sum(x => x.Count)).Select(d => 0.2 + (d * 37 % 11) / 4.0)] })];
        BassBasisFit resampled = BassComponents.Assemble(fit.Decomposition, predictors, weighted, 3);
        Assert.True(resampled.JointDraws.Resampled);
        var effects = new EffectSpec { Subsets = [["z0"], ["z0", "z1"]], GridPoints = 4 };
        BasisEffectAnalysis analysis = resampled.Effects(effects: effects);
        EffectFunctions[] own = [.. Enumerable.Range(0, weighted.Length).Select(k => weighted[k].Effects(resampled.ComponentData(k), effects: effects).At[0])];
        LinearReconstruction map = fit.Decomposition.Map;
        int l = 5;
        EffectFunctions at = analysis.AtLocation(l);
        for (int d = 0; d < analysis.Count; d++)
        {
            double mean = map.Offsets[l];
            for (int k = 0; k < own.Length; k++) mean += map.Scales[l] * map.Basis(l, k) * own[k].Mean[resampled.JointDraws.Draw(k, d)];
            Assert.Equal(mean, at.Mean[d], 12);
            for (int e = 0; e < 2; e++)
                for (int c = 0; c < at.Effects[e].Values.GetLength(1); c++)
                {
                    double expected = 0.0;
                    for (int k = 0; k < own.Length; k++)
                        expected += map.Scales[l] * map.Basis(l, k) * own[k].Effects[e].Values[resampled.JointDraws.Draw(k, d), c];
                    Assert.Equal(expected, at.Effects[e].Values[d, c], 12);
                }
        }
    }

    [Fact]
    public void DefaultSet_IsEveryMainAndTheComponentsPairs()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BasisEffectAnalysis analysis = fit.Effects(effects: new EffectSpec { GridPoints = 3 });
        var expected = new SortedSet<(int, int)>();
        for (int k = 0; k < fit.Components.Count; k++)
            foreach (Effect e in fit.Components[k].Effects(fit.ComponentData(k), effects: new EffectSpec { GridPoints = 3 }).At[0].Effects)
                if (e.Inputs.Length == 2) expected.Add((e.Inputs[0], e.Inputs[1]));
        EffectFunctions at = analysis.AtLocation(0);
        Assert.Equal(new[] { new[] { 0 }, new[] { 1 } }, at.Effects.Take(2).Select(e => e.Inputs));
        Assert.Equal(expected.Select(p => new[] { p.Item1, p.Item2 }), at.Effects.Skip(2).Select(e => e.Inputs));
        Assert.Equal(analysis.Locations, analysis.PerLocation().Count);
        _out.WriteLine($"pairs: {string.Join(", ", expected)}");
    }

    [Fact]
    public void Refusals_CarryTheirReasons()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        Assert.Contains("training locations", Assert.Throws<ArgumentException>(() => fit.Effects(new SensitivitySpec { Coordinates = [[0.5]] })).Message);
        Assert.Contains("output-variable measures", Assert.Throws<ArgumentException>(() =>
            fit.Effects(new SensitivitySpec { OutputMeasures = new Dictionary<string, InputMeasure> { ["t"] = new UniformMeasure() } })).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.Effects(effects: new EffectSpec { GridPoints = 1 }));
        Assert.Contains("No input is named", Assert.Throws<ArgumentException>(() => fit.Effects(effects: new EffectSpec { Subsets = [["w"]] })).Message);

        // The budget counts the components' state; the integrated form, built on demand, checks its own cells.
        int r = fit.Components.Count, s = fit.JointDraws.Count;
        var two = new EffectSpec { Subsets = [["z0"], ["z1"]], GridPoints = 5 };
        Assert.Contains("budget", Assert.Throws<ArgumentException>(() => fit.Effects(effects: two with { CellBudget = 10L * s * r - 1 })).Message);
        BasisEffectAnalysis tight = fit.Effects(effects: two with { CellBudget = 10L * s * r });
        Assert.Equal(fit.Decomposition.Locations, tight.PerLocation().Count);
        Assert.Contains("budget", Assert.Throws<ArgumentException>(() => tight.Integrated()).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => tight.AtLocation(fit.Decomposition.Locations));
        Assert.Throws<ArgumentOutOfRangeException>(() => tight.ComponentEffects(r));
    }
}
