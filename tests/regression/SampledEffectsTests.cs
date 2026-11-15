using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled-sensitivity brief's slice 2 gates: effects on the sampled route for a scalar fit (<see cref="BassSensitivity.Effects"/>
/// with <see cref="SensitivityRoute"/>, decision 5), each cell's ANOVA component by pinning its inputs over the 2N shared base rows.
/// <list type="bullet">
/// <item>Separable fits read on both routes: every draw's sampled cells and mean against the exact ones, replicated over seeds, within
/// their stated errors, under every kind of product law, categorical inputs and two-way effects included.</item>
/// <item>A hand-built directional fit against tensor quadrature of its surface; a chain-built ridge-block fit (E36's slice 3) against
/// quadrature of its own predicted surface: methods apart from the estimator.</item>
/// <item>The grid's mean against the surface's mean; the route, the shared points, the budgets and the refusals.</item>
/// </list>
/// </summary>
public sealed class SampledEffectsTests(ITestOutputHelper output)
{
    private static readonly Dictionary<string, InputMeasure> MixedLaws = new()
    {
        ["x0"] = new MixtureMeasure([new NormalComponent(0.5, 0.35, 0.15), new NormalComponent(0.5, 0.75, 0.25)]),
        ["x1"] = new UniformMeasure(0.1, 0.9),
        ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }),
    };

    private static readonly Dictionary<string, InputMeasure> EmpiricalAndStudent = new()
    {
        ["x0"] = new EmpiricalMeasure(),
        ["x1"] = new MixtureMeasure([new StudentComponent(1.0, 4.0, 0.5, 0.2)]),
    };

    private static readonly Lazy<BassFit> Legendre = new(() => Bass.Fit(KhaosReadoutTests.Config(5), KhaosReadoutTests.Data(120, 20261110)));

    private static (BassFit Fit, BassData Data, Dictionary<string, InputMeasure>? Laws) Separable(string name) => name switch
    {
        "hinges, mixed laws" => (SensitivityFixture.FitScalar(), SensitivityFixture.Scalar(), MixedLaws),
        "hinges, the boxes" => (SensitivityFixture.FitScalar(), SensitivityFixture.Scalar(), null),
        "degree 2, empirical and Student" => (SensitivityFixture.FitScalar(2), SensitivityFixture.Scalar(), EmpiricalAndStudent),
        "Legendre" => (Legendre.Value, KhaosReadoutTests.Data(120, 20261110), null),
        _ => throw new ArgumentException(name),
    };

    public static TheoryData<string> SeparableNames() => ["hinges, mixed laws", "hinges, the boxes", "degree 2, empirical and Student", "Legendre"];

    // Every draw's mean and every effect's cells, against the truth's.
    private static Dictionary<string, SampledSensitivityTests.Piece> Cells(EffectFunctions sampled, Func<int, double> mean, Func<int, int, int, double> truth)
    {
        var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
        for (int d = 0; d < sampled.Mean.Length; d++)
        {
            cells[$"draw {d} mean"] = new("mean", mean(d), sampled.Mean[d], sampled.MeanError![d]);
            for (int k = 0; k < sampled.Effects.Count; k++)
            {
                Effect effect = sampled.Effects[k];
                string label = string.Join(",", effect.Inputs);
                for (int c = 0; c < effect.Values.GetLength(1); c++)
                    cells[$"draw {d} f_{label}[{c}]"] = new($"f_{label}[{c}]", truth(d, k, c), effect.Values[d, c], effect.Errors![d, c]);
            }
        }
        return cells;
    }

    private void Report(string name, SampledSensitivityTests.Tally t) =>
        output.WriteLine($"{name}: {t.Cells} cells, {t.Exact} exact zeros; worst √R·|mean z| {t.WorstBias:F2} ({t.WorstBiasCell}); pooled sd of z "
                         + $"{t.PooledSd:F3}, per cell {t.LowestSd:F2} to {t.HighestSd:F2}; beyond 3: {t.Beyond3:P2}; worst z {t.WorstZ:F2} ({t.WorstZCell})");

    [Theory]
    [MemberData(nameof(SeparableNames))]
    public void SeparableFits_SampledMeetsExact_WithinTheStatedError(string name)
    {
        const int Seeds = 24, N = 1024;
        (BassFit fit, BassData data, Dictionary<string, InputMeasure>? laws) = Separable(name);
        var grid = new EffectSpec { GridPoints = 7 };
        EffectFunctions exact = BassSensitivity.Effects(fit, data, new SensitivitySpec { Measures = laws, Route = SensitivityRoute.Exact }, grid).At[0];
        double scale = exact.Effects.Max(e => Enumerable.Range(0, e.Values.GetLength(0)).Max(d => Enumerable.Range(0, e.Values.GetLength(1)).Max(c => Math.Abs(e.Values[d, c]))));
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(20261116, Seeds))
        {
            EffectAnalysis analysis = BassSensitivity.Effects(fit, data, new SensitivitySpec { Measures = laws, Route = SensitivityRoute.Sampled, Points = N, Seed = seed }, grid);
            Assert.Equal(SensitivityRoute.Sampled, analysis.Route);
            Assert.Equal(N, analysis.Points);
            EffectFunctions sampled = analysis.At[0];
            Assert.Equal(exact.Effects.Select(e => string.Join(",", e.Inputs)), sampled.Effects.Select(e => string.Join(",", e.Inputs)));
            for (int k = 0; k < exact.Effects.Count; k++)
                for (int v = 0; v < exact.Effects[k].Grids.Count; v++) Assert.Equal(exact.Effects[k].Grids[v].Values, sampled.Effects[k].Grids[v].Values);
            runs.Add(Cells(sampled, d => exact.Mean[d], (d, k, c) => exact.Effects[k].Values[d, c]));
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, cell => cell.EndsWith(" mean") ? Math.Abs(exact.Mean[int.Parse(cell.Split(' ')[1])]) + scale : scale);
        Report(name, t);
        Assert.True(t.WorstBias <= SampledSensitivityTests.BiasCut(t.Cells), $"{t.WorstBiasCell}: √R·|mean z| {t.WorstBias:F2} against {SampledSensitivityTests.BiasCut(t.Cells):F2}");
        Assert.True(t.PooledSd >= SampledSensitivityTests.PooledLow && t.PooledSd <= SampledSensitivityTests.PooledHigh, $"pooled sd of z {t.PooledSd:F3}");
        Assert.True(t.Beyond3 <= 0.01, $"{t.Beyond3:P2} beyond 3");
    }

    // ── A directional fit by hand, and a ridge-block fit from the chain ────────────────────────────────────────────────────

    // The hand fit's draws and boxes, its surfaces (SampledSensitivityTests' hand computation) and their quadrature.
    private static (BassFit Fit, BassData Data, QuadratureReference[] Truth) HandDirectional(int degree, int panels)
    {
        BassData data = SampledSensitivityTests.DirectionalData();
        BassFit fit = SampledSensitivityTests.DirectionalFit(data, degree);
        var box = new ReferenceBox([.. fit.Transform.Predictors.Select(p => ((RangeScale)p).Min)], [.. fit.Transform.Predictors.Select(p => ((RangeScale)p).Max)]);
        return (fit, data, [.. Enumerable.Range(0, 2).Select(d => new QuadratureReference($"draw {d}", SampledSensitivityTests.HandSurface(data, degree, d), box, panels))]);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void DirectionalFit_SampledMeetsQuadrature_WithinTheStatedError(int degree)
    {
        const int Seeds = 40, N = 2048;
        (BassFit fit, BassData data, QuadratureReference[] truth) = HandDirectional(degree, 48);
        var grid = new EffectSpec { GridPoints = 6 };
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        EffectFunctions? first = null;
        double[][][] truthCells = [];
        foreach (int seed in SeedTree.Derive(37, Seeds))
        {
            EffectFunctions sampled = BassSensitivity.Effects(fit, data, new SensitivitySpec { Points = N, Seed = seed }, grid).At[0];
            if (first is null)
            {
                first = sampled;
                Assert.Equal(6, sampled.Effects.Count);   // three mains and the three pairs the hinges hold
                truthCells = [.. truth.Select(q => sampled.Effects.Select(e => CellsOf(e).Select(at => q.Effect(e.Inputs, at)).ToArray()).ToArray())];
            }
            runs.Add(Cells(sampled, d => truth[d].Mean, (d, k, c) => truthCells[d][k][c]));
        }
        double scale = truthCells.SelectMany(x => x).SelectMany(x => x).Max(Math.Abs);
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, cell => cell.EndsWith(" mean") ? Math.Abs(truth[int.Parse(cell.Split(' ')[1])].Mean) + scale : scale);
        Report($"directional, degree {degree}", t);
        Assert.True(t.WorstBias <= SampledSensitivityTests.BiasCut(t.Cells), $"{t.WorstBiasCell}: √R·|mean z| {t.WorstBias:F2} against {SampledSensitivityTests.BiasCut(t.Cells):F2}");
        Assert.True(t.PooledSd >= SampledSensitivityTests.PooledLow && t.PooledSd <= SampledSensitivityTests.PooledHigh, $"pooled sd of z {t.PooledSd:F3}");
        Assert.True(t.Beyond3 <= 0.01, $"{t.Beyond3:P2} beyond 3");
    }

    // Two inputs on [0, 1], a single index the chain follows with ridge blocks of K = 4 columns.
    private static readonly Lazy<(BassFit Fit, BassData Data)> Ridge = new(() =>
    {
        const int n = 120;
        var rng = new Xoshiro256PlusPlus(36);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = Math.Exp(0.6 * x0[i] + 0.8 * x1[i]) + Math.Sin(3.0 * x0[i]) + 0.05 * Normal.SampleStandard(rng);
        }
        var data = new BassData([new ContinuousPredictor("x0", x0, new InputDomain(0.0, 1.0)), new ContinuousPredictor("x1", x1, new InputDomain(0.0, 1.0))],
                                new ScalarResponse(y));
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6, Directional = new DirectionalBasis { Columns = 4 } },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 3655, BurnIn = 600, Samples = 6, Thin = 20 },
        };
        return (Bass.Fit(config, data), data);
    });

    // The ridge fit's mean and effects per draw by tensor Gauss–Legendre quadrature of its predicted surface: 3 nodes on each of
    // `panels` panels per input, the cells' pinned values against every node of the other input.
    private static (double[] Mean, double[][][] Cells) RidgeTruth(BassFit fit, BassData data, EffectFunctions shape, int panels)
    {
        QuadratureRule rule = QuadratureRule.Uniform(0.0, 1.0, 3, panels);
        int m = rule.Nodes.Length, draws = shape.Mean.Length;
        double[,,] Predict(double[] a, double[] b) => fit.Predict(data, PredictionCases.Of([new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)]),
                                                                  new BassPredictOptions { Predictive = false }).Latent;
        double[] x0 = [.. Enumerable.Range(0, m * m).Select(k => rule.Nodes[k / m])], x1 = [.. Enumerable.Range(0, m * m).Select(k => rule.Nodes[k % m])];
        double[,,] full = Predict(x0, x1);
        var mean = new double[draws];
        for (int d = 0; d < draws; d++)
            for (int k = 0; k < m * m; k++) mean[d] += rule.Weights[k / m] * rule.Weights[k % m] * full[d, k, 0];

        // E[f | x_v = g] at each grid value against the other input's nodes; the pair's cells at their points.
        double[] grid = shape.Effects[0].Grids[0].Values;
        int gm = grid.Length;
        var given = new double[2][,];
        for (int v = 0; v < 2; v++)
        {
            double[] a = [.. Enumerable.Range(0, gm * m).Select(k => v == 0 ? grid[k / m] : rule.Nodes[k % m])];
            double[] b = [.. Enumerable.Range(0, gm * m).Select(k => v == 0 ? rule.Nodes[k % m] : grid[k / m])];
            double[,,] at = Predict(a, b);
            given[v] = new double[draws, gm];
            for (int d = 0; d < draws; d++)
                for (int k = 0; k < gm * m; k++) given[v][d, k / m] += rule.Weights[k % m] * at[d, k, 0];
        }
        double[,,] pair = Predict([.. Enumerable.Range(0, gm * gm).Select(k => grid[k / gm])], [.. Enumerable.Range(0, gm * gm).Select(k => grid[k % gm])]);
        var cells = new double[draws][][];
        for (int d = 0; d < draws; d++)
        {
            cells[d] = new double[shape.Effects.Count][];
            for (int k = 0; k < shape.Effects.Count; k++)
            {
                int[] u = shape.Effects[k].Inputs;
                cells[d][k] = u.Length == 1
                    ? [.. Enumerable.Range(0, gm).Select(g => given[u[0]][d, g] - mean[d])]
                    : [.. Enumerable.Range(0, gm * gm).Select(c => pair[d, c, 0] - given[0][d, c / gm] - given[1][d, c % gm] + mean[d])];
            }
        }
        return (mean, cells);
    }

    [Fact]
    public void RidgeBlockFit_SampledMeetsQuadratureOfItsSurface_WithinTheStatedError()
    {
        const int Seeds = 40, N = 2048;
        (BassFit fit, BassData data) = Ridge.Value;
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws.SelectMany(d => d.Structure.TermIds), id => Assert.NotNull(fit.Terms[id].Ridge));
        Assert.Contains(draws, d => d.Structure.TermIds.Length > 0);
        var grid = new EffectSpec { GridPoints = 7, Subsets = [["x0"], ["x1"], ["x0", "x1"]] };
        EffectFunctions shape = BassSensitivity.Effects(fit, data, new SensitivitySpec { Points = 64 }, grid).At[0];
        (double[] mean, double[][][] truth) = RidgeTruth(fit, data, shape, 64);
        (double[] coarseMean, double[][][] coarse) = RidgeTruth(fit, data, shape, 32);
        double quadrature = truth.Zip(coarse).SelectMany(t => t.First.Zip(t.Second).SelectMany(e => e.First.Zip(e.Second, (a, b) => Math.Abs(a - b)))).Max();
        quadrature = Math.Max(quadrature, mean.Zip(coarseMean, (a, b) => Math.Abs(a - b)).Max());

        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(38, Seeds))
        {
            EffectAnalysis analysis = BassSensitivity.Effects(fit, data, new SensitivitySpec { Points = N, Seed = seed }, grid);
            Assert.Equal(SensitivityRoute.Sampled, analysis.Route);
            runs.Add(Cells(analysis.At[0], d => mean[d], (d, k, c) => truth[d][k][c]));
        }
        double scale = truth.SelectMany(x => x).SelectMany(x => x).Max(Math.Abs);
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, cell => cell.EndsWith(" mean") ? Math.Abs(mean[int.Parse(cell.Split(' ')[1])]) + scale : scale);
        Report($"ridge blocks, {draws.Length} draws, up to {draws.Max(d => d.Structure.TermIds.Length)} columns", t);
        output.WriteLine($"the quadrature moves by at most {quadrature:E2} on halving its panels, against effects up to {scale:E2}");
        Assert.True(quadrature <= 1e-3 * scale, $"quadrature at {quadrature:E2}");
        Assert.True(t.WorstBias <= SampledSensitivityTests.BiasCut(t.Cells), $"{t.WorstBiasCell}: √R·|mean z| {t.WorstBias:F2} against {SampledSensitivityTests.BiasCut(t.Cells):F2}");
        Assert.True(t.PooledSd >= SampledSensitivityTests.PooledLow && t.PooledSd <= SampledSensitivityTests.PooledHigh, $"pooled sd of z {t.PooledSd:F3}");
        Assert.True(t.Beyond3 <= 0.01, $"{t.Beyond3:P2} beyond 3");
    }

    // ── The grid's mean, the route and the points ───────────────────────────────────────────────────────────────────────────

    [Fact]
    public void GridMean_MeetsTheSurfacesMean()
    {
        // Under a uniform law on a continuous input the grid's trapezoid weights are the law's, so the main effect's weighted mean over
        // the grid is 0 but for the grid's own quadrature error; under level probabilities it is 0 exactly. The sampled effect's weighted
        // mean is read against the truth's, within four of its weighted errors (an upper bound on its own error), on the separable fit
        // under the boxes, both routes, and on the hand directional fit; and each draw's mean against the truth's within its error.
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var grid = new EffectSpec { GridPoints = 25, Subsets = [["x0"], ["x1"], ["z"]] };
        EffectFunctions exact = BassSensitivity.Effects(fit, data, new SensitivitySpec { Route = SensitivityRoute.Exact }, grid).At[0];
        EffectFunctions sampled = BassSensitivity.Effects(fit, data, new SensitivitySpec { Route = SensitivityRoute.Sampled, Points = 4096 }, grid).At[0];
        double worst = 0.0;
        for (int d = 0; d < sampled.Mean.Length; d++)
        {
            for (int k = 0; k < 3; k++)
            {
                double[] w = Weights(sampled.Effects[k]);
                double est = 0.0, truth = 0.0, bound = 0.0;
                for (int c = 0; c < w.Length; c++)
                {
                    est += w[c] * sampled.Effects[k].Values[d, c];
                    truth += w[c] * exact.Effects[k].Values[d, c];
                    bound += w[c] * sampled.Effects[k].Errors![d, c];
                }
                if (bound == 0.0) Assert.Equal(truth, est, 12);
                else worst = Math.Max(worst, Math.Abs(est - truth) / bound);
                if (sampled.Effects[k].Grids[0].Levels is not null) Assert.True(Math.Abs(truth) <= 1e-12 * (1.0 + Math.Abs(exact.Mean[d])), $"draw {d}: {truth:E2}");
            }
            Assert.True(Math.Abs(sampled.Mean[d] - exact.Mean[d]) <= 4.0 * sampled.MeanError![d] + 1e-12, $"draw {d}: mean {sampled.Mean[d]:R} against {exact.Mean[d]:R} ± {sampled.MeanError[d]:E2}");
        }
        output.WriteLine($"the separable fit: the grid mean's worst |sampled − exact| over its bound {worst:F3}");
        Assert.True(worst <= 4.0, $"{worst:F3}");

        (BassFit hand, BassData handData, QuadratureReference[] q) = HandDirectional(2, 48);
        EffectFunctions directional = BassSensitivity.Effects(hand, handData, new SensitivitySpec { Points = 4096 }, new EffectSpec { GridPoints = 25, Subsets = [["x0"], ["x1"], ["x2"]] }).At[0];
        worst = 0.0;
        for (int d = 0; d < 2; d++)
            for (int k = 0; k < 3; k++)
            {
                Effect e = directional.Effects[k];
                double[] w = Weights(e);
                double est = 0.0, truth = 0.0, bound = 0.0;
                double[][] at = [.. CellsOf(e)];
                for (int c = 0; c < w.Length; c++)
                {
                    est += w[c] * e.Values[d, c];
                    truth += w[c] * q[d].Effect(e.Inputs, at[c]);
                    bound += w[c] * e.Errors![d, c];
                }
                worst = Math.Max(worst, Math.Abs(est - truth) / bound);
            }
        output.WriteLine($"the directional fit: the grid mean's worst |sampled − quadrature| over its bound {worst:F3}");
        Assert.True(worst <= 4.0, $"{worst:F3}");
    }

    [Fact]
    public void Auto_IsTheExactRoute_BitForBit_OnASeparableFit()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var grid = new EffectSpec { GridPoints = 9 };
        EffectAnalysis auto = BassSensitivity.Effects(fit, data, new SensitivitySpec { Measures = MixedLaws }, grid);
        EffectAnalysis exact = BassSensitivity.Effects(fit, data, new SensitivitySpec { Measures = MixedLaws, Route = SensitivityRoute.Exact }, grid);
        Assert.Equal(SensitivityRoute.Exact, auto.Route);
        Assert.Equal(0, auto.Points);
        Assert.Equal(exact.At[0].Mean, auto.At[0].Mean);
        Assert.Null(auto.At[0].MeanError);
        for (int k = 0; k < exact.At[0].Effects.Count; k++)
        {
            Assert.Equal(exact.At[0].Effects[k].Values, auto.At[0].Effects[k].Values);
            Assert.Null(auto.At[0].Effects[k].Errors);
        }
    }

    [Fact]
    public void Points_AreTheSobolRoutes_SharedByEveryDraw()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var spec = new SensitivitySpec { Measures = MixedLaws, Route = SensitivityRoute.Sampled, Points = 256, Seed = 9 };
        // Named, since the default pairs follow the selected draws' terms.
        var grid = new EffectSpec { GridPoints = 5, Subsets = [["x0"], ["z"], ["x0", "x1"], ["x1", "z"]] };
        EffectFunctions all = BassSensitivity.Effects(fit, data, spec, grid).At[0];
        Assert.Equal(all.Mean, BassSensitivity.Effects(fit, data, spec, grid).At[0].Mean);
        Assert.NotEqual(all.Mean, BassSensitivity.Effects(fit, data, spec with { Seed = 10 }, grid).At[0].Mean);
        foreach (int d in new[] { 0, 31, all.Mean.Length - 1 })
        {
            EffectFunctions one = BassSensitivity.Effects(fit, data, spec with { Draws = [d] }, grid).At[0];
            Assert.Equal(all.Mean[d], one.Mean[0]);
            for (int k = 0; k < all.Effects.Count; k++)
                for (int c = 0; c < all.Effects[k].Values.GetLength(1); c++)
                {
                    Assert.Equal(all.Effects[k].Values[d, c], one.Effects[k].Values[0, c]);
                    Assert.Equal(all.Effects[k].Errors![d, c], one.Effects[k].Errors![0, c]);
                }
        }
    }

    [Fact]
    public void Refusals_NameTheirReason()
    {
        (BassFit hand, BassData data, _) = HandDirectional(1, 4);
        ArgumentException exact = Assert.Throws<ArgumentException>(() => hand.Effects(data, new SensitivitySpec { Route = SensitivityRoute.Exact }));
        Assert.Contains("not a product of one-variable factors", exact.Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => hand.Effects(data, new SensitivitySpec { Points = 1 }));
        ArgumentException budget = Assert.Throws<ArgumentException>(() => hand.Effects(data, effects: new EffectSpec { EvaluationBudget = 1_000_000 }));
        Assert.Contains("term values", budget.Message);
        ArgumentException later = Assert.Throws<ArgumentException>(() =>
            BassSensitivity.Effects(SensitivityFixture.FitFunctional(), SensitivityFixture.Functional(), new SensitivitySpec { Route = SensitivityRoute.Sampled }));
        Assert.Contains("slice 3", later.Message);
    }

    // An effect's cells' pinned values, the last input fastest.
    private static IEnumerable<double[]> CellsOf(Effect e)
    {
        int cells = e.Values.GetLength(1);
        for (int c = 0; c < cells; c++)
        {
            var at = new double[e.Inputs.Length];
            int rest = c;
            for (int t = e.Inputs.Length - 1; t >= 0; t--)
            {
                at[t] = e.Grids[t].Values[rest % e.Grids[t].Values.Length];
                rest /= e.Grids[t].Values.Length;
            }
            yield return at;
        }
    }

    // A main effect's grid weights under its law: the trapezoid rule's on a continuous grid (the uniform law's), the level probabilities
    // the fixture's default law gives a categorical one (uniform).
    private static double[] Weights(Effect e)
    {
        int m = e.Grids[0].Values.Length;
        if (e.Grids[0].Levels is not null) return [.. Enumerable.Repeat(1.0 / m, m)];
        return [.. Enumerable.Range(0, m).Select(g => (g == 0 || g == m - 1 ? 0.5 : 1.0) / (m - 1))];
    }
}
