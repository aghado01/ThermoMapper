using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 5, categorical inputs, through the readouts: a gated fit's term counts, its prediction against its draws' columns by hand,
/// the surface and E8's direction against the function, from outside the model and beside BASS's hinges on the same data; the calibration
/// emulator's split on a categorical parameter that gates a direction; Sobol's sampled route on f(C, X) = C·X, where C enters only as a gate,
/// against its closed form, and E8's Poincaré check on that fit.
/// </summary>
public sealed class CategoricalReadoutTests(ITestOutputHelper output)
{
    // y = exp(0.6x₀ + 0.8x₁) where c ≠ 1, plus 1 where c = 2, and 0.5 where c = 1, noise sd 0.05; three inputs uniform on [0, 1] and c of
    // three levels: a ridge gated to {a, c} and a level shift on {c} that no intercept absorbs.
    private static (BassData Data, double[] Truth) Gated(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[3][];
        for (int j = 0; j < 3; j++) x[j] = new double[n];
        int[] c = new int[n];
        var y = new double[n];
        var truth = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 3; j++) x[j][i] = rng.NextDouble();
            c[i] = rng.NextInt(3);
            truth[i] = c[i] != 1 ? Math.Exp(0.6 * x[0][i] + 0.8 * x[1][i]) + (c[i] == 2 ? 1.0 : 0.0) : 0.5;
            y[i] = truth[i] + 0.05 * Normal.SampleStandard(rng);
        }
        return (new BassData([.. Enumerable.Range(0, 3).Select(j => (BassPredictor)new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0))),
                              new CategoricalPredictor("c", c, ["a", "b", "c"])], new ScalarResponse(y)), truth);
    }

    private static double Rmse(BassFit fit, BassData data, double[] truth)
    {
        double[,,] latent = fit.Predict(data, PredictionCases.Of(data)).Latent;
        int draws = latent.GetLength(0);
        return Math.Sqrt(Enumerable.Range(0, truth.Length).Average(i =>
        {
            double f = Enumerable.Range(0, draws).Average(d => latent[d, i, 0]);
            return (f - truth[i]) * (f - truth[i]);
        }));
    }

    // A column by hand: its cells' mask times slice 1's hinge, or slice 3's paper column on the hand projection.
    private static double[] HandColumn(BassData data, BassTerm term)
    {
        if (term.Ridge is not { } r) return CategoricalRidgeTests.HandColumn(data, term);
        double[] mask = CategoricalRidgeTests.Mask(CategoricalRidgeTests.Levels(data), term.Categorical);
        double[] projection = DirectionalChainTests.Projections(CategoricalRidgeTests.Continuous(data), r.Direction);
        return [.. projection.Select((u, i) => mask[i] * RidgeColumnTests.Paper([.. r.Knots], r.Column, u).Value)];
    }

    [Theory]
    [InlineData("chain", 4)]
    [InlineData("smc", 4)]
    [InlineData("chain", 1)]
    public void AGatedFit_CountsItsTerms_PredictsAsItsColumns_AndFollowsTheFunction(string sampler, int columns)
    {
        (BassData data, double[] truth) = Gated(300, 3731);
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 10, Directional = new DirectionalBasis { Columns = columns } },
            Sampling = new BassSampling { Chains = 2, MasterSeed = 3732, BurnIn = 6000, Samples = 300, Thin = 2 },
        };
        if (sampler == "smc")
            config = config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 200, Replicates = 2, MutationSteps = 12, TargetConditionalEss = 0.9 } } };
        BassFit fit = Bass.Fit(config, data);
        BassModel model = BassModel.Create(config, data);

        // M counts structure terms through the table: a ridge block's K columns once, a term of categorical factors alone once.
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        int gatedDirections = 0, categoricalAlone = 0;
        foreach (BassDraw draw in pooled)
        {
            BassTerm[] terms = [.. draw.Structure.TermIds.Select(id => fit.Terms[id])];
            IReadOnlyList<BassTerm> units = model.Expansion is { } expansion ? expansion.Units(terms) : terms;
            Assert.Equal(units.Count, draw.Structure.TermCount(fit.Terms));
            gatedDirections += units.Count(u => u.Directional is not null && u.Categorical.Count > 0);
            categoricalAlone += units.Count(u => u.Directional is null);
        }
        double[] weights = fit.DrawWeights is { } w ? [.. w] : [.. Enumerable.Repeat(1.0 / pooled.Length, pooled.Length)];
        double meanM = pooled.Select((d, i) => weights[i] * d.Structure.TermCount(fit.Terms)).Sum() / weights.Sum();
        Assert.Equal(meanM, fit.Summary().Posterior.Count.Mean, 10);
        if (sampler != "smc")
        {
            BassTraces traces = fit.Traces();
            for (int c = 0; c < fit.Draws.Count; c++)
                for (int i = 0; i < fit.Draws[c].Count; i++) Assert.Equal(fit.Draws[c][i].Structure.TermCount(fit.Terms), traces.Chains[c].Count[i]);
        }

        // The prediction at the training cases against each draw's columns by hand, for every tenth draw.
        double[,,] latent = fit.Predict(data, PredictionCases.Of(data)).Latent;
        double worst = 0.0;
        for (int d = 0; d < pooled.Length; d += 10)
        {
            int[] ids = pooled[d].Structure.TermIds;
            double[][] hand = [.. ids.Select(id => HandColumn(data, fit.Terms[id]))];
            for (int i = 0; i < truth.Length; i++)
            {
                double f = pooled[d].Coefficients[0];
                for (int j = 0; j < ids.Length; j++) f += pooled[d].Coefficients[j + 1] * hand[j][i];
                worst = Math.Max(worst, Math.Abs(f - latent[d, i, 0]) / Math.Max(1.0, Math.Abs(f)));
            }
        }

        // The surface and E8's leading direction against the function's, beside BASS's hinges on the same data.
        double rmse = Rmse(fit, data, truth);
        double spread = Math.Sqrt(truth.Average(t => (t - truth.Average()) * (t - truth.Average())));
        ActiveSubspaceAnalysis analysis = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Points = 1024 });
        double[,] vectors = analysis.Posterior.MeanCVectors;
        double[,] leading = new double[3, 1], target = new double[3, 1];
        double[] direction = [0.6, 0.8, 0.0];
        for (int j = 0; j < 3; j++) { leading[j, 0] = vectors[j, 0]; target[j, 0] = direction[j]; }
        double angle = SubspacePosterior.PrincipalAngles(leading, target)[0];
        BassFit hinges = Bass.Fit(config with { Structure = config.Structure with { Directional = null } }, data);
        output.WriteLine($"{sampler}, K = {columns}: mean M {meanM:F2}; gated directions {gatedDirections}, categorical terms alone {categoricalAlone} over "
                         + $"{pooled.Length} draws; prediction against the hand columns {worst:E2}; surface RMSE {rmse:F4} of a spread {spread:F3} "
                         + $"(BASS's hinges {Rmse(hinges, data, truth):F4}); E8's angle to (0.6, 0.8, 0) {angle:F4} rad");
        Assert.True(gatedDirections > 0 && categoricalAlone > 0, "the fit drew no gated direction or no categorical term alone");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.True(rmse <= 0.05 * spread, $"{rmse:F4}");
        Assert.True(angle <= 0.1, $"{angle:F4}");
    }

    // ---- the calibration emulator's split --------------------------------------------------------------------------------------------

    [Fact]
    public void TheEmulatorsSplitOnACategoricalParameter_KeepsTheDirectionItGates()
    {
        (BassData data, _) = Gated(200, 3733);
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 8, Directional = new DirectionalBasis() },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 3734, BurnIn = 3000, Samples = 40 },
        };
        BassFit fit = Bass.Fit(config, data);
        int gated = fit.Draws.SelectMany(c => c).Sum(d => d.Structure.TermIds.Count(id => fit.Terms[id] is { Directional: not null, Categorical.Count: > 0 }));
        Assert.True(gated > 0, "no draw holds a gated direction, so the split is not reached");
        CategoricalParameter c = new("c", ["b", "c", "a"]);   // declared in another order than the training levels
        var emulator = new BassFitEmulator(fit, data, [c, new ContinuousParameter("x0", 0.0, 1.0), new ContinuousParameter("x1", 0.0, 1.0), new ContinuousParameter("x2", 0.0, 1.0)]);
        var rng = new Xoshiro256PlusPlus(3735);
        var outputs = new double[1];
        var enumerated = new double[3, 1];
        int[][] all = [[0], [1], [2]];
        int compared = 0;
        for (int trial = 0; trial < 20; trial++)
        {
            int s = rng.NextInt(emulator.DrawCount);
            double[] t = [rng.NextDouble(), rng.NextDouble(), rng.NextDouble()];
            emulator.PredictLevels(s, t, all, enumerated);
            for (int level = 0; level < 3; level++)
            {
                emulator.Predict(s, t, [level], outputs);
                var cases = new PredictionCases([new ContinuousColumn("x0", [t[0]]), new ContinuousColumn("x1", [t[1]]), new ContinuousColumn("x2", [t[2]]),
                                                 new CategoricalColumn("c", [c.Levels[level]])]);
                BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Draws = [s], Predictive = false });
                Assert.Equal(p.Latent[0, 0, 0], outputs[0]);
                Assert.Equal(outputs[0], enumerated[level, 0]);
                compared++;
            }
        }
        output.WriteLine($"gated directions held across the draws {gated}; {compared} emulator readings equal the fit's prediction bit for bit");
    }

    // ---- Sobol's sampled route on a gate ------------------------------------------------------------------------------------------------

    // f(C, X) = C·X with C ~ Bernoulli(½) and X ~ U(0, 1), independent: Var f = 5/48, S_C = 3/5, T_C = 4/5, S_X = 1/5, T_X = 2/5. C enters
    // only as a gate: a term reading it holds a direction too, which no separable factor carries.
    private static BassData GateOnly(int n = 200)
    {
        var rng = new Xoshiro256PlusPlus(3736);
        double[] x = [.. Enumerable.Range(0, n).Select(i => (0.5 + i * 0.6180339887498949) % 1.0)];
        int[] c = [.. Enumerable.Range(0, n).Select(i => i % 2)];
        double[] y = [.. Enumerable.Range(0, n).Select(i => c[i] * x[i] + 0.01 * Normal.SampleStandard(rng))];
        return new BassData([new ContinuousPredictor("x", x, new InputDomain(0.0, 1.0)), new CategoricalPredictor("c", c, ["off", "on"])], new ScalarResponse(y));
    }

    [Fact]
    public void Sobol_OnAGateOnlyCategory_VariesTheGate_AndMeetsTheClosedForm_AndPoincareHolds()
    {
        BassData data = GateOnly();
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 6, Directional = new DirectionalBasis() },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 3737, BurnIn = 3000, Samples = 40 },
        };
        BassFit fit = Bass.Fit(config, data);
        SobolAnalysis analysis = fit.Sensitivity(data, new SensitivitySpec { Route = SensitivityRoute.Auto, Points = 4096 });
        Assert.Equal(SensitivityRoute.Sampled, analysis.Route);
        SobolDecomposition s = analysis.Decompositions[0];
        int xv = analysis.Measure.OfContinuous(0), cv = analysis.Measure.OfCategorical(0);
        int sx = Enumerable.Range(0, s.Subsets.Count).Single(k => s.Subsets[k].SequenceEqual([xv]));
        int sc = Enumerable.Range(0, s.Subsets.Count).Single(k => s.Subsets[k].SequenceEqual([cv]));
        int[] defined = [.. Enumerable.Range(0, s.Variance.Length).Where(s.IsDefined)];
        double Mean(Func<int, double> f) => defined.Average(f);
        (string Name, double Engine, double Truth)[] rows =
        [
            ("S_C", Mean(d => s.MainIndex(d, sc)), 0.6), ("T_C", Mean(d => s.TotalIndex(d, cv)), 0.8),
            ("S_X", Mean(d => s.MainIndex(d, sx)), 0.2), ("T_X", Mean(d => s.TotalIndex(d, xv)), 0.4),
        ];
        foreach ((string name, double engine, double truth) in rows) output.WriteLine($"{name}: sampled {engine:F4}, closed form {truth:F4}");
        // Before the gate's variable was read, C went unvaried and S_C read 0.
        Assert.All(rows, r => Assert.True(Math.Abs(r.Engine - r.Truth) <= 0.03, $"{r.Name} {r.Engine:F4} against {r.Truth:F4}"));

        // E8's Poincaré check on the same fit: totals with their error present, no violation.
        PoincareCheck check = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Points = 2048 }).Poincare!;
        Assert.NotNull(check.TotalError);
        output.WriteLine($"Poincaré: violations {check.Violations}");
        Assert.Equal(0, check.Violations);
    }
}
