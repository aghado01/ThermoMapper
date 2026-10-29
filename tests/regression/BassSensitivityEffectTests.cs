using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Effect functions (O8, <see cref="BassSensitivity.Effects"/>) on a degree-2 fit under a non-uniform product measure: each
/// main effect is E[f | x_v] − f₀ and each two-way effect E[f | x_v, x_w] − f_v − f_w − f₀, against a tensor quadrature of
/// the basis's own columns over the other inputs, with f₀ = E[f] likewise; the default effects, the grids, and the refusals.
/// R hard-codes the uniform measure and degree 1 here (the O-O2 defect).
/// </summary>
public sealed class BassSensitivityEffectTests
{
    private static readonly Dictionary<string, InputMeasure> Measures = new()
    {
        ["x0"] = new MixtureMeasure([new NormalComponent(0.5, 0.35, 0.15), new StudentComponent(0.5, 3.0, 0.75, 0.2)]),
        ["x1"] = new UniformMeasure(0.1, 0.9),
        ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }),
    };

    [Fact]
    public void Effects_AreConditionalExpectationsLessTheirLowerOrders()
    {
        BassFit fit = SensitivityFixture.FitScalar(degree: 2);
        BassData data = SensitivityFixture.Scalar();
        var spec = new SensitivitySpec { Measures = Measures };
        EffectAnalysis analysis = BassSensitivity.Effects(fit, data, spec, new EffectSpec
        {
            Subsets = [["x0"], ["x1"], ["z"], ["x0", "x1"], ["z", "x1"]], GridPoints = 6,
        });
        EffectFunctions e = Assert.Single(analysis.At);
        Assert.Equal(new[] { 1, 2 }, e.Effects[4].Inputs);   // sorted
        Assert.Equal(new[] { 0.1, 0.26, 0.42, 0.58, 0.74, 0.9 }, e.Effects[1].Grids[0].Values.Select(v => Math.Round(v, 12)));
        Assert.Equal(new[] { "a", "b", "c" }, e.Effects[2].Grids[0].Levels!);

        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var basis = new BassBasis(inputs, 2);
        var reference = new SensitivityReference(inputs, analysis.Measure, 2);
        var nodes = Enumerable.Range(0, 3).Select(reference.Nodes).ToArray();
        int count = fit.Draws.Sum(c => c.Count);
        foreach (int d in new[] { 0, count / 3, count - 1 })
        {
            BassDraw draw = Draw(fit, d);
            double Mean(params (int Variable, double Value)[] fixedAt) => ConditionalMean(fit, basis, inputs, draw, nodes, fixedAt);
            double f0 = Mean();
            Assert.True(Math.Abs(e.Mean[d] - f0) <= 1e-11 * Math.Max(1.0, Math.Abs(f0)), $"f0 of draw {d}: {e.Mean[d]} against {f0}");
            var main = new Func<double, double>[3];
            for (int v = 0; v < 3; v++)
            {
                int variable = v;
                main[v] = x => Mean((variable, x)) - f0;
                Effect effect = e.Effects[v];
                double[] grid = Scaled(effect.Grids[0]);
                for (int c = 0; c < grid.Length; c++)
                {
                    double expected = main[v](grid[c]);
                    Assert.True(Math.Abs(effect.Values[d, c] - expected) <= 1e-11 * Math.Max(1.0, Math.Abs(f0)),
                        $"main effect of {v} at {grid[c]}, draw {d}: {effect.Values[d, c]} against {expected}");
                }
            }
            foreach (Effect effect in e.Effects.Skip(3))
            {
                double[] g0 = Scaled(effect.Grids[0]), g1 = Scaled(effect.Grids[1]);
                int v0 = effect.Inputs[0], v1 = effect.Inputs[1];
                for (int c0 = 0; c0 < g0.Length; c0++)
                    for (int c1 = 0; c1 < g1.Length; c1++)
                    {
                        double expected = Mean((v0, g0[c0]), (v1, g1[c1])) - main[v0](g0[c0]) - main[v1](g1[c1]) - f0;
                        Assert.True(Math.Abs(effect.Values[d, c0 * g1.Length + c1] - expected) <= 1e-11 * Math.Max(1.0, Math.Abs(f0)),
                            $"two-way effect of ({v0}, {v1}) at ({g0[c0]}, {g1[c1]}), draw {d}");
                    }
            }
        }
    }

    [Fact]
    public void DefaultEffects_AreTheMainsAndThePosteriorsPairs()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        EffectFunctions e = BassSensitivity.Effects(fit, data).At[0];
        Assert.Equal(new[] { 0, 1, 2 }, e.Effects.Take(3).Select(x => x.Inputs.Single()));
        var pairs = new HashSet<(int, int)>();
        var moments = BassMoments.On(BassInputs.Prepare(data, fit.Transform), 1, BassSensitivity.Of(fit, data).Measure);
        foreach (var chain in fit.Draws)
            foreach (BassDraw draw in chain)
                foreach (int id in draw.Structure.TermIds)
                {
                    IReadOnlyList<int> v = moments.Variables(fit.Terms[id]);
                    for (int i = 0; i < v.Count; i++)
                        for (int k = i + 1; k < v.Count; k++) pairs.Add((v[i], v[k]));
                }
        Assert.Equal(pairs.Order(), e.Effects.Skip(3).Select(x => (x.Inputs[0], x.Inputs[1])));
        Assert.Equal(100, e.Effects[0].Grids[0].Values.Length);
        Assert.Equal(fit.Draws.Sum(c => c.Count), e.Mean.Length);
    }

    [Fact]
    public void Refusals_NameTheirReason()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        Assert.Contains("main or two-way", Assert.Throws<ArgumentException>(() =>
            BassSensitivity.Effects(fit, data, effects: new EffectSpec { Subsets = [["x0", "x1", "z"]] })).Message);
        Assert.Contains("'w'", Assert.Throws<ArgumentException>(() =>
            BassSensitivity.Effects(fit, data, effects: new EffectSpec { Subsets = [["w"]] })).Message);
        Assert.Contains("twice", Assert.Throws<ArgumentException>(() =>
            BassSensitivity.Effects(fit, data, effects: new EffectSpec { Subsets = [["x0", "x0"]] })).Message);
        Assert.Contains("budget", Assert.Throws<ArgumentException>(() =>
            BassSensitivity.Effects(fit, data, effects: new EffectSpec { CellBudget = 1000 })).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => BassSensitivity.Effects(fit, data, effects: new EffectSpec { GridPoints = 1 }));
    }

    private static BassDraw Draw(BassFit fit, int d)
    {
        foreach (var chain in fit.Draws)
        {
            if (d < chain.Count) return chain[d];
            d -= chain.Count;
        }
        throw new ArgumentOutOfRangeException(nameof(d));
    }

    // A grid in the fit's scaled units: the identity transform here, and level indices for a categorical input.
    private static double[] Scaled(EffectGrid grid) => grid.Values;

    // E[f | the fixed inputs] by the tensor rule over the others.
    private static double ConditionalMean(BassFit fit, BassBasis basis, BassInputs inputs, BassDraw draw, (double[] X, double[] W)[] nodes,
                                          (int Variable, double Value)[] fixedAt)
    {
        var each = nodes.Select((n, v) => fixedAt.Any(f => f.Variable == v) ? (X: new[] { fixedAt.First(f => f.Variable == v).Value }, W: new[] { 1.0 }) : n).ToArray();
        int total = each.Aggregate(1, (acc, n) => acc * n.X.Length);
        var x0 = new double[total];
        var x1 = new double[total];
        var z = new int[total];
        var w = new double[total];
        for (int i = 0; i < total; i++)
        {
            int a = i / (each[1].X.Length * each[2].X.Length), b = i / each[2].X.Length % each[1].X.Length, c = i % each[2].X.Length;
            x0[i] = each[0].X[a];
            x1[i] = each[1].X[b];
            z[i] = (int)each[2].X[c];
            w[i] = each[0].W[a] * each[1].W[b] * each[2].W[c];
        }
        var cases = new CaseInputs([x0, x1], [z]);
        var f = Enumerable.Repeat(draw.Coefficients[0], total).ToArray();
        var column = new double[total];
        for (int j = 0; j < draw.Structure.TermIds.Length; j++)
        {
            basis.CaseColumn(fit.Terms[draw.Structure.TermIds[j]], cases, column);
            for (int i = 0; i < total; i++) f[i] += draw.Coefficients[1 + j] * column[i];
        }
        return f.Zip(w, (u, v) => u * v).Sum();
    }
}
