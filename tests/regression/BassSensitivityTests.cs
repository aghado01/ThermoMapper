using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Direct sensitivity (O6, <see cref="BassSensitivity.Of"/>) on fitted draws:
/// <list type="bullet">
/// <item>the Sobol identities on every draw: the main variances over all orders sum to the total, each input's total sums
/// the subsets that hold it, and totals are at least the mains;</item>
/// <item>each index of a fitted draw against a pick-freeze Monte Carlo estimate under a non-uniform product measure, within
/// five standard errors;</item>
/// <item>the two functional forms against each other by the law of total variance (the integrated form under the empirical
/// output grid against the conditioned form on that grid), and the conditioned variance against a tensor quadrature of the
/// basis's columns;</item>
/// <item>draw order, weights, zero variance (undefined, never 0 or 1), clamping, the order, the subset budget and the
/// refusals.</item>
/// </list>
/// </summary>
public sealed class BassSensitivityTests
{
    private static readonly Dictionary<string, InputMeasure> Measures = new()
    {
        ["x0"] = new MixtureMeasure([new NormalComponent(0.5, 0.35, 0.15), new NormalComponent(0.5, 0.75, 0.25)]),
        ["x1"] = new UniformMeasure(0.1, 0.9),
        ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }),
    };

    [Fact]
    public void SobolIdentities_HoldOnEveryDraw()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        SobolAnalysis analysis = BassSensitivity.Of(fit, SensitivityFixture.Scalar(), new SensitivitySpec { Measures = Measures });
        SobolDecomposition s = Assert.Single(analysis.Decompositions);
        Assert.Null(analysis.Form);
        Assert.Equal(new[] { "x0", "x1", "z" }, s.Inputs);
        Assert.Equal(s.Subsets.Max(u => u.Length), analysis.MaxOrder);
        int defined = 0;
        for (int d = 0; d < s.Count; d++)
        {
            if (!s.IsDefined(d)) continue;
            defined++;
            double v = s.Variance[d], sum = 0.0;
            for (int u = 0; u < s.Subsets.Count; u++) sum += s.MainVariance[d, u];
            Assert.True(Math.Abs(sum - v) <= 1e-12 * v, $"draw {d}: mains sum to {sum}, total {v}");
            for (int i = 0; i < s.Inputs.Count; i++)
            {
                double holding = Enumerable.Range(0, s.Subsets.Count).Where(u => s.Subsets[u].Contains(i)).Sum(u => s.MainVariance[d, u]);
                Assert.True(Math.Abs(s.TotalVariance[d, i] - holding) <= 1e-12 * v, $"draw {d}, input {i}");
                int single = IndexOf(s, [i]);
                if (single >= 0) Assert.True(s.TotalIndex(d, i) >= s.MainIndex(d, single) - 1e-12);
            }
            Assert.Equal(DecompositionStatus.Defined, s.Status[d] & ~DecompositionStatus.Clamped);
        }
        Assert.True(defined > s.Count / 2);
        double[] means = s.MeanMainIndices();
        Assert.Equal(1.0, means.Sum(), 10);
    }

    [Fact]
    public void Indices_MatchPickFreezeMonteCarloOnAFittedDraw()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        SobolDecomposition s = BassSensitivity.Of(fit, data, new SensitivitySpec { Measures = Measures }).Decompositions[0];
        // The draw with the most terms, an interaction among them.
        int d = Enumerable.Range(0, s.Count).OrderByDescending(k => DrawOf(fit, s, k).Structure.TermIds.Length).First();
        BassDraw draw = DrawOf(fit, s, d);
        var basis = new BassBasis(BassInputs.Prepare(data, fit.Transform), 1);

        const int n = 200_000;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(2026, 1)[0]);
        double[][] a = Sample(rng, n), b = Sample(rng, n);
        double[] fa = Surface(fit, basis, draw, a), fb = Surface(fit, basis, draw, b);
        // V from the two independent samples.
        double mean = (fa.Sum() + fb.Sum()) / (2 * n);
        double[] sq = [.. fa.Concat(fb).Select(f => (f - mean) * (f - mean))];
        (double v, double vSe) = MeanSe(sq);
        Assert.True(Math.Abs(v - s.Variance[d]) <= 5 * vSe, $"V: Monte Carlo {v} ± {vSe}, exact {s.Variance[d]}");
        for (int i = 0; i < 3; i++)
        {
            double[][] ab = [.. a.Select((col, k) => k == i ? b[k] : col)];
            double[] fab = Surface(fit, basis, draw, ab);
            (double main, double mainSe) = MeanSe([.. fb.Select((f, r) => f * (fab[r] - fa[r]))]);           // Saltelli 2010
            (double total, double totalSe) = MeanSe([.. fa.Select((f, r) => 0.5 * (f - fab[r]) * (f - fab[r]))]);   // Jansen
            int single = IndexOf(s, [i]);
            double exactMain = single >= 0 ? s.MainVariance[d, single] : 0.0;
            Assert.True(Math.Abs(main - exactMain) <= 5 * mainSe + 1e-12, $"main {i}: Monte Carlo {main} ± {mainSe}, exact {exactMain}");
            Assert.True(Math.Abs(total - s.TotalVariance[d, i]) <= 5 * totalSe + 1e-12, $"total {i}: Monte Carlo {total} ± {totalSe}, exact {s.TotalVariance[d, i]}");
        }
        // Each pair's closed variance Var(E[f | x_i, x_k]) = V_i + V_k + V_ik.
        for (int i = 0; i < 3; i++)
            for (int k = i + 1; k < 3; k++)
            {
                double[][] ab = [.. a.Select((col, c) => c == i || c == k ? b[c] : col)];
                double[] fab = Surface(fit, basis, draw, ab);
                (double closed, double se) = MeanSe([.. fb.Select((f, r) => f * (fab[r] - fa[r]))]);
                double exact = new[] { new[] { i }, new[] { k }, new[] { i, k } }.Select(u => IndexOf(s, u)).Where(u => u >= 0).Sum(u => s.MainVariance[d, u]);
                Assert.True(Math.Abs(closed - exact) <= 5 * se + 1e-12, $"closed ({i}, {k}): Monte Carlo {closed} ± {se}, exact {exact}");
            }
        Assert.Contains(s.Subsets, u => u.Length == 2 && IndexOf(s, u) >= 0 && s.MainVariance[d, IndexOf(s, u)] > 0.0);
    }

    [Fact]
    public void FunctionalForms_ObeyTheLawOfTotalVariance()
    {
        BassFit fit = SensitivityFixture.FitFunctional();
        BassData data = SensitivityFixture.Functional();
        var grid = new Dictionary<string, InputMeasure> { ["t"] = new EmpiricalMeasure() };
        SobolAnalysis integrated = BassSensitivity.Of(fit, data, new SensitivitySpec { OutputMeasures = grid });
        SobolAnalysis conditioned = BassSensitivity.Of(fit, data, new SensitivitySpec { Functional = FunctionalForm.Condition });
        EffectAnalysis means = BassSensitivity.Effects(fit, data, new SensitivitySpec { Functional = FunctionalForm.Condition }, new EffectSpec { Subsets = [] });
        SobolDecomposition whole = Assert.Single(integrated.Decompositions);
        Assert.Equal(FunctionalForm.Integrate, integrated.Form);
        Assert.Equal(new[] { "x0", "x1", "t" }, whole.Inputs);
        Assert.Equal(7, conditioned.Decompositions.Count);
        Assert.Equal(7, conditioned.Coordinates!.Count);
        Assert.All(conditioned.CoordinatesOutside!, Assert.False);
        Assert.Equal(new[] { "x0", "x1" }, conditioned.Decompositions[0].Inputs);
        for (int d = 0; d < whole.Count; d++)
        {
            double within = conditioned.Decompositions.Average(c => c.Variance[d]);
            double grand = means.At.Average(e => e.Mean[d]);
            double between = means.At.Average(e => (e.Mean[d] - grand) * (e.Mean[d] - grand));
            Assert.True(Math.Abs(whole.Variance[d] - (within + between)) <= 1e-11 * Math.Max(whole.Variance[d], 1e-12),
                $"draw {d}: integrated {whole.Variance[d]}, within {within} + between {between}");
            foreach (SobolDecomposition part in conditioned.Decompositions.Append(whole))
            {
                if (!part.IsDefined(d)) continue;
                double sum = 0.0;
                for (int u = 0; u < part.Subsets.Count; u++) sum += part.MainVariance[d, u];
                Assert.True(Math.Abs(sum - part.Variance[d]) <= 1e-12 * part.Variance[d]);
                for (int i = 0; i < part.Inputs.Count; i++)
                {
                    int single = IndexOf(part, [i]);
                    if (single >= 0) Assert.True(part.TotalVariance[d, i] >= part.MainVariance[d, single] - 1e-12 * part.Variance[d]);
                }
            }
        }
        Assert.Contains(whole.Subsets, u => u.SequenceEqual(new[] { 2 }));   // t is an input with main and total indices
    }

    [Fact]
    public void ConditionedVariance_MatchesTensorQuadratureAtEachCoordinate()
    {
        BassFit fit = SensitivityFixture.FitFunctional();
        BassData data = SensitivityFixture.Functional();
        double[][] coordinates = [[0.0], [0.4], [1.0], [1.3]];
        SobolAnalysis conditioned = BassSensitivity.Of(fit, data, new SensitivitySpec
        {
            Functional = FunctionalForm.Condition, Coordinates = coordinates,
            Measures = new Dictionary<string, InputMeasure> { ["x1"] = new UniformMeasure(0.2, 0.8) },
        });
        Assert.Equal(new[] { false, false, false, true }, conditioned.CoordinatesOutside!);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var basis = new BassBasis(inputs, 1);
        var reference = new SensitivityReference(inputs, conditioned.Measure, 1);
        (double[] x0, double[] w0) = reference.Nodes(0);
        (double[] x1, double[] w1) = reference.Nodes(1);
        var cases = new CaseInputs([[.. x0.SelectMany(u => x1.Select(_ => u))], [.. x0.SelectMany(_ => x1)]], []);
        double[] weight = [.. w0.SelectMany(u => w1.Select(v => u * v))];
        var outputs = new OutputInputs([[.. coordinates.Select(c => fit.Transform.OutputVariables[0].Scale(c[0]))]]);
        SobolDecomposition first = conditioned.Decompositions[0];
        foreach (int d in new[] { 0, first.Count / 2, first.Count - 1 })
        {
            BassDraw draw = DrawOf(fit, first, d);
            var f = new double[coordinates.Length][];
            for (int r = 0; r < coordinates.Length; r++) f[r] = Enumerable.Repeat(draw.Coefficients[0], cases.Cases).ToArray();
            var column = new double[cases.Cases];
            var side = new double[coordinates.Length];
            for (int j = 0; j < draw.Structure.TermIds.Length; j++)
            {
                BassTerm term = fit.Terms[draw.Structure.TermIds[j]];
                basis.CaseColumn(term, cases, column);
                basis.OutputColumn(term, outputs, side);
                for (int r = 0; r < coordinates.Length; r++)
                    for (int i = 0; i < cases.Cases; i++) f[r][i] += draw.Coefficients[1 + j] * side[r] * column[i];
            }
            for (int r = 0; r < coordinates.Length; r++)
            {
                double m = f[r].Zip(weight, (u, w) => u * w).Sum();
                double v = f[r].Zip(weight, (u, w) => (u - m) * (u - m) * w).Sum();
                Assert.True(Math.Abs(conditioned.Decompositions[r].Variance[d] - v) <= 1e-11 * Math.Max(v, 1e-9),
                    $"draw {d}, coordinate {coordinates[r][0]}: {conditioned.Decompositions[r].Variance[d]} against {v}");
            }
        }
    }

    [Fact]
    public void Draws_KeepChainOrder_AndTheirWeights()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        SobolDecomposition s = BassSensitivity.Of(fit, data).Decompositions[0];
        int at = 0;
        for (int c = 0; c < fit.Draws.Count; c++)
            for (int i = 0; i < fit.Draws[c].Count; i++, at++)
            {
                Assert.Equal(c, s.Chains[at]);
                Assert.Equal(i, s.Draws[at]);
            }
        Assert.Equal(at, s.Count);
        Assert.Null(s.Weights);

        double[] weights = [.. Enumerable.Range(0, at).Select(k => 1.0 + k % 3)];
        SobolDecomposition w = BassSensitivity.Of(fit with { DrawWeights = weights }, data).Decompositions[0];
        double[] mean = w.MeanTotalIndices();
        for (int i = 0; i < w.Inputs.Count; i++)
        {
            double num = 0.0, den = 0.0;
            for (int d = 0; d < at; d++)
            {
                if (!w.IsDefined(d)) continue;
                num += weights[d] * w.TotalIndex(d, i);
                den += weights[d];
            }
            Assert.Equal(num / den, mean[i], 13);
        }
        Assert.Throws<ArgumentException>(() => BassSensitivity.Of(fit with { DrawWeights = [1.0] }, data));
    }

    [Fact]
    public void ZeroVariance_LeavesTheIndicesUndefined()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var table = new BassTermTable();
        int dead = table.Intern(new BassTerm([new HingeFactor(0, 1, 6)]));   // (x0 − 2/3)₊: zero on [0, 0.5]
        int live = table.Intern(new BassTerm([new HingeFactor(1, 1, 2)]));
        var store = new DrawStore<BassDraw>();
        store.Add(new BassDraw(new BassStructure([]), 1.0, 1.0, [2.0], 1.0));
        store.Add(new BassDraw(new BassStructure([dead]), 1.0, 1.0, [2.0, 5.0], 1.0));
        store.Add(new BassDraw(new BassStructure([dead, live]), 1.0, 1.0, [2.0, 5.0, 1.0], 1.0));
        BassFit synthetic = fit with { Terms = table, Draws = [store], DrawWeights = null };
        var spec = new SensitivitySpec { Measures = new Dictionary<string, InputMeasure> { ["x0"] = new UniformMeasure(0.0, 0.5) } };
        SobolDecomposition s = BassSensitivity.Of(synthetic, data, spec).Decompositions[0];
        Assert.Equal(new[] { DecompositionStatus.ZeroVariance, DecompositionStatus.ZeroVariance }, s.Status.Take(2));
        Assert.True(double.IsNaN(s.MainIndex(0, 0)) && double.IsNaN(s.TotalIndex(1, 0)));
        Assert.True(s.IsDefined(2));
        Assert.Equal(0.0, s.TotalIndex(2, 0));
        Assert.Equal(1.0, s.TotalIndex(2, 1), 14);
        Assert.Equal(1.0, s.MeanTotalIndices()[1], 14);   // the undefined draws are left out, not counted as 0
        double[] none = BassSensitivity.Of(synthetic with { Draws = [Only(store, 1)] }, data, spec).Decompositions[0].MeanMainIndices();
        Assert.True(double.IsNaN(Assert.Single(none)));   // no draw defined: undefined, not 0
    }

    [Fact]
    public void Status_ClampsWithinTheTolerance_AndFlagsBeyondIt()
    {
        var main = new double[,] { { 1.0, -1e-12, 0.5 } };
        var total = new double[,] { { 1.0, 0.5 } };
        Assert.Equal(DecompositionStatus.Clamped, BassSensitivity.Status(1.5, 3.0, 1e-10, main, total, 0));
        Assert.Equal(0.0, main[0, 1]);
        main[0, 1] = -1e-3;
        Assert.Equal(DecompositionStatus.NegativeBeyondTolerance, BassSensitivity.Status(1.5, 3.0, 1e-10, main, total, 0));
        Assert.Equal(-1e-3, main[0, 1]);
        Assert.Equal(DecompositionStatus.ZeroVariance, BassSensitivity.Status(1e-14, 3.0, 1e-10, main, total, 0));
        Assert.Equal(DecompositionStatus.ZeroVariance, BassSensitivity.Status(0.0, 0.0, 1e-10, main, total, 0));
    }

    [Fact]
    public void OrderAndBudget_BoundTheSubsets()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        SobolDecomposition full = BassSensitivity.Of(fit, data).Decompositions[0];
        Assert.True(full.Subsets.Count > 3);
        SobolAnalysis first = BassSensitivity.Of(fit, data, new SensitivitySpec { MaxOrder = 1 });
        Assert.Equal(1, first.MaxOrder);
        Assert.All(first.Decompositions[0].Subsets, u => Assert.Single(u));
        // The totals do not depend on the reported order.
        Assert.Equal(full.TotalVariance, first.Decompositions[0].TotalVariance);
        var ex = Assert.Throws<ArgumentException>(() => BassSensitivity.Of(fit, data, new SensitivitySpec { SubsetBudget = 2 }));
        Assert.Contains("budget", ex.Message);
    }

    [Fact]
    public void Refusals_NameTheirReason()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        Assert.Contains("other data", Assert.Throws<ArgumentException>(() => BassSensitivity.Of(fit, SensitivityFixture.Scalar(seed: 8))).Message);
        Assert.Throws<ArgumentException>(() => BassSensitivity.Of(fit, data, new SensitivitySpec { Coordinates = [[0.5]] }));
        Assert.Throws<ArgumentException>(() => BassSensitivity.Of(fit, data, new SensitivitySpec
        {
            OutputMeasures = new Dictionary<string, InputMeasure> { ["t"] = new UniformMeasure() },
        }));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassSensitivity.Of(fit, data, new SensitivitySpec { RoundingTolerance = -1.0 }));
        BassFit functional = SensitivityFixture.FitFunctional();
        Assert.Throws<ArgumentException>(() => BassSensitivity.Of(functional, SensitivityFixture.Functional(), new SensitivitySpec { Coordinates = [[0.5]] }));
        Assert.Throws<ArgumentException>(() => BassSensitivity.Of(functional, SensitivityFixture.Functional(), new SensitivitySpec
        {
            Functional = FunctionalForm.Condition, Coordinates = [[0.5, 0.1]],
        }));
    }

    internal static int IndexOf(SobolDecomposition s, int[] subset)
    {
        for (int u = 0; u < s.Subsets.Count; u++)
            if (s.Subsets[u].SequenceEqual(subset)) return u;
        return -1;
    }

    internal static BassDraw DrawOf(BassFit fit, SobolDecomposition s, int d) => fit.Draws[s.Chains[d]][s.Draws[d]];

    private static DrawStore<BassDraw> Only(DrawStore<BassDraw> store, int index)
    {
        var one = new DrawStore<BassDraw>();
        one.Add(store[index]);
        return one;
    }

    // N draws of (x0, x1, z) from Measures: the normal mixture truncated to [0, 1] by rejection (weights before
    // truncation), the uniform on [0.1, 0.9], the levels by their probabilities.
    private static double[][] Sample(Xoshiro256PlusPlus rng, int n)
    {
        var x0 = new double[n];
        var x1 = new double[n];
        var z = new double[n];
        Normal c1 = Normal.MeanSd(0.35, 0.15), c2 = Normal.MeanSd(0.75, 0.25);
        for (int i = 0; i < n; i++)
        {
            double x;
            do x = rng.NextDouble() < 0.5 ? c1.Sample(rng) : c2.Sample(rng);
            while (x < 0.0 || x > 1.0);
            x0[i] = x;
            x1[i] = 0.1 + 0.8 * rng.NextDouble();
            double u = rng.NextDouble();
            z[i] = u < 0.2 ? 0 : u < 0.7 ? 1 : 2;
        }
        return [x0, x1, z];
    }

    private static double[] Surface(BassFit fit, BassBasis basis, BassDraw draw, double[][] x)
    {
        var cases = new CaseInputs([x[0], x[1]], [[.. x[2].Select(v => (int)v)]]);
        var f = Enumerable.Repeat(draw.Coefficients[0], cases.Cases).ToArray();
        var column = new double[cases.Cases];
        for (int j = 0; j < draw.Structure.TermIds.Length; j++)
        {
            basis.CaseColumn(fit.Terms[draw.Structure.TermIds[j]], cases, column);
            for (int i = 0; i < f.Length; i++) f[i] += draw.Coefficients[1 + j] * column[i];
        }
        return f;
    }

    private static (double Mean, double Se) MeanSe(double[] values)
    {
        double mean = values.Average();
        double var = values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1);
        return (mean, Math.Sqrt(var / values.Length));
    }
}
