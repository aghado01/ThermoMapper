using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The completion audit's probes (issues/regression/bass/completion-audit-brief.md, part 1's table and part 4's edge
/// cases): what the hub's contract names and no earlier test exercised end to end.
/// <list type="bullet">
/// <item>Every cell of the response representations by the input kinds: prediction at new inputs (and output coordinates)
/// against a hand evaluation of the draws, sensitivity and effects against enumeration and a midpoint grid, reconstruction
/// through K, and restart against the uninterrupted run, bit for bit.</item>
/// <item>The edge cases: more than 64 levels, the intercept-only posterior against its closed form, one predictor, one
/// draw, two cases, full PCA rank through the fit, a zero-probability level and a truncated measure through a fit.</item>
/// <item>What the earlier fixtures left out: degrees 2 and 3 and two output variables through a direct fit, a measure wider
/// than the training range, and the residual R calls <c>other</c>.</item>
/// </list>
/// References stand outside the code under test: the hinge is written out here, and the moments are a brute-force sum.
/// </summary>
public sealed class BassCompletionAuditProbes
{
    private readonly ITestOutputHelper _out;
    public BassCompletionAuditProbes(ITestOutputHelper output) => _out = output;

    public enum Inputs { Continuous, Categorical, Mixed }

    public enum Response { Scalar, Direct, Pca, Supplied }

    private const int Cases = 30, Locations = 6, Outputs = 10;

    // ── The cells ──────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double[] X0() => [.. Enumerable.Range(0, Cases).Select(i => (i % 6) / 5.0)];

    private static double[] X1() => [.. Enumerable.Range(0, Cases).Select(i => (i * 7 % 5) / 4.0)];

    private static int[] Z0() => [.. Enumerable.Range(0, Cases).Select(i => i % 3)];

    private static int[] Z1() => [.. Enumerable.Range(0, Cases).Select(i => (i / 3) % 2)];

    private static BassPredictor[] Predictors(Inputs inputs) => inputs switch
    {
        Inputs.Continuous => [new ContinuousPredictor("x0", X0()), new ContinuousPredictor("x1", X1())],
        Inputs.Categorical => [new CategoricalPredictor("z0", Z0(), ["a", "b", "c"]), new CategoricalPredictor("z1", Z1(), ["p", "q"])],
        _ => [new ContinuousPredictor("x0", X0()), new CategoricalPredictor("z0", Z0(), ["a", "b", "c"])],
    };

    private static double Surface(Inputs inputs, int i)
    {
        double x0 = X0()[i], x1 = X1()[i];
        int z0 = Z0()[i], z1 = Z1()[i];
        return inputs switch
        {
            Inputs.Continuous => 1 + 2 * x0 + 1.5 * Math.Max(x1 - 0.4, 0) + 2 * x0 * x1,
            Inputs.Categorical => 1 + (z0 == 1 ? 1.5 : 0) - (z1 == 1 ? 0.7 : 0) + (z0 == 2 && z1 == 0 ? 1 : 0),
            _ => 1 + 2 * x0 + (z0 == 1 ? 1.5 * x0 : 0),
        };
    }

    private static double[,] Curves(Inputs inputs, int m, int seed = 8)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var y = new double[Cases, m];
        for (int i = 0; i < Cases; i++)
            for (int l = 0; l < m; l++)
            {
                double t = l / (m - 1.0);
                y[i, l] = Surface(inputs, i) * (1 + t) + Math.Sin(3 * t) * (i % 4) / 3.0 + 0.1 * Normal.SampleStandard(rng);
            }
        return y;
    }

    private static BassData Data(Inputs inputs, Response response)
    {
        if (response == Response.Scalar)
        {
            var rng = new Xoshiro256PlusPlus(5);
            return new BassData(Predictors(inputs), new ScalarResponse([.. Enumerable.Range(0, Cases).Select(i => Surface(inputs, i) + 0.1 * Normal.SampleStandard(rng))]));
        }
        double[] t = [.. Enumerable.Range(0, Locations).Select(l => l / (Locations - 1.0))];
        return new BassData(Predictors(inputs), new FunctionalResponse(Curves(inputs, Locations), [new ContinuousPredictor("t", t)]));
    }

    private static BassConfig Config(int samples = 24, int chains = 2, int burn = 120, int maxBasis = 6)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxCategoricalInteraction = 2, MaxBasis = maxBasis },
            Sampling = new BassSampling { Chains = chains, MasterSeed = 11, BurnIn = burn, Samples = samples },
        };

    private static BassDecomposition Decomposition(Inputs inputs, Response response)
    {
        double[,] y = Curves(inputs, Outputs);
        if (response == Response.Pca) return BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2 });

        // A nonorthogonal basis fixed in advance, its coefficients the caller's least squares: two curves with a shared part.
        var k = new double[Outputs, 2];
        for (int l = 0; l < Outputs; l++)
        {
            double t = l / (Outputs - 1.0);
            k[l, 0] = 1 + t;
            k[l, 1] = 1 + t + Math.Sin(3 * t);
        }
        double a = 0, b = 0, c = 0;
        for (int l = 0; l < Outputs; l++)
        {
            a += k[l, 0] * k[l, 0];
            b += k[l, 0] * k[l, 1];
            c += k[l, 1] * k[l, 1];
        }
        double det = a * c - b * b;
        var w = new double[2, Cases];
        for (int i = 0; i < Cases; i++)
        {
            double u = 0, v = 0;
            for (int l = 0; l < Outputs; l++)
            {
                u += k[l, 0] * y[i, l];
                v += k[l, 1] * y[i, l];
            }
            w[0, i] = (c * u - b * v) / det;
            w[1, i] = (a * v - b * u) / det;
        }
        return BassDecomposition.FromBasis(k, w, new double[Outputs], [.. Enumerable.Repeat(1.0, Outputs)], responses: y);
    }

    private static PredictionCases NewCases(Inputs inputs) => inputs switch
    {
        Inputs.Continuous => new PredictionCases([new ContinuousColumn("x1", [0.13, 0.62, 0.95]), new ContinuousColumn("x0", [0.21, 0.5, 0.77])]),
        Inputs.Categorical => new PredictionCases([new CategoricalColumn("z1", ["q", "p", "q"]), new CategoricalColumn("z0", ["c", "a", "b"])]),
        _ => new PredictionCases([new CategoricalColumn("z0", ["c", "a", "b"]), new ContinuousColumn("x0", [0.21, 0.5, 0.77])]),
    };

    // The new cases, as each kind's positions read them: scaled by hand, since these predictors' ranges are [0, 1].
    private static (double[][] X, int[][] Z) NewScaled(Inputs inputs) => inputs switch
    {
        Inputs.Continuous => ([[0.21, 0.5, 0.77], [0.13, 0.62, 0.95]], []),
        Inputs.Categorical => ([], [[2, 0, 1], [1, 0, 1]]),
        _ => ([[0.21, 0.5, 0.77]], [[2, 0, 1]]),
    };

    // ── The hand evaluation: R's normalized hinge, written out ─────────────────────────────────────────────────────────

    private static double Hinges(IEnumerable<HingeFactor> factors, Func<int, double> value, Func<int, int, double> knot, int degree)
    {
        double f = 1.0, peak = 1.0;
        foreach (HingeFactor h in factors)
        {
            double t = knot(h.Variable, h.Knot);
            f *= Math.Pow(Math.Max(h.Sign * (value(h.Variable) - t), 0.0), degree);
            peak *= Math.Pow((h.Sign + 1) / 2.0 - h.Sign * t, degree);
        }
        return peak != 0.0 ? f / peak : f;
    }

    // f(x, z, t) of one draw at scaled inputs. With several output variables, `outputs` holds them; `t` is the one's.
    private static double Hand(BassFit fit, BassInputs inputs, BassDraw draw, double[] x, int[] z, double t, double[]? outputs = null)
    {
        int degree = fit.Config.Structure.Degree;
        double f = draw.Coefficients[0];
        int[] ids = draw.Structure.TermIds;
        for (int j = 0; j < ids.Length; j++)
        {
            BassTerm term = fit.Terms[ids[j]];
            double b = Hinges(term.Continuous, v => x[v], (v, k) => inputs.ContinuousKnots(v)[k], degree);
            foreach (CategoricalFactor c in term.Categorical)
                if (!c.Levels.Contains(z[c.Variable])) b = 0.0;
            if (term.Functional.Count > 0) b *= Hinges(term.Functional, v => outputs?[v] ?? t, (v, k) => inputs.FunctionalKnots(v)[k], degree);
            f += draw.Coefficients[j + 1] * b;
        }
        return f;
    }

    private static List<BassDraw> Pooled(BassFit fit) => [.. fit.Draws.SelectMany(c => c)];

    private static int Richest(IEnumerable<int> draws, Func<int, int> terms) => draws.OrderByDescending(terms).First();

    // ── The grid: every input by name, its nodes and weights, each input's marginal summing to 1 ───────────────────────

    private sealed record Node(double Value, int Level, double Weight);

    private static Node[] Continuous(int points, double lower = 0.0, double upper = 1.0)
        => [.. Enumerable.Range(0, points).Select(g => new Node(lower + (g + 0.5) / points * (upper - lower), -1, 1.0 / points))];

    private static Node[] Levels(params double[] probabilities) => [.. probabilities.Select((p, l) => new Node(0.0, l, p))];

    // Var f, f's mean, and each input's Var E[f | input] over the product of `dims`, f taking one node per input.
    private static (double V, double Mean, double[] Main) Moments(Node[][] dims, Func<Node[], double> f)
    {
        int count = dims.Aggregate(1, (a, d) => a * d.Length);
        var nodes = new Node[dims.Length];
        var at = new int[dims.Length];
        var conditional = dims.Select(d => new double[d.Length]).ToArray();
        double sum = 0.0, square = 0.0;
        for (int cell = 0; cell < count; cell++)
        {
            int rest = cell;
            double w = 1.0;
            for (int v = 0; v < dims.Length; v++)
            {
                at[v] = rest % dims[v].Length;
                rest /= dims[v].Length;
                nodes[v] = dims[v][at[v]];
                w *= nodes[v].Weight;
            }
            double value = f(nodes);
            sum += w * value;
            square += w * value * value;
            for (int v = 0; v < dims.Length; v++) conditional[v][at[v]] += w * value;
        }
        var main = new double[dims.Length];
        for (int v = 0; v < dims.Length; v++)
            for (int a = 0; a < dims[v].Length; a++)
                if (dims[v][a].Weight > 0.0)
                {
                    double given = conditional[v][a] / dims[v][a].Weight;
                    main[v] += dims[v][a].Weight * (given - sum) * (given - sum);
                }
        return (square - sum * sum, sum, main);
    }

    // The grid over `predictors`, uniform on each, the continuous midpoints `points` a side.
    private static Node[][] Grid(IReadOnlyList<BassPredictor> predictors, int points)
        => [.. predictors.Select(p => p is CategoricalPredictor c ? Levels([.. c.Levels.Select(_ => 1.0 / c.Levels.Count)]) : Continuous(points))];

    private static int SubsetIndex(SobolDecomposition s, int input)
    {
        for (int u = 0; u < s.Subsets.Count; u++)
            if (s.Subsets[u].Length == 1 && s.Subsets[u][0] == input) return u;
        return -1;
    }

    private static void AgreeWithGrid(string where, SobolDecomposition s, int draw, (double V, double Mean, double[] Main) grid, double tolerance)
    {
        double slack = tolerance * Math.Max(grid.V, 1e-12);
        Assert.True(Math.Abs(s.Variance[draw] - grid.V) <= slack, $"{where}: variance {s.Variance[draw]} against the grid's {grid.V}");
        for (int v = 0; v < grid.Main.Length; v++)
        {
            int u = SubsetIndex(s, v);
            double main = u >= 0 ? s.MainVariance[draw, u] : 0.0;
            Assert.True(Math.Abs(main - grid.Main[v]) <= slack, $"{where}: input {v}'s main variance {main} against the grid's {grid.Main[v]}");
        }
    }

    // ── A cell, direct: scalar or functional ───────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(Inputs.Continuous, Response.Scalar)]
    [InlineData(Inputs.Categorical, Response.Scalar)]
    [InlineData(Inputs.Mixed, Response.Scalar)]
    [InlineData(Inputs.Continuous, Response.Direct)]
    [InlineData(Inputs.Categorical, Response.Direct)]
    [InlineData(Inputs.Mixed, Response.Direct)]
    public void Cell_Direct_PredictsDecomposesAndRestarts(Inputs kind, Response response)
    {
        BassData data = Data(kind, response);
        BassFit fit = Bass.Fit(Config(), data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        bool direct = response == Response.Direct;

        // Prediction at new inputs and, for a direct fit, new output coordinates, against the hand evaluation.
        double[] t = [0.1, 0.45, 0.8];
        (double[][] nx, int[][] nz) = NewScaled(kind);
        BassPrediction p = fit.Predict(data, NewCases(kind), new BassPredictOptions
        {
            Predictive = false,
            OutputLocations = direct ? [new ContinuousColumn("t", t)] : null,
        });
        Assert.Equal((pooled.Count, 3, direct ? 3 : 1), (p.DrawCount, p.Cases, p.Locations));
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 3; i++)
                for (int l = 0; l < p.Locations; l++)
                {
                    double hand = Hand(fit, inputs, pooled[d], [.. nx.Select(c => c[i])], [.. nz.Select(c => c[i])], direct ? t[l] : 0.0);
                    Assert.True(Math.Abs(hand - p.Latent[d, i, l]) <= 1e-12, $"draw {d}, case {i}, location {l}: {p.Latent[d, i, l]} against {hand}");
                }

        // Sensitivity on the richest draw: its variance and main variances against the grid under the default measures;
        // the integrated form takes the output coordinate as an input.
        int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
        Assert.True(pooled[richest].Structure.TermIds.Length > 0, "the chains leave the intercept");
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data).Decompositions);
        Node[][] grid = direct ? [.. Grid(data.Predictors, 30), Continuous(30)] : Grid(data.Predictors, 200);
        int pc = data.Predictors.Count(q => q is ContinuousPredictor);
        double Eval(Node[] nodes)
        {
            double[] x = [.. Enumerable.Range(0, data.Predictors.Count).Where(v => data.Predictors[v] is ContinuousPredictor).Select(v => nodes[v].Value)];
            int[] z = [.. Enumerable.Range(0, data.Predictors.Count).Where(v => data.Predictors[v] is CategoricalPredictor).Select(v => nodes[v].Level)];
            return Hand(fit, inputs, pooled[richest], x, z, direct ? nodes[^1].Value : 0.0);
        }
        (double V, double Mean, double[] Main) moments = Moments(grid, Eval);
        bool exact = pc == 0 && !direct;
        AgreeWithGrid($"{kind}, {response}", s, richest, moments, exact ? 1e-11 : 0.02);
        for (int d = 0; d < s.Count; d++)
        {
            if (!s.IsDefined(d)) continue;
            double total = 0.0;
            for (int u = 0; u < s.Subsets.Count; u++) total += s.MainVariance[d, u];
            Assert.True(Math.Abs(total - s.Variance[d]) <= 1e-10 * s.Variance[d], "the components sum to the variance");
        }

        // Effects: each main effect at its grid, against the conditional mean from the same grid, and the draw's mean.
        EffectFunctions effects = fit.Effects(data, null, new EffectSpec { GridPoints = 5 }).At[0];
        Assert.True(Math.Abs(effects.Mean[richest] - moments.Mean) <= (exact ? 1e-11 : 0.02) * Math.Max(Math.Abs(moments.Mean), 1.0));
        foreach (Effect e in effects.Effects.Where(e => e.Inputs.Length == 1))
        {
            int v = e.Inputs[0];
            for (int c = 0; c < e.Grids[0].Values.Length; c++)
            {
                double at = e.Grids[0].Values[c];
                Node[][] given = [.. grid];
                given[v] = e.Grids[0].Levels is null ? [new Node(at, -1, 1.0)] : [new Node(0.0, (int)at, 1.0)];
                double conditional = Moments(given, Eval).Mean;
                double expected = conditional - moments.Mean;
                Assert.True(Math.Abs(e.Values[richest, c] - expected) <= (exact ? 1e-11 : 0.03) * Math.Max(Math.Abs(moments.Mean), 1.0),
                            $"{kind}, {response}: effect of input {v} at {at}: {e.Values[richest, c]} against {expected}");
            }
        }

        // Restart: the checkpoint continues the run bit for bit.
        BassConfig half = Config(samples: 12);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassRunTests.AssertSameFit(Bass.Fit(Config(samples: 24), data), Bass.Resume(first.Checkpoint!, half, data));
        _out.WriteLine($"{kind}, {response}: {fit.Terms.Count} terms; richest draw {pooled[richest].Structure.TermIds.Length} terms, V {moments.V:G5}");
    }

    // ── A cell, through a basis ────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(Inputs.Continuous, Response.Pca)]
    [InlineData(Inputs.Categorical, Response.Pca)]
    [InlineData(Inputs.Mixed, Response.Pca)]
    [InlineData(Inputs.Continuous, Response.Supplied)]
    [InlineData(Inputs.Categorical, Response.Supplied)]
    [InlineData(Inputs.Mixed, Response.Supplied)]
    public void Cell_Basis_ReconstructsDecomposesAndRestarts(Inputs kind, Response response)
    {
        BassDecomposition decomposition = Decomposition(kind, response);
        BassPredictor[] predictors = Predictors(kind);
        BassBasisFit fit = BassComponents.Fit(decomposition, predictors, Config());
        LinearReconstruction map = decomposition.Map;
        int r = decomposition.Rank, m = decomposition.Locations;
        BassInputs inputs = BassInputs.Prepare(fit.ComponentData(0), fit.Components[0].Transform);
        List<BassDraw>[] pooled = [.. fit.Components.Select(Pooled)];

        // f at a location for joint draw j: y_m + y_s Σ_k K_lk w_k(x), each w_k by hand.
        double At(int j, int location, double[] x, int[] z)
        {
            double value = 0.0;
            for (int k = 0; k < r; k++)
                value += map.Basis(location, k) * Hand(fit.Components[k], inputs, pooled[k][fit.JointDraws.Draw(k, j)], x, z, 0.0);
            return map.Offsets[location] + map.Scales[location] * value;
        }

        // Reconstruction at new inputs.
        (double[][] nx, int[][] nz) = NewScaled(kind);
        BassBasisPrediction p = fit.Predict(NewCases(kind), new BassPredictOptions { Predictive = false });
        Assert.Equal((fit.JointDraws.Count, 3, m), (p.DrawCount, p.Cases, p.Locations));
        for (int j = 0; j < p.DrawCount; j++)
            for (int i = 0; i < 3; i++)
                for (int l = 0; l < m; l++)
                {
                    double hand = At(j, l, [.. nx.Select(c => c[i])], [.. nz.Select(c => c[i])]);
                    Assert.True(Math.Abs(hand - p.Latent[j, i, l]) <= 1e-11, $"joint draw {j}, case {i}, location {l}: {p.Latent[j, i, l]} against {hand}");
                }

        // Sensitivity at a location and integrated over the locations, on the richest joint draw, against the grid.
        int richest = Richest(Enumerable.Range(0, p.DrawCount), j => Enumerable.Range(0, r).Sum(k => pooled[k][fit.JointDraws.Draw(k, j)].Structure.TermIds.Length));
        Node[][] grid = Grid(predictors, 200);
        double Eval(Node[] nodes, int location)
        {
            double[] x = [.. Enumerable.Range(0, predictors.Length).Where(v => predictors[v] is ContinuousPredictor).Select(v => nodes[v].Value)];
            int[] z = [.. Enumerable.Range(0, predictors.Length).Where(v => predictors[v] is CategoricalPredictor).Select(v => nodes[v].Level)];
            return At(richest, location, x, z);
        }
        bool exact = predictors.All(q => q is CategoricalPredictor);
        BasisSobolAnalysis analysis = fit.Sensitivity();
        int middle = m / 2;
        (double V, double Mean, double[] Main) moments = Moments(grid, n => Eval(n, middle));
        AgreeWithGrid($"{kind}, {response}, location {middle}", analysis.AtLocation(middle), richest, moments, exact ? 1e-10 : 0.02);

        Node[][] withLocation = [.. Grid(predictors, exact ? 2 : 60), .. new[] { Levels([.. Enumerable.Repeat(1.0 / m, m)]) }];
        (double V, double Mean, double[] Main) integrated = Moments(withLocation, n => Eval(n, n[^1].Level));
        AgreeWithGrid($"{kind}, {response}, integrated", analysis.Integrated, richest, integrated, exact ? 1e-10 : 0.02);

        // Effects at the same location.
        EffectFunctions effects = fit.Effects(null, new EffectSpec { GridPoints = 5 }).AtLocation(middle);
        Assert.True(Math.Abs(effects.Mean[richest] - moments.Mean) <= (exact ? 1e-10 : 0.02) * Math.Max(Math.Abs(moments.Mean), 1.0));
        foreach (Effect e in effects.Effects.Where(e => e.Inputs.Length == 1))
        {
            int v = e.Inputs[0];
            for (int c = 0; c < e.Grids[0].Values.Length; c++)
            {
                double at = e.Grids[0].Values[c];
                Node[][] given = [.. grid];
                given[v] = e.Grids[0].Levels is null ? [new Node(at, -1, 1.0)] : [new Node(0.0, (int)at, 1.0)];
                double expected = Moments(given, n => Eval(n, middle)).Mean - moments.Mean;
                Assert.True(Math.Abs(e.Values[richest, c] - expected) <= (exact ? 1e-10 : 0.03) * Math.Max(Math.Abs(moments.Mean), 1.0),
                            $"{kind}, {response}: effect of input {v} at {at}: {e.Values[richest, c]} against {expected}");
            }
        }

        // Restart: each component's checkpoint continues it, and the fits assemble into the uninterrupted basis fit.
        BassConfig half = Config(samples: 12);
        BassBasisFit first = BassComponents.Fit(decomposition, predictors, half, new BassRunOptions { Checkpoint = true });
        BassFit[] resumed = [.. Enumerable.Range(0, r).Select(k => Bass.Resume(first.Components[k].Checkpoint!, first.Components[k].Config, first.ComponentData(k)))];
        BassBasisFit continued = BassComponents.Assemble(decomposition, predictors, resumed, half.Sampling.MasterSeed);
        BassBasisPrediction again = continued.Predict(NewCases(kind), new BassPredictOptions { Predictive = false });
        for (int k = 0; k < r; k++) BassRunTests.AssertSameFit(fit.Components[k], continued.Components[k]);
        Assert.Equal(p.Latent.Cast<double>(), again.Latent.Cast<double>());
        _out.WriteLine($"{kind}, {response}: r = {r}; richest joint draw V at location {middle} {moments.V:G5}; integrated {integrated.V:G5}");
    }

    // ── More than 64 levels ────────────────────────────────────────────────────────────────────────────────────────────

    private static BassData ManyLevels(int levels, int seed = 3)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        string[] names = [.. Enumerable.Range(0, levels).Select(l => $"l{l}")];
        int n = 2 * levels;
        int[] codes = [.. Enumerable.Range(0, n).Select(i => i % levels)];
        double[] x = [.. Enumerable.Range(0, n).Select(i => (i % 7) / 6.0)];
        double[] y = [.. Enumerable.Range(0, n).Select(i => (i % levels >= levels - 3 ? 4.0 : 1.0) + x[i] + 0.1 * Normal.SampleStandard(rng))];
        return new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", codes, names)], new ScalarResponse(y));
    }

    [Theory]
    [InlineData(70)]
    [InlineData(300)]
    public void Levels_AboveSixtyFour_DrawAndScoreUnderBothMeasures(int levels)
    {
        BassData data = ManyLevels(levels);
        foreach (CategoricalSubsetMeasure measure in new[] { CategoricalSubsetMeasure.SizeUniform, CategoricalSubsetMeasure.Uniform })
        {
            BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
            var prior = new BassStructurePrior(inputs, 1, 1, 0, measure, new BassCountPrior());
            var rng = new Xoshiro256PlusPlus(19);
            var seen = new bool[levels];
            double sizes = 0.0;
            int subsets = 0, above = 0;
            for (int draw = 0; draw < 4000; draw++)
            {
                BassTerm? term = prior.DrawTerm(rng);
                if (term is null || term.Categorical.Count == 0) continue;
                CategoricalFactor f = term.Categorical[0];
                Assert.InRange(f.Levels.Count, 1, levels - 1);
                subsets++;
                sizes += f.Levels.Count;
                foreach (int l in f.Levels)
                {
                    seen[l] = true;
                    if (l >= 64) above++;
                }

                // The term measure's log probability: the part's order, then the subset's own.
                double logSubset = measure == CategoricalSubsetMeasure.Uniform
                    ? -(levels * Math.Log(2.0) + Math.Log(1.0 - Math.Pow(2.0, 1 - levels)))
                    : -Math.Log(levels - 1) - LogChoose(levels, f.Levels.Count);
                double logRest = term.Continuous.Count == 0 ? -Math.Log(2.0) : -2 * Math.Log(2.0) - Math.Log(7.0);
                Assert.Equal(-Math.Log(2.0) + logSubset + logRest, prior.LogPrior(term), 9);
            }
            Assert.All(seen, Assert.True);
            Assert.True(above > 0, "a level above the 64th is drawn");
            double expected = levels / 2.0, spread = measure == CategoricalSubsetMeasure.Uniform ? Math.Sqrt(levels) / 2.0 : levels / Math.Sqrt(12.0);
            Assert.InRange(sizes / subsets, expected - 6 * spread / Math.Sqrt(subsets), expected + 6 * spread / Math.Sqrt(subsets));
        }
    }

    private static double LogChoose(int n, int k)
        => SpecialFunctions.LogGamma(n + 1.0) - SpecialFunctions.LogGamma(k + 1.0) - SpecialFunctions.LogGamma(n - k + 1.0);

    [Fact]
    public void Levels_AboveSixtyFour_FitPredictMeasureAndRestart()
    {
        const int levels = 70;
        BassData data = ManyLevels(levels);
        BassConfig config = Config(samples: 30, burn: 400, maxBasis: 8);
        BassFit fit = Bass.Fit(config, data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        Assert.True(pooled.Any(d => d.Structure.TermIds.Any(id => fit.Terms[id].Categorical.Any(c => c.Levels.Any(l => l >= 64)))),
                    "a retained term's subset holds a level above the 64th");

        // Prediction at levels on both sides of the 64th, against the hand evaluation.
        string[] named = ["l3", "l63", "l64", "l69"];
        int[] indices = [3, 63, 64, 69];
        double[] x = [0.2, 0.4, 0.6, 0.8];
        BassPrediction p = fit.Predict(data, new PredictionCases([new ContinuousColumn("x", x), new CategoricalColumn("z", named)]),
                                       new BassPredictOptions { Predictive = false });
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 4; i++)
                Assert.Equal(Hand(fit, inputs, pooled[d], [x[i]], [indices[i]], 0.0), p.Latent[d, i, 0], 11);

        // The moments under a declared measure: the mass on four levels, two above the 64th, one of them at zero.
        var probabilities = new Dictionary<string, double> { ["l3"] = 0.25, ["l64"] = 0.5, ["l69"] = 0.25, ["l65"] = 0.0 };
        var spec = new SensitivitySpec { Measures = new Dictionary<string, InputMeasure> { ["z"] = new CategoricalMeasure(probabilities) } };
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data, spec).Decompositions);
        int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
        var dist = new double[levels];
        dist[3] = 0.25;
        dist[64] = 0.5;
        dist[69] = 0.25;
        (double V, double Mean, double[] Main) grid = Moments([Continuous(400), Levels(dist)], n => Hand(fit, inputs, pooled[richest], [n[0].Value], [n[1].Level], 0.0));
        AgreeWithGrid("70 levels", s, richest, grid, 0.01);

        BassConfig half = config with { Sampling = config.Sampling with { Samples = 15 } };
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassRunTests.AssertSameFit(fit, Bass.Resume(first.Checkpoint!, half, data));
        _out.WriteLine($"70 levels: {fit.Terms.Count} terms; richest draw {pooled[richest].Structure.TermIds.Length}, V {grid.V:G5}");
    }

    // ── The intercept-only posterior ───────────────────────────────────────────────────────────────────────────────────

    // A support threshold of every case: no proper subset reaches it, so every birth is refused and M stays 0.
    private static (BassConfig Config, BassData Data) InterceptOnly(int samples = 4000)
    {
        BassData data = Data(Inputs.Categorical, Response.Scalar);
        BassConfig config = Config(samples: samples, burn: 200) with
        {
            Structure = new BassStructureSpec { MaxCategoricalInteraction = 1, MaxBasis = 4, Support = new SupportRule(SupportKind.Count, Cases) },
        };
        return (config, data);
    }

    [Fact]
    public void InterceptOnly_PosteriorIsItsClosedForm()
    {
        (BassConfig config, BassData data) = InterceptOnly();
        BassFit fit = Bass.Fit(config, data);
        List<BassDraw> pooled = Pooled(fit);
        Assert.All(pooled, d => Assert.Empty(d.Structure.TermIds));
        Assert.Equal(0, fit.Terms.Count);
        Assert.True(fit.Diagnostics.Refusals.Count(BassRefusal.Support) > 0);

        // With M = 0 the posterior is closed: σ² | y ~ IG(c + N/2, c·s² + RSS/2), RSS about ȳ under level centring and
        // s² = RSS/N, and a₀ | σ², τ, y centred on ȳ.
        double[] y = ((ScalarResponse)data.Response).Y.ToArray();
        double mean = y.Average(), rss = y.Sum(v => (v - mean) * (v - mean));
        const double c = 1e-3;
        double shape = c + y.Length / 2.0, scale = c * rss / y.Length + rss / 2.0;
        double expected = scale / (shape - 1.0);
        double sd = expected / Math.Sqrt(shape - 2.0);
        // Chains and independent: a₀ and σ² are drawn fresh at each retained draw, given τ.
        double se = sd / Math.Sqrt(pooled.Count);
        double sigma2 = pooled.Average(d => d.Sigma2);
        _out.WriteLine($"E[σ²]: {sigma2:F5} against {expected:F5} (± {se:F5}); E[a₀]: {pooled.Average(d => d.Coefficients[0]):F5} against {mean:F5}");
        Assert.InRange(sigma2, expected - 6 * se, expected + 6 * se);
        double level = pooled.Average(d => d.Coefficients[0]);
        Assert.InRange(level, mean - 6 * Math.Sqrt(expected / y.Length / pooled.Count), mean + 6 * Math.Sqrt(expected / y.Length / pooled.Count));
    }

    [Fact]
    public void InterceptOnly_EveryViewReadsIt()
    {
        (BassConfig config, BassData data) = InterceptOnly(samples: 40);
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { Checkpoint = true });
        List<BassDraw> pooled = Pooled(fit);
        Assert.All(pooled, d => Assert.Empty(d.Structure.TermIds));

        // Prediction: each draw's latent is its intercept at every case; the noise is its σ².
        BassPrediction p = fit.Predict(data, NewCases(Inputs.Categorical));
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 3; i++) Assert.Equal(pooled[d].Coefficients[0], p.Latent[d, i, 0]);

        // Sensitivity: every draw's variance is zero, so its indices are undefined, not zero.
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data).Decompositions);
        Assert.All(s.Status, status => Assert.Equal(DecompositionStatus.ZeroVariance, status & DecompositionStatus.ZeroVariance));
        Assert.All(s.MeanTotalIndices(), v => Assert.True(double.IsNaN(v)));
        EffectFunctions effects = fit.Effects(data, null, new EffectSpec { GridPoints = 3 }).At[0];
        for (int d = 0; d < pooled.Count; d++) Assert.Equal(pooled[d].Coefficients[0], effects.Mean[d], 12);
        Assert.All(effects.Effects.Where(e => e.Inputs.Length == 1), e => Assert.All(e.Values.Cast<double>(), v => Assert.Equal(0.0, v, 12)));

        // Presentation.
        BassSummary summary = fit.Summary();
        Assert.Equal((0.0, 0.0), (summary.Posterior.Count.Min, summary.Posterior.Count.Max));
        BassTrainingFit training = fit.TrainingFit(data);
        Assert.Equal(Cases, training.Cases);
        Assert.Equal(new[] { pooled.Count / 2 }, new[] { fit.Traces().Chains[0].Count.Length });

        // The structural diagnostics read a constant: not a claim of convergence.
        Assert.NotEqual(ConvergenceOutcome.Converged, fit.Diagnostics.Verdict.Group("structural").Outcome);

        // Restart.
        BassFit resumed = Bass.Resume(fit.Checkpoint!, config, data);
        Assert.Equal(2 * 40, resumed.Draws[0].Count);
    }

    // ── One predictor, one draw, two cases, one case ───────────────────────────────────────────────────────────────────

    [Fact]
    public void OnePredictor_FitPredictAndSensitivity()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double[] x = [.. Enumerable.Range(0, Cases).Select(i => i / (Cases - 1.0))];
        var data = new BassData([new ContinuousPredictor("x", x)], new ScalarResponse([.. x.Select(v => 1 + 2 * Math.Max(v - 0.5, 0) + 0.05 * Normal.SampleStandard(rng))]));
        BassFit fit = Bass.Fit(Config(samples: 30, burn: 300), data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        double[] at = [0.1, 0.9];
        BassPrediction p = fit.Predict(data, new PredictionCases([new ContinuousColumn("x", at)]), new BassPredictOptions { Predictive = false });
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 2; i++) Assert.Equal(Hand(fit, inputs, pooled[d], [at[i]], [], 0.0), p.Latent[d, i, 0], 12);

        // R refuses `sobol` with one input; here the one input holds the whole variance.
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data).Decompositions);
        Assert.Equal(new[] { "x" }, s.Inputs);
        int defined = 0;
        for (int d = 0; d < s.Count; d++)
        {
            if (!s.IsDefined(d)) continue;
            defined++;
            Assert.Equal(1.0, s.TotalIndex(d, 0), 12);
            Assert.Equal(1.0, s.MainIndex(d, 0), 12);
        }
        Assert.True(defined > 0);
    }

    [Fact]
    public void OneDraw_EveryViewReadsIt()
    {
        BassData data = Data(Inputs.Mixed, Response.Scalar);
        BassConfig config = Config(samples: 1, chains: 1, burn: 80);
        BassFit fit = Bass.Fit(config, data);
        Assert.Equal(1, Assert.Single(fit.Draws).Count);
        BassPrediction p = fit.Predict(data, NewCases(Inputs.Mixed));
        Assert.Equal((1, 3, 1), (p.DrawCount, p.Cases, p.Locations));
        Assert.Equal(p.Latent[0, 1, 0], p.Mean()[1, 0]);
        Assert.Equal(p.Latent[0, 1, 0], p.Quantile(0.5)[1, 0], 9);    // one Gaussian: its median is its mean
        Assert.NotNull(p.Predictive);
        Assert.Single(fit.Sensitivity(data).Decompositions[0].Variance);
        Assert.Single(fit.Effects(data, null, new EffectSpec { GridPoints = 3 }).At);
        Assert.Equal(1, fit.Summary().Metadata.Lengths.Draws);
        Assert.Single(fit.Traces().Chains[0].Count);
        Assert.Equal(Cases, fit.TrainingFit(data).Cases);
        // One chain of one draw has no R̂ to read: undetermined or failing, never converged.
        Assert.NotEqual(ConvergenceOutcome.Converged, fit.Diagnostics.Verdict.Outcome);
        BassPointwiseLogDensityShape(fit.PointwiseLogDensity(data));

        // And through a basis.
        BassBasisFit basis = BassComponents.Fit(BassDecomposition.FromPca(Curves(Inputs.Mixed, Outputs), new BassPcaSpec { Rank = 2 }), Predictors(Inputs.Mixed), config);
        Assert.Equal(1, basis.JointDraws.Count);
        Assert.Equal(1, basis.Predict(NewCases(Inputs.Mixed)).DrawCount);
        Assert.Equal(1, basis.Sensitivity().Count);
        Assert.Equal(Outputs, basis.TrainingFit().Locations);
    }

    private static void BassPointwiseLogDensityShape(PointwiseLogDensity density)
        => Assert.Equal((1, Cases), (density.Draws, density.Values.GetLength(1)));

    [Fact]
    public void OneCase_IsRefusedWithItsReason_AndTwoCasesFit()
    {
        Assert.Contains("no knot", Assert.Throws<ArgumentException>(() => new ContinuousPredictor("x", [1.0])).Message);
        Assert.Contains("no information", Assert.Throws<ArgumentException>(() => new CategoricalPredictor("z", [0], ["a", "b"])).Message);
        Assert.Contains("two cases", Assert.Throws<ArgumentException>(() => BassDecomposition.FromPca(new double[1, 4])).Message);

        var data = new BassData([new ContinuousPredictor("x", [0.0, 1.0])], new ScalarResponse([1.0, 2.0]));
        BassFit fit = Bass.Fit(Config(samples: 20, burn: 50), data);
        Assert.Equal(1, fit.Resolved.MaxBasis);    // the intercept counts toward the dimension bound
        Assert.All(Pooled(fit), d => Assert.InRange(d.Structure.TermIds.Length, 0, 1));
        BassPrediction p = fit.Predict(data, new PredictionCases([new ContinuousColumn("x", [0.5])]));
        Assert.Equal(1, p.Cases);
        Assert.All(p.Latent.Cast<double>(), v => Assert.True(double.IsFinite(v)));
    }

    // ── Full PCA rank through a fit ────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void FullPcaRank_FitsPredictsAndDecomposes()
    {
        double[,] y = Curves(Inputs.Continuous, 4);    // m = 4 < n, so the centred rank is 4
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 1.0 });
        Assert.Equal(4, d.Rank);
        Assert.Equal(d.Pca!.NumericalRank, d.Rank);
        BassPredictor[] predictors = Predictors(Inputs.Continuous);
        BassBasisFit fit = BassComponents.Fit(d, predictors, Config(samples: 12, burn: 60));
        BassBasisPrediction p = fit.Predict(NewCases(Inputs.Continuous));
        Assert.Equal((fit.JointDraws.Count, 3, 4), (p.DrawCount, p.Cases, p.Locations));
        Assert.All(p.Predictive!.Cast<double>(), v => Assert.True(double.IsFinite(v)));
        BasisSobolAnalysis analysis = fit.Sensitivity();
        Assert.Equal(4, analysis.Components);
        BassTrainingFit training = fit.TrainingFit();
        Assert.Equal((Cases, 4), (training.Cases, training.Locations));
    }

    // ── Zero-probability levels and truncated measures through a fit ───────────────────────────────────────────────────

    [Fact]
    public void ZeroProbabilityLevel_SensitivityMatchesEnumeration()
    {
        BassData data = Data(Inputs.Categorical, Response.Scalar);
        BassFit fit = Bass.Fit(Config(samples: 30, burn: 200), data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        var spec = new SensitivitySpec
        {
            Measures = new Dictionary<string, InputMeasure>
            {
                ["z0"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.0, ["b"] = 0.3, ["c"] = 0.7 }),
                ["z1"] = new CategoricalMeasure(new Dictionary<string, double> { ["p"] = 1.0, ["q"] = 0.0 }),
            },
        };
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data, spec).Decompositions);
        for (int d = 0; d < pooled.Count; d++)
        {
            (double V, double Mean, double[] Main) grid = Moments([Levels(0.0, 0.3, 0.7), Levels(1.0, 0.0)], n => Hand(fit, inputs, pooled[d], [], [n[0].Level, n[1].Level], 0.0));
            if (grid.V < 1e-20)
            {
                Assert.False(s.IsDefined(d));
                continue;
            }
            AgreeWithGrid($"draw {d}", s, d, grid, 1e-10);
        }
    }

    [Theory]
    [InlineData(Response.Scalar)]
    [InlineData(Response.Pca)]
    public void TruncatedMeasure_SensitivityMatchesQuadrature(Response response)
    {
        var uniform = new UniformMeasure(0.25, 0.75);
        var mixture = new MixtureMeasure([new NormalComponent(0.5, 0.2, 0.15), new StudentComponent(0.5, 3.0, 0.7, 0.1)]) { Lower = 0.1, Upper = 0.9 };
        var measures = new Dictionary<string, InputMeasure> { ["x0"] = uniform, ["x1"] = mixture };
        var spec = new SensitivitySpec { Measures = measures };

        // The truncated mixture's nodes, the midpoint rule on its density; weights before truncation renormalized.
        const int points = 400;
        double Density(double v) => 0.5 * Math.Exp(-0.5 * (v - 0.2) * (v - 0.2) / (0.15 * 0.15)) / (0.15 * Math.Sqrt(2 * Math.PI))
            + 0.5 * Math.Pow(1 + (v - 0.7) * (v - 0.7) / (0.1 * 0.1) / 3.0, -2.0) * Math.Exp(SpecialFunctions.LogGamma(2.0) - SpecialFunctions.LogGamma(1.5)) / (0.1 * Math.Sqrt(3.0 * Math.PI));
        double[] xs = [.. Enumerable.Range(0, points).Select(g => 0.1 + (g + 0.5) / points * 0.8)];
        double total = xs.Sum(Density);
        Node[] x1Nodes = [.. xs.Select(v => new Node(v, -1, Density(v) / total))];
        Node[][] grid = [Continuous(points, 0.25, 0.75), x1Nodes];

        if (response == Response.Scalar)
        {
            BassData data = Data(Inputs.Continuous, Response.Scalar);
            BassFit fit = Bass.Fit(Config(samples: 16, burn: 200), data);
            BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
            List<BassDraw> pooled = Pooled(fit);
            SobolDecomposition s = Assert.Single(fit.Sensitivity(data, spec).Decompositions);
            int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
            AgreeWithGrid("scalar", s, richest, Moments(grid, n => Hand(fit, inputs, pooled[richest], [n[0].Value, n[1].Value], [], 0.0)), 0.01);
        }
        else
        {
            BassDecomposition d = BassDecomposition.FromPca(Curves(Inputs.Continuous, Outputs), new BassPcaSpec { Rank = 2 });
            BassBasisFit fit = BassComponents.Fit(d, Predictors(Inputs.Continuous), Config(samples: 16, burn: 200));
            BassInputs inputs = BassInputs.Prepare(fit.ComponentData(0), fit.Components[0].Transform);
            List<BassDraw>[] pooled = [.. fit.Components.Select(Pooled)];
            int joint = Richest(Enumerable.Range(0, fit.JointDraws.Count), j => Enumerable.Range(0, 2).Sum(k => pooled[k][fit.JointDraws.Draw(k, j)].Structure.TermIds.Length));
            int location = d.Locations / 2;
            double At(Node[] n)
            {
                double value = 0.0;
                for (int k = 0; k < 2; k++)
                    value += d.Map.Basis(location, k) * Hand(fit.Components[k], inputs, pooled[k][fit.JointDraws.Draw(k, joint)], [n[0].Value, n[1].Value], [], 0.0);
                return d.Map.Offsets[location] + d.Map.Scales[location] * value;
            }
            AgreeWithGrid("basis", fit.Sensitivity(spec).AtLocation(location), joint, Moments(grid, At), 0.01);
        }
    }

    // A measure wider than the training range: the hinges extrapolate past [0, 1], and the moments follow them.
    [Fact]
    public void WideMeasure_ExtrapolatesTheHinges()
    {
        BassData data = Data(Inputs.Continuous, Response.Scalar);
        BassFit fit = Bass.Fit(Config(samples: 16, burn: 200), data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        var measures = new Dictionary<string, InputMeasure> { ["x0"] = new UniformMeasure(-0.5, 1.5), ["x1"] = new UniformMeasure(0.0, 1.25) };
        SobolAnalysis analysis = fit.Sensitivity(data, new SensitivitySpec { Measures = measures });
        Assert.True(analysis.Measure.Continuous(0)!.Extrapolates);
        int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
        Node[][] grid = [Continuous(300, -0.5, 1.5), Continuous(300, 0.0, 1.25)];
        AgreeWithGrid("wide", analysis.Decompositions[0], richest, Moments(grid, n => Hand(fit, inputs, pooled[richest], [n[0].Value, n[1].Value], [], 0.0)), 0.01);
    }

    // ── Higher degree and several output variables through a fit ───────────────────────────────────────────────────────

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void HigherDegree_PredictsDecomposesAndHasEffects(int degree)
    {
        BassData data = Data(Inputs.Mixed, Response.Direct);
        BassConfig config = Config(samples: 16, burn: 200);
        config = config with { Structure = config.Structure with { Degree = degree } };
        BassFit fit = Bass.Fit(config, data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        double[] t = [0.1, 0.8];
        BassPrediction p = fit.Predict(data, NewCases(Inputs.Mixed), new BassPredictOptions { Predictive = false, OutputLocations = [new ContinuousColumn("t", t)] });
        (double[][] nx, int[][] nz) = NewScaled(Inputs.Mixed);
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 3; i++)
                for (int l = 0; l < 2; l++)
                    Assert.Equal(Hand(fit, inputs, pooled[d], [nx[0][i]], [nz[0][i]], t[l]), p.Latent[d, i, l], 11);
        int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data).Decompositions);
        Node[][] grid = [Continuous(60), Levels(1.0 / 3, 1.0 / 3, 1.0 / 3), Continuous(60)];
        double Eval(Node[] n) => Hand(fit, inputs, pooled[richest], [n[0].Value], [n[1].Level], n[2].Value);
        (double V, double Mean, double[] Main) moments = Moments(grid, Eval);
        AgreeWithGrid($"degree {degree}", s, richest, moments, 0.02);
        EffectFunctions effects = fit.Effects(data, null, new EffectSpec { GridPoints = 5 }).At[0];
        Effect main = effects.Effects.Single(e => e.Inputs.SequenceEqual(new[] { 0 }));
        for (int c = 0; c < main.Grids[0].Values.Length; c++)
        {
            Node[][] given = [[new Node(main.Grids[0].Values[c], -1, 1.0)], grid[1], grid[2]];
            Assert.True(Math.Abs(main.Values[richest, c] - (Moments(given, Eval).Mean - moments.Mean)) <= 0.03 * Math.Max(Math.Abs(moments.Mean), 1.0));
        }
    }

    private static BassData TwoOutputVariables(int locations = 12)
    {
        var rng = new Xoshiro256PlusPlus(21);
        double[] t0 = [.. Enumerable.Range(0, locations).Select(l => (l % 4) / 3.0)];
        double[] t1 = [.. Enumerable.Range(0, locations).Select(l => (l * 5 % locations) / (locations - 1.0))];   // no grid
        double[] x = X0();
        var y = new double[Cases, locations];
        for (int i = 0; i < Cases; i++)
            for (int l = 0; l < locations; l++)
                y[i, l] = 1 + 2 * x[i] * t0[l] + Math.Max(t1[l] - 0.4, 0) + 0.05 * Normal.SampleStandard(rng);
        return new BassData([new ContinuousPredictor("x0", x), new ContinuousPredictor("x1", X1())],
                            new FunctionalResponse(y, [new ContinuousPredictor("t0", t0), new ContinuousPredictor("t1", t1)]));
    }

    [Fact]
    public void TwoOutputVariables_PredictAndIntegrate()
    {
        BassData data = TwoOutputVariables();
        BassFit fit = Bass.Fit(Config(samples: 16, burn: 200), data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        double[] t0 = [0.2, 0.7], t1 = [0.9, 0.3];
        BassPrediction p = fit.Predict(data, NewCases(Inputs.Continuous), new BassPredictOptions
        {
            Predictive = false,
            OutputLocations = [new ContinuousColumn("t1", t1), new ContinuousColumn("t0", t0)],
        });
        (double[][] nx, _) = NewScaled(Inputs.Continuous);
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 3; i++)
                for (int l = 0; l < 2; l++)
                    Assert.Equal(Hand(fit, inputs, pooled[d], [nx[0][i], nx[1][i]], [], 0.0, [t0[l], t1[l]]), p.Latent[d, i, l], 11);

        SobolAnalysis integrated = fit.Sensitivity(data);
        SobolDecomposition s = Assert.Single(integrated.Decompositions);
        Assert.Equal(new[] { "x0", "x1", "t0", "t1" }, s.Inputs);
        int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
        Node[][] grid = [Continuous(16), Continuous(16), Continuous(16), Continuous(16)];
        AgreeWithGrid("two output variables", s, richest, Moments(grid, n => Hand(fit, inputs, pooled[richest], [n[0].Value, n[1].Value], [], 0.0, [n[2].Value, n[3].Value])), 0.04);

        // Conditioning on a whole coordinate, every output variable named, reads as the grid over the inputs alone.
        double[] at = [0.5, 0.25];
        SobolDecomposition c = Assert.Single(fit.Sensitivity(data, new SensitivitySpec { Functional = FunctionalForm.Condition, Coordinates = [at] }).Decompositions);
        Assert.Equal(new[] { "x0", "x1" }, c.Inputs);
        AgreeWithGrid("a whole coordinate", c, richest, Moments([Continuous(200), Continuous(200)], n => Hand(fit, inputs, pooled[richest], [n[0].Value, n[1].Value], [], 0.0, at)), 0.02);
    }

    // R's func.var: indices as a function of one output variable, the others integrated. Audit gap G-2, closed.
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void ConditionOnOneOfTwoOutputVariables_IntegratesTheOther(int degree)
    {
        BassData data = TwoOutputVariables();
        BassConfig config = Config(samples: 16, burn: 200);
        BassFit fit = Bass.Fit(config with { Structure = config.Structure with { Degree = degree } }, data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        List<BassDraw> pooled = Pooled(fit);
        var spec = new SensitivitySpec { Functional = FunctionalForm.Condition, ConditionOn = ["t0"], Coordinates = [[0.5]] };
        SobolDecomposition s = Assert.Single(fit.Sensitivity(data, spec).Decompositions);
        Assert.Equal(new[] { "x0", "x1", "t1" }, s.Inputs);
        int richest = Richest(Enumerable.Range(0, pooled.Count), d => pooled[d].Structure.TermIds.Length);
        Node[][] grid = [Continuous(40), Continuous(40), Continuous(40)];
        double Eval(Node[] n) => Hand(fit, inputs, pooled[richest], [n[0].Value, n[1].Value], [], 0.0, [0.5, n[2].Value]);
        (double V, double Mean, double[] Main) moments = Moments(grid, Eval);
        AgreeWithGrid("t0 = 0.5", s, richest, moments, 0.03);

        // The effects share the split: the mean, and t1's main effect at its grid.
        EffectFunctions effects = fit.Effects(data, spec, new EffectSpec { GridPoints = 5 }).At[0];
        Assert.True(Math.Abs(effects.Mean[richest] - moments.Mean) <= 0.01 * Math.Max(Math.Abs(moments.Mean), 1.0));
        Effect t1 = effects.Effects.Single(e => e.Inputs.SequenceEqual(new[] { 2 }));
        for (int c = 0; c < t1.Grids[0].Values.Length; c++)
        {
            Node[][] given = [grid[0], grid[1], [new Node(t1.Grids[0].Values[c], -1, 1.0)]];
            double expected = Moments(given, Eval).Mean - moments.Mean;
            Assert.True(Math.Abs(t1.Values[richest, c] - expected) <= 0.01 * Math.Max(Math.Abs(moments.Mean), 1.0),
                        $"t1's effect at {t1.Grids[0].Values[c]}: {t1.Values[richest, c]} against {expected}");
        }
    }

    // Terms where the split between held and integrated factors matters: one held whole, which folds into the level; one
    // integrated whole, with no case factor, which still varies; a case factor with both; and a hinge facing its end of the
    // range, whose part R leaves unnormalized, at a coordinate past that end where it is not zero. The hand evaluation
    // divides the whole part by one constant, so a normalizer split or counted twice fails here.
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void ConditionOnOne_HandBuiltTerms_CarryTheNormalizerOnce(int degree)
    {
        BassData data = TwoOutputVariables();
        BassConfig config = Config(samples: 2, chains: 1, burn: 5);
        BassFit run = Bass.Fit(config with { Structure = config.Structure with { Degree = degree } }, data);
        BassInputs inputs = BassInputs.Prepare(data, run.Transform);
        Assert.Equal(1.0, inputs.FunctionalKnots(0)[3]);
        var table = new BassTermTable();
        int[] ids =
        [
            table.Intern(new BassTerm(functional: [new HingeFactor(0, 1, 1)])),
            table.Intern(new BassTerm(functional: [new HingeFactor(1, -1, 6)])),
            table.Intern(new BassTerm([new HingeFactor(0, 1, 2)], functional: [new HingeFactor(0, -1, 2), new HingeFactor(1, 1, 3)])),
            table.Intern(new BassTerm([new HingeFactor(1, -1, 3)], functional: [new HingeFactor(0, 1, 3), new HingeFactor(1, 1, 1)])),
        ];
        var store = new DrawStore<BassDraw>();
        store.Add(new BassDraw(new BassStructure(ids), 1.0, 1.0, [0.5, 1.5, -0.8, 2.0, 3.0], 1.0));
        BassFit fit = run with { Terms = table, Draws = [store], DrawWeights = null };

        double[] coordinates = [0.5, 1.2];
        var spec = new SensitivitySpec { Functional = FunctionalForm.Condition, ConditionOn = ["t0"], Coordinates = [[coordinates[0]], [coordinates[1]]] };
        SobolAnalysis analysis = fit.Sensitivity(data, spec);
        EffectAnalysis effects = fit.Effects(data, spec, new EffectSpec { Subsets = [] });
        Assert.Equal(new[] { false, true }, analysis.CoordinatesOutside);
        Node[][] grid = [Continuous(60), Continuous(60), Continuous(60)];
        for (int c = 0; c < 2; c++)
        {
            double at = coordinates[c];
            (double V, double Mean, double[] Main) moments = Moments(grid, n => Hand(fit, inputs, store[0], [n[0].Value, n[1].Value], [], 0.0, [at, n[2].Value]));
            AgreeWithGrid($"t0 = {at}", analysis.Decompositions[c], 0, moments, 0.01);
            Assert.True(moments.V > 0.0);
            Assert.True(Math.Abs(effects.At[c].Mean[0] - moments.Mean) <= 0.005 * Math.Max(Math.Abs(moments.Mean), 1.0), $"t0 = {at}: mean {effects.At[c].Mean[0]} against {moments.Mean}");
        }
    }

    // ── Declared controls ──────────────────────────────────────────────────────────────────────────────────────────────

    // The audit's re-trace found these controls read at their defaults only, or built by hand past the config: a declared
    // value must reach the model and bind the retained terms.
    [Fact]
    public void DeclaredControls_ReachTheModelAndTheRetainedTerms()
    {
        // w₁, w₂, h₁ and h₂, and caps below what the variables allow.
        BassData two = TwoOutputVariables();
        BassConfig declared = Config(samples: 16, burn: 200) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 1, MaxFunctionalInteraction = 1, MaxBasis = 6, OrderWeight = 2, VariableWeight = 0.5 },
            Count = new BassCountPrior(3, 7),
        };
        BassModel model = BassModel.Create(declared, two);
        Assert.Equal((1, 1), (model.Resolved.MaxInteraction, model.Resolved.MaxFunctionalInteraction));
        Assert.Equal((2.0, 0.5), (model.Proposal.OrderWeight, model.Proposal.VariableWeight));
        Assert.Equal(Math.Log(3.0 / (1 * 8)), model.Prior.LogCountRatio(0), 14);   // (M + h₁)/((M + 1)(h₂ + 1))
        Assert.Equal(Math.Log(4.0 / (2 * 8)), model.Prior.LogCountRatio(1), 14);
        BassFit capped = Bass.Fit(declared, two);
        BassTerm[] terms = [.. Enumerable.Range(0, capped.Terms.Count).Select(id => capped.Terms[id])];
        Assert.NotEmpty(terms);
        Assert.All(terms, t => Assert.True(t.Continuous.Count <= 1 && t.Functional.Count <= 1, t.ToString()));

        BassData factors = Data(Inputs.Categorical, Response.Scalar);
        BassFit oneFactor = Bass.Fit(Config(samples: 16, burn: 200) with { Structure = new BassStructureSpec { MaxCategoricalInteraction = 1, MaxBasis = 6 } }, factors);
        Assert.All(Enumerable.Range(0, oneFactor.Terms.Count), id => Assert.True(oneFactor.Terms[id].Categorical.Count <= 1));

        // A declared npart.func every functional factor falls short of: t's hinges are zero at their knot's end of the six
        // locations, so at six no term keeps a functional part, where the default's threshold of 0.6 keeps them.
        BassData direct = Data(Inputs.Continuous, Response.Direct);
        BassFit open = Bass.Fit(Config(samples: 16, burn: 200), direct);
        Assert.Contains(Enumerable.Range(0, open.Terms.Count), id => open.Terms[id].Functional.Count > 0);
        BassConfig strict = Config(samples: 16, burn: 200) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6, Support = new SupportRule(SupportKind.Count, null, Locations) },
        };
        BassFit closed = Bass.Fit(strict, direct);
        Assert.Equal(Locations, closed.Resolved.SupportMinimumFunctional);
        Assert.NotEqual(0, closed.Diagnostics.Refusals.Count(BassRefusal.Support));
        Assert.All(Enumerable.Range(0, closed.Terms.Count), id => Assert.Empty(closed.Terms[id].Functional));
    }

    // ── What the hub names and the members leave to the reader ─────────────────────────────────────────────────────────

    // R's `other`: the variance outside the reported subsets, at a location and over the locations, the higher orders' mains
    // read off the full decomposition. Audit gap G-1, closed.
    [Fact]
    public void BasisOther_IsTheVarianceLessTheReportedMains()
    {
        BassDecomposition d = Decomposition(Inputs.Continuous, Response.Pca);
        BassBasisFit fit = BassComponents.Fit(d, Predictors(Inputs.Continuous), Config());
        BasisSobolAnalysis wholeFit = fit.Sensitivity(), firstOrder = fit.Sensitivity(new SensitivitySpec { MaxOrder = 1 });
        foreach ((SobolDecomposition full, SobolDecomposition first) in new[]
                 {
                     (wholeFit.AtLocation(d.Locations / 2), firstOrder.AtLocation(d.Locations / 2)),
                     (wholeFit.Integrated, firstOrder.Integrated),
                 })
        {
            Assert.Contains(full.Subsets, u => u.Length == 2);
            Assert.All(full.OtherVariance, v => Assert.Equal(0.0, v));
            int seen = 0;
            for (int j = 0; j < first.Count; j++)
            {
                if (!first.IsDefined(j)) continue;
                double reported = Enumerable.Range(0, first.Subsets.Count).Sum(u => first.MainVariance[j, u]);
                double higher = Enumerable.Range(0, full.Subsets.Count).Where(u => full.Subsets[u].Length > 1).Sum(u => full.MainVariance[j, u]);
                Assert.Equal(higher, first.Variance[j] - reported, 1e-9 * first.Variance[j]);
                Assert.Equal(higher, first.OtherVariance[j], 1e-9 * first.Variance[j]);
                seen++;
            }
            Assert.True(seen > 0);
        }
    }

    // A direct fit needs its output variable to vary: a single location has no knot, and a scalar fit stands in.
    [Fact]
    public void OneLocation_DirectIsRefused_AScalarFitStandsIn()
    {
        Assert.Contains("no knot", Assert.Throws<ArgumentException>(() => new ContinuousPredictor("t", [0.5])).Message);
        Assert.Contains("no knot", Assert.Throws<ArgumentException>(() =>
            new FunctionalResponse(new double[3, 1], [new ContinuousPredictor("t", [0.5])])).Message);
    }
}
