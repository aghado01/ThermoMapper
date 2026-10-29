using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The residual-dependence check (stage 4, S9 and S10; F11): a hand-worked case, ties and an unordered index; the BASS check
/// against a direct computation (the posterior mean from explicit columns, equal and given weights, an output variable
/// stored out of order), which the per-term route must match; a deterministic smooth misfit that a case-only structure
/// cannot reach reading high; independent noise reading about −1/m within its sampling error; null where no order is
/// decided; Moran's I hand-worked, and on two output coordinates under independence and under a smooth misfit; each
/// declared index against its direct reading; and the refusals.
/// </summary>
public sealed class BassResidualDependenceTests
{
    private readonly ITestOutputHelper _out;
    public BassResidualDependenceTests(ITestOutputHelper output) => _out = output;

    // x₀ on four values and x₁ on n/4 (cases i), the output variable t at the given locations, y = f(x₀, x₁, t) + σε.
    private static BassData Functional(int n, double[] t, Func<double, double, double, double> f, double sigma, int seed,
                                       ContinuousPredictor? covariate = null)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int m = t.Length, rows = n / 4;
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i % 4) / 3.0;
            x1[i] = (i / 4) / (rows - 1.0);
            for (int r = 0; r < m; r++) y[i, r] = f(x0[i], x1[i], t[r]) + sigma * Normal.SampleStandard(rng);
        }
        List<ContinuousPredictor> outputs = [new ContinuousPredictor("t", t)];
        if (covariate is not null) outputs.Add(covariate);
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new FunctionalResponse(y, outputs));
    }

    private static double[] Grid(int m) => [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];

    // r₁ of a sequence by its definition.
    private static double LagOne(double[] e)
    {
        double mean = e.Average();
        double numerator = 0, denominator = 0;
        for (int k = 0; k < e.Length; k++) denominator += (e[k] - mean) * (e[k] - mean);
        for (int k = 0; k + 1 < e.Length; k++) numerator += (e[k] - mean) * (e[k + 1] - mean);
        return numerator / denominator;
    }

    // The type-7 quantile: position p(n − 1), interpolated between its neighbours.
    private static double Quantile(double[] values, double p)
    {
        double[] sorted = [.. values.OrderBy(v => v)];
        double h = p * (sorted.Length - 1);
        int lo = (int)Math.Floor(h);
        return lo + 1 < sorted.Length ? sorted[lo] + (h - lo) * (sorted[lo + 1] - sorted[lo]) : sorted[lo];
    }

    [Fact]
    public void Of_HandWorkedTiesAndOrder()
    {
        // e = 1, 2, 3, 4: ē = 2.5, Σd² = 5, Σdₖdₖ₊₁ = 0.75 − 0.25 + 0.75, so r₁ = 0.25; reversed, the same.
        Assert.Equal(0.25, ResidualDependence.LagOne(new double[,] { { 1, 2, 3, 4 } }, [0, 1, 2, 3]).PerCase[0], 15);
        Assert.Equal(0.25, ResidualDependence.LagOne(new double[,] { { 1, 2, 3, 4 } }, [3, 2, 1, 0]).PerCase[0], 15);
        // Index 0, 2, 1, 3 reads the locations as 1, 3, 2, 4: Σdₖdₖ₊₁ = −0.75 − 0.25 − 0.75, r₁ = −0.35.
        Assert.Equal(-0.35, ResidualDependence.LagOne(new double[,] { { 1, 2, 3, 4 } }, [0, 2, 1, 3]).PerCase[0], 15);
        // Ties keep the locations' order: 0, 0, 1, 1 reads 1, 2, 3, 4; and 1, 1, 0, 0 reads 3, 4, 1, 2.
        Assert.Equal(0.25, ResidualDependence.LagOne(new double[,] { { 1, 2, 3, 4 } }, [0, 0, 1, 1]).PerCase[0], 15);
        Assert.Equal(LagOne([3, 4, 1, 2]), ResidualDependence.LagOne(new double[,] { { 1, 2, 3, 4 } }, [1, 1, 0, 0]).PerCase[0], 15);

        // A case whose residual does not vary has no r₁, and the summary leaves it out.
        ResidualDependence check = ResidualDependence.LagOne(new double[,] { { 1, 2, 3, 4 }, { 5, 5, 5, 5 }, { 1, 3, 2, 4 } }, [0, 1, 2, 3], [0.0, 0.5, 1.0]);
        Assert.True(double.IsNaN(check.PerCase[1]));
        Assert.Equal((0.25 - 0.35) / 2, check.Mean, 15);
        Assert.Equal([0.0, 0.5, 1.0], check.Probabilities);
        Assert.Equal(-0.35, check.Quantiles[0], 15);
        Assert.Equal(-0.05, check.Quantiles[1], 15);
        Assert.Equal(0.25, check.Quantiles[2], 15);
        Assert.Equal(ResidualDependence.DefaultProbabilities, ResidualDependence.LagOne(new double[,] { { 1, 2 } }, [0, 1]).Probabilities);

        Assert.Empty(ResidualDependence.LagOne(new double[,] { { 1, 2 } }, [0, 1], []).Quantiles);

        ResidualDependence none = ResidualDependence.LagOne(new double[,] { { 2, 2, 2 } }, [0, 1, 2]);
        Assert.True(double.IsNaN(none.Mean));
        Assert.All(none.Quantiles, q => Assert.True(double.IsNaN(q)));
    }

    [Fact]
    public void Bass_MatchesADirectComputation()
    {
        // t stored out of order, so the check must sort the locations; a smooth signal with noise, fitted by a chain.
        double[] t = [.. new[] { 4, 0, 7, 1, 6, 2, 5, 3 }.Select(r => r / 7.0)];
        BassData data = Functional(12, t, (x0, x1, s) => 1 + 2 * x0 + 1.5 * s * x1 + Math.Sin(3 * s), 0.2, seed: 7);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), data, seed: 61, kept: 15);
        BassDraw[] pooled = [.. draws.SelectMany(c => c)];
        var rng = new Xoshiro256PlusPlus(62);
        double[] weights = [.. pooled.Select(_ => rng.NextDouble())];
        weights[3] = 0.0;
        double[] y = model.Data.Response.Values.ToArray();
        int n = 12, m = t.Length;
        int[] order = [.. Enumerable.Range(0, m).OrderBy(r => t[r])];
        double[] probabilities = [0.0, 0.1, 0.5, 0.9, 1.0];

        foreach (double[]? w in new[] { null, weights })
        {
            ResidualDependence check = BassResidualDependence.Of(model, table, draws, weights: w, probabilities: probabilities)!;
            // The posterior mean from each draw's explicit uncentred columns.
            var mean = new double[y.Length];
            double total = 0;
            for (int s = 0; s < pooled.Length; s++)
            {
                double ws = w?[s] ?? 1.0;
                double[][] columns = [.. pooled[s].Structure.TermIds.Select(id => BassPointwiseTests.Column(model, table[id]))];
                for (int k = 0; k < y.Length; k++)
                {
                    double f = pooled[s].Coefficients[0];
                    for (int q = 0; q < columns.Length; q++) f += pooled[s].Coefficients[q + 1] * columns[q][k];
                    mean[k] += ws * f;
                }
                total += ws;
            }
            var expected = new double[n];
            bool orderMatters = false;
            for (int i = 0; i < n; i++)
            {
                double[] e = [.. Enumerable.Range(0, m).Select(r => y[i * m + r] - mean[i * m + r] / total)];
                expected[i] = LagOne([.. order.Select(r => e[r])]);
                orderMatters |= Math.Abs(LagOne(e) - expected[i]) > 0.05;
                Assert.True(Math.Abs(expected[i] - check.PerCase[i]) <= 1e-12, $"case {i}: {check.PerCase[i]:R} against {expected[i]:R}");
            }
            Assert.True(orderMatters, "the stored order would read the same values");
            Assert.Equal(expected.Average(), check.Mean, 12);
            for (int q = 0; q < probabilities.Length; q++) Assert.Equal(Quantile(expected, probabilities[q]), check.Quantiles[q], 12);
            _out.WriteLine($"{(w is null ? "equal" : "given")} weights: mean r₁ {check.Mean:F3}, quantiles {string.Join(", ", check.Quantiles.Select(v => v.ToString("F3")))}");
        }
        // The weights move the posterior mean, and so the values.
        Assert.NotEqual(BassResidualDependence.Of(model, table, draws)!.PerCase, BassResidualDependence.Of(model, table, draws, weights: weights)!.PerCase);
    }

    [Fact]
    public void SmoothMisfit_ReadsHigh()
    {
        // A deterministic simulator, y = 1 + 2x₀ + sin(2πt) on 20 locations, fitted by a structure with no functional
        // factor: the fit is flat in t in every case, so each case's residual is the sine about its own level and reads
        // the sine's lag-one value.
        double[] t = Grid(20);
        BassData data = Functional(12, t, (x0, _, s) => 1 + 2 * x0 + Math.Sin(2 * Math.PI * s), 0.0, seed: 1);
        var model = BassModel.Create(BassPointwiseTests.Config(maxBasis: 4), data);
        BassTerm[] caseOnly = [new BassTerm([new HingeFactor(0, 1, 0)]), new BassTerm([new HingeFactor(1, -1, 2)])];
        var moves = new BassStructureMoves(model, caseOnly, Math.Log(0.01));
        var table = new BassTermTable();
        var rng = new Xoshiro256PlusPlus(2);
        List<List<BassDraw>> draws = [[.. Enumerable.Range(0, 20).Select(_ => moves.Compose(moves.Start, rng, table))]];
        ResidualDependence check = BassResidualDependence.Of(model, table, draws)!;
        double sine = LagOne([.. t.Select(s => Math.Sin(2 * Math.PI * s))]);
        _out.WriteLine($"case-only structure: mean r₁ {check.Mean:F4}, the sine's {sine:F4}");
        Assert.True(sine > 0.85);
        Assert.All(check.PerCase, r => Assert.True(Math.Abs(r - sine) <= 1e-9, $"{r:R} against {sine:R}"));

        // The same simulator fitted by chains capped at two terms, which cannot reach the sine either.
        (BassModel capped, BassTermTable ct, List<List<BassDraw>> cd) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(maxBasis: 2), data, seed: 63, kept: 20);
        ResidualDependence chains = BassResidualDependence.Of(capped, ct, cd)!;
        _out.WriteLine($"capped chains: mean r₁ {chains.Mean:F3}, 5% quantile {chains.Quantiles[0]:F3}");
        Assert.True(chains.Quantiles[0] > 0.5, $"5% quantile {chains.Quantiles[0]:F3}");
    }

    [Fact]
    public void IndependentNoise_ReadsAboutMinusOneOverM()
    {
        // y = 1 + 2x₀ + 1.5·t·x₁ + 0.3ε, within the model's span, on 16 locations: under independence r₁ is about −1/m with
        // a standard deviation about 1/√m per case, so the mean over 32 cases sits within its sampling error of −1/m.
        double[] t = Grid(16);
        BassData data = Functional(32, t, (x0, x1, s) => 1 + 2 * x0 + 1.5 * s * x1, 0.3, seed: 9);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), data, seed: 64, burn: 300, kept: 40, thin: 4);
        ResidualDependence check = BassResidualDependence.Of(model, table, draws)!;
        double m = t.Length, sd = Math.Sqrt(check.PerCase.Select(r => (r - check.Mean) * (r - check.Mean)).Sum() / (check.PerCase.Length - 1));
        double error = sd / Math.Sqrt(check.PerCase.Length);
        _out.WriteLine($"independent noise: mean r₁ {check.Mean:F3} ± {error:F3} against {-1 / m:F3}; sd {sd:F3} against {1 / Math.Sqrt(m):F3}");
        Assert.True(Math.Abs(check.Mean + 1 / m) <= 4 * error, $"mean r₁ {check.Mean:F3} ± {error:F3}");
        Assert.True(sd < 2 / Math.Sqrt(m), $"sd {sd:F3}");
        Assert.All(check.PerCase, r => Assert.True(double.IsFinite(r)));
    }

    [Fact]
    public void NullWhereNoOrderIsDecided()
    {
        // A scalar response, and a direct response with two output variables (t and a covariate known at each location).
        (BassModel scalar, BassTermTable st, List<List<BassDraw>> sd) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), BassPointwiseTests.Data(BassPointwiseTests.Kind.Scalar), seed: 65, kept: 3);
        Assert.Null(BassResidualDependence.Of(scalar, st, sd));

        double[] t = Grid(6);
        var elevation = new ContinuousPredictor("elevation", [0.3, 0.1, 0.4, 0.1, 0.5, 0.9]);
        BassData two = Functional(8, t, (x0, _, s) => x0 + s, 0.1, seed: 3, covariate: elevation);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = BassPointwiseTests.Run(BassPointwiseTests.Config(), two, seed: 66, kept: 3);
        Assert.Equal(2, model.Inputs.Outputs!.Count);
        Assert.Null(BassResidualDependence.Of(model, table, draws));
        Assert.Equal(ResidualStatistic.LagOne, BassResidualDependence.Of(model, table, draws, new ResidualIndex.AlongVariable(1))!.Statistic);
    }

    [Fact]
    public void Refusals()
    {
        Assert.Throws<ArgumentNullException>("residuals", () => ResidualDependence.LagOne(null!, [0, 1]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.LagOne(new double[0, 2], [0, 1]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.LagOne(new double[2, 1], [0]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.LagOne(new double[,] { { 1, double.NaN } }, [0, 1]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.LagOne(new double[,] { { 1, double.PositiveInfinity } }, [0, 1]));
        Assert.Throws<ArgumentException>("index", () => ResidualDependence.LagOne(new double[2, 3], [0, 1]));
        Assert.Throws<ArgumentException>("index", () => ResidualDependence.LagOne(new double[2, 3], [0, double.NaN, 1]));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => ResidualDependence.LagOne(new double[2, 3], [0, 1, 2], [0.5, 1.5]));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => ResidualDependence.LagOne(new double[2, 3], [0, 1, 2], [double.NaN]));

        BassData data = Functional(8, Grid(5), (x0, _, s) => x0 + s, 0.1, seed: 4);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = BassPointwiseTests.Run(BassPointwiseTests.Config(), data, seed: 67, kept: 3);
        int count = draws.Sum(c => c.Count);
        Assert.Throws<ArgumentNullException>("model", () => BassResidualDependence.Of(null!, table, draws));
        Assert.Throws<ArgumentNullException>("terms", () => BassResidualDependence.Of(model, null!, draws));
        Assert.Throws<ArgumentNullException>("draws", () => BassResidualDependence.Of(model, table, null!));
        Assert.Throws<ArgumentException>("draws", () => BassResidualDependence.Of(model, table, []));
        Assert.Throws<ArgumentException>("draws", () => BassResidualDependence.Of(model, table, [draws[0], []]));
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, weights: new double[count - 1]));
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, weights: new double[count]));
        double[] negative = [.. Enumerable.Repeat(1.0, count)];
        negative[0] = -1;
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, weights: negative));
        negative[0] = double.NaN;
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, weights: negative));
        negative[0] = double.PositiveInfinity;
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, weights: negative));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => BassResidualDependence.Of(model, table, draws, probabilities: [-0.1]));

        // The arguments are checked where the check is null too.
        (BassModel scalar, BassTermTable st, List<List<BassDraw>> sd) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), BassPointwiseTests.Data(BassPointwiseTests.Kind.Scalar), seed: 68, kept: 2);
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(scalar, st, sd, weights: [1.0]));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => BassResidualDependence.Of(scalar, st, sd, probabilities: [2.0]));
    }

    [Fact]
    public void MoransI_HandWorked()
    {
        // On a line at 0, 1, 2, 3 with K = 1 (ties to the lower location): 0→1, 1→0, 2→1, 3→2. e = 1, 2, 3, 4 has
        // d = −1.5, −0.5, 0.5, 1.5 and Σd² = 5, so Σ dₖd_N(k) = 0.75 + 0.75 − 0.25 + 0.75 = 2 and I = 2/5.
        double[,] rising = { { 1, 2, 3, 4 } };
        ResidualDependence line = ResidualDependence.MoransI(rising, [[0, 1, 2, 3]], neighbours: 1);
        Assert.Equal(ResidualStatistic.MoransI, line.Statistic);
        Assert.Equal(0.4, line.PerCase[0], 15);
        Assert.Equal(-1.0 / 3, line.Expected, 15);

        // A 2 × 2 grid with K = 2 reads each location's rook neighbours: a checkerboard is I = −1 and a ramp across the
        // rows is I = 0, its neighbours' deviations cancelling.
        double[][] grid = [[0, 1, 0, 1], [0, 0, 1, 1]];
        Assert.Equal(-1.0, ResidualDependence.MoransI(new double[,] { { 1, 0, 0, 1 } }, grid, neighbours: 2).PerCase[0], 15);
        Assert.Equal(0.0, ResidualDependence.MoransI(new double[,] { { 1, 2, 3, 4 } }, grid, neighbours: 2).PerCase[0], 15);
        // Two rows that are constant along x: each location's neighbours are its row-mate and its column-mate, I = 0.
        Assert.Equal(0.0, ResidualDependence.MoransI(new double[,] { { 1, 1, 2, 2 } }, grid, neighbours: 2).PerCase[0], 15);

        // A constant case has no value, and the summary leaves it out.
        ResidualDependence mixed = ResidualDependence.MoransI(new double[,] { { 1, 0, 0, 1 }, { 3, 3, 3, 3 } }, grid, 2, [0.5]);
        Assert.True(double.IsNaN(mixed.PerCase[1]));
        Assert.Equal(-1.0, mixed.Mean, 15);
        Assert.Equal(-1.0, mixed.Quantiles[0], 15);

        // r₁'s expectation is −1/m, to first order.
        Assert.Equal(-0.25, ResidualDependence.LagOne(rising, [0, 1, 2, 3]).Expected, 15);
    }

    [Fact]
    public void MoransI_OnTwoCoordinates_ReadsIndependenceAndMisfit()
    {
        // A direct response over a 5 × 5 grid of output locations (s, u): y = 1 + 2x₀ + 1.5·s·x₁ + 0.3ε lies in the model's
        // span, so Moran's I over the rook neighbours sits about −1/(m − 1); adding sin(πs)·sin(πu), which no structure
        // capped at two terms reaches, reads high.
        double[] s = [.. Enumerable.Range(0, 25).Select(r => (r % 5) / 4.0)];
        double[] u = [.. Enumerable.Range(0, 25).Select(r => (r / 5) / 4.0)];
        var second = new ContinuousPredictor("u", u);
        var near = new ResidualIndex.NearestNeighbours([0, 1]);

        BassData independent = Functional(32, s, (x0, x1, a) => 1 + 2 * x0 + 1.5 * a * x1, 0.3, seed: 19, covariate: second);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), independent, seed: 69, burn: 300, kept: 40, thin: 4);
        ResidualDependence check = BassResidualDependence.Of(model, table, draws, near)!;
        double sd = Math.Sqrt(check.PerCase.Select(v => (v - check.Mean) * (v - check.Mean)).Sum() / (check.PerCase.Length - 1));
        double error = sd / Math.Sqrt(check.PerCase.Length);
        _out.WriteLine($"independent noise: mean I {check.Mean:F3} ± {error:F3} against {check.Expected:F3}");
        Assert.Equal(-1.0 / 24, check.Expected, 15);
        Assert.True(Math.Abs(check.Mean - check.Expected) <= 4 * error, $"mean I {check.Mean:F3} ± {error:F3}");

        BassData misfit = Functional(32, s, (x0, _, a) => 1 + 2 * x0 + 0.0 * a, 0.05, seed: 20, covariate: second);
        double[] y = misfit.Response.Values.ToArray();
        var values = new double[32, 25];
        for (int i = 0; i < 32; i++)
            for (int r = 0; r < 25; r++) values[i, r] = y[i * 25 + r] + Math.Sin(Math.PI * s[r]) * Math.Sin(Math.PI * u[r]);
        BassData smooth = new(misfit.Predictors, new FunctionalResponse(values, [new ContinuousPredictor("s", s), second]));
        (BassModel capped, BassTermTable ct, List<List<BassDraw>> cd) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(maxBasis: 2), smooth, seed: 70, kept: 20);
        ResidualDependence high = BassResidualDependence.Of(capped, ct, cd, near)!;
        _out.WriteLine($"smooth misfit: mean I {high.Mean:F3}, 5% quantile {high.Quantiles[0]:F3}");
        Assert.True(high.Quantiles[0] > 0.3, $"5% quantile {high.Quantiles[0]:F3}");
    }

    [Fact]
    public void DeclaredIndices_MatchTheirDirectReadings()
    {
        // Two output variables, t and a covariate; each declared index reads the posterior-mean residuals as the
        // record's own constructors read them.
        double[] t = Grid(6);
        var elevation = new ContinuousPredictor("elevation", [0.3, 0.1, 0.4, 0.2, 0.5, 0.9]);
        BassData two = Functional(8, t, (x0, x1, s) => x0 + s * x1, 0.1, seed: 5, covariate: elevation);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = BassPointwiseTests.Run(BassPointwiseTests.Config(), two, seed: 71, kept: 6);
        OutputInputs outputs = model.Inputs.Outputs!;
        double[,] residuals = Residuals(model, table, draws);

        ResidualDependence alongT = BassResidualDependence.Of(model, table, draws, new ResidualIndex.AlongVariable(0))!;
        ResidualDependence alongElevation = BassResidualDependence.Of(model, table, draws, new ResidualIndex.AlongVariable(1))!;
        Assert.Equal(ResidualDependence.LagOne(residuals, outputs.Variable(0)).PerCase, alongT.PerCase, Close);
        Assert.Equal(ResidualDependence.LagOne(residuals, outputs.Variable(1)).PerCase, alongElevation.PerCase, Close);
        Assert.NotEqual(alongT.PerCase, alongElevation.PerCase);

        int[] order = [5, 3, 0, 2, 4, 1];
        double[] keys = new double[6];
        for (int k = 0; k < 6; k++) keys[order[k]] = k;
        ResidualDependence alongOrder = BassResidualDependence.Of(model, table, draws, new ResidualIndex.AlongOrder(order))!;
        Assert.Equal(ResidualDependence.LagOne(residuals, keys).PerCase, alongOrder.PerCase, Close);

        ResidualDependence near = BassResidualDependence.Of(model, table, draws, new ResidualIndex.NearestNeighbours([1, 0], 2))!;
        Assert.Equal(ResidualDependence.MoransI(residuals, [outputs.Variable(1).ToArray(), outputs.Variable(0).ToArray()], 2).PerCase, near.PerCase, Close);

        // The default on one output variable is the declared index along it.
        BassData one = Functional(8, t, (x0, x1, s) => x0 + s * x1, 0.1, seed: 5);
        (BassModel m1, BassTermTable t1, List<List<BassDraw>> d1) = BassPointwiseTests.Run(BassPointwiseTests.Config(), one, seed: 72, kept: 4);
        Assert.Equal(BassResidualDependence.Of(m1, t1, d1)!.PerCase, BassResidualDependence.Of(m1, t1, d1, new ResidualIndex.AlongVariable(0))!.PerCase);
    }

    [Fact]
    public void Fit_ReadsADeclaredIndex()
    {
        double[] t = Grid(6);
        var elevation = new ContinuousPredictor("elevation", [0.3, 0.1, 0.4, 0.2, 0.5, 0.9]);
        BassData two = Functional(12, t, (x0, x1, s) => x0 + s * x1, 0.1, seed: 6, covariate: elevation);
        BassConfig config = BassRunTests.Config(chains: 2, samples: 20, maxBasis: 3);
        BassFit fit = Bass.Fit(config, two);
        Assert.Null(fit.ResidualDependence);
        var model = BassModel.Create(config, two);
        var index = new ResidualIndex.NearestNeighbours([0, 1], 3);
        Assert.Equal(BassResidualDependence.Of(model, fit.Terms, fit.Draws, index)!.PerCase, fit.ResidualDependenceFor(two, index).PerCase);
        Assert.Throws<ArgumentNullException>("index", () => fit.ResidualDependenceFor(two, null!));
        BassData other = Functional(12, t, (x0, x1, s) => x0 + s * x1, 0.1, seed: 7, covariate: elevation);
        Assert.Throws<ArgumentException>("data", () => fit.ResidualDependenceFor(other, index));
    }

    [Fact]
    public void DeclaredIndex_Refusals()
    {
        Assert.Throws<ArgumentNullException>("residuals", () => ResidualDependence.MoransI(null!, [[0, 1]]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.MoransI(new double[1, 1], [[0]]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.MoransI(new double[,] { { 1, double.NaN } }, [[0, 1]]));
        Assert.Throws<ArgumentNullException>("coordinates", () => ResidualDependence.MoransI(new double[1, 2], null!));
        Assert.Throws<ArgumentException>("coordinates", () => ResidualDependence.MoransI(new double[1, 2], []));
        Assert.Throws<ArgumentException>("coordinates", () => ResidualDependence.MoransI(new double[1, 2], [[0, 1, 2]]));
        Assert.Throws<ArgumentException>("coordinates", () => ResidualDependence.MoransI(new double[1, 2], [[0, double.NaN]]));
        Assert.Throws<ArgumentOutOfRangeException>("neighbours", () => ResidualDependence.MoransI(new double[1, 3], [[0, 1, 2]], 0));
        Assert.Throws<ArgumentOutOfRangeException>("neighbours", () => ResidualDependence.MoransI(new double[1, 3], [[0, 1, 2]], 3));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => ResidualDependence.MoransI(new double[1, 3], [[0, 1, 2]], 1, [1.5]));

        double[] t = Grid(5);
        var elevation = new ContinuousPredictor("elevation", [0.3, 0.1, 0.4, 0.2, 0.5]);
        BassData two = Functional(8, t, (x0, _, s) => x0 + s, 0.1, seed: 8, covariate: elevation);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = BassPointwiseTests.Run(BassPointwiseTests.Config(), two, seed: 73, kept: 3);
        foreach (ResidualIndex bad in new ResidualIndex[]
                 {
                     new ResidualIndex.AlongVariable(2), new ResidualIndex.AlongVariable(-1),
                     new ResidualIndex.AlongOrder([0, 1, 2, 3]), new ResidualIndex.AlongOrder([0, 1, 2, 3, 3]), new ResidualIndex.AlongOrder([0, 1, 2, 3, 5]),
                     new ResidualIndex.AlongOrder(null!),
                     new ResidualIndex.NearestNeighbours([]), new ResidualIndex.NearestNeighbours(null!), new ResidualIndex.NearestNeighbours([0, 0]),
                     new ResidualIndex.NearestNeighbours([0, 2]), new ResidualIndex.NearestNeighbours([0, 1], 0), new ResidualIndex.NearestNeighbours([0, 1], 5),
                 })
            Assert.Throws<ArgumentException>("index", () => BassResidualDependence.Of(model, table, draws, bad));

        // A scalar response has no output index to declare.
        (BassModel scalar, BassTermTable st, List<List<BassDraw>> sd) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), BassPointwiseTests.Data(BassPointwiseTests.Kind.Scalar), seed: 74, kept: 2);
        Assert.Throws<ArgumentException>("index", () => BassResidualDependence.Of(scalar, st, sd, new ResidualIndex.AlongVariable(0)));
    }

    private static readonly IEqualityComparer<double> Close = EqualityComparer<double>.Create((a, b) => Math.Abs(a - b) <= 1e-12 || (double.IsNaN(a) && double.IsNaN(b)), _ => 0);

    // The posterior-mean residuals of the draws, each draw's fit from its explicit uncentred columns.
    private static double[,] Residuals(BassModel model, BassTermTable table, List<List<BassDraw>> draws)
    {
        BassDraw[] pooled = [.. draws.SelectMany(c => c)];
        double[] y = model.Data.Response.Values.ToArray();
        int n = model.Data.Response.Cases, m = model.Data.Response.Locations;
        var mean = new double[y.Length];
        foreach (BassDraw draw in pooled)
        {
            double[][] columns = [.. draw.Structure.TermIds.Select(id => BassPointwiseTests.Column(model, table[id]))];
            for (int k = 0; k < y.Length; k++)
            {
                double f = draw.Coefficients[0];
                for (int q = 0; q < columns.Length; q++) f += draw.Coefficients[q + 1] * columns[q][k];
                mean[k] += f / pooled.Length;
            }
        }
        var residuals = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int r = 0; r < m; r++) residuals[i, r] = y[i * m + r] - mean[i * m + r];
        return residuals;
    }
}
