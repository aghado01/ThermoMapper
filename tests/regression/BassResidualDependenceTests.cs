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
/// The residual-dependence check (stage 4, S9; F11): a hand-worked case, ties and an unordered index; the BASS check
/// against a direct computation (the posterior mean from explicit columns, equal and given weights, an output variable
/// stored out of order); a deterministic smooth misfit that a case-only structure cannot reach reading high; independent
/// noise reading about −1/m within its sampling error; null where no order is decided; and the refusals.
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
        Assert.Equal(0.25, ResidualDependence.Of(new double[,] { { 1, 2, 3, 4 } }, [0, 1, 2, 3]).LagOne[0], 15);
        Assert.Equal(0.25, ResidualDependence.Of(new double[,] { { 1, 2, 3, 4 } }, [3, 2, 1, 0]).LagOne[0], 15);
        // Index 0, 2, 1, 3 reads the locations as 1, 3, 2, 4: Σdₖdₖ₊₁ = −0.75 − 0.25 − 0.75, r₁ = −0.35.
        Assert.Equal(-0.35, ResidualDependence.Of(new double[,] { { 1, 2, 3, 4 } }, [0, 2, 1, 3]).LagOne[0], 15);
        // Ties keep the locations' order: 0, 0, 1, 1 reads 1, 2, 3, 4; and 1, 1, 0, 0 reads 3, 4, 1, 2.
        Assert.Equal(0.25, ResidualDependence.Of(new double[,] { { 1, 2, 3, 4 } }, [0, 0, 1, 1]).LagOne[0], 15);
        Assert.Equal(LagOne([3, 4, 1, 2]), ResidualDependence.Of(new double[,] { { 1, 2, 3, 4 } }, [1, 1, 0, 0]).LagOne[0], 15);

        // A case whose residual does not vary has no r₁, and the summary leaves it out.
        ResidualDependence check = ResidualDependence.Of(new double[,] { { 1, 2, 3, 4 }, { 5, 5, 5, 5 }, { 1, 3, 2, 4 } }, [0, 1, 2, 3], [0.0, 0.5, 1.0]);
        Assert.True(double.IsNaN(check.LagOne[1]));
        Assert.Equal((0.25 - 0.35) / 2, check.Mean, 15);
        Assert.Equal([0.0, 0.5, 1.0], check.Probabilities);
        Assert.Equal(-0.35, check.Quantiles[0], 15);
        Assert.Equal(-0.05, check.Quantiles[1], 15);
        Assert.Equal(0.25, check.Quantiles[2], 15);
        Assert.Equal(ResidualDependence.DefaultProbabilities, ResidualDependence.Of(new double[,] { { 1, 2 } }, [0, 1]).Probabilities);

        ResidualDependence none = ResidualDependence.Of(new double[,] { { 2, 2, 2 } }, [0, 1, 2]);
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
            ResidualDependence check = BassResidualDependence.Of(model, table, draws, w, probabilities)!;
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
                Assert.True(Math.Abs(expected[i] - check.LagOne[i]) <= 1e-12, $"case {i}: {check.LagOne[i]:R} against {expected[i]:R}");
            }
            Assert.True(orderMatters, "the stored order would read the same values");
            Assert.Equal(expected.Average(), check.Mean, 12);
            for (int q = 0; q < probabilities.Length; q++) Assert.Equal(Quantile(expected, probabilities[q]), check.Quantiles[q], 12);
            _out.WriteLine($"{(w is null ? "equal" : "given")} weights: mean r₁ {check.Mean:F3}, quantiles {string.Join(", ", check.Quantiles.Select(v => v.ToString("F3")))}");
        }
        // The weights move the posterior mean, and so the values.
        Assert.NotEqual(BassResidualDependence.Of(model, table, draws)!.LagOne, BassResidualDependence.Of(model, table, draws, weights)!.LagOne);
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
        Assert.All(check.LagOne, r => Assert.True(Math.Abs(r - sine) <= 1e-9, $"{r:R} against {sine:R}"));

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
        double m = t.Length, sd = Math.Sqrt(check.LagOne.Select(r => (r - check.Mean) * (r - check.Mean)).Sum() / (check.LagOne.Length - 1));
        double error = sd / Math.Sqrt(check.LagOne.Length);
        _out.WriteLine($"independent noise: mean r₁ {check.Mean:F3} ± {error:F3} against {-1 / m:F3}; sd {sd:F3} against {1 / Math.Sqrt(m):F3}");
        Assert.True(Math.Abs(check.Mean + 1 / m) <= 4 * error, $"mean r₁ {check.Mean:F3} ± {error:F3}");
        Assert.True(sd < 2 / Math.Sqrt(m), $"sd {sd:F3}");
        Assert.All(check.LagOne, r => Assert.True(double.IsFinite(r)));
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
    }

    [Fact]
    public void Refusals()
    {
        Assert.Throws<ArgumentNullException>("residuals", () => ResidualDependence.Of(null!, [0, 1]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.Of(new double[0, 2], [0, 1]));
        Assert.Throws<ArgumentException>("residuals", () => ResidualDependence.Of(new double[2, 1], [0]));
        Assert.Throws<ArgumentException>("index", () => ResidualDependence.Of(new double[2, 3], [0, 1]));
        Assert.Throws<ArgumentException>("index", () => ResidualDependence.Of(new double[2, 3], [0, double.NaN, 1]));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => ResidualDependence.Of(new double[2, 3], [0, 1, 2], [0.5, 1.5]));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => ResidualDependence.Of(new double[2, 3], [0, 1, 2], [double.NaN]));

        BassData data = Functional(8, Grid(5), (x0, _, s) => x0 + s, 0.1, seed: 4);
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = BassPointwiseTests.Run(BassPointwiseTests.Config(), data, seed: 67, kept: 3);
        int count = draws.Sum(c => c.Count);
        Assert.Throws<ArgumentNullException>("model", () => BassResidualDependence.Of(null!, table, draws));
        Assert.Throws<ArgumentNullException>("terms", () => BassResidualDependence.Of(model, null!, draws));
        Assert.Throws<ArgumentNullException>("draws", () => BassResidualDependence.Of(model, table, null!));
        Assert.Throws<ArgumentException>("draws", () => BassResidualDependence.Of(model, table, []));
        Assert.Throws<ArgumentException>("draws", () => BassResidualDependence.Of(model, table, [draws[0], []]));
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, new double[count - 1]));
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, new double[count]));
        double[] negative = [.. Enumerable.Repeat(1.0, count)];
        negative[0] = -1;
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, negative));
        negative[0] = double.NaN;
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, negative));
        negative[0] = double.PositiveInfinity;
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(model, table, draws, negative));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => BassResidualDependence.Of(model, table, draws, probabilities: [-0.1]));

        // The arguments are checked where the check is null too.
        (BassModel scalar, BassTermTable st, List<List<BassDraw>> sd) =
            BassPointwiseTests.Run(BassPointwiseTests.Config(), BassPointwiseTests.Data(BassPointwiseTests.Kind.Scalar), seed: 68, kept: 2);
        Assert.Throws<ArgumentException>("weights", () => BassResidualDependence.Of(scalar, st, sd, [1.0]));
        Assert.Throws<ArgumentOutOfRangeException>("probabilities", () => BassResidualDependence.Of(scalar, st, sd, probabilities: [2.0]));
    }
}
