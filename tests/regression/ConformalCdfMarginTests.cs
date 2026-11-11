using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class ConformalCdfMarginTests(ITestOutputHelper output)
{
    private static readonly ConformalOptions Kept = new() { MarginBounds = MarginBounds.Kept };

    internal static PointwiseLogDensity Training(double[] x, double[] y, double[] beta)
    {
        var ll = new double[beta.Length, x.Length]; var cdf = new double[beta.Length, x.Length];
        for (int d = 0; d < beta.Length; d++) for (int i = 0; i < x.Length; i++)
        {
            ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
            cdf[d, i] = Normal.StandardLogCdf(y[i] - beta[d] * x[i]);
        }
        return new(ll, UnitPartition.EachObservation(x.Length), ConformalBayesTests.Theta) { LogCdf = cdf };
    }

    internal static CandidateLogDensity Candidates(double[] beta, double x, double[] grid)
    {
        var ll = new double[beta.Length, grid.Length]; var cdf = new double[beta.Length, grid.Length];
        for (int d = 0; d < beta.Length; d++) for (int g = 0; g < grid.Length; g++)
        {
            ll[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * x);
            cdf[d, g] = Normal.StandardLogCdf(grid[g] - beta[d] * x);
        }
        return new(ll, [grid], ConformalBayesTests.Theta) { LogCdf = cdf };
    }

    [Fact]
    public void CountsMatchLiteralTwoTailUnionAndPitAtEveryLevel()
    {
        var rng = new Xoshiro256PlusPlus(2026111110);
        for (int m = 3; m <= 7; m++) for (int fixture = 0; fixture < 40; fixture++)
        {
            double[] s = Enumerable.Range(0, m).Select(_ => Math.Floor(8 * rng.NextDouble()) / 8).ToArray();
            double[] e = Enumerable.Range(0, m).Select(_ => Math.Floor(4 * rng.NextDouble()) / 16).ToArray();
            bool[] flags = Enumerable.Range(0, m).Select(i => i < m - 1 && rng.NextDouble() < .3).ToArray();
            var r = CdfMarginCounts.Read(s.Select(v => new CdfScore(v, 1 - v)).ToArray(), e, e, false, flags, MarginBounds.Kept, InflationBounds.PerPosition);
            for (int k = 0; k <= m; k++)
            {
                int band = 0, accepted = 0, lower = 0, upper = 0;
                for (int u = 0; u < m; u++)
                {
                    bool Low(bool negate)
                    {
                        double[] a = negate ? s.Select(v => -v).ToArray() : s;
                        int inflated = Enumerable.Range(0, m).Count(j => a[j] - e[j] <= a[u] + e[u]);
                        int shift = Enumerable.Range(0, m).Count(j => j == u || a[j] + 2 * e[j] + 2 * e[u] <= a[u]);
                        return inflated > k && shift <= k;
                    }
                    if (Low(false) || Low(true)) band++;
                    bool inf = Enumerable.Range(0, m).Count(j => s[j] - e[j] <= s[u] + e[u]) > k
                        && Enumerable.Range(0, m).Count(j => -s[j] - e[j] <= -s[u] + e[u]) > k;
                    if (flags[u] && inf) accepted++;
                    if (!flags[u] && Enumerable.Range(0, m).Count(j => s[j] <= s[u] + 2 * (e[j] + e[u])) <= k) lower++;
                    if (flags[u] || Enumerable.Range(0, m).Count(j => s[j] + 2 * (e[j] + e[u]) < s[u]) < k) upper++;
                }
                Assert.Equal(band, r.Fragile![k]); Assert.Equal(accepted, r.DegenerateAccepted![k]);
                Assert.Equal(lower, r.PitLowerCounts![k]); Assert.Equal(upper, r.PitUpperCounts![k]);
            }
        }
    }

    [Fact]
    public void DoubledPitCountsDominateTheExactShiftCountsUnderBoundedError()
    {
        var rng = new Xoshiro256PlusPlus(2026111111);
        for (int fixture = 0; fixture < 100; fixture++)
        {
            const int m = 9;
            double[] s = Enumerable.Range(0, m).Select(_ => .2 + .6 * rng.NextDouble()).ToArray();
            double[] e = Enumerable.Range(0, m).Select(_ => .1 * rng.NextDouble()).ToArray();
            double[] estimate = s.Select((v, i) => v + (2 * rng.NextDouble() - 1) * e[i]).ToArray();
            var r = CdfMarginCounts.Read(estimate.Select(v => new CdfScore(v, 1 - v)).ToArray(), e, e, false, null, MarginBounds.Carried, InflationBounds.PerPosition);
            for (int k = 0; k <= m; k++)
            {
                int upper = Enumerable.Range(0, m).Count(u => Enumerable.Range(0, m).Count(j => s[j] + e[j] + e[u] < s[u]) < k);
                int lower = Enumerable.Range(0, m).Count(u => Enumerable.Range(0, m).Count(j => s[j] <= s[u] + e[j] + e[u]) <= k);
                Assert.True(r.PitUpperCounts![k] >= upper); Assert.True(r.PitLowerCounts![k] <= lower);
            }
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ProbabilityErrorsMatchCenteredInfluencesAndCarriedTurns(bool weighted)
    {
        double[] x = [-2, -1, 0, 1, 2, 3], y = [-1, -1, 0, 1, 1, 4];
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 128, new Xoshiro256PlusPlus(43));
        var training = Training(x, y, beta); var candidates = Candidates(beta, 2, [-3, 0, 3]);
        if (weighted)
        {
            double[] weights = Enumerable.Range(0, beta.Length).Select(i => 1.0 + i % 7).ToArray();
            training = training with { Weights = weights }; candidates = candidates with { Weights = weights };
        }
        var d = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Distributions[0];
        double[] baseWeights = weighted ? Enumerable.Range(0, beta.Length).Select(i => 1.0 + i % 7).ToArray() : Enumerable.Repeat(1.0, beta.Length).ToArray();
        double[] W(int column, bool candidate)
        {
            double[] w = Enumerable.Range(0, beta.Length).Select(i => baseWeights[i] * Math.Exp(candidate ? candidates.Values[i, column] : training.Values[i, column])).ToArray();
            double sum = w.Sum(); return w.Select(v => v / sum).ToArray();
        }
        double Error(double[,] field, int column, double[] w)
        {
            double[] f = Enumerable.Range(0, beta.Length).Select(i => Math.Exp(field[i, column])).ToArray();
            double mean = f.Select((v, i) => w[i] * v).Sum();
            return d.Margin!.Multiple * Math.Sqrt(f.Select((v, i) => w[i] * w[i] * (v - mean) * (v - mean)).Sum());
        }
        for (int g = 0; g < 3; g++)
        {
            var r = d.Margin!.Readings[g]; var w = W(g, true);
            for (int j = 0; j <= x.Length; j++)
            {
                double[,] f = j == x.Length ? candidates.LogCdf! : training.LogCdf!; int column = j == x.Length ? g : j;
                double held = Error(f, column, w), carried = held;
                for (int i = 0; i < x.Length; i++) carried = Math.Max(carried, Error(f, column, W(i, false)));
                Assert.Equal(held, r.PositionErrors[j], 11); Assert.Equal(carried, r.CarriedErrors![j], 11);
            }
        }
    }

    [Fact]
    public void TinyUpperTailsRemainDistinctAtZeroErrorAndTiesKeepTheirMass()
    {
        CdfScore[] scores = [new(1, 4e-200), new(1, 2e-200), new(1, 1e-200), new(1, 3e-200)];
        var r = CdfMarginCounts.Read(scores, new double[4], new double[4], false, null, MarginBounds.Carried, InflationBounds.PerPosition);
        Assert.Equal(2, r.InflatedLowerRank); Assert.Equal(3, r.InflatedUpperRank);
        for (int k = 0; k <= 4; k++) { Assert.Equal(k, r.PitLowerCounts![k]); Assert.Equal(k, r.PitUpperCounts![k]); }
        var tied = CdfMarginCounts.Read(Enumerable.Repeat(new CdfScore(.5, .5), 4).ToArray(), new double[4], new double[4], false, null, MarginBounds.Carried, InflationBounds.PerPosition);
        Assert.Equal(4, tied.InflatedLowerRank); Assert.Equal(4, tied.InflatedUpperRank);
        Assert.Equal(0, tied.PitLowerCounts![2]); Assert.Equal(4, tied.PitUpperCounts![2]);

        var ll = new double[16, 3]; var fc = new double[16, 3]; var lc = new double[16, 1]; var cc = new double[16, 1];
        for (int d = 0; d < 16; d++)
        {
            for (int i = 0; i < 3; i++) { ll[d, i] = -1; fc[d, i] = -scores[i].Upper; }
            lc[d, 0] = -1; cc[d, 0] = -3e-200;
        }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(3), ConformalBayesTests.Theta) { LogCdf = fc };
        var candidates = new CandidateLogDensity(lc, [[0.0]], ConformalBayesTests.Theta) { LogCdf = cc };
        var distribution = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Distributions[0];
        Assert.Equal(1, distribution.Below[0]); Assert.Equal(0, distribution.Tied[0]);
        Assert.Equal(2, distribution.Margin!.InflatedLowerRanks[0]); Assert.Equal(3, distribution.Margin.InflatedUpperRanks[0]);
        Assert.All(distribution.Margin.Readings[0].PositionErrors, e => Assert.Equal(0, e));
        for (int k = 0; k <= 4; k++) { Assert.Equal(k, distribution.Margin.PitLowerCounts![k]); Assert.Equal(k, distribution.Margin.PitUpperCounts![k]); }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void CdfAndDensityUseTheSameWeightFlagsForChainsAndWeightedRows(bool weighted)
    {
        double[] beta = Enumerable.Range(0, 128).Select(i => .4 * Math.Sin(.2 * i)).ToArray();
        var training = Training([-2, -1, 0, 1, 2, 3, 4], [-1, -1, 0, 1, 1, 4, 2], beta);
        var candidates = Candidates(beta, 2, [-20, 0, 20]);
        if (weighted)
        {
            double[] weights = Enumerable.Range(0, 128).Select(i => i % 7 == 0 ? 0.0 : 1.0 + i % 5).ToArray();
            training = training with { Weights = weights }; candidates = candidates with { Weights = weights };
        }
        else { training = training with { ChainLengths = new[] { 64, 64 } }; candidates = candidates with { ChainLengths = new[] { 64, 64 } }; }
        foreach (double floor in new[] { 10.0, 100.0, 400.0 })
        {
            var options = Kept with { Alpha = .4, DegenerateEss = floor };
            var density = ConformalBayes.FullStep(training, candidates, options).Sets[0];
            var cdf = ConformalBayes.Distribution(training, candidates, options).Distributions[0];
            for (int g = 0; g < 3; g++)
            {
                Assert.Equal(density.Margin!.Flagged![g], cdf.Margin!.Readings[g].Flagged);
                Assert.Equal(density.Margin.CandidateFlagged![g], cdf.Margin.Readings[g].CandidateFlagged);
            }
        }
    }

    [Theory]
    [InlineData(1e-16)]
    [InlineData(double.Epsilon)]
    public void TinyErrorProbabilityKeepsBothConformalPathsFinite(double nu)
    {
        double[] beta = new double[16];
        var training = Training([1, 1, 1], [-1, 0, 1], beta); var candidates = Candidates(beta, 1, [0]);
        var options = new ConformalOptions { Alpha = .5, ErrorLevel = nu, MarginBounds = MarginBounds.Carried };
        var density = ConformalBayes.FullStep(training, candidates, options).Sets[0];
        var cdf = ConformalBayes.Distribution(training, candidates, options).Distributions[0];
        Assert.True(double.IsFinite(density.Margin!.Multiple)); Assert.Equal(density.Margin.Multiple, cdf.Margin!.Multiple);
        Assert.All(density.Margin.CandidateErrors, e => Assert.Equal(0, e));
        Assert.All(cdf.Margin.Readings[0].PositionErrors, e => Assert.Equal(0, e));
        Assert.True(double.IsFinite(density.CoverageUnderError)); Assert.True(double.IsFinite(cdf.CoverageUnderError(.5)));
    }

    [Fact]
    public void CarriedTurnBlocksPreserveScalarMaxima()
    {
        const int n = 259, s = 8;
        double[] x = Enumerable.Range(0, n).Select(i => .01 * (i - 100)).ToArray();
        double[] y = x.Select(v => Math.Sin(v)).ToArray(), beta = Enumerable.Range(0, s).Select(i => .05 * i).ToArray();
        var training = Training(x, y, beta); var candidates = Candidates(beta, 1, [.2]);
        var result = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Distributions[0];
        foreach (int j in new[] { 0, 128, 258 })
        {
            double largest = 0;
            for (int i = 0; i < n; i++)
            {
                double[] w = Enumerable.Range(0, s).Select(d => Math.Exp(training.Values[d, i])).ToArray(); double total = w.Sum();
                for (int d = 0; d < s; d++) w[d] /= total;
                double[] f = Enumerable.Range(0, s).Select(d => Math.Exp(training.LogCdf![d, j])).ToArray();
                double mean = f.Select((v, d) => w[d] * v).Sum();
                largest = Math.Max(largest, result.Margin!.Multiple * Math.Sqrt(f.Select((v, d) => w[d] * w[d] * (v - mean) * (v - mean)).Sum()));
            }
            Assert.Equal(largest, result.Margin!.UnitCarriedErrors[j], 11);
        }
    }

    [Fact]
    public void SubnormalWeightedTurnsNormalizeCombinedLogsBeforeComputingCarriedErrors()
    {
        double[,] ll = { { -745, 0, 0 }, { 0, 0, 0 } };
        double[,] fc = { { Math.Log(.25), Math.Log(.25), Math.Log(.25) }, { Math.Log(.75), Math.Log(.75), Math.Log(.75) } };
        double[] weights = [1, double.Epsilon];
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(3), ConformalBayesTests.Theta) { LogCdf = fc, Weights = weights };
        var candidates = new CandidateLogDensity(new double[2, 1], [[0.0]], ConformalBayesTests.Theta)
            { LogCdf = new double[,] { { Math.Log(.5) }, { Math.Log(.5) } }, Weights = weights };
        var d = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Distributions[0];
        double a = Math.Exp(-745 - Math.Log(double.Epsilon)), w0 = a / (1 + a), w1 = 1 / (1 + a);
        double mean = w0 * .25 + w1 * .75;
        double expected = d.Margin!.Multiple * Math.Sqrt(w0 * w0 * (.25 - mean) * (.25 - mean) + w1 * w1 * (.75 - mean) * (.75 - mean));
        Assert.Equal(expected, d.Margin.UnitCarriedErrors[1], 12);
        Assert.InRange(w0, .36, .37);
    }

    [Fact]
    public void CenteredInfluencesRetainTheDominantDrawsTinyDeviation()
    {
        var field = new CdfMomentField(new double[,] { { Math.Log(.25) }, { Math.Log(.75) } }, [0, 1], [0]);
        var result = field.Measure([1, 1e-200], 1, [1]);
        Assert.InRange(result.Errors[0] / 1e-200, .707, .708);
    }

    [Fact]
    public void BothTinyTailsKeepNonzeroStandardErrorsBelowTheSquaredDoubleRange()
    {
        foreach (bool upper in new[] { false, true })
        {
            var cdf = new double[4, 1];
            for (int d = 0; d < 4; d++) cdf[d, 0] = upper ? -(d + 1) * 1e-200 : Math.Log((d + 1) * 1e-200);
            var field = new CdfMomentField(cdf, [0, 1, 2, 3], [0]);
            var reading = field.Measure([.25, .25, .25, .25], 1, [1]);
            Assert.True(reading.Errors[0] > 0); Assert.True(double.IsFinite(reading.Errors[0]));
            Assert.InRange(reading.Errors[0] / 1e-200, .55, .57);
        }
    }

    [Fact]
    public void MarginOptionsAndUnsupportedCarriedTurnsAreRefused()
    {
        var training = Training([1, 1, 1], [-1, 0, 1], new double[16]);
        var candidates = Candidates(new double[16], 1, [0]);
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training, candidates, Kept with { Score = FullStepScore.LeaveOneOut }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.Distribution(training, candidates, Kept with { ErrorLevel = 0 }));
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginRefinement = new((t, y) => candidates, .1) }));
        var unsupported = (double[,])training.Values.Clone();
        for (int d = 0; d < 16; d++) unsupported[d, 0] = double.NegativeInfinity;
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training with { Values = unsupported }, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }));
    }

    [Fact]
    public void OptionsKeepRanksAndFlagsAndCanChargeEveryTurn()
    {
        double[] beta = Enumerable.Range(0, 128).Select(i => .1 * Math.Sin(i)).ToArray();
        var training = Training(Enumerable.Repeat(1.0, 9).ToArray(), new double[9], beta);
        var candidates = Candidates(beta, 1, [-3, 0, 3]);
        var old = ConformalBayes.Distribution(training, candidates).Distributions[0];
        var d = ConformalBayes.Distribution(training, candidates, Kept with { DegenerateEss = 10_000 }).Distributions[0];
        Assert.Equal(old.Below, d.Below); Assert.Equal(old.Tied, d.Tied); Assert.Equal(old.Fall, d.Fall); Assert.Equal(old.ScoreFall, d.ScoreFall);
        Assert.Null(old.Margin); Assert.Equal(10, d.Margin!.FlaggedCount);
        Assert.All(d.Margin.Readings, r => { Assert.True(r.CandidateFlagged); Assert.Null(r.Fragile); });
        Assert.Equal(-.42, d.CoverageUnderError(.4), 12);
        Assert.Equal(-.41, (d with { Inflated = true }).CoverageUnderError(.4), 12);
        Assert.Equal(new ConformalPitBounds(0, 1), d.Margin.Pit(5));
        Assert.Equal(1, d.CoverageUnderError(.1));
        var inflationOnly = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Inflate = true }).Distributions[0];
        Assert.Null(inflationOnly.Margin!.Counts); Assert.True(double.IsNaN((inflationOnly with { Inflated = false }).CoverageUnderError(.4)));
        for (int g = 0; g < 3; g++) Assert.True(!d.AcceptsPlain(g, .4) || d.AcceptsInflated(g, .4));
    }

    [Fact]
    public void RefinementReusesTheScorerAndChecksTheSeamAtExtremeFiniteEnds()
    {
        var ll = new double[128, 9]; var fc = new double[128, 9]; var lc = new double[128, 2]; var cc = new double[128, 2];
        for (int d = 0; d < 128; d++)
        {
            for (int i = 0; i < 9; i++) { ll[d, i] = -1; fc[d, i] = Math.Log(.5); }
            lc[d, 0] = d == 0 ? 0 : -100; lc[d, 1] = -1; cc[d, 0] = cc[d, 1] = Math.Log(.5);
        }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(9), ConformalBayesTests.Theta) { LogCdf = fc };
        var candidates = new CandidateLogDensity(lc, [[-1e308, 1e308]], ConformalBayesTests.Theta) { LogCdf = cc };
        var plain = ConformalBayes.Distribution(training, candidates, Kept).Distributions[0];
        int calls = 0;
        var refined = ConformalBayes.Distribution(training, candidates, Kept with { MarginRefinement = new((t, y) =>
        {
            calls++; Assert.Equal(0, y); var a = new double[128, 1]; var b = new double[128, 1];
            for (int d = 0; d < 128; d++) { a[d, 0] = -1; b[d, 0] = Math.Log(.5); }
            return new(a, [[y]], ConformalBayesTests.Theta) { LogCdf = b };
        }, 1e308) }).Distributions[0];
        Assert.Equal(1, calls); Assert.Equal(1, refined.Margin!.RefinementProbes);
        Assert.Equal(plain.Below, refined.Below); Assert.Equal(plain.Tied, refined.Tied);
        Assert.True(refined.Margin.FlaggedCount >= plain.Margin!.FlaggedCount);
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training, candidates, Kept with { MarginRefinement = new((t, y) => new(lc, [[-1.0, 1.0]], ConformalBayesTests.Theta), 1e308) }));
    }

    [Theory]
    [InlineData(MarginBounds.Uniform)]
    [InlineData(MarginBounds.Carried)]
    [InlineData(MarginBounds.Kept)]
    public void PitInterpolatesAndBothInflationModesContainThePlainSet(MarginBounds mode)
    {
        double[] x = [-2, -1, 0, 1, 2, 3], y = [-1, -1, 0, 1, 1, 4];
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 256, new Xoshiro256PlusPlus(72));
        var training = Training(x, y, beta); var candidates = Candidates(beta, 2, [-5, 0, 5]);
        foreach (var inflation in new[] { InflationBounds.PerPosition, InflationBounds.Uniform })
        {
            var d = ConformalBayes.Distribution(training, candidates, Kept with { MarginBounds = mode, InflationBounds = inflation, DegenerateEss = 10 }).Distributions[0];
            for (int g = 0; g < 3; g++) Assert.True(!d.AcceptsPlain(g, .6) || d.AcceptsInflated(g, .6));
            for (int k = 0; k < d.Units; k++)
            {
                var a = d.Margin!.Pit(k); var b = d.Margin.Pit(k + 1); var mid = d.Margin.PitAt((k + .5) / d.Units);
                Assert.Equal((a.Lower + b.Lower) / 2, mid.Lower, 12); Assert.Equal((a.Upper + b.Upper) / 2, mid.Upper, 12);
                Assert.True(a.Lower <= a.Upper);
            }
        }
    }

    [Theory]
    [InlineData(64)]
    [InlineData(512)]
    public void ExactDrawCoverageAndPitReading(int draws)
    {
        if (!BassQualification.Enabled) return;
        const int n = 19, repetitions = 1000;
        int[] seeds = SeedTree.Derive(2026111112, repetitions);
        foreach (double floor in new[] { 10.0, 100.0, 400.0 })
        {
            int hits = 0, inflated = 0; double diagnostic = 0, inflatedDiagnostic = 0, departure = 0;
            var bins = new int[10];
            foreach (int seed in seeds)
            {
                var rng = new Xoshiro256PlusPlus(seed);
                double[] x = Enumerable.Range(0, n).Select(_ => 4 * rng.NextDouble() - 2).ToArray();
                double[] y = x.Select(v => ConformalBayesTests.Truth(v, .5, rng)).ToArray();
                double xt = 4 * rng.NextDouble() - 2, yt = ConformalBayesTests.Truth(xt, .5, rng);
                double[] beta = ConformalBayesTests.ExactDraws(x, y, draws, rng);
                var d = ConformalBayes.Distribution(Training(x, y, beta), Candidates(beta, xt, [yt]), Kept with { DegenerateEss = floor }).Distributions[0];
                if (d.AcceptsPlain(0, .1)) hits++; if (d.AcceptsInflated(0, .1)) inflated++;
                diagnostic += d.CoverageUnderError(.1); inflatedDiagnostic += (d with { Inflated = true }).CoverageUnderError(.1);
                departure += d.Margin!.PitDeparture;
                bins[Math.Min(9, (int)(10 * d.At(0, rng.NextDouble())))]++;
            }
            output.WriteLine($"Exact Gaussian draws S{draws}, floor {floor}, {repetitions} repetitions: coverage {hits / (double)repetitions:F4}, inflated {inflated / (double)repetitions:F4}, nominal .9 SE {Math.Sqrt(.09 / repetitions):F4}; singleton diagnostics {diagnostic / repetitions:F4}/{inflatedDiagnostic / repetitions:F4}, mean PIT departure diagnostic {departure / repetitions:F4}, PIT bins {string.Join(',', bins)}. Singleton response grids do not certify response-uniform coverage.");
        }
    }

    [Theory]
    [InlineData("gaussian")]
    [InlineData("t")]
    [InlineData("wrong")]
    public void BassCdfKeptBesideCarried(string kind)
    {
        if (!BassQualification.Enabled) return;
        const int n = 99, fits = 5, tests = 30, widths = 4;
        int[] seeds = SeedTree.Derive(2026111113, fits);
        double[] floors = [10, 100, 400];
        var hits = new int[3]; var inflated = new int[3]; var diagnostic = new double[3]; var old = new double[3];
        var width = new double[3]; var infWidth = new double[3]; var pit = new double[3]; var bins = new int[10]; int nonmonotone = 0;
        double[] grid = Enumerable.Range(0, 161).Select(g => -6 + .1 * g).ToArray();
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            var (x0, x1, y) = BassCandidateTests.Draw(n, rng); var (t0, t1, yt) = BassCandidateTests.Draw(tests, rng);
            BassPredictor[] Predictors(double[] a, double[] b) => kind == "wrong" ? [new ContinuousPredictor("x1", b)] : [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];
            var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
            var fit = Bass.Fit(BassStudentTTests.Config(kind == "t" ? BassNoise.StudentT(3, 1.0 / 3) : null, burnIn: 1000, samples: 500), data);
            var training = fit.PointwiseLogDensity(data, cdf: true);
            var candidates = fit.CandidateLogDensity(data, PredictionCases.Of(Predictors([.. t0, .. t0[..widths]], [.. t1, .. t1[..widths]])),
                [.. yt.Select(v => new[] { v }), .. Enumerable.Repeat(grid, widths)], cdf: true);
            var carried = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried });
            for (int f = 0; f < 3; f++)
            {
                var result = ConformalBayes.Distribution(training, candidates, Kept with { DegenerateEss = floors[f] });
                for (int t = 0; t < tests; t++)
                {
                    var d = result.Distributions[t]; if (d.AcceptsPlain(0, .1)) hits[f]++; if (d.AcceptsInflated(0, .1)) inflated[f]++;
                    if (f == 0) bins[Math.Min(9, (int)(10 * d.At(0, rng.NextDouble())))]++;
                    Assert.Equal(carried.Distributions[t].Below, d.Below); Assert.Equal(carried.Distributions[t].Tied, d.Tied);
                }
                for (int t = tests; t < tests + widths; t++)
                {
                    var d = result.Distributions[t]; diagnostic[f] += d.CoverageUnderError(.1); old[f] += carried.Distributions[t].CoverageUnderError(.1);
                    width[f] += d.Width(.1); infWidth[f] += (d with { Inflated = true }).Width(.1); pit[f] += d.Margin!.PitDeparture;
                    if (f == 0 && !d.Monotone) nonmonotone++;
                }
            }
        }
        for (int f = 0; f < 3; f++) output.WriteLine($"BASS CDF {kind}, floor {floors[f]}, {fits} fits/{fits * tests} truths/{fits * widths} fixed grids: actual coverage {hits[f] / (double)(fits * tests):F4}/{inflated[f] / (double)(fits * tests):F4}, widths {width[f] / (fits * widths):F3}/{infWidth[f] / (fits * widths):F3}; kept/carried diagnostics {diagnostic[f] / (fits * widths):F4}/{old[f] / (fits * widths):F4}, PIT departure diagnostic {pit[f] / (fits * widths):F4}; nonmonotone grids {nonmonotone}, PIT bins {string.Join(',', bins)}. Conditional proxy/grid diagnostics, not guarantees.");
    }
}
