using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The tournament correction (Bhattacharyya, Zhang and Barber, arXiv 2605.29200; E16, the conformal layer's slice 13). The rule is the
/// cross-fitted count at k = ⌊αm⌋ (<c>tournament_count_iff</c>, <c>AcceptsPairwiseBy</c>); the pair fit from leave-one-out draws cannot tell
/// which member is the candidate (<c>tournFit_comp</c>); coverage over every ordering holds at 1 − (2k − 1)/m with an order-reading pair fit
/// that is not symmetric (<c>tournament_orderReading</c>), and fails with one that marks the candidate; the inflated set covers its stable
/// bound where the pair scores lean against the candidate by up to ε (<c>coverage_tournament_stable</c>). With the qualification switch
/// (<c>BASS_QUALIFICATION</c> = 1): coverage with chains stopped before mixing, beside the full step from the same chains, and the cost row.
/// </summary>
public sealed class TournamentTests(ITestOutputHelper output)
{
    [Fact]
    public void TheCountRule_IsTheDirectCount()
    {
        // Random scores rounded to tenths, so ties occur: the candidate is accepted iff its wins, games where its score plus δ falls below the
        // unit's, plus k fall below m.
        var rng = new Xoshiro256PlusPlus(5);
        const int n = 13, columns = 40;
        var units = new double[n, columns];
        var candidates = new double[n, columns];
        for (int i = 0; i < n; i++)
            for (int c = 0; c < columns; c++)
            {
                units[i, c] = Math.Round(rng.NextDouble() * 2.0, 1);
                candidates[i, c] = Math.Round(rng.NextDouble() * 2.0, 1);
            }
        double[] grid = [.. Enumerable.Range(0, columns).Select(c => 0.1 * c)];
        foreach (double alpha in new[] { 0.1, 0.2, 0.35 })
            foreach (double delta in new[] { 0.0, 0.3 })
            {
                ConformalSet set = CrossConformal.Tournament(units, candidates, [grid], alpha, delta).Sets[0];
                int m = n + 1, k = (int)Math.Floor(alpha * m);
                Assert.Equal(k, set.Level);
                Assert.Equal(1.0 - (2.0 * k - 1) / m, set.Coverage, 12);
                for (int c = 0; c < columns; c++)
                {
                    int wins = Enumerable.Range(0, n).Count(i => candidates[i, c] + delta < units[i, c]);
                    Assert.Equal(wins + k < m, set.AcceptsPlain(c));
                    Assert.Equal(set.AcceptsPlain(c), set.Accepts(c));
                }
            }
    }

    [Fact]
    public void ThePairFit_CannotTellWhichMemberIsTheCandidate()
    {
        // Swapping the pair swaps the scores bit for bit, each meets a direct sum over the draws, and the density form reads the scores' sets.
        var rng = new Xoshiro256PlusPlus(17);
        const int s = 300;
        double[] lw = [.. Enumerable.Range(0, s).Select(_ => Math.Log(0.2 + rng.NextDouble()))];
        double total = Math.Log(lw.Sum(Math.Exp));
        lw = [.. lw.Select(v => v - total)];
        double[] h = [.. Enumerable.Range(0, s).Select(_ => -3.0 * rng.NextDouble() - 1.0)];
        double[] c = [.. Enumerable.Range(0, s).Select(d => d % 17 == 0 ? double.NegativeInfinity : -2.0 * rng.NextDouble())];
        double unit = CrossConformal.PairScore(lw, h, c), candidate = CrossConformal.PairScore(lw, c, h);
        double Direct(double[] scored, double[] other)
        {
            double num = 0.0, den = 0.0;
            for (int d = 0; d < s; d++)
            {
                double w = Math.Exp(lw[d] + scored[d] + other[d]);
                den += w;
                num += w * Math.Exp(scored[d]);
            }
            return Math.Log(num / den);
        }
        Assert.Equal(Direct(h, c), unit, 12);
        Assert.Equal(Direct(c, h), candidate, 12);
        // The same pair handed over the other way round.
        Assert.Equal(unit, CrossConformal.PairScore(lw, h, c));
        double[] none = [.. Enumerable.Repeat(double.NegativeInfinity, s)];
        // No draw gives the pair positive density: each member's own predictive under the fit to the rest.
        Assert.Equal(Math.Log(Enumerable.Range(0, s).Sum(d => Math.Exp(lw[d] + h[d]))), CrossConformal.PairScore(lw, h, none), 12);
        Assert.Equal(double.NegativeInfinity, CrossConformal.PairScore(lw, none, h));

        // The density form against the scores' form, over four fits and a grid of five values.
        const int fits = 4;
        double[] grid = [-1.0, -0.5, 0.0, 0.5, 1.0];
        var held = new CandidateLogDensity[fits];
        var tests = new CandidateLogDensity[fits];
        var units = new double[fits, grid.Length];
        var candidates = new double[fits, grid.Length];
        for (int i = 0; i < fits; i++)
        {
            double[] mu = [.. Enumerable.Range(0, 50).Select(_ => rng.NextDouble() - 0.5)];
            double y = rng.NextDouble() - 0.5;
            var hv = new double[50, 1];
            var cv = new double[50, grid.Length];
            for (int d = 0; d < 50; d++)
            {
                hv[d, 0] = -0.5 * (y - mu[d]) * (y - mu[d]) / 0.1;
                for (int g = 0; g < grid.Length; g++) cv[d, g] = -0.5 * (grid[g] - mu[d]) * (grid[g] - mu[d]) / 0.1;
            }
            held[i] = new CandidateLogDensity(hv, [[y]], []);
            tests[i] = new CandidateLogDensity(cv, [grid], []);
            double[] zero = new double[50];
            double[] hcol = [.. Enumerable.Range(0, 50).Select(d => hv[d, 0])];
            for (int g = 0; g < grid.Length; g++)
            {
                double[] ccol = [.. Enumerable.Range(0, 50).Select(d => cv[d, g])];
                units[i, g] = CrossConformal.PairScore(zero, hcol, ccol);
                candidates[i, g] = CrossConformal.PairScore(zero, ccol, hcol);
            }
        }
        ConformalSet a = CrossConformal.Tournament(fits, i => (held[i], tests[i]), alpha: 0.4).Sets[0];
        ConformalSet b = CrossConformal.Tournament(units, candidates, [grid], alpha: 0.4).Sets[0];
        Assert.Equal(b.Below, a.Below);
        Assert.Equal(b.Tied, a.Tied);
        Assert.Throws<ArgumentException>(() => CrossConformal.Tournament(fits, i => (tests[i], tests[i]), alpha: 0.4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrossConformal.Tournament(units, candidates, [grid], 0.4, inflation: -0.1));
    }

    // Coverage over every ordering of each multiset of m values from {0, 1, 2}: the candidate at the last place, each game's pair fit read
    // from the rest in order with the left-out place closed up, the pair handed over as values (or, marked, with roles).
    private static (double Least, double Bound) EveryOrdering(int m, double alpha, bool marked)
    {
        int n = m - 1, k = (int)Math.Floor(alpha * m);
        double bound = 1.0 - (2.0 * k - 1) / m, least = double.PositiveInfinity;
        IEnumerable<int[]> Multisets(int length, int from)
        {
            if (length == 0) { yield return []; yield break; }
            for (int v = from; v < 3; v++)
                foreach (int[] rest in Multisets(length - 1, v)) yield return [v, .. rest];
        }
        IEnumerable<int[]> Permutations(int[] items)
        {
            if (items.Length <= 1) { yield return items; yield break; }
            for (int i = 0; i < items.Length; i++)
                foreach (int[] rest in Permutations([.. items[..i], .. items[(i + 1)..]])) yield return [items[i], .. rest];
        }
        // An order-reading fit to the rest, not symmetric: the rest's values weighted by their places, 1, 2, 3, ….
        double Fit(IReadOnlyList<int> rest)
        {
            double num = 0.0, den = 0.0;
            for (int j = 0; j < rest.Count; j++) { num += (j + 1) * rest[j]; den += j + 1; }
            return den > 0.0 ? num / den : 0.0;
        }
        foreach (int[] multiset in Multisets(m, 0))
        {
            int covered = 0, orderings = 0;
            foreach (int[] order in Permutations(multiset))
            {
                orderings++;
                int test = order[n];
                var units = new double[n, 1];
                var candidates = new double[n, 1];
                for (int i = 0; i < n; i++)
                {
                    int[] rest = [.. order[..i], .. order[(i + 1)..n]];
                    double fit = Fit(rest);
                    int sum = order[i] + test;   // the pair as a multiset of values
                    double Score(int a) => -Math.Abs(a - (0.6 * fit + 0.2 * sum));
                    units[i, 0] = Score(order[i]);
                    candidates[i, 0] = marked ? Score(test) - 10.0 : Score(test);
                }
                if (CrossConformal.Tournament(units, candidates, [[test + 0.0]], alpha).Sets[0].Accepts(0)) covered++;
            }
            least = Math.Min(least, covered / (double)orderings);
        }
        return (least, bound);
    }

    [Theory]
    [InlineData(5, 0.2)]
    [InlineData(5, 0.4)]
    [InlineData(6, 0.34)]
    public void Coverage_OverEveryOrdering_WithAnOrderReadingFit(int m, double alpha)
    {
        (double least, double bound) = EveryOrdering(m, alpha, marked: false);
        output.WriteLine($"m {m}, α {alpha}: least coverage over the multisets {least:F4} against 1 − (2k − 1)/m = {bound:F4}");
        Assert.True(least >= bound - 1e-12, $"coverage {least} against {bound}");
    }

    [Fact]
    public void AFitThatMarksTheCandidate_Fails()
    {
        (double least, double bound) = EveryOrdering(5, 0.2, marked: true);
        Assert.Equal(0.0, least);
        Assert.True(bound > 0.5);
    }

    [Fact]
    public void TheStableBound_OptimizesTheMarkovThreshold()
    {
        // The ledger's figure: at m 20, k 2 and ν 0.01 the excess over k/m is 0.227 against their 0.300.
        Assert.Equal(1.0 - 0.1 - 0.1 - 2.0 * 19 * 0.01 / 3.0, CrossConformal.TournamentStableCoverage(20, 2, 0.01), 12);
        Assert.Equal(0.9, CrossConformal.TournamentStableCoverage(20, 2, 0.0), 12);
        Assert.True(CrossConformal.TournamentStableCoverage(1000, 100, 1e-4) >= 1.0 - 0.1 - 3.0 * 0.01);
        Assert.Throws<ArgumentOutOfRangeException>(() => CrossConformal.TournamentStableCoverage(20, 2, 1.5));
    }

    [Fact]
    public void TheInflatedSet_CoversItsStableBound_WherePairScoresLeanAgainstTheCandidate()
    {
        // m 20, k 2. Each game's scores lean against the candidate by up to ε about the reference R(z) = −|z|, and with chance ν the
        // candidate's side strays far below: no symmetry, so Theorem 1 does not apply, and the plain set undercovers. The 2ε-inflated set
        // covers at least its stable bound.
        const int m = 20, n = m - 1, replications = 40000;
        const double eps = 0.05, nu = 0.01, alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(23);
        var units = new double[n, 1];
        var candidates = new double[n, 1];
        int plain = 0, inflated = 0;
        for (int r = 0; r < replications; r++)
        {
            double test = Gaussian(rng);
            for (int i = 0; i < n; i++)
            {
                double z = Gaussian(rng);
                units[i, 0] = -Math.Abs(z) + eps * rng.NextDouble();
                candidates[i, 0] = rng.NextDouble() < nu ? -Math.Abs(test) - 5.0 : -Math.Abs(test) - eps * rng.NextDouble();
            }
            if (CrossConformal.Tournament(units, candidates, [[test]], alpha).Sets[0].Accepts(0)) plain++;
            if (CrossConformal.Tournament(units, candidates, [[test]], alpha, inflation: 2.0 * eps).Sets[0].Accepts(0)) inflated++;
        }
        double p = plain / (double)replications, q = inflated / (double)replications;
        double bound = CrossConformal.TournamentStableCoverage(m, 2, nu);
        double se = Math.Sqrt(q * (1.0 - q) / replications);
        output.WriteLine($"plain {p:F4}, inflated {q:F4} (se {se:F4}); the stable bound {bound:F4}, 1 − 2α {1.0 - (2.0 * 2 - 1) / m:F4}");
        Assert.True(q >= bound - 4.0 * se);
        Assert.True(q > p);
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // ---- Behind the qualification switch ----

    // Bayesian linear regression y = β₀ + β₁x + N(0, σ²), σ² = 0.25, β ~ N(0, 16 I), by random-walk Metropolis from β = (4, −4) with steps of
    // 0.1, from one seed, the log likelihood summed over the units in their order: a fit that reads its units in order. Returns the kept draws.
    private static double[][] Metropolis(IReadOnlyList<double> x, IReadOnlyList<double> y, int steps, int keep, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double LogTarget(double b0, double b1)
        {
            double l = -(b0 * b0 + b1 * b1) / 32.0;
            for (int i = 0; i < x.Count; i++) { double e = y[i] - b0 - b1 * x[i]; l -= e * e / 0.5; }
            return l;
        }
        double c0 = 4.0, c1 = -4.0, current = LogTarget(c0, c1);
        var kept = new List<double[]>();
        for (int t = 0; t < steps; t++)
        {
            double p0 = c0 + 0.1 * Gaussian(rng), p1 = c1 + 0.1 * Gaussian(rng), proposed = LogTarget(p0, p1);
            if (Math.Log(1.0 - rng.NextDouble()) < proposed - current) { c0 = p0; c1 = p1; current = proposed; }
            if (t >= steps - keep) kept.Add([c0, c1]);
        }
        return [.. kept];
    }

    private static double LogNormal(double y, double mean) => -0.5 * Math.Log(2.0 * Math.PI * 0.25) - (y - mean) * (y - mean) / 0.5;

    private static CandidateLogDensity Seam(double[][] draws, double x, double[] grid)
    {
        var v = new double[draws.Length, grid.Length];
        for (int d = 0; d < draws.Length; d++)
            for (int g = 0; g < grid.Length; g++) v[d, g] = LogNormal(grid[g], draws[d][0] + draws[d][1] * x);
        return new CandidateLogDensity(v, [grid], []);
    }

    [Theory]
    [InlineData("stopped before mixing", 30, 15)]
    [InlineData("mixed", 3000, 1000)]
    public void Coverage_WithMetropolisChains_BesideTheFullStep(string chains, int steps, int keep)
    {
        // n 39 (m 40, k 4, 1 − (2k − 1)/m = 0.825), 1,000 replications of one held-out case at its truth from the model, and 200 more on a grid
        // for the width; every chain from the replication's seed, the full step's on all 39 units and the tournament's on each D₋ᵢ.
        if (!BassQualification.Enabled) return;
        const int n = 39, replications = 1000, widths = 200;
        double[] grid = [.. Enumerable.Range(0, 301).Select(g => -2.0 + 0.02 * g)];
        int fullCovered = 0, tournamentCovered = 0;
        double fullWidth = 0.0, tournamentWidth = 0.0;
        long fullTicks = 0, tournamentTicks = 0;
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(2026_1113 + 31 * r);
            double[] x = [.. Enumerable.Range(0, n + 1).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => 1.0 + 2.0 * v + 0.5 * Gaussian(rng))];
            double xt = x[n], yt = y[n];
            bool wide = r < widths;
            double[] cases = wide ? grid : [yt];
            long t0 = System.Diagnostics.Stopwatch.GetTimestamp();
            double[][] all = Metropolis(x[..n], y[..n], steps, keep, seed: r);
            var training = new double[keep, n];
            for (int d = 0; d < keep; d++)
                for (int i = 0; i < n; i++) training[d, i] = LogNormal(y[i], all[d][0] + all[d][1] * x[i]);
            ConformalSet full = ConformalBayes.FullStep(new PointwiseLogDensity(training, UnitPartition.EachObservation(n), []), Seam(all, xt, cases)).Sets[0];
            long t1 = System.Diagnostics.Stopwatch.GetTimestamp();
            ConformalSet tournament = CrossConformal.Tournament(n, i =>
            {
                double[] xi = [.. x[..i], .. x[(i + 1)..n]], yi = [.. y[..i], .. y[(i + 1)..n]];
                double[][] draws = Metropolis(xi, yi, steps, keep, seed: r);
                return (Seam(draws, x[i], [y[i]]), Seam(draws, xt, cases));
            }).Sets[0];
            long t2 = System.Diagnostics.Stopwatch.GetTimestamp();
            fullTicks += t1 - t0;
            tournamentTicks += t2 - t1;
            if (wide)
            {
                fullWidth += full.Width;
                tournamentWidth += tournament.Width;
                continue;
            }
            if (full.Accepts(0)) fullCovered++;
            if (tournament.Accepts(0)) tournamentCovered++;
        }
        int count = replications - widths;
        double f = fullCovered / (double)count, t = tournamentCovered / (double)count, se = Math.Sqrt(0.09 / count);
        output.WriteLine($"{chains} ({steps} steps, {keep} kept): coverage the full step {f:F4}, the tournament {t:F4} (se about {se:F4}; 1 − α 0.9000, "
                         + $"1 − (2k − 1)/m 0.8250); width {fullWidth / widths:F3} and {tournamentWidth / widths:F3}; "
                         + $"{fullTicks / (double)System.Diagnostics.Stopwatch.Frequency:F1} s and {tournamentTicks / (double)System.Diagnostics.Stopwatch.Frequency:F1} s");
        Assert.True(t >= 0.825 - 3.5 * Math.Sqrt(0.825 * 0.175 / count), $"the tournament's coverage {t:F4}");
    }

    [Fact]
    public void Coverage_WithBassChainsStoppedEarly_BesideTheFullStep()
    {
        // Slice 4's chain stopped early (no burn-in, 20 draws, one chain), n 39 (m 40, k 4), 100 replications of 8 held-out cases from the
        // truth outside every model; the full step on all 39 units and the tournament on each D₋ᵢ, every fit from the config's one seed.
        if (!BassQualification.Enabled) return;
        const int n = 39, replications = 100, tests = 8;
        BassConfig config = BassConformalDistributionTests.ConfigFor("a chain stopped early");
        int fullCovered = 0, tournamentCovered = 0;
        double fullWidth = 0.0, tournamentWidth = 0.0;
        double[] grid = [.. Enumerable.Range(0, 201).Select(g => -4.0 + 0.04 * g)];
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(2026_1114 + 37 * r);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(n, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(tests, rng);
            BassPredictor[] P(double[] a, double[] b) => [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];
            double[][] grids = [.. yt.Select(v => new[] { v }), grid];
            var cases = PredictionCases.Of(P([.. t0, t0[0]], [.. t1, t1[0]]));
            var data = new BassData(P(x0, x1), new ScalarResponse(y));
            BassFit fit = Bass.Fit(config, data);
            ConformalPrediction full = ConformalBayes.FullStep(fit.PointwiseLogDensity(data), fit.CandidateLogDensity(data, cases, grids));
            ConformalPrediction tournament = CrossConformal.Tournament(n, i =>
            {
                double[] a = [.. x0[..i], .. x0[(i + 1)..]], b = [.. x1[..i], .. x1[(i + 1)..]], c = [.. y[..i], .. y[(i + 1)..]];
                var rest = new BassData(P(a, b), new ScalarResponse(c));
                BassFit fi = Bass.Fit(config, rest);
                var heldOut = new PredictionCases([new ContinuousColumn("x0", [x0[i]]), new ContinuousColumn("x1", [x1[i]])]);
                return (fi.CandidateLogDensity(rest, heldOut, [[y[i]]]), fi.CandidateLogDensity(rest, cases, grids));
            });
            for (int t = 0; t < tests; t++)
            {
                if (full.Sets[t].Accepts(0)) fullCovered++;
                if (tournament.Sets[t].Accepts(0)) tournamentCovered++;
            }
            fullWidth += full.Sets[tests].Width;
            tournamentWidth += tournament.Sets[tests].Width;
        }
        int count = replications * tests;
        output.WriteLine($"BASS, a chain stopped early: coverage the full step {fullCovered / (double)count:F4}, the tournament {tournamentCovered / (double)count:F4} "
                         + $"over {count} cases; width {fullWidth / replications:F3} and {tournamentWidth / replications:F3}");
        Assert.True(tournamentCovered / (double)count >= 0.825 - 3.5 * Math.Sqrt(0.825 * 0.175 / count));
    }

    [Fact]
    public void CostRow()
    {
        // n 399 fits of S 1,000 draws and one grid of 200: the tournament's reduction from given seams, beside the full step's on one fit of
        // the same size. The n chains themselves are the engine's, n times one fit.
        if (!BassQualification.Enabled) return;
        const int n = 399, s = 1000, g = 200;
        var rng = new Xoshiro256PlusPlus(9);
        double[] grid = [.. Enumerable.Range(0, g).Select(v => -2.0 + 0.02 * v)];
        double[] mu = [.. Enumerable.Range(0, s).Select(_ => 0.2 * Gaussian(rng))];
        var cv = new double[s, g];
        for (int d = 0; d < s; d++)
            for (int k = 0; k < g; k++) cv[d, k] = LogNormal(grid[k], mu[d]);
        var candidates = new CandidateLogDensity(cv, [grid], []);
        double[] y = [.. Enumerable.Range(0, n).Select(_ => 0.5 * Gaussian(rng))];
        var held = new CandidateLogDensity[n];
        var training = new double[s, n];
        for (int i = 0; i < n; i++)
        {
            var hv = new double[s, 1];
            for (int d = 0; d < s; d++) training[d, i] = hv[d, 0] = LogNormal(y[i], mu[d]);
            held[i] = new CandidateLogDensity(hv, [[y[i]]], []);
        }
        long t0 = BassCost.Tick();
        CrossConformal.Tournament(n, i => (held[i], candidates));
        double tournament = (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e3;
        t0 = BassCost.Tick();
        ConformalBayes.FullStep(new PointwiseLogDensity(training, UnitPartition.EachObservation(n), []), candidates);
        double full = (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e3;
        output.WriteLine($"n {n}, S {s}, G {g} ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): the tournament's reduction {tournament:F0} ms, "
                         + $"the full step {full:F0} ms");
    }
}
