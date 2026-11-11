using System;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class ConformalKeptMarginTests(ITestOutputHelper output)
{
    private static readonly ConformalOptions Kept = new() { Alpha = 0.25, MarginBounds = MarginBounds.Kept };

    private static CandidateLogDensity Candidates(double[] beta, double x, double[] grid)
    {
        var ll = new double[beta.Length, grid.Length];
        for (int d = 0; d < beta.Length; d++)
            for (int g = 0; g < grid.Length; g++) ll[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * x);
        return new CandidateLogDensity(ll, [grid], ConformalBayesTests.Theta);
    }

    private static PointwiseLogDensity Training(double[] x, double[] y, double[] beta)
    {
        var ll = new double[beta.Length, x.Length];
        for (int d = 0; d < beta.Length; d++)
            for (int i = 0; i < x.Length; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
        return new PointwiseLogDensity(ll, UnitPartition.EachObservation(x.Length), ConformalBayesTests.Theta);
    }

    [Fact]
    public void EveryTurnIsChargedWhenTheFloorExceedsTheDrawCount()
    {
        // coverage_close_excluding, coverage_grid_excluding: no kept reading, every turn charged.
        double[] beta = Enumerable.Range(0, 128).Select(d => 0.1 * Math.Sin(d)).ToArray();
        double[] x = Enumerable.Repeat(1.0, 9).ToArray(), y = new double[9];
        var set = ConformalBayes.FullStep(Training(x, y, beta), Candidates(beta, 1, [-1, 0, 1]),
                                        Kept with { DegenerateEss = 10_000 }).Sets[0];
        Assert.Equal(10, set.Margin!.FlaggedCount);
        Assert.Equal(0, set.Margin.Count);
        Assert.All(set.Margin.CandidateFlagged!, Assert.True);
        Assert.All(set.Margin.Fragile!, f => Assert.Equal(-1, f));
        Assert.Equal(-0.22, set.CoverageUnderError, 12);
    }

    [Fact]
    public void KeptCountsMatchTheLiteralProxyAndEstimatedReading()
    {
        // card_excluded_le_estimated: the band plus D's inflated-accepted units, read on a kept ordering.
        var rng = new Xoshiro256PlusPlus(20261111);
        double[] x = [-2, -1, 0, 0.5, 1, 2, 3], y = [-1, 0, 0, 1, 1, 1, 6];
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 512, rng);
        var training = Training(x, y, beta);
        double[] grid = [-8, -2, 0, 2, 8];
        var candidates = Candidates(beta, 2, grid);
        foreach (double floor in new[] { 1.0, 50.0, 400.0 })
        {
            var set = ConformalBayes.FullStep(training, candidates, Kept with { DegenerateEss = floor }).Sets[0];
            var scorer = new ConformalBayes.StratumScorer(training.Values, Enumerable.Range(0, beta.Length).ToArray(), null,
                Enumerable.Range(0, x.Length).ToArray(), null, set.Level, set.Margin!.Multiple, MarginBounds.Kept, degenerateEss: floor);
            for (int g = 0; g < grid.Length; g++)
            {
                var sc = ConformalBayes.Scores(training, candidates, g);
                double[] weights = beta.Select(b => Math.Exp(ConformalBayesTests.LogNormal(grid[g], b * 2))).ToArray();
                double total = weights.Sum();
                for (int d = 0; d < weights.Length; d++) weights[d] /= total;
                bool[] flags = Enumerable.Range(0, x.Length).Select(i =>
                {
                    double ep = 0, inv = 0;
                    for (int d = 0; d < beta.Length; d++)
                    {
                        double p = Math.Exp(training.Values[d, i]);
                        ep += weights[d] * p;
                        inv += weights[d] / p;
                    }
                    return scorer.FreeFlagged[i] || (sc.EffectiveSampleSize >= floor && !(beta.Length / (ep * inv) >= floor));
                }).ToArray();
                bool ownFlag = !(sc.EffectiveSampleSize >= floor);
                Assert.Equal(flags.Count(f => f) + (ownFlag ? 1 : 0), set.Margin.Flagged![g]);
                Assert.Equal(ownFlag, set.Margin.CandidateFlagged![g]);
                if (ownFlag) continue;
                double[] scores = [.. sc.Scores, sc.Candidate];
                double[] errors = [.. Enumerable.Range(0, x.Length).Select(i => Math.Max(set.Margin.Multiple * sc.Errors[i], set.Margin.UnitCarriedErrors![i])),
                                    set.Margin.CandidateCarriedErrors![g]];
                Assert.Equal(ConformalBayes.EstimatedBand(scores, errors, set.Level), set.Margin.Fragile![g]);
                int accepted = Enumerable.Range(0, x.Length).Count(i => flags[i] &&
                    Enumerable.Range(0, scores.Length).Count(j => scores[j] - errors[j] <= scores[i] + errors[i]) > set.Level);
                Assert.Equal(accepted, set.Margin.DegenerateAccepted![g]);
            }
        }
    }

    [Fact]
    public void RefinementScoresTransitionsWithTheSameStratumScorer()
    {
        double[] x = [-1, 0, 1, 2, 3, 4, 5], y = [-1, 0, 1, 2, 3, 4, 5];
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 256, new Xoshiro256PlusPlus(72));
        var training = Training(x, y, beta);
        var candidates = Candidates(beta, 4, [-20, 4]);
        int calls = 0;
        var plain = ConformalBayes.FullStep(training, candidates, Kept with { DegenerateEss = 100 }).Sets[0];
        var refined = ConformalBayes.FullStep(training, candidates, Kept with
        {
            DegenerateEss = 100,
            MarginRefinement = new MarginRefinement((t, y0) => { calls++; return Candidates(beta, 4, [y0]); }, 0.01)
        }).Sets[0];
        Assert.True(calls > 0);
        Assert.Equal(calls, refined.Margin!.Refined!.Probes);
        Assert.True(refined.Margin.Count >= plain.Margin!.Count);
        Assert.True(refined.Margin.FlaggedCount >= plain.Margin.FlaggedCount);
        Assert.Equal(plain.PValues, refined.PValues);
    }

    [Fact]
    public void RefinementProbesOppositeSignFiniteExtremeEndpoints()
    {
        var ll = new double[128, 9];
        var lc = new double[128, 2];
        for (int d = 0; d < 128; d++)
        {
            for (int i = 0; i < 9; i++) ll[d, i] = -1;
            lc[d, 0] = d == 0 ? 0 : -100;
            lc[d, 1] = -1;
        }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(9), ConformalBayesTests.Theta);
        var candidates = new CandidateLogDensity(lc, [[-1e308, 1e308]], ConformalBayesTests.Theta);
        int calls = 0;
        var set = ConformalBayes.FullStep(training, candidates, Kept with
        {
            MarginRefinement = new MarginRefinement((t, y) =>
            {
                Assert.Equal(0, t);
                Assert.Equal(0.0, y);
                calls++;
                var at = new double[128, 1];
                for (int d = 0; d < 128; d++) at[d, 0] = -1;
                return new CandidateLogDensity(at, [[y]], ConformalBayesTests.Theta);
            }, 1e308)
        }).Sets[0];
        Assert.NotEqual(set.Margin!.CandidateFlagged![0], set.Margin.CandidateFlagged[1]);
        Assert.Equal(1, calls);
        Assert.Equal(1, set.Margin.Refined!.Probes);
    }

    [Fact]
    public void CarriedBoundsDominateEveryRepresentedKeptOrdering()
    {
        // coverage_close_excluding's hfail: dominance of the computed standard errors, not a sure bound on exact score error.
        const int n = 7, s = 512;
        double[] x = [-2, -1, 0, 0.5, 1, 2, 3], y = [-1, 0, 0, 1, 1, 1, 6];
        double[] beta = ConformalBayesTests.ExactDraws(x, y, s, new Xoshiro256PlusPlus(311));
        var training = Training(x, y, beta);
        var candidates = Candidates(beta, 2, [-4, 0, 4]);
        var set = ConformalBayes.FullStep(training, candidates, Kept with { DegenerateEss = 100 }).Sets[0];
        var scorer = new ConformalBayes.StratumScorer(training.Values, Enumerable.Range(0, s).ToArray(), null,
            Enumerable.Range(0, n).ToArray(), null, set.Level, set.Margin!.Multiple, MarginBounds.Kept, degenerateEss: 100);
        var aug = new double[s, n + 3];
        for (int d = 0; d < s; d++)
        {
            for (int i = 0; i < n; i++) aug[d, i] = training.Values[d, i];
            for (int g = 0; g < 3; g++) aug[d, n + g] = candidates.Values[d, g];
        }
        var augmented = new PointwiseLogDensity(aug, UnitPartition.EachObservation(n + 3), ConformalBayesTests.Theta);
        for (int i = 0; i < n; i++)
        {
            if (scorer.FreeFlagged[i]) continue;
            var column = new double[s, 1];
            for (int d = 0; d < s; d++) column[d, 0] = training.Values[d, i];
            var sc = ConformalBayes.Scores(augmented, new CandidateLogDensity(column, [[0.0]], ConformalBayesTests.Theta), 0);
            for (int j = 0; j < n; j++)
                Assert.True(set.Margin.UnitCarriedErrors![j] + 1e-10 >= set.Margin.Multiple * sc.Errors[j]);
            for (int g = 0; g < 3; g++)
                if (!set.Margin.CandidateFlagged![g])
                    Assert.True(set.Margin.CandidateCarriedErrors![g] + 1e-10 >= set.Margin.Multiple * sc.Errors[n + g]);
        }
    }

    [Fact]
    public void ProxyAgainstExactConjugateOwnTurnRefits()
    {
        if (!BassQualification.Enabled) return;
        const int n = 19, draws = 512, refitDraws = 16384, fixtures = 24;
        int[] seeds = SeedTree.Derive(2026111101, fixtures);
        foreach (double floor in new[] { 10.0, 100.0, 400.0 })
        {
            int flagged = 0, exactFlagged = 0, missed = 0, conservative = 0, turns = 0, empiricalMissed = 0, refitDisagreements = 0;
            int flaggedYTurns = 0, flaggedYMissed = 0, flaggedYExact = 0, unionMissed = 0;
            double worstRatio = 0;
            foreach (int seed in seeds)
            {
                var rng = new Xoshiro256PlusPlus(seed);
                double[] x = Enumerable.Range(0, n).Select(_ => 4 * rng.NextDouble() - 2).ToArray();
                double[] y = x.Select(v => ConformalBayesTests.Truth(v, 0.5, rng)).ToArray();
                double xt = 4 * rng.NextDouble() - 2;
                double[] grid = [-6, -2, 0, 2, 6];
                double[] beta = ConformalBayesTests.ExactDraws(x, y, draws, rng);
                var training = Training(x, y, beta);
                var candidate = Candidates(beta, xt, grid);
                var set = ConformalBayes.FullStep(training, candidate, Kept with { DegenerateEss = floor }).Sets[0];
                var scorer = new ConformalBayes.StratumScorer(training.Values, Enumerable.Range(0, draws).ToArray(), null,
                    Enumerable.Range(0, n).ToArray(), null, set.Level, set.Margin!.Multiple, MarginBounds.Kept, degenerateEss: floor);
                for (int i = 0; i < n; i++)
                {
                    double[] rx = x.Where((_, j) => j != i).ToArray(), ry = y.Where((_, j) => j != i).ToArray();
                    double[] free = ConformalBayesTests.ExactDraws(rx, ry, refitDraws, rng);
                    double Ess(double[] b)
                    {
                        double[] log = b.Select(v => ConformalBayesTests.LogNormal(y[i], v * x[i])).ToArray();
                        double top = log.Max(), sum = 0, squares = 0;
                        foreach (double l in log) { double w = Math.Exp(l - top); sum += w; squares += w * w; }
                        return draws * sum * sum / (squares * log.Length);
                    }
                    double freeEss = Ess(free);
                    for (int g = 0; g < grid.Length; g++)
                    {
                        // Independent posterior refit of D minus i plus candidate; p_i is now the add-one-in ratio.
                        double[] own = ConformalBayesTests.ExactDraws([.. rx, xt], [.. ry, grid[g]], refitDraws, rng);
                        double empirical = Ess(own), precision = 0.1, xy = 0;
                        for (int j = 0; j < rx.Length; j++) { precision += rx[j] * rx[j]; xy += rx[j] * ry[j]; }
                        precision += xt * xt; xy += xt * grid[g];
                        double variance = 1 / precision, mean = xy / precision, residual = y[i] - x[i] * mean;
                        double Moment(int q) => -q * Math.Log(2 * Math.PI) / 2 - Math.Log(1 + q * x[i] * x[i] * variance) / 2
                            - q * residual * residual / (2 * (1 + q * x[i] * x[i] * variance));
                        double exact = draws * Math.Exp(2 * Moment(1) - Moment(2));
                        double[] w = beta.Select(v => Math.Exp(ConformalBayesTests.LogNormal(grid[g], v * xt))).ToArray();
                        double total = w.Sum(), ep = 0, inv = 0, squares = 0;
                        for (int d = 0; d < draws; d++)
                        {
                            w[d] /= total; squares += w[d] * w[d];
                            double p = Math.Exp(training.Values[d, i]); ep += w[d] * p; inv += w[d] / p;
                        }
                        double proxyEss = draws / (ep * inv);
                        bool ownFlag = 1 / squares < floor;
                        bool proxy = scorer.FreeFlagged[i] || (!ownFlag && !(proxyEss >= floor)), actual = exact < floor;
                        if (ownFlag)
                        {
                            flaggedYTurns++; if (actual) flaggedYExact++; if (actual && !proxy) flaggedYMissed++;
                            continue;
                        }
                        turns++; if (proxy) flagged++; if (actual) exactFlagged++;
                        if (actual && !proxy) missed++; if (proxy && !actual) conservative++;
                        if (empirical < floor && !proxy) empiricalMissed++;
                        if ((empirical < floor) != actual) refitDisagreements++;
                        if (Math.Min(freeEss, empirical) < floor && !proxy) unionMissed++;
                        worstRatio = Math.Max(worstRatio, proxyEss / exact);
                    }
                }
            }
            output.WriteLine($"exact own-turn refits: floor {floor}, S {draws}, refit S {refitDraws}, {turns} turns at kept y; proxy flags {flagged}, analytic own ESS flags {exactFlagged}, misses {missed}, conservative flags {conservative}, empirical-own misses {empiricalMissed}, empirical/analytic classifications disagree {refitDisagreements}, free/own union misses {unionMissed}, largest proxy/analytic ESS {worstRatio:G5}; at flagged y {flaggedYTurns} turns, analytic own flags {flaggedYExact}, y-free proxy misses {flaggedYMissed}. This is a measurement, not hflag or hcount.");
        }
    }

    [Theory]
    [InlineData(64)]
    [InlineData(512)]
    public void CoverageOnExactDrawsWithTheKeptReading(int draws)
    {
        if (!BassQualification.Enabled) return;
        const int n = 19, replications = 1500;
        int[] seeds = SeedTree.Derive(2026111102, replications);
        foreach (double floor in new[] { 10.0, 100.0, 400.0 })
        {
            int covered = 0, inflated = 0;
            double reading = 0, inflatedReading = 0, flags = 0, band = 0;
            foreach (int seed in seeds)
            {
                var rng = new Xoshiro256PlusPlus(seed);
                double[] x = Enumerable.Range(0, n).Select(_ => 4 * rng.NextDouble() - 2).ToArray();
                double[] y = x.Select(v => ConformalBayesTests.Truth(v, 0.5, rng)).ToArray();
                double xt = 4 * rng.NextDouble() - 2, yt = ConformalBayesTests.Truth(xt, 0.5, rng);
                double[] beta = ConformalBayesTests.ExactDraws(x, y, draws, rng);
                var set = ConformalBayes.FullStep(Training(x, y, beta), Candidates(beta, xt, [yt]),
                    Kept with { Alpha = 0.1, DegenerateEss = floor }).Sets[0];
                if (set.AcceptsPlain(0)) covered++; if (set.AcceptsInflated(0)) inflated++;
                reading += set.CoverageUnderError; inflatedReading += (set with { Inflated = true }).CoverageUnderError;
                flags += set.Margin!.FlaggedCount!.Value; band += set.Margin.Count!.Value;
            }
            double se = Math.Sqrt(0.09 / replications);
            output.WriteLine($"S {draws}, floor {floor}: coverage {covered / (double)replications:F4}, inflated {inflated / (double)replications:F4}, se {se:F4}; singleton true-value readings plain {reading / replications:F4}, inflated {inflatedReading / replications:F4}, flags {flags / replications:F2}, band {band / replications:F2}. Singleton readings are not a response-uniform grid bound.");
            Assert.True(covered / (double)replications >= reading / replications - 3.5 * se);
            Assert.True(inflated / (double)replications >= inflatedReading / replications - 3.5 * se);
        }
    }

    [Theory]
    [InlineData(InflationBounds.PerPosition)]
    [InlineData(InflationBounds.Uniform)]
    public void AFarCandidateDoesNotLendItsLowEssToConstantTrainingTurns(InflationBounds inflation)
    {
        var ll = new double[128, 9]; var lc = new double[128, 2];
        for (int d = 0; d < 128; d++)
        {
            for (int i = 0; i < 9; i++) ll[d, i] = -1;
            lc[d, 0] = d == 0 ? 0 : -100;
            lc[d, 1] = -1;
        }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(9), ConformalBayesTests.Theta);
        var candidate = new CandidateLogDensity(lc, [[-20.0, 0.0]], ConformalBayesTests.Theta);
        var options = Kept with { DegenerateEss = 100, InflationBounds = inflation };
        var set = ConformalBayes.FullStep(training, candidate, options).Sets[0];
        var original = ConformalBayes.FullStep(training, candidate, options with { MarginBounds = null }).Sets[0];
        Assert.Equal(1, set.Margin!.Flagged![0]);
        Assert.True(set.Margin.CandidateFlagged![0]);
        Assert.Equal(0, set.Margin.Flagged[1]);
        Assert.Equal(original.Margin!.InflatedRanks, set.Margin.InflatedRanks);
    }

    [Theory]
    [InlineData(1000)]
    [InlineData(4000)]
    public void KeptCostRow(int s)
    {
        if (!BassQualification.Enabled) return;
        const int n = 1000, g = 200;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(2026111104, 2)[s == 1000 ? 0 : 1]);
        var ll = new double[s, n]; var lc = new double[s, g];
        for (int d = 0; d < s; d++)
        {
            for (int i = 0; i < n; i++) ll[d, i] = -1 - 4 * rng.NextDouble();
            for (int j = 0; j < g; j++) lc[d, j] = -1 - 6 * rng.NextDouble() * (1 + j / 50.0);
        }
        int[] rows = Enumerable.Range(0, s).ToArray(), members = Enumerable.Range(0, n).ToArray();
        ConformalBayes.StratumScorer scorer = null!;
        var reading = new ConformalBayes.Reading(g);
        double Best(Action call)
        {
            double best = double.PositiveInfinity;
            call();
            for (int repeat = 0; repeat < 3; repeat++)
            {
                long start = BassCost.Tick(); call();
                best = Math.Min(best, 1000 * (BassCost.Tick() - start) / BassCost.TicksPerSecond);
            }
            return best;
        }
        double build = Best(() => scorer = new ConformalBayes.StratumScorer(ll, rows, null, members, null, 100, 3.5, MarginBounds.Kept));
        double weights = Best(() => scorer.Weigh(lc, 0, g, reading));
        double product = Best(() => scorer.Multiply(g));
        double count = Best(() => { scorer.Multiply(g); scorer.Count(lc, 0, g, reading); }) - product;
        output.WriteLine($"Kept cost, n {n}, S {s}, G {g}, {(BassCost.Cycles ? "cycle" : "wall")} clock, best of three indicative: once/stratum {build:F2} ms; per grid weights {weights:F2}, products {product:F2}, counts {count:F2} ms; per candidate {weights * 1000 / g:F2}, {product * 1000 / g:F2}, {count * 1000 / g:F2} us.");
    }
}
