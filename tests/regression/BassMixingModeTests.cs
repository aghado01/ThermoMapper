using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's M4 (the local-move measurements brief): the chains' modes. S10's untempered fits of its three fixtures, each
/// chain run as the run runs it (5,000 sweeps of burn-in, 2,000 sampled), and over the sampling phase each state's
/// unnormalized log posterior, the log marginal plus the log structure prior (τ's prior density left out), its M, and
/// its structure by the variables and orders of its terms. Per chain: the log posterior's quantiles, M's band and the
/// structures held; across chains: their medians' differences and whether one chain's states reach another's range.
/// It reads, it does not judge; it asserts that the harness's chains are S10's (the run's draws hold the same M, draw
/// for draw).
/// </summary>
public sealed class BassMixingModeTests
{
    private readonly ITestOutputHelper _out;
    public BassMixingModeTests(ITestOutputHelper output) => _out = output;

    // One sampled state: its M, log posterior and structure by variables and orders.
    private sealed record Held(int M, double LogPosterior, string Structure);

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void ChainsModes_OverTheSamplingPhase(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassMixing.Fixture fixture = BassMixing.Get(name);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassTarget target = model.Target;
        int chains = fixture.Config.Sampling.Chains;
        var held = new Held[chains][];
        Parallel.For(0, chains, c =>
        {
            BassChain chain = BassMixing.Chain(model, c);
            for (int s = 0; s < BassMixing.BurnIn; s++) chain.Step();
            held[c] = new Held[BassMixing.Samples];
            for (int s = 0; s < BassMixing.Samples; s++)
            {
                BassState state = chain.Step();
                var evaluation = (BassEvaluation)chain.Moves.CurrentEvaluation;
                double logPosterior = evaluation.LogMarginal + target.LogStructurePrior(state.Fit);
                Assert.True(double.IsFinite(logPosterior));
                held[c][s] = new Held(state.Count, logPosterior, string.Join(" + ", state.Terms.Select(fixture.Variables).Order(StringComparer.Ordinal)));
            }
        });
        double seconds = clock.Elapsed.TotalSeconds;

        // The harness's chains are S10's: the run's retained draws hold the same M, draw for draw.
        BassFit fit = Bass.Fit(fixture.Config, fixture.Data);
        for (int c = 0; c < chains; c++)
            Assert.Equal(held[c].Select(h => h.M), fit.Draws[c].Select(d => d.Structure.TermIds.Length));
        _out.WriteLine($"## M4, {name}: {chains} chains, {BassMixing.Samples:N0} sampled states each ({seconds:F1} s; the run's R̂ of M {fit.Diagnostics.Count.RHat.Value:F2})");

        _out.WriteLine("");
        _out.WriteLine("chain | mean M | M range | M 5%–95% | log post min | 5% | 25% | median | 75% | 95% | max | structures | top three's share");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < chains; c++)
        {
            double[] lp = [.. held[c].Select(h => h.LogPosterior)];
            double[] m = [.. held[c].Select(h => (double)h.M)];
            int top = held[c].GroupBy(h => h.Structure).Select(g => g.Count()).OrderDescending().Take(3).Sum();
            _out.WriteLine($"{c} | {m.Average():F2} | {m.Min()}–{m.Max()} | {BassMixing.Quantile(m, 0.05)}–{BassMixing.Quantile(m, 0.95)} | {lp.Min():F1} | " +
                           $"{string.Join(" | ", new[] { 0.05, 0.25, 0.5, 0.75, 0.95 }.Select(q => BassMixing.Quantile(lp, q).ToString("F1")))} | {lp.Max():F1} | " +
                           $"{held[c].Select(h => h.Structure).Distinct().Count()} | {(double)top / held[c].Length:P1}");
        }

        // Each chain's three most frequent structures, their share in every chain and their median log posterior there.
        _out.WriteLine("");
        _out.WriteLine($"chain | rank | structure | M | {string.Join(" | ", Enumerable.Range(0, chains).Select(j => $"share in {j}"))} | median log post");
        _out.WriteLine("--- | --- | --- | ---" + string.Concat(Enumerable.Repeat(" | ---", chains + 1)));
        for (int c = 0; c < chains; c++)
            foreach ((IGrouping<string, Held> g, int rank) in held[c].GroupBy(h => h.Structure).OrderByDescending(g => g.Count()).Take(3).Select((g, r) => (g, r + 1)))
                _out.WriteLine($"{c} | {rank} | {g.Key} | {g.First().M} | " +
                               string.Join(" | ", Enumerable.Range(0, chains).Select(j => $"{(double)held[j].Count(h => h.Structure == g.Key) / held[j].Length:P1}")) +
                               $" | {BassMixing.Quantile(g.Select(h => h.LogPosterior), 0.5):F1}");

        // Across chains: each median against the highest, and the share of chain i's states inside chain j's range, then at
        // or above chain j's median.
        double[] medians = [.. held.Select(h => BassMixing.Quantile(h.Select(x => x.LogPosterior), 0.5))];
        _out.WriteLine("");
        _out.WriteLine($"median log post less the highest chain's: {string.Join(", ", medians.Select((x, c) => $"chain {c} {x - medians.Max():F1}"))}");
        _out.WriteLine("");
        _out.WriteLine($"chain i's states in chain j's range / at or above j's median | {string.Join(" | ", Enumerable.Range(0, chains).Select(j => $"j = {j}"))}");
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", chains)));
        for (int i = 0; i < chains; i++)
            _out.WriteLine($"i = {i} | " + string.Join(" | ", Enumerable.Range(0, chains).Select(j =>
            {
                double lo = held[j].Min(h => h.LogPosterior), hi = held[j].Max(h => h.LogPosterior);
                double inRange = held[i].Count(h => h.LogPosterior >= lo && h.LogPosterior <= hi) / (double)held[i].Length;
                double above = held[i].Count(h => h.LogPosterior >= medians[j]) / (double)held[i].Length;
                return $"{inRange:P1} / {above:P1}";
            })));
    }
}
