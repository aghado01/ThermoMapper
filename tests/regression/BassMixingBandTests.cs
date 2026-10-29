using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's M5 (the local-move measurements brief): the steep band of S10's tempered Friedman fit.
/// <list type="bullet">
///   <item><b>Part 1:</b> S10's run, each chain's ladder self-sized during its 5,000 scans of burn-in as the run sizes it,
///   then 2,000 scans sampled on the frozen ladder. Each rung's mean and sd of M and of ℓ, the untempered log marginal
///   ln p₁(y | structure, τ), over the sampling phase, read from the replicas after each scan's swaps, beside the
///   readout's per-pair rejections.</item>
///   <item><b>Part 2:</b> per chain, a dense fixed ladder through the band, adaptation off: the chain's own frozen rungs
///   above 0.40, then 31 rungs from 0.40 to 0.10 at Δb = 0.01, then its own rungs below 0.10 down to its hot end 1/N. The
///   same fixture, start, seed and lengths; the same rows and rejections.</item>
///   <item><b>The Gaussian reading</b> (the design draft's lever 4a): each pair's increment h(r) = √(2/π)·Φ⁻¹((1 + r)/2),
///   with 1 − r floored at 1/(2·scans); the band's Σh against its Σr on both ladders, over the window [0.10, 0.40] and
///   over the span of Part 1's pairs above 0.9, a pair straddling a window's edge counted by the share of its Δb
///   inside.</item>
/// </list>
/// It reads, it does not judge; it asserts that each rung's rows cover every sampled scan and that the readout's rejections
/// lie in [0, 1].
/// </summary>
public sealed class BassMixingBandTests
{
    private const double Top = 0.40, Bottom = 0.10, Step = 0.01;
    private readonly ITestOutputHelper _out;
    public BassMixingBandTests(ITestOutputHelper output) => _out = output;

    // One ladder's sampling phase: its rungs, each rung's running moments of M and ℓ, and each pair's rejection.
    private sealed record Rungs(double[] InverseTemperatures, RunningMoments[] M, RunningMoments[] LogLikelihood, double[] Rejection, double RoundTrips, double Ideal);

    [Fact]
    public void SteepBand_FrozenAndDenseLadders()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman, new TemperingSpec());
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        int chains = fixture.Config.Sampling.Chains;

        var clock = System.Diagnostics.Stopwatch.StartNew();
        var coarse = new Rungs[chains];
        Parallel.For(0, chains, c =>
        {
            BassLadder ladder = BassMixing.Ladder(model, c, fixture.Config.Sampling.Tempering!);
            BassMixing.Burn(ladder);
            coarse[c] = Sample(ladder);
        });
        double coarseSeconds = clock.Elapsed.TotalSeconds;

        clock.Restart();
        var dense = new Rungs[chains];
        Parallel.For(0, chains, c =>
        {
            double[] inverseTemperatures = Dense(coarse[c].InverseTemperatures);
            var spec = new TemperingSpec { Levels = inverseTemperatures.Length, Adaptive = false, MaxLevels = Math.Max(64, inverseTemperatures.Length) };
            BassLadder ladder = BassMixing.Ladder(model, c, spec, inverseTemperatures);
            BassMixing.Burn(ladder);
            dense[c] = Sample(ladder);
        });
        double denseSeconds = clock.Elapsed.TotalSeconds;

        _out.WriteLine($"## M5, Friedman tempered: Part 1 self-sized ({coarseSeconds:F1} s), Part 2 dense through [{Bottom}, {Top}] at Δb = {Step} ({denseSeconds:F1} s); " +
                       $"{BassMixing.BurnIn:N0} scans of burn-in, {BassMixing.Samples:N0} sampled");
        for (int c = 0; c < chains; c++)
        {
            Table($"Part 1, chain {c}", coarse[c]);
            Table($"Part 2, chain {c}", dense[c]);
        }
        Band(coarse, dense);
    }

    // The sampling phase on a ladder past its burn-in: each scan's rungs read after its swaps, and the readout's rejections.
    private static Rungs Sample(BassLadder ladder)
    {
        int levels = ladder.Tempering.Replicas.Count;
        var m = new RunningMoments[levels];
        var ell = new RunningMoments[levels];
        for (int s = 0; s < BassMixing.Samples; s++)
        {
            ladder.Step();
            for (int i = 0; i < levels; i++)
            {
                var replica = ladder.Tempering.Replicas[i];
                m[i].Add(replica.Current.Count);
                ell[i].Add(((BassEvaluation)replica.CurrentEvaluation).LogMarginal);
            }
        }
        BassTemperingReadout readout = BassTemperingReadout.Of([ladder]);
        double[] rejection = readout.PairRejection[0];
        Assert.All(rejection, r => Assert.InRange(r, 0.0, 1.0 + 1e-12));
        Assert.All(m, x => Assert.Equal(BassMixing.Samples, x.Count));
        return new Rungs(ladder.Tempering.InverseTemperatures, m, ell, rejection, readout.RoundTripsPerScan, readout.IdealRoundTripRate);
    }

    // The chain's own rungs above the window, the window at Δb = 0.01, then its own rungs below it.
    private static double[] Dense(double[] frozen)
    {
        const double slack = 1e-9;
        int count = (int)Math.Round((Top - Bottom) / Step) + 1;
        return [.. frozen.Where(b => b > Top + slack), .. Enumerable.Range(0, count).Select(k => Math.Round(Top - k * Step, 2)), .. frozen.Where(b => b < Bottom - slack)];
    }

    // Each rung's row, and the pair below it (to the next hotter rung): its rejection and its Gaussian increment.
    private void Table(string label, Rungs rungs)
    {
        _out.WriteLine("");
        _out.WriteLine($"{label}: {rungs.InverseTemperatures.Length} rungs, round trips per scan {rungs.RoundTrips:G3} against the ideal {rungs.Ideal:G3}, Λ̂ = Σr {rungs.Rejection.Sum():F2}, Σh {rungs.Rejection.Sum(H):F2}");
        _out.WriteLine("rung | b | M mean | M sd | ℓ mean | ℓ sd | pair to the next | r | h(r)");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int i = 0; i < rungs.InverseTemperatures.Length; i++)
        {
            string pair = i + 1 < rungs.InverseTemperatures.Length ? $"[{rungs.InverseTemperatures[i + 1]:G3}, {rungs.InverseTemperatures[i]:G3}] | {rungs.Rejection[i]:F3} | {H(rungs.Rejection[i]):F2}" : "— | — | —";
            _out.WriteLine($"{i} | {rungs.InverseTemperatures[i]:G4} | {rungs.M[i].Mean:F2} | {rungs.M[i].Sd:F2} | {rungs.LogLikelihood[i].Mean:F1} | {rungs.LogLikelihood[i].Sd:F1} | {pair}");
        }
    }

    // The band per chain: Σr and Σh over the window and over Part 1's steep span, on both ladders.
    private void Band(Rungs[] coarse, Rungs[] dense)
    {
        _out.WriteLine("");
        _out.WriteLine("chain | Part 1's pairs above 0.9 | steep span | Part 1 Σr, Σh over the span | Part 2 Σr, Σh over the span | Part 2 pairs above 0.9 | " +
                       $"Part 1 Σr, Σh over [{Bottom}, {Top}] | Part 2 Σr, Σh over [{Bottom}, {Top}]");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < coarse.Length; c++)
        {
            Rungs a = coarse[c], d = dense[c];
            int[] steep = [.. Enumerable.Range(0, a.Rejection.Length).Where(i => a.Rejection[i] > 0.9)];
            string pairs = steep.Length == 0 ? "none" : string.Join(", ", steep.Select(i => $"[{a.InverseTemperatures[i + 1]:G3}, {a.InverseTemperatures[i]:G3}] {a.Rejection[i]:F3}"));
            (double lo, double hi) = steep.Length == 0 ? (double.NaN, double.NaN) : (steep.Min(i => a.InverseTemperatures[i + 1]), steep.Max(i => a.InverseTemperatures[i]));
            string over(Rungs r, double l, double u) => double.IsNaN(l) ? "—" : $"{Sum(r, l, u, x => x):F2}, {Sum(r, l, u, H):F2}";
            int[] denseSteep = [.. Enumerable.Range(0, d.Rejection.Length).Where(i => d.Rejection[i] > 0.9)];
            string densePairs = denseSteep.Length == 0 ? "none" : string.Join(", ", denseSteep.Select(i => $"[{d.InverseTemperatures[i + 1]:G3}, {d.InverseTemperatures[i]:G3}] {d.Rejection[i]:F3}"));
            _out.WriteLine($"{c} | {pairs} | {(steep.Length == 0 ? "—" : $"[{lo:G3}, {hi:G3}]")} | {over(a, lo, hi)} | {over(d, lo, hi)} | {densePairs} | " +
                           $"{over(a, Bottom, Top)} | {over(d, Bottom, Top)}");
        }
    }

    // Σ f(r) over the ladder's pairs within [lo, hi], each pair weighted by the share of its Δb inside.
    private static double Sum(Rungs rungs, double lo, double hi, Func<double, double> f)
    {
        double sum = 0.0;
        for (int i = 0; i < rungs.Rejection.Length; i++)
        {
            double upper = rungs.InverseTemperatures[i], lower = rungs.InverseTemperatures[i + 1];
            double inside = Math.Min(upper, hi) - Math.Max(lower, lo);
            if (inside > 0.0) sum += inside / (upper - lower) * f(rungs.Rejection[i]);
        }
        return sum;
    }

    // The Gaussian reading's increment, √(2/π)·Φ⁻¹((1 + r)/2), its 1 − r floored at 1/(2·scans).
    private static double H(double r)
    {
        double floored = Math.Clamp(r, 0.0, 1.0 - 1.0 / (2.0 * BassMixing.Samples));
        return Math.Sqrt(2.0 / Math.PI) * Normal.StandardQuantile((1.0 + floored) / 2.0);
    }

    // Welford's running mean and sd.
    private struct RunningMoments
    {
        private double _mean, _squares;
        public int Count { get; private set; }
        public readonly double Mean => _mean;
        public readonly double Sd => Count > 1 ? Math.Sqrt(_squares / (Count - 1)) : double.NaN;

        public void Add(double x)
        {
            Count++;
            double delta = x - _mean;
            _mean += delta / Count;
            _squares += delta * (x - _mean);
        }
    }
}
