using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Tempering for BARS, measured. Every run uses BARS's default moves: births from the uniform base (ε = 1, step 3.9)
/// with relocations by the local kernel.
/// <list type="bullet">
/// <item><c>TEMPER_BENCH</c>: round trips per unit work (step 3.3, Q5's revisit trigger), by schedule, ladder adaptation
/// and path. Round trips compare like with like only between even-odd schedules on one path: they are journeys between
/// a ladder's own ends, which differ by path, and the sequential sweep's cascades let a state cross several rungs with no
/// local step between.</item>
/// <item><c>TEMPER_BENCH_ESS</c>: the decision study for OQ1 and OQ4. Adapted DEO ladders on the likelihood path with
/// the hot end at β_min = 0.05, 0.2, 0.35 and 0.5, and on the powered path at 0.05, judged by the cold chain's ESS of the
/// knot count, the energy (its log likelihood) and the peak location (the argmax of the conditional-mean curve). Each
/// estimate pools 8 ladders, and two independent groups show its noise.</item>
/// <item><c>TEMPER_BENCH_BASE</c>: the equal-work baseline, untempered chains as many as the ladders' replicas
/// (8 × L), each run as long.</item>
/// <item><c>TEMPER_BENCH_OQ2</c>: the OQ1 rerun with exact reference rungs (OQ2), pre-registered in the preparation
/// brief: the Normal model's integrated and linear tempered forms at the four hot ends, each path ending at its exact
/// reference rung, and the powered path at 0.05, by the same cold-chain ESS, with each arm's communication barrier and
/// round trips beside it.</item>
/// </list>
/// Work is counted in replica-steps of the measured phase, which tempered and untempered runs spend alike; seconds are
/// reported beside it. Each variable names the TSV file written; <c>TEMPER_BENCH_FITS</c> (comma-separated) restricts the
/// fits.
/// </summary>
public sealed class TemperingBenchmark
{
    private const int Burn = 4095;      // 1 + 2 + … + 2048: Algorithm 4's rounds exactly
    private const int Measured = 20_000;

    private sealed record Fit(string Name, IObservationModel Model, double[] X, double[] Y, double PriorMean);

    private static Fit[] Fits(params string[] names)
    {
        (double[] sx, double[] sy) = Sine();
        (double[] px, double[] py) = Peak();
        (double[] cx, double[] cy) = Counts();
        (double[] ox, double[] oy) = Oscillatory();
        Fit[] all =
        [
            new("peak", new WeightedNormalModel(), px, py, 5.0),
            new("sine", new WeightedNormalModel(), sx, sy, 5.0),
            new("counts", new PoissonModel(), cx, cy, 5.0),
            new("oscillatory", new WeightedNormalModel(), ox, oy, 15.0),
        ];
        string[]? only = Environment.GetEnvironmentVariable("TEMPER_BENCH_FITS")?.Split(',');
        return all.Where(f => names.Contains(f.Name) && (only is null || only.Contains(f.Name))).ToArray();
    }

    // BARS's default moves: one birth density shared by birth and death, from the uniform base.
    private static IRjMove<KnotConfig>[] Moves()
    {
        var kernel = new LocalBetaKernel(50.0);
        var density = new KnotBirthDensity(kernel, exploration: 1.0);
        return [new KnotBirthMove(density), new KnotDeathMove(density), new KnotRelocateMove(kernel)];
    }

    [Fact]
    public void RoundTrips_PerUnitWork()
    {
        if (Environment.GetEnvironmentVariable("TEMPER_BENCH") is not string path) return;   // heavy; set TEMPER_BENCH to run

        var sb = new StringBuilder();
        sb.AppendLine("data\tpath\tschedule\tadapted\tlevels\tseed\ttrips_per_scan\ttrips_per_replica_step\ttrips_per_second\t" +
                      "barrier\tideal_trips_per_scan\tcold_knots\thot_knots\tswap_rate\tladder");
        (TemperingPath Path, SwapSchedule Schedule, bool Adapt)[] configs =
        [
            (TemperingPath.Likelihood, SwapSchedule.DeterministicEvenOdd, true),
            (TemperingPath.Powered, SwapSchedule.DeterministicEvenOdd, true),
            (TemperingPath.Likelihood, SwapSchedule.DeterministicEvenOdd, false),
            (TemperingPath.Likelihood, SwapSchedule.Sequential, false),
            (TemperingPath.Powered, SwapSchedule.Sequential, false),
            (TemperingPath.Likelihood, SwapSchedule.Sequential, true),
            (TemperingPath.Powered, SwapSchedule.Sequential, true),
        ];
        foreach (Fit fit in Fits("sine", "peak", "counts"))
            foreach ((TemperingPath tp, SwapSchedule schedule, bool adapt) in configs)
                foreach (int levels in new[] { 4, 8 })
                    for (int seed = 1; seed <= 3; seed++)
                        sb.AppendLine(Run(fit, tp, schedule, adapt, levels, 0.05, seed).Row);
        File.WriteAllText(path, sb.ToString());
    }

    /// <summary>
    /// The decision study for OQ1 and OQ4: adapted DEO ladders on the likelihood path with the hot end at four β_min,
    /// and on the powered path, by the cold chain's ESS per unit work. Guarded by <c>TEMPER_BENCH_ESS</c>.
    /// </summary>
    [Fact]
    public void HotEnds_ByColdChainEss()
    {
        if (Environment.GetEnvironmentVariable("TEMPER_BENCH_ESS") is not string path) return;   // heavy; set TEMPER_BENCH_ESS to run

        const int chains = 8;
        (TemperingPath Path, double BetaMin)[] arms =
        [
            (TemperingPath.Likelihood, 0.05), (TemperingPath.Likelihood, 0.2), (TemperingPath.Likelihood, 0.35),
            (TemperingPath.Likelihood, 0.5), (TemperingPath.Powered, 0.05),
        ];
        var ess = new StringBuilder(EssHeader + "\n");
        foreach (Fit fit in Fits("peak", "sine", "counts", "oscillatory"))
            foreach ((TemperingPath tp, double betaMin) in arms)
                foreach (int levels in new[] { 4, 8 })
                    for (int group = 0; group < 2; group++)
                    {
                        var series = new Series(chains);
                        for (int c = 0; c < chains; c++)
                        {
                            RunResult r = Run(fit, tp, SwapSchedule.DeterministicEvenOdd, adapt: true, levels, betaMin,
                                              seed: 1000 * (group + 1) + c);
                            series.Add(c, fit, r, levels);
                        }
                        ess.AppendLine(series.Row(fit.Name, tp.ToString(), betaMin, levels, group));
                        File.WriteAllText(path, ess.ToString());   // partial results survive an interrupted run
                    }
    }

    /// <summary>
    /// The equal-work baseline for <see cref="HotEnds_ByColdChainEss"/>: untempered chains, as many as the ladders'
    /// replicas (8 ladders × L rungs), each run as long. Guarded by <c>TEMPER_BENCH_BASE</c>.
    /// </summary>
    [Fact]
    public void Untempered_EqualWorkBaseline()
    {
        if (Environment.GetEnvironmentVariable("TEMPER_BENCH_BASE") is not string path) return;   // heavy; set TEMPER_BENCH_BASE to run

        var ess = new StringBuilder(EssHeader + "\n");
        foreach (Fit fit in Fits("peak", "sine", "counts", "oscillatory"))
            foreach (int levels in new[] { 4, 8 })
            {
                int chains = 8 * levels;
                var series = new Series(chains);
                for (int c = 0; c < chains; c++)
                {
                    var target = new SplineTarget(new SplineBasis(3), fit.Model, new PoissonPrior(fit.PriorMean), fit.X, fit.Y);
                    var chain = new ReversibleJumpChain<KnotConfig>(Moves(), target, new KnotConfig([]), new Xoshiro256PlusPlus(5000 + c));
                    for (int i = 0; i < Burn; i++) chain.Step();
                    var configs = new KnotConfig[Measured];
                    var energy = new double[Measured];
                    var clock = Stopwatch.StartNew();
                    for (int i = 0; i < Measured; i++)
                    {
                        configs[i] = chain.Step();
                        energy[i] = chain.CurrentEvaluation.LogLikelihood(1.0);
                    }
                    series.Add(c, fit, new RunResult("", configs, energy, clock.Elapsed.TotalSeconds), 1);
                }
                ess.AppendLine(series.Row(fit.Name, "Untempered", double.NaN, levels, -1));
                File.WriteAllText(path, ess.ToString());
            }
    }

    private const string EssHeader = "data\tpath\tbeta_min\tlevels\tgroup\tchains\tess_knots_bulk\tess_knots_tail\t" +
                                     "ess_energy_bulk\tess_peak_bulk\tess_peak_tail\trhat_knots\trhat_peak\tmean_knots\t" +
                                     "mean_peak\tseconds\twork\tess_knots_per_kwork\tess_peak_per_kwork\tess_knots_per_second\t" +
                                     "ess_peak_per_second\tbarrier\ttrips_per_scan";

    /// <summary>
    /// The OQ1 rerun with exact reference rungs (OQ2), pre-registered in the preparation brief. Adapted DEO ladders of
    /// 4 and 8 rungs, 8 ladders × 2 groups per arm, judged by the cold chain's ESS per unit work as
    /// <see cref="HotEnds_ByColdChainEss"/> is, with the arms' barriers and round trips beside it. Guarded by
    /// <c>TEMPER_BENCH_OQ2</c>; <c>TEMPER_BENCH_FITS</c> splits the fits across processes.
    /// </summary>
    [Fact]
    public void ReferenceRungs_ByColdChainEss()
    {
        if (Environment.GetEnvironmentVariable("TEMPER_BENCH_OQ2") is not string path) return;   // heavy; set TEMPER_BENCH_OQ2 to run

        const int chains = 8;
        double[] hotEnds = [0.05, 0.2, 0.35, 0.5];
        var ess = new StringBuilder(EssHeader + "\tarm\n");
        foreach (Fit fit in Fits("peak", "sine", "counts", "oscillatory"))
        {
            bool normal = fit.Model is WeightedNormalModel;
            var arms = new List<(string Name, TemperingPath Path, IObservationModel Model, bool Reference, double BetaMin)>();
            IObservationModel linear = normal ? new WeightedNormalModel(form: TemperedForm.Linear) : fit.Model;
            foreach (double b in hotEnds) arms.Add((normal ? "L-int" : "L", TemperingPath.Likelihood, fit.Model, false, b));
            if (normal) foreach (double b in hotEnds) arms.Add(("L-lin", TemperingPath.Likelihood, linear, false, b));
            arms.Add((normal ? "L-lin-ref" : "L-ref", TemperingPath.Likelihood, linear, true, 0.05));
            arms.Add(("P", TemperingPath.Powered, fit.Model, false, 0.05));
            arms.Add(("P-ref", TemperingPath.Powered, fit.Model, true, 0.05));
            foreach ((string name, TemperingPath tp, IObservationModel model, bool reference, double betaMin) in arms)
                foreach (int levels in new[] { 4, 8 })
                    for (int group = 0; group < 2; group++)
                    {
                        var series = new Series(chains);
                        for (int c = 0; c < chains; c++)
                        {
                            RunResult r = Run(fit with { Model = model }, tp, SwapSchedule.DeterministicEvenOdd, adapt: true, levels,
                                              betaMin, seed: 1000 * (group + 1) + c, reference);
                            series.Add(c, fit, r, levels);
                        }
                        ess.AppendLine(series.Row(fit.Name, tp.ToString(), reference ? 0.0 : betaMin, levels, group) + "\t" + name);
                        File.WriteAllText(path, ess.ToString());   // partial results survive an interrupted run
                    }
        }
    }

    // Work is replica-steps of the measured phase, a reference rung's scan counted at its draw's attempts; NaN takes
    // levels × scans.
    private sealed record RunResult(string Row, KnotConfig[] Configs, double[] Energy, double Seconds, double Barrier = double.NaN,
                                    double TripsPerScan = double.NaN, double Work = double.NaN);

    // The cold chains' series, pooled into one set of rank-normalized diagnostics with the work they took.
    private sealed class Series(int chains)
    {
        private readonly double[][] _knots = new double[chains][], _energy = new double[chains][], _peak = new double[chains][];
        private double _seconds, _work, _barrier, _trips;

        public void Add(int c, Fit fit, RunResult r, int levels)
        {
            _knots[c] = r.Configs.Select(k => (double)k.Count).ToArray();
            _energy[c] = r.Energy;
            _peak[c] = PeakLocations(fit, r.Configs);
            _seconds += r.Seconds;
            _work += double.IsNaN(r.Work) ? (double)levels * r.Configs.Length : r.Work;
            _barrier += r.Barrier / chains;
            _trips += r.TripsPerScan / chains;
        }

        public string Row(string name, string path, double betaMin, int levels, int group)
        {
            ConvergenceDiagnostic k = ChainDiagnostics.Diagnose(_knots), e = ChainDiagnostics.Diagnose(_energy),
                                  p = ChainDiagnostics.Diagnose(_peak);
            string G(double v) => v.ToString("G5", CultureInfo.InvariantCulture);
            double kilo = _work / 1000.0;
            return string.Join('\t', name, path, G(betaMin), levels, group, _knots.Length, G(k.BulkEss.Value), G(k.TailEss.Value),
                               G(e.BulkEss.Value), G(p.BulkEss.Value), G(p.TailEss.Value), G(k.RHat.Value), G(p.RHat.Value),
                               G(_knots.SelectMany(s => s).Average()), G(_peak.SelectMany(s => s).Average()), G(_seconds), G(_work),
                               G(k.Ess.Value / kilo), G(p.Ess.Value / kilo), G(k.Ess.Value / _seconds), G(p.Ess.Value / _seconds),
                               G(_barrier), G(_trips));
        }
    }

    // The argmax of each draw's conditional-mean curve, computed after the timed loop; repeated configurations reuse it.
    private static double[] PeakLocations(Fit fit, KnotConfig[] configs)
    {
        var basis = new SplineBasis(3);
        var peaks = new double[configs.Length];
        KnotConfig? last = null;
        double lastPeak = double.NaN;
        for (int i = 0; i < configs.Length; i++)
        {
            if (!ReferenceEquals(configs[i], last))
            {
                double[] coef = fit.Model.PosteriorMeanCoefficients(basis.Design(configs[i], fit.X), fit.Y, null);
                lastPeak = SplineExtrema.Argmax(configs[i], coef, basis).Location;
                last = configs[i];
            }
            peaks[i] = lastPeak;
        }
        return peaks;
    }

    private static RunResult Run(Fit fit, TemperingPath tp, SwapSchedule schedule, bool adapt, int levels, double betaMin, int seed,
                                 bool reference = false)
    {
        var target = new SplineTarget(new SplineBasis(3), fit.Model, new PoissonPrior(fit.PriorMean), fit.X, fit.Y);
        double[] betas = reference
            ? ParallelTempering<KnotConfig>.ReferenceLadder(levels, betaMin)
            : ParallelTempering<KnotConfig>.GeometricLadder(levels, betaMin);
        IReferenceSampler<KnotConfig>? referenceSampler = reference ? target.Reference(tp) : null;
        var pt = new ParallelTempering<KnotConfig>(
            _ => Moves(), target, new KnotConfig([]), betas, tp, seed, schedule, referenceSampler);

        if (adapt)
            for (int round = 1; round <= 2048; round *= 2)
            {
                for (int i = 0; i < round; i++) pt.Step();
                pt.AdaptLadder();
            }
        else
            for (int i = 0; i < Burn; i++) pt.Step();

        long trips0 = pt.RoundTrips, scans0 = pt.Scans, swaps0 = pt.SwapAccepts, attempts0 = pt.SwapAttempts;
        long referenceAttempts0 = referenceSampler?.Attempts ?? 0;
        var before = new double[levels - 1];
        for (int i = 0; i < before.Length; i++) before[i] = pt.PairStats[i].MeanRejection * scans0;
        double coldKnots = 0.0, hotKnots = 0.0;
        var configs = new KnotConfig[Measured];
        var energySeries = new double[Measured];
        var clock = Stopwatch.StartNew();
        for (int i = 0; i < Measured; i++)
        {
            pt.Step();
            coldKnots += pt.Replicas[0].Current.Count;
            hotKnots += pt.Replicas[levels - 1].Current.Count;
            configs[i] = pt.Cold;
            energySeries[i] = pt.ColdChain.CurrentEvaluation.LogLikelihood(1.0);
        }
        double seconds = clock.Elapsed.TotalSeconds;

        var rejection = new double[levels - 1];
        for (int i = 0; i < rejection.Length; i++) rejection[i] = (pt.PairStats[i].MeanRejection * pt.Scans - before[i]) / Measured;
        double trips = pt.RoundTrips - trips0;
        double ideal = schedule == SwapSchedule.Sequential ? double.NaN : TemperingLadder.IdealRoundTripRate(rejection, schedule);
        string F(double v) => v.ToString("G5", CultureInfo.InvariantCulture);
        string row = string.Join('\t', fit.Name, tp, schedule, adapt, levels, seed, F(trips / Measured), F(trips / Measured / levels),
                           F(trips / seconds), F(TemperingLadder.Barrier(rejection)), F(ideal), F(coldKnots / Measured),
                           F(hotKnots / Measured), F((double)(pt.SwapAccepts - swaps0) / (pt.SwapAttempts - attempts0)),
                           string.Join(',', Array.ConvertAll(pt.Betas, b => b.ToString("F4", CultureInfo.InvariantCulture))));
        double work = referenceSampler is null
            ? (double)levels * Measured
            : (double)(levels - 1) * Measured + (referenceSampler.Attempts - referenceAttempts0);
        return new RunResult(row, configs, energySeries, seconds, TemperingLadder.Barrier(rejection), trips / Measured, work);
    }

    // The refactor probe's fit: a sine with 0.15 noise and an outlier every 17th point.
    private static (double[] X, double[] Y) Sine()
    {
        var rng = new Xoshiro256PlusPlus(seed: 77);
        int n = 120;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Gaussian(rng) + (i % 17 == 0 ? 2.0 : 0.0);
        }
        return (x, y);
    }

    // TemperingTests' peak: a narrow bump at 0.4 under 0.05 noise.
    private static (double[] X, double[] Y) Peak()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        int n = 150;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-120.0 * (x[i] - 0.4) * (x[i] - 0.4)) + 0.05 * Gaussian(rng);
        }
        return (x, y);
    }

    // The refactor probe's counts: rounded exp(1 + sin 2πx) plus a squared-noise bump.
    private static (double[] X, double[] Y) Counts()
    {
        var rng = new Xoshiro256PlusPlus(seed: 77);
        int n = 120;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double z = Gaussian(rng);
            y[i] = Math.Round(Math.Exp(1.0 + Math.Sin(2.0 * Math.PI * x[i])) + 0.5 * z * z);
        }
        return (x, y);
    }

    // BirthBenchmark's oscillation: sin 12πx under 0.15 noise, n = 200, about 13 knots; the slowest-mixing fit there.
    private static (double[] X, double[] Y) Oscillatory()
    {
        var rng = new Xoshiro256PlusPlus(seed: 61);
        int n = 200;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(12.0 * Math.PI * x[i]) + 0.15 * Gaussian(rng);
        }
        return (x, y);
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
