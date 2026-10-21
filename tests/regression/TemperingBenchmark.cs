using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Round trips per unit work for BARS ladders (step 3.3, and Q5's revisit trigger): the swap schedule, ladder
/// adaptation and tempering path, measured by full round trips across the ladder (Syed et al. 2022) rather than by
/// swap acceptance, on three fits. Round trips compare like with like only between even-odd schedules on one path:
/// they are journeys between a ladder's own ends, which differ by path, and the sequential sweep's cascades let a
/// state cross several rungs with no local step between. So the cold chain's effective sample size per unit work (knot
/// count, and the cold state's log likelihood as an energy) is measured beside them, across the three seeds of each
/// configuration; at this budget it is the noisier but neutral measure. Guarded by the <c>TEMPER_BENCH</c> environment
/// variable, which names the TSV file the per-run table is written to (the ESS table goes beside it, <c>.ess.tsv</c>), so
/// the default suite stays fast.
/// </summary>
public sealed class TemperingBenchmark
{
    private const int Burn = 4095;      // 1 + 2 + … + 2048: Algorithm 4's rounds exactly
    private const int Measured = 20_000;
    private const int Thin = 10;        // ESS is O(n²) in the series length; the chains are far more correlated than this

    [Fact]
    public void RoundTrips_PerUnitWork()
    {
        if (Environment.GetEnvironmentVariable("TEMPER_BENCH") is not string path) return;   // heavy; set TEMPER_BENCH to run

        var sb = new StringBuilder();
        sb.AppendLine("data\tpath\tschedule\tadapted\tlevels\tseed\ttrips_per_scan\ttrips_per_replica_step\ttrips_per_second\t" +
                      "barrier\tideal_trips_per_scan\tcold_knots\thot_knots\tswap_rate\tladder");
        (double[] sx, double[] sy) = Sine();
        (double[] px, double[] py) = Peak();
        (double[] cx, double[] cy) = Counts();
        (string Name, IObservationModel Model, double[] X, double[] Y)[] fits =
            [("sine", new WeightedNormalModel(), sx, sy), ("peak", new WeightedNormalModel(), px, py), ("counts", new PoissonModel(), cx, cy)];
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
        var ess = new StringBuilder();
        ess.AppendLine("data\tpath\tschedule\tadapted\tlevels\tess_knots\tess_energy\tseconds\tess_knots_per_second\t" +
                       "ess_energy_per_second\tess_knots_per_1k_replica_steps");

        foreach ((string name, IObservationModel model, double[] x, double[] y) in fits)
            foreach ((TemperingPath tp, SwapSchedule schedule, bool adapt) in configs)
                foreach (int levels in new[] { 4, 8 })
                {
                    var knots = new double[3][];
                    var energy = new double[3][];
                    double seconds = 0.0;
                    for (int seed = 1; seed <= 3; seed++)
                    {
                        (string row, knots[seed - 1], energy[seed - 1], double sec) = Run(name, model, x, y, tp, schedule, adapt, levels, seed);
                        sb.AppendLine(row);
                        seconds += sec;
                    }
                    double essK = Maths.Samplers.Ensemble.ChainDiagnostics.Ess(knots) * Thin;
                    double essE = Maths.Samplers.Ensemble.ChainDiagnostics.Ess(energy) * Thin;
                    string G(double v) => v.ToString("G5", CultureInfo.InvariantCulture);
                    ess.AppendLine(string.Join('\t', name, tp, schedule, adapt, levels, G(essK), G(essE), G(seconds),
                                               G(essK / seconds), G(essE / seconds), G(essK / (3.0 * Measured * levels / 1000.0))));
                }
        File.WriteAllText(path, sb.ToString());
        File.WriteAllText(Path.ChangeExtension(path, ".ess.tsv"), ess.ToString());
    }

    private static (string Row, double[] Knots, double[] Energy, double Seconds) Run(
        string name, IObservationModel model, double[] x, double[] y, TemperingPath tp, SwapSchedule schedule, bool adapt,
        int levels, int seed)
    {
        var kernel = new LocalBetaKernel(50.0);
        var target = new SplineTarget(new SplineBasis(3), model, new PoissonPrior(5.0), x, y);
        var pt = new ParallelTempering<KnotConfig>(
            _ => [new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel)],
            target, new KnotConfig([]), ParallelTempering<KnotConfig>.GeometricLadder(levels, 0.05), tp, seed, schedule);

        if (adapt)
            for (int round = 1; round <= 2048; round *= 2)
            {
                for (int i = 0; i < round; i++) pt.Step();
                pt.AdaptLadder();
            }
        else
            for (int i = 0; i < Burn; i++) pt.Step();

        long trips0 = pt.RoundTrips, scans0 = pt.Scans, swaps0 = pt.SwapAccepts, attempts0 = pt.SwapAttempts;
        var before = new double[levels - 1];
        for (int i = 0; i < before.Length; i++) before[i] = pt.PairStats[i].MeanRejection * scans0;
        double coldKnots = 0.0, hotKnots = 0.0;
        var knotSeries = new double[Measured / Thin];
        var energySeries = new double[Measured / Thin];
        var clock = Stopwatch.StartNew();
        for (int i = 0; i < Measured; i++)
        {
            pt.Step();
            coldKnots += pt.Replicas[0].Current.Count;
            hotKnots += pt.Replicas[levels - 1].Current.Count;
            if (i % Thin == 0)
            {
                knotSeries[i / Thin] = pt.Cold.Count;
                energySeries[i / Thin] = pt.ColdChain.CurrentEvaluation.LogLikelihood(1.0);
            }
        }
        double seconds = clock.Elapsed.TotalSeconds;

        var rejection = new double[levels - 1];
        for (int i = 0; i < rejection.Length; i++) rejection[i] = (pt.PairStats[i].MeanRejection * pt.Scans - before[i]) / Measured;
        double trips = pt.RoundTrips - trips0;
        double ideal = schedule == SwapSchedule.Sequential ? double.NaN : TemperingLadder.IdealRoundTripRate(rejection, schedule);
        string F(double v) => v.ToString("G5", CultureInfo.InvariantCulture);
        string row = string.Join('\t', name, tp, schedule, adapt, levels, seed, F(trips / Measured), F(trips / Measured / levels),
                           F(trips / seconds), F(TemperingLadder.Barrier(rejection)), F(ideal), F(coldKnots / Measured),
                           F(hotKnots / Measured), F((double)(pt.SwapAccepts - swaps0) / (pt.SwapAttempts - attempts0)),
                           string.Join(',', Array.ConvertAll(pt.Betas, b => b.ToString("F4", CultureInfo.InvariantCulture))));
        return (row, knotSeries, energySeries, seconds);
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

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
