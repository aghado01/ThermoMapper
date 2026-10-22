using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The exploratory share of births (step 3.9, OQ11): knot-count and peak-location ESS per second for ε from 0 (local
/// births only) to 1 (births from the base only), with the uniform base on every fit and the design base on the two uneven designs. The fits' MSE and
/// mean knot count check that the posterior does not move with ε. Guarded by <c>BIRTH_BENCH</c> (seeds per arm);
/// <c>BIRTH_BENCH_FITS</c> (comma-separated names) restricts the fits and <c>BIRTH_BENCH_EPS</c> (comma-separated values)
/// the exploratory shares. The table and summary go under
/// <c>artifacts/test-runs/BirthBenchmark/</c>.
/// </summary>
public sealed class BirthBenchmark
{
    private readonly ITestOutputHelper _out;
    public BirthBenchmark(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    internal sealed record Fit(string Name, double[] X, double[] Y, Func<double, double> Truth, double PriorMean, bool Uneven);

    private static Fit Make(string name, double[] x, Func<double, double> truth, double noise, int seed, double priorMean, bool uneven,
                            Func<int, double>? extra = null)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] y = x.Select((t, i) => truth(t) + noise * Gaussian(rng) + (extra?.Invoke(i) ?? 0.0)).ToArray();
        return new Fit(name, x, y, truth, priorMean, uneven);
    }

    // Shared with MixingBenchmark's relocation study.
    internal static Fit[] Fits()
    {
        double[] Even(int n) => Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
        double[] gap = Enumerable.Range(0, 60).Select(i => 0.35 * (i + 0.5) / 60)
            .Concat(Enumerable.Range(0, 60).Select(i => 0.65 + 0.35 * (i + 0.5) / 60)).ToArray();
        double[] skewed = Enumerable.Range(0, 150).Select(i => Math.Pow((i + 0.5) / 150, 2.0)).ToArray();   // dense near 0
        Func<double, double> sine = t => Math.Sin(2.0 * Math.PI * t);
        Func<double, double> bumpSine = t => Math.Exp(-Math.Pow((t - 0.5) / 0.06, 2)) + Math.Sin(3.0 * t);
        return
        [
            Make("adaptive-stop", Even(100), sine, 0.15, 12, 4.0, false),                                   // EnsembleTests
            Make("sine-outliers", Even(120), sine, 0.15, 77, 5.0, false, i => i % 17 == 0 ? 2.0 : 0.0),        // TemperingBenchmark
            Make("peak", Even(150), t => Math.Exp(-120.0 * (t - 0.4) * (t - 0.4)), 0.05, 11, 5.0, false),    // TemperingBenchmark
            Make("gap", gap, bumpSine, 0.1, 37, 5.0, true),
            Make("skewed", skewed, bumpSine, 0.1, 41, 5.0, true),
            // A jump, which needs a cluster of nearly adjacent knots: the case DMGK's local births were designed for.
            Make("jump", Even(150), t => 0.5 * Math.Sin(2.0 * Math.PI * t) + (t > 0.5 ? 1.5 : 0.0), 0.1, 53, 5.0, false),
            // Many densely packed knots, about fifteen to twenty.
            Make("oscillatory", Even(200), t => Math.Sin(12.0 * Math.PI * t), 0.15, 61, 15.0, false),
        ];
    }

    [Fact]
    public void ExploratoryBirths_EssPerSecond()
    {
        if (!int.TryParse(Environment.GetEnvironmentVariable("BIRTH_BENCH"), out int seeds) || seeds <= 0) return;   // heavy
        static string F(double v) => v.ToString("R", CultureInfo.InvariantCulture);
        double[] epsilons = Environment.GetEnvironmentVariable("BIRTH_BENCH_EPS") is string shares
            ? shares.Split(',').Select(e => double.Parse(e, CultureInfo.InvariantCulture)).ToArray()
            : [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0];
        double[] grid = Enumerable.Range(0, 51).Select(i => i / 50.0).ToArray();
        var table = new StringBuilder("fit\tbase\texploration\tseed\tseconds\tknots_ess\tpeak_ess\tknots_ess_per_s\tpeak_ess_per_s\t" +
                                      "mean_knots\tmse\tacceptance\tthin_support_share\n");
        var rows = new List<(string Fit, BirthBase Base, double Eps, double Seconds, double KnotsEss, double PeakEss, double Knots, double Mse)>();

        string[]? only = Environment.GetEnvironmentVariable("BIRTH_BENCH_FITS")?.Split(',');
        foreach (Fit fit in Fits().Where(f => only is null || only.Contains(f.Name)))
            foreach (BirthBase birthBase in fit.Uneven ? new[] { BirthBase.Uniform, BirthBase.Design } : new[] { BirthBase.Uniform })
                foreach (double eps in epsilons)
                    for (int s = 0; s < seeds; s++)
                    {
                        var ensemble = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(fit.PriorMean),
                                                        new LocalBetaKernel(50.0), birthExploration: eps, birthBase: birthBase);
                        var clock = Stopwatch.StartNew();
                        BarsResult r = ensemble.Run(fit.X, fit.Y, grid: grid, chains: 4, masterSeed: 100 + s, burn: 2000, samples: 10_000);
                        double seconds = clock.Elapsed.TotalSeconds;
                        double mse = grid.Select((t, g) => Math.Pow(r.Fit[g] - fit.Truth(t), 2)).Where((_, g) =>
                            fit.Name != "gap" || grid[g] < 0.35 || grid[g] > 0.65).Average();
                        double knotsEss = r.Knots.BulkEss.Value, peakEss = r.Peak.LocationConvergence.BulkEss.Value;
                        double thin = (double)r.Evaluations.Count(DesignRejection.ThinSupport) / r.Evaluations.Total;
                        table.Append(fit.Name).Append('\t').Append(birthBase).Append('\t').Append(F(eps)).Append('\t').Append(s).Append('\t')
                             .Append(F(seconds)).Append('\t').Append(F(knotsEss)).Append('\t').Append(F(peakEss)).Append('\t')
                             .Append(F(knotsEss / seconds)).Append('\t').Append(F(peakEss / seconds)).Append('\t').Append(F(r.MeanKnots))
                             .Append('\t').Append(F(mse)).Append('\t').Append(F(r.AcceptanceRate)).Append('\t').Append(F(thin)).Append('\n');
                        rows.Add((fit.Name, birthBase, eps, seconds, knotsEss, peakEss, r.MeanKnots, mse));
                    }

        // Per fit, base and ε: medians over seeds.
        var summary = new StringBuilder("fit\tbase\texploration\tknots_ess_per_s\tpeak_ess_per_s\tmean_knots\tmse\n");
        static double Median(IEnumerable<double> v) { double[] a = v.Order().ToArray(); return a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]); }
        foreach (var g in rows.GroupBy(r => (r.Fit, r.Base, r.Eps)))
        {
            double kps = Median(g.Select(r => r.KnotsEss / r.Seconds)), pps = Median(g.Select(r => r.PeakEss / r.Seconds));
            double knots = Median(g.Select(r => r.Knots)), mse = Median(g.Select(r => r.Mse));
            summary.Append(g.Key.Fit).Append('\t').Append(g.Key.Base).Append('\t').Append(F(g.Key.Eps)).Append('\t').Append(F(kps)).Append('\t')
                   .Append(F(pps)).Append('\t').Append(F(knots)).Append('\t').Append(F(mse)).Append('\n');
            _out.WriteLine($"[birth-bench] {g.Key.Fit,-13} {g.Key.Base,-7} ε={g.Key.Eps:F2}: knots ESS/s {kps,7:F1}  peak ESS/s {pps,7:F1}  " +
                           $"mean k {knots:F2}  mse {mse:E2}");
        }

        ArtifactRun run = HarnessArtifacts.Create("test-runs", "BirthBenchmark", "Maths.Regression.Tests",
                                                  new Dictionary<string, object?> { ["Seeds"] = seeds });
        _out.WriteLine($"[birth-bench] table {run.WriteText("birth-benchmark.tsv", table.ToString())}");
        _out.WriteLine($"[birth-bench] summary {run.WriteText("summary.tsv", summary.ToString())}");
    }
}
