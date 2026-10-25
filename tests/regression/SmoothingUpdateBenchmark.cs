using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The anisotropic smoothing variances' two updates, Newton (with its random-walk escape steps) and slice, compared by the
/// worst component's bulk ESS of its log smoothing per second of wall time, over five regimes: BAPS on an 18 × 18
/// surface under the inverse-gamma prior (R1) and a scaled Weibull prior (R2), BAPS on a 30 × 30 surface with 196
/// coefficients (R3), <see cref="HierarchicalPSpline"/> on three replicate surfaces (R4), and R1 at five times the noise,
/// where the variances are weakly identified (R5). Each arm runs 4 chains of 1000 burn-in and 2000 draws per seed.
/// Pre-registered rule: slice becomes the default only if its ratio to Newton is at least 1.2 in some regime and at least
/// 0.8 in every one. Guarded by <c>SMOOTHING_BENCH</c> (seeds per arm); the table and summary go under
/// <c>artifacts/test-runs/SmoothingUpdateBenchmark/</c>.
/// </summary>
public sealed class SmoothingUpdateBenchmark
{
    private readonly ITestOutputHelper _out;
    public SmoothingUpdateBenchmark(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (TensorDesign Td, double[][] Responses) Surface(int grid, int knots, double noise, int replicates, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var xs = new double[grid * grid];
        var ys = new double[grid * grid];
        for (int a = 0; a < grid; a++)
            for (int b = 0; b < grid; b++) { xs[a * grid + b] = (a + 0.5) / grid; ys[a * grid + b] = (b + 0.5) / grid; }
        var responses = new double[replicates][];
        for (int r = 0; r < replicates; r++)
            responses[r] = xs.Select((x, i) => Math.Sin(1.5 * Math.PI * x) * Math.Sin(5.0 * Math.PI * ys[i]) + 0.05 * r + noise * Gaussian(rng))
                             .ToArray();
        var config = new KnotConfig(Enumerable.Range(1, knots).Select(k => k / (knots + 1.0)).ToArray());
        return (TensorProductDesign.Build(new SplineBasis(3), config, xs, new SplineBasis(3), config, ys), responses);
    }

    private static double WorstEss(IReadOnlyList<ConvergenceDiagnostic> components) => components.Min(c => c.BulkEss.Value);

    [Fact]
    public void Benchmark()
    {
        if (!int.TryParse(Environment.GetEnvironmentVariable("SMOOTHING_BENCH"), out int seeds) || seeds <= 0) return;   // heavy

        var regimes = new List<(string Name, Func<SmoothingVarianceUpdate, int, double> EssPerSecond)>();
        (TensorDesign r1, double[][] y1) = Surface(18, 6, 0.1, 1, seed: 31);
        (TensorDesign r3, double[][] y3) = Surface(30, 10, 0.1, 1, seed: 33);
        (TensorDesign r4, double[][] y4) = Surface(16, 5, 0.1, 3, seed: 41);
        (TensorDesign r5, double[][] y5) = Surface(18, 6, 0.5, 1, seed: 35);
        var p1 = new TensorPenalty(r1.NuX, r1.NuY);
        WeibullVariancePrior weibull = WeibullVariancePrior.ScaledTo(r1.Design, p1, bound: 3.0);

        double Baps(TensorDesign td, double[] y, ISmoothingVariancePrior? prior, SmoothingVarianceUpdate update, int seed)
        {
            var model = new Baps(td.Design, y, new TensorPenalty(td.NuX, td.NuY), componentPrior: prior, componentUpdate: update);
            var clock = Stopwatch.StartNew();
            BapsResult r = model.Run(BapsScheme.AnisotropicGibbs, chains: 4, masterSeed: seed, burn: 1000, samples: 2000);
            return WorstEss(r.ComponentLogLambdaConvergence!) / clock.Elapsed.TotalSeconds;
        }
        double Hierarchical(SmoothingVarianceUpdate update, int seed)
        {
            var model = new HierarchicalPSpline(r4.Design, y4, new TensorPenalty(r4.NuX, r4.NuY), anisotropicSmoothing: true,
                                                componentUpdate: update);
            var clock = Stopwatch.StartNew();
            HierarchicalResult h = model.Run(chains: 4, masterSeed: seed, burn: 1000, samples: 2000);
            return WorstEss(h.ComponentSmoothingVarianceConvergence!) / clock.Elapsed.TotalSeconds;
        }
        regimes.Add(("R1 BAPS IG", (u, s) => Baps(r1, y1[0], null, u, s)));
        regimes.Add(("R2 BAPS Weibull", (u, s) => Baps(r1, y1[0], weibull, u, s)));
        regimes.Add(("R3 BAPS nu=196", (u, s) => Baps(r3, y3[0], null, u, s)));
        regimes.Add(("R4 HPS centred", Hierarchical));
        regimes.Add(("R5 BAPS noisy", (u, s) => Baps(r5, y5[0], null, u, s)));

        static string F(double v) => v.ToString("G6", CultureInfo.InvariantCulture);
        var table = new StringBuilder("regime\tupdate\tseed\tworst_ess_per_s\n");
        var summary = new StringBuilder("regime\tnewton_ess_per_s\tslice_ess_per_s\tratio\n");
        var ratios = new List<double>();
        foreach ((string name, Func<SmoothingVarianceUpdate, int, double> essPerSecond) in regimes)
        {
            essPerSecond(SmoothingVarianceUpdate.Newton, 999);   // warm-up
            var means = new Dictionary<SmoothingVarianceUpdate, double>();
            foreach (SmoothingVarianceUpdate update in new[] { SmoothingVarianceUpdate.Newton, SmoothingVarianceUpdate.Slice })
            {
                var values = new List<double>();
                for (int s = 1; s <= seeds; s++)
                {
                    double v = essPerSecond(update, s);
                    values.Add(v);
                    table.Append(name).Append('\t').Append(update).Append('\t').Append(s).Append('\t').Append(F(v)).Append('\n');
                }
                means[update] = values.Average();
            }
            double ratio = means[SmoothingVarianceUpdate.Slice] / means[SmoothingVarianceUpdate.Newton];
            ratios.Add(ratio);
            summary.Append(name).Append('\t').Append(F(means[SmoothingVarianceUpdate.Newton])).Append('\t')
                   .Append(F(means[SmoothingVarianceUpdate.Slice])).Append('\t').Append(F(ratio)).Append('\n');
            _out.WriteLine($"[smoothing-bench] {name,-16} Newton {means[SmoothingVarianceUpdate.Newton],8:F1}  slice {means[SmoothingVarianceUpdate.Slice],8:F1}  ratio {ratio:F2}");
        }
        bool switchDefault = ratios.Max() >= 1.2 && ratios.Min() >= 0.8;
        _out.WriteLine($"[smoothing-bench] rule: slice {(switchDefault ? "becomes" : "does not become")} the default (ratios {ratios.Min():F2}..{ratios.Max():F2})");

        ArtifactRun run = HarnessArtifacts.Create("test-runs", "SmoothingUpdateBenchmark", "Maths.Regression.Tests",
                                                  new Dictionary<string, object?> { ["Seeds"] = seeds });
        _out.WriteLine($"[smoothing-bench] table {run.WriteText("smoothing-benchmark.tsv", table.ToString())}");
        _out.WriteLine($"[smoothing-bench] summary {run.WriteText("summary.tsv", summary.ToString())}");
    }
}
