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
/// <c>HPS_MIXING_BENCH</c> (seeds per arm): how fast <see cref="HierarchicalPSpline"/>'s schemes mix, OQ18 phase C's
/// exploration. Five regimes, each a family of offset sine replicates: the 6-replicate base fit; near-identical replicates
/// (no offsets, so τ_u is near 0); few (R = 3) and many (R = 30) replicates; and weak data per replicate (12 points at
/// noise 0.5). Five arms: the centred sweep, the collapsed population block and the scale interweaving under the
/// inverse-gamma prior, and the parameter-expanded sweep under the half-t prior (ν = 1, scale the responses' sd), which
/// changes the posterior as well as the sampler; and the collapsed block with the interweaving in one sweep. Each run is 4 chains × 4,000 draws after 1,000 burn-in sweeps. The
/// readouts are bulk ESS per 1,000 draws of log τ_u², log σ², log τ_0² and the least favourable population coefficient,
/// τ_u²'s ESS per second, the worst R̂, and τ_u's posterior median. Tables go under
/// <c>artifacts/test-runs/HierarchicalMixingBenchmark/</c>. Exploration on single data sets, not a study for defaults.
/// </summary>
public sealed class HierarchicalMixingBenchmark
{
    private readonly ITestOutputHelper _out;
    public HierarchicalMixingBenchmark(ITestOutputHelper output) => _out = output;

    private static string F(double v) => v.ToString("G4", CultureInfo.InvariantCulture);

    private sealed record Regime(string Name, int Replicates, int Points, double Noise, double Spread);

    private static readonly Regime[] Regimes =
    {
        new("base", 6, 40, 0.2, 0.3),
        new("near-identical", 6, 40, 0.2, 0.0),
        new("few (R=3)", 3, 40, 0.2, 0.3),
        new("many (R=30)", 30, 40, 0.2, 0.3),
        new("weak data", 6, 12, 0.5, 0.3),
    };

    private static (double[,] Z, double[][] Y) Data(Regime regime)
    {
        var rng = new Xoshiro256PlusPlus(seed: 41);
        int n = regime.Points;
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = (i + 0.5) / n;
        var y = new double[regime.Replicates][];
        for (int r = 0; r < regime.Replicates; r++)
        {
            y[r] = new double[n];
            double shift = regime.Spread * (rng.NextDouble() - 0.5);
            for (int i = 0; i < n; i++)
            {
                double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
                y[r][i] = Math.Sin(2.0 * Math.PI * x[i]) + shift + regime.Noise * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            }
        }
        var knots = new double[10];
        for (int k = 0; k < knots.Length; k++) knots[k] = (k + 1.0) / (knots.Length + 1);
        return (new SplineBasis(3).Design(new KnotConfig(knots), x), y);
    }

    private static double Median(IEnumerable<double> values)
    {
        double[] a = values.Order().ToArray();
        return a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
    }

    [Fact]
    public void Schemes_AcrossRegimes()
    {
        if (!int.TryParse(Environment.GetEnvironmentVariable("HPS_MIXING_BENCH"), out int seeds) || seeds <= 0) return;   // heavy
        const int chains = 4, burn = 1000, samples = 4000;
        var arms = new (string Name, HierarchicalScheme Scheme)[]
        {
            ("centred", HierarchicalScheme.Centred), ("collapsed", HierarchicalScheme.CollapsedPopulation),
            ("interweaved", HierarchicalScheme.Interweaved), ("collapsed+interweaved", HierarchicalScheme.CollapsedInterweaved),
            ("half-t PX", HierarchicalScheme.ParameterExpanded),
        };
        var table = new StringBuilder("regime\tarm\tseed\ttauU2_ess_per_1000\tsigma2_ess_per_1000\ttau02_ess_per_1000\tbeta0_min_ess_per_1000\ttauU2_ess_per_second\tworst_rhat\ttauU_median\tseconds\n");
        var summary = new StringBuilder("regime\tarm\ttauU2_ess_per_1000\tsigma2_ess_per_1000\ttau02_ess_per_1000\tbeta0_min_ess_per_1000\ttauU2_ess_per_second\tworst_rhat\ttauU_median\n");
        foreach (Regime regime in Regimes)
        {
            (double[,] z, double[][] y) = Data(regime);
            int nu = z.GetLength(1);
            foreach ((string name, HierarchicalScheme scheme) in arms)
            {
                var sampler = scheme == HierarchicalScheme.ParameterExpanded
                    ? new HierarchicalPSpline(z, y, new DifferencePenalty(2), deviationPrior: new HalfTDeviationPrior())
                    : new HierarchicalPSpline(z, y, new DifferencePenalty(2));
                var rows = new List<double[]>();
                for (int s = 0; s < seeds; s++)
                {
                    var clock = Stopwatch.StartNew();
                    HierarchicalResult h = sampler.Run(chains: chains, masterSeed: 100 + s, burn: burn, samples: samples, scheme: scheme);
                    double seconds = clock.Elapsed.TotalSeconds;
                    double perThousand = 1000.0 / h.Draws;
                    double beta0 = double.PositiveInfinity;
                    for (int j = 0; j < nu; j++)
                    {
                        int jj = j;
                        double[][] traces = h.PosteriorDraws.Select(c => c.Select(d => d.Population[jj]).ToArray()).ToArray();
                        beta0 = Math.Min(beta0, ChainDiagnostics.Diagnose(traces).BulkEss.Value);
                    }
                    double worst = new[] { h.NoiseVarianceConvergence.RHat.Value, h.DeviationVarianceConvergence.RHat.Value,
                                           h.PopulationSmoothingVarianceConvergence.RHat.Value, h.RHatPopulationMax.Value }.Max();
                    double tauU = Math.Sqrt(h.DeviationVarianceInterval.Median);
                    double[] row =
                    {
                        h.DeviationVarianceConvergence.BulkEss.Value * perThousand, h.NoiseVarianceConvergence.BulkEss.Value * perThousand,
                        h.PopulationSmoothingVarianceConvergence.BulkEss.Value * perThousand, beta0 * perThousand,
                        h.DeviationVarianceConvergence.BulkEss.Value / seconds, worst, tauU, seconds,
                    };
                    rows.Add(row);
                    table.Append($"{regime.Name}\t{name}\t{s}\t{string.Join('\t', row.Select(F))}\n");
                }
                double[] med = Enumerable.Range(0, 7).Select(k => Median(rows.Select(r => r[k]))).ToArray();
                summary.Append($"{regime.Name}\t{name}\t{string.Join('\t', med.Select(F))}\n");
                _out.WriteLine($"[hps-mixing] {regime.Name,-15} {name,-21}: τ_u² ESS/1000 {F(med[0]),6}  σ² {F(med[1]),6}  τ_0² {F(med[2]),6}  " +
                               $"β_0 min {F(med[3]),6}  τ_u² ESS/s {F(med[4]),7}  worst R̂ {F(med[5])}  τ_u median {F(med[6])}");
            }
        }
        ArtifactRun run = HarnessArtifacts.Create("test-runs", "HierarchicalMixingBenchmark", "Maths.Regression.Tests",
                                                  new Dictionary<string, object?> { ["Seeds"] = seeds, ["Chains"] = chains, ["Burn"] = burn, ["Samples"] = samples });
        _out.WriteLine($"[hps-mixing] table {run.WriteText("schemes.tsv", table.ToString())}");
        _out.WriteLine($"[hps-mixing] summary {run.WriteText("schemes-summary.tsv", summary.ToString())}");
    }
}
