using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Two mixing studies on BARS's default kernel (the preparation brief's "After Phase 3: BARS mixing"); births come from
/// the uniform base throughout, so the local kernel drives relocations only.
/// <list type="bullet">
/// <item><c>RELOC_BENCH</c> (seeds per arm): the relocation kernel's reach, on <see cref="BirthBenchmark"/>'s fits. The
/// local Beta kernel fixed at τ = 10, 50 (the default) or 500, step SDs of about 0.15, 0.07 and 0.02; tuned per chain in
/// burn-in toward 44% relocation acceptance and then frozen (<c>adaptKernelTau</c> = 50); the uniform kernel, a global
/// relocation as in ABMRS; and an even mixture of the uniform and the default.</item>
/// <item><c>SCHEDULE_BENCH</c> (seeds per arm): the move schedule. Constant weights 0.4/0.4/0.2 (the default) and
/// 0.3/0.3/0.4 (He, Yang and Kang's); DMGK's prior-aware schedule at c = 0.4 and at Green's budget-scaled c; and a
/// constant schedule at DMGK's budget, masked where the prior has no mass. Three regimes: the prior at the posterior's
/// mode, a prior in conflict with the likelihood, and an active hard cap. The readouts are ESJD(k), the reversal rate
/// (dimension changes undone by the next one), and ESS of the knot count, the energy (log likelihood), the fit at
/// x = 0.5 and the peak location.</item>
/// </list>
/// <c>MIXING_BENCH_FITS</c> (comma-separated) restricts the relocation study's fits. Tables and summaries go under
/// <c>artifacts/test-runs/MixingBenchmark/</c>.
/// </summary>
public sealed class MixingBenchmark
{
    private readonly ITestOutputHelper _out;
    public MixingBenchmark(ITestOutputHelper output) => _out = output;

    private static string F(double v) => v.ToString("G5", CultureInfo.InvariantCulture);

    private static double Median(IEnumerable<double> values) => Quantile(values, 0.5);

    private static double Quantile(IEnumerable<double> values, double q)
    {
        double[] a = values.Order().ToArray();
        double h = q * (a.Length - 1);
        int lo = (int)Math.Floor(h);
        return lo + 1 < a.Length ? a[lo] + (h - lo) * (a[lo + 1] - a[lo]) : a[lo];
    }

    [Fact]
    public void Relocation_KernelSpread()
    {
        if (!int.TryParse(Environment.GetEnvironmentVariable("RELOC_BENCH"), out int seeds) || seeds <= 0) return;   // heavy
        (string Name, Func<IKnotKernel> Kernel, double? Adapt)[] arms =
        [
            ("fixed-10", () => new LocalBetaKernel(10.0), null), ("fixed-50", () => new LocalBetaKernel(50.0), null),
            ("fixed-500", () => new LocalBetaKernel(500.0), null), ("adaptive-50", () => new LocalBetaKernel(50.0), 50.0),
            ("uniform", () => new UniformKernel(), null), ("mixture-0.5", () => new MixtureKernel(0.5, new LocalBetaKernel(50.0)), null),
        ];
        double[] grid = Enumerable.Range(0, 51).Select(i => i / 50.0).ToArray();
        string[]? only = Environment.GetEnvironmentVariable("MIXING_BENCH_FITS")?.Split(',');
        var table = new StringBuilder("fit\tarm\tseed\tseconds\tknots_ess\tpeak_ess\tmean_knots\tmse\trelocate_acceptance\tacceptance\n");
        var rows = new List<(string Fit, string Arm, double KnotsEss, double PeakEss, double Knots, double Mse, double Relocate, double Seconds)>();
        foreach (BirthBenchmark.Fit fit in BirthBenchmark.Fits().Where(f => only is null || only.Contains(f.Name)))
            foreach ((string arm, Func<IKnotKernel> kernel, double? adapt) in arms)
                for (int s = 0; s < seeds; s++)
                {
                    var ensemble = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(fit.PriorMean),
                                                    kernel(), adaptKernelTau: adapt);
                    var clock = Stopwatch.StartNew();
                    BarsResult r = ensemble.Run(fit.X, fit.Y, grid: grid, chains: 4, masterSeed: 100 + s, burn: 2000, samples: 10_000);
                    double seconds = clock.Elapsed.TotalSeconds;
                    double mse = grid.Select((t, g) => Math.Pow(r.Fit[g] - fit.Truth(t), 2)).Where((_, g) =>
                        fit.Name != "gap" || grid[g] < 0.35 || grid[g] > 0.65).Average();
                    MoveStat relocate = r.MoveStats.First(m => m.Key == "relocate");
                    double relocateRate = (double)relocate.Accepted / relocate.Attempts;
                    double knotsEss = r.Knots.BulkEss.Value, peakEss = r.Peak.LocationConvergence.BulkEss.Value;
                    table.Append(string.Join('\t', fit.Name, arm, s, F(seconds), F(knotsEss), F(peakEss), F(r.MeanKnots), F(mse),
                                             F(relocateRate), F(r.AcceptanceRate))).Append('\n');
                    rows.Add((fit.Name, arm, knotsEss, peakEss, r.MeanKnots, mse, relocateRate, seconds));
                }

        var summary = new StringBuilder("fit\tarm\tknots_ess\tpeak_ess\tknots_ess_vs_default\tpeak_ess_vs_default\tmean_knots\tmse\trelocate_acceptance\tseconds\n");
        foreach (var g in rows.GroupBy(r => r.Fit))
        {
            double k0 = Median(g.Where(r => r.Arm == "fixed-50").Select(r => r.KnotsEss));
            double p0 = Median(g.Where(r => r.Arm == "fixed-50").Select(r => r.PeakEss));
            foreach (var a in g.GroupBy(r => r.Arm))
            {
                double k = Median(a.Select(r => r.KnotsEss)), p = Median(a.Select(r => r.PeakEss));
                summary.Append(string.Join('\t', g.Key, a.Key, F(k), F(p), F(k / k0), F(p / p0), F(Median(a.Select(r => r.Knots))),
                                           F(Median(a.Select(r => r.Mse))), F(Median(a.Select(r => r.Relocate))), F(Median(a.Select(r => r.Seconds))))).Append('\n');
                _out.WriteLine($"[reloc] {g.Key,-13} {a.Key,-11}: knots ESS {k,7:F0} ({k / k0:F2}×)  peak ESS {p,7:F0} ({p / p0:F2}×)  " +
                               $"relocate acc {Median(a.Select(r => r.Relocate)):F3}  mean k {Median(a.Select(r => r.Knots)):F2}");
            }
        }
        ArtifactRun run = HarnessArtifacts.Create("test-runs", "MixingBenchmark", "Maths.Regression.Tests",
                                                  new Dictionary<string, object?> { ["Study"] = "relocation", ["Seeds"] = seeds });
        _out.WriteLine($"[reloc] table {run.WriteText("relocation.tsv", table.ToString())}");
        _out.WriteLine($"[reloc] summary {run.WriteText("relocation-summary.tsv", summary.ToString())}");
    }

    // Relocation with reach: from the uniform with probability share, otherwise from the local kernel.
    private sealed class MixtureKernel(double share, IKnotKernel local) : IKnotKernel
    {
        public double Sample(double center, Xoshiro256PlusPlus rng) => rng.NextDouble() < share ? rng.NextDouble() : local.Sample(center, rng);

        public double LogDensity(double x, double center)
        {
            double a = Math.Log(share), b = Math.Log(1.0 - share) + local.LogDensity(x, center);
            double hi = Math.Max(a, b);
            return hi + Math.Log(Math.Exp(a - hi) + Math.Exp(b - hi));
        }
    }

    private sealed record Regime(string Name, double[] X, double[] Y, IComplexityPrior Prior);

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static Regime[] Regimes()
    {
        double[] Even(int n) => Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
        double[] Noisy(double[] x, Func<double, double> f, int seed)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            return x.Select(t => f(t) + 0.15 * Gaussian(rng)).ToArray();
        }
        double[] sx = Even(100), sy = Noisy(sx, t => Math.Sin(2.0 * Math.PI * t), 12);   // BirthBenchmark's adaptive-stop sine
        double[] wx = Even(150), wy = Noisy(wx, t => Math.Sin(7.0 * Math.PI * t), 71);   // wants about 8 knots
        return
        [
            new("mode", sx, sy, new PoissonPrior(2.0)),       // the posterior's 1.8 knots sit at the prior's mode
            new("conflict", wx, wy, new PoissonPrior(1.0)),   // the likelihood wants about 8 knots, the prior about 1
            new("cap", wx, wy, new UniformPrior(8)),          // the same data, capped where they want to be
        ];
    }

    // The schedule arms, each building one chain's moves; births come from the uniform base.
    private static (string Name, Func<IComplexityPrior, IRjMove<KnotConfig>[]> Moves)[] ScheduleArms()
    {
        static (KnotBirthDensity Density, LocalBetaKernel Kernel) Parts()
        {
            var kernel = new LocalBetaKernel(50.0);
            return (new KnotBirthDensity(kernel, exploration: 1.0), kernel);
        }
        static IRjMove<KnotConfig>[] Constant(double b, double d, double r)
        {
            (KnotBirthDensity density, LocalBetaKernel kernel) = Parts();
            return [new KnotBirthMove(density, b), new KnotDeathMove(density, d), new KnotRelocateMove(kernel, r)];
        }
        static IRjMove<KnotConfig>[] Dmgk(IComplexityPrior prior, double c)
        {
            (KnotBirthDensity density, LocalBetaKernel kernel) = Parts();
            return DmgkSchedule.Wrap(new KnotBirthMove(density), new KnotDeathMove(density), new KnotRelocateMove(kernel), prior, c);
        }
        // Constant birth and death at half DMGK(0.4)'s largest b_k + d_k, with no mass on moves into states the prior
        // excludes; relocation takes the rest.
        static IRjMove<KnotConfig>[] Masked(IComplexityPrior prior)
        {
            IRjMove<KnotConfig>[] dmgk = Dmgk(prior, 0.4);
            double budget = 0.0;
            for (int k = 0; k <= 40; k++)
            {
                var state = new KnotConfig(Enumerable.Range(0, k).Select(i => (i + 0.5) / (k + 1)).ToArray());
                if (!double.IsNegativeInfinity(prior.LogPrior(k))) budget = Math.Max(budget, dmgk[0].Weight(state) + dmgk[1].Weight(state));
            }
            double half = 0.5 * budget;
            (KnotBirthDensity density, LocalBetaKernel kernel) = Parts();
            double Birth(KnotConfig s) => double.IsNegativeInfinity(prior.LogPrior(s.Count + 1)) ? 0.0 : half;
            double Death(KnotConfig s) => s.Count == 0 ? 0.0 : half;
            return
            [
                new ScheduledMove<KnotConfig>(new KnotBirthMove(density), Birth),
                new ScheduledMove<KnotConfig>(new KnotDeathMove(density), Death),
                new ScheduledMove<KnotConfig>(new KnotRelocateMove(kernel), s => 1.0 - Birth(s) - Death(s)),
            ];
        }
        return
        [
            ("const-0.4", _ => Constant(0.4, 0.4, 0.2)),
            ("const-0.3", _ => Constant(0.3, 0.3, 0.4)),
            ("dmgk-0.4", p => Dmgk(p, 0.4)),
            ("dmgk-green", p => Dmgk(p, DmgkSchedule.GreenConstant(p, maxKnots: 40))),
            ("masked-equal", Masked),
        ];
    }

    [Fact]
    public void MoveSchedules_ByRegime()
    {
        if (!int.TryParse(Environment.GetEnvironmentVariable("SCHEDULE_BENCH"), out int seeds) || seeds <= 0) return;   // heavy
        const int chains = 4, burn = 2000, measured = 10_000;
        var basis = new SplineBasis(3);
        var model = new WeightedNormalModel();
        var table = new StringBuilder("regime\tarm\tseed\tseconds\tesjd_k\treversal_rate\tdimension_changes\tess_knots\t" +
                                      "ess_energy\tess_fit_mid\tess_peak\tmean_knots\tacceptance\n");
        var rows = new List<(string Regime, string Arm, double[] Values)>();
        Regime[] regimes = Regimes();
        var arms = ScheduleArms();
        int[] masters = SeedTree.Derive(2026_10_22, regimes.Length * arms.Length);
        for (int ri = 0; ri < regimes.Length; ri++)
            for (int ai = 0; ai < arms.Length; ai++)
            {
                Regime regime = regimes[ri];
                (string arm, Func<IComplexityPrior, IRjMove<KnotConfig>[]> moves) = arms[ai];
                int[] chainSeeds = SeedTree.Derive(masters[ri * arms.Length + ai], seeds * chains);
                for (int s = 0; s < seeds; s++)
                {
                    var knots = new double[chains][];
                    var energy = new double[chains][];
                    var fitMid = new double[chains][];
                    var peak = new double[chains][];
                    double seconds = 0.0, jumps = 0.0, accepted = 0.0;
                    long changes = 0, reversals = 0;
                    for (int c = 0; c < chains; c++)
                    {
                        var target = new SplineTarget(basis, model, regime.Prior, regime.X, regime.Y);
                        var chain = new ReversibleJumpChain<KnotConfig>(moves(regime.Prior), target, new KnotConfig([]),
                                                                        new Xoshiro256PlusPlus(chainSeeds[s * chains + c]));
                        for (int i = 0; i < burn; i++) chain.Step();
                        var configs = new KnotConfig[measured];
                        energy[c] = new double[measured];
                        KnotConfig previous = chain.Current;
                        int lastDirection = 0;
                        var clock = Stopwatch.StartNew();
                        for (int i = 0; i < measured; i++)
                        {
                            KnotConfig k = chain.Step();
                            configs[i] = k;
                            energy[c][i] = chain.CurrentEvaluation.LogLikelihood(1.0);
                            int delta = k.Count - previous.Count;
                            jumps += delta * delta;
                            if (!ReferenceEquals(k, previous)) accepted++;
                            if (delta != 0)
                            {
                                int direction = Math.Sign(delta);
                                if (lastDirection != 0 && direction == -lastDirection) reversals++;
                                if (lastDirection != 0) changes++;
                                lastDirection = direction;
                            }
                            previous = k;
                        }
                        seconds += clock.Elapsed.TotalSeconds;
                        knots[c] = configs.Select(k => (double)k.Count).ToArray();
                        (fitMid[c], peak[c]) = Readouts(basis, model, regime, configs);
                    }
                    ConvergenceDiagnostic dk = ChainDiagnostics.Diagnose(knots), de = ChainDiagnostics.Diagnose(energy),
                                          df = ChainDiagnostics.Diagnose(fitMid), dp = ChainDiagnostics.Diagnose(peak);
                    double[] values =
                    [
                        seconds, jumps / (chains * measured), changes == 0 ? double.NaN : (double)reversals / changes, changes,
                        dk.BulkEss.Value, de.BulkEss.Value, df.BulkEss.Value, dp.BulkEss.Value,
                        knots.SelectMany(v => v).Average(), accepted / (chains * measured),
                    ];
                    table.Append(regime.Name).Append('\t').Append(arm).Append('\t').Append(s).Append('\t')
                         .Append(string.Join('\t', values.Select(F))).Append('\n');
                    rows.Add((regime.Name, arm, values));
                }
            }

        string[] names = ["seconds", "esjd_k", "reversal_rate", "dimension_changes", "ess_knots", "ess_energy", "ess_fit_mid", "ess_peak", "mean_knots", "acceptance"];
        var summary = new StringBuilder("regime\tarm\t" + string.Join('\t', names.Select(n => n + "_median")) + "\tess_knots_q25\tess_knots_q75\tess_peak_q25\tess_peak_q75\n");
        foreach (var g in rows.GroupBy(r => (r.Regime, r.Arm)))
        {
            double[] med = Enumerable.Range(0, names.Length).Select(i => Median(g.Select(r => r.Values[i]))).ToArray();
            summary.Append(g.Key.Regime).Append('\t').Append(g.Key.Arm).Append('\t').Append(string.Join('\t', med.Select(F)))
                   .Append('\t').Append(F(Quantile(g.Select(r => r.Values[4]), 0.25))).Append('\t').Append(F(Quantile(g.Select(r => r.Values[4]), 0.75)))
                   .Append('\t').Append(F(Quantile(g.Select(r => r.Values[7]), 0.25))).Append('\t').Append(F(Quantile(g.Select(r => r.Values[7]), 0.75)))
                   .Append('\n');
            _out.WriteLine($"[schedule] {g.Key.Regime,-8} {g.Key.Arm,-12}: ESJD {med[1]:F4}  reversal {med[2]:F3}  ESS k {med[4],6:F0}  " +
                           $"energy {med[5],6:F0}  fit {med[6],6:F0}  peak {med[7],6:F0}  k {med[8]:F2}  {med[0]:F2} s");
        }
        ArtifactRun run = HarnessArtifacts.Create("test-runs", "MixingBenchmark", "Maths.Regression.Tests",
                                                  new Dictionary<string, object?> { ["Study"] = "schedule", ["Seeds"] = seeds });
        _out.WriteLine($"[schedule] table {run.WriteText("schedule.tsv", table.ToString())}");
        _out.WriteLine($"[schedule] summary {run.WriteText("schedule-summary.tsv", summary.ToString())}");
    }

    // The conditional-mean fit at x = 0.5 and its argmax, for each draw; repeated configurations reuse the last result.
    private static (double[] FitMid, double[] Peak) Readouts(SplineBasis basis, IObservationModel model, Regime regime, KnotConfig[] configs)
    {
        var fitMid = new double[configs.Length];
        var peak = new double[configs.Length];
        KnotConfig? last = null;
        double f = double.NaN, p = double.NaN;
        for (int i = 0; i < configs.Length; i++)
        {
            if (!ReferenceEquals(configs[i], last))
            {
                double[] coef = model.PosteriorMeanCoefficients(basis.Design(configs[i], regime.X), regime.Y, null);
                f = basis.Evaluate(configs[i], coef, 0.5);
                p = SplineExtrema.Argmax(configs[i], coef, basis).Location;
                last = configs[i];
            }
            fitMid[i] = f;
            peak[i] = p;
        }
        return (fitMid, peak);
    }
}
