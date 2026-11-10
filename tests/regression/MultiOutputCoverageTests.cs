using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputCoverageTests(ITestOutputHelper output)
{
    internal sealed record Reading(double CorrelationCoverage, double JointCoverage, double IntervalWidth, double BandWidth,
                                   double CorrelationFlip, double JointFlip);

    [Fact]
    public void MatchedCovariancePriors_AtThreeToThirtyCasesPerOutput()
    {
        if (!BassQualification.Enabled) return;
        int replicates = int.TryParse(Environment.GetEnvironmentVariable("E34_COVERAGE_REPLICATES"), out int r) ? r : 400;
        const int retained = 4000, burn = 1000;
        var cells = (from k in new[] { 2, 5 } from ratio in new[] { 3, 5, 10, 30 } from rho in new[] { 0.0, 0.7 } select (k, ratio, rho)).ToArray();
        int[] seeds = SeedTree.Derive(20261215, cells.Length);
        var csv = new List<string> { "K,N,rho,replicate,prior,correlation_coverage,joint_coverage,interval_width,band_width,correlation_half_flip,joint_half_flip" };
        for (int cell = 0; cell < cells.Length; cell++)
        {
            var (k, ratio, rho) = cells[cell]; int n = k * ratio;
            var truth = new double[k, k];
            for (int a = 0; a < k; a++) for (int b = 0; b < k; b++)
                truth[a, b] = (a == b ? 1 : rho * (a % 2 == b % 2 ? 1 : -1)) * (0.7 + 0.4 * a) * (0.7 + 0.4 * b);
            double[,] root = ReferenceAlgebra.Cholesky(truth)!;
            int[] dataSeeds = SeedTree.Derive(seeds[cell], replicates);
            var results = new Reading[2][] { new Reading[replicates], new Reading[replicates] };
            for (int replicate = 0; replicate < replicates; replicate++)
            {
                int[] streams = SeedTree.Derive(dataSeeds[replicate], 5);
                var dataRng = new Xoshiro256PlusPlus(streams[0]);
                var y = new double[n, k]; var mean = new double[k];
                for (int i = 0; i < n; i++) { double[] row = Gaussian(root, dataRng); for (int a = 0; a < k; a++) { y[i, a] = row[a]; mean[a] += row[a] / n; } }
                double[] future = Gaussian(root, dataRng);
                var residual = new double[k, k];
                for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) for (int i = 0; i < n; i++) residual[a, b] += (y[i, a] - mean[a]) * (y[i, b] - mean[b]);
                double[] scales = [.. Enumerable.Range(0, k).Select(a => Math.Sqrt(residual[a, a] / n))];
                // Mean-centred intercept inside the prior, held tau=1/N: the actual E34 count is N.
                double h = 1 + 1.0 / (n * (1 + 1.0 / n));
                for (int prior = 0; prior < 2; prior++)
                {
                    var covariance = new double[2 * retained][,]; var noise = new double[2 * retained][];
                    for (int chain = 0; chain < 2; chain++)
                    {
                        var rng = new Xoshiro256PlusPlus(streams[1 + 2 * prior + chain]);
                        var oracle = new CovarianceGibbsReference(residual, n, scales, prior == 1);
                        for (int step = 0; step < burn + retained; step++)
                        {
                            double[,] sigma = oracle.Step(rng);
                            if (step < burn) continue;
                            int at = chain * retained + step - burn;
                            covariance[at] = sigma;
                            noise[at] = Gaussian(ReferenceAlgebra.Cholesky(sigma)!, rng);
                            for (int a = 0; a < k; a++) noise[at][a] *= Math.Sqrt(h);
                        }
                    }
                    var full = Summarize(covariance, noise, 0, 2 * retained, truth, mean, future, h);
                    var first = Summarize(covariance, noise, 0, retained, truth, mean, future, h);
                    var second = Summarize(covariance, noise, retained, retained, truth, mean, future, h);
                    double flip = 0;
                    for (int pair = 0; pair < full.Covered.Length; pair++) if (first.Covered[pair] != second.Covered[pair]) flip += 1.0 / full.Covered.Length;
                    var reading = new Reading(full.Covered.Count(v => v) / (double)full.Covered.Length, full.Joint ? 1 : 0,
                        full.IntervalWidth, full.BandWidth, flip, first.Joint != second.Joint ? 1 : 0);
                    results[prior][replicate] = reading;
                    csv.Add(FormattableString.Invariant($"{k},{n},{rho},{replicate},{(prior == 0 ? "scaled" : "huang-wand")},{reading.CorrelationCoverage},{reading.JointCoverage},{reading.IntervalWidth},{reading.BandWidth},{reading.CorrelationFlip},{reading.JointFlip}"));
                }
            }
            foreach (int prior in new[] { 0, 1 })
            {
                var a = results[prior];
                var corr = MeanError(a.Select(t => t.CorrelationCoverage)); var joint = MeanError(a.Select(t => t.JointCoverage));
                var interval = Wilson(joint.Mean, replicates);
                output.WriteLine($"K={k}, N={n}, rho={rho}, {(prior == 0 ? "scaled" : "HW")}: correlation {corr.Mean:F4} ± {1.96 * corr.Se:F4}, joint {joint.Mean:F4} [{interval.Lo:F4},{interval.Hi:F4}]; widths {a.Average(t => t.IntervalWidth):F4}/{a.Average(t => t.BandWidth):F4}; half-chain decision flips {a.Average(t => t.CorrelationFlip):F4}/{a.Average(t => t.JointFlip):F4}");
            }
            var dc = MeanError(results[1].Zip(results[0], (a, b) => a.CorrelationCoverage - b.CorrelationCoverage));
            var dj = MeanError(results[1].Zip(results[0], (a, b) => a.JointCoverage - b.JointCoverage));
            double cHalf = dc.Se > 0 ? 1.96 * dc.Se : 1 - Math.Pow(0.05, 1.0 / replicates);
            double jHalf = dj.Se > 0 ? 1.96 * dj.Se : 1 - Math.Pow(0.05, 1.0 / replicates);
            output.WriteLine($"paired HW-minus-scaled: correlation {dc.Mean:F4} ± {cHalf:F4}, joint {dj.Mean:F4} ± {jHalf:F4}; datasets={replicates}, posterior draws={2 * retained}; zero-discordance cells use their 95% upper probability bound");
            if (Environment.GetEnvironmentVariable("E34_COVERAGE_OUT") is { Length: > 0 } path) File.WriteAllLines(path, csv);
        }
    }

    internal static (bool[] Covered, bool Joint, double IntervalWidth, double BandWidth) Summarize(double[][,] cov, double[][] noise, int offset, int count, double[,] truth, double[] mean, double[] future, double h)
    {
        int k = mean.Length; var sd = new double[k]; var covered = new List<bool>(); double interval = 0;
        for (int a = 0; a < k; a++) sd[a] = Math.Sqrt(h * Enumerable.Range(offset, count).Average(i => cov[i][a, a]));
        for (int a = 0; a < k; a++) for (int b = 0; b < a; b++)
        {
            double actual = truth[a, b] / Math.Sqrt(truth[a, a] * truth[b, b]);
            double[] rho = [.. Enumerable.Range(offset, count).Select(i => cov[i][a, b] / Math.Sqrt(cov[i][a, a] * cov[i][b, b])).Order()];
            double lo = QuantileSorted(rho, 0.025), hi = QuantileSorted(rho, 0.975);
            covered.Add(lo <= actual && actual <= hi); interval += (hi - lo) / (k * (k - 1) / 2);
        }
        double[] maxima = [.. Enumerable.Range(offset, count).Select(i => Enumerable.Range(0, k).Max(a => Math.Abs(noise[i][a]) / sd[a])).Order()];
        double q = QuantileSorted(maxima, 0.95);
        bool joint = Enumerable.Range(0, k).All(a => Math.Abs(future[a] - mean[a]) <= q * sd[a]);
        return ([.. covered], joint, interval, sd.Average() * 2 * q);
    }

    internal static double[] Gaussian(double[,] root, Xoshiro256PlusPlus rng)
    {
        int k = root.GetLength(0); var z = new double[k]; var y = new double[k];
        for (int i = 0; i < k; i++) z[i] = Normal.SampleStandard(rng);
        for (int i = 0; i < k; i++) for (int j = 0; j <= i; j++) y[i] += root[i, j] * z[j];
        return y;
    }

    internal static double QuantileSorted(double[] values, double p)
    {
        double at = p * (values.Length - 1); int lo = (int)at, hi = Math.Min(lo + 1, values.Length - 1);
        return values[lo] + (at - lo) * (values[hi] - values[lo]);
    }

    private static (double Mean, double Se) MeanError(IEnumerable<double> sequence)
    {
        double[] x = [.. sequence]; double mean = x.Average();
        return (mean, Math.Sqrt(x.Sum(v => (v - mean) * (v - mean)) / (x.Length * (x.Length - 1))));
    }

    private static (double Lo, double Hi) Wilson(double p, int n)
    {
        const double z = 1.959963984540054;
        double denominator = 1 + z * z / n, centre = (p + z * z / (2 * n)) / denominator;
        double half = z * Math.Sqrt(p * (1 - p) / n + z * z / (4.0 * n * n)) / denominator;
        return (centre - half, centre + half);
    }
}
