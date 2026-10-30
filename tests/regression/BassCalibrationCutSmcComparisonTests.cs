using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Stage 8's C5 comparison on stage 7's recovery fixture: exploration, which reads and recommends nothing. The tour
/// against a random order; cut SMC against the nested reference (C3) at matched error on E[θ]; how many of C3's inner
/// runs converge under the default gate (R̂ 1.01, ESS 100 per chain, 4,000 draws per chain), stage 7's open item; and
/// θ*'s coverage over five seeds, as stage 7 read it. Heavy: set CUT_SMC_COMPARISON=1 to run it.
/// </summary>
public sealed class BassCalibrationCutSmcComparisonTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationCutSmcComparisonTests(ITestOutputHelper output) => _out = output;

    private static bool Enabled => Environment.GetEnvironmentVariable("CUT_SMC_COMPARISON") is "1";

    private static readonly ContinuousParameter[] Parameters = [BassCalibrationRecoveryTests.Theta1, BassCalibrationRecoveryTests.Theta2];
    private static readonly double[] Sd = [.. Enumerable.Repeat(0.05, BassCalibrationRecoveryTests.M)];

    private static BassBasisEmulator Basis() => new(BassCalibrationRecoveryTests.BasisFit.Value, Parameters);

    private static BassFitEmulator Direct() => new(BassCalibrationRecoveryTests.DirectFit.Value.Fit, BassCalibrationRecoveryTests.DirectFit.Value.Data, Parameters);

    private static CalibrationConfig Config(int seed, CutSmcSpec? smc, GaussianMeasurement? groups = null, CovarianceFormulation? covariance = null) => new()
    {
        Parameters = Parameters,
        Measurement = groups ?? BassCalibrationRecoveryTests.OneGroup(),
        Covariance = covariance ?? new FunctionalFormulation(),
        Sampling = new CalibrationSampling { MasterSeed = seed },
        Smc = smc,
    };

    /// <summary>A fixed subset of an emulator's draws, every one taken once, so replicates share their outer draws and
    /// their spread is the inner error alone.</summary>
    private sealed class Subset(ICalibrationEmulator inner, int[] draws) : ICalibrationEmulator
    {
        public IReadOnlyList<CalibrationParameter> Parameters => inner.Parameters;
        public int DrawCount => draws.Length;
        public IReadOnlyList<double>? DrawWeights => null;
        public int Outputs => inner.Outputs;
        public EmulatorNoise Noise(int draw) => inner.Noise(draws[draw]);
        public double[,]? TruncationResiduals() => inner.TruncationResiduals();
        public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output) => inner.Predict(draws[draw], continuous, levels, output);
        public void PredictLevels(int draw, ReadOnlySpan<double> continuous, IReadOnlyList<int[]> combinations, double[,] output)
            => inner.PredictLevels(draws[draw], continuous, combinations, output);
    }

    private static double[] Means(IReadOnlyList<CalibrationDraw> draws, IReadOnlyList<double> weights)
        => [.. Enumerable.Range(0, 2).Select(j => Enumerable.Range(0, draws.Count).Sum(i => weights[i] * draws[i].Theta[j]) / weights.Sum())];

    private static (double Mean, double Sd) Spread(IReadOnlyList<double> x)
    {
        double mean = x.Average();
        return (mean, Math.Sqrt(x.Sum(v => (v - mean) * (v - mean)) / (x.Count - 1)));
    }

    // Wall and processor seconds of `run`.
    private static (T Result, double Wall, double Cpu) Timed<T>(Func<T> run)
    {
        var process = Process.GetCurrentProcess();
        TimeSpan cpu = process.TotalProcessorTime;
        var watch = Stopwatch.StartNew();
        T result = run();
        watch.Stop();
        process.Refresh();
        return (result, watch.Elapsed.TotalSeconds, (process.TotalProcessorTime - cpu).TotalSeconds);
    }

    [Fact]
    public void TourAgainstRandomOrder()
    {
        if (!Enabled) return;
        // Twenty outer draws of the basis fit's 200, the same field data; each seed's tour and random runs share the
        // outer draws and the first leg, so the tour's distances measure both paths.
        BassBasisEmulator emulator = Basis();
        double[] y = BassCalibrationRecoveryTests.Field(1, Sd);
        _out.WriteLine("seed | order | path length | bridge stages | least ESS (bridges) | families left | wall s");
        var ratios = new List<double>();
        var rows = new List<(int Seed, string Order, double Length, int Stages, double LeastEss, int Families)>();
        foreach (int seed in Enumerable.Range(1, 5))
        {
            CutSmcCalibration? tour = null;
            foreach (CutSmcOrder order in new[] { CutSmcOrder.Tour, CutSmcOrder.Random })
            {
                var model = new CalibrationModel(Config(seed, new CutSmcSpec { OuterDraws = 20, Order = order }), y, emulator);
                (CutSmcCalibration result, double wall, _) = Timed(() => CutSmcSampler.Run(model));
                tour ??= result;
                double length = CutSmcOrdering.Length(tour.Distances!, [.. result.Order]);
                CutSmcPopulation[] bridges = [.. result.Populations.Skip(1)];
                int stages = bridges.Sum(p => p.Leg.StageCount);
                rows.Add((seed, order.ToString(), length, stages, bridges.Min(p => p.Leg.MinEss), result.Families));
                _out.WriteLine($"{seed} | {order} | {length:F2} | {stages} | {bridges.Min(p => p.Leg.MinEss):F0} | {result.Families} | {wall:F1}");
                foreach (CutSmcPopulation p in bridges)
                    if (tour.Distances![result.Order[p.Visit - 1], result.Order[p.Visit]] > 0)
                        ratios.Add(p.Leg.StageCount / tour.Distances[result.Order[p.Visit - 1], result.Order[p.Visit]]);
            }
        }
        foreach (string order in new[] { "Tour", "Random" })
        {
            var mine = rows.Where(r => r.Order == order).ToArray();
            _out.WriteLine($"{order}: mean path {mine.Average(r => r.Length):F2}, mean bridge stages {mine.Average(r => r.Stages):F1}, " +
                           $"mean least ESS {mine.Average(r => r.LeastEss):F0}, mean families {mine.Average(r => r.Families):F0}");
        }
        _out.WriteLine($"bridge stages per unit of distance: median {ratios.Order().ElementAt(ratios.Count / 2):F1} " +
                       $"(the conditional-ESS rule's √(2(1 − 0.99))⁻¹ = {1 / Math.Sqrt(0.02):F1} for small Gaussian steps)");
    }

    [Fact]
    public void CutSmcAgainstNested_AtMatchedError_AndTheDefaultGate()
    {
        if (!Enabled) return;
        // Six fixed draws of the basis fit, the same field data: five replicates of each method, so the spread of the
        // pooled E[θ] is each one's inner Monte Carlo error. C3 runs under the default gate, its every inner run pooled
        // for the error (its verdicts read beside); the cut SMC at its defaults.
        int[] six = NestedCalibrationSampler.OuterDraws(Basis(), 6, new Xoshiro256PlusPlus(41));
        foreach ((string label, ICalibrationEmulator emulator) in new (string, ICalibrationEmulator)[] { ("basis", new Subset(Basis(), six)), ("direct", new Subset(Direct(), six)) })
        {
            double[] y = BassCalibrationRecoveryTests.Field(1, Sd);
            const int r = 5;
            var smc = new List<(double[] Means, double Wall, double Cpu, int Stages)>();
            var nested = new List<(double[] Means, double Wall, double Cpu, NestedCalibration Result)>();
            for (int k = 0; k < r; k++)
            {
                var model = new CalibrationModel(Config(600 + k, new CutSmcSpec()), y, emulator);
                (CutSmcCalibration a, double wa, double ca) = Timed(() => CutSmcSampler.Run(model));
                (IReadOnlyList<CalibrationDraw> da, double[] pa) = a.Pooled();
                smc.Add((Means(da, pa), wa, ca, a.TotalStages));
                (NestedCalibration b, double wb, double cb) = Timed(() => NestedCalibrationSampler.Run(model));
                (IReadOnlyList<CalibrationDraw> db, double[] pb) = b.Pooled(includeNotConverged: true);
                nested.Add((Means(db, pb), wb, cb, b));
                _out.WriteLine($"{label}, replicate {k}: SMC E[θ] ({smc[^1].Means[0]:F5}, {smc[^1].Means[1]:F5}) in {wa:F1} s wall, {ca:F0} s CPU, {a.TotalStages} stages; " +
                               $"C3 ({nested[^1].Means[0]:F5}, {nested[^1].Means[1]:F5}) in {wb:F1} s wall, {cb:F0} s CPU, " +
                               $"{b.Runs.Count(run => run.Converged)}/{b.Runs.Count} converged, draws per chain {string.Join(",", b.Runs.Select(run => run.SamplesPerChain))}");
            }
            double smcWall = smc.Average(x => x.Wall), smcCpu = smc.Average(x => x.Cpu), refWall = nested.Average(x => x.Wall), refCpu = nested.Average(x => x.Cpu);
            for (int j = 0; j < 2; j++)
            {
                (double ms, double ss) = Spread([.. smc.Select(x => x.Means[j])]);
                (double mn, double sn) = Spread([.. nested.Select(x => x.Means[j])]);
                double ratio = (sn / ss) * (sn / ss);   // C3's work over SMC's, per run, to a matched error
                _out.WriteLine($"{label} θ{j + 1}: SMC {ms:F5} ± {ss:G2} per run, C3 {mn:F5} ± {sn:G2} per run (difference {ms - mn:G2}); " +
                               $"to SMC's error C3 takes {refWall * ratio:F1} s wall, {refCpu * ratio:F0} s CPU, against SMC's {smcWall:F1} s wall, {smcCpu:F0} s CPU");
            }
            InnerRun[] runs = [.. nested.SelectMany(x => x.Result.Runs)];
            _out.WriteLine($"{label}, C3 under the default gate: {runs.Count(x => x.Converged)}/{runs.Length} inner runs converged; " +
                           $"draws per chain: median {runs.Select(x => x.SamplesPerChain).Order().ElementAt(runs.Length / 2)}, max {runs.Max(x => x.SamplesPerChain)}");
            foreach (var failing in runs.Where(x => !x.Converged).SelectMany(x => x.Verdict.Groups.Where(g => g.Outcome != ConvergenceOutcome.Converged)
                                                                         .SelectMany(g => g.Failing.Select(f => (Group: g.Name, f.Name)))).GroupBy(f => f).OrderByDescending(g => g.Count()))
                _out.WriteLine($"  failing {failing.Key.Group}/{failing.Key.Name}: {failing.Count()} runs");
        }
    }

    [Fact]
    public void ThetaStarCoverage_OverFiveSeeds()
    {
        if (!Enabled) return;
        // Stage 7's four recovery cases and seeds, calibrated by the cut SMC at its defaults over six outer draws.
        var cases = new (string Label, Func<ICalibrationEmulator> Emulator, Func<int, double[]> Field, GaussianMeasurement Groups, CovarianceFormulation Covariance)[]
        {
            ("basis", Basis, seed => BassCalibrationRecoveryTests.Field(seed, Sd), BassCalibrationRecoveryTests.OneGroup(), new FunctionalFormulation()),
            ("direct", Direct, seed => BassCalibrationRecoveryTests.Field(seed, Sd), BassCalibrationRecoveryTests.OneGroup(), new FunctionalFormulation()),
            ("discrepancy", Basis, seed => BassCalibrationRecoveryTests.Field(seed, Sd, t => 0.1 + 0.1 * (t - 0.5)), BassCalibrationRecoveryTests.OneGroup(),
                new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy([.. BassCalibrationRecoveryTests.T.Select(t => new[] { 1.0, t - 0.5 })], 0.2) }),
            ("two groups", Basis, seed => BassCalibrationRecoveryTests.Field(seed, [.. Enumerable.Range(0, BassCalibrationRecoveryTests.M).Select(l => l < BassCalibrationRecoveryTests.M / 2 ? 0.03 : 0.08)]),
                new GaussianMeasurement([new VarianceGroup([.. Enumerable.Range(0, BassCalibrationRecoveryTests.M / 2)], new InverseGammaVariance(4, 0.05)),
                                         new VarianceGroup([.. Enumerable.Range(BassCalibrationRecoveryTests.M / 2, BassCalibrationRecoveryTests.M / 2)], new HalfCauchyDeviation(0.05))]),
                new FunctionalFormulation()),
        };
        foreach (var c in cases)
        {
            ICalibrationEmulator emulator = c.Emulator();
            int[] covered = new int[2];
            for (int seed = 1; seed <= 5; seed++)
            {
                var model = new CalibrationModel(Config(seed, new CutSmcSpec { OuterDraws = 6 }, c.Groups, c.Covariance), c.Field(seed), emulator);
                (CutSmcCalibration result, double wall, _) = Timed(() => CutSmcSampler.Run(model));
                (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
                var line = new List<string>();
                for (int j = 0; j < 2; j++)
                {
                    double[] values = [.. draws.Select(d => d.Theta[j])];
                    double lo = BassCalibrationRecoveryTests.WeightedQuantile(values, weights, 0.025), hi = BassCalibrationRecoveryTests.WeightedQuantile(values, weights, 0.975);
                    if (BassCalibrationRecoveryTests.ThetaStar[j] >= lo && BassCalibrationRecoveryTests.ThetaStar[j] <= hi) covered[j]++;
                    line.Add($"θ{j + 1} {Means(draws, weights)[j]:F3} [{lo:F3}, {hi:F3}]");
                }
                _out.WriteLine($"{c.Label}, seed {seed}: {string.Join("; ", line)}; {result.TotalStages} stages, families {result.Families}, {wall:F1} s");
            }
            _out.WriteLine($"{c.Label}: θ* covered {covered[0]}/5 and {covered[1]}/5");
        }
    }
}
