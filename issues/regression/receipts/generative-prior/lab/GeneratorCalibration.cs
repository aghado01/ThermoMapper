using System.Diagnostics;
using System.Text.Json;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace GenerativePrior;

internal sealed record InversionSettings(int MasterSeed, int[] Dimensions, int[] Ratios, int Repetitions,
    double Box, double TotalNoiseDeviation, double VarianceDegreesOfFreedom, int BurnIn, int DecorrelationEvery,
    int Chains, int MaxSamples, int BatchSize, double RHatTarget, double MinEssPerChain, double[] Ladder,
    int Particles, int SmcPopulations, double TargetConditionalEss, int MutationSteps, double ResampleBelow,
    int MaxStagesPerLeg, double SmcEssFraction, double PopulationSpreadTarget, int RequiredGoodProblems,
    double RecoveryTarget, double ExpansionSeconds, int AlignmentStates, int[] GridSizes,
    double GridMeanTolerance, double GridQuantileTolerance, double SamplerGridMeanTolerance, double ReferenceTolerance);

internal sealed class GeneratorEmulator(Settings settings, Network network, IReadOnlyList<CalibrationParameter> parameters) : ICalibrationEmulator
{
    public IReadOnlyList<CalibrationParameter> Parameters => parameters;
    public int DrawCount => 1;
    public IReadOnlyList<double>? DrawWeights => null;
    public int Outputs => settings.Measurements;
    public EmulatorNoise Noise(int draw) => new IndependentNoise(0);
    public double[,]? TruncationResiduals() => null;
    public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
    {
        var state = new Forward(settings);
        network.Predict(continuous.ToArray(), state);
        state.Residual.AsSpan().CopyTo(output);
    }
}

internal static class GeneratorCalibration
{
    private sealed record Reading(int K, int Ratio, int Repetition, int DataSeed, int SamplerSeed, string Method,
        double[] Truth, double[] Mean, double[] Q05, double[] Q95, int Covered, double RelativeError,
        bool DiagnosticGate, string Pooling, double? MaximumRHat, double? MinimumEssPerChain, double? MinimumParticleEss,
        int? MinimumFamilies, double? PopulationSpread, double NoiseDeviation, double NegativeMass, double Seconds,
        double[] AlignmentClean, double[] AlignmentNoisy, double? GridMeanError);

    private static double Dot(double[] a, double[] b) => PairwiseSum.Dot(a, b);
    private static double Norm(double[] a) => Math.Sqrt(Dot(a, a));
    private static double? Finite(double x) => double.IsFinite(x) ? x : null;

    private static double Quantile(double[] values, double[] weights, double probability)
    {
        int[] order = Enumerable.Range(0, values.Length).OrderBy(i => values[i]).ToArray();
        double total = weights.Sum(), current = 0;
        foreach (int index in order)
        {
            current += weights[index];
            if (current >= total * probability) return values[index];
        }
        return values[order[^1]];
    }

    private static double[] BoxCode(int k, double box, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        while (true)
        {
            double[] code = Enumerable.Range(0, k).Select(_ => Normal.SampleStandard(rng)).ToArray();
            if (code.All(x => Math.Abs(x) <= box)) return code;
        }
    }

    private static (double[] Mean, double[] Lower, double[] Upper) Grid(Settings settings, Network net, double[] y,
        double deviation, InversionSettings config, int size)
    {
        var points = new List<double[]>(size * size);
        var logs = new List<double>(size * size);
        var state = new Forward(settings);
        double beta = config.VarianceDegreesOfFreedom * deviation * deviation / 2;
        double power = (config.VarianceDegreesOfFreedom + y.Length) / 2;
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
            {
                double[] code = [-config.Box + (i + 0.5) * 2 * config.Box / size,
                                 -config.Box + (j + 0.5) * 2 * config.Box / size];
                net.Evaluate(code, y, state);
                points.Add(code);
                logs.Add(-power * Math.Log(beta + state.Objective));
            }
        double maximum = logs.Max();
        double[] weights = logs.Select(x => Math.Exp(x - maximum)).ToArray();
        double total = weights.Sum();
        var mean = new double[2]; var lower = new double[2]; var upper = new double[2];
        for (int j = 0; j < 2; j++)
        {
            double[] values = points.Select(p => p[j]).ToArray();
            mean[j] = values.Zip(weights, (v, w) => v * w).Sum() / total;
            lower[j] = Quantile(values, weights, 0.05);
            upper[j] = Quantile(values, weights, 0.95);
        }
        return (mean, lower, upper);
    }

    private static double Spread(IReadOnlyList<CalibrationDraw[]> populations, IReadOnlyList<double[]> weights, int k)
    {
        double maximum = 0;
        for (int j = 0; j < k; j++)
        {
            double[] means = populations.Select((p, i) => p.Select((d, n) => d.Theta[j] * weights[i][n]).Sum() / weights[i].Sum()).ToArray();
            double grand = means.Average();
            double variance = populations.Select((p, i) => p.Select((d, n) => weights[i][n] * Math.Pow(d.Theta[j] - grand, 2)).Sum() / weights[i].Sum()).Average();
            maximum = Math.Max(maximum, variance > 0 ? (means.Max() - means.Min()) / Math.Sqrt(variance) : 0);
        }
        return maximum;
    }

    private static double Cosine(double[] gradient, double[] truth, double[] theta)
    {
        double[] direction = truth.Zip(theta, (a, b) => a - b).ToArray();
        double denominator = Norm(gradient) * Norm(direction);
        return denominator > 0 ? -Dot(gradient, direction) / denominator : 0;
    }

    public static int Run(InversionSettings config, string output, bool replayCase = false)
    {
        var readings = new List<Reading>();
        using var rowsFile = new StreamWriter(Path.Combine(output, "calibrations.jsonl")) { AutoFlush = true };
        using var drawsFile = new StreamWriter(Path.Combine(output, "alignment-states.jsonl")) { AutoFlush = true };
        using var fullDrawsFile = new StreamWriter(Path.Combine(output, "posterior-draws.jsonl"));
        using var dataFile = new StreamWriter(Path.Combine(output, "problems.jsonl")) { AutoFlush = true };
        int[] masters = SeedTree.Derive(config.MasterSeed, config.Dimensions.Length);
        bool expanded = false;
        for (int dimension = 0; dimension < config.Dimensions.Length; dimension++)
        {
            int k = config.Dimensions[dimension];
            if (replayCase && k != 2) continue;
            if (k == 16)
            {
                expanded = Groups(readings, config).Where(g => g.K == 8 && g.Ratio >= 4).All(g =>
                    g.Passed && g.MedianSeconds <= config.ExpansionSeconds);
                if (!expanded) break;
            }
            int[] seeds = SeedTree.Derive(masters[dimension], config.Repetitions);
            for (int rep = 0; rep < seeds.Length; rep++)
            {
                if (replayCase && rep != 0) continue;
                int[] branches = SeedTree.Derive(seeds[rep], 10);
                double[] truth = BoxCode(k, config.Box, branches[3]), start = BoxCode(k, config.Box, branches[4]);
                var noiseRng = new Xoshiro256PlusPlus(branches[5]);
                double[] noise = Enumerable.Range(0, k * config.Ratios.Max()).Select(_ => Normal.SampleStandard(noiseRng)).ToArray();
                foreach (int ratio in config.Ratios)
                {
                    if (replayCase && ratio != 4) continue;
                    int m = k * ratio;
                    var settings = new Settings(k, 250, 600, m, 1, seeds[rep], 1, Math.ScaleB(1, -52),
                        50000, 0.001, config.ReferenceTolerance, true, true, false);
                    var net = new Network(settings, branches);
                    var state = new Forward(settings);
                    net.Predict(truth, state);
                    double[] clean = (double[])state.Residual.Clone();
                    double deviation = config.TotalNoiseDeviation / Math.Sqrt(m);
                    double[] y = clean.Select((v, j) => v + deviation * noise[j]).ToArray();
                    net.Evaluate(start, y, state);
                    var gradient = new double[k];
                    net.Gradient(state, gradient);
                    double check = Reference.Check(net, start, y, state.Objective, gradient);
                    if (!(check <= config.ReferenceTolerance)) throw new ArithmeticException("Scalar generator control failed.");
                    var descentClock = Stopwatch.StartNew();
                    Program.Endpoint descent = Program.Recover(settings, net, y, start);
                    double descentSeconds = descentClock.Elapsed.TotalSeconds;
                    dataFile.WriteLine(JsonSerializer.Serialize(new { k, ratio, Repetition = rep + 1, DataSeed = seeds[rep],
                        Branches = branches, Truth = truth, Start = start, Observations = y, Clean = clean, deviation,
                        Descent = descent, DescentSeconds = descentSeconds, ReferenceDiscrepancy = check }));
                    (double[] Mean, double[] Lower, double[] Upper)? fineGrid = null;
                    if (k == 2 && rep == 0 && ratio == 4)
                    {
                        var coarse = Grid(settings, net, y, deviation, config, config.GridSizes[0]);
                        var fine = Grid(settings, net, y, deviation, config, config.GridSizes[1]);
                        bool passed = coarse.Mean.Zip(fine.Mean, (a, b) => Math.Abs(a - b)).Max() <= config.GridMeanTolerance &&
                            coarse.Lower.Zip(fine.Lower, (a, b) => Math.Abs(a - b)).Max() <= config.GridQuantileTolerance &&
                            coarse.Upper.Zip(fine.Upper, (a, b) => Math.Abs(a - b)).Max() <= config.GridQuantileTolerance;
                        LabProgram.Write(Path.Combine(output, "grid-control.json"), new { Coarse = new { coarse.Mean, coarse.Lower, coarse.Upper },
                            Fine = new { fine.Mean, fine.Lower, fine.Upper }, Passed = passed });
                        if (!passed) throw new ArithmeticException("Independent grid resolution control failed.");
                        fineGrid = fine;
                    }
                    IReadOnlyList<CalibrationParameter> parameters = Enumerable.Range(0, k)
                        .Select(j => (CalibrationParameter)new ContinuousParameter($"x{j}", -config.Box, config.Box)).ToArray();
                    var emulator = new GeneratorEmulator(settings, net, parameters);
                    foreach (string method in new[] { "nested", "tempered", "smc" })
                    {
                        int samplerSeed = branches[method == "nested" ? 6 : method == "tempered" ? 7 : 8];
                        CalibrationConfig MakeConfig(int seed) => new() {
                            Parameters = parameters,
                            Measurement = new GaussianMeasurement([new VarianceGroup(Enumerable.Range(0, m).ToArray(),
                                new InverseGammaVariance(config.VarianceDegreesOfFreedom, deviation))]),
                            Covariance = new FunctionalFormulation(),
                            Sampling = new CalibrationSampling { MasterSeed = seed, BurnIn = config.BurnIn,
                                DecorrelationEvery = config.DecorrelationEvery, Tempering = method == "tempered" ? config.Ladder : null },
                            Reference = new ReferenceSampling { OuterDraws = 1, Chains = config.Chains, MaxSamples = config.MaxSamples,
                                BatchSize = config.BatchSize, RHatTarget = config.RHatTarget, MinEssPerChain = config.MinEssPerChain },
                            Smc = method == "smc" ? new CutSmcSpec { Particles = config.Particles, OuterDraws = 1,
                                TargetConditionalEss = config.TargetConditionalEss, MutationSteps = config.MutationSteps,
                                ResampleBelow = config.ResampleBelow, MaxStagesPerLeg = config.MaxStagesPerLeg } : null
                        };
                        var clock = Stopwatch.StartNew();
                        var populations = new List<CalibrationDraw[]>();
                        var populationWeights = new List<double[]>();
                        bool diagnostic;
                        double? maximumRhat = null, minimumEss = null;
                        int families = 0;
                        double spread = 0;
                        if (method == "smc")
                        {
                            int[] populationSeeds = SeedTree.Derive(samplerSeed, config.SmcPopulations);
                            double smallestEss = double.PositiveInfinity;
                            families = int.MaxValue;
                            foreach (int seed in populationSeeds)
                            {
                                var model = new CalibrationModel(MakeConfig(seed), y, emulator);
                                CutSmcCalibration result = CutSmcSampler.Run(model, 1);
                                var pooled = result.Pooled();
                                populations.Add(pooled.Draws.ToArray());
                                populationWeights.Add(pooled.Weights);
                                smallestEss = Math.Min(smallestEss, 1 / pooled.Weights.Sum(w => w * w));
                                families = Math.Min(families, result.Populations.Min(p => p.Leg.Families));
                            }
                            spread = Spread(populations, populationWeights, k);
                            minimumEss = smallestEss;
                            diagnostic = smallestEss >= config.SmcEssFraction * config.Particles && spread <= config.PopulationSpreadTarget;
                        }
                        else
                        {
                            var model = new CalibrationModel(MakeConfig(samplerSeed), y, emulator);
                            NestedCalibration result = NestedCalibrationSampler.Run(model);
                            var pooled = result.Pooled(includeNotConverged: true);
                            populations.Add(pooled.Draws.ToArray()); populationWeights.Add(pooled.Weights);
                            diagnostic = result.Runs.All(run => run.Converged);
                            var diagnostics = new List<ConvergenceDiagnostic>();
                            foreach (InnerRun run in result.Runs)
                                for (int j = 0; j <= k; j++)
                                    diagnostics.Add(ChainDiagnostics.Diagnose(run.Chains.Select(chain => chain.Select(d =>
                                        j == k ? Math.Log(d.Variances[0]) : d.Theta[j]).ToArray()).ToArray()));
                            maximumRhat = Finite(diagnostics.Where(d => double.IsFinite(d.RHat.Value)).Select(d => d.RHat.Value).DefaultIfEmpty(double.NaN).Max());
                            minimumEss = Finite(diagnostics.Where(d => double.IsFinite(d.Ess.Value)).Select(d => d.Ess.Value / config.Chains).DefaultIfEmpty(double.NaN).Min());
                        }
                        double seconds = clock.Elapsed.TotalSeconds;
                        CalibrationDraw[] draws = populations.SelectMany(p => p).ToArray();
                        double[] weights = populationWeights.SelectMany(w => w.Select(x => x / populations.Count)).ToArray();
                        double total = weights.Sum();
                        for (int i = 0; i < draws.Length; i++)
                            fullDrawsFile.WriteLine(JsonSerializer.Serialize(new {
                                k, ratio, Repetition = rep + 1, Method = method, SamplerSeed = samplerSeed, Index = i,
                                Population = method == "smc" ? i / config.Particles : 0,
                                Chain = method == "smc" ? (int?)null : i / (draws.Length / config.Chains),
                                Weight = weights[i], draws[i].Theta, draws[i].Variances, draws[i].LogLikelihood
                            }));
                        fullDrawsFile.Flush();
                        var mean = new double[k]; var lower = new double[k]; var upper = new double[k];
                        int covered = 0;
                        for (int j = 0; j < k; j++)
                        {
                            double[] values = draws.Select(d => d.Theta[j]).ToArray();
                            mean[j] = values.Zip(weights, (v, w) => v * w).Sum() / total;
                            lower[j] = Quantile(values, weights, 0.05); upper[j] = Quantile(values, weights, 0.95);
                            if (truth[j] >= lower[j] && truth[j] <= upper[j]) covered++;
                        }
                        double error = Norm(mean.Zip(truth, (a, b) => a - b).ToArray()) / Norm(truth);
                        double negativeMass = 0;
                        for (int i = 0; i < draws.Length; i++)
                        {
                            double norm = Norm(draws[i].Theta);
                            double cosine = norm > 0 ? Dot(draws[i].Theta, truth) / (norm * Norm(truth)) : 1;
                            double radius = norm / Norm(truth);
                            if (cosine <= -Math.Cos(Math.PI / 9) && Math.Abs(radius - 1 / Math.PI) <= 1 / (4 * Math.PI))
                                negativeMass += weights[i] / total;
                        }
                        double[] alignmentClean = new double[Math.Min(config.AlignmentStates, draws.Length)];
                        double[] alignmentNoisy = new double[alignmentClean.Length];
                        for (int i = 0; i < alignmentClean.Length; i++)
                        {
                            int index = i * draws.Length / alignmentClean.Length;
                            double[] theta = draws[index].Theta;
                            net.Evaluate(theta, clean, state); net.Gradient(state, gradient);
                            alignmentClean[i] = Cosine(gradient, truth, theta);
                            net.Evaluate(theta, y, state); net.Gradient(state, gradient);
                            alignmentNoisy[i] = Cosine(gradient, truth, theta);
                            drawsFile.WriteLine(JsonSerializer.Serialize(new { k, ratio, Repetition = rep + 1, Method = method,
                                Index = index, Weight = weights[index], Theta = theta, draws[index].Variances,
                                RelativeError = Norm(theta.Zip(truth, (a, b) => a - b).ToArray()) / Norm(truth),
                                CleanCosine = alignmentClean[i], NoisyCosine = alignmentNoisy[i] }));
                        }
                        double sigma = Math.Sqrt(draws.Select((d, i) => weights[i] * d.Variances[0]).Sum() / total);
                        double? gridError = fineGrid.HasValue ? mean.Zip(fineGrid.Value.Mean, (a, b) => Math.Abs(a - b)).Max() : null;
                        var reading = new Reading(k, ratio, rep + 1, seeds[rep], samplerSeed, method, truth, mean, lower, upper,
                            covered, error, diagnostic, method == "smc" ? "three independent particle populations" :
                                diagnostic ? "converged inner ensemble" : "diagnostic pool including unconverged ensemble",
                            maximumRhat, method == "smc" ? null : minimumEss,
                            method == "smc" ? minimumEss : null, method == "smc" ? families : null,
                            method == "smc" ? spread : null, sigma, negativeMass,
                            seconds, alignmentClean, alignmentNoisy, gridError);
                        readings.Add(reading); rowsFile.WriteLine(JsonSerializer.Serialize(reading));
                        Console.WriteLine($"G2 k={k} ratio={ratio} rep={rep + 1} {method} error={error:G4} diagnostic={diagnostic} seconds={seconds:F2}");
                    }
                }
            }
        }
        LabProgram.Write(Path.Combine(output, "summary.json"), new { LabProgram.Revision, LabProgram.SettingsSha256,
            Groups = Groups(readings, config), ExpandedTo16 = expanded,
            GridSamplerMeanGates = readings.Where(r => r.K == 2 && r.Ratio == 4 && r.Repetition == 1)
                .Select(r => new { r.Method, r.GridMeanError, Passed = r.GridMeanError <= config.SamplerGridMeanTolerance }) });
        return 0;
    }

    private sealed record Group(int K, int Ratio, string Method, double MedianError, double Coverage,
        double? DiagnosticQualifiedCoverage, bool IncludesFailedDiagnostics,
        int DiagnosticPasses, double MedianSeconds, double MeanNegativeMass, bool Passed);

    private static Group[] Groups(IEnumerable<Reading> readings, InversionSettings s) =>
        readings.GroupBy(r => new { r.K, r.Ratio, r.Method }).Select(g => {
            double[] errors = g.Select(r => r.RelativeError).Order().ToArray();
            double[] seconds = g.Select(r => r.Seconds).Order().ToArray();
            double median = 0.5 * (errors[(errors.Length - 1) / 2] + errors[errors.Length / 2]);
            Reading[] good = g.Where(r => r.DiagnosticGate).ToArray();
            return new Group(g.Key.K, g.Key.Ratio, g.Key.Method, median, g.Sum(r => r.Covered) / (double)(g.Count() * g.Key.K),
                good.Length > 0 ? good.Sum(r => r.Covered) / (double)(good.Length * g.Key.K) : null, good.Length < g.Count(),
                g.Count(r => r.DiagnosticGate), 0.5 * (seconds[(seconds.Length - 1) / 2] + seconds[seconds.Length / 2]),
                g.Average(r => r.NegativeMass), median <= s.RecoveryTarget && g.Count(r => r.DiagnosticGate) >= s.RequiredGoodProblems);
        }).ToArray();
}
