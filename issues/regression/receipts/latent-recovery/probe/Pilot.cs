using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;

namespace LatentRecoveryProbe;

/// <summary>
/// Part A: the synthesizer emulated by a BASS fit (code as inputs, the feature over frequency as the response), then
/// inverted from m kept bins by the calibration stack's nested reference sampler, beside the same inversion through the
/// synthesizer itself. Every stream derives from the settings' master seed by <see cref="SeedTree"/>.
/// </summary>
public static partial class Pilot
{
    public sealed record Settings(
        string Stage, int MasterSeed, int Replicates, double[] NoiseVariances, int[] Ratios, int TrainingBase, int TrainingPerK, int TestRuns,
        string[] Layouts, int MaxBasis, int MaxInteraction, int FitChains, int FitBurnIn, int FitSamples, int FitThin,
        int CalBurnIn, int CalDecorrelation, int OuterDraws, int CalChains, int MaxSamples, int BatchSize, double RHat, double MinEss,
        double PriorDf, double PriorSd, double DiscrepancyScale, int DiscrepancyNoise, double[] Ladder, int SmcParticles, int[] ModalityKs,
        int[] ModalityRatios, int ModalityNoise, int ModalityReplicates, int ReconstructionDraws, int Threads, int DirectionalColumns = 1);

    public static Settings Stage1 => new(
        Stage: "s1", MasterSeed: 2026111601, Replicates: 30, NoiseVariances: [1e-6, 1e-5, 1e-4], Ratios: [1, 2, 4, 8],
        TrainingBase: 80, TrainingPerK: 80, TestRuns: 100, Layouts: ["functional"], MaxBasis: 80, MaxInteraction: 3,
        FitChains: 2, FitBurnIn: 8000, FitSamples: 500, FitThin: 4,
        CalBurnIn: 1000, CalDecorrelation: 100, OuterDraws: 8, CalChains: 3, MaxSamples: 2000, BatchSize: 250, RHat: 1.05, MinEss: 50,
        PriorDf: 2, PriorSd: 0.3, DiscrepancyScale: 1.0, DiscrepancyNoise: 1, Ladder: [1.0, 0.6, 0.35, 0.2, 0.1, 0.05], SmcParticles: 400,
        ModalityKs: [2, 5], ModalityRatios: [1, 4], ModalityNoise: 1, ModalityReplicates: 10, ReconstructionDraws: 200, Threads: 6);

    public static bool Dispatch(string command, string dir, string[] args)
    {
        switch (command)
        {
            case "s1-peek": Peek(new ResonanceSynthesizer(5), Stage1, dir, "s1"); return true;
            case "s1-run": Run(Stage1, k => new ResonanceSynthesizer(k), int.Parse(args[0]), dir, args.Contains("rerun"), args.Contains("fit-only")); return true;
            case "s1-gate": File.WriteAllText(Path.Combine(dir, "s1-gate.txt"), PilotGates.Stage1(dir, Stage1)); Console.Write(File.ReadAllText(Path.Combine(dir, "s1-gate.txt"))); return true;
            default: return DispatchStage2(command, dir, args);
        }
    }

    private static string F(double v, string f = "G6") => v.ToString(f, CultureInfo.InvariantCulture);
    private static string J(IEnumerable<double> v, string f = "F5") => string.Join(";", v.Select(x => F(x, f)));

    /// <summary>Synthesizer levels at a few codes: design facts for the floor and the noise, read before any fit.</summary>
    public static void Peek(Synthesizer synth, Settings s, string dir, string stage)
    {
        var rng = new Xoshiro256PlusPlus(7);
        var sb = new StringBuilder();
        sb.AppendLine($"{synth.GetType().Name}: frame {synth.Frame}, fs {synth.SampleRate}, floor {synth.Floor}, L {synth.L}, f {F(synth.Frequencies[0])}..{F(synth.Frequencies[^1])}");
        for (int i = 0; i < 6; i++)
        {
            double[] theta = i == 0 ? synth.Lower.Zip(synth.Upper, (a, b) => 0.5 * (a + b)).ToArray() : synth.Draw(rng);
            double[] clean = synth.Clean(theta);
            double energy = synth.Energy(theta);
            sb.AppendLine($"theta {J(theta, "F3")}: energy {F(energy, "G4")}, feature min {F(clean.Min(), "F2")} max {F(clean.Max(), "F2")} " +
                          $"(ln floor {F(Math.Log(synth.Floor), "F2")}); SNR dB at noise levels: " +
                          string.Join(", ", s.NoiseVariances.Select(v => F(10 * Math.Log10(energy / (synth.Frame * v)), "F1"))));
            foreach (double v in s.NoiseVariances)
            {
                double[] noisy = synth.Noisy(theta, v, rng);
                var d = noisy.Zip(clean, (a, b) => a - b).ToArray();
                sb.AppendLine($"   noise {v}: noisy − clean mean {F(d.Average(), "F3")} rms {F(Math.Sqrt(d.Average(x => x * x)), "F3")} max {F(d.Max(), "F2")}");
            }
        }
        File.WriteAllText(Path.Combine(dir, $"{stage}-peek.txt"), sb.ToString());
        Console.Write(sb);
    }

    // Seeds per k: SeedTree.Derive(master, 5)[k − 1] → Derive(·, 16): [design, test, fit, truth, noise, subsets, calibration,
    // modality] for the main replicates, the same eight again for a rerun's fresh replicates (design, test and fit unused there).
    private static int[] KSeeds(Settings s, int k) => SeedTree.Derive(SeedTree.Derive(s.MasterSeed, 5)[k - 1], 16);

    public sealed record Layout(string Name, BassData Data, Func<IReadOnlyList<double>, BassFitEmulator> Emulator);

    private static Layout BuildLayout(string name, Synthesizer synth, double[][] design, IReadOnlyList<CalibrationParameter> parameters, double[] axis)
    {
        int n = design.Length, L = synth.L, k = synth.K;
        double[][] features = design.Select(t => synth.Clean(t)).ToArray();
        if (name == "functional")
        {
            var y = new double[n, L];
            for (int i = 0; i < n; i++) for (int l = 0; l < L; l++) y[i, l] = features[i][l];
            var predictors = Enumerable.Range(0, k).Select(j => (BassPredictor)new ContinuousPredictor(synth.Names[j], design.Select(t => t[j]).ToArray())).ToList();
            var data = new BassData(predictors, new FunctionalResponse(y, [new ContinuousPredictor("lnf", axis)]));
            return new Layout(name, data, null!);
        }
        if (name is "scalar-directional" or "scalar-axis")
        {
            var cols = Enumerable.Range(0, k + 1).Select(_ => new double[n * L]).ToArray();
            var y = new double[n * L];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < L; l++)
                {
                    for (int j = 0; j < k; j++) cols[j][i * L + l] = design[i][j];
                    cols[k][i * L + l] = axis[l];
                    y[i * L + l] = features[i][l];
                }
            var predictors = Enumerable.Range(0, k).Select(j => (BassPredictor)new ContinuousPredictor(synth.Names[j], cols[j])).Append(new ContinuousPredictor("lnf", cols[k])).ToList();
            return new Layout(name, new BassData(predictors, new ScalarResponse(y)), null!);
        }
        throw new ArgumentException($"Unknown layout '{name}'.");
    }

    private static BassFitEmulator EmulatorAt(Layout layout, BassFit fit, IReadOnlyList<CalibrationParameter> parameters, IReadOnlyList<double> keptAxis)
        => layout.Name == "functional"
            ? new BassFitEmulator(fit, layout.Data, parameters, null, [new ContinuousColumn("lnf", keptAxis.ToArray())])
            : new BassFitEmulator(fit, layout.Data, parameters, new PredictionCases([new ContinuousColumn("lnf", keptAxis.ToArray())]));

    private static BassConfig FitConfig(Settings s, string layout, int seed) => new()
    {
        Structure = layout == "scalar-directional"
            ? new BassStructureSpec { MaxInteraction = s.MaxInteraction, MaxBasis = s.MaxBasis, Directional = new DirectionalBasis { Columns = s.DirectionalColumns } }
            : new BassStructureSpec { MaxInteraction = s.MaxInteraction, MaxBasis = s.MaxBasis },
        Sampling = new BassSampling { Chains = s.FitChains, MasterSeed = seed, BurnIn = s.FitBurnIn, Samples = s.FitSamples, Thin = s.FitThin },
    };

    internal static double WeightedQuantile(double[] values, double[] weights, double p)
    {
        int[] order = [.. Enumerable.Range(0, values.Length).OrderBy(i => values[i])];
        double total = weights.Sum(), cumulative = 0.0;
        foreach (int i in order)
            if ((cumulative += weights[i] / total) >= p) return values[i];
        return values[order[^1]];
    }

    /// <summary>Systematic resampling of weighted draws to n equally weighted ones, offset ½ (no randomness).</summary>
    internal static List<T> Resample<T>(IReadOnlyList<T> draws, double[] weights, int n)
    {
        double total = weights.Sum(), step = total / n, u = 0.5 * step, cumulative = 0;
        var picked = new List<T>(n);
        int i = 0;
        for (int j = 0; j < n; j++)
        {
            double target = u + j * step;
            while (i < draws.Count - 1 && cumulative + weights[i] < target) cumulative += weights[i++];
            picked.Add(draws[i]);
        }
        return picked;
    }

    public sealed record FitReport(string Layout, int K, int TrainingRuns, int Rows, double Seconds, int Draws, double SigmaMean, double HeldOutRmse,
                                   double HeldOutR2, double HeldOutMaxBinRmse, double ResponseSd, string ResidualDependence);

    private static (BassFit Fit, FitReport Report) FitAndCheck(Settings s, Synthesizer synth, Layout layout, IReadOnlyList<CalibrationParameter> parameters,
                                                               double[] axis, int fitSeed, int testSeed, int n)
    {
        var sw = Stopwatch.StartNew();
        BassFit fit = Bass.Fit(FitConfig(s, layout.Name, fitSeed), layout.Data, new BassRunOptions { MaxDegreeOfParallelism = s.FitChains });
        double seconds = sw.Elapsed.TotalSeconds;
        BassFitEmulator full = EmulatorAt(layout, fit, parameters, axis);
        int draws = full.DrawCount;
        double sigma = Enumerable.Range(0, draws).Average(d => ((IndependentNoise)full.Noise(d)).Variance);
        var rng = new Xoshiro256PlusPlus(testSeed);
        int L = synth.L;
        double sse = 0, sst = 0;
        var binSse = new double[L];
        var all = new List<double>();
        var truths = new List<double[]>();
        var preds = new List<double[]>();
        for (int t = 0; t < s.TestRuns; t++)
        {
            double[] theta = synth.Draw(rng);
            double[] truth = synth.Clean(theta);
            var mean = new double[L];
            var buf = new double[L];
            for (int d = 0; d < draws; d++)
            {
                full.Predict(d, theta, [], buf);
                for (int l = 0; l < L; l++) mean[l] += buf[l] / draws;
            }
            truths.Add(truth);
            preds.Add(mean);
            all.AddRange(truth);
        }
        double grand = all.Average();
        for (int t = 0; t < truths.Count; t++)
            for (int l = 0; l < L; l++)
            {
                double e = preds[t][l] - truths[t][l];
                sse += e * e;
                binSse[l] += e * e;
                sst += (truths[t][l] - grand) * (truths[t][l] - grand);
            }
        int count = truths.Count * L;
        // Held-out lag-1 autocorrelation of the residuals along the frequency index, per test code, for every layout alike.
        var r1 = new List<double>();
        for (int t = 0; t < truths.Count; t++)
        {
            double[] e = Enumerable.Range(0, L).Select(l => preds[t][l] - truths[t][l]).ToArray();
            double me = e.Average(), num = 0, den = 0;
            for (int l = 0; l < L; l++) den += (e[l] - me) * (e[l] - me);
            for (int l = 1; l < L; l++) num += (e[l] - me) * (e[l - 1] - me);
            if (den > 0) r1.Add(num / den);
        }
        r1.Sort();
        string heldOut = $"held-out lag-1 r along the bins: median {F(r1[r1.Count / 2], "F3")}, 5%–95% {F(r1[(int)(0.05 * r1.Count)], "F3")}–{F(r1[(int)(0.95 * r1.Count)], "F3")}";
        string rd = (fit.ResidualDependence is { } r
            ? $"training lag-1 r along lnf: mean {F(r.Mean, "F3")} (independent noise {F(r.Expected, "F3")}), quantiles {J(r.Quantiles, "F3")} at {J(r.Probabilities, "F2")}"
            : "training check: none for a scalar layout") + "; " + heldOut;
        var report = new FitReport(layout.Name, synth.K, n, layout.Data.Response is FunctionalResponse ? n * L : n * L, seconds, draws, Math.Sqrt(sigma),
                                   Math.Sqrt(sse / count), 1 - sse / sst, Math.Sqrt(binSse.Max() / truths.Count), Math.Sqrt(sst / count), rd);
        return (fit, report);
    }

    public sealed record Job(int Rep, string Emulator, string Mode, int Noise, int Ratio, bool Discrepancy);

    /// <summary>Stage 1 (and Stage 2's shell): fit, held-out check, then the replicate grid and the modality subset.</summary>
    public static void Run(Settings s, Func<int, Synthesizer> make, int k, string dir, bool rerun, bool fitOnly, Func<Job, bool>? filter = null)
    {
        ThreadPool.SetMaxThreads(s.Threads, s.Threads);
        File.WriteAllText(Path.Combine(dir, $"{s.Stage}-settings.json"), JsonSerializer.Serialize(s, new JsonSerializerOptions { WriteIndented = true }));
        Synthesizer synth = make(k);
        int[] seeds = KSeeds(s, k);
        int off = rerun ? 8 : 0;
        double[] axis = synth.Frequencies.Select(f => Math.Log(f)).ToArray();
        IReadOnlyList<CalibrationParameter> parameters = Enumerable.Range(0, k).Select(j => (CalibrationParameter)new ContinuousParameter(synth.Names[j], synth.Lower[j], synth.Upper[j])).ToList();
        int n = s.TrainingBase + s.TrainingPerK * k;
        double[][] design = synth.LatinHypercube(n, new Xoshiro256PlusPlus(seeds[0]));
        var fits = new Dictionary<string, (Layout Layout, BassFit Fit)>();
        var reportLines = new List<string>();
        foreach (string layoutName in s.Layouts)
        {
            Layout layout = BuildLayout(layoutName, synth, design, parameters, axis);
            (BassFit fit, FitReport report) = FitAndCheck(s, synth, layout, parameters, axis, seeds[2], seeds[1], n);
            fits[layoutName] = (layout, fit);
            reportLines.Add(JsonSerializer.Serialize(report));
            Console.WriteLine(JsonSerializer.Serialize(report));
        }
        File.WriteAllLines(Path.Combine(dir, $"{s.Stage}-k{k}-fit{(rerun ? "-rerun" : "")}.jsonl"), reportLines);
        if (fitOnly) return;

        // Replicates: θ* uniform on the box, one noise realization per level, one random order of the bins (the kept set is
        // its first m), shared across ratios, emulators and modes so their comparisons are paired.
        var truthRng = new Xoshiro256PlusPlus(seeds[3 + off]);
        double[][] truths = Enumerable.Range(0, s.Replicates).Select(_ => synth.Draw(truthRng)).ToArray();
        int[] noiseSeeds = SeedTree.Derive(seeds[4 + off], s.Replicates * s.NoiseVariances.Length);
        double[][][] observed = truths.Select((t, r) => s.NoiseVariances.Select((v, i) => synth.Noisy(t, v, new Xoshiro256PlusPlus(noiseSeeds[r * s.NoiseVariances.Length + i]))).ToArray()).ToArray();
        double[][] clean = truths.Select(t => synth.Clean(t)).ToArray();
        var subsetRng = new Xoshiro256PlusPlus(seeds[5 + off]);
        int[][] orders = truths.Select(_ =>
        {
            int[] p = Enumerable.Range(0, synth.L).ToArray();
            for (int i = p.Length - 1; i > 0; i--) { int j = subsetRng.NextInt(i + 1); (p[i], p[j]) = (p[j], p[i]); }
            return p;
        }).ToArray();

        var jobs = new List<Job>();
        string[] emulators = [.. s.Layouts.Select(l => "bass-" + l), "oracle"];
        for (int r = 0; r < s.Replicates; r++)
            for (int v = 0; v < s.NoiseVariances.Length; v++)
                foreach (int ratio in s.Ratios)
                {
                    if (ratio * k > synth.L) continue;
                    foreach (string e in emulators) jobs.Add(new Job(r, e, "nested", v, ratio, false));
                    if (v == s.DiscrepancyNoise) foreach (string e in emulators.Where(e => e != "oracle")) jobs.Add(new Job(r, e, "nested", v, ratio, true));
                }
        if (s.ModalityKs.Contains(k))
            for (int r = 0; r < Math.Min(s.ModalityReplicates, s.Replicates); r++)
                foreach (int ratio in s.ModalityRatios)
                    foreach (string e in emulators)
                        foreach (string mode in new[] { "tempered", "smc" })
                            jobs.Add(new Job(r, e, mode, s.ModalityNoise, ratio, false));
        if (filter is not null) jobs = jobs.Where(filter).ToList();
        int[] jobSeeds = SeedTree.Derive(seeds[6 + off], jobs.Count);

        // The discrepancy basis over all bins: 1 and the standardized ln f.
        double mu = axis.Average(), sd = Math.Sqrt(axis.Average(a => (a - mu) * (a - mu)));
        double[][] discrepancyRows = axis.Select(a => new[] { 1.0, (a - mu) / sd }).ToArray();

        string file = Path.Combine(dir, $"{s.Stage}-k{k}-replicates{(rerun ? "-rerun" : "")}.csv");
        var lines = new ConcurrentBag<(int, string)>();
        string header = "stage,k,rep,emulator,mode,noise,ratio,m,disc,seed,snr_db,truth_unit,mean_unit,q05_unit,q95_unit,sd_unit,err_rms,covered," +
                        "converged,runs,sigma_post,recon_truth_held,recon_obs_held,fit_obs_kept,mass_truth,mass_double,mass_half,seconds,chain_spread";
        var emulatorCache = new ConcurrentDictionary<(string, int, int), ICalibrationEmulator>();
        int done = 0;
        var clock = Stopwatch.StartNew();
        Parallel.For(0, jobs.Count, new ParallelOptions { MaxDegreeOfParallelism = s.Threads }, ji =>
        {
            Job job = jobs[ji];
            int m = job.Ratio * k;
            int[] kept = orders[job.Rep][..m].OrderBy(i => i).ToArray();
            int[] held = orders[job.Rep][m..].OrderBy(i => i).ToArray();
            double[] keptAxis = kept.Select(i => axis[i]).ToArray();
            ICalibrationEmulator emulator = emulatorCache.GetOrAdd((job.Emulator, job.Rep, m), key => key.Item1 == "oracle"
                ? new SynthesizerEmulator(synth, parameters, kept.Select(i => synth.Frequencies[i]).ToArray())
                : EmulatorAt(fits[key.Item1[5..]].Layout, fits[key.Item1[5..]].Fit, parameters, keptAxis));
            double[] y = kept.Select(i => observed[job.Rep][job.Noise][i]).ToArray();
            var config = new CalibrationConfig
            {
                Parameters = parameters,
                Measurement = new GaussianMeasurement([new VarianceGroup(Enumerable.Range(0, m).ToArray(), new InverseGammaVariance(s.PriorDf, s.PriorSd))]),
                Covariance = job.Discrepancy
                    ? new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(kept.Select(i => discrepancyRows[i]).ToArray(), s.DiscrepancyScale) }
                    : new FunctionalFormulation(),
                Sampling = new CalibrationSampling { MasterSeed = jobSeeds[ji], BurnIn = s.CalBurnIn, DecorrelationEvery = s.CalDecorrelation, Tempering = job.Mode == "tempered" ? s.Ladder : null },
                Reference = new ReferenceSampling { OuterDraws = s.OuterDraws, Chains = s.CalChains, MaxSamples = s.MaxSamples, BatchSize = s.BatchSize, RHatTarget = s.RHat, MinEssPerChain = s.MinEss },
                Smc = job.Mode == "smc" ? new CutSmcSpec { Particles = s.SmcParticles, OuterDraws = s.OuterDraws } : null,
            };
            var model = new CalibrationModel(config, y, emulator);
            var sw = Stopwatch.StartNew();
            IReadOnlyList<CalibrationDraw> draws;
            double[] weights;
            int converged, runs;
            double spread = double.NaN;
            if (job.Mode == "smc")
            {
                CutSmcCalibration smc = CutSmcSampler.Run(model, 1);
                (draws, weights) = smc.Pooled();
                converged = runs = smc.Populations.Count;
            }
            else
            {
                NestedCalibration nested = NestedCalibrationSampler.Run(model);
                (draws, weights) = nested.Pooled(includeNotConverged: true);
                converged = nested.Runs.Count(x => x.Converged);
                runs = nested.Runs.Count;
                // The chains' disagreement: per inner run and coordinate, the range of the chains' means over the pooled sd.
                spread = 0;
                foreach (InnerRun run in nested.Runs)
                    for (int j = 0; j < k; j++)
                    {
                        double[][] chainValues = run.Chains.Select(c => c.Select(d => d.Theta[j]).ToArray()).ToArray();
                        double[] all = chainValues.SelectMany(v => v).ToArray();
                        double mean = all.Average(), sdAll = Math.Sqrt(all.Average(v => (v - mean) * (v - mean)));
                        double[] means = chainValues.Select(v => v.Average()).ToArray();
                        if (sdAll > 0) spread = Math.Max(spread, (means.Max() - means.Min()) / sdAll);
                    }
            }
            double seconds = sw.Elapsed.TotalSeconds;
            double[] truthUnit = Enumerable.Range(0, k).Select(j => (truths[job.Rep][j] - synth.Lower[j]) / (synth.Upper[j] - synth.Lower[j])).ToArray();
            var meanU = new double[k];
            var q05 = new double[k];
            var q95 = new double[k];
            var sdU = new double[k];
            int covered = 0;
            for (int j = 0; j < k; j++)
            {
                double[] u = draws.Select(d => (d.Theta[j] - synth.Lower[j]) / (synth.Upper[j] - synth.Lower[j])).ToArray();
                double wsum = weights.Sum();
                meanU[j] = u.Zip(weights, (a, b) => a * b).Sum() / wsum;
                sdU[j] = Math.Sqrt(u.Zip(weights, (a, b) => b * (a - meanU[j]) * (a - meanU[j])).Sum() / wsum);
                q05[j] = WeightedQuantile(u, weights, 0.05);
                q95[j] = WeightedQuantile(u, weights, 0.95);
                if (truthUnit[j] >= q05[j] && truthUnit[j] <= q95[j]) covered++;
            }
            double errRms = Math.Sqrt(Enumerable.Range(0, k).Average(j => (meanU[j] - truthUnit[j]) * (meanU[j] - truthUnit[j])));
            double sigmaPost = Math.Sqrt(draws.Zip(weights, (d, w) => d.Variances[0] * w).Sum() / weights.Sum());

            // Octave masses on the first coordinate (Stage 2's f₀): posterior mass within a semitone of θ*, 2θ* and θ*/2.
            double massTruth = 0, massDouble = 0, massHalf = 0;
            if (s.Stage == "s2")
            {
                double wsum = weights.Sum();
                for (int i = 0; i < draws.Count; i++)
                {
                    double octaves = (draws[i].Theta[0] - truths[job.Rep][0]) / Math.Log(2);
                    if (Math.Abs(octaves) < 1.0 / 12) massTruth += weights[i] / wsum;
                    if (Math.Abs(octaves - 1) < 1.0 / 12) massDouble += weights[i] / wsum;
                    if (Math.Abs(octaves + 1) < 1.0 / 12) massHalf += weights[i] / wsum;
                }
            }

            // Reconstruction over all bins from resampled posterior draws: the emulator draw each conditioned on, at its θ,
            // plus the discrepancy D·v where declared.
            List<CalibrationDraw> picked = Resample(draws, weights, s.ReconstructionDraws);
            var recon = new double[synth.L];
            var buf = new double[synth.L];
            ICalibrationEmulator fullEmulator = emulatorCache.GetOrAdd((job.Emulator, -1, synth.L), key => key.Item1 == "oracle"
                ? new SynthesizerEmulator(synth, parameters, synth.Frequencies)
                : EmulatorAt(fits[key.Item1[5..]].Layout, fits[key.Item1[5..]].Fit, parameters, axis));
            foreach (CalibrationDraw d in picked)
            {
                fullEmulator.Predict(d.EmulatorDraw, d.Theta, [], buf);
                for (int l = 0; l < synth.L; l++)
                {
                    double delta = job.Discrepancy ? discrepancyRows[l][0] * d.Coefficients[0] + discrepancyRows[l][1] * d.Coefficients[1] : 0.0;
                    recon[l] += (buf[l] + delta) / picked.Count;
                }
            }
            double Rmse(int[] idx, double[] reference) => idx.Length == 0 ? double.NaN : Math.Sqrt(idx.Average(i => (recon[i] - reference[i]) * (recon[i] - reference[i])));
            double snr = 10 * Math.Log10(synth.Energy(truths[job.Rep]) / (synth.Frame * s.NoiseVariances[job.Noise]));
            string line = string.Join(",", s.Stage, k, job.Rep, job.Emulator, job.Mode, F(s.NoiseVariances[job.Noise]), job.Ratio, m, job.Discrepancy ? 1 : 0, jobSeeds[ji], F(snr, "F2"),
                J(truthUnit), J(meanU), J(q05), J(q95), J(sdU), F(errRms, "F5"), covered, converged, runs, F(sigmaPost, "F4"),
                F(Rmse(held, clean[job.Rep]), "F4"), F(Rmse(held, observed[job.Rep][job.Noise]), "F4"), F(Rmse(kept, observed[job.Rep][job.Noise]), "F4"),
                F(massTruth, "F3"), F(massDouble, "F3"), F(massHalf, "F3"), F(seconds, "F2"), F(spread, "F3"));
            lines.Add((ji, line));
            int now = System.Threading.Interlocked.Increment(ref done);
            if (now % 50 == 0 || now == jobs.Count)
                Console.WriteLine($"[{s.Stage} k={k}] {now}/{jobs.Count} calibrations, {clock.Elapsed.TotalMinutes:F1} min");
        });
        File.WriteAllLines(file, [header, .. lines.OrderBy(t => t.Item1).Select(t => t.Item2)]);
        Console.WriteLine($"{file}: {jobs.Count} rows");
    }
}
