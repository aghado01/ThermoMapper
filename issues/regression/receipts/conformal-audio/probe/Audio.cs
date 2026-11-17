using System.Diagnostics;
using System.Text.Json;
using LatentRecoveryProbe;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;

namespace ConformalAudio;

internal sealed record Turn(double[] Short, double[] Plus, double[] Target);
internal sealed record Conversation(int Seed, int Mark, double NoiseSd, Turn[] Turns);
internal sealed record ScoredConversation(Conversation Conversation, double[][] Centers, double[] Scores);

internal sealed class FrozenPredictor(BassFit fit, BassData data, int features, int[] draws, double scale)
{
    internal double Scale => scale;
    internal ScoredConversation[] Score(Conversation[] groups)
    {
        Turn[] turns = groups.SelectMany(c => c.Turns).ToArray();
        BassPrediction prediction = fit.Predict(data, Audio.Cases(turns, features), new BassPredictOptions { Draws = draws, Predictive = false });
        double[,] mean = prediction.Mean(); int at = 0, bins = Audio.TargetHz.Length;
        return groups.Select(c => {
            var centers = new double[c.Turns.Length][]; var scores = new double[c.Turns.Length];
            for (int t = 0; t < c.Turns.Length; t++)
            {
                centers[t] = new double[bins];
                for (int b = 0; b < bins; b++) { centers[t][b] = mean[at++, 0]; scores[t] = Math.Max(scores[t], Math.Abs(c.Turns[t].Target[b] - centers[t][b]) / scale); }
                App.Require(double.IsFinite(scores[t]), "Finite audio score.");
            }
            return new ScoredConversation(c, centers, scores);
        }).ToArray();
    }
}

internal static class Audio
{
    internal static double[] TargetHz = [];
    private static JsonElement config;
    private static readonly ResonanceSynthesizer resonance = new(2);
    private sealed class Transform(int width, double[] frequencies, double rate)
    {
        private readonly double[] taper = Windows.Hann(width);
        private readonly double[][] cos = frequencies.Select(f => Enumerable.Range(0, width).Select(n => Math.Cos(2 * Math.PI * f * n / rate)).ToArray()).ToArray();
        private readonly double[][] sin = frequencies.Select(f => Enumerable.Range(0, width).Select(n => Math.Sin(2 * Math.PI * f * n / rate)).ToArray()).ToArray();
        internal double[] Power(double[] signal, int start)
        {
            var powers = new double[frequencies.Length];
            for (int k = 0; k < powers.Length; k++)
            {
                double re = 0, im = 0;
                for (int j = 0; j < width; j++) { double x = signal[start + j] * taper[j]; re += x * cos[k][j]; im += x * sin[k][j]; }
                powers[k] = re * re + im * im;
            }
            return powers;
        }
    }
    private static Transform shortTransform = null!, longTransform = null!, targetTransform = null!;
    private static double[] Features(double[] waveform)
    {
        int width = App.I(config, "ShortWindow"), hop = App.I(config, "Hop"), count = (waveform.Length - width) / hop + 1;
        int bins = App.Doubles(config, "ShortHz").Length;
        var sum = new double[bins]; var first = new double[bins]; var second = new double[bins];
        for (int f = 0; f < count; f++)
        {
            double location = (f * hop + 0.5 * width) / waveform.Length;
            double[] p = shortTransform.Power(waveform, f * hop);
            for (int b = 0; b < bins; b++) { sum[b] += p[b]; first[b] += p[b] * location; second[b] += p[b] * location * location; }
        }
        var result = new double[3 * bins];
        for (int b = 0; b < bins; b++)
        {
            double denominator = Math.Max(1e-12, sum[b]), center = first[b] / denominator;
            result[3 * b] = Math.Log(sum[b] / count + App.D(config, "Floor")); result[3 * b + 1] = center;
            result[3 * b + 2] = Math.Sqrt(Math.Max(0, second[b] / denominator - center * center));
        }
        return result;
    }
    private static Conversation[] Generate(int seed, int count, int gap)
    {
        int[] seeds = SeedTree.Derive(seed, count); int frame = App.I(config, "Frame");
        return seeds.Select(s => {
            int[] branches = SeedTree.Derive(s, 3); var source = new Xoshiro256PlusPlus(branches[0]);
            double amplitude = 0.6 + 0.8 * source.NextDouble(); int[] lengths = App.Ints(config, "Lengths");
            int length = lengths[source.NextInt(lengths.Length)]; double[] noiseLevels = App.Doubles(config, "NoiseSd");
            double noise = noiseLevels[source.NextInt(noiseLevels.Length)]; double[] code = resonance.Draw(source);
            int[] turnSeeds = SeedTree.Derive(branches[1], length); var turns = new Turn[length];
            for (int t = 0; t < length; t++)
            {
                int[] child = SeedTree.Derive(turnSeeds[t], 2); var rng = new Xoshiro256PlusPlus(child[0]); var noiseRng = new Xoshiro256PlusPlus(child[1]);
                double gain = 0.9 + 0.2 * rng.NextDouble(); int sign = rng.NextInt(2) == 0 ? -1 : 1, jitter = rng.NextInt(17) - 8;
                var waveform = new double[frame];
                if (gap > 0) { int first = (frame - gap) / 2 + jitter; waveform[first] = amplitude * gain; waveform[first + gap] = 0.8 * amplitude * gain * sign; }
                else { double[] clean = resonance.Waveform(code); for (int j = 0; j < clean.Length; j++) waveform[(frame - clean.Length) / 2 + jitter + j] = amplitude * gain * clean[j]; }
                for (int j = 0; j < frame; j++) waveform[j] += noise * Normal.SampleStandard(noiseRng);
                double[] small = Features(waveform), large = longTransform.Power(waveform, (frame - App.I(config, "LongWindow")) / 2)
                    .Select(p => Math.Log(p + App.D(config, "Floor"))).ToArray();
                double[] target = targetTransform.Power(waveform, 0).Select(p => Math.Log(p + App.D(config, "Floor"))).ToArray();
                turns[t] = new(small, small.Concat(large).ToArray(), target);
            }
            return new Conversation(s, new Xoshiro256PlusPlus(branches[2]).NextInt(length), noise, turns);
        }).ToArray();
    }
    private static double[] X(Turn t, int features) => features == t.Short.Length ? t.Short : t.Plus;
    internal static PredictionCases Cases(Turn[] turns, int features)
    {
        var columns = new List<PredictionColumn>();
        for (int j = 0; j < features; j++) { int column = j; columns.Add(new ContinuousColumn($"v{j}", turns.SelectMany(t => TargetHz.Select(_ => X(t, features)[column])).ToArray())); }
        columns.Add(new ContinuousColumn("frequency", turns.SelectMany(_ => TargetHz).ToArray()));
        return new(columns);
    }
    private static FrozenPredictor Fit(Conversation[] training, int features, int seed, string name)
    {
        Turn[] turns = training.SelectMany(c => c.Turns).ToArray();
        var predictors = new List<BassPredictor>();
        for (int j = 0; j < features; j++) { int column = j; predictors.Add(new ContinuousPredictor($"v{j}", turns.SelectMany(t => TargetHz.Select(_ => X(t, features)[column])).ToArray())); }
        predictors.Add(new ContinuousPredictor("frequency", turns.SelectMany(_ => TargetHz).ToArray()));
        var data = new BassData(predictors, new ScalarResponse(turns.SelectMany(t => t.Target).ToArray()));
        var settings = new BassConfig {
            Structure = new BassStructureSpec { MaxBasis = App.I(config, "MaxBasis"), MaxInteraction = App.I(config, "MaxInteraction") },
            Sampling = new BassSampling { Chains = App.I(config, "Chains"), MasterSeed = seed, BurnIn = App.I(config, "BurnIn"), Samples = App.I(config, "Samples"), Thin = App.I(config, "Thin") }
        };
        var clock = Stopwatch.StartNew(); BassFit fit = Bass.Fit(settings, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        double seconds = clock.Elapsed.TotalSeconds; int total = fit.Draws.Sum(d => d.Count), used = Math.Min(total, App.I(config, "PredictionDraws"));
        int[] draws = Enumerable.Range(0, used).Select(i => i * total / used).ToArray();
        BassPrediction first = fit.Predict(data, Cases([turns[0]], features), new BassPredictOptions { Draws = draws, Predictive = false });
        double scale = Math.Max(1e-6, Math.Sqrt(first.NoiseVariances.Select((v, i) => v * first.Weights[i]).Sum()));
        App.Require(double.IsFinite(scale) && scale > 0 && total > 0, "Frozen fit scale/draws.");
        App.Write($"fit-{name}.json", new { Name = name, Seed = seed, Settings = settings, TrainingHash = App.Hash(training),
            Identity = fit.Identity, TotalDraws = total, SelectedDraws = draws, Scale = scale, Seconds = seconds,
            Verdict = fit.Diagnostics.Verdict.Outcome.ToString(), Groups = fit.Diagnostics.Verdict.Groups.Select(g => new { g.Name, Outcome = g.Outcome.ToString(),
                Failing = g.Failing.Select(f => new { f.Name, RHat = f.RHat.ToString(), Ess = f.Ess?.ToString() }).ToArray() }).ToArray(),
            fit.Diagnostics.AcceptanceRate, fit.Diagnostics.MassAtCap });
        return new(fit, data, features, draws, scale);
    }
    private sealed record SplitReading(int Gap, int Split, string Predictor, string GroupEvent, double Threshold, double Width,
        double MarkedCoverage, double AllCoverage, double Rmse, int CalibrationGroups, int TestGroups, bool FullSpace, string PredictionHash);
    internal static void Run(int seed, bool firstSplit)
    {
        config = App.Settings.GetProperty("Audio"); TargetHz = App.Doubles(config, "TargetHz");
        int frame = App.I(config, "Frame"); double rate = App.D(config, "Rate");
        shortTransform = new(App.I(config, "ShortWindow"), App.Doubles(config, "ShortHz"), rate);
        longTransform = new(App.I(config, "LongWindow"), App.Doubles(config, "LongHz"), rate); targetTransform = new(frame, TargetHz, rate);
        var readings = new List<SplitReading>();
        using var splitFile = new StreamWriter(Path.Combine(App.Output, "audio-splits.jsonl")) { AutoFlush = true };
        using var testFile = new StreamWriter(Path.Combine(App.Output, "audio-tests.jsonl"));
        using var predictionFile = new StreamWriter(Path.Combine(App.Output, "audio-predictions.jsonl"));
        int[] gaps = App.Ints(config, "Gaps").Concat(new[] { 0 }).ToArray(), geometrySeeds = SeedTree.Derive(seed, gaps.Length);
        for (int g = 0; g < gaps.Length; g++)
        {
            if (firstSplit && g > 0) break;
            int gap = gaps[g]; int[] branches = SeedTree.Derive(geometrySeeds[g], 4);
            Conversation[] training = Generate(branches[0], App.I(config, "TrainingGroups"), gap); App.Write($"training-gap{gap}.json", training);
            FrozenPredictor shortFit = Fit(training, training[0].Turns[0].Short.Length, branches[1], $"gap{gap}-short");
            FrozenPredictor plusFit = Fit(training, training[0].Turns[0].Plus.Length, branches[2], $"gap{gap}-plus");
            double[] baselineCenter = Enumerable.Range(0, TargetHz.Length).Select(b => training.SelectMany(c => c.Turns).Average(t => t.Target[b])).ToArray();
            int[] splits = SeedTree.Derive(branches[3], App.I(config, "Splits"));
            for (int split = 0; split < splits.Length; split++)
            {
                if (firstSplit && split > 0) break;
                int[] children = SeedTree.Derive(splits[split], 2);
                Conversation[] calibration = Generate(children[0], App.I(config, "CalibrationGroups"), gap), test = Generate(children[1], App.I(config, "TestGroups"), gap);
                ScoredConversation[] Base(Conversation[] groups) => groups.Select(c => new ScoredConversation(c, c.Turns.Select(_ => baselineCenter).ToArray(),
                    c.Turns.Select(t => t.Target.Zip(baselineCenter, (y, m) => Math.Abs(y - m)).Max()).ToArray())).ToArray();
                foreach (var (name, cal, tst, scale) in new[] {
                    ("short", shortFit.Score(calibration), shortFit.Score(test), shortFit.Scale),
                    ("plus", plusFit.Score(calibration), plusFit.Score(test), plusFit.Scale),
                    ("baseline", Base(calibration), Base(test), 1.0) })
                {
                    for (int c = 0; c < tst.Length; c++) App.Line(predictionFile, new { Gap = gap, Split = split + 1, Predictor = name,
                        tst[c].Conversation.Seed, tst[c].Conversation.Mark, tst[c].Conversation.NoiseSd, Targets = tst[c].Conversation.Turns.Select(t => t.Target).ToArray(), tst[c].Centers, tst[c].Scores, Scale = scale });
                    string hash = App.Hash(tst);
                    foreach (GroupScore groupEvent in new[] { GroupScore.OneCase, GroupScore.WorstCase })
                    {
                        double[] scores = cal.Select(c => -(groupEvent == GroupScore.OneCase ? c.Scores[c.Conversation.Mark] : c.Scores.Max())).ToArray();
                        SplitCalibration fit = SplitConformal.Calibrate(scores, alpha: App.D(App.Settings, "Alpha")); double q = -fit.Threshold(0);
                        App.Require(double.IsFinite(q) && q >= 0, "Finite nonnegative audio threshold.");
                        double[] flatCal = groupEvent == GroupScore.OneCase ? scores : cal.SelectMany(c => c.Scores.Select(s => -s)).ToArray();
                        int[] groupOf = groupEvent == GroupScore.OneCase ? Enumerable.Range(0, cal.Length).ToArray() : cal.SelectMany((c, i) => Enumerable.Repeat(i, c.Scores.Length)).ToArray();
                        double[] flatTest = tst.SelectMany(c => c.Scores.Select(s => -s)).ToArray();
                        var engine = GroupConformal.Density(App.Seam(flatCal), groupOf, App.Seam(flatTest), alpha: App.D(App.Settings, "Alpha"), score: groupEvent);
                        int at = 0, marked = 0, all = 0; double squared = 0; int values = 0;
                        for (int c = 0; c < tst.Length; c++)
                        {
                            bool mc = tst[c].Scores[tst[c].Conversation.Mark] <= q, ac = tst[c].Scores.Max() <= q;
                            for (int t = 0; t < tst[c].Scores.Length; t++) { App.Require(engine.Sets[at++].AcceptsPlain(0) == (tst[c].Scores[t] <= q), "Public group rank agreement.");
                                for (int b = 0; b < TargetHz.Length; b++) { double e = tst[c].Conversation.Turns[t].Target[b] - tst[c].Centers[t][b]; squared += e * e; values++; } }
                            marked += mc ? 1 : 0; all += ac ? 1 : 0;
                            App.Line(testFile, new { Gap = gap, Split = split + 1, Predictor = name, GroupEvent = groupEvent.ToString(), Conversation = c,
                                Length = tst[c].Scores.Length, tst[c].Conversation.NoiseSd, MarkedCovered = mc, AllCovered = ac, Width = 2 * q * scale });
                        }
                        var row = new SplitReading(gap, split + 1, name, groupEvent.ToString(), q, 2 * q * scale, marked / (double)tst.Length,
                            all / (double)tst.Length, Math.Sqrt(squared / values), cal.Length, tst.Length, fit.Level(0) == 0, hash);
                        readings.Add(row); App.Line(splitFile, row);
                    }
                }
                Console.WriteLine($"Audio gap={gap}, split {split + 1}/{splits.Length} complete.");
            }
        }
        var summaries = readings.GroupBy(r => new { r.Gap, r.Predictor, r.GroupEvent }).Select(g => new { g.Key.Gap, g.Key.Predictor, g.Key.GroupEvent,
            MarkedCoverage = App.MeanInterval(g.Select(r => r.MarkedCoverage).ToArray()), AllCoverage = App.MeanInterval(g.Select(r => r.AllCoverage).ToArray()),
            Width = App.MeanInterval(g.Select(r => r.Width).ToArray()), MeanRmse = g.Average(r => r.Rmse), FullSpaceFrequency = g.Count(r => r.FullSpace) / (double)g.Count() }).ToArray();
        var changes = readings.Where(r => r.Predictor == "short").Select(r => new { r.Gap, r.GroupEvent, r.Split,
            Change = readings.Single(p => p.Gap == r.Gap && p.Split == r.Split && p.GroupEvent == r.GroupEvent && p.Predictor == "plus").Width - r.Width });
        App.Write("audio-summary.json", new { App.Revision, App.SettingsHash, NumericalControlsPassed = true, FirstSplitOnly = firstSplit,
            Groups = summaries, PairedWidthChanges = changes.GroupBy(c => new { c.Gap, c.GroupEvent }).Select(g => new { g.Key.Gap, g.Key.GroupEvent, Change = App.MeanInterval(g.Select(c => c.Change).ToArray()) }).ToArray() });
    }
}
