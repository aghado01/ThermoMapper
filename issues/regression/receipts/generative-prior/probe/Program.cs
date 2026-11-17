using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text.Json;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace GenerativePrior;

internal sealed record Settings(int CodeDimension, int HiddenWidth, int OutputWidth, int Measurements,
    int Repetitions, int MasterSeed, double Step, double GradientTolerance, int MaxUpdates,
    double SuccessTolerance, double ReferenceTolerance, bool SignCheck = true, bool FastForward = false,
    bool RequireAllRecovered = true, int? IntermediateWidth = null, bool Magnitude = false)
{
    public int Depth => IntermediateWidth.HasValue ? 3 : 2;
    public double NegativeRadius
    {
        get
        {
            if (Depth == 2 && !Magnitude) return 1 / Math.PI;
            double angle = Math.PI;
            for (int i = 0; i < Depth; i++)
                angle = Math.Acos(Math.Clamp(((Math.PI - angle) * Math.Cos(angle) + Math.Sin(angle)) / Math.PI, -1, 1));
            return Magnitude ? 2 * Math.Sin(angle) / Math.PI + (Math.PI - 2 * angle) / Math.PI * Math.Cos(angle) :
                Math.Cos(angle);
        }
    }
}

internal sealed record Trial(int Repetition, int Seed, int[] ChildSeeds, double[] Truth, double[] Start,
    double[] Estimate, int Updates, string Stop, int SignChecks, int[] FlipIterations,
    long SelectedPreactivationZeros, double InitialObjectivePlus, double InitialObjectiveMinus,
    double FinalObjective, double FinalGradientNorm, double RelativeCodeError, bool Success,
    double MaxReferenceDiscrepancy, int EvaluatedUpdates, int SkippedUpdates, int CycleLength,
    double CosineToTruth, double? NegativeAngleDegrees, double RadiusRatio, double AxialProjection,
    double PopulationRadiusDeviation, bool NegativeOrientation, bool PopulationRadiusMatch,
    string EndpointClass, bool FlatNegativeRay, double RayInnerProduct, double RaySquaredNorm,
    double? RayOptimalRadius, double? FiniteRayRadiusDeviation, long SelectedMeasurementZeros = 0,
    double PredictedNegativeRadius = 0.3183098861837907);

internal static class Program
{
    private static double Norm(double[] v) => Math.Sqrt(PairwiseSum.Dot(v, v));

    internal static double[] Code(int dimension, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return Enumerable.Range(0, dimension).Select(_ => Normal.SampleStandard(rng)).ToArray();
    }

    internal sealed record Endpoint(double[] X, int Updates, int Checks, string Stop, int[] Flips, long Zeros,
        double InitialPlus, double InitialMinus, double Objective, double GradientNorm, double Discrepancy,
        int SkippedUpdates, int CycleLength, long MeasurementZeros);

    private sealed record State(double[] X, int Updates, long Zeros, int Flips, long MeasurementZeros);

    private static bool SameBits(double[] a, double[] b)
    {
        for (int i = 0; i < a.Length; i++)
            if (BitConverter.DoubleToInt64Bits(a[i]) != BitConverter.DoubleToInt64Bits(b[i])) return false;
        return true;
    }

    internal static Endpoint Recover(Settings s, Network net, double[] y, double[] start)
    {
        double[] x = (double[])start.Clone();
        var negative = new double[x.Length];
        var gradient = new double[x.Length];
        var plus = new Forward(s);
        var minus = new Forward(s);
        var flips = new List<int>();
        var history = new Queue<State>();
        int checks = 0, updates = 0;
        int skipped = 0, cycleLength = 0;
        long zeros = 0;
        long measurementZeros = 0;
        double discrepancy = 0, initialPlus = 0, initialMinus = 0;
        string stop = "update-limit";

        while (updates < s.MaxUpdates)
        {
            if (s.FastForward)
            {
                if (history.Count == 128) history.Dequeue();
                history.Enqueue(new((double[])x.Clone(), updates, zeros, flips.Count, measurementZeros));
            }
            for (int i = 0; i < x.Length; i++) negative[i] = -x[i];
            net.Evaluate(x, y, plus);
            net.Evaluate(negative, y, minus);
            if (!double.IsFinite(plus.Objective) || !double.IsFinite(minus.Objective))
                throw new ArithmeticException("Nonfinite signed objective.");
            bool flip = s.SignCheck && minus.Objective < plus.Objective;
            checks++;
            if (flip) flips.Add(checks);
            var z = flip ? negative : x;
            var state = flip ? minus : plus;
            zeros += state.Zeros;
            measurementZeros += state.MeasurementZeros;
            net.Gradient(state, gradient);

            if (checks == 1)
            {
                initialPlus = plus.Objective;
                initialMinus = minus.Objective;
                discrepancy = Math.Max(Reference.Check(net, x, y, initialPlus, flip ? null : gradient),
                    Reference.Check(net, negative, y, initialMinus, flip ? gradient : null));
            }

            double norm = Norm(gradient);
            if (!double.IsFinite(norm)) throw new ArithmeticException("Nonfinite gradient.");
            if (norm < s.GradientTolerance)
            {
                Array.Copy(z, x, x.Length);
                stop = "gradient";
                break;
            }

            for (int j = 0; j < x.Length; j++) x[j] = z[j] - s.Step * gradient[j];
            updates++;

            if (s.FastForward)
            {
                State? repeated = history.LastOrDefault(state => SameBits(state.X, x));
                if (repeated is not null)
                {
                    int period = updates - repeated.Updates;
                    int copies = (s.MaxUpdates - updates) / period;
                    if (copies > 0)
                    {
                        long cycleZeros = zeros - repeated.Zeros;
                        long cycleMeasurementZeros = measurementZeros - repeated.MeasurementZeros;
                        int[] cycleFlips = flips.Skip(repeated.Flips).Select(i => i - repeated.Updates).ToArray();
                        for (int copy = 0; copy < copies; copy++)
                            foreach (int offset in cycleFlips) flips.Add(updates + copy * period + offset);
                        int count = copies * period;
                        zeros += copies * cycleZeros;
                        measurementZeros += copies * cycleMeasurementZeros;
                        updates += count;
                        checks += count;
                        skipped += count;
                        cycleLength = period;
                        history.Clear();
                    }
                }
            }
        }

        net.Evaluate(x, y, plus);
        net.Gradient(plus, gradient);
        for (int i = 0; i < x.Length; i++) negative[i] = -x[i];
        net.Evaluate(negative, y, minus);
        discrepancy = Math.Max(discrepancy, Reference.Check(net, x, y, plus.Objective, gradient));
        discrepancy = Math.Max(discrepancy, Reference.Check(net, negative, y, minus.Objective));
        return new(x, updates, checks, stop, flips.ToArray(), zeros, initialPlus, initialMinus,
            plus.Objective, Norm(gradient), discrepancy, skipped, cycleLength, measurementZeros);
    }

    private static Trial Run(Settings s, int repetition, int seed)
    {
        var children = SeedTree.Derive(seed, s.IntermediateWidth.HasValue ? 6 : 5);
        var net = new Network(s, children);
        var truth = Code(s.CodeDimension, children[3]);
        var start = Code(s.CodeDimension, children[4]);
        if (Norm(start) == 0 || Norm(truth) == 0) throw new ArithmeticException("Zero code.");
        var observation = new Forward(s);
        net.Predict(truth, observation);
        double[] target = s.Magnitude ? observation.Residual.Select(Math.Abs).ToArray() : observation.Residual;
        var oppositeTruth = truth.Select(v => -v).ToArray();
        var ray = new Forward(s);
        net.Predict(oppositeTruth, ray);
        double rayNorm = PairwiseSum.Dot(ray.Residual, ray.Residual);
        double[] rayResponse = s.Magnitude ? ray.Residual.Select(Math.Abs).ToArray() : ray.Residual;
        double inner = PairwiseSum.Dot(rayResponse, target);
        double? rayRadius = rayNorm == 0 ? null : Math.Max(0, inner / rayNorm);
        Endpoint end = Recover(s, net, target, start);
        double error = Norm(end.X.Zip(truth, (a, b) => a - b).ToArray()) / Norm(truth);
        double truthNorm = Norm(truth), endpointNorm = Norm(end.X);
        double cosine = endpointNorm == 0 ? 0 : Math.Clamp(PairwiseSum.Dot(end.X, truth) / (endpointNorm * truthNorm), -1, 1);
        double? angle = endpointNorm == 0 ? null : Math.Acos(-cosine) * 180 / Math.PI;
        double radius = endpointNorm / truthNorm, rho = s.NegativeRadius;
        bool orientation = angle.HasValue && angle.Value <= 20;
        bool radiusMatch = Math.Abs(radius - rho) <= rho / 4;
        string classification = error < s.SuccessTolerance ? "recovered" :
            orientation && radiusMatch ? "negative-neighbourhood" : "elsewhere";
        return new(repetition, seed, children, truth, start, end.X, end.Updates, end.Stop, end.Checks,
            end.Flips, end.Zeros, end.InitialPlus, end.InitialMinus, end.Objective, end.GradientNorm,
            error, error < s.SuccessTolerance, end.Discrepancy, end.Updates - end.SkippedUpdates,
            end.SkippedUpdates, end.CycleLength, cosine, angle, radius,
            PairwiseSum.Dot(end.X, truth) / (truthNorm * truthNorm), radius - rho,
            orientation, radiusMatch, classification, rayNorm == 0, inner, rayNorm,
            rayRadius, rayRadius.HasValue ? radius - rayRadius.Value : null, end.MeasurementZeros, rho);
    }

    public static int Main(string[] args)
    {
        if (args.Length < 3) throw new ArgumentException("Expected settings.json, output directory, source revision, optional trial indices.");
        using var process = Process.GetCurrentProcess();
        process.PriorityClass = ProcessPriorityClass.BelowNormal;
        if (process.PriorityClass != ProcessPriorityClass.BelowNormal || Environment.ProcessorCount != 1)
            throw new InvalidOperationException("Run at below-normal priority with DOTNET_PROCESSOR_COUNT=1.");

        byte[] settingsBytes = File.ReadAllBytes(args[0]);
        if (JsonDocument.Parse(settingsBytes).RootElement.TryGetProperty("ExtendedCells", out _))
            return Extensions.Run(args, settingsBytes);
        if (JsonDocument.Parse(settingsBytes).RootElement.TryGetProperty("MeasurementCounts", out _))
            return Measurements.Run(args, settingsBytes);
        if (JsonDocument.Parse(settingsBytes).RootElement.TryGetProperty("Cells", out _))
            return Curve.Run(args, settingsBytes);
        var settings = JsonSerializer.Deserialize<Settings>(settingsBytes)!;
        return Execute(settings, settingsBytes, args[1], args[2], args.Length == 4 ? args[3] : null);
    }

    internal static int Execute(Settings settings, byte[] settingsBytes, string results, string revision, string? selected)
    {
        using var process = Process.GetCurrentProcess();
        Directory.CreateDirectory(results);
        int[] seeds = SeedTree.Derive(settings.MasterSeed, settings.Repetitions);
        int[] indices = selected is not null ? selected.Split(',').Select(int.Parse).Select(i => i - 1).ToArray() :
            Enumerable.Range(0, seeds.Length).ToArray();
        var trials = new List<Trial>(seeds.Length);
        using (var rows = new StreamWriter(Path.Combine(results, "trials.jsonl")))
        {
            rows.AutoFlush = true;
            foreach (int i in indices)
            {
                Trial trial = Run(settings, i + 1, seeds[i]);
                trials.Add(trial);
                rows.WriteLine(JsonSerializer.Serialize(trial));
                Console.WriteLine($"{i + 1}/{seeds.Length} {trial.Stop} updates={trial.Updates} relative-error={trial.RelativeCodeError:G6}");
            }
        }

        int successes = trials.Count(t => t.Success);
        long zeros = trials.Sum(t => t.SelectedPreactivationZeros);
        double maxDiscrepancy = trials.Max(t => t.MaxReferenceDiscrepancy);
        bool gate = zeros == 0 && trials.Sum(t => t.SelectedMeasurementZeros) == 0 && maxDiscrepancy <= settings.ReferenceTolerance &&
            (!settings.RequireAllRecovered || successes == trials.Count);
        var summary = new
        {
            SourceRevision = revision, SettingsSha256 = Convert.ToHexStringLower(SHA256.HashData(settingsBytes)),
            Runtime = RuntimeInformation.FrameworkDescription, Architecture = RuntimeInformation.ProcessArchitecture.ToString(),
            Priority = process.PriorityClass.ToString(), ProcessorCount = Environment.ProcessorCount,
            settings.CodeDimension, settings.Measurements, settings.MasterSeed, settings.SignCheck, settings.FastForward,
            settings.IntermediateWidth, settings.Magnitude, settings.Depth, settings.NegativeRadius,
            Repetitions = trials.Count, Successes = successes, GradientStops = trials.Count(t => t.Stop == "gradient"),
            UpdateLimitStops = trials.Count(t => t.Stop == "update-limit"),
            MinUpdates = trials.Min(t => t.Updates), MaxUpdates = trials.Max(t => t.Updates),
            SignChecks = trials.Sum(t => t.SignChecks), Flips = trials.Sum(t => t.FlipIterations.Length),
            RepetitionsWithFlips = trials.Count(t => t.FlipIterations.Length > 0),
            SelectedPreactivationZeros = zeros, MinRelativeCodeError = trials.Min(t => t.RelativeCodeError),
            SelectedMeasurementZeros = trials.Sum(t => t.SelectedMeasurementZeros),
            MaxRelativeCodeError = trials.Max(t => t.RelativeCodeError), MaxReferenceDiscrepancy = maxDiscrepancy,
            EvaluatedUpdates = trials.Sum(t => t.EvaluatedUpdates), SkippedUpdates = trials.Sum(t => t.SkippedUpdates),
            CycleLengths = trials.Where(t => t.CycleLength > 0).Select(t => t.CycleLength).Distinct().Order().ToArray(),
            EndpointClasses = trials.GroupBy(t => t.EndpointClass).ToDictionary(g => g.Key, g => g.Count()),
            GatePassed = gate
        };
        File.WriteAllText(Path.Combine(results, "summary.json"),
            JsonSerializer.Serialize(summary, new JsonSerializerOptions { WriteIndented = true }) + "\n");
        return gate ? 0 : 1;
    }
}
