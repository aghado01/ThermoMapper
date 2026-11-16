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
    double SuccessTolerance, double ReferenceTolerance);

internal sealed record Trial(int Repetition, int Seed, int[] ChildSeeds, double[] Truth, double[] Start,
    double[] Estimate, int Updates, string Stop, int SignChecks, int[] FlipIterations,
    long SelectedPreactivationZeros, double InitialObjectivePlus, double InitialObjectiveMinus,
    double FinalObjective, double FinalGradientNorm, double RelativeCodeError, bool Success,
    double MaxReferenceDiscrepancy);

internal static class Program
{
    private static double Norm(double[] v) => Math.Sqrt(PairwiseSum.Dot(v, v));

    private static double[] Code(int dimension, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return Enumerable.Range(0, dimension).Select(_ => Normal.SampleStandard(rng)).ToArray();
    }

    private sealed record Endpoint(double[] X, int Updates, int Checks, string Stop, int[] Flips, long Zeros,
        double InitialPlus, double InitialMinus, double Objective, double GradientNorm, double Discrepancy);

    private static Endpoint Recover(Settings s, Network net, double[] y, double[] start)
    {
        double[] x = (double[])start.Clone();
        var negative = new double[x.Length];
        var gradient = new double[x.Length];
        var plus = new Forward(s);
        var minus = new Forward(s);
        var flips = new List<int>();
        int checks = 0, updates = 0;
        long zeros = 0;
        double discrepancy = 0, initialPlus = 0, initialMinus = 0;
        string stop = "update-limit";

        while (updates < s.MaxUpdates)
        {
            for (int i = 0; i < x.Length; i++) negative[i] = -x[i];
            net.Evaluate(x, y, plus);
            net.Evaluate(negative, y, minus);
            if (!double.IsFinite(plus.Objective) || !double.IsFinite(minus.Objective))
                throw new ArithmeticException("Nonfinite signed objective.");
            bool flip = minus.Objective < plus.Objective;
            checks++;
            if (flip) flips.Add(checks);
            var z = flip ? negative : x;
            var state = flip ? minus : plus;
            zeros += state.Zeros;
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
        }

        net.Evaluate(x, y, plus);
        net.Gradient(plus, gradient);
        for (int i = 0; i < x.Length; i++) negative[i] = -x[i];
        net.Evaluate(negative, y, minus);
        discrepancy = Math.Max(discrepancy, Reference.Check(net, x, y, plus.Objective, gradient));
        discrepancy = Math.Max(discrepancy, Reference.Check(net, negative, y, minus.Objective));
        return new(x, updates, checks, stop, flips.ToArray(), zeros, initialPlus, initialMinus,
            plus.Objective, Norm(gradient), discrepancy);
    }

    private static Trial Run(Settings s, int repetition, int seed)
    {
        var children = SeedTree.Derive(seed, 5);
        var net = new Network(s, children);
        var truth = Code(s.CodeDimension, children[3]);
        var start = Code(s.CodeDimension, children[4]);
        if (Norm(start) == 0 || Norm(truth) == 0) throw new ArithmeticException("Zero code.");
        var observation = new Forward(s);
        net.Predict(truth, observation);
        Endpoint end = Recover(s, net, observation.Residual, start);
        double error = Norm(end.X.Zip(truth, (a, b) => a - b).ToArray()) / Norm(truth);
        return new(repetition, seed, children, truth, start, end.X, end.Updates, end.Stop, end.Checks,
            end.Flips, end.Zeros, end.InitialPlus, end.InitialMinus, end.Objective, end.GradientNorm,
            error, error < s.SuccessTolerance, end.Discrepancy);
    }

    public static int Main(string[] args)
    {
        if (args.Length != 3) throw new ArgumentException("Expected settings.json, output directory, source revision.");
        using var process = Process.GetCurrentProcess();
        process.PriorityClass = ProcessPriorityClass.BelowNormal;
        if (process.PriorityClass != ProcessPriorityClass.BelowNormal || Environment.ProcessorCount != 1)
            throw new InvalidOperationException("Run at below-normal priority with DOTNET_PROCESSOR_COUNT=1.");

        byte[] settingsBytes = File.ReadAllBytes(args[0]);
        var settings = JsonSerializer.Deserialize<Settings>(settingsBytes)!;
        string results = args[1];
        Directory.CreateDirectory(results);
        int[] seeds = SeedTree.Derive(settings.MasterSeed, settings.Repetitions);
        var trials = new List<Trial>(seeds.Length);
        using (var rows = new StreamWriter(Path.Combine(results, "trials.jsonl")))
        {
            rows.AutoFlush = true;
            for (int i = 0; i < seeds.Length; i++)
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
        bool gate = successes == trials.Count && zeros == 0 && maxDiscrepancy <= settings.ReferenceTolerance;
        var summary = new
        {
            SourceRevision = args[2], SettingsSha256 = Convert.ToHexStringLower(SHA256.HashData(settingsBytes)),
            Runtime = RuntimeInformation.FrameworkDescription, Architecture = RuntimeInformation.ProcessArchitecture.ToString(),
            Priority = process.PriorityClass.ToString(), ProcessorCount = Environment.ProcessorCount,
            Repetitions = seeds.Length, Successes = successes, GradientStops = trials.Count(t => t.Stop == "gradient"),
            UpdateLimitStops = trials.Count(t => t.Stop == "update-limit"),
            MinUpdates = trials.Min(t => t.Updates), MaxUpdates = trials.Max(t => t.Updates),
            SignChecks = trials.Sum(t => t.SignChecks), Flips = trials.Sum(t => t.FlipIterations.Length),
            RepetitionsWithFlips = trials.Count(t => t.FlipIterations.Length > 0),
            SelectedPreactivationZeros = zeros, MinRelativeCodeError = trials.Min(t => t.RelativeCodeError),
            MaxRelativeCodeError = trials.Max(t => t.RelativeCodeError), MaxReferenceDiscrepancy = maxDiscrepancy,
            GatePassed = gate
        };
        File.WriteAllText(Path.Combine(results, "summary.json"),
            JsonSerializer.Serialize(summary, new JsonSerializerOptions { WriteIndented = true }) + "\n");
        return gate ? 0 : 1;
    }
}
