using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text.Json;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace GenerativePrior;

internal sealed record MeasurementSettings(Settings Solver, int[] MeasurementCounts, string BaselineReceipt);

internal static class Measurements
{
    private static bool Same(double x, double y) =>
        BitConverter.DoubleToInt64Bits(x) == BitConverter.DoubleToInt64Bits(y);

    private static string Hash(double[][] rows)
    {
        if (!BitConverter.IsLittleEndian) throw new PlatformNotSupportedException("Hashes use little-endian binary64.");
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        foreach (double[] row in rows) hash.AppendData(MemoryMarshal.AsBytes(row.AsSpan()));
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    private static object CheckData(MeasurementSettings config, Trial baseline)
    {
        Settings s = config.Solver;
        int[] children = SeedTree.Derive(baseline.Seed, 5);
        double[] truth = Program.Code(s.CodeDimension, children[3]);
        double[] start = Program.Code(s.CodeDimension, children[4]);
        if (!children.SequenceEqual(baseline.ChildSeeds) ||
            !truth.Zip(baseline.Truth, Same).All(v => v) || !start.Zip(baseline.Start, Same).All(v => v))
            throw new InvalidOperationException("Baseline streams or vectors differ.");
        var net = new Network(s, children);
        var observation = new Forward(s);
        net.Predict(truth, observation);
        var plus = new Forward(s);
        var minus = new Forward(s);
        net.Evaluate(start, observation.Residual, plus);
        net.Evaluate(start.Select(v => -v).ToArray(), observation.Residual, minus);
        if (!Same(plus.Objective, baseline.InitialObjectivePlus) || !Same(minus.Objective, baseline.InitialObjectiveMinus))
            throw new InvalidOperationException("Baseline initial objectives differ.");
        var gradient = new double[s.CodeDimension];
        Forward chosen = minus.Objective < plus.Objective ? minus : plus;
        double[] initial = minus.Objective < plus.Objective ? start.Select(v => -v).ToArray() : start;
        net.Gradient(chosen, gradient);
        double reference = Reference.Check(net, initial, observation.Residual, chosen.Objective, gradient);
        net.Evaluate(baseline.Estimate, observation.Residual, plus);
        net.Gradient(plus, gradient);
        double norm = Math.Sqrt(PairwiseSum.Dot(gradient, gradient));
        reference = Math.Max(reference, Reference.Check(net, baseline.Estimate, observation.Residual, plus.Objective, gradient));
        var ray = new Forward(s);
        net.Predict(truth.Select(v => -v).ToArray(), ray);
        if (!Same(plus.Objective, baseline.FinalObjective) || !Same(norm, baseline.FinalGradientNorm) ||
            !Same(PairwiseSum.Dot(ray.Residual, ray.Residual), baseline.RaySquaredNorm) ||
            !Same(PairwiseSum.Dot(ray.Residual, observation.Residual), baseline.RayInnerProduct) ||
            !(reference <= s.ReferenceTolerance))
            throw new InvalidOperationException("Baseline endpoint or ray control differs.");

        var rng = new Xoshiro256PlusPlus(children[2]);
        var standard = new double[config.MeasurementCounts.Max()][];
        for (int i = 0; i < standard.Length; i++)
        {
            standard[i] = new double[s.OutputWidth];
            for (int j = 0; j < s.OutputWidth; j++) standard[i][j] = Normal.SampleStandard(rng);
        }
        string w1 = Hash(net.W1.Rows), w2 = Hash(net.W2.Rows);
        var matrices = new List<object>();
        foreach (int m in config.MeasurementCounts)
        {
            var paired = new Network(s with { Measurements = m }, children);
            if (Hash(paired.W1.Rows) != w1 || Hash(paired.W2.Rows) != w2)
                throw new InvalidOperationException("Paired generator differs.");
            double scale = 1.0 / Math.Sqrt(m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < s.OutputWidth; j++)
                    if (!Same(paired.A.Rows[i][j], scale * standard[i][j]))
                        throw new InvalidOperationException("Scaled measurement prefix differs.");
            matrices.Add(new { Measurements = m, StandardPrefixSha256 = Hash(standard[..m]),
                ScaledMatrixSha256 = Hash(paired.A.Rows), EntryBitMismatches = 0 });
        }
        return new { baseline.Repetition, baseline.Seed, ChildSeeds = children, W1Sha256 = w1, W2Sha256 = w2,
            TruthSha256 = Hash([truth]), StartSha256 = Hash([start]), Matrices = matrices,
            BaselineScalarAndVectorBitsAgree = true, MaxReferenceDiscrepancy = reference };
    }

    internal static int Run(string[] args, byte[] bytes)
    {
        var config = JsonSerializer.Deserialize<MeasurementSettings>(bytes)!;
        string baselineDir = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(args[0])!, "..", config.BaselineReceipt));
        Trial[] baseline = File.ReadLines(Path.Combine(baselineDir, "trials.jsonl"))
            .Select(line => JsonSerializer.Deserialize<Trial>(line)!).ToArray();
        int[] seeds = SeedTree.Derive(config.Solver.MasterSeed, config.Solver.Repetitions);
        if (!baseline.Select(row => row.Seed).SequenceEqual(seeds)) throw new InvalidOperationException("Baseline seeds differ.");
        string? selected = args.Length > 4 && args[4] != "all" ? args[4] : null;
        int[] indices = selected is null ? Enumerable.Range(0, seeds.Length).ToArray() :
            selected.Split(',').Select(int.Parse).Select(i => i - 1).ToArray();
        string[]? cells = args.Length > 3 && args[3] != "all" ? args[3].Split(',') : null;
        Directory.CreateDirectory(args[1]);
        using (var controls = new StreamWriter(Path.Combine(args[1], "data-checks.jsonl")))
            foreach (int i in indices) controls.WriteLine(JsonSerializer.Serialize(CheckData(config, baseline[i])));
        Console.WriteLine($"Nested rows and archived m=150 data checks passed for {indices.Length} problems.");
        bool accelerate = args.Length <= 5 || args[5] != "full";
        int exit = 0;
        foreach (int m in config.MeasurementCounts)
        {
            string name = $"m{m}";
            if (cells is not null && !cells.Contains(name)) continue;
            string output = Path.Combine(args[1], name);
            if (m == config.Solver.Measurements)
            {
                Directory.CreateDirectory(output);
                File.WriteAllLines(Path.Combine(output, "trials.jsonl"),
                    indices.Select(i => JsonSerializer.Serialize(baseline[i])));
                File.WriteAllText(Path.Combine(output, "summary.json"), JsonSerializer.Serialize(new {
                    ReusedSourceRevision = "8a77d0006ac8e50ca3539f3e3fd26a0097d57bb2",
                    ReusedFrom = config.BaselineReceipt, Measurements = m, Repetitions = indices.Length,
                    Successes = indices.Count(i => baseline[i].Success),
                    GradientStops = indices.Count(i => baseline[i].Stop == "gradient"),
                    UpdateLimitStops = indices.Count(i => baseline[i].Stop == "update-limit"), GatePassed = true
                }, new JsonSerializerOptions { WriteIndented = true }) + "\n");
                Console.WriteLine($"Reused {indices.Length} archived signed m={m} endpoints.");
                continue;
            }
            Settings settings = config.Solver with { Measurements = m, FastForward = accelerate };
            Console.WriteLine($"Starting m={m}, master={settings.MasterSeed}, fast-forward={accelerate}");
            exit |= Program.Execute(settings, bytes, output, args[2], selected);
        }
        return exit;
    }
}
