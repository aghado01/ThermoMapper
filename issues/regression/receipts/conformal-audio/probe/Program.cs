using System.Diagnostics;
using System.Security.Cryptography;
using System.Text.Json;
using Maths.Rng;
using Maths.Regression;
using Maths.Regression.Conformal;

namespace ConformalAudio;

internal static class App
{
    internal static JsonElement Settings;
    internal static string Output = "", Revision = "", SettingsHash = "";
    internal static readonly JsonSerializerOptions Json = new() { WriteIndented = true };
    internal static int I(JsonElement x, string key) => x.GetProperty(key).GetInt32();
    internal static double D(JsonElement x, string key) => x.GetProperty(key).GetDouble();
    internal static int[] Ints(JsonElement x, string key) => x.GetProperty(key).EnumerateArray().Select(v => v.GetInt32()).ToArray();
    internal static double[] Doubles(JsonElement x, string key) => x.GetProperty(key).EnumerateArray().Select(v => v.GetDouble()).ToArray();
    internal static void Require(bool condition, string message) { if (!condition) throw new ArithmeticException(message); }
    internal static string Hash(object value) => Convert.ToHexStringLower(SHA256.HashData(JsonSerializer.SerializeToUtf8Bytes(value)));
    internal static void Write(string name, object value) => File.WriteAllText(Path.Combine(Output, name), JsonSerializer.Serialize(value, Json) + "\n");
    internal static void Line(StreamWriter writer, object value) => writer.WriteLine(JsonSerializer.Serialize(value));
    internal static object Interval(int count, int total)
    {
        double p = count / (double)total, z = 1.959963984540054, den = 1 + z * z / total;
        double center = (p + z * z / (2 * total)) / den;
        double half = z * Math.Sqrt(p * (1 - p) / total + z * z / (4.0 * total * total)) / den;
        return new { Count = count, Total = total, Rate = p, Lower = center - half, Upper = center + half };
    }
    internal static object MeanInterval(double[] values)
    {
        double mean = values.Average();
        double se = values.Length > 1 ? Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1) / values.Length) : 0;
        return new { Mean = mean, StandardError = se, Lower = mean - 1.959963984540054 * se, Upper = mean + 1.959963984540054 * se, Splits = values.Length };
    }
    internal static CandidateLogDensity Seam(double[] conformity) => new(ToMatrix(conformity), conformity.Select(_ => new[] { 0.0 }).ToArray(), []) { Forms = ConformalForms.Split };
    private static double[,] ToMatrix(double[] x) { var m = new double[1, x.Length]; for (int i = 0; i < x.Length; i++) m[0, i] = x[i]; return m; }
    public static int Main(string[] args)
    {
        Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
        Require(Environment.ProcessorCount == 1, "Expose one processor.");
        byte[] bytes = File.ReadAllBytes(args[1]);
        Settings = JsonDocument.Parse(bytes).RootElement.Clone();
        SettingsHash = Convert.ToHexStringLower(SHA256.HashData(bytes)); Output = args[2]; Revision = args[3];
        Directory.CreateDirectory(Output);
        int index = Array.IndexOf(new[] { "controls", "audio", "transport", "risk", "alarm" }, args[0]);
        Require(index >= 0, "Unknown block.");
        int seed = SeedTree.Derive(I(Settings, "MasterSeed"), 5)[index];
        switch (index) { case 0: Controls.Run(); break; case 1: Audio.Run(seed, args.Length > 4 && args[4] == "first-split"); break;
            case 2: Transport.Run(); break; case 3: Risk.Run(seed); break; case 4: Alarms.Run(seed); break; }
        Write("execution.json", new { Revision, SettingsHash, Block = args[0], Seed = seed, Workers = Environment.ProcessorCount,
            Priority = Process.GetCurrentProcess().PriorityClass.ToString() });
        return 0;
    }
}
