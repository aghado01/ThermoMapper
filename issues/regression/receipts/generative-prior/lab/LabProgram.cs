using System.Diagnostics;
using System.Security.Cryptography;
using System.Text.Json;

namespace GenerativePrior;

internal static class LabProgram
{
    internal static string Revision = "";
    internal static string SettingsSha256 = "";

    public static int Main(string[] args)
    {
        using var process = Process.GetCurrentProcess();
        process.PriorityClass = ProcessPriorityClass.BelowNormal;
        if (Environment.ProcessorCount != 1 || process.PriorityClass != ProcessPriorityClass.BelowNormal)
            throw new InvalidOperationException("Expose one .NET processor and use below-normal priority.");
        byte[] bytes = File.ReadAllBytes(args[1]);
        Revision = args[3];
        SettingsSha256 = Convert.ToHexStringLower(SHA256.HashData(bytes));
        Directory.CreateDirectory(args[2]);
        return args[0] switch {
            "g3" => WeightConditions.Run(JsonSerializer.Deserialize<WeightSettings>(bytes)!, args[2]),
            "g4" => AudioGeometry.Run(JsonSerializer.Deserialize<AudioSettings>(bytes)!, args[2]),
            "g2" => GeneratorCalibration.Run(JsonSerializer.Deserialize<InversionSettings>(bytes)!, args[2]),
            _ => throw new ArgumentException("Unknown lab experiment.")
        };
    }

    internal static void Write(string path, object value) =>
        File.WriteAllText(path, JsonSerializer.Serialize(value, new JsonSerializerOptions { WriteIndented = true }) + "\n");
}
