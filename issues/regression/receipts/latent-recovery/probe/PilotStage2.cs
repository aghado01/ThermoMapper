using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace LatentRecoveryProbe;

/// <summary>Stage 2: the harmonic source, the two layouts, octave masses.</summary>
public static partial class Pilot
{
    public static Settings Stage2 => Stage1 with
    {
        Stage = "s2", MasterSeed = 2026111602, Replicates = 20, NoiseVariances = [1e-3], Ratios = [2, 4, 8, 16],
        TrainingBase = 100, TrainingPerK = 100, Layouts = ["functional", "scalar-directional"], MaxBasis = 120,
        DiscrepancyNoise = -1, ModalityKs = [2], ModalityRatios = [4, 16], ModalityNoise = 0, ModalityReplicates = 10,
    };

    private static bool DispatchStage2(string command, string dir, string[] args)
    {
        switch (command)
        {
            case "s2-peek": Peek(new HarmonicSynthesizer(3), Stage2, dir, "s2"); return true;
            case "s2-run": Run(Stage2, k => new HarmonicSynthesizer(k), int.Parse(args[0]), dir, args.Contains("rerun"), args.Contains("fit-only")); return true;
            case "s1-fit-only": Run(Stage1, k => new ResonanceSynthesizer(k), int.Parse(args[0]), dir, false, true); return true;
            case "s2-gate": File.WriteAllText(Path.Combine(dir, "s2-readings.txt"), Stage2Readings(dir)); Console.Write(File.ReadAllText(Path.Combine(dir, "s2-readings.txt"))); return true;
            default: return false;
        }
    }

    private static string Stage2Readings(string dir)
    {
        var sb = new StringBuilder("## Stage 2 readings\n\n### Fits\n\n");
        var rows = new List<PilotGates.Row>();
        foreach (int k in new[] { 2, 3 })
        {
            string fit = Path.Combine(dir, $"s2-k{k}-fit.jsonl");
            if (File.Exists(fit)) foreach (string l in File.ReadAllLines(fit)) sb.AppendLine("- " + l);
            string path = Path.Combine(dir, $"s2-k{k}-replicates.csv");
            if (File.Exists(path)) rows.AddRange(PilotGates.Load(path));
        }
        sb.AppendLine("\n### Cells (nested sampler)\n");
        sb.AppendLine(PilotGates.Table(PilotGates.Cells(rows)));
        sb.AppendLine("### f₀: error and octave masses\n");
        sb.AppendLine("| k | emulator | mode | ratio | n | median |f₀ err| (unit) | f₀ coverage | median mass at θ* | replicates with mass > 0.5 at 2f₀ or f₀/2 | mean mass at 2f₀ | mean mass at f₀/2 |\n|---|---|---|---|---|---|---|---|---|---|---|");
        foreach (var g in rows.GroupBy(r => (r.K, r.Emulator, r.Mode, r.Ratio)).OrderBy(g => g.Key.K).ThenBy(g => g.Key.Emulator).ThenBy(g => g.Key.Mode).ThenBy(g => g.Key.Ratio))
        {
            double cover = g.Count(r => r.Truth[0] >= r.Q05[0] && r.Truth[0] <= r.Q95[0]) / (double)g.Count();
            sb.AppendLine($"| {g.Key.K} | {g.Key.Emulator} | {g.Key.Mode} | {g.Key.Ratio} | {g.Count()} | {PilotGates.Median(g.Select(r => Math.Abs(r.Mean[0] - r.Truth[0]))).ToString("F3", CultureInfo.InvariantCulture)} | " +
                          $"{cover.ToString("F2", CultureInfo.InvariantCulture)} | {PilotGates.Median(g.Select(r => r.MassTruth)).ToString("F2", CultureInfo.InvariantCulture)} | " +
                          $"{g.Count(r => r.MassDouble > 0.5 || r.MassHalf > 0.5)} | {g.Average(r => r.MassDouble).ToString("F2", CultureInfo.InvariantCulture)} | {g.Average(r => r.MassHalf).ToString("F2", CultureInfo.InvariantCulture)} |");
        }
        return sb.ToString();
    }
}
