using System.Text.Json;
using Maths.Rng;

namespace GenerativePrior;

internal sealed record Cell(int K, int? MasterIndex, bool Paired, int PaperSuccesses);
internal sealed record CurveSettings(int NewMasterSeed, int Repetitions, int HiddenWidth, int OutputWidth,
    int Measurements, double Step, double GradientTolerance, int MaxUpdates, double SuccessTolerance,
    double ReferenceTolerance, Cell[] Cells);

internal static class Curve
{
    internal static int Run(string[] args, byte[] bytes)
    {
        var config = JsonSerializer.Deserialize<CurveSettings>(bytes)!;
        int[] masters = SeedTree.Derive(config.NewMasterSeed, 14);
        string[]? selectedCells = args.Length > 3 && args[3] != "all" ? args[3].Split(',') : null;
        string? trials = args.Length > 4 && args[4] != "all" ? args[4] : null;
        bool accelerate = args.Length <= 5 || args[5] != "full";
        int exit = 0;
        foreach (Cell cell in config.Cells)
        {
            foreach (bool sign in cell.Paired ? new[] { true, false } : new[] { true })
            {
                string name = $"k{cell.K}-{(sign ? "on" : "off")}";
                if (selectedCells is not null && !selectedCells.Contains(name)) continue;
                int master = cell.MasterIndex.HasValue ? masters[cell.MasterIndex.Value] : 181204176;
                var settings = new Settings(cell.K, config.HiddenWidth, config.OutputWidth,
                    config.Measurements, config.Repetitions, master, config.Step, config.GradientTolerance,
                    config.MaxUpdates, config.SuccessTolerance, config.ReferenceTolerance, sign, accelerate, false);
                Console.WriteLine($"Starting {name}, master={master}, fast-forward={accelerate}");
                exit |= Program.Execute(settings, bytes, Path.Combine(args[1], name), args[2], trials);
            }
        }
        return exit;
    }
}
