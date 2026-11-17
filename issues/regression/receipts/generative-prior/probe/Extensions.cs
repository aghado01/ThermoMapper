using System.Text.Json;

namespace GenerativePrior;

internal sealed record ExtendedCell(string Name, int K, int MasterSeed, int Depth, bool Magnitude, bool Paired, bool Bridge = false);
internal sealed record ExtendedSettings(Settings Solver, ExtendedCell[] ExtendedCells);

internal static class Extensions
{
    internal static int Run(string[] args, byte[] bytes)
    {
        var config = JsonSerializer.Deserialize<ExtendedSettings>(bytes)!;
        string[]? selected = args.Length > 3 && args[3] != "all" ? args[3].Split(',') : null;
        string? trials = args.Length > 4 && args[4] != "all" ? args[4] : null;
        bool accelerate = args.Length <= 5 || args[5] != "full";
        int exit = 0;
        foreach (ExtendedCell cell in config.ExtendedCells)
        {
            if (selected is null && cell.Bridge) continue;
            if (selected is not null && !selected.Contains(cell.Name)) continue;
            foreach (bool sign in cell.Paired ? new[] { true, false } : new[] { true })
            {
                var settings = config.Solver with {
                    CodeDimension = cell.K, MasterSeed = cell.MasterSeed,
                    IntermediateWidth = cell.Depth == 3 ? 400 : null, Magnitude = cell.Magnitude,
                    SignCheck = sign, FastForward = accelerate, Step = Math.Pow(2, cell.Depth) / (cell.Depth * cell.Depth)
                };
                string name = $"{cell.Name}-{(sign ? "on" : "off")}";
                Console.WriteLine($"Starting {name}, predicted negative radius={settings.NegativeRadius:R}");
                exit |= Program.Execute(settings, bytes, Path.Combine(args[1], name), args[2], trials);
            }
        }
        return exit;
    }
}
