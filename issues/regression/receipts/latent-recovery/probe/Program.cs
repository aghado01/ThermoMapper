using System;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace LatentRecoveryProbe;

/// <summary>
/// Entry point. Every command runs at below-normal priority with its own thread cap (the machine is shared) and writes
/// under the results folder given as the second argument.
/// </summary>
public static class Program
{
    public static int Main(string[] args)
    {
        try { Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal; } catch { }
        if (args.Length < 2)
        {
            Console.Error.WriteLine("usage: LatentRecoveryProbe <command> <results dir> [options]");
            return 2;
        }
        string command = args[0], dir = args[1];
        Directory.CreateDirectory(dir);
        var sw = Stopwatch.StartNew();
        switch (command)
        {
            case "b1":
                LogPowerCovarianceSweep.RunPairs(LogPowerCovarianceSweep.Default, dir);
                LogPowerCovarianceSweep.RunCells(LogPowerCovarianceSweep.Default, dir, "complex");
                LogPowerCovarianceSweep.RunCells(LogPowerCovarianceSweep.Default, dir, "real");
                break;
            case "b1-rerun-pairs":
                LogPowerCovarianceSweep.RunPairs(LogPowerCovarianceSweep.Default, dir, args.Skip(2).Select(int.Parse).ToArray(), rerun: true);
                break;
            case "b1-rerun-cells":
                LogPowerCovarianceSweep.RunCells(LogPowerCovarianceSweep.Default, dir, args[2], rerun: true);
                break;
            case "b1-gate":
                File.WriteAllText(Path.Combine(dir, "b1-gate.txt"), SweepGates.B1(dir));
                Console.Write(File.ReadAllText(Path.Combine(dir, "b1-gate.txt")));
                break;
            case "b2":
                WindowScalingSweep.Run(WindowScalingSweep.Default, dir);
                break;
            case "b2-rerun-noise":
                WindowScalingSweep.Run(WindowScalingSweep.Default, dir, rerunNoise: true, noiseOnly: args.Skip(2).Select(int.Parse).ToArray());
                break;
            case "b2-supplement":
                File.WriteAllText(Path.Combine(dir, "b2-supplement.txt"), SweepGates.B2Supplement(dir));
                Console.Write(File.ReadAllText(Path.Combine(dir, "b2-supplement.txt")));
                break;
            case "b2-gate":
                File.WriteAllText(Path.Combine(dir, "b2-gate.txt"), SweepGates.B2(dir));
                Console.Write(File.ReadAllText(Path.Combine(dir, "b2-gate.txt")));
                break;
            default:
                if (!Pilot.Dispatch(command, dir, args.Skip(2).ToArray()))
                {
                    Console.Error.WriteLine($"unknown command '{command}'");
                    return 2;
                }
                break;
        }
        Console.WriteLine($"[{command}] {sw.Elapsed.TotalSeconds:F1} s");
        return 0;
    }
}
