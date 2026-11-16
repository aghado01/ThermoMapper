using System;
using System.Diagnostics;
using System.IO;
using System.Threading;

namespace SplineSpaceProbe;

/// <summary>Entry point; every command runs at below-normal priority on a capped thread pool (the machine is shared).</summary>
public static class Program
{
    public const int Threads = 6;

    public static int Main(string[] args)
    {
        try { Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal; } catch { }
        ThreadPool.SetMaxThreads(Threads, Threads);
        if (args.Length < 2) { Console.Error.WriteLine("usage: SplineSpaceProbe <p1|p1-summary|p2|p2-rerun|p2-summary|p3|p3-sweep|p3-summary> <dir>"); return 2; }
        string dir = args[1];
        Directory.CreateDirectory(dir);
        var sw = Stopwatch.StartNew();
        switch (args[0])
        {
            case "p1": P1.Run(dir); break;
            case "p1-summary": Write(dir, "p1-summary.md", P1.Summary(dir)); break;
            case "p2": P2.Run(dir, rerun: false); break;
            case "p2-rerun": P2.Run(dir, rerun: true); break;
            case "p2-summary": Write(dir, "p2-summary.md", P2.Summary(dir)); break;
            case "p2-diagnostic": P2.Diagnose(dir); break;
            case "p3": P3.Run(dir); break;
            case "p3-sweep": P3.Sweep(dir); break;
            case "p3-uniform": P3.UniformCubic(dir); break;
            case "p3-summary": Write(dir, "p3-summary.md", P3.Summary(dir)); break;
            default: Console.Error.WriteLine($"unknown command {args[0]}"); return 2;
        }
        Console.WriteLine($"[{args[0]}] {sw.Elapsed.TotalSeconds:F1} s");
        return 0;
    }

    private static void Write(string dir, string name, string text)
    {
        File.WriteAllText(Path.Combine(dir, name), text);
        Console.Write(text);
    }
}
