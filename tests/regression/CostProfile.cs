using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime;
using System.Runtime.InteropServices;
using System.Text;

namespace Maths.Regression.Tests;

/// <summary>
/// The cost profile's meters and rows (the brief, <c>issues/regression/performance-cost-regression-samplers-brief.md</c>).
/// <list type="bullet">
///   <item><b>Method</b> (M1's, as <see cref="BassCost"/> takes it): each part is a batch of calls, calibrated to a few
///   milliseconds; a warm-up batch, then the median over nine batches of a batch's mean per call, the parts' batches
///   round-robin with the whole's so that load falls on all alike. A call over 100 ms is a heavy point: five
///   single-call batches.</item>
///   <item><b>Per batch</b> it also reads the thread's allocated bytes, the GC pause and the collection counts.</item>
///   <item><b>Rows</b> go to <c>COST_PROFILE_OUT</c> as block, point, part, quantity, value, flag. A block already in the
///   file is skipped on a rerun, so a restart loses nothing. <c>COST_PROFILE_SCALE</c> = 0.1 is the trial; a block's
///   own first pass runs at a tenth with its rows discarded (the tier-0 guard).</item>
///   <item><b>Load</b> is read around each block (<see cref="Load"/>): the other test hosts and their priorities, the other
///   processes' cores busy, the free memory. A block is tier Q when no other test host ran, B when only below-normal
///   ones did, L otherwise.</item>
/// </list>
/// </summary>
internal static class CostProfile
{
    internal const int Repeats = 9, HeavyRepeats = 5;
    internal const double BatchMs = 5.0, WholeBatchMs = 20.0, HeavyMs = 100.0;

    internal static bool Enabled => Environment.GetEnvironmentVariable("BASS_QUALIFICATION") == "1";

    /// <summary>Whether <c>COST_PROFILE_BLOCK</c> (semicolon-separated prefixes) selects the block; every block when unset.</summary>
    internal static bool Selected(string block)
    {
        string? filter = Environment.GetEnvironmentVariable("COST_PROFILE_BLOCK");
        if (string.IsNullOrWhiteSpace(filter)) return true;
        return filter.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                     .Any(f => block.StartsWith(f, StringComparison.OrdinalIgnoreCase));
    }

    internal static double Scale => ScaleOverride ?? (double.TryParse(Environment.GetEnvironmentVariable("COST_PROFILE_SCALE"),
        NumberStyles.Float, CultureInfo.InvariantCulture, out double s) && s > 0 && s <= 1 ? s : 1.0);

    internal static double? ScaleOverride { get; set; }

    internal static string OutDirectory
    {
        get
        {
            string dir = Environment.GetEnvironmentVariable("COST_PROFILE_OUT") ?? Path.Combine(Path.GetTempPath(), "cost-profile");
            Directory.CreateDirectory(dir);
            return dir;
        }
    }

    /// <summary>The setting tag a row carries: the <c>DOTNET_*</c> variables in force, or <c>default</c>.</summary>
    internal static string Setting
    {
        get
        {
            var vars = new List<string>();
            foreach (string name in new[] { "DOTNET_TieredPGO", "DOTNET_TieredCompilation", "DOTNET_TC_QuickJitForLoops", "DOTNET_gcConcurrent",
                                             "DOTNET_gcServer", "DOTNET_GCgen0size", "DOTNET_ReadyToRun", "DOTNET_TC_OnStackReplacement" })
            {
                string? v = Environment.GetEnvironmentVariable(name);
                if (!string.IsNullOrEmpty(v)) vars.Add($"{name.Replace("DOTNET_", "", StringComparison.Ordinal)}={v}");
            }
            return vars.Count == 0 ? "default" : string.Join(' ', vars);
        }
    }

    // ---- Rows ----

    private static readonly object Gate = new();
    private static HashSet<string>? _done;
    internal static bool Discard { get; set; }

    private static string RowsPath => Path.Combine(OutDirectory, "cost-profile-rows.csv");

    private static HashSet<string> Done()
    {
        if (_done is not null) return _done;
        _done = [];
        if (File.Exists(RowsPath))
            foreach (string line in File.ReadLines(RowsPath).Skip(1))
            {
                string[] f = line.Split(',');
                if (f.Length >= 8 && f[2] == "(load)" && f[4] == "tier") _done.Add(f[0] + "|" + f[1]);
            }
        return _done;
    }

    internal static bool BlockDone(string block) { lock (Gate) return Done().Contains(Setting + "|" + block); }

    internal static void Row(string block, string point, string part, string quantity, double value, string flag = "")
        => Row(block, point, part, quantity, value.ToString("R", CultureInfo.InvariantCulture), flag);

    internal static void Row(string block, string point, string part, string quantity, string value, string flag = "")
    {
        if (Discard) return;
        string line = string.Join(',', Clean(Setting), Clean(block), Clean(point), Clean(part), Clean(quantity), Clean(value), Clean(flag),
                                  Scale.ToString("R", CultureInfo.InvariantCulture)) + "\n";
        lock (Gate)
        {
            if (_buffer is not null) { _buffer.Append(line); return; }
            Append(line);
        }
    }

    // A block's rows are held until it finishes, so a run that dies midway leaves none of them and the block reruns whole.
    private static StringBuilder? _buffer;

    private static void Append(string text)
    {
        bool fresh = !File.Exists(RowsPath);
        using var w = new StreamWriter(RowsPath, append: true, new UTF8Encoding(false));
        if (fresh) w.Write("setting,block,point,part,quantity,value,flag,scale\n");
        w.Write(text);
    }

    private static string Clean(string s) => s.Replace(',', ';').Replace('\n', ' ');

    // ---- Parts and their reading ----

    /// <summary>A timed part: a batch of <c>count</c> calls, its counted units per call, and the calls per whole.</summary>
    internal sealed class Part(string name, Action<int> run, double units = 0, string unitName = "", double perWhole = double.NaN)
    {
        public string Name { get; } = name;
        public Action<int> Run { get; } = run;
        public double Units { get; } = units;
        public string UnitName { get; } = unitName;
        public double PerWhole { get; } = perWhole;
    }

    internal readonly record struct Reading(string Name, double Median, double Q1, double Q3, double Bytes, double PauseNs, double Gen0, double Gen1,
                                            double Gen2, int Count, bool Heavy, double Units, string UnitName, double PerWhole)
    {
        public double NsPerUnit => Units > 0 ? Median / Units : double.NaN;
    }

    internal static double Median(IReadOnlyList<double> v) => Quantile(v, 0.5);

    internal static double Quantile(IReadOnlyList<double> v, double q)
    {
        if (v.Count == 0) return double.NaN;
        double[] a = v.Order().ToArray();
        double h = q * (a.Length - 1);
        int lo = (int)Math.Floor(h);
        return lo + 1 < a.Length ? a[lo] + (h - lo) * (a[lo + 1] - a[lo]) : a[lo];
    }

    private static double Seconds(long ticks) => ticks / (double)Stopwatch.Frequency;

    /// <summary>Reads the parts round-robin: median ns per call over the batches, with bytes, pause and collections per call.</summary>
    internal static Reading[] Measure(IReadOnlyList<Part> parts, double batchMs = BatchMs, int repeats = Repeats)
    {
        int m = parts.Count;
        double scaled = Math.Max(0.5, batchMs * Scale);
        var counts = new int[m];
        var heavy = new bool[m];
        for (int p = 0; p < m; p++)
        {
            if (Warmed.Add(CurrentBlock + "/" + parts[p].Name)) WarmUp(parts[p]);
            (counts[p], heavy[p]) = Calibrate(parts[p], scaled);
        }
        for (int p = 0; p < m; p++) parts[p].Run(counts[p]);   // the warm-up batch, discarded
        var times = new List<double>[m];
        var bytes = new List<double>[m];
        var pause = new double[m];
        var calls = new long[m];
        var gen = new double[m, 3];
        for (int p = 0; p < m; p++) { times[p] = []; bytes[p] = []; }
        bool anyHeavy = heavy.Any(h => h);
        int rounds = anyHeavy ? HeavyRepeats : repeats;
        for (int r = 0; r < rounds; r++)
            for (int p = 0; p < m; p++)
            {
                int c = counts[p];
                long a0 = GC.GetAllocatedBytesForCurrentThread();
                TimeSpan g0 = GC.GetTotalPauseDuration();
                int c0 = GC.CollectionCount(0), c1 = GC.CollectionCount(1), c2 = GC.CollectionCount(2);
                long t0 = Stopwatch.GetTimestamp();
                parts[p].Run(c);
                long t1 = Stopwatch.GetTimestamp();
                times[p].Add(Seconds(t1 - t0) * 1e9 / c);
                bytes[p].Add((GC.GetAllocatedBytesForCurrentThread() - a0) / (double)c);
                pause[p] += (GC.GetTotalPauseDuration() - g0).TotalMilliseconds * 1e6;
                calls[p] += c;
                gen[p, 0] += GC.CollectionCount(0) - c0;
                gen[p, 1] += GC.CollectionCount(1) - c1;
                gen[p, 2] += GC.CollectionCount(2) - c2;
            }
        var readings = new Reading[m];
        for (int p = 0; p < m; p++)
            readings[p] = new Reading(parts[p].Name, Median(times[p]), Quantile(times[p], 0.25), Quantile(times[p], 0.75), Median(bytes[p]),
                                      pause[p] / calls[p], gen[p, 0] / calls[p], gen[p, 1] / calls[p], gen[p, 2] / calls[p], counts[p], heavy[p],
                                      parts[p].Units, parts[p].UnitName, parts[p].PerWhole);
        return readings;
    }

    /// <summary>The block being run; a part is warmed on its first sight in a block, in this process.</summary>
    internal static string CurrentBlock { get; set; } = "";

    private static readonly HashSet<string> Warmed = [];

    /// <summary>Milliseconds a part runs the first time it is seen before it is calibrated. The runtime promotes a method from
    /// tier 0 only after a 100 ms call-counting delay and then compiles it in the background (twice, with dynamic PGO), so a
    /// short warm-up reads unoptimized code: a vector <c>Exp</c> read about 2 µs per element in the trial.</summary>
    internal const double WarmMs = 400.0;

    private static void WarmUp(Part part)
    {
        var watch = Stopwatch.StartNew();
        int c = 1;
        while (watch.Elapsed.TotalMilliseconds < WarmMs)
        {
            long t0 = Stopwatch.GetTimestamp();
            part.Run(c);
            double ms = Seconds(Stopwatch.GetTimestamp() - t0) * 1e3;
            if (ms / c > HeavyMs) break;
            if (ms < 2.0 && c < 1 << 24) c *= 2;
        }
    }

    private static (int Count, bool Heavy) Calibrate(Part part, double batchMs)
    {
        part.Run(1);
        int c = 1;
        double elapsed;
        while (true)
        {
            long t0 = Stopwatch.GetTimestamp();
            part.Run(c);
            elapsed = Seconds(Stopwatch.GetTimestamp() - t0) * 1e3;
            if (elapsed / c > HeavyMs) return (1, true);
            if (elapsed >= batchMs / 4 || c >= 1 << 24) break;
            double grow = Math.Clamp(batchMs / 4 / Math.Max(elapsed, 1e-4), 2.0, 64.0);
            c = (int)Math.Min(c * grow, 1 << 24);
        }
        return (Math.Max(1, (int)Math.Round(batchMs / (elapsed / c))), false);
    }

    /// <summary>Writes a reading as rows: the median and its quartiles in ns per call, bytes and pause per call, the counted
    /// units, ns per unit and the calls per whole.</summary>
    internal static void Write(string block, string point, Reading r, string flag = "")
    {
        Row(block, point, r.Name, "ns", r.Median, flag);
        Row(block, point, r.Name, "ns_q1", r.Q1, flag);
        Row(block, point, r.Name, "ns_q3", r.Q3, flag);
        Row(block, point, r.Name, "bytes", r.Bytes, flag);
        Row(block, point, r.Name, "pause_ns", r.PauseNs, flag);
        Row(block, point, r.Name, "gen0", r.Gen0, flag);
        Row(block, point, r.Name, "gen2", r.Gen2, flag);
        if (r.Units > 0)
        {
            Row(block, point, r.Name, "units", r.Units, flag);
            Row(block, point, r.Name, "ns_per_unit", r.NsPerUnit, flag);
            if (!string.IsNullOrEmpty(r.UnitName)) Row(block, point, r.Name, "unit:" + r.UnitName, r.Units, flag);
        }
        if (!double.IsNaN(r.PerWhole)) Row(block, point, r.Name, "per_whole", r.PerWhole, flag);
        if (r.Heavy) Row(block, point, r.Name, "heavy", 1, flag);
    }

    /// <summary>Times a point's parts and writes them; a part called <c>(empty)</c> is the timer's overhead, added here.</summary>
    internal static Reading[] Point(string block, string point, IReadOnlyList<Part> parts, string flag = "", double batchMs = BatchMs)
    {
        var all = new List<Part>(parts) { new("(empty)", c => { for (int i = 0; i < c; i++) { } }) };
        Reading[] readings = Measure(all, batchMs);
        foreach (Reading r in readings) Write(block, point, r, flag);
        return readings;
    }

    // ---- Blocks and the load around them ----

    internal sealed record LoadSnapshot(DateTime At, Dictionary<int, double> Cpu, double OwnCpu, int OtherTestHosts, string Priorities, double FreeMemoryGb);

    [StructLayout(LayoutKind.Sequential, CharSet = CharSet.Auto)]
    private struct MemoryStatus
    {
        public uint Length, MemoryLoad;
        public ulong TotalPhys, AvailPhys, TotalPageFile, AvailPageFile, TotalVirtual, AvailVirtual, AvailExtendedVirtual;
    }

    [DllImport("kernel32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool GlobalMemoryStatusEx(ref MemoryStatus status);

    internal static double FreeMemoryGb()
    {
        var status = new MemoryStatus { Length = (uint)Marshal.SizeOf<MemoryStatus>() };
        return GlobalMemoryStatusEx(ref status) ? status.AvailPhys / 1073741824.0 : double.NaN;
    }

    internal static LoadSnapshot Snapshot()
    {
        var cpu = new Dictionary<int, double>();
        int hosts = 0;
        var priorities = new List<string>();
        int self = Environment.ProcessId;
        foreach (Process p in Process.GetProcesses())
        {
            try
            {
                cpu[p.Id] = p.TotalProcessorTime.TotalSeconds;
                if (p.Id != self && p.ProcessName.Contains("testhost", StringComparison.OrdinalIgnoreCase))
                {
                    hosts++;
                    priorities.Add(p.PriorityClass.ToString());
                }
            }
            catch (Exception) { /* protected or exited */ }
            finally { p.Dispose(); }
        }
        double own = Process.GetCurrentProcess().TotalProcessorTime.TotalSeconds;
        return new LoadSnapshot(DateTime.UtcNow, cpu, own, hosts, string.Join('/', priorities), FreeMemoryGb());
    }

    /// <summary>The tier of a block from the test hosts seen at its start and end: Q, B (only below-normal hosts) or L.</summary>
    internal static string Tier(LoadSnapshot start, LoadSnapshot end)
    {
        int hosts = Math.Max(start.OtherTestHosts, end.OtherTestHosts);
        if (hosts == 0) return "Q";
        string seen = start.Priorities + "/" + end.Priorities;
        bool onlyBelow = seen.Split('/', StringSplitOptions.RemoveEmptyEntries).All(p => p is "BelowNormal" or "Idle");
        return onlyBelow ? "B" : "L";
    }

    /// <summary>Runs a block: skipped when its rows exist, else its body once at a tenth with rows discarded, then for the
    /// record, with the load read before and after. Returns whether it ran.</summary>
    internal static bool Run(string block, Action body, bool discardPass = true)
    {
        if (!Enabled || !Selected(block)) return false;
        if (BlockDone(block)) return false;
        CurrentBlock = block;
        double free = FreeMemoryGb();
        LoadSnapshot start = Snapshot();
        if (discardPass)
        {
            Discard = true;
            ScaleOverride = 0.1;
            try { body(); }
            finally { Discard = false; ScaleOverride = null; }
        }
        var watch = Stopwatch.StartNew();
        LoadSnapshot mid = Snapshot();
        lock (Gate) _buffer = new StringBuilder();
        try { body(); }
        catch { lock (Gate) _buffer = null; throw; }
        watch.Stop();
        LoadSnapshot end = Snapshot();
        string tier = Tier(start, end);
        double wall = (end.At - mid.At).TotalSeconds;
        double others = 0.0;
        foreach (KeyValuePair<int, double> kv in end.Cpu)
            if (kv.Key != Environment.ProcessId)
                others += kv.Value - (mid.Cpu.TryGetValue(kv.Key, out double c0) ? c0 : 0.0);
        double busy = wall > 0 ? others / wall : double.NaN;
        double own = wall > 0 ? (end.OwnCpu - mid.OwnCpu) / wall : double.NaN;
        Row(block, "(load)", "", "tier", tier);
        Row(block, "(load)", "", "other_test_hosts_start", start.OtherTestHosts);
        Row(block, "(load)", "", "other_test_hosts_end", end.OtherTestHosts);
        Row(block, "(load)", "", "priorities", end.Priorities.Length == 0 ? "-" : end.Priorities);
        Row(block, "(load)", "", "other_cores_busy", busy);
        Row(block, "(load)", "", "own_cpu_over_wall", own);
        Row(block, "(load)", "", "free_memory_gb_start", free);
        Row(block, "(load)", "", "wall_s", watch.Elapsed.TotalSeconds);
        lock (Gate)
        {
            Append(_buffer!.ToString());
            _buffer = null;
            Done().Add(Setting + "|" + block);
        }
        return true;
    }

    // ---- The machine ----

    internal static void MachineRows()
    {
        const string b = "(machine)";
        Row(b, "", "", "runtime", System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription);
        Row(b, "", "", "os", System.Runtime.InteropServices.RuntimeInformation.OSDescription);
        Row(b, "", "", "logical_processors", Environment.ProcessorCount);
        Row(b, "", "", "server_gc", GCSettings.IsServerGC ? "true" : "false");
        Row(b, "", "", "gc_latency", GCSettings.LatencyMode.ToString());
        Row(b, "", "", "setting", Setting);
        Row(b, "", "", "tiered_pgo_env", Environment.GetEnvironmentVariable("DOTNET_TieredPGO") ?? "unset");
        Row(b, "", "", "free_memory_gb", FreeMemoryGb());
        Row(b, "", "", "vector_hw", $"V128={System.Runtime.Intrinsics.Vector128.IsHardwareAccelerated} V256={System.Runtime.Intrinsics.Vector256.IsHardwareAccelerated} V512={System.Runtime.Intrinsics.Vector512.IsHardwareAccelerated}");
    }
}
