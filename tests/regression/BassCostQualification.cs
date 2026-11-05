using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The cost qualification's grid and its reader (the cost brief, <c>issues/regression/bass/performance-cost-brief.md</c>).
/// Every value the brief fixes before a run is declared here: the surfaces, the base point, the axes, the burn-in rule
/// and M1's timing method (<see cref="BassMixingCostTests"/>), so a change after a timed run shows in the diff.
/// <list type="bullet">
///   <item><b>Surfaces.</b> S10's (x₁, x₂ carry the signal, the rest inert); a rich one, every predictor active through a
///   one-way sinusoid and a two-way one with its neighbour, noise sd 0.05, so the cap sets M; a categorical one, five
///   continuous predictors on the rich surface and two categoricals of L balanced levels with level effects; a direct
///   functional one over p = 5 and an output grid of m points; a basis one of 40 smooth components over m = 200.</item>
///   <item><b>A point</b> is read on chain 0 of its config, started and burned as the run does it: 2,000 sweeps, then
///   blocks of 1,000 while M is below the cap, until three blocks in a row leave M no higher than it has been, to 20,000
///   at most. Counts and times read the same state, since the burn-in is deterministic.</item>
///   <item><b>Its reading</b> is M1's: each part of a sweep at that state, the median over nine batches of a batch's mean
///   per call after a warm-up batch, the parts' batches round-robin with the sweep's; or, in the counts pass, each part's
///   and the sweep's bytes allocated per call, by <see cref="GC.GetTotalAllocatedBytes(bool)"/> deltas beside the
///   thread's own counter.</item>
/// </list>
/// Rows go to <c>BASS_COST_OUT</c> (a folder; the temp folder's <c>bass-cost</c> otherwise) as block, point, quantity,
/// value. <c>BASS_COST_SCALE</c> = 0.1 is the trial: a tenth of every timed batch, the burn-in whole.
/// </summary>
internal static class BassCost
{
    // M1's method, and the brief's base point.
    internal const int Repeats = 9, Pool = 256, SweepBatch = 100, PartBatch = 200, SegmentBatch = 100;
    internal const int BurnIn = 2000, BurnBlock = 1000, BurnMax = 20000, BurnPatience = 3;
    internal const int BaseN = 1000, BaseP = 10, BaseCap = 20, BaseOrder = 2;
    internal const double Noise = 0.05;

    /// <summary>The fraction of every timed batch: 1, or the trial's 0.1; a block's warm-up pass overrides it.</summary>
    internal static double Scale => ScaleOverride ?? (double.TryParse(Environment.GetEnvironmentVariable("BASS_COST_SCALE"), NumberStyles.Float,
                                                      CultureInfo.InvariantCulture, out double s) && s > 0 && s <= 1 ? s : 1.0);

    /// <summary>The scale a block's warm-up pass runs at, while it runs.</summary>
    internal static double? ScaleOverride { get; set; }

    internal static int Scaled(int count) => Math.Max(2, (int)Math.Round(count * Scale));

    internal enum Surface { S10, Rich, Categorical, Functional, Basis }

    /// <summary>A grid point: its surface and sizes. L is the categoricals' level count; Nx and Locations the direct
    /// functional's cases and output grid; Rank a basis fit's components.</summary>
    internal sealed record Point(Surface Surface, int N, int P, int Cap, int Order, int L = 0, int Locations = 0, int Rank = 0)
    {
        public string Key => Surface switch
        {
            Surface.Categorical => $"categorical n={N} p=5+2 L={L} M_max={Cap} K={Order}",
            Surface.Functional => $"functional n_x={N} m={Locations} p={P} M_max={Cap} K={Order}",
            Surface.Basis => $"basis n={N} m={Locations} r={Rank} p={P} M_max={Cap} K={Order}",
            _ => $"{(Surface == Surface.S10 ? "S10" : "rich")} n={N} p={P} M_max={Cap} K={Order}",
        };

        /// <summary>The rows of the design: n, or n_x·m for a direct functional.</summary>
        public int Rows => Surface == Surface.Functional ? N * Locations : N;
    }

    internal static Point Base => new(Surface.Rich, BaseN, BaseP, BaseCap, BaseOrder);

    /// <summary>The grid points <c>BASS_COST_POINTS</c> names (their <see cref="Point.Key"/>s, separated by semicolons, in the
    /// order given), read as one block "points" instead of the axes; none when it is unset. A before-and-after reading of a
    /// change reads the handful of points it prices, not whole axes.</summary>
    internal static Point[] NamedPoints()
    {
        if (Environment.GetEnvironmentVariable("BASS_COST_POINTS") is not { Length: > 0 } list) return [];
        Point[] all = [.. Axes.SelectMany(a => a.Points)];
        return [.. list.Split(';', StringSplitOptions.TrimEntries | StringSplitOptions.RemoveEmptyEntries)
                       .Select(k => all.FirstOrDefault(p => p.Key == k) ?? throw new ArgumentException($"No grid point is named '{k}'."))];
    }

    /// <summary>The axes, each varied with the others at the base, and S10's crossed n × M extended; a point on two lists
    /// is read once.</summary>
    internal static readonly (string Axis, Point[] Points)[] Axes =
    [
        ("n", [.. new[] { 200, 500, 1000, 2000, 5000, 10000, 20000, 50000 }.Select(n => Base with { N = n })]),
        ("M_max at n = 1,000", [.. new[] { 5, 10, 20, 50, 100, 200 }.Select(c => Base with { Cap = c })]),
        ("M_max at n = 10,000", [.. new[] { 5, 10, 20, 50, 100, 200 }.Select(c => Base with { N = 10000, Cap = c })]),
        ("p, rich", [.. new[] { 2, 5, 10, 20, 50, 100, 200 }.Select(p => Base with { P = p })]),
        ("p, S10", [.. new[] { 2, 5, 10, 20, 50, 100, 200 }.Select(p => Base with { Surface = Surface.S10, P = p })]),
        ("order", [.. new[] { 1, 2, 3, 4, 5 }.Select(k => Base with { Order = k })]),
        ("L", [.. new[] { 2, 4, 16, 64, 256, 1024 }.Select(l => new Point(Surface.Categorical, 4096, 7, BaseCap, BaseOrder, L: l))]),
        ("direct functional", [.. from nx in new[] { 100, 500 } from m in new[] { 10, 50, 200, 1000 }
                                  select new Point(Surface.Functional, nx, 5, BaseCap, BaseOrder, Locations: m)]),
        ("basis component", [new Point(Surface.Basis, 500, 5, BaseCap, BaseOrder, Locations: 200, Rank: 1)]),
        ("n × M_max", [.. from n in new[] { 200, 1000, 5000, 20000 } from c in new[] { 5, 20, 50, 100 } select Base with { N = n, Cap = c }]),
    ];

    /// <summary>The largest point of each axis, whose peak working set is read in its own process.</summary>
    internal static readonly (string Axis, Point Point)[] Largest =
    [
        ("n", Base with { N = 50000 }),
        ("M_max", Base with { N = 10000, Cap = 200 }),
        ("p", Base with { P = 200 }),
        ("order", Base with { Order = 5 }),
        ("L", new Point(Surface.Categorical, 4096, 7, BaseCap, BaseOrder, L: 1024)),
        ("functional", new Point(Surface.Functional, 500, 5, BaseCap, BaseOrder, Locations: 1000)),
        ("basis rank", new Point(Surface.Basis, 500, 5, BaseCap, BaseOrder, Locations: 200, Rank: 30)),
    ];

    // ---- The surfaces ----

    /// <summary>The point's config: D01–D06's defaults, the cap and the interaction caps, two chains.</summary>
    internal static BassConfig Config(Point x, int burnIn = BurnIn, int samples = 500, int chains = 2) => new()
    {
        Structure = new BassStructureSpec { MaxBasis = x.Cap, MaxInteraction = x.Order, MaxCategoricalInteraction = x.Order, MaxFunctionalInteraction = x.Order },
        Sampling = new BassSampling { Chains = chains, MasterSeed = 71, BurnIn = burnIn, Samples = samples },
    };

    internal static BassData Data(Point x) => x.Surface switch
    {
        Surface.S10 => S10(x.N, x.P),
        Surface.Rich => Rich(x.N, x.P),
        Surface.Categorical => Categorical(x.N, x.L),
        Surface.Functional => Functional(x.N, x.Locations, x.P),
        _ => BasisComponent(x),
    };

    // S10's surface, BassQualificationCostTests' generator.
    internal static BassData S10(int n, int p)
    {
        var rng = new Xoshiro256PlusPlus(seed: 73);
        double[][] x = Columns(p, n);
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int v = 0; v < p; v++) x[v][i] = rng.NextDouble();
            y[i] = 10 * Math.Sin(Math.PI * x[0][i] * x[1][i]) + 10 * x[1][i] * x[1][i] + 3 * Math.Sin(6 * Math.PI * x[0][i])
                   + Noise * BassQualification.Gaussian(rng);
        }
        return new BassData(Continuous(x), new ScalarResponse(y));
    }

    // The rich surface: Σ_v w_v [a_v sin(2π f_v x_v + φ_v) + 0.8 sin(2π x_v x_{v+1})], every predictor in a one-way and a
    // two-way term, the inputs uniform. w_v = 1 for the first ten predictors and 10/(v + 1) beyond, so that past p = 10 a
    // term's share of the variance does not fall as 1/p and the cap, not the Occam cost, sets M (revised after the first
    // counts pass, before any timed run: with w_v = 1 throughout, p = 100 held 5 terms of 20 and p = 200 none).
    internal static BassData Rich(int n, int p)
    {
        var rng = new Xoshiro256PlusPlus(seed: 79);
        double[][] x = Columns(p, n);
        for (int v = 0; v < p; v++)
            for (int i = 0; i < n; i++) x[v][i] = rng.NextDouble();
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = RichMean(x, i, p) + Noise * BassQualification.Gaussian(rng);
        return new BassData(Continuous(x), new ScalarResponse(y));
    }

    // The rich mean at case i over the first p columns of x.
    private static double RichMean(double[][] x, int i, int p)
    {
        double f = 0.0;
        for (int v = 0; v < p; v++)
        {
            double a = 1.0 + 0.5 * Math.Cos(v), frequency = 1.0 + v % 3, phase = 0.7 * v, weight = v < 10 ? 1.0 : 10.0 / (v + 1);
            double term = a * Math.Sin(2 * Math.PI * frequency * x[v][i] + phase);
            if (p > 1) term += 0.8 * Math.Sin(2 * Math.PI * x[v][i] * x[(v + 1) % p][i]);
            f += weight * term;
        }
        return f;
    }

    /// <summary>The categorical surface: five continuous predictors on the rich surface and two categoricals of
    /// <paramref name="levels"/> levels, balanced (case i takes level i mod L in the first and a shuffled one in the
    /// second, so every level has ⌊n/L⌋ or more cases), each level an effect, and the first's effect sloped in x₁.</summary>
    internal static BassData Categorical(int n, int levels)
    {
        var rng = new Xoshiro256PlusPlus(seed: 83);
        double[][] x = Columns(5, n);
        for (int v = 0; v < 5; v++)
            for (int i = 0; i < n; i++) x[v][i] = rng.NextDouble();
        int[] z1 = [.. Enumerable.Range(0, n).Select(i => i % levels)];
        int[] z2 = [.. z1];
        for (int i = n - 1; i > 0; i--)
        {
            int j = (int)(rng.NextUInt64() % (ulong)(i + 1));
            (z2[i], z2[j]) = (z2[j], z2[i]);
        }
        double[] e1 = [.. Enumerable.Range(0, levels).Select(_ => BassQualification.Gaussian(rng))];
        double[] e2 = [.. Enumerable.Range(0, levels).Select(_ => BassQualification.Gaussian(rng))];
        var y = new double[n];
        for (int i = 0; i < n; i++)
            y[i] = RichMean(x, i, 5) + e1[z1[i]] * (1 + x[0][i]) + e2[z2[i]] + Noise * BassQualification.Gaussian(rng);
        string[] names = [.. Enumerable.Range(0, levels).Select(l => $"l{l}")];
        return new BassData([.. Continuous(x), new CategoricalPredictor("z1", z1, names), new CategoricalPredictor("z2", z2, names)], new ScalarResponse(y));
    }

    /// <summary>The direct functional surface: y(x, t) = f(x)(1 + t) + sin(2π t x₁) over the rich f on p inputs, t on m
    /// equally spaced locations in [0, 1], noise sd 0.05.</summary>
    internal static BassData Functional(int cases, int locations, int p)
    {
        var rng = new Xoshiro256PlusPlus(seed: 89);
        double[][] x = Columns(p, cases);
        for (int v = 0; v < p; v++)
            for (int i = 0; i < cases; i++) x[v][i] = rng.NextDouble();
        double[] t = [.. Enumerable.Range(0, locations).Select(r => locations == 1 ? 0.0 : r / (locations - 1.0))];
        var y = new double[cases, locations];
        for (int i = 0; i < cases; i++)
        {
            double f = RichMean(x, i, p);
            for (int r = 0; r < locations; r++) y[i, r] = f * (1 + t[r]) + Math.Sin(2 * Math.PI * t[r] * x[0][i]) + Noise * BassQualification.Gaussian(rng);
        }
        return new BassData(Continuous(x), new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    /// <summary>The basis surface's responses, n × m: Σ_k g_k(x)·sin(kπt)/k over 40 components, g_k a rich-surface
    /// sinusoid in the inputs, noise sd 0.05; and its predictors.</summary>
    internal static (double[,] Y, BassPredictor[] Predictors) BasisSurface(int n, int locations, int p)
    {
        var rng = new Xoshiro256PlusPlus(seed: 97);
        double[][] x = Columns(p, n);
        for (int v = 0; v < p; v++)
            for (int i = 0; i < n; i++) x[v][i] = rng.NextDouble();
        var y = new double[n, locations];
        for (int i = 0; i < n; i++)
            for (int k = 1; k <= 40; k++)
            {
                double g = Math.Sin(2 * Math.PI * (x[k % p][i] + 0.3 * k * x[(k + 1) % p][i]) + k);
                for (int r = 0; r < locations; r++) y[i, r] += g * Math.Sin(k * Math.PI * r / (locations - 1.0)) / k;
            }
        for (int i = 0; i < n; i++)
            for (int r = 0; r < locations; r++) y[i, r] += Noise * BassQualification.Gaussian(rng);
        return (y, Continuous(x));
    }

    /// <summary>The basis point's decomposition at its rank.</summary>
    internal static BassDecomposition Decomposition(Point x)
    {
        (double[,] y, _) = BasisSurface(x.N, x.Locations, x.P);
        return BassDecomposition.FromPca(y, new BassPcaSpec { Rank = x.Rank });
    }

    // A basis fit's first component, the scalar fit its sweeps run on.
    private static BassData BasisComponent(Point x)
    {
        (double[,] y, BassPredictor[] predictors) = BasisSurface(x.N, x.Locations, x.P);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 1, Residual = TruncationResidual.InSample });
        return new BassData(predictors, new ScalarResponse(d.CoefficientResponse(0)));
    }

    private static double[][] Columns(int p, int n)
    {
        var x = new double[p][];
        for (int v = 0; v < p; v++) x[v] = new double[n];
        return x;
    }

    private static BassPredictor[] Continuous(double[][] x) => [.. x.Select((c, v) => (BassPredictor)new ContinuousPredictor($"x{v + 1}", c))];

    // ---- Burn-in ----

    /// <summary>Chain 0 of the point, started as the run starts it and burned by the declared rule: the sweeps taken.</summary>
    internal static (BassModel Model, BassChain Chain, int Sweeps) Burned(Point x)
    {
        BassModel model = BassModel.Create(Config(x), Data(x));
        BassChain chain = BassMixing.Chain(model, 0);
        int sweeps = 0, flat = 0, best = 0;
        for (; sweeps < BurnIn; sweeps++) chain.Step();
        while (chain.Current.Count < x.Cap && sweeps < BurnMax && flat < BurnPatience)
        {
            best = Math.Max(best, chain.Current.Count);
            for (int s = 0; s < BurnBlock; s++) chain.Step();
            sweeps += BurnBlock;
            flat = chain.Current.Count > best ? 0 : flat + 1;
        }
        chain.Structure.ResetTallies();
        return (model, chain, sweeps);
    }

    /// <summary>A term's order: its factors over every part.</summary>
    internal static int OrderOf(BassTerm t) => t.Continuous.Count + t.Categorical.Count + t.Functional.Count;

    // ---- The clock ----

    /// <summary>Whether a time is read on the thread's own cycle counter (<c>BASS_COST_CLOCK</c> = cycles) rather than the wall clock. A
    /// thread starved by other processes' load keeps its cycles, so the counter reads the work a part does when the machine is not quiet; it
    /// is an indicative reading for a comparison of builds in one window, in µs at the rate the counter ran at over a busy loop (the best of
    /// several, the least starved), and a timing for the record still waits for a quiet machine and the wall clock.</summary>
    internal static bool Cycles { get; } = Environment.GetEnvironmentVariable("BASS_COST_CLOCK") == "cycles";

    [DllImport("kernel32.dll")]
    private static extern IntPtr GetCurrentThread();

    [DllImport("kernel32.dll")]
    private static extern bool QueryThreadCycleTime(IntPtr thread, out ulong cycles);

    /// <summary>The clock's reading: a wall-clock tick, or the thread's cycles.</summary>
    internal static long Tick()
    {
        if (!Cycles) return Stopwatch.GetTimestamp();
        QueryThreadCycleTime(GetCurrentThread(), out ulong cycles);
        return (long)cycles;
    }

    /// <summary>The clock's ticks per second: the wall clock's, or the cycle counter's calibrated rate.</summary>
    internal static double TicksPerSecond { get; } = Cycles ? CalibratedCycleRate() : Stopwatch.Frequency;

    // The best rate over eight 25 ms busy loops: a starved loop reads low, an unstarved one reads the core's clock.
    private static double CalibratedCycleRate()
    {
        double best = 0.0, sink = 0.0;
        for (int trial = 0; trial < 8; trial++)
        {
            QueryThreadCycleTime(GetCurrentThread(), out ulong c0);
            long w0 = Stopwatch.GetTimestamp();
            while (Stopwatch.GetTimestamp() - w0 < Stopwatch.Frequency / 40) sink += 1.0;
            QueryThreadCycleTime(GetCurrentThread(), out ulong c1);
            best = Math.Max(best, (c1 - c0) / ((Stopwatch.GetTimestamp() - w0) / (double)Stopwatch.Frequency));
        }
        return sink >= 0 ? best : 0.0;
    }

    // ---- Meters ----

    /// <summary>What a batch reads: µs per call, or bytes allocated per call (the process's precise counter, and the
    /// thread's own beside it).</summary>
    internal enum Meter { Time, Bytes }

    /// <summary>The process's allocated bytes and the thread's, read at once.</summary>
    internal static (long Process, long Thread) Allocated() => (GC.GetTotalAllocatedBytes(precise: true), GC.GetAllocatedBytesForCurrentThread());

    /// <summary>A batch of <paramref name="batch"/> calls, read whole: µs, or bytes, per call.</summary>
    internal static Func<int, double> Batch(Meter meter, Action call, int batch = PartBatch)
    {
        int size = Scaled(batch);
        return meter == Meter.Time
            ? [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
            {
                long t0 = Tick();
                for (int k = 0; k < size; k++) call();
                return (Tick() - t0) * 1e6 / TicksPerSecond / size;
            }
            : (int _) =>
            {
                long b0 = GC.GetTotalAllocatedBytes(precise: true);
                for (int k = 0; k < size; k++) call();
                return (double)(GC.GetTotalAllocatedBytes(precise: true) - b0) / size;
            };
    }

    /// <summary>A batch of calls that each read their own segment (ticks, or bytes), so their setup is left out, or skip
    /// themselves with −1: the mean per counted call, NaN when none was.</summary>
    internal static Func<int, double> EachBatch(Meter meter, Func<Func<long>, long> call, int batch = SegmentBatch)
    {
        int size = Scaled(batch);
        Func<long> read = meter == Meter.Time ? Tick : () => GC.GetTotalAllocatedBytes(precise: true);
        double unit = meter == Meter.Time ? 1e6 / TicksPerSecond : 1.0;
        return [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
        {
            long total = 0, counted = 0;
            for (int k = 0; k < size; k++)
                if (call(read) is var t && t >= 0)
                {
                    total += t;
                    counted++;
                }
            return counted > 0 ? total * unit / counted : double.NaN;
        };
    }

    /// <summary>Each timer's batch round-robin with the others' (a warm-up round, numbered −1, then nine), so load falls on
    /// every part alike: each timer's nine batch values.</summary>
    internal static Dictionary<string, double[]> Interleave(IReadOnlyList<(string Name, Func<int, double> Batch)> timers, int repeats = Repeats)
    {
        var values = timers.ToDictionary(t => t.Name, _ => new double[repeats]);
        for (int r = -1; r < repeats; r++)
            foreach ((string name, Func<int, double> batch) in timers)
            {
                double v = batch(r);
                if (r >= 0) values[name][r] = v;
            }
        return values;
    }

    internal static double Median(double[] values) => BassMixing.Quantile(values.Where(double.IsFinite), 0.5);

    internal static double Iqr(double[] values)
    {
        double[] finite = [.. values.Where(double.IsFinite)];
        return BassMixing.Quantile(finite, 0.75) - BassMixing.Quantile(finite, 0.25);
    }

    // ---- Load and the machine ----

    /// <summary>The other test hosts on the machine: their ids and names.</summary>
    internal static string[] OtherTestHosts()
    {
        int self = Environment.ProcessId;
        return [.. Process.GetProcesses().Where(p => p.Id != self && p.ProcessName.StartsWith("testhost", StringComparison.OrdinalIgnoreCase))
                                         .Select(p => $"{p.ProcessName}:{p.Id}")];
    }

    /// <summary>The CPU seconds of every other process this one may read, by id.</summary>
    internal static Dictionary<int, double> OtherCpu()
    {
        int self = Environment.ProcessId;
        var cpu = new Dictionary<int, double>();
        foreach (Process p in Process.GetProcesses())
            using (p)
            {
                if (p.Id == self) continue;
                try { cpu[p.Id] = p.TotalProcessorTime.TotalSeconds; }
                catch (Exception e) when (e is InvalidOperationException or System.ComponentModel.Win32Exception or NotSupportedException) { }
            }
        return cpu;
    }

    /// <summary>A block's load: the other test hosts at its start and end, this process's CPU over its wall time, and the
    /// cores the other processes kept busy over it (their CPU seconds over the wall, for those alive at both ends).</summary>
    internal sealed class Load
    {
        private readonly Stopwatch _clock = Stopwatch.StartNew();
        private readonly TimeSpan _cpu0 = Process.GetCurrentProcess().TotalProcessorTime;
        private readonly Dictionary<int, double> _others0 = OtherCpu();
        public string[] AtStart { get; } = OtherTestHosts();

        public (string[] AtEnd, double CpuOverWall, double Seconds, double OtherCores) End()
        {
            double wall = _clock.Elapsed.TotalSeconds;
            double cpu = (Process.GetCurrentProcess().TotalProcessorTime - _cpu0).TotalSeconds;
            Dictionary<int, double> others = OtherCpu();
            double busy = others.Where(kv => _others0.ContainsKey(kv.Key)).Sum(kv => kv.Value - _others0[kv.Key]);
            return (OtherTestHosts(), cpu / wall, wall, busy / wall);
        }
    }

    /// <summary>The machine as the brief records it, from inside the process; the CPU model and the power plan are read
    /// beside the run.</summary>
    internal static IEnumerable<(string Quantity, string Value)> Machine()
    {
#if DEBUG
        yield return ("build", "Debug");
#else
        yield return ("build", "Release");
#endif
        yield return ("commit", Commit());
        yield return ("runtime", RuntimeInformation.FrameworkDescription);
        yield return ("os", RuntimeInformation.OSDescription);
        yield return ("processor identifier", Environment.GetEnvironmentVariable("PROCESSOR_IDENTIFIER") ?? "");
        yield return ("logical processors", Environment.ProcessorCount.ToString(CultureInfo.InvariantCulture));
        yield return ("server GC", GCSettings.IsServerGC.ToString());
        yield return ("GC latency mode", GCSettings.LatencyMode.ToString());
        yield return ("concurrent GC", (AppContext.GetData("System.GC.Concurrent") ?? "default").ToString()!);
        yield return ("tiered compilation", (AppContext.GetData("System.Runtime.TieredCompilation") ?? "default").ToString()!);
        yield return ("tiered PGO", (AppContext.GetData("System.Runtime.TieredPGO") ?? "default").ToString()!);
        yield return ("DOTNET_TieredCompilation", Environment.GetEnvironmentVariable("DOTNET_TieredCompilation") ?? "unset");
        yield return ("scale", Scale.ToString(CultureInfo.InvariantCulture));
        yield return ("clock", Cycles ? $"thread cycles, calibrated at {TicksPerSecond / 1e9:F2} GHz (indicative)" : "wall");
    }

    // The checkout's commit, read from .git above the test's folder.
    private static string Commit()
    {
        for (string? dir = AppContext.BaseDirectory; dir is not null; dir = Path.GetDirectoryName(dir))
        {
            string git = Path.Combine(dir, ".git");
            if (!Directory.Exists(git)) continue;
            string head = File.ReadAllText(Path.Combine(git, "HEAD")).Trim();
            if (!head.StartsWith("ref: ", StringComparison.Ordinal)) return head;
            string refPath = Path.Combine(git, head[5..].Replace('/', Path.DirectorySeparatorChar));
            return File.Exists(refPath) ? File.ReadAllText(refPath).Trim() : head;
        }
        return "unknown";
    }

    // ---- Rows ----

    /// <summary>A block's rows, appended to its file under the output folder as block, point, quantity, value.</summary>
    internal sealed class Rows(string pass, string block, bool discard = false)
    {
        private readonly List<string> _lines = [];
        public string Block => block;

        public void Add(string point, string quantity, double value) => Add(point, quantity, value.ToString("R", CultureInfo.InvariantCulture));

        public void Add(string point, string quantity, string value) => _lines.Add($"{Csv(block)},{Csv(point)},{Csv(quantity)},{Csv(value)}");

        public string Flush()
        {
            if (discard)
            {
                _lines.Clear();
                return "(discarded)";
            }
            string dir = Environment.GetEnvironmentVariable("BASS_COST_OUT") is { Length: > 0 } d ? d : Path.Combine(Path.GetTempPath(), "bass-cost");
            Directory.CreateDirectory(dir);
            string path = Path.Combine(dir, $"{pass}-{Slug(block)}.csv");
            bool fresh = !File.Exists(path);
            using var writer = new StreamWriter(path, append: true);
            if (fresh) writer.WriteLine("block,point,quantity,value");
            foreach (string line in _lines) writer.WriteLine(line);
            _lines.Clear();
            return path;
        }

        private static string Csv(string s) => s.Contains(',') || s.Contains('"') ? $"\"{s.Replace("\"", "\"\"")}\"" : s;

        private static string Slug(string s) => new([.. s.Select(c => char.IsLetterOrDigit(c) ? char.ToLowerInvariant(c) : '-')]);
    }
}
