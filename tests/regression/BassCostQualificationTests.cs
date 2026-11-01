using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Performance qualification's cost half (the cost brief): what a sweep, a draw and an output cost, by case count,
/// predictors, active terms, interaction order, category cardinality, functional size and basis rank, around the base
/// point (<see cref="BassCost"/>). It reads, it does not judge, and it changes nothing in production.
/// <list type="bullet">
///   <item><b>Counts</b> (<see cref="Counts_EveryPoint"/>, <see cref="Counts_RetainedMemory"/>,
///   <see cref="Counts_PeakWorkingSet"/>): each point's realized M, frequencies and analytic counts, its bytes allocated
///   per sweep and per part, retained bytes per fit and per draw, and the axes' peak working sets, each in its own
///   process. Load does not move them, so they run at once.</item>
///   <item><b>Times</b> (<see cref="Timed_Blocks"/>): each block, named by <c>BASS_COST_BLOCK</c>, run on a quiet
///   machine, its load recorded: the grid's axes, a repeat of the base point, the outputs, the other backends,
///   calibration's N15 and cut SMC, and the parallel fit.</item>
/// </list>
/// Run it alone, in Release, behind <c>BASS_QUALIFICATION</c> = 1; rows go to <c>BASS_COST_OUT</c>.
/// </summary>
public sealed class BassCostQualificationTests
{
    private readonly ITestOutputHelper _out;
    public BassCostQualificationTests(ITestOutputHelper output) => _out = output;

    /// <summary>The grid's blocks: each axis's points not on an earlier axis.</summary>
    internal static (string Block, BassCost.Point[] Points)[] GridBlocks()
    {
        var seen = new HashSet<string>();
        return [.. BassCost.Axes.Select(a => (a.Axis, a.Points.Where(x => seen.Add(x.Key)).ToArray()))];
    }

    internal static readonly string[] OtherBlocks = ["repeat", "outputs", "backends", "calibration", "parallel"];

    // ---- Counts ----

    [Fact]
    public void Counts_EveryPoint()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var rows = new BassCost.Rows("counts", "grid");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        Checks(rows);
        ExportForR(rows);
        rows.Flush();
        BassCostSweep.Read(BassCost.Base with { N = 200 }, BassCost.Meter.Bytes);   // warm-up: JIT and caches
        foreach ((string block, BassCost.Point[] points) in GridBlocks())
        {
            foreach (BassCost.Point x in points)
            {
                BassCostSweep.Reading r = BassCostSweep.Read(x, BassCost.Meter.Bytes);
                Write(rows, r, BassCost.Meter.Bytes, block);
                _out.WriteLine($"{x.Key}: burn {r.BurnSweeps}, M {r.M}/{r.Cap}, {r.Of("sweep"):F0} B/sweep ({r.ThreadBytesPerSweep:F0} on the thread)");
            }
            rows.Flush();
        }
        string path = rows.Flush();
        _out.WriteLine($"rows: {path}");
    }

    // The checks before the runs that are counts: the allocation counter, the surfaces.
    private void Checks(BassCost.Rows rows)
    {
        // A known allocation reads its size: a byte[1000] is 16 bytes of header, 8 of length and 1,000 of data.
        static byte[] Allocate() => new byte[1000];
        Allocate();
        long b0 = GC.GetTotalAllocatedBytes(precise: true);
        byte[] held = Allocate();
        long known = GC.GetTotalAllocatedBytes(precise: true) - b0;
        GC.KeepAlive(held);
        rows.Add("check", "allocation counter: byte[1000], bytes", known);
        Assert.InRange(known, 1024, 1024 + 64);

        // N12's MarginalAt allocates nothing, at the base point's burned design.
        (BassModel model, BassChain chain, _) = BassCost.Burned(BassCost.Base with { N = 200 });
        BassDesign design = chain.Structure.DesignAt(chain.Current);
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                        design.FittedSquares, design.Observations, design.ResponseMean, GPriorCentre.Mean);
        double sink = block.MarginalAt(1.0).LogMarginal(null);
        b0 = GC.GetTotalAllocatedBytes(precise: true);
        for (int k = 0; k < 1000; k++) sink += block.MarginalAt(1.0 + k * 1e-3).LogMarginal(null);
        long marginal = GC.GetTotalAllocatedBytes(precise: true) - b0;
        rows.Add("check", "allocation counter: MarginalAt × 1,000, bytes", marginal);
        Assert.True(double.IsFinite(sink));
        Assert.Equal(0, marginal);

        // L's designs give every level at least 4 cases.
        foreach (BassCost.Point x in BassCost.Axes.Single(a => a.Axis == "L").Points)
        {
            BassData data = BassCost.Data(x);
            int least = int.MaxValue;
            foreach (CategoricalPredictor z in data.Predictors.OfType<CategoricalPredictor>())
                least = Math.Min(least, z.Codes.ToArray().GroupBy(c => c).Min(g => g.Count()));
            rows.Add(x.Key, "check: fewest cases at a level", least);
            Assert.True(least >= 4, $"{x.Key}: a level has {least} cases");
        }
        _ = model;
    }

    // R's context points' data, for r/scripts/bass_cost_context.R: Friedman at n = 500, p = 10; the rich surface at
    // n = 5,000, p = 10; the direct functional at n_x = 100, m = 50.
    private static void ExportForR(BassCost.Rows rows)
    {
        string dir = Path.Combine(OutputFolder(), "r");
        Directory.CreateDirectory(dir);
        Scalar(Path.Combine(dir, "friedman.csv"), BassQualificationRecoveryTests.Friedman().Data);
        Scalar(Path.Combine(dir, "rich.csv"), BassCost.Rich(5000, 10));
        BassData f = BassCost.Functional(100, 50, 5);
        Scalar(Path.Combine(dir, "functional-x.csv"), new BassData(f.Predictors, new ScalarResponse(new double[100])));
        var functional = (FunctionalResponse)f.Response;
        using (var w = new StreamWriter(Path.Combine(dir, "functional-y.csv")))
        {
            w.WriteLine(string.Join(",", Enumerable.Range(0, 50).Select(r => $"t{r}")));
            for (int i = 0; i < 100; i++) w.WriteLine(string.Join(",", Enumerable.Range(0, 50).Select(r => functional[i, r].ToString("R", CultureInfo.InvariantCulture))));
        }
        File.WriteAllLines(Path.Combine(dir, "functional-t.csv"),
            ["t", .. functional.OutputVariables[0].Values.ToArray().Select(v => v.ToString("R", CultureInfo.InvariantCulture))]);
        rows.Add("R", "data folder", dir);

        static void Scalar(string path, BassData data)
        {
            ContinuousPredictor[] x = [.. data.Predictors.Cast<ContinuousPredictor>()];
            double[][] columns = [.. x.Select(c => c.Values.ToArray())];
            double[] y = data.Response.Values.ToArray();
            using var w = new StreamWriter(path);
            w.WriteLine(string.Join(",", x.Select(c => c.Name).Append("y")));
            for (int i = 0; i < y.Length; i++)
                w.WriteLine(string.Join(",", columns.Select(c => c[i].ToString("R", CultureInfo.InvariantCulture)).Append(y[i].ToString("R", CultureInfo.InvariantCulture))));
        }
    }

    [Fact]
    public void Counts_RetainedMemory()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var rows = new BassCost.Rows("counts", "retained");
        int[] lengths = [BassCost.Scaled(250), BassCost.Scaled(1000)];
        BassCost.Point[] points =
        [
            BassCost.Base, BassCost.Base with { N = 10000 }, BassCost.Base with { N = 50000 },
            BassCost.Base with { Cap = 5 }, BassCost.Base with { Cap = 100 }, BassCost.Base with { Cap = 200 },
            new(BassCost.Surface.Functional, 100, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: 10),
            new(BassCost.Surface.Functional, 100, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: 200),
            new(BassCost.Surface.Functional, 100, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: 1000),
        ];
        var run = new BassRunOptions { MaxDegreeOfParallelism = 2 };
        foreach (BassCost.Point x in points)
        {
            BassData data = BassCost.Data(x);
            var retained = new double[2];
            var draws = new int[2];
            for (int k = 0; k < 2; k++)
            {
                long before = GC.GetTotalMemory(forceFullCollection: true);
                BassFit fit = Bass.Fit(BassCost.Config(x, samples: lengths[k]), data, run);
                long after = GC.GetTotalMemory(forceFullCollection: true);
                retained[k] = after - before;
                draws[k] = fit.Draws.Sum(c => c.Count);
                rows.Add(x.Key, $"retained bytes, {draws[k]} draws", retained[k]);
                rows.Add(x.Key, $"mean M over the draws, {draws[k]} draws", fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length));
                GC.KeepAlive(fit);
            }
            double perDraw = (retained[1] - retained[0]) / (draws[1] - draws[0]);
            rows.Add(x.Key, "retained bytes per draw", perDraw);
            rows.Add(x.Key, "retained bytes per fit (no draws)", retained[0] - perDraw * draws[0]);
            _out.WriteLine($"{x.Key}: {perDraw:F0} B per draw, {retained[0] - perDraw * draws[0]:F0} B per fit");
        }
        rows.Flush();
    }

    // The peak's fits burn 10,000 sweeps, past what the burn-in rule took at every largest point but M_max = 200's.
    private const int PeakBurnIn = 10000;

    /// <summary>The peak working set at the axis <c>BASS_COST_PEAK</c> names (<see cref="BassCost.Largest"/>), read in
    /// this process alone: the axis's largest point fitted, two chains, one thread.</summary>
    [Fact]
    public void Counts_PeakWorkingSet()
    {
        if (!BassQualification.Enabled || Environment.GetEnvironmentVariable("BASS_COST_PEAK") is not { Length: > 0 } axis) return;
        BassCost.Point x = BassCost.Largest.Single(a => a.Axis == axis).Point;
        var rows = new BassCost.Rows("counts", $"peak {axis}");
        using var self = Process.GetCurrentProcess();
        self.Refresh();
        rows.Add(x.Key, "peak working set before the fit, bytes", self.PeakWorkingSet64);
        var run = new BassRunOptions { MaxDegreeOfParallelism = 1 };
        BassConfig config = BassCost.Config(x, burnIn: PeakBurnIn, samples: BassCost.Scaled(1000));
        double meanM;
        if (x.Surface == BassCost.Surface.Basis)
        {
            (_, BassPredictor[] predictors) = BassCost.BasisSurface(x.N, x.Locations, x.P);
            BassBasisFit fit = BassComponents.Fit(BassCost.Decomposition(x), predictors, config, run);
            meanM = fit.Components.SelectMany(c => c.Draws.SelectMany(d => d)).Average(d => d.Structure.TermIds.Length);
            GC.KeepAlive(fit);
        }
        else
        {
            BassFit fit = Bass.Fit(config, BassCost.Data(x), run);
            meanM = fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
            GC.KeepAlive(fit);
        }
        self.Refresh();
        rows.Add(x.Key, "peak working set, bytes", self.PeakWorkingSet64);
        rows.Add(x.Key, "mean M over the draws", meanM);
        rows.Add(x.Key, "GC heap after a full collection, bytes", GC.GetTotalMemory(forceFullCollection: true));
        rows.Flush();
        _out.WriteLine($"{axis}: {x.Key} peak {self.PeakWorkingSet64 / 1048576.0:F0} MiB, mean M {meanM:F1}");
    }

    // ---- Times ----

    /// <summary>The timed blocks <c>BASS_COST_BLOCK</c> names (separated by semicolons, since axis names hold commas;
    /// every block when unset), each with its load recorded at its start and end.</summary>
    [Fact]
    public void Timed_Blocks()
    {
        if (!BassQualification.Enabled) return;   // heavy
        (string Block, BassCost.Point[] Points)[] grid = GridBlocks();
        string[] all = [.. grid.Select(g => g.Block), .. OtherBlocks];
        string[] asked = Environment.GetEnvironmentVariable("BASS_COST_BLOCK") is { Length: > 0 } list
            ? [.. list.Split(';', StringSplitOptions.TrimEntries | StringSplitOptions.RemoveEmptyEntries)]
            : all;
        foreach (string block in asked)
        {
            Assert.Contains(block, all);
            var rows = new BassCost.Rows("timed", block);
            foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
            var load = new BassCost.Load();
            BassCostSweep.Read(BassCost.Base with { N = 200 }, BassCost.Meter.Time);   // warm-up: JIT and caches
            switch (block)
            {
                case "repeat": Repeat(rows); break;
                case "outputs": Outputs(rows); break;
                case "backends": Backends(rows); break;
                case "calibration": Calibration(rows); break;
                case "parallel": Parallel(rows); break;
                default:
                    foreach (BassCost.Point x in grid.Single(g => g.Block == block).Points)
                    {
                        BassCostSweep.Reading r = BassCostSweep.Read(x, BassCost.Meter.Time);
                        Write(rows, r, BassCost.Meter.Time, block);
                        _out.WriteLine($"{x.Key}: M {r.M}/{r.Cap}, {r.Of("sweep"):F1} µs per sweep, rest {r.PerSweep()[^1].Value:F1}");
                    }
                    break;
            }
            (string[] end, double cpu, double seconds, double others) = load.End();
            bool flagged = load.AtStart.Length > 0 || end.Length > 0;
            rows.Add("load", "other test hosts at the start", string.Join(" ", load.AtStart));
            rows.Add("load", "other test hosts at the end", string.Join(" ", end));
            rows.Add("load", "own CPU over wall", cpu);
            rows.Add("load", "other processes' cores busy", others);
            rows.Add("load", "block seconds", seconds);
            rows.Add("load", "flagged", flagged ? "yes" : "no");
            _out.WriteLine($"{block}: {seconds:F0} s, CPU/wall {cpu:F2}, others {others:F2} cores, {(flagged ? "FLAGGED" : "quiet")}; rows {rows.Flush()}");
        }
    }

    // A reading's rows: the state, the frequencies, every part with its spread, per sweep, the analytic counts and, in
    // time, each part's time per counted unit.
    private static void Write(BassCost.Rows rows, BassCostSweep.Reading r, BassCost.Meter meter, string axis)
    {
        string p = r.Point.Key, unit = meter == BassCost.Meter.Time ? "µs" : "bytes";
        rows.Add(p, "axis", axis);
        rows.Add(p, "rows N", r.Point.Rows);
        rows.Add(p, "burn-in sweeps", r.BurnSweeps);
        rows.Add(p, "burn-in seconds", r.BurnSeconds);
        rows.Add(p, "read seconds", r.ReadSeconds);
        rows.Add(p, "M at the state", r.M);
        rows.Add(p, "M_max resolved", r.Cap);
        rows.Add(p, "at the cap", r.AtCap ? "yes" : "no");
        rows.Add(p, "mean M at the sweep batches' ends", r.SweepEndM);
        rows.Add(p, "mean order held", r.HeldOrder);
        rows.Add(p, "mean order proposed", r.PoolOrder);
        rows.Add(p, "empty birth share", r.EmptyShare);
        rows.Add(p, "rebuild interval", r.Rebuild);
        rows.Add(p, "parts unread at the state", string.Join(" ", r.Unread));
        (string, double)[] frequencies =
        [
            ("births per sweep", r.BirthsPerSweep), ("births with a term per sweep", r.NonEmptyPerSweep), ("deaths per sweep", r.DeathsPerSweep),
            ("changes per sweep", r.ChangesPerSweep), ("admitted per sweep", r.AdmittedPerSweep), ("refused per sweep", r.RefusedPerSweep),
            ("birth accepts per sweep", r.BirthAccepts), ("death accepts per sweep", r.DeathAccepts), ("change accepts per sweep", r.ChangeAccepts),
            ("rebuilds per sweep", r.RebuildsPerSweep), ("τ evaluations per step", r.TauEvaluations), ("columns evaluated per sweep", r.ColumnsPerSweep),
        ];
        foreach ((string q, double v) in frequencies) rows.Add(p, q, v);
        foreach (string name in r.Parts.Keys)
        {
            rows.Add(p, $"{name} {unit}", r.Of(name));
            rows.Add(p, $"{name} {unit} IQR", r.IqrOf(name));
        }
        rows.Add(p, $"commit birth {unit}, rebuild amortized", r.Commit("birth"));
        rows.Add(p, $"commit death {unit}, rebuild amortized", r.Commit("death"));
        rows.Add(p, $"commit change {unit}, rebuild amortized", r.Commit("change"));
        foreach ((string name, double v) in r.PerSweep()) rows.Add(p, $"per sweep {unit}: {name}", v);
        if (meter == BassCost.Meter.Bytes) rows.Add(p, "bytes per sweep on the thread", r.ThreadBytesPerSweep);
        rows.Add(p, "GC pause µs per sweep", r.GcPausePerSweep);
        rows.Add(p, "gen0 collections per 1,000 sweeps", r.Gen0PerThousand);
        rows.Add(p, "gen1 collections per 1,000 sweeps", r.Gen1PerThousand);
        rows.Add(p, "gen2 collections per 1,000 sweeps", r.Gen2PerThousand);
        foreach ((string name, double units) in r.Units)
        {
            rows.Add(p, $"units: {name}", units);
            double value = name == "commit" ? r.Commit("birth") : r.Of(name);
            if (meter == BassCost.Meter.Time && units > 0 && double.IsFinite(value)) rows.Add(p, $"ns per unit: {name}", value * 1000 / units);
        }
        if (double.IsFinite(r.MembershipCases) && meter == BassCost.Meter.Time) rows.Add(p, "membership ns per case", r.Of("membership") * 1000 / r.MembershipCases);
    }

    // The repeatability check: the base point read twice in one process; its counts agree exactly, its times within the
    // batches' interquartile range.
    private void Repeat(BassCost.Rows rows)
    {
        BassCostSweep.Reading a = BassCostSweep.Read(BassCost.Base, BassCost.Meter.Time), b = BassCostSweep.Read(BassCost.Base, BassCost.Meter.Time);
        Write(rows, a, BassCost.Meter.Time, "repeat 1");
        Write(rows, b, BassCost.Meter.Time, "repeat 2");
        Assert.Equal(a.M, b.M);
        Assert.Equal(a.BirthsPerSweep, b.BirthsPerSweep);
        Assert.Equal(a.AdmittedPerSweep, b.AdmittedPerSweep);
        Assert.Equal(a.BirthAccepts + a.DeathAccepts + a.ChangeAccepts, b.BirthAccepts + b.DeathAccepts + b.ChangeAccepts);
        int within = 0;
        foreach (string name in a.Parts.Keys)
        {
            double gap = Math.Abs(a.Of(name) - b.Of(name)), iqr = Math.Max(a.IqrOf(name), b.IqrOf(name));
            bool agrees = gap <= iqr;
            if (agrees) within++;
            rows.Add("repeat", $"{name}: |difference| µs, IQR µs, within", $"{gap:G4}; {iqr:G4}; {(agrees ? "yes" : "no")}");
        }
        rows.Add("repeat", "parts within their IQR", $"{within} of {a.Parts.Count}");
        _out.WriteLine($"repeat: {within} of {a.Parts.Count} parts within their IQR");
    }

    // ---- The outputs ----

    // Prediction per case per draw, scalar and functional; Sobol indices and effects per draw; for basis fits against r.
    private void Outputs(BassCost.Rows rows)
    {
        var run = new BassRunOptions { MaxDegreeOfParallelism = 1 };
        int samples = BassCost.Scaled(250);
        BassCost.Point b = BassCost.Base;
        BassCost.Point[] scalar =
        [
            b, .. new[] { 5, 10, 50, 100, 200 }.Select(c => b with { Cap = c }), .. new[] { 2, 5, 20, 50, 100, 200 }.Select(p => b with { P = p }),
            .. new[] { 1, 3, 4, 5 }.Select(k => b with { Order = k }),
            .. new[] { 10, 50, 200, 1000 }.Select(m => new BassCost.Point(BassCost.Surface.Functional, 100, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: m)),
        ];
        foreach (BassCost.Point x in scalar)
        {
            BassData data = BassCost.Data(x);
            BassFit fit = Bass.Fit(BassCost.Config(x, samples: samples), data, run);
            int draws = fit.Draws.Sum(c => c.Count);
            int cases = x.Surface == BassCost.Surface.Functional ? 50 : 200;
            PredictionCases at = PredictionCases.Of(x.Surface == BassCost.Surface.Functional ? BassCost.Functional(cases, x.Locations, x.P).Predictors : BassCost.Rich(cases, x.P).Predictors);
            // The effects on a 20-point grid: R's 100 holds 100² cells per held pair, past the default budget at the full
            // draw count (61,000 cells × 500 draws at the base point) and near 2 GB at the cap points; the cost per cell
            // carries to any grid.
            var effectSpec = new EffectSpec { GridPoints = 20, CellBudget = 50_000_000 };
            EffectAnalysis probe = BassSensitivity.Effects(fit, data, effects: effectSpec);
            long cells = probe.At.Sum(f => f.Effects.Sum(e => (long)e.Values.Length)) / draws;
            Dictionary<string, double[]> times = BassCost.Interleave(
            [
                ("predict", BassCost.Batch(BassCost.Meter.Time, () => fit.Predict(data, at), batch: 1)),
                ("sobol", BassCost.Batch(BassCost.Meter.Time, () => BassSensitivity.Of(fit, data), batch: 1)),
                ("effects", BassCost.Batch(BassCost.Meter.Time, () => BassSensitivity.Effects(fit, data, effects: effectSpec), batch: 1)),
            ]);
            double locations = Math.Max(1, x.Locations);
            rows.Add(x.Key, "draws", draws);
            rows.Add(x.Key, "mean M over the draws", fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length));
            rows.Add(x.Key, "prediction cases", cases);
            rows.Add(x.Key, "predict µs per case per location per draw", BassCost.Median(times["predict"]) / (cases * locations * draws));
            rows.Add(x.Key, "sobol µs per draw", BassCost.Median(times["sobol"]) / draws);
            rows.Add(x.Key, "effects µs per draw (20-point grid)", BassCost.Median(times["effects"]) / draws);
            rows.Add(x.Key, "effect cells (20-point grid)", cells);
            rows.Add(x.Key, "effects ns per cell per draw", BassCost.Median(times["effects"]) * 1000 / (draws * (double)cells));
            foreach ((string name, double[] v) in times) rows.Add(x.Key, $"{name} µs per call IQR", BassCost.Iqr(v));
            _out.WriteLine($"outputs {x.Key}: {draws} draws, predict {BassCost.Median(times["predict"]) / (cases * locations * draws):G3} µs/case/draw, " +
                           $"sobol {BassCost.Median(times["sobol"]) / draws:G3} µs/draw, effects {BassCost.Median(times["effects"]) / draws:G3} µs/draw");
            rows.Flush();
        }

        // Basis fits against r: the components' fit, prediction at 50 new cases over the m locations, and the Sobol
        // decomposition, per joint draw.
        foreach (int r in new[] { 1, 3, 10, 30 })
        {
            var x = new BassCost.Point(BassCost.Surface.Basis, 500, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: 200, Rank: r);
            (_, BassPredictor[] predictors) = BassCost.BasisSurface(x.N, x.Locations, x.P);
            BassDecomposition decomposition = BassCost.Decomposition(x);
            BassConfig config = BassCost.Config(x, samples: samples);
            var clock = Stopwatch.StartNew();
            BassBasisFit fit = BassComponents.Fit(decomposition, predictors, config, run);
            double fitSeconds = clock.Elapsed.TotalSeconds;
            int draws = fit.Components[0].Draws.Sum(c => c.Count);
            PredictionCases at = PredictionCases.Of(BassCost.BasisSurface(50, x.Locations, x.P).Predictors);
            Dictionary<string, double[]> times = BassCost.Interleave(
            [
                ("predict", BassCost.Batch(BassCost.Meter.Time, () => fit.Predict(at), batch: 1)),
                ("sobol", BassCost.Batch(BassCost.Meter.Time, () => BassBasisSensitivity.Of(fit), batch: 1)),
            ]);
            long sweeps = (long)r * config.Sampling.Chains * (config.Sampling.BurnIn + config.Sampling.Samples);
            rows.Add(x.Key, "draws", draws);
            rows.Add(x.Key, "fit seconds", fitSeconds);
            rows.Add(x.Key, "fit µs per component sweep", fitSeconds * 1e6 / sweeps);
            rows.Add(x.Key, "predict µs per case per location per draw", BassCost.Median(times["predict"]) / (50.0 * x.Locations * draws));
            rows.Add(x.Key, "basis sobol µs per draw", BassCost.Median(times["sobol"]) / draws);
            foreach ((string name, double[] v) in times) rows.Add(x.Key, $"{name} µs per call IQR", BassCost.Iqr(v));
            _out.WriteLine($"outputs {x.Key}: fit {fitSeconds:F1} s, sobol {BassCost.Median(times["sobol"]) / draws:G3} µs/draw");
            rows.Flush();
        }
    }

    // ---- The other backends ----

    private void Backends(BassCost.Rows rows)
    {
        (BassModel model, BassChain chain, _) = BassCost.Burned(BassCost.Base);
        string p = BassCost.Base.Key;
        // The untempered sweep (moves and τ, no composition) at the state, for scale.
        BassChainSnapshot fixedChain = chain.Snapshot();
        BassChain plain = BassChain.Restore(model, fixedChain);
        Dictionary<string, double[]> untempered = BassCost.Interleave([("step", BassCost.Batch(BassCost.Meter.Time, () => plain.Step(), batch: 100))]);
        rows.Add(p, "untempered step µs (moves and τ)", BassCost.Median(untempered["step"]));

        Tempering(rows, model);
        Smc(rows, model);
        FlatHistogram(rows, model, chain);
    }

    // Tempering at 5, 10 and 20 rungs, geometric from b = 1 to the self-sized ladder's hot end: per rung-sweep and per
    // swap round (every pair's ratio, and the hand-overs of the accepted swaps), M1's tempered decomposition.
    private void Tempering(BassCost.Rows rows, BassModel model)
    {
        (BassState start, int[] seeds) = BassMixing.Start(model, 0);
        double hot = new BassLadder(model, new TemperingSpec(), start.Terms, start.LogTau, seeds[3]).Tempering.InverseTemperatures.Min();
        rows.Add("tempering", "hot end b", hot);
        foreach (int levels in new[] { 5, 10, 20 })
        {
            double[] betas = [.. Enumerable.Range(0, levels).Select(k => Math.Pow(hot, k / (levels - 1.0)))];
            var spec = new TemperingSpec { Levels = levels, Adaptive = false };
            var ladder = new BassLadder(model, spec, betas, start.Terms, start.LogTau, seeds[3]);
            BassMixing.Burn(ladder, BassCost.BurnIn);
            ParallelTempering<BassState> tempering = ladder.Tempering;
            int pairs = levels - 1;
            long swaps = 0, scans = 0;
            double sink = 0.0;
            int i = 0;
            int[] givers = [0, levels / 2, levels - 1];
            int scanBatch = BassCost.Scaled(20);
            var timers = new List<(string Name, Func<int, double> Batch)>
            {
                ("scan", round =>
                {
                    long before = tempering.SwapAccepts;
                    double micro = BassCost.Batch(BassCost.Meter.Time, () => ladder.Step(), batch: 20)(round);
                    if (round >= 0)
                    {
                        swaps += tempering.SwapAccepts - before;
                        scans += scanBatch;
                    }
                    return micro;
                }),
                ("ratio", BassCost.Batch(BassCost.Meter.Time, () =>
                {
                    var a = tempering.Replicas[i % pairs];
                    var b = tempering.Replicas[i++ % pairs + 1];
                    sink += a.Temperature.LogDensity(b.CurrentEvaluation) + b.Temperature.LogDensity(a.CurrentEvaluation) - a.CurrentLogDensity - b.CurrentLogDensity;
                })),
            };
            foreach (int j in givers)
                timers.Add(($"hand-over {j}", round =>
                {
                    BassState state = tempering.Replicas[j].Current;
                    BassStructureMoves giver = ladder.Structures[j];
                    var receiver = new BassStructureMoves(model, state, giver.Snapshot(state));
                    return BassCost.Batch(BassCost.Meter.Time, () => receiver.Adopt(state, giver.DesignFor(state)))(round);
                }));
            Dictionary<string, double[]> times = BassCost.Interleave(timers);
            double scan = BassCost.Median(times["scan"]), ratio = BassCost.Median(times["ratio"]), accepted = (double)swaps / scans;
            double handOver = givers.Average(j => BassCost.Median(times[$"hand-over {j}"]));
            double swapRound = pairs * ratio + 2 * accepted * handOver;
            string key = $"tempering, {levels} rungs";
            rows.Add(key, "scan µs", scan);
            rows.Add(key, "scan µs IQR", BassCost.Iqr(times["scan"]));
            rows.Add(key, "pair ratio µs", ratio);
            rows.Add(key, "hand-over µs", handOver);
            rows.Add(key, "accepted swaps per scan", accepted);
            rows.Add(key, "swap round µs", swapRound);
            rows.Add(key, "rung-sweep µs", (scan - swapRound) / levels);
            rows.Add(key, "M by rung", string.Join(" ", tempering.Replicas.Select(r => r.Current.Count)));
            Assert.True(double.IsFinite(sink));
            _out.WriteLine($"{key}: scan {scan:F0} µs, rung-sweep {(scan - swapRound) / levels:F1}, swap round {swapRound:F1}");
        }
    }

    // SMC for fitting at 1,000 and 4,000 particles, one thread, D05's defaults otherwise: each stage's reweighting search
    // and its advance read whole; then at the run's end, a particle's mutation step, a reweighting, the ancestors' draw and
    // a fork, read in batches.
    private void Smc(BassCost.Rows rows, BassModel model)
    {
        var defaults = new SmcSpec();
        foreach (int particles in new[] { 1000, 4000 })
        {
            var target = new BassParticleTarget(model.Target);
            var source = new BassParticles(model);
            var kernel = new BassSmcKernel(model);
            var options = new SmcOptions
            {
                Particles = particles, Resampling = defaults.Resampling, ResampleBelow = defaults.ResampleBelow, MutationSteps = defaults.MutationSteps,
                MaxStages = defaults.MaxStages, MaxDegreeOfParallelism = 1,
            };
            SmcRun<BassParticle> run = SmcRun<BassParticle>.Start(new SmcStage<BassParticle>(target, new Temperature(0.0, TemperingPath.Likelihood)),
                                                                  source, kernel, options, 211);
            int stageCap = BassCost.Scaled(defaults.MaxStages);
            var search = new List<double>();
            var advance = new List<double>();
            while (run.Stage.Temperature.InverseTemperature < 1.0 && run.StagesDone < stageCap)
            {
                long t0 = Stopwatch.GetTimestamp();
                double next = SmcSampler.NextInverseTemperature(run, defaults.TargetConditionalEss);
                long t1 = Stopwatch.GetTimestamp();
                run.Advance(new SmcStage<BassParticle>(target, new Temperature(next, TemperingPath.Likelihood)));
                long t2 = Stopwatch.GetTimestamp();
                search.Add((t1 - t0) * 1e6 / Stopwatch.Frequency);
                advance.Add((t2 - t1) * 1e6 / Stopwatch.Frequency);
            }
            SmcResult<BassParticle> result = run.Result();
            SmcStage<BassParticle> last = run.Stage;
            IReadOnlyList<BassParticle> held = result.Particles;
            ITargetEvaluation[] evaluations = [.. held.Select(target.Evaluate)];
            var rng = new Xoshiro256PlusPlus(223);
            var ancestors = new int[particles];
            double[] uniform = [.. Enumerable.Repeat(1.0 / particles, particles)];
            int i = 0;
            Dictionary<string, double[]> times = BassCost.Interleave(
            [
                ("mutation step", BassCost.Batch(BassCost.Meter.Time, () =>
                {
                    int k = i++ % particles;
                    evaluations[k] = kernel.Move(held[k], evaluations[k], last, rng).Evaluation;
                }, batch: 200)),
                ("reweighting", BassCost.Batch(BassCost.Meter.Time, () => run.ConditionalEss(last), batch: 20)),
                ("ancestors", BassCost.Batch(BassCost.Meter.Time, () => Resampling.Ancestors(defaults.Resampling, uniform, rng, ancestors), batch: 20)),
                ("fork", BassCost.Batch(BassCost.Meter.Time, () => source.Fork(held[i++ % particles]), batch: 200)),
            ]);
            IReadOnlyList<SmcStageReport> stages = result.Stages;
            SmcStageReport[] resampled = [.. stages.Skip(1).Where(s => s.Resampled)];
            double forks = resampled.Length == 0 ? 0.0 : resampled.Average(s => particles - (double)s.DistinctAncestors);
            string key = $"SMC, {particles} particles";
            rows.Add(key, "stages", stages.Count - 1);
            rows.Add(key, "reached the posterior", result.ReachedPosterior ? "yes" : "no");
            rows.Add(key, "resampled stages", resampled.Length);
            rows.Add(key, "mean forks per resampling", forks);
            rows.Add(key, "mean M at the end", held.Average(q => q.State.Count));
            rows.Add(key, "stage search µs, median", BassCost.Median([.. search]));
            rows.Add(key, "stage advance µs, median", BassCost.Median([.. advance]));
            rows.Add(key, "run seconds (search + advance)", (search.Sum() + advance.Sum()) / 1e6);
            rows.Add(key, "advance µs per particle mutation step, run mean", advance.Sum() / ((stages.Count - 1.0) * particles * defaults.MutationSteps));
            rows.Add(key, "mutation step µs at the end", BassCost.Median(times["mutation step"]));
            rows.Add(key, "reweighting µs at the end", BassCost.Median(times["reweighting"]));
            rows.Add(key, "ancestors µs", BassCost.Median(times["ancestors"]));
            rows.Add(key, "fork µs", BassCost.Median(times["fork"]));
            rows.Add(key, "resampling µs (ancestors + mean forks)", BassCost.Median(times["ancestors"]) + forks * BassCost.Median(times["fork"]));
            foreach ((string name, double[] v) in times) rows.Add(key, $"{name} µs IQR", BassCost.Iqr(v));
            _out.WriteLine($"{key}: {stages.Count - 1} stages, mutation {BassCost.Median(times["mutation step"]):F1} µs, fork {BassCost.Median(times["fork"]):F1} µs");
            rows.Flush();
        }
    }

    // The flat histogram at the base point: its learning iterations whole (four chains stepped one thread, the visits'
    // weight update and the chains' reweighing), from a flat start over a window around the burned state's ℓ; and its
    // weight lookup, update and reweighing alone.
    private void FlatHistogram(BassCost.Rows rows, BassModel model, BassChain chain)
    {
        double top = model.Target.Evaluate(chain.Current).LogLikelihood(1.0);
        double low = model.Target.Evaluate(BassMixing.Start(model, 0).Start).LogLikelihood(1.0);
        var window = new FlatHistogramWindow(Math.Min(low, top) - 50.0, top + 50.0);
        int iterations = BassCost.Scaled(2000);
        var spec = new FlatHistogramSpec
        {
            Chains = 4, LearningIterations = 100_000, ProductionDraws = 1, Start = FlatHistogramStart.Flat, Window = window, LowSlope = 0.0, Nodes = 64,
        };
        string fingerprint = Bass.Fingerprint(model.Config, model.Data);
        var client = new BassMulticanonicalModel(model, fingerprint, spec.LearningIterations, new BassSweepCount(new BassRunOptions(), spec.LearningIterations * spec.Chains));
        var parallel = new System.Threading.Tasks.ParallelOptions { MaxDegreeOfParallelism = 1 };
        MulticanonicalRun.Run(spec, model.Target, client, 229, iterations: BassCost.Scaled(100), parallel: parallel);   // warm-up
        long t0 = Stopwatch.GetTimestamp();
        MulticanonicalResult<BassState> result = MulticanonicalRun.Run(spec, model.Target, client, 229, iterations: iterations, parallel: parallel);
        double perIteration = (Stopwatch.GetTimestamp() - t0) * 1e6 / Stopwatch.Frequency / iterations;

        MulticanonicalWeights weights = result.Weights;
        var learning = new SamcLearning(weights, 10.0 * weights.Nodes);
        double[] levels = [.. Enumerable.Range(0, spec.Chains).Select(c => window.Low + (c + 0.5) / spec.Chains * (window.High - window.Low))];
        var evaluation = new MulticanonicalEvaluation(model.Target.Evaluate(chain.Current), weights);
        double sink = 0.0;
        int i = 0;
        Dictionary<string, double[]> times = BassCost.Interleave(
        [
            ("lookup", BassCost.Batch(BassCost.Meter.Time, () => sink += weights.LogWeight(levels[i++ % levels.Length]))),
            ("update", BassCost.Batch(BassCost.Meter.Time, () => learning.Update(levels))),
            ("reweigh", BassCost.Batch(BassCost.Meter.Time, () => sink += evaluation.Reweighed(weights).LogWeight)),
        ]);
        string key = "flat histogram, 4 chains";
        rows.Add(key, "window", $"{window.Low:G6} to {window.High:G6}, 64 nodes");
        rows.Add(key, "iterations read", iterations);
        rows.Add(key, "iteration µs (4 chain steps, update, reweighs)", perIteration);
        rows.Add(key, "weight lookup µs", BassCost.Median(times["lookup"]));
        rows.Add(key, "weight update µs (one iteration's)", BassCost.Median(times["update"]));
        rows.Add(key, "reweigh µs (per chain)", BassCost.Median(times["reweigh"]));
        rows.Add(key, "chain step µs (the iteration less the update and reweighs, per chain)",
                 (perIteration - BassCost.Median(times["update"]) - spec.Chains * BassCost.Median(times["reweigh"])) / spec.Chains);
        rows.Add(key, "mean M at the end", result.Chains.Average(c => ((BassMulticanonicalChain)c).Chain.Current.Count));
        Assert.True(double.IsFinite(sink));
        _out.WriteLine($"{key}: {perIteration:F1} µs per iteration");
    }

    // ---- Calibration ----

    // N15: DiagonalPlusLowRank's construction and solve per σ² proposal, m ∈ {20, 50, 200} × k ∈ {10, 39, 100}; and the
    // cut SMC per leg at stage 7's sizes.
    private void Calibration(BassCost.Rows rows)
    {
        var rng = new Xoshiro256PlusPlus(233);
        foreach (int m in new[] { 20, 50, 200 })
            foreach (int k in new[] { 10, 39, 100 })
            {
                double[] d = [.. Enumerable.Range(0, m).Select(_ => 0.5 + rng.NextDouble())];
                var u = new double[m, k];
                for (int r = 0; r < m; r++)
                    for (int c = 0; c < k; c++) u[r, c] = BassQualification.Gaussian(rng);
                double[] s = [.. Enumerable.Range(0, k).Select(_ => 0.1 + rng.NextDouble())];
                double[] b = [.. Enumerable.Range(0, m).Select(_ => BassQualification.Gaussian(rng))];
                var fixedSigma = new DiagonalPlusLowRank(d, u, s);
                double sink = 0.0;
                int batch = Math.Clamp(2_000_000 / (m * k * k), 4, 200);
                Dictionary<string, double[]> times = BassCost.Interleave(
                [
                    ("proposal", BassCost.Batch(BassCost.Meter.Time, () =>
                    {
                        var sigma = new DiagonalPlusLowRank(d, u, s);
                        sink += sigma.InverseQuadraticForm(b) + sigma.LogDet;
                    }, batch)),
                    ("solve", BassCost.Batch(BassCost.Meter.Time, () => sink += fixedSigma.InverseQuadraticForm(b), batch: 200)),
                ]);
                int t = Math.Min(m, k);
                // One-sided Jacobi on W (m × k): per sweep k(k − 1)/2 pairs, each two norms and a cosine (3m) and, rotated,
                // the two columns and V's (4m + 4k); the capacitance Cholesky: WᵀW, its factor and two solves.
                double jacobiSweep = k * (k - 1) / 2.0 * (7.0 * m + 4.0 * k);
                double capacitance = m * k * (k + 1) / 2.0 + k * (double)k * k / 6.0 + 2.0 * m * k + k * (double)k;
                string key = $"N15 m={m} k={k}";
                rows.Add(key, "proposal µs (construction and solve)", BassCost.Median(times["proposal"]));
                rows.Add(key, "proposal µs IQR", BassCost.Iqr(times["proposal"]));
                rows.Add(key, "solve µs", BassCost.Median(times["solve"]));
                rows.Add(key, "rank t", t);
                rows.Add(key, "Jacobi multiply-adds per sweep", jacobiSweep);
                rows.Add(key, "capacitance Cholesky multiply-adds per proposal", capacitance);
                Assert.True(double.IsFinite(sink));
                _out.WriteLine($"{key}: {BassCost.Median(times["proposal"]):F1} µs per proposal");
            }
        rows.Flush();

        // The cut SMC at stage 7's pooled-posterior fixture and its default spec, one thread.
        var walls = new List<double>();
        CutSmcCalibration? last = null;
        for (int rep = 0; rep < Math.Max(1, (int)Math.Round(3 * BassCost.Scale)); rep++)
        {
            CalibrationModel model = CutSmcModel(100 + rep, new CutSmcSpec());
            var clock = Stopwatch.StartNew();
            last = CutSmcSampler.Run(model, 1);
            walls.Add(clock.Elapsed.TotalSeconds);
        }
        rows.Add("cut SMC, stage 7", "legs", last!.Populations.Count);
        rows.Add("cut SMC, stage 7", "stages", last.TotalStages);
        rows.Add("cut SMC, stage 7", "particles", new CutSmcSpec().Particles);
        rows.Add("cut SMC, stage 7", "run seconds, median", BassCost.Median([.. walls]));
        rows.Add("cut SMC, stage 7", "seconds per leg", BassCost.Median([.. walls]) / last.Populations.Count);
        rows.Add("cut SMC, stage 7", "ms per stage", 1000 * BassCost.Median([.. walls]) / last.TotalStages);
        rows.Flush();
    }

    // BassCalibrationCutSmcTests' stage-7 fixture: one θ, a categorical of 3 levels, one inverse-gamma group, the
    // discrepancy, three draws.
    private static CalibrationModel CutSmcModel(int seed, CutSmcSpec smc)
    {
        TestEmulator emulator = CalibrationFixture.Emulator(1);
        CalibrationConfig config = CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true)) with
        {
            Sampling = new CalibrationSampling { MasterSeed = seed },
            Smc = smc,
        };
        return new CalibrationModel(config, CalibrationFixture.Observations(emulator), emulator);
    }

    // ---- Parallel fits ----

    // At three points, the default fit (four chains, default lengths) with its chains in parallel, against the same fit
    // one chain at a time: one chain's time × the chain count.
    private void Parallel(BassCost.Rows rows)
    {
        BassCost.Point[] points =
        [
            BassCost.Base, BassCost.Base with { N = 10000 },
            new(BassCost.Surface.Functional, 100, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: 200),
        ];
        foreach (BassCost.Point x in points)
        {
            BassData data = BassCost.Data(x);
            var sampling = new BassSampling();
            BassConfig config = BassCost.Config(x, burnIn: BassCost.Scaled(sampling.BurnIn), samples: BassCost.Scaled(sampling.Samples), chains: sampling.Chains);
            var clock = Stopwatch.StartNew();
            BassFit serial = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
            double one = clock.Elapsed.TotalSeconds;
            clock.Restart();
            BassFit together = Bass.Fit(config, data, new BassRunOptions());
            double all = clock.Elapsed.TotalSeconds;
            int chains = config.Sampling.Chains, cores = Math.Min(chains, Environment.ProcessorCount);
            rows.Add(x.Key, "chains", chains);
            rows.Add(x.Key, "serial fit seconds (one chain's time × the chains)", one);
            rows.Add(x.Key, "parallel fit seconds", all);
            rows.Add(x.Key, "speedup", one / all);
            rows.Add(x.Key, "parallel efficiency over min(chains, processors)", one / all / cores);
            rows.Add(x.Key, "draws identical", serial.Draws.SelectMany(c => c).Zip(together.Draws.SelectMany(c => c)).All(z => z.First.Sigma2 == z.Second.Sigma2) ? "yes" : "no");
            _out.WriteLine($"parallel {x.Key}: {one:F1} s serial, {all:F1} s parallel, speedup {one / all:F2}");
            rows.Flush();
        }
    }

    private static string OutputFolder() => Environment.GetEnvironmentVariable("BASS_COST_OUT") is { Length: > 0 } d ? d : Path.Combine(Path.GetTempPath(), "bass-cost");
}
