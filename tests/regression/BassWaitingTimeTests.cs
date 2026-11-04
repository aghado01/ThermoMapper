using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using F = Maths.Regression.Tests.BassPathProfile;

namespace Maths.Regression.Tests;

/// <summary>
/// The route waiting times against exact hitting times and first passages (the waiting-time brief; gap G1, reading (d)).
/// <list type="bullet">
/// <item>Always run: the reduction against closed forms and a dense solve, the series and its cut against the path profile's reading, Rayleigh on an enumerated kernel, and the chain's record on one small passage.</item>
/// <item>Qualification runs (<see cref="BassQualification.Enabled"/>): the parts, each its own test, heavy ones in the brief's order. Rows go beside the brief as <c>waiting-time-*.csv</c>, appended as reached and skipped on a rerun.</item>
/// </list>
/// Choices the brief leaves open:
/// <list type="bullet">
/// <item>The chains' 99% interval on a mean of at least 1,000 passages, by the mean's normal law: about 5% across the checks together.</item>
/// <item>"The series cut at the route's first state in the basin" keeps a slide the route takes first, whose every state is the start's basin, since the cut says where the series ends, not how it begins. <c>partial</c> and <c>exact</c> beside <c>aligned</c> read other starts of the same series.</item>
/// <item>The mass from a chain reads the highest structure of the start's basin seen before first entry, its dwell counted from when first reached (<see cref="PassageTally"/>): an estimator of the same ratio that does not charge the climb to it.</item>
/// <item>The ladder's grid in b is 0.005 apart; a rung is the grid's b nearest the brief's ln. The rungs run lowest ln first, so a budget that stops leaves the costliest ones.</item>
/// <item>At M ≤ 4 the dense check runs on the exact targets alone: its cost is about a minute a target, and the set targets share the same arithmetic.</item>
/// </list>
/// </summary>
public sealed class BassWaitingTimeTests
{
    private const double Rung10 = 10.0;
    private static readonly double[] LadderNats = [Rung10, 12.0, 14.0, 16.0];
    private const long Million = 1_000_000;

    private static readonly string[] Shallow = ["a2", "a1x"], Deep = ["A", "B", "chain0", "chain1", "chain2", "chain3"];
    private static readonly Dictionary<string, string> PairName = new()
    {
        ["A"] = "A-a1", ["B"] = "B-a1", ["a2"] = "a2-a1", ["a1x"] = "a1+x6x8-a1", ["chain0"] = "chain0-a1", ["chain1"] = "chain1-a1", ["chain2"] = "chain2-a1", ["chain3"] = "chain3-a1",
    };

    private static readonly string[] ExactHeader =
    [
        "cap", "seed", "from", "to", "target", "m_from", "m_to", "states", "target_states", "ln_h", "ln_w_eff", "ln_m_ratio", "ln_excursion", "ln_basin_ratio", "v_basin", "dense_ln_h",
        "exact_ln_w", "exact_ln_w_class", "exact_depth", "exact_highest_m", "exact_hops", "exact_ln_w_lattice",
        "aligned_ln_w", "aligned_ln_w_class", "aligned_depth", "aligned_highest_m", "aligned_hops", "aligned_slides", "aligned_ln_w_lattice", "aligned_position", "note",
    ];

    private static readonly string[] LatticeHeader =
    [
        "n", "data", "pair", "exemplar", "alignment", "direction", "d", "position", "slides", "states_admitted", "pi_from", "ln_w_route", "ln_w_lattice", "ln_credit", "ln_w_route_cut", "ln_w_lattice_cut", "ln_credit_cut", "stored_ln_w", "seconds",
    ];

    private static readonly string[] ReadingHeader =
    [
        "part", "start", "target", "alignment", "pi_from", "ln_w_full", "ln_w_class_full", "ln_w_cut", "ln_w_class_cut", "hops_cut", "slides_cut", "low_cut", "depth_cut", "highest_m_cut", "stored_ln_w", "stored_ln_w_class",
    ];

    private static readonly string[] ValidationHeader = ["check", "case", "count", "worst", "value_a", "value_b", "value_c", "value_d", "note"];

    private static readonly string[] LadderHeader = ["pair", "b", "position", "pi_from", "ln_w_full", "ln_w_class_full", "ln_w_cut", "ln_w_class_cut", "low_cut", "low_full", "hops_cut", "highest_m_cut"];

    private static readonly string[] RungHeader = ["pair", "nat", "b", "ln_w_cut", "ln_w_class_cut", "off", "cap_sweeps", "duplicate_of"];

    private static readonly string[] FlatHeader = ["label", "basin", "share", "stored_share", "ln_mass_over_top", "structure_draws", "basin_draws", "pi_top", "ln_mass", "levels_note"];

    private readonly ITestOutputHelper _out;
    public BassWaitingTimeTests(ITestOutputHelper output) => _out = output;

    private static PathCsv Csv(string file, string[] header, int keys) => new(file, header, keys, "waiting-time");
    private static PathCsv Validation() => Csv("validation", ValidationHeader, 2);
    private static void Quietly() => Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
    private static double D(string v) => double.Parse(v, CultureInfo.InvariantCulture);
    private static string Stamp => DateTime.UtcNow.ToString("yyyy-MM-ddTHH:mm", CultureInfo.InvariantCulture);

    // ---- Always run ----

    [Theory]
    [InlineData(0.3, 0.1)]
    [InlineData(1e-40, 1e-60)]
    public void TwoStates_ReadOneOverTheRate(double a, double b)
    {
        double logTotal = LogMath.Add(Math.Log(a), Math.Log(b));
        var flows = new Dictionary<(int, int), double> { [(0, 1)] = Math.Log(a), [(1, 0)] = Math.Log(b) };
        var chain = ReversibleChain.Of([Math.Log(b) - logTotal, Math.Log(a) - logTotal], [flows]);
        (HittingReading r, double[] v) = new HittingReduction(chain, [0, 1]).Toward([1], [0]).Read(0)!.Value;
        Assert.Equal(-Math.Log(a), r.LnHit, 10);
        Assert.Equal(r.LnHit, r.LnHitFromVoltage, 10);
        Assert.Equal(1.0, v[0]);
        Assert.Equal(0.0, v[1]);
    }

    // The lazy reflecting path, π uniform: R_eff = 4n(n − 1), M = ½, so h = 2n(n − 1).
    [Fact]
    public void Path_ReadsItsClosedForm()
    {
        const int n = 40;
        var flows = new Dictionary<(int, int), double>();
        for (int i = 0; i + 1 < n; i++) flows[(i, i + 1)] = flows[(i + 1, i)] = Math.Log(0.25);
        var chain = ReversibleChain.Of(new double[n], [flows]);
        (HittingReading r, double[] v) = new HittingReduction(chain, [0, n - 1]).Toward([n - 1], [0]).Read(0)!.Value;
        Assert.Equal(Math.Log(2.0 * n * (n - 1)), r.LnHit, 10);
        Assert.Equal(r.LnHit, r.LnHitFromVoltage, 10);
        for (int i = 0; i < n; i++) Assert.Equal(1.0 - (double)i / (n - 1), v[i], 12);
    }

    // Two wells of two states joined by a rate ε near 10⁻⁴⁰, far past a plain solve: R_eff = 16 + 4/ε, M = ½.
    [Fact]
    public void TwoWells_ReadATinyRate()
    {
        const double epsilon = 1e-40;
        var flows = new Dictionary<(int, int), double>
        {
            [(0, 1)] = Math.Log(0.5), [(1, 0)] = Math.Log(0.5), [(2, 3)] = Math.Log(0.5), [(3, 2)] = Math.Log(0.5),
            [(1, 2)] = Math.Log(epsilon), [(2, 1)] = Math.Log(epsilon),
        };
        var chain = ReversibleChain.Of(new double[4], [flows]);
        HittingReading r = new HittingReduction(chain, [0, 3]).Toward([3], [0]).Read(0)!.Value.Reading;
        Assert.Equal(Math.Log(8.0 + 2.0 / epsilon), r.LnHit, 10);
        Assert.Equal(r.LnHit, r.LnHitFromVoltage, 10);
        Assert.Equal(Math.Log(0.25 * (16.0 + 4.0 / epsilon)), r.LnWaitEffective, 10);
    }

    // A random reversible chain whose π spans 30 nats: exact and set targets, shared stages, against the dense solve; the kernel's walker
    // against both; and Rayleigh on the chain without a share of its edges.
    [Fact]
    public void Reduction_AgreesWithTheDenseSolve()
    {
        const int n = 60;
        var rng = new Xoshiro256PlusPlus(171);
        double[] logPi = [.. Enumerable.Range(0, n).Select(_ => 30 * rng.NextDouble())];
        var flows = new Dictionary<(int, int), double>();
        var thinned = new Dictionary<(int, int), double>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                bool backbone = j == i + 1;
                if (!backbone && rng.NextDouble() > 0.15) continue;
                double q = Math.Log((0.1 + 0.9 * rng.NextDouble()) / n);
                double up = q + Math.Min(0.0, logPi[j] - logPi[i]), down = q + Math.Min(0.0, logPi[i] - logPi[j]);
                (flows[(i, j)], flows[(j, i)]) = (up, down);
                if (!backbone) continue;
                (thinned[(i, j)], thinned[(j, i)]) = (up, down);
            }
        var chain = ReversibleChain.Of(logPi, [flows]);
        int[] kept = [3, 11, 27, 40, 41, 42, 59];
        var reduction = new HittingReduction(chain, kept, threads: 2);
        int compared = 0;
        foreach (int[] target in new[] { new[] { 59 }, new[] { 3 }, new[] { 40, 41, 42 } })
        {
            double[] dense = HittingReduction.Dense(chain, target);
            int[] starts = [.. kept.Except(target)];
            HittingReduction.Target stage = reduction.Toward(target, starts);
            foreach (int s in starts)
            {
                (HittingReading r, double[] v) = stage.Read(s)!.Value;
                // A plain solve holds six digits to about e²⁰; past it the voltage's sum is the reference.
                if (r.LnHit <= 20.0)
                {
                    Assert.True(Math.Abs(r.LnHit - Math.Log(dense[s])) <= 1e-7, $"h from {s}: {r.LnHit} against {Math.Log(dense[s])}");
                    compared++;
                }
                Assert.True(Math.Abs(r.LnHit - r.LnHitFromVoltage) <= 1e-9, $"R_eff · M from {s}: {r.LnHitFromVoltage} against {r.LnHit}");
                Assert.True(r.LnMassRatio >= 0.0);
                Assert.All(target, t => Assert.Equal(0.0, v[t]));
                // Harmonic off the ends.
                for (int x = 0; x < n; x++)
                {
                    if (x == s || target.Contains(x)) continue;
                    double next = chain.Hold[x] * v[x];
                    for (int e = chain.Starts[x]; e < chain.Starts[x + 1]; e++) next += chain.Rate[e] * v[chain.To[e]];
                    Assert.True(Math.Abs(next - v[x]) <= 1e-9 * Math.Max(v[x], 1e-100) + 1e-30, $"v is not harmonic at {x}: {next} against {v[x]}");
                }
            }
        }
        Assert.True(compared > 0, "No start read under e²⁰.");
        var thin = ReversibleChain.Of(logPi, [thinned]);
        double whole = new HittingReduction(chain, [0, 59]).Toward([59], [0]).Read(0)!.Value.Reading.LnWaitEffective;
        double sparse = new HittingReduction(thin, [0, 59]).Toward([59], [0]).Read(0)!.Value.Reading.LnWaitEffective;
        Assert.True(whole <= sparse + 1e-12, $"Rayleigh: the whole chain reads {whole}, the thinned {sparse}.");
    }

    [Fact]
    public void Walker_ReadsTheExactMean()
    {
        var rng = new Xoshiro256PlusPlus(173);
        const int n = 12;
        double[] logPi = [.. Enumerable.Range(0, n).Select(i => -0.3 * Math.Abs(i - 5.0))];
        var flows = new Dictionary<(int, int), double>();
        for (int i = 0; i + 1 < n; i++)
        {
            flows[(i, i + 1)] = Math.Log(0.2) + Math.Min(0.0, logPi[i + 1] - logPi[i]);
            flows[(i + 1, i)] = Math.Log(0.2) + Math.Min(0.0, logPi[i] - logPi[i + 1]);
        }
        var chain = ReversibleChain.Of(logPi, [flows]);
        double h = Math.Exp(new HittingReduction(chain, [0, n - 1]).Toward([n - 1], [0]).Read(0)!.Value.Reading.LnHit);
        var target = new bool[n];
        target[n - 1] = true;
        var walker = new KernelWalker(chain);
        const int passages = 40_000;
        double[] t = new double[passages];
        for (int i = 0; i < passages; i++) t[i] = walker.Passage(0, target, rng);
        double mean = t.Average(), se = Math.Sqrt(t.Sum(x => (x - mean) * (x - mean)) / (passages - 1) / passages);
        Assert.True(Math.Abs(mean - h) <= 4.0 * se, $"The walker's mean {mean} ± {se} against {h}.");
    }

    [Fact]
    public void Series_ReadsThePathProfilesArithmetic()
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s);
        var fixture = new BassMixing.Fixture("grid", data, config, null);
        var surface = new PathSurface(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        BassKernelEnumeration.State[] top = [.. space.States.Where(s => s.Terms.Length == 3).OrderByDescending(s => s.LogPi).Take(40)];
        BassTerm[] a = top[0].Terms, b = top.First(s => !s.Terms.Intersect(a).Any() && fixture.Basin(s.Terms) != fixture.Basin(a)).Terms;
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        ReversibleChain chain = ReversibleChain.Of([.. space.States.Select(s => s.LogPi)], flows.Values);
        int ia = space.Find(a), ib = space.Find(b);
        Assert.Single(chain.Classes());
        var reduction = new HittingReduction(chain, [ia, ib]);
        foreach (bool exact in new[] { true, false })
        {
            PathRoute route = F.Build(surface, a, b, exact, PathReading.Held, [])!;
            foreach (bool forward in new[] { true, false })
            {
                PathDirection read = F.Read(route, forward);
                RouteSeries whole = BassWaitingTime.Series(BassWaitingTime.Directed(route, forward));
                Assert.Equal(read.LnWaitExact, whole.LnWait, 9);
                Assert.Equal(read.LnWaitClass, whole.LnWaitClass, 9);
                Assert.Equal(read.DepthFrom, whole.Depth, 9);
                int[] end = BassWaitingTime.Code(forward ? b : a);
                PathHop[] kept = BassWaitingTime.CutAt(route, forward, end)!;
                RouteSeries cut = BassWaitingTime.Series(kept);
                Assert.True(cut.LnWait <= whole.LnWait + 1e-12);
                Assert.True(cut.Hops <= whole.Hops);
                int target = forward ? ib : ia, start = forward ? ia : ib;
                HittingReading r = reduction.Toward([target], [start]).Read(start)!.Value.Reading;
                double lnWait = whole.LnWait, identity = lnWait - r.LnHit - (lnWait - r.LnWaitEffective - r.LnMassRatio);
                Assert.True(Math.Abs(identity) <= 1e-9, $"The two corrections do not sum to the error: {identity}.");
                Assert.True(r.LnWaitEffective <= lnWait + 1e-9, $"Rayleigh: the whole network reads {r.LnWaitEffective}, the route {lnWait}.");
                double? lattice = BassWaitingTime.LatticeWait(surface, route.Lattice, forward ? route.Masks[0] : route.Masks[^1], [forward ? route.Masks[^1] : route.Masks[0]], whole.PiFrom, whole.LnSlideResistance);
                Assert.NotNull(lattice);
                Assert.True(r.LnWaitEffective <= lattice!.Value + 1e-9 && lattice.Value <= lnWait + 1e-9, $"Rayleigh: {r.LnWaitEffective} ≤ {lattice} ≤ {lnWait}.");
            }
        }
    }

    [Fact]
    public void Tally_ReadsOnePassage()
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s);
        var fixture = new BassMixing.Fixture("grid", data, config, null);
        var surface = new PathSurface(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        BassKernelEnumeration.State[] top = [.. space.States.Where(s => s.Terms.Length == 3).OrderByDescending(s => s.LogPi).Take(40)];
        BassTerm[] a = top[0].Terms, b = top.First(s => fixture.Basin(s.Terms) != fixture.Basin(a)).Terms;
        Assert.Equal(fixture.Basin(a) == fixture.Basin(b), BassWaitingTime.Code(a).SequenceEqual(BassWaitingTime.Code(b)));
        var spec = new PassageSpec
        {
            Part = "test", Start = "a", Arm = "H", Rung = "-", Chain = 0, B = 1.0, Surface = surface, Model = surface.Held, Terms = a, LogTau = Math.Log(surface.HeldTau),
            TargetExact = b, Basin = fixture.Basin, StayWindow = 50, Trace = true,
        };
        var tally = new PassageTally();
        var traced = new List<string>();
        tally.Advance(spec, 5_000_000, e => traced.Add(e.Basin));
        Assert.NotNull(tally.FirstEntry);
        Assert.True(tally.Crossed(spec));
        Assert.Equal(tally.FirstEntry!.Value + 50, tally.Sweeps);
        Assert.Equal(tally.FirstEntry.Value - 1, tally.PreSweeps);
        Assert.True(tally.InStartPre <= tally.PreSweeps && tally.TopDwell <= tally.TopSinceBasin && tally.TopSinceBasin <= tally.TopSinceAll);
        Assert.True(tally.StaySweeps <= tally.StaySeen);
        if (traced.Count > 0) Assert.Equal(fixture.Basin(b), traced[^1]);

        // The row reads back, and a resumed chain carries on from it.
        PassageTally back = PassageTally.Parse(tally.Row(spec));
        Assert.Equal(tally.FirstEntry, back.FirstEntry);
        Assert.Equal(tally.Sweeps, back.Sweeps);
        Assert.Equal(tally.BestPi, back.BestPi, 5);
        Assert.Equal(tally.Structure, BassMixing.Key(BassWaitingTime.ParseStructure(tally.Structure)));
        Assert.True(tally.BestPi >= surface.Score(a).Pi(PathReading.Held) - 1e-9);
    }

    // ---- The parts' shared fixtures ----

    private sealed record S10(BassMixing.Fixture Fixture, PathSurface Surface, Dictionary<(string Label, string Rank), BassPathProfileTests.Exemplar> Modes)
    {
        public BassPathProfileTests.Exemplar A1 => Modes[("a1", "1")];
    }

    private static S10 Load()
    {
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        return new S10(fixture, new PathSurface(fixture.Config, fixture.Data), BassPathProfileTests.ReadModes());
    }

    private static PassageSpec Chain(S10 s, string part, string label, string arm, string rung, int chain, double b, bool trace)
    {
        BassPathProfileTests.Exemplar e = s.Modes[(label, "1")];
        return new PassageSpec
        {
            Part = part, Start = label, Arm = arm, Rung = rung, Chain = chain, B = b, Surface = s.Surface, Model = arm == "H" ? s.Surface.Held : s.Surface.Model, Terms = e.Terms,
            LogTau = arm == "H" ? Math.Log(s.Surface.HeldTau) : s.Surface.Score(e.Terms).BestLogTau, TargetExact = s.A1.Terms, Basin = s.Fixture.Basin, Trace = trace,
        };
    }

    // ---- Part 1 ----

    [Fact]
    public void Part1_ExactHittingTimes()
    {
        if (!BassQualification.Enabled) return;   // heavy: M ≤ 4 holds about 12,000 states
        Quietly();
        PathCsv csv = Csv("exact", ExactHeader, 5);
        var pairs = new PathCsv("groundtruth", ["fixture"], 5).Rows().Where(r => r[0] == "grid").Select(r => (Cap: int.Parse(r[1], CultureInfo.InvariantCulture), Seed: int.Parse(r[2], CultureInfo.InvariantCulture), A: int.Parse(r[3], CultureInfo.InvariantCulture), B: int.Parse(r[4], CultureInfo.InvariantCulture), PiA: D(r[7]))).ToArray();
        foreach (var fixture in pairs.GroupBy(p => (p.Cap, p.Seed)).OrderBy(g => g.Key.Cap).ThenBy(g => g.Key.Seed))
            ExactFixture(csv, fixture.Key.Cap, fixture.Key.Seed, [.. fixture.Select(p => (p.A, p.B, p.PiA))]);
    }

    private void ExactFixture(PathCsv csv, int cap, int seed, (int A, int B, double PiA)[] pairs)
    {
        string[] Key(int from, int to, string target) => [F.I(cap), F.I(seed), F.I(from), F.I(to), target];
        if (pairs.All(p => new[] { "structure", "basin" }.All(t => csv.Has(Key(p.A, p.B, t)) && csv.Has(Key(p.B, p.A, t))))) return;
        var clock = Stopwatch.StartNew();
        BassConfig config = BassLocalStepTests.Continuous();
        config = config with { Structure = config.Structure with { MaxBasis = cap } };
        BassData data = BassSweepTests.Data(false, seed);
        var surface = new PathSurface(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        ReversibleChain chain = ReversibleChain.Of([.. space.States.Select(s => s.LogPi)], flows.Values);
        int[][] classes = chain.Classes();
        if (classes.Length != 1) throw new InvalidOperationException($"cap {cap} seed {seed}: {classes.Length} communicating classes, not read here.");
        foreach ((int a, int _, double piA) in pairs) Assert.True(Math.Abs(space.States[a].LogPi - piA) <= 2e-6, "The stored pair's state moved.");

        int[] modes = [.. pairs.SelectMany(p => new[] { p.A, p.B }).Distinct()];
        int[][] codes = [.. space.States.Select(s => BassWaitingTime.Code(s.Terms))];
        Dictionary<int, int[]> basins = modes.ToDictionary(m => m, m => Enumerable.Range(0, space.States.Count).Where(i => codes[i].SequenceEqual(codes[m])).ToArray());
        var reduction = new HittingReduction(chain, modes.Concat(basins.Values.SelectMany(b => b)).Distinct(), PassageBook.Threads);
        _out.WriteLine($"cap {cap} seed {seed}: {space.States.Count} states, {modes.Length} modes, kept {modes.Concat(basins.Values.SelectMany(b => b)).Distinct().Count()}; {clock.Elapsed.TotalSeconds:F0} s to the first stage.");

        var routes = new Dictionary<(int, int), (PathRoute? Exact, PathRoute? Aligned)>();
        foreach ((int a, int b, double _) in pairs)
            routes[(a, b)] = (F.Build(surface, space.States[a].Terms, space.States[b].Terms, true, PathReading.Held, []), F.Build(surface, space.States[a].Terms, space.States[b].Terms, false, PathReading.Held, []));

        foreach (int to in modes)
            foreach (string kind in new[] { "structure", "basin" })
            {
                int[] target = kind == "basin" ? basins[to] : [to];
                int[] froms = [.. pairs.Where(p => p.A == to || p.B == to).Select(p => p.A == to ? p.B : p.A)];
                int[] open = [.. froms.Where(f => !csv.Has(Key(f, to, kind)))];
                if (open.Length == 0) continue;
                int[] starts = [.. open.Where(f => !target.Contains(f))];
                foreach (int f in open.Except(starts))
                {
                    var row = new string[ExactHeader.Length];
                    Array.Fill(row, "");
                    Key(f, to, kind).CopyTo(row, 0);
                    row[^1] = "start in the target's basin";
                    csv.Append(row);
                }
                if (starts.Length == 0) continue;
                HittingReduction.Target stage = reduction.Toward(target, starts);
                var read = new Dictionary<int, (HittingReading R, double[] V)>();
                foreach (int f in starts)
                    if (stage.Read(f) is { } one) read[f] = one;
                double[]? dense = (cap <= 3 || kind == "structure") && read.Values.Any(r => r.R.LnHit <= 20.0) ? HittingReduction.Dense(chain, target, PassageBook.Threads) : null;
                foreach (int f in starts)
                {
                    if (!read.TryGetValue(f, out (HittingReading R, double[] V) r)) throw new InvalidOperationException($"cap {cap} seed {seed}: {f} cannot reach {to}.");
                    double denseH = dense is null || r.R.LnHit > 20.0 ? double.NaN : Math.Log(dense[f]);
                    if (!double.IsNaN(denseH) && Math.Abs(denseH - r.R.LnHit) > 5e-6) throw new InvalidOperationException($"cap {cap} seed {seed} {f} → {to} {kind}: the dense solve reads {denseH}, the reduction {r.R.LnHit}.");
                    if (Math.Abs(r.R.LnHit - r.R.LnHitFromVoltage) > 5e-6) throw new InvalidOperationException($"cap {cap} seed {seed} {f} → {to} {kind}: R_eff · M reads {r.R.LnHitFromVoltage}, h {r.R.LnHit}.");
                    (PathRoute? exact, PathRoute? aligned) = routes.TryGetValue((f, to), out var pair) ? pair : routes[(to, f)];
                    bool forward = routes.ContainsKey((f, to));
                    double massStart = 0.0, weighted = 0.0;
                    foreach (int x in basins[f])
                    {
                        double w = Math.Exp(chain.LogPi[x] - chain.LogPi[f]);
                        massStart += w;
                        weighted += w * r.V[x];
                    }
                    int[] targetCode = kind == "basin" ? codes[to] : [];
                    string[] Route(PathRoute? route, bool slides)
                    {
                        if (route is null || route.Hops.Length == 0) return slides ? ["", "", "", "", "", "", "", ""] : ["", "", "", "", "", ""];
                        PathHop[] hops = kind == "basin" ? BassWaitingTime.CutAt(route, forward, targetCode)! : BassWaitingTime.Directed(route, forward);
                        RouteSeries series = BassWaitingTime.Series(hops);
                        int fromMask = forward ? route.Masks[0] : route.Masks[^1];
                        int[] ends = kind == "basin" ? BassWaitingTime.MasksIn(route.Lattice, targetCode) : [forward ? route.Masks[^1] : route.Masks[0]];
                        double lattice = BassWaitingTime.LatticeWait(surface, route.Lattice, fromMask, ends, series.PiFrom, series.LnSlideResistance) ?? throw new InvalidOperationException("The route's own lattice joins no ends.");
                        const double slack = 1e-8;
                        if (r.R.LnWaitEffective > lattice + slack || lattice > series.LnWait + slack)
                            throw new InvalidOperationException($"Rayleigh fails: cap {cap} seed {seed} {f} → {to} {kind}: R_eff reads {r.R.LnWaitEffective}, the lattice {lattice}, the route {series.LnWait}.");
                        return slides
                            ? [F.F(series.LnWait), F.F(series.LnWaitClass), F.F(series.Depth), F.I(series.HighestM), F.I(series.Hops), F.I(series.SlideHops), F.F(lattice), route.Position]
                            : [F.F(series.LnWait), F.F(series.LnWaitClass), F.F(series.Depth), F.I(series.HighestM), F.I(series.Hops), F.F(lattice)];
                    }
                    csv.Append([.. Key(f, to, kind), F.I(space.States[f].Terms.Length), F.I(space.States[to].Terms.Length), F.I(space.States.Count), F.I(target.Length), F.F(r.R.LnHit), F.F(r.R.LnWaitEffective),
                                F.F(r.R.LnMassRatio), F.F(r.R.LnExcursion), F.F(Math.Log(massStart)), F.F(weighted / massStart), F.F(denseH), .. Route(exact, false), .. Route(aligned, true), ""]);
                }
            }
        _out.WriteLine($"cap {cap} seed {seed}: {clock.Elapsed.TotalSeconds:F0} s.");
    }

    // ---- Part 2 ----

    [Fact]
    public void Part2_LatticeResistance()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        PathCsv csv = Csv("lattice", LatticeHeader, 6);
        var pairs = new List<(string Pair, string Rank, BassTerm[] S, BassTerm[] T)>();
        foreach (string label in Shallow.Concat(Deep)) pairs.Add((PairName[label], "1/1", s.Modes[(label, "1")].Terms, s.A1.Terms));
        pairs.Add(("A-a1", "2/2", s.Modes[("A", "2")].Terms, s.Modes[("a1", "2")].Terms));
        pairs.Add(("A-a1", "d/1", s.Modes[("A", "d")].Terms, s.A1.Terms));
        pairs.Add(("A-a1", "1/d", s.Modes[("A", "1")].Terms, s.Modes[("a1", "d")].Terms));
        pairs.Add(("A-a1", "d/d", s.Modes[("A", "d")].Terms, s.Modes[("a1", "d")].Terms));
        var stored = new PathCsv("routes", ["n"], 7);
        foreach ((string pair, string rank, BassTerm[] a, BassTerm[] b) in pairs)
        {
            if (csv.Has("500", "S10", pair, rank, "aligned", "from S") && csv.Has("500", "S10", pair, rank, "aligned", "from T")) continue;
            var clock = Stopwatch.StartNew();
            PathRoute route = F.Build(s.Surface, a, b, false, PathReading.Held, [], priced: false)!;
            Assert.True(route.Lattice.Dim <= 11, $"{pair} {rank}: |D| = {route.Lattice.Dim}.");
            foreach (bool forward in new[] { true, false })
            {
                string direction = forward ? "from S" : "from T";
                double storedWait = D(stored.Rows().First(r => r[0] == "500" && r[1] == "S10" && r[2] == pair && r[3] == rank && r[4] == "aligned" && r[5] == nameof(PathReading.Held) && r[6] == direction)[29]);
                RouteSeries whole = BassWaitingTime.Series(BassWaitingTime.Directed(route, forward));
                if (Math.Abs(whole.LnWait - storedWait) > 1e-4) throw new InvalidOperationException($"{pair} {rank} {direction}: the rebuilt route reads {whole.LnWait}, the stored {storedWait}.");
                int from = forward ? route.Masks[0] : route.Masks[^1], to = forward ? route.Masks[^1] : route.Masks[0];
                double lattice = BassWaitingTime.LatticeWait(s.Surface, route.Lattice, from, [to], whole.PiFrom, whole.LnSlideResistance, PassageBook.Threads)!.Value;
                int[] code = BassWaitingTime.Code(forward ? b : a);
                PathHop[] kept = BassWaitingTime.CutAt(route, forward, code)!;
                RouteSeries cut = BassWaitingTime.Series(kept);
                double latticeCut = BassWaitingTime.LatticeWait(s.Surface, route.Lattice, from, BassWaitingTime.MasksIn(route.Lattice, code), cut.PiFrom, cut.LnSlideResistance, PassageBook.Threads)!.Value;
                if (lattice > whole.LnWait + 1e-8 || latticeCut > cut.LnWait + 1e-8) throw new InvalidOperationException($"Rayleigh fails: {pair} {rank} {direction}: the lattice reads {lattice}, the route {whole.LnWait}; cut, {latticeCut} and {cut.LnWait}.");
                int admitted = route.Lattice.Held.Count(p => !double.IsNegativeInfinity(p));
                csv.Append("500", "S10", pair, rank, "aligned", direction, F.I(route.Lattice.Dim), route.Position, F.I(whole.SlideHops), F.I(admitted), F.F(whole.PiFrom), F.F(whole.LnWait), F.F(lattice), F.F(whole.LnWait - lattice),
                           F.F(cut.LnWait), F.F(latticeCut), F.F(cut.LnWait - latticeCut), F.F(storedWait), clock.Elapsed.TotalSeconds.ToString("F1", CultureInfo.InvariantCulture));
                _out.WriteLine($"{pair} {rank} {direction}: |D| {route.Lattice.Dim}, ln W route {whole.LnWait:F2}, lattice {lattice:F2}; cut {cut.LnWait:F2}, lattice {latticeCut:F2}.");
            }
        }
    }

    // ---- The readings the chains are read against, from the stored steps ----

    private static (string Part, string Start, string Target, string N, string Data, string Pair, string Exemplar, string Alignment, bool Forward, BassTerm[] From, BassTerm[] To)[] ReadingPlan(S10 s)
    {
        var plan = new List<(string, string, string, string, string, string, string, string, bool, BassTerm[], BassTerm[])>();
        foreach (string label in Shallow.Concat(Deep))
            foreach (string alignment in new[] { "aligned", "partial", "exact" })
                plan.Add((Shallow.Contains(label) ? "3s" : "3d", label, "a1", "500", "S10", PairName[label], "1/1", alignment, true, s.Modes[(label, "1")].Terms, s.A1.Terms));
        var carried = new PathCsv("carried", ["n"], 4).Rows().ToArray();
        BassTerm[] C(string data, string label) => BassWaitingTime.ParseStructure(carried.Single(r => r[0] == "125" && r[1] == data && r[2] == label && r[3] == "1")[21]);
        plan.Add(("5", "a1", "A", "125", "stream", "A-a1", "1", "aligned", false, C("stream", "a1"), C("stream", "A")));
        plan.Add(("5", "B", "a1", "125", "fresh", "B-a1", "1", "partial+redraw", true, C("fresh", "B"), C("fresh", "a1")));
        return [.. plan];
    }

    [Fact]
    public void Readings_FromTheStoredSteps()
    {
        if (!BassQualification.Enabled) return;
        S10 s = Load();
        PathCsv csv = Csv("readings", ReadingHeader, 4);
        var stored = new PathCsv("routes", ["n"], 7).Rows().ToArray();
        foreach (var r in ReadingPlan(s))
        {
            StepRow[] rows = StepSeries.Read(r.N, r.Data, r.Pair, r.Exemplar, r.Alignment);
            Assert.NotEmpty(rows);
            RouteSeries full = StepSeries.Walk(rows, r.Forward, r.From, null), cut = StepSeries.Walk(rows, r.Forward, r.From, BassWaitingTime.Code(r.To));
            string[] row = stored.Single(x => x[0] == r.N && x[1] == r.Data && x[2] == r.Pair && x[3] == r.Exemplar && x[4] == r.Alignment && x[5] == nameof(PathReading.Held) && x[6] == (r.Forward ? "from S" : "from T"));
            double storedWait = D(row[29]), storedClass = D(row[30]);
            if (Math.Abs(full.LnWait - storedWait) > 2e-4 || Math.Abs(full.LnWaitClass - storedClass) > 2e-4)
                throw new InvalidOperationException($"{r.Part} {r.Start} {r.Alignment}: the rows' series reads {full.LnWait} and {full.LnWaitClass}, the stored {storedWait} and {storedClass}.");
            csv.Append(r.Part, r.Start, r.Target, r.Alignment, F.F(full.PiFrom), F.F(full.LnWait), F.F(full.LnWaitClass), F.F(cut.LnWait), F.F(cut.LnWaitClass), F.I(cut.Hops), F.I(cut.SlideHops), F.F(cut.Low), F.F(cut.Depth), F.I(cut.HighestM), F.F(storedWait), F.F(storedClass));
            _out.WriteLine($"{r.Part} {r.Start} → {r.Target} {r.Alignment}: full {full.LnWait:F2} (class {full.LnWaitClass:F2}); cut {cut.LnWait:F2} (class {cut.LnWaitClass:F2}), {cut.Hops} hops, {cut.SlideHops} slides, depth {cut.Depth:F2}.");
        }
    }

    // ---- The checks ----

    private static (double Mean, double Low, double High, double Cv, double Max) Interval(double[] t)
    {
        double mean = t.Average(), sd = Math.Sqrt(t.Sum(x => (x - mean) * (x - mean)) / (t.Length - 1)), half = 2.5758 * sd / Math.Sqrt(t.Length);
        return (mean, mean - half, mean + half, sd / mean, t.Max());
    }

    [Fact]
    public void Checks_KernelAndProductionChain()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        PathCsv validation = Validation();
        const int passages = 1000;
        string[][] rows = [.. Csv("exact", ExactHeader, 5).Rows().Where(r => r[9].Length > 0 && r[0] == "3")];
        // At least 3 pairs in the brief's range: the rows nearest ln 7, 9 and 11 to the exact structure, one to a basin, and one at b = ½ for the ladder.
        var chosen = new List<(string[] Row, double B)>();
        foreach (double wanted in new[] { 7.0, 9.0, 11.0 })
            chosen.Add((rows.Where(r => r[4] == "structure" && D(r[9]) is >= 6.0 and <= 12.0).OrderBy(r => Math.Abs(D(r[9]) - wanted)).First(), 1.0));
        chosen.Add((rows.Where(r => r[4] == "basin" && D(r[9]) is >= 6.0 and <= 12.0).OrderBy(r => Math.Abs(D(r[9]) - 8.0)).First(), 1.0));
        chosen.Add((chosen[1].Row, 0.5));
        chosen.Add((chosen[3].Row, 0.5));
        int index = 0;
        foreach ((string[] row, double b) in chosen)
        {
            int seed = int.Parse(row[1], CultureInfo.InvariantCulture), from = int.Parse(row[2], CultureInfo.InvariantCulture), to = int.Parse(row[3], CultureInfo.InvariantCulture);
            string kind = row[4], name = $"seed {seed} {from} to {to} {kind} b {b:F2}";
            index++;
            if (validation.Has("production chain", name)) continue;
            var clock = Stopwatch.StartNew();
            BassConfig config = BassLocalStepTests.Continuous();
            config = config with { Structure = config.Structure with { MaxBasis = 3 } };
            BassData data = BassSweepTests.Data(false, seed);
            var fixture = new BassMixing.Fixture("grid", data, config, null);
            var surface = new PathSurface(config, data);
            BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model, b);
            var flows = BassKernelEnumeration.Kernel(space, (key, st) => BassKernelEnumeration.D01(space, key, st, null));
            ReversibleChain chain = ReversibleChain.Of([.. space.States.Select(st => st.LogPi)], flows.Values);
            int[] code = BassWaitingTime.Code(space.States[to].Terms);
            int[] target = kind == "basin" ? [.. Enumerable.Range(0, space.States.Count).Where(i => BassWaitingTime.Code(space.States[i].Terms).SequenceEqual(code))] : [to];
            double lnH = new HittingReduction(chain, [from, .. target]).Toward(target, [from]).Read(from)!.Value.Reading.LnHit;
            if (lnH > 13.0) throw new InvalidOperationException($"{name}: ln h = {lnH:F1} is past the check's reach.");
            var inside = new bool[chain.Size];
            foreach (int t in target) inside[t] = true;
            var walker = new KernelWalker(chain);
            int[] seeds = SeedTree.Derive(9301 + index, passages);
            var kernel = new double[passages];
            var production = new double[passages];
            var options = new ParallelOptions { MaxDegreeOfParallelism = PassageBook.Threads };
            Parallel.For(0, passages, options, i => kernel[i] = walker.Passage(from, inside, new Xoshiro256PlusPlus(seeds[i])));
            Parallel.For(0, passages, options, i =>
            {
                var spec = new PassageSpec
                {
                    Part = "check", Start = name, Arm = "H", Rung = "-", Chain = i, B = b, Surface = surface, Model = surface.Held, Terms = space.States[from].Terms, LogTau = Math.Log(surface.HeldTau),
                    TargetExact = space.States[to].Terms, Basin = fixture.Basin, ToExact = kind == "structure", StayWindow = 0,
                };
                var tally = new PassageTally();
                tally.Advance(spec, long.MaxValue, null);
                production[i] = tally.FirstEntry!.Value;
            });
            double h = Math.Exp(lnH);
            var k = Interval(kernel);
            var p = Interval(production);
            bool agrees = h >= k.Low && h <= k.High, chainAgrees = h >= p.Low && h <= p.High;
            validation.Append("kernel walk", name, F.I(passages), F.F(Math.Log(k.Mean) - lnH), F.F(lnH), F.F(Math.Log(k.Mean)), F.F(Math.Log(k.Low)), F.F(Math.Log(k.High)), $"{(agrees ? "agrees" : "OUTSIDE")}; cv {k.Cv:F3}");
            validation.Append("production chain", name, F.I(passages), F.F(Math.Log(p.Mean) - lnH), F.F(lnH), F.F(Math.Log(p.Mean)), F.F(Math.Log(p.Low)), F.F(Math.Log(p.High)), $"{(chainAgrees ? "agrees" : "OUTSIDE")}; cv {p.Cv:F3}");
            _out.WriteLine($"{name}: ln h {lnH:F3}; kernel {Math.Log(k.Mean):F3} [{Math.Log(k.Low):F3}, {Math.Log(k.High):F3}], production {Math.Log(p.Mean):F3} [{Math.Log(p.Low):F3}, {Math.Log(p.High):F3}]; {clock.Elapsed.TotalSeconds:F0} s.");
            Assert.True(agrees && chainAgrees, $"{name}: h = {h:F1}, kernel [{k.Low:F1}, {k.High:F1}], production [{p.Low:F1}, {p.High:F1}].");
        }
    }

    // ---- The starts, against the stored π ----

    [Fact]
    public void Starts_RecheckedOnAFreshDesign()
    {
        if (!BassQualification.Enabled) return;
        S10 s = Load();
        double worst = 0.0;
        int count = 0;
        foreach (string label in Shallow.Concat(Deep).Append("a1"))
        {
            BassPathProfileTests.Exemplar e = s.Modes[(label, "1")];
            worst = Math.Max(worst, Math.Abs(s.Surface.Score(e.Terms).Pi(PathReading.Held) - e.PiHeld));
            Assert.Equal(e.Basin, s.Fixture.Basin(e.Terms));
            count++;
        }
        foreach (var part5 in Part5Plan())
        {
            worst = Math.Max(worst, Math.Abs(part5.Surface.Score(part5.From).Pi(PathReading.Held) - part5.PiFrom));
            worst = Math.Max(worst, Math.Abs(part5.Surface.Score(part5.To).Pi(PathReading.Held) - part5.PiTo));
            count += 2;
        }
        Validation().Append("starts", $"fresh design {Stamp}", F.I(count), F.F(worst), "", "", "", "", "each chain's start and target at the held τ against the stored π");
        _out.WriteLine($"{count} structures: worst {worst:E2}.");
        Assert.True(worst <= 2e-6);
    }

    [Fact]
    public void Rate_SweepsPerSecond()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        var part5 = Part5Plan()[1];
        PassageSpec[] specs =
        [
            Chain(s, "rate", "A", "H", "-", 0, 1.0, false), Chain(s, "rate", "a2", "H", "-", 1, 1.0, false), Chain(s, "rate", "A", "S", "-", 2, 1.0, false),
            Spec5(part5, "rate", 3),
        ];
        const long sweeps = 400_000;
        PathCsv validation = Validation();
        Parallel.ForEach(specs, new ParallelOptions { MaxDegreeOfParallelism = PassageBook.Threads }, spec =>
        {
            var tally = new PassageTally();
            tally.Advance(spec, sweeps, null);
            double rate = tally.Sweeps / tally.Seconds;
            validation.Append("rate", $"{spec.Start} arm {spec.Arm} n {spec.Surface.Observations} {Stamp}", F.I(tally.Sweeps), F.F(rate), F.F(tally.Accepted / (double)tally.Sweeps), F.F(tally.Seconds), "", "", "sweeps a second under four concurrent chains at below-normal priority; the acceptance beside it");
            _out.WriteLine($"{spec.Start} arm {spec.Arm} N = {spec.Surface.Observations}: {rate:N0} sweeps a second, acceptance {tally.Accepted / (double)tally.Sweeps:P2}.");
        });
    }

    // ---- Part 3 ----

    [Fact]
    public void Part3_ShallowStarts()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        var book = new PassageBook();
        PassageSpec[] specs = [.. Shallow.SelectMany(l => Enumerable.Range(0, 8).Select(c => Chain(s, "3s", l, "H", "-", c, 1.0, true)))];
        foreach (long limit in new[] { Million, 10 * Million, 20 * Million, 30 * Million, 40 * Million, 50 * Million, 60 * Million, 70 * Million, 80 * Million, 90 * Million, 100 * Million })
            if (!book.Advance(specs, limit, _out)) break;
    }

    [Fact]
    public void Part3_DeepStartsToTenMillion()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        var book = new PassageBook();
        foreach (long limit in new[] { Million, 10 * Million })
        {
            PassageSpec[] specs = [.. new[] { "H", "S" }.SelectMany(arm => Deep.SelectMany(l => Enumerable.Range(0, 4).Select(c => Chain(s, "3d", l, arm, "-", c, 1.0, true))))];
            if (!book.Advance(specs, limit, _out)) break;
        }
    }

    [Fact]
    public void Part3_DeepStartsToHundredMillion()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        var book = new PassageBook();
        PassageSpec[] specs = [.. Deep.SelectMany(l => Enumerable.Range(0, 4).Select(c => Chain(s, "3d", l, "H", "-", c, 1.0, true)))];
        for (int tens = 2; tens <= 10; tens++)
            if (!book.Advance(specs, tens * 10 * Million, _out)) break;
    }

    // Added after registration: twelve more chains on each deep start, for the exponential check, the departure and the gateway at b = 1.
    [Fact]
    public void Part3_MoreDeepChains()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        var book = new PassageBook();
        PassageSpec[] specs = [.. Deep.SelectMany(l => Enumerable.Range(0, 12).Select(c => Chain(s, "3x", l, "H", "-", c, 1.0, false)))];
        foreach (long limit in new[] { Million, 10 * Million, 20 * Million, 30 * Million, 40 * Million, 50 * Million, 60 * Million, 70 * Million, 80 * Million, 90 * Million, 100 * Million })
            if (!book.Advance(specs, limit, _out)) break;
    }

    // ---- Part 4 ----

    [Fact]
    public void Part4_ReadingsAndRungs()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        PathCsv ladder = Csv("ladder", LadderHeader, 2), rungs = Csv("rungs", RungHeader, 2);
        var stored = new PathCsv("routes", ["n"], 7).Rows().ToArray();
        foreach (string label in Deep)
        {
            string pair = PairName[label];
            if (LadderNats.All(n => rungs.Has(pair, F.F(n)))) continue;
            var clock = Stopwatch.StartNew();
            var at = new LadderPair(s.Surface, s.Modes[(label, "1")].Terms, s.A1.Terms);
            var read = new List<LadderReading>();
            for (int k = 1; k <= 200; k++)
                if (at.Read(k * 0.005) is { } one)
                {
                    read.Add(one);
                    ladder.Append(pair, F.F(one.B), one.Position, F.F(one.PiFrom), F.F(one.LnWaitFull), F.F(one.LnWaitClassFull), F.F(one.LnWait), F.F(one.LnWaitClass), F.F(one.Low), F.F(one.LowFull), F.I(one.Hops), F.I(one.HighestM));
                }
            // b = 1 is the stored route's reading.
            LadderReading top = read[^1];
            string[] row = stored.Single(x => x[0] == "500" && x[1] == "S10" && x[2] == pair && x[3] == "1/1" && x[4] == "aligned" && x[5] == nameof(PathReading.Held) && x[6] == "from S");
            if (Math.Abs(top.LnWaitFull - D(row[29])) > 1e-4 || Math.Abs(top.LnWaitClassFull - D(row[30])) > 1e-4)
                throw new InvalidOperationException($"{pair}: b = 1 reads {top.LnWaitFull} and {top.LnWaitClassFull}, the stored route {row[29]} and {row[30]}.");
            var taken = new Dictionary<double, double>();
            foreach (double nat in LadderNats)
            {
                LadderReading best = read.OrderBy(r => Math.Abs(r.LnWaitClass - nat)).ThenByDescending(r => r.B).First();
                string duplicate = taken.TryGetValue(best.B, out double first) ? F.F(first) : "";
                taken.TryAdd(best.B, nat);
                rungs.Append(pair, F.F(nat), F.F(best.B), F.F(best.LnWait), F.F(best.LnWaitClass), F.F(best.LnWaitClass - nat), F.F(Math.Ceiling(10.0 * Math.Exp(best.LnWaitClass))), duplicate);
                _out.WriteLine($"{pair}, ln {nat:F0}: b = {best.B:F3}, class {best.LnWaitClass:F2}, route {best.LnWait:F2}, {best.Hops} hops, slide {best.Position}.");
            }
            _out.WriteLine($"{pair}: {clock.Elapsed.TotalSeconds:F0} s.");
        }
    }

    [Fact]
    public void Part4_Chains()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        var book = new PassageBook();
        var rungs = Csv("rungs", RungHeader, 2).Rows().ToArray();
        Assert.NotEmpty(rungs);
        foreach (double nat in LadderNats)
        {
            var specs = new List<(PassageSpec Spec, long Cap)>();
            foreach (string[] r in rungs.Where(r => Math.Abs(D(r[1]) - nat) < 1e-9 && r[7].Length == 0))
            {
                string label = PairName.Single(p => p.Value == r[0]).Key;
                for (int c = 0; c < 8; c++) specs.Add((Chain(s, "4", label, "H", $"ln{nat:F0}", c, D(r[2]), false), (long)D(r[6])));
            }
            long top = specs.Max(x => x.Cap);
            for (long limit = Math.Min(top, 10 * Million); ; limit = Math.Min(top, limit + 10 * Million))
            {
                foreach (var group in specs.GroupBy(x => Math.Min(x.Cap, limit)))
                    if (!book.Advance([.. group.Select(x => x.Spec)], group.Key, _out)) return;
                if (limit >= top) break;
            }
        }
    }

    // ---- Part 5 ----

    private sealed record Plan5(string Start, string TargetLabel, PathSurface Surface, BassMixing.Fixture Fixture, BassTerm[] From, BassTerm[] To, double PiFrom, double PiTo);

    private static Plan5[] Part5Plan()
    {
        var carried = new PathCsv("carried", ["n"], 4).Rows().ToArray();
        Plan5 Make(string data, int seed, string start, string target)
        {
            BassData on = BassQualificationRecoveryTests.Friedman(125, seed).Data;
            BassConfig config = BassQualificationRecoveryTests.Config(53);
            string[] from = carried.Single(r => r[0] == "125" && r[1] == data && r[2] == start && r[3] == "1"), to = carried.Single(r => r[0] == "125" && r[1] == data && r[2] == target && r[3] == "1");
            return new Plan5(start, target, new PathSurface(config, on), new BassMixing.Fixture("N = 125", on, config, null), BassWaitingTime.ParseStructure(from[21]), BassWaitingTime.ParseStructure(to[21]), D(from[7]), D(to[7]));
        }
        return [Make("stream", 51, "a1", "A"), Make("fresh", SeedTree.Derive(9101, 5)[0], "B", "a1")];
    }

    private static PassageSpec Spec5(Plan5 p, string part, int chain) => new()
    {
        Part = part, Start = p.Start, Arm = "H", Rung = "-", Chain = chain, B = 1.0, Surface = p.Surface, Model = p.Surface.Held, Terms = p.From, LogTau = Math.Log(p.Surface.HeldTau),
        TargetExact = p.To, Basin = p.Fixture.Basin, Trace = part == "5",
    };

    [Fact]
    public void Part5_ChainsAtN125()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        var book = new PassageBook();
        PassageSpec[] specs = [.. Part5Plan().SelectMany(p => Enumerable.Range(0, 8).Select(c => Spec5(p, "5", c)))];
        foreach (long limit in new[] { Million, 10 * Million, 20 * Million, 30 * Million, 40 * Million, 50 * Million, 60 * Million, 70 * Million, 80 * Million, 90 * Million, 100 * Million })
            if (!book.Advance(specs, limit, _out)) break;
    }

    // ---- The flat histogram's mass ----

    [Fact]
    public void FlatHistogram_TheBasinsMass()
    {
        if (!BassQualification.Enabled) return;   // heavy: arm (i) again, about half an hour
        Quietly();
        S10 s = Load();
        PathCsv csv = Csv("flathistogram", FlatHeader, 1);
        string[] labels = [.. s.Modes.Keys.Where(k => k.Rank == "1").Select(k => k.Label)];
        if (labels.All(l => csv.Has(l)) && csv.Has("largest 1")) return;
        BassConfig config = BassFlatHistogramQualificationTests.ArmConfig(true, null, BassFlatHistogramQualificationTests.ArmOneSeed);
        BassFit fit = Bass.Fit(config, s.Fixture.Data);
        var draws = fit.Diagnostics.FlatHistogram!.Draws!;
        double[] weights = fit.DrawWeights!;
        var byBasin = new Dictionary<string, (double W, long N)>();
        var byStructure = new Dictionary<string, (double W, long N)>();
        for (int c = 0; c < fit.Draws.Count; c++)
            for (int d = 0; d < draws.Draws; d++)
            {
                BassTerm[] terms = [.. fit.Draws[c][d].Structure.TermIds.Select(id => fit.Terms[id])];
                string basin = s.Fixture.Basin(terms), key = BassMixing.Key(terms);
                double w = weights[c * draws.Draws + d];
                (double W, long N) b = byBasin.GetValueOrDefault(basin), st = byStructure.GetValueOrDefault(key);
                byBasin[basin] = (b.W + w, b.N + 1);
                byStructure[key] = (st.W + w, st.N + 1);
            }
        string[][] modeRows = [.. new PathCsv("modes", ["label"], 4).Rows()];
        foreach (string label in labels)
        {
            BassPathProfileTests.Exemplar e = s.Modes[(label, "1")];
            (double W, long N) basin = byBasin.GetValueOrDefault(e.Basin), top = byStructure.GetValueOrDefault(BassMixing.Key(e.Terms));
            double stored = modeRows.Where(r => r[0] == label && r[1] == "1").Select(r => r[18].Length == 0 ? double.NaN : D(r[18])).FirstOrDefault(double.NaN);
            double ratio = top.N > 0 ? Math.Log(basin.W) - Math.Log(top.W) : double.NaN;
            csv.Append(label, e.Basin, F.F(basin.W), F.F(stored), F.F(ratio), F.I(top.N), F.I(basin.N), F.F(e.PiHeld), F.F(e.PiHeld + ratio), "draws of the structure beside it: a single draw holds a factor of e or more");
            _out.WriteLine($"{label}: share {basin.W:P3}, {basin.N} draws; the structure {top.N} draws; ln(mass / π(top)) {ratio:F2}.");
        }
        int place = 0;
        foreach ((string name, (double W, long N) held) in byBasin.OrderByDescending(b => b.Value.W).Take(12))
        {
            place++;
            csv.Append($"largest {place}", name, F.F(held.W), "", "", "", F.I(held.N), "", "", "the run's largest basins by share");
            _out.WriteLine($"{held.W:P3} | {held.N:N0} draws | {name}");
        }
    }
}
