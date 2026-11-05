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
/// The route series credited with a term's whole measure, and the wait to hold a1's basin (the basin-hold brief; gap G1, readings (f) and (g)).
/// <list type="bullet">
/// <item>Always run: the hold tally against a brute-force read of random streams, stopped and resumed, and on one enumerated fixture the class sums against the kernel's own flows and Rayleigh over the lumped kernel.</item>
/// <item>Qualification runs (<see cref="BassQualification.Enabled"/>): the parts, each its own test. Rows go beside the brief as <c>basin-hold-*.csv</c>, appended as reached and skipped on a rerun.</item>
/// </list>
/// Choices the brief leaves open:
/// <list type="bullet">
/// <item>Every class is summed whole. No route holds a three-factor term, so no class passes about 10⁶ members and none needs the brief's sample; the sampled sum's gate has nothing to read.</item>
/// <item>R_short shorts the target's whole basin for both targets, the structure's included, so on a structure target ln(R_eff / R_short) also holds finding the target's own knots. A structure target in the start's basin has no R_short, and its row says so.</item>
/// <item>R_lump reads the route's basin sequence as the route takes it, a basin it returns to included. The rows count those returns.</item>
/// <item>A term's credit in W_held is read at the edge's lower state, and an addition and a deletion take the same c_e (<see cref="BasinSeries"/>).</item>
/// <item>W_held,25 is the arithmetic behind the path profile's table: the class flow on an addition, ln of its class's size on every later step, a slide's included, whether or not the term has since gone.</item>
/// <item>"a1's side" of the two-well balance is every structure not holding both an x1 term and an x2 term, the set the share is conditioned on. The narrower sides are written beside it. The balance reads the first return past T_hold + W, and its exposure runs from there.</item>
/// <item>A visit ends at 1,000 sweeps outside the target, the departure's length; a touch is a visit closed before the hold starts.</item>
/// <item>The checks run on the kernel with windows of 50, 200 and 1,000, primary (1,000, ½), and visits ended at 20 sweeps outside, since these spaces mix in thousands of sweeps and a narrower window leaves no touch to count. Their intervals are 99%, as the waiting-time brief's checks; the chains' waits take the exact Poisson 95% the outcome rules name.</item>
/// </list>
/// </summary>
public sealed class BassBasinHoldTests
{
    private const long Million = 1_000_000;

    // The brief's 4 core-hours for Parts 1–3, less what Part 1 and the shapes take.
    private const double ClassBudgetHours = 3.9;

    private static readonly int[] CheckWindows = [50, 200, 1_000];

    private static readonly string[] Shallow = ["a2", "a1x"], Deep = ["A", "B", "chain0", "chain1", "chain2", "chain3"], Alignments = ["aligned", "partial", "exact"];
    private static readonly Dictionary<string, string> PairName = new()
    {
        ["A"] = "A-a1", ["B"] = "B-a1", ["a2"] = "a2-a1", ["a1x"] = "a1+x6x8-a1", ["chain0"] = "chain0-a1", ["chain1"] = "chain1-a1", ["chain2"] = "chain2-a1", ["chain3"] = "chain3-a1",
    };

    // The path profile's column "every addition's class credited", to one decimal, to a1's structure.
    private static readonly Dictionary<(string Start, string Alignment), double> CreditedTable = new()
    {
        [("chain0", "aligned")] = 35.4, [("chain0", "partial")] = 35.0, [("chain0", "exact")] = 25.5,
        [("chain1", "aligned")] = 18.8, [("chain1", "partial")] = 13.7, [("chain1", "exact")] = 11.0,
        [("chain2", "aligned")] = 23.5, [("chain2", "partial")] = 21.6, [("chain2", "exact")] = 13.9,
        [("chain3", "aligned")] = 34.7, [("chain3", "partial")] = 33.2, [("chain3", "exact")] = 22.8,
    };

    private static readonly string[] ExactRoute = ["ln_w", "ln_w_class", "depth", "highest_m", "basin_hops", "returns", "ln_w_lump", "ln_w_own", "ln_w_held", "ln_w_held25"];

    private static readonly string[] ExactHeader =
    [
        "cap", "seed", "from", "to", "target", "m_from", "m_to", "states", "dictionary", "basins", "ln_h", "ln_w_eff", "ln_w_short",
        .. ExactRoute.Select(c => $"exact_{c}"), .. ExactRoute.Select(c => $"aligned_{c}"), "note",
    ];

    private static readonly string[] ClassHeader = ["kind", "structure", "mask", "ln", "members", "admitted", "seconds"];

    private static readonly string[] HopHeader = ["start", "alignment", "index", "kind", "term", "m_lower", "from", "to", "ln_flow", "ln_own", "members", "admitted", "earlier", "ln_credit", "credits", "ln_w_edge"];

    private static readonly string[] SeriesHeader =
    [
        "start", "alignment", "pi_from", "edges", "slides", "ln_w_route", "ln_w_class", "ln_w_held25", "ln_w_own", "ln_w_held", "ln_w_route_full", "ln_w_class_full", "ln_w_held25_full", "seconds",
    ];

    private static readonly string[] ShapeHeader =
    [
        "start", "alignment", "sequence_basins", "gateways", "gateways_on", "gateways_last_door", "doors_on", "doors_off", "crossings", "moves", "moves_on",
    ];

    private static readonly string[] ValidationHeader = ["check", "case", "count", "worst", "value_a", "value_b", "value_c", "value_d", "note"];

    private static readonly string[] FlatHeader = ["kind", "label", "basin", "holds_both", "share", "conditional", "draws", "by_chain", "note"];

    private readonly ITestOutputHelper _out;
    public BassBasinHoldTests(ITestOutputHelper output) => _out = output;

    private static PathCsv Csv(string file, string[] header, int keys) => new(file, header, keys, "basin-hold");
    private static PathCsv Validation() => Csv("validation", ValidationHeader, 2);
    private static void Quietly() => Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
    private static double D(string v) => double.Parse(v, CultureInfo.InvariantCulture);
    private static string Stamp => DateTime.UtcNow.ToString("yyyy-MM-ddTHH:mm", CultureInfo.InvariantCulture);
    private static string Text(double v) => v.ToString("F1", CultureInfo.InvariantCulture);

    // The processor time this host has taken, for the budget's ledger.
    private static void Ledger(string name) => Validation().Append("ledger", $"{name} {Stamp}", "", F.F(Process.GetCurrentProcess().TotalProcessorTime.TotalHours), "", "", "", "", "core-hours this test host took");

    // ---- Always run ----

    // Random streams through the tally, against a read of the whole stream by brute force, once straight and once stopped and resumed from its row.
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Tally_AgreesWithABruteForceRead(int seed)
    {
        const int n = 40_000, gap = 20;
        static BassTerm T(params int[] v) => new([.. v.Select(x => new HingeFactor(x, 1, 0))], null, null);
        BassTerm[] both = [T(0), T(1)], middle = [T(2)], other = [T(3), T(2)], target = [T(0, 1), T(2)];
        BassTerm[][] states = [both, middle, other, target];
        double[] stay = [0.999, 0.99, 0.995, 0.995 + 0.004 * (seed % 2)];
        string Basin(IEnumerable<BassTerm> terms) => string.Join(" + ", terms.Select(t => string.Join("·", t.Continuous.Select(h => $"x{h.Variable + 1}"))).Order(StringComparer.Ordinal));
        var spec = new HoldSpec
        {
            Start = "s", Arm = "H", Chain = seed, Model = null!, Terms = both, LogTau = 0.0, Target = BassWaitingTime.Code(target), Basin = Basin, Largest = [BassWaitingTime.Code(other), BassWaitingTime.Code(target)],
            Windows = [10, 50, 200], Gap = gap, Block = 1_000,
        };
        var rng = new Xoshiro256PlusPlus(600 + seed);
        var stream = new int[n + 1];
        for (int t = 1; t <= n; t++) stream[t] = rng.NextDouble() < stay[stream[t - 1]] ? stream[t - 1] : (stream[t - 1] + 1 + (int)(3 * rng.NextDouble())) % 4;

        HoldTally Run(int? stopAt)
        {
            var tally = new HoldTally();
            tally.Begin(spec, states[0]);
            var visits = new List<string[]>();
            var blocks = new List<string[]>();
            for (int t = 1; t <= n; t++)
            {
                tally.Tick(stream[t] == stream[t - 1] ? null : states[stream[t]]);
                if (t != stopAt) continue;
                visits.AddRange(tally.PendingVisits);
                blocks.AddRange(tally.PendingBlocks);
                tally = HoldTally.Parse(tally.Row(spec));
                tally.Begin(spec, states[stream[t]]);
            }
            tally.PendingVisits.InsertRange(0, visits);
            tally.PendingBlocks.InsertRange(0, blocks);
            return tally;
        }

        HoldTally straight = Run(null);
        bool[] inTarget = [.. stream.Select(s => s == 3)];
        int first = Array.FindIndex(inTarget, 1, b => b);
        Assert.True(first > 0);
        Assert.Equal(first, straight.FirstEntry);
        Assert.Equal(Basin(states[stream[first - 1]]), straight.EntryFrom);
        Assert.Equal(Array.FindIndex(stream, 1, s => s != 0), straight.FirstExit);
        Assert.Equal(first - 1, straight.PreSweeps);
        Assert.Equal(first - 1, straight.PreFamily.Sum());
        long window = 0;
        for (int w = 0; w < 3; w++)
            for (int c = 0; c < 2; c++)
            {
                int width = spec.Windows[w];
                int need = (int)Math.Ceiling(spec.Levels[c] * width - 1e-9);
                long? expected = null;
                window = 0;
                for (int t = 1; t <= n && expected is null; t++)
                {
                    if (inTarget[t]) window++;
                    if (t > width && inTarget[t - width]) window--;
                    if (t >= width && window >= need) expected = t - width + 1;
                }
                Assert.Equal(expected, straight.Held[w * 2 + c]);
            }
        Assert.NotNull(straight.Hold);
        long hold = straight.Hold!.Value, post = hold + 200;
        long detected = post - 1;
        Assert.Equal(n - post + 1, straight.PostSweeps);
        Assert.Equal(n - post + 1, straight.PostFamily.Sum());
        long Count(Func<int, bool> test) => Enumerable.Range((int)post, (int)(n - post + 1)).Count(test);
        Assert.Equal(Count(t => stream[t] != 0), straight.PostNotBoth);
        Assert.Equal(Count(t => inTarget[t]), straight.PostTarget);
        Assert.Equal(Count(t => stream[t] == 2), straight.PostLargest[0]);
        Assert.Equal(Count(t => inTarget[t]), straight.PostLargest[1]);
        int returning = Array.FindIndex(stream, (int)post, s => s == 0);
        Assert.Equal(returning < 0 ? null : returning, straight.Return);

        // Entries and visits to the hold's detection: runs in the target, merged across at most a gap's sweeps outside.
        var runs = new List<(int First, int Last)>();
        int entries = 0;
        for (int t = 1; t <= detected; t++)
        {
            if (!inTarget[t]) continue;
            if (!inTarget[t - 1]) entries++;
            if (runs.Count > 0 && t - runs[^1].Last <= gap) runs[^1] = (runs[^1].First, t);
            else runs.Add((t, t));
        }
        Assert.Equal(entries, straight.Entries);
        Assert.Equal(runs.Count, straight.Visits);
        Assert.Equal(runs.Count(r => r.Last < hold), straight.Touches);
        Assert.Equal(runs.Count, straight.PendingVisits.Count);

        // Whole blocks hold the share's sweeps.
        Assert.Equal(straight.PostSweeps, straight.PendingBlocks.Sum(b => long.Parse(b[4], CultureInfo.InvariantCulture)));
        Assert.Equal(straight.PostTarget, straight.PendingBlocks.Sum(b => long.Parse(b[6], CultureInfo.InvariantCulture)));
        Assert.Equal(straight.PostNotBoth, straight.PendingBlocks.Sum(b => long.Parse(b[5], CultureInfo.InvariantCulture)));

        // A departure is the first stretch of more than a gap's sweeps outside the start's basin.
        int departure = -1;
        for (int t = 1; t <= n - gap && departure < 0; t++)
            if (stream[t] != 0 && stream[t - 1] == 0 && Enumerable.Range(t, gap + 1).All(u => stream[u] != 0)) departure = t;
        Assert.Equal(departure < 0 ? null : departure, straight.Departure);
        long region = departure < 0 || departure + gap >= first ? 0 : first - departure - gap;
        Assert.Equal(region, straight.RegionSweeps);
        Assert.Equal(region, straight.RegionBasins(8).Sum(c => c.Sweeps));

        // Stopped anywhere a block ends, the tally reads as if it had not, its rows the same but for the restart's stored stretches.
        foreach (int stop in new[] { (int)(first / 1000 + 1) * 1000, (int)(post / 1000) * 1000, (int)(post / 1000 + 2) * 1000 })
        {
            HoldTally resumed = Run(stop);
            string[] a = straight.Row(spec), b = resumed.Row(spec);
            for (int k = 0; k < a.Length; k++)
            {
                if (HoldTally.Header[k] is "region_basins") continue;
                Assert.True(a[k] == b[k], $"stopped at {stop}: {HoldTally.Header[k]} reads {b[k]}, against {a[k]}");
            }
            Assert.Equal(straight.RegionBasins(8), resumed.RegionBasins(8));
            Assert.Equal(straight.PendingVisits.Count, resumed.PendingVisits.Count);
            Assert.Equal(straight.PendingBlocks.Count, resumed.PendingBlocks.Count);
        }
    }

    [Fact]
    public void CodeOf_ReadsWhatBasinWrote()
    {
        foreach ((string basin, int[] code) in new (string, int[])[] { ("4 × x1·x2 + 2 × x3 + x4 + x5", [3, 3, 3, 3, 4, 4, 8, 16]), ("x10·x9 + x1", [1, 768]), ("∅", []) })
            Assert.Equal(code, BasinSeries.CodeOf(basin));
    }

    // The class sums against the kernel's flows, window 0 and ±25 against the route's own reading, and Rayleigh over the lumped kernel.
    [Fact]
    public void Classes_ReadTheKernel()
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s);
        var surface = new PathSurface(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        ReversibleChain chain = ReversibleChain.Of([.. space.States.Select(s => s.LogPi)], flows.Values);
        Dictionary<(int From, int To), double> births = flows[BassKernelEnumeration.BirthDeath];
        int[][] codes = [.. space.States.Select(s => BassWaitingTime.Code(s.Terms))];
        Lumped lumped = Lumped.Of(space, chain, codes);

        // A term's whole class from a state is the kernel's flows to the structures holding one more term on those variables.
        BassKernelEnumeration.State lower = space.States.First(s => s.Terms.Length == 2);
        int at = space.States.IndexOf(lower);
        foreach (BassTerm t in space.Terms.Where(u => u.Continuous.Count <= 2 && !lower.Terms.Contains(u)).DistinctBy(BassWaitingTime.Mask))
        {
            ClassRead own = ClassSums.Own(surface, lower.Terms, t, null);
            double kernel = double.NegativeInfinity;
            int count = 0;
            foreach (((int from, int to), double logK) in births)
            {
                if (from != at || space.States[to].Terms.Length != 3) continue;
                BassTerm added = space.States[to].Terms.Single(u => !lower.Terms.Contains(u));
                if (BassWaitingTime.Mask(added) != BassWaitingTime.Mask(t)) continue;
                kernel = LogMath.Add(kernel, lower.LogPi + logK);
                count++;
            }
            if (count == 0) continue;
            Assert.Equal(count, own.Admitted);
            Assert.True(Math.Abs(own.Ln - kernel) <= 1e-9, $"{BassMixing.Key([t])}: the sum reads {own.Ln}, the kernel {kernel}.");
        }

        BassKernelEnumeration.State[] top = [.. space.States.Where(s => s.Terms.Length == 3).OrderByDescending(s => s.LogPi).Take(40)];
        BassTerm[] a = top[0].Terms, b = top.First(s => !s.Terms.Intersect(a).Any() && !BassWaitingTime.Code(s.Terms).SequenceEqual(BassWaitingTime.Code(a))).Terms;
        int ia = space.Find(a), ib = space.Find(b);
        foreach (bool exact in new[] { true, false })
        {
            PathRoute route = F.Build(surface, a, b, exact, PathReading.Held, [])!;
            foreach (bool forward in new[] { true, false })
            {
                BasinRoute whole = BasinSeries.FromRoute(route, forward, null);
                RouteSeries series = BassWaitingTime.Series(BassWaitingTime.Directed(route, forward));
                Assert.Equal(series.LnWait, BasinSeries.Route(whole), 9);
                Assert.Equal(series.LnWait, BasinSeries.CreditedToLater(whole, e => (ClassSums.Own(surface, e.Lower, e.Term, 0).Ln, ClassSums.Own(surface, e.Lower, e.Term, 0).Members)), 9);
                Assert.Equal(series.LnWaitClass, BasinSeries.OnAdditions(whole, e => ClassSums.Own(surface, e.Lower, e.Term, F.ClassWindow).Ln), 9);
                Assert.Equal(series.LnWaitClass, BasinSeries.OnAdditions(whole, e => e.LnFlowClass), 9);
                foreach (BasinEdge e in whole.Edges)
                {
                    Assert.Equal(e.LnFlow, ClassSums.Own(surface, e.Lower, e.Term, 0).Ln, 9);
                    foreach (BassTerm held in e.Earlier)
                        Assert.Equal(0.0, ClassSums.Swap(surface, e.Lower, held, 0).Ln, 9);
                }
                int[] end = BassWaitingTime.Code(forward ? b : a);
                BasinRoute cut = BasinSeries.FromRoute(route, forward, end);
                Assert.Equal(BassWaitingTime.Series(BassWaitingTime.CutAt(route, forward, end)!).LnWait, BasinSeries.Route(cut), 9);
                double own = BasinSeries.Own(cut, e => ClassSums.Own(surface, e.Lower, e.Term, null).Ln);
                double whole2 = BasinSeries.Held(cut, e => ClassSums.Own(surface, e.Lower, e.Term, null).Ln, (e, t) => ClassSums.Swap(surface, e.Lower, t, null).Ln);
                double lump = BasinSeries.Lumped(cut, lumped.LnFlow);
                int from = forward ? ia : ib;
                double shorted = lumped.LnWaitShort(chain, from, end);
                HittingReading r = new HittingReduction(chain, [from, .. space.States.Select((s, i) => i).Where(i => codes[i].SequenceEqual(end))]).Toward([.. space.States.Select((s, i) => i).Where(i => codes[i].SequenceEqual(end))], [from]).Read(from)!.Value.Reading;
                Assert.True(shorted <= r.LnWaitEffective + 1e-9, $"Rayleigh: shorted {shorted}, R_eff {r.LnWaitEffective}.");
                Assert.True(r.LnWaitEffective <= BasinSeries.Route(cut) + 1e-9);
                Assert.True(shorted <= lump + 1e-9 && lump <= own + 1e-9 && own <= BasinSeries.Route(cut) + 1e-9, $"R_short {shorted}, R_lump {lump}, R_own {own}, route {BasinSeries.Route(cut)}.");
                Assert.True(whole2 <= own + 1e-9, $"W_held {whole2} above W_own {own}: a class holds its own member, so a credit is at least 1.");
            }
        }
    }

    // ---- The lumped kernel ----

    /// <summary>A kernel lumped by basin: each basin's π and the flows between basins, summed.</summary>
    private sealed class Lumped
    {
        private readonly Dictionary<string, int> _ids;
        private readonly Dictionary<(int, int), double> _flow;
        private readonly double _total;
        private readonly int[] _basinOf;
        private readonly Dictionary<(int, int), double> _escape = [];
        private readonly HittingReduction _reduction;

        private Lumped(Dictionary<string, int> ids, Dictionary<(int, int), double> flow, double total, int[] basinOf, ReversibleChain chain)
        {
            _ids = ids;
            _flow = flow;
            _total = total;
            _basinOf = basinOf;
            Chain = chain;
            _reduction = new HittingReduction(chain, Enumerable.Range(0, chain.Size), 1);
        }

        public ReversibleChain Chain { get; }

        public int Basins => Chain.Size;

        private static string K(int[] code) => string.Join(".", code);

        public static Lumped Of(BassKernelEnumeration.Space space, ReversibleChain chain, int[][] codes)
        {
            var ids = new Dictionary<string, int>();
            var basinOf = new int[codes.Length];
            var pi = new List<double>();
            for (int i = 0; i < codes.Length; i++)
            {
                string key = K(codes[i]);
                if (!ids.TryGetValue(key, out int id))
                {
                    ids[key] = id = pi.Count;
                    pi.Add(double.NegativeInfinity);
                }
                basinOf[i] = id;
                pi[id] = LogMath.Add(pi[id], chain.LogPi[i]);
            }
            var flow = new Dictionary<(int, int), double>();
            for (int i = 0; i < chain.Size; i++)
                for (int e = chain.Starts[i]; e < chain.Starts[i + 1]; e++)
                {
                    int a = basinOf[i], b = basinOf[chain.To[e]];
                    if (a == b) continue;
                    flow[(a, b)] = LogMath.Add(flow.GetValueOrDefault((a, b), double.NegativeInfinity), chain.LogPi[i] + chain.LogRate[e]);
                }
            foreach (((int a, int b), double c) in flow)
                if (Math.Abs(c - flow[(b, a)]) > 1e-8) throw new InvalidOperationException("The lumped flows are not balanced.");
            var rates = flow.ToDictionary(kv => kv.Key, kv => kv.Value - pi[kv.Key.Item1]);
            var lumped = ReversibleChain.Of([.. pi], [rates.ToDictionary(kv => (kv.Key.Item1, kv.Key.Item2), kv => kv.Value)]);
            return new Lumped(ids, flow, space.States[0].LogPi - chain.LogPi[0], basinOf, lumped);
        }

        /// <summary>ln c between two basins, in the surface's units.</summary>
        public double LnFlow(int[] from, int[] to) => _flow[(_ids[K(from)], _ids[K(to)])] + _total;

        /// <summary>ln π(S)R with each basin shorted, from a state of <paramref name="chain"/> to a basin; null within one basin.</summary>
        public double LnWaitShort(ReversibleChain chain, int from, int[] target)
        {
            int a = _basinOf[from], b = _ids[K(target)];
            if (a == b) return double.NaN;
            if (!_escape.TryGetValue((a, b), out double lnEscape))
                _escape[(a, b)] = lnEscape = _reduction.Toward([b], [a]).Read(a)!.Value.Reading.LnEscape;
            return chain.LogPi[from] - Chain.LogPi[a] - lnEscape;
        }
    }

    // ---- The parts' shared fixture ----

    private sealed record S10(BassMixing.Fixture Fixture, PathSurface Surface, Dictionary<(string Label, string Rank), BassPathProfileTests.Exemplar> Modes)
    {
        public BassPathProfileTests.Exemplar A1 => Modes[("a1", "1")];
    }

    private static S10 Load()
    {
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        return new S10(fixture, new PathSurface(fixture.Config, fixture.Data), BassPathProfileTests.ReadModes());
    }

    private static StepRow[] Rows(string label, string alignment) => StepSeries.Read("500", "S10", PairName[label], "1/1", alignment);

    // ---- Part 1 ----

    [Fact]
    public void Part1_ExactResistances()
    {
        if (!BassQualification.Enabled) return;   // heavy: M ≤ 4 holds about 12,000 states
        Quietly();
        PathCsv csv = Csv("exact", ExactHeader, 5);
        string[][] stored = [.. new PathCsv("exact", ["cap"], 5, "waiting-time").Rows()];
        var pairs = new PathCsv("groundtruth", ["fixture"], 5).Rows().Where(r => r[0] == "grid").Select(r => (Cap: r[1], Seed: r[2], A: int.Parse(r[3], CultureInfo.InvariantCulture), B: int.Parse(r[4], CultureInfo.InvariantCulture))).ToArray();
        var worst = new double[3];
        Array.Fill(worst, double.PositiveInfinity);
        foreach (var fixture in stored.GroupBy(r => (Cap: r[0], Seed: r[1])).OrderBy(g => int.Parse(g.Key.Cap, CultureInfo.InvariantCulture)).ThenBy(g => int.Parse(g.Key.Seed, CultureInfo.InvariantCulture)))
            ExactFixture(csv, fixture.Key.Cap, fixture.Key.Seed, [.. pairs.Where(p => p.Cap == fixture.Key.Cap && p.Seed == fixture.Key.Seed).Select(p => (p.A, p.B))], [.. fixture], worst);
        Validation().Append("rayleigh", $"part 1 {Stamp}", F.I(csv.Rows().Count()), F.F(Math.Min(worst[0], worst[1])), F.F(worst[0]), F.F(worst[1]), F.F(worst[2]), "", "least margins, in nats: R_eff over R_short; R_route over R_eff; R_lump over R_short. Any negative stops the run");
        Ledger("part 1");
    }

    private void ExactFixture(PathCsv csv, string cap, string seed, (int A, int B)[] pairs, string[][] rows, double[] worst)
    {
        if (rows.All(r => csv.Has(r[..5]))) return;
        var clock = Stopwatch.StartNew();
        BassConfig config = BassLocalStepTests.Continuous();
        config = config with { Structure = config.Structure with { MaxBasis = int.Parse(cap, CultureInfo.InvariantCulture) } };
        BassData data = BassSweepTests.Data(false, int.Parse(seed, CultureInfo.InvariantCulture));
        var surface = new PathSurface(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        ReversibleChain chain = ReversibleChain.Of([.. space.States.Select(s => s.LogPi)], flows.Values);
        int[][] codes = [.. space.States.Select(s => BassWaitingTime.Code(s.Terms))];
        Lumped lumped = Lumped.Of(space, chain, codes);
        var routes = new Dictionary<(int, int), (PathRoute? Exact, PathRoute? Aligned)>();
        foreach ((int a, int b) in pairs)
            routes[(a, b)] = (F.Build(surface, space.States[a].Terms, space.States[b].Terms, true, PathReading.Held, []), F.Build(surface, space.States[a].Terms, space.States[b].Terms, false, PathReading.Held, []));

        var own = new Dictionary<string, double>();
        var swap = new Dictionary<string, double>();
        double Own(BasinEdge e)
        {
            string key = $"{BassMixing.Key(e.Lower)}|{BassWaitingTime.Mask(e.Term)}";
            if (!own.TryGetValue(key, out double v)) own[key] = v = ClassSums.Own(surface, e.Lower, e.Term, null).Ln;
            return v;
        }
        double Swap(BasinEdge e, BassTerm t)
        {
            string key = $"{BassMixing.Key(e.Lower)}|{BassMixing.Key([t])}";
            if (!swap.TryGetValue(key, out double v)) swap[key] = v = ClassSums.Swap(surface, e.Lower, t, null).Ln;
            return v;
        }

        foreach (string[] r in rows)
        {
            if (csv.Has(r[..5])) continue;
            int from = int.Parse(r[2], CultureInfo.InvariantCulture), to = int.Parse(r[3], CultureInfo.InvariantCulture);
            string kind = r[4];
            var row = new string[ExactHeader.Length];
            Array.Fill(row, "");
            r[..5].CopyTo(row, 0);
            row[5] = F.I(space.States[from].Terms.Length);
            row[6] = F.I(space.States[to].Terms.Length);
            row[7] = F.I(space.States.Count);
            row[8] = F.I(space.Terms.Count);
            row[9] = F.I(lumped.Basins);
            if (r[9].Length == 0)
            {
                row[^1] = "start in the target's basin";
                csv.Append(row);
                continue;
            }
            double lnH = D(r[9]), eff = D(r[10]);
            row[10] = F.F(lnH);
            row[11] = F.F(eff);
            if (codes[from].SequenceEqual(codes[to]))
            {
                row[^1] = "structure target in the start's basin: no R_short";
                csv.Append(row);
                continue;
            }
            double shorted = lumped.LnWaitShort(chain, from, codes[to]);
            row[12] = F.F(shorted);
            if (shorted > eff + 2e-5) throw new InvalidOperationException($"Rayleigh fails: cap {cap} seed {seed} {from} → {to} {kind}: R_short reads {shorted}, R_eff {eff}.");
            worst[0] = Math.Min(worst[0], eff - shorted);
            bool forward = routes.ContainsKey((from, to));
            (PathRoute? exact, PathRoute? aligned) = forward ? routes[(from, to)] : routes[(to, from)];
            int[]? cut = kind == "basin" ? codes[to] : null;
            void Read(PathRoute? route, int column, int stored, int width)
            {
                if (route is null || route.Hops.Length == 0 || r[stored].Length == 0) return;
                BasinRoute walk = BasinSeries.FromRoute(route, forward, cut);
                double lnRoute = D(r[stored]);
                if (Math.Abs(BasinSeries.Route(walk) - lnRoute) > 2e-5) throw new InvalidOperationException($"cap {cap} seed {seed} {from} → {to} {kind}: the walk reads {BasinSeries.Route(walk)}, the stored route {lnRoute}.");
                double lump = BasinSeries.Lumped(walk, lumped.LnFlow), lnOwn = BasinSeries.Own(walk, Own), held = BasinSeries.Held(walk, Own, Swap);
                double held25 = BasinSeries.CreditedToLater(walk, e => (e.LnFlowClass, e.ClassSize));
                if (eff > lnRoute + 2e-5 || shorted > lump + 1e-8 || lump > lnOwn + 1e-8 || lnOwn > lnRoute + 2e-5)
                    throw new InvalidOperationException($"Rayleigh fails: cap {cap} seed {seed} {from} → {to} {kind}: R_short {shorted}, R_eff {eff}, R_lump {lump}, R_own {lnOwn}, the route {lnRoute}.");
                worst[1] = Math.Min(worst[1], lnRoute - eff);
                worst[2] = Math.Min(worst[2], lump - shorted);
                int returns = walk.Along.Length - walk.Along.Select(a => string.Join(".", BassWaitingTime.Code(a))).Distinct().Count();
                string[] cells = [F.F(lnRoute), r[stored + 1], r[stored + 2], r[stored + 3], F.I(walk.Edges.Length), F.I(returns), F.F(lump), F.F(lnOwn), F.F(held), F.F(held25)];
                if (cells.Length != width) throw new InvalidOperationException("A route's cells moved.");
                cells.CopyTo(row, column);
            }
            Read(exact, 13, 16, ExactRoute.Length);
            Read(aligned, 13 + ExactRoute.Length, 22, ExactRoute.Length);
            csv.Append(row);
        }
        _out.WriteLine($"cap {cap} seed {seed}: {space.States.Count} states, {space.Terms.Count} terms, {lumped.Basins} basins; {clock.Elapsed.TotalSeconds:F0} s.");
    }

    // ---- Part 2 ----

    // How fast a class sums at N = 500, on the deepest route's largest structure.
    [Fact]
    public void Rate_ClassSum()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        BasinRoute walk = BasinSeries.FromSteps(s.Surface, Rows("chain0", "exact"), s.Modes[("chain0", "1")].Terms, s.A1.Terms, BassWaitingTime.Code(s.A1.Terms));
        BasinEdge edge = walk.Edges.Where(e => e.Term.Continuous.Count == 2).MaxBy(e => e.Lower.Length)!;
        var clock = Stopwatch.StartNew();
        ClassRead read = ClassSums.Own(s.Surface, edge.Lower, edge.Term, null, PassageBook.Threads);
        double seconds = clock.Elapsed.TotalSeconds;
        Validation().Append("class rate", $"M {edge.Lower.Length} {Stamp}", F.I(read.Members), F.F(read.Members / seconds), F.F(read.Ln), F.F(read.Admitted), F.F(seconds), F.F(edge.LnFlow), "members a second on four threads; ln of the sum; those admitted; seconds; the route's exact flow");
        _out.WriteLine($"{read.Members:N0} members, {read.Admitted:N0} admitted, {seconds:F1} s: {read.Members / seconds:N0} a second. ln own {read.Ln:F3}, exact flow {edge.LnFlow:F3}.");
    }

    [Fact]
    public void Part2_WholeMeasureSeries()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        PathCsv classes = Csv("classes", ClassHeader, 3), hops = Csv("hops", HopHeader, 3), series = Csv("series", SeriesHeader, 2);
        var held = new Dictionary<string, ClassRead>();
        foreach (string[] r in classes.Rows()) held[$"{r[0]}|{r[1]}|{r[2]}"] = new ClassRead(D(r[3]), long.Parse(r[4], CultureInfo.InvariantCulture), long.Parse(r[5], CultureInfo.InvariantCulture));
        ClassRead Class(string kind, BassTerm[] structure, BassTerm t)
        {
            string key = BassMixing.Key(structure), mask = F.I(BassWaitingTime.Mask(t));
            if (held.TryGetValue($"{kind}|{key}|{mask}", out ClassRead read)) return read;
            var clock = Stopwatch.StartNew();
            read = kind == "own" ? ClassSums.Own(s.Surface, structure, t, null, PassageBook.Threads) : ClassSums.Read(s.Surface, structure, t, null, p => p, PassageBook.Threads);
            classes.Append(kind, key, mask, F.F(read.Ln), F.I(read.Members), F.I(read.Admitted), clock.Elapsed.TotalSeconds.ToString("F1", CultureInfo.InvariantCulture));
            return held[$"{kind}|{key}|{mask}"] = read;
        }
        int[] target = BassWaitingTime.Code(s.A1.Terms);
        // The cheapest routes first, so a budget that stops leaves the costliest: an alignment at a time, exact last.
        foreach (string alignment in Alignments)
            foreach (string label in Shallow.Concat(Deep))
            {
                if (series.Has(label, alignment)) continue;
                double spent = classes.Rows().Sum(r => D(r[6])) * PassageBook.Threads / 3600.0;
                if (spent >= ClassBudgetHours)
                {
                    Validation().Append("budget", $"part 2 {Stamp}", "", F.F(spent), "", "", "", "", $"the class sums' core-hours reached {ClassBudgetHours}: {label} {alignment} and what follows not run");
                    return;
                }
                var clock = Stopwatch.StartNew();
                BassTerm[] start = s.Modes[(label, "1")].Terms;
                StepRow[] rows = Rows(label, alignment);
                BasinRoute cut = BasinSeries.FromSteps(s.Surface, rows, start, s.A1.Terms, target), whole = BasinSeries.FromSteps(s.Surface, rows, start, s.A1.Terms, null);
                double Own(BasinEdge e) => Class("own", e.Lower, e.Term).Ln;
                double Swap(BasinEdge e, BassTerm t) => Class("swap", [.. e.Lower.Where(u => !u.Equals(t))], t).Ln - s.Surface.Score(e.Lower, profile: false).Pi(PathReading.Held);
                int index = 0;
                foreach (BasinEdge e in cut.Edges)
                {
                    ClassRead own = Class("own", e.Lower, e.Term);
                    double[] credits = [.. e.Earlier.Select(t => Swap(e, t))];
                    hops.Append(label, alignment, F.I(index++), e.Kind.ToString().ToLowerInvariant(), BassMixing.Key([e.Term]), F.I(e.Lower.Length), s.Fixture.Basin(e.Kind == PathHopKind.Add ? e.Lower : [.. e.Lower, e.Term]), s.Fixture.Basin(e.Kind == PathHopKind.Add ? [.. e.Lower, e.Term] : e.Lower),
                                F.F(e.LnFlow), F.F(own.Ln), F.I(own.Members), F.I(own.Admitted), F.I(e.Earlier.Length), F.F(credits.Sum()), string.Join(";", e.Earlier.Zip(credits, (t, c) => $"{BassMixing.Key([t])}={F.F(c)}")), F.F(cut.PiFrom - own.Ln - credits.Sum()));
                }
                double lnOwn = BasinSeries.Own(cut, Own), lnHeld = BasinSeries.Held(cut, Own, Swap);
                static (double, long) Stored(BasinEdge e) => (e.LnFlowClass, e.ClassSize);
                series.Append(label, alignment, F.F(cut.PiFrom), F.I(cut.Edges.Length), F.I(cut.Steps.Length - cut.Edges.Length), F.F(BasinSeries.Route(cut)), F.F(BasinSeries.OnAdditions(cut, e => e.LnFlowClass)), F.F(BasinSeries.CreditedToLater(cut, Stored)),
                              F.F(lnOwn), F.F(lnHeld), F.F(BasinSeries.Route(whole)), F.F(BasinSeries.OnAdditions(whole, e => e.LnFlowClass)), F.F(BasinSeries.CreditedToLater(whole, Stored)), clock.Elapsed.TotalSeconds.ToString("F0", CultureInfo.InvariantCulture));
                _out.WriteLine($"{label} {alignment}: route {BasinSeries.Route(cut):F2}, held25 {BasinSeries.CreditedToLater(cut, Stored):F2}, own {lnOwn:F2}, held {lnHeld:F2}; {clock.Elapsed.TotalSeconds:F0} s.");
            }
        Ledger("part 2");
    }

    // ---- Part 3 ----

    [Fact]
    public void Part3_RouteShape()
    {
        if (!BassQualification.Enabled) return;
        S10 s = Load();
        PathCsv csv = Csv("shape", ShapeHeader, 2);
        var latest = new Dictionary<string, string[]>();
        foreach (string[] r in new PathCsv("chains", PassageTally.Header, 6, "waiting-time").Rows()) latest[string.Join("|", r.Take(5))] = r;
        string[][] traces = [.. new PathCsv("traces", ["part"], 6, "waiting-time").Rows().Where(r => r[0] == "3d")];
        int[] target = BassWaitingTime.Code(s.A1.Terms);
        string a1 = s.Fixture.Basin(s.A1.Terms);
        // The basin before a trace's first row, from the row's own terms.
        string Before(string after, string gone, string come)
        {
            var count = new Dictionary<string, int>();
            foreach (string token in after.Split(" + "))
            {
                string[] parts = token.Split(" × ");
                count[parts[^1]] = parts.Length == 2 ? int.Parse(parts[0], CultureInfo.InvariantCulture) : 1;
            }
            if (come.Length > 0)
                foreach (BassTerm t in BassWaitingTime.ParseStructure(come)) count[s.Fixture.Variables(t)]--;
            if (gone.Length > 0)
                foreach (BassTerm t in BassWaitingTime.ParseStructure(gone)) count[s.Fixture.Variables(t)] = count.GetValueOrDefault(s.Fixture.Variables(t)) + 1;
            return string.Join(" + ", count.Where(kv => kv.Value > 0).OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => kv.Value > 1 ? $"{kv.Value} × {kv.Key}" : kv.Key));
        }
        foreach (string label in Deep)
        {
            string[] doors = [.. latest.Values.Where(r => r[0] == "3x" && r[1] == label && r[31].Length > 0).Select(r => r[31])];
            var crossings = traces.Where(r => r[1] == label).GroupBy(r => (Arm: r[2], Chain: r[4])).ToArray();
            foreach (string alignment in Alignments)
            {
                BasinRoute walk = BasinSeries.FromSteps(s.Surface, Rows(label, alignment), s.Modes[(label, "1")].Terms, s.A1.Terms, target);
                string[] along = [.. walk.Along.Select(x => s.Fixture.Basin(x))];
                var on = new HashSet<string>(along);
                string lastDoor = along[^2];
                if (along[^1] != a1) throw new InvalidOperationException($"{label} {alignment}: the cut ends in {along[^1]}.");
                int moves = 0, movesOn = 0;
                foreach (var crossing in crossings)
                {
                    string? before = null;
                    foreach (string[] row in crossing.OrderBy(r => long.Parse(r[5], CultureInfo.InvariantCulture)))
                    {
                        before ??= Before(row[10], row[7], row[8]);
                        moves++;
                        if (on.Contains(before) && on.Contains(row[10])) movesOn++;
                        before = row[10];
                    }
                }
                string Doors(bool wanted) => string.Join(";", doors.Where(d => on.Contains(d) == wanted).GroupBy(d => d).OrderByDescending(g => g.Count()).Select(g => $"{g.Key}={g.Count()}"));
                csv.Append(label, alignment, F.I(on.Count), F.I(doors.Length), F.I(doors.Count(on.Contains)), F.I(doors.Count(d => d == lastDoor)), Doors(true), Doors(false), F.I(crossings.Length), F.I(moves), F.I(movesOn));
                _out.WriteLine($"{label} {alignment}: {doors.Count(on.Contains)} of {doors.Length} gateways on its {on.Count} basins, {doors.Count(d => d == lastDoor)} at its last door; {movesOn} of {moves} moves on it.");
            }
        }
    }

    // ---- The gates ----

    // The credit reduces to what exists: window 0 is the route's reading, ±25 the path profile's class reading and its credited column.
    [Fact]
    public void Gate_CreditReduces()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        PathCsv validation = Validation();
        string[][] readings = [.. new PathCsv("readings", ["part"], 4, "waiting-time").Rows()];
        int[] target = BassWaitingTime.Code(s.A1.Terms);
        double worstZero = 0.0, worstClass = 0.0, worstCredited = 0.0, worstStored = 0.0, worstTable = 0.0;
        int count = 0;
        foreach (string label in Deep.Concat(Shallow))
            foreach (string alignment in Alignments)
            {
                string[] read = readings.Single(r => r[0] is "3s" or "3d" && r[1] == label && r[2] == "a1" && r[3] == alignment);
                StepRow[] rows = Rows(label, alignment);
                BassTerm[] start = s.Modes[(label, "1")].Terms;
                foreach ((BasinRoute walk, double lnRoute, double lnClass) in new[] { (BasinSeries.FromSteps(s.Surface, rows, start, s.A1.Terms, null), D(read[5]), D(read[6])), (BasinSeries.FromSteps(s.Surface, rows, start, s.A1.Terms, target), D(read[7]), D(read[8])) })
                {
                    var at = new Dictionary<(BassTerm, int), ClassRead>();
                    ClassRead Window(BasinEdge e, int w)
                    {
                        if (!at.TryGetValue((e.Term, w), out ClassRead c)) at[(e.Term, w)] = c = ClassSums.Own(s.Surface, e.Lower, e.Term, w, PassageBook.Threads);
                        return c;
                    }
                    (double, long) Live(BasinEdge e, int w) => (Window(e, w).Ln, Window(e, w).Members);
                    foreach (BasinEdge e in walk.Edges)
                    {
                        worstZero = Math.Max(worstZero, Math.Abs(Window(e, 0).Ln - e.LnFlow));
                        worstClass = Math.Max(worstClass, double.IsNaN(e.LnFlowClass) ? 0.0 : Math.Abs(Window(e, F.ClassWindow).Ln - e.LnFlowClass));
                        Assert.Equal(1, Window(e, 0).Members);
                    }
                    worstZero = Math.Max(worstZero, Math.Abs(BasinSeries.CreditedToLater(walk, e => Live(e, 0)) - lnRoute));
                    worstClass = Math.Max(worstClass, Math.Abs(BasinSeries.OnAdditions(walk, e => Window(e, F.ClassWindow).Ln) - lnClass));
                    worstCredited = Math.Max(worstCredited, Math.Abs(BasinSeries.CreditedToLater(walk, e => Live(e, F.ClassWindow)) - BasinSeries.CreditedToLater(walk, e => (e.LnFlowClass, e.ClassSize))));
                    count++;
                }
                BasinRoute full = BasinSeries.FromSteps(s.Surface, rows, start, s.A1.Terms, null);
                if (CreditedTable.TryGetValue((label, alignment), out double table))
                    worstTable = Math.Max(worstTable, Math.Abs(BasinSeries.CreditedToLater(full, e => (e.LnFlowClass, e.ClassSize)) - table));
                worstStored = Math.Max(worstStored, Math.Abs(BasinSeries.Route(full) - D(read[5])));
            }
        validation.Append("credit reduces", $"N = 500 {Stamp}", F.I(count), F.F(Math.Max(Math.Max(worstZero, worstClass), worstCredited)), F.F(worstZero), F.F(worstClass), F.F(worstCredited), F.F(worstTable), "against the stored readings: window 0 and the route; ±25 and the class reading; ±25 credited live and from the stored rows; and that against the path profile's column, printed to one decimal");
        _out.WriteLine($"{count} series: window 0 {worstZero:E2}; ±25 class {worstClass:E2}; credited {worstCredited:E2}; table {worstTable:F3}; replay {worstStored:E2}.");
        Assert.True(worstZero <= 2e-4 && worstClass <= 2e-4 && worstCredited <= 2e-4 && worstTable <= 0.0501 && worstStored <= 2e-4);
    }

    // The fixtures' kernels, the case each check reads: seed, a basin of at least the share given, its share, and its states.
    private sealed record Enumerated(int Seed, PathSurface Surface, BassMixing.Fixture Fixture, BassKernelEnumeration.Space Space, ReversibleChain Chain, int[][] Codes)
    {
        public static Enumerated Of(int seed, int cap)
        {
            BassConfig config = BassLocalStepTests.Continuous();
            config = config with { Structure = config.Structure with { MaxBasis = cap } };
            BassData data = BassSweepTests.Data(false, seed);
            var surface = new PathSurface(config, data);
            BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
            var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
            ReversibleChain chain = ReversibleChain.Of([.. space.States.Select(s => s.LogPi)], flows.Values);
            return new Enumerated(seed, surface, new BassMixing.Fixture("grid", data, config, null), space, chain, [.. space.States.Select(s => BassWaitingTime.Code(s.Terms))]);
        }

        /// <summary>The basins by π, largest first, as (code, share).</summary>
        public (int[] Code, double Share)[] Basins()
        {
            double[] pi = Chain.Pi();
            return [.. Enumerable.Range(0, pi.Length).GroupBy(i => string.Join(".", Codes[i])).Select(g => (Codes[g.First()], g.Sum(i => pi[i]))).OrderByDescending(b => b.Item2)];
        }
    }

    private static HoldSpec CheckSpec(Enumerated e, BassModel model, BassTerm[] start, double logTau, int[] target, int[] returns, int chain, string arm, double[]? levels = null) => new()
    {
        Start = $"seed {e.Seed}", Arm = arm, Chain = chain, Model = model, Terms = start, LogTau = logTau, Target = target, Basin = e.Fixture.Basin,
        Windows = CheckWindows, Levels = levels ?? [0.25, 0.5], Gap = 20, Returns = returns,
    };

    // One walk through the tally, on the enumerated kernel, until its first return.
    private static HoldTally Walk(Enumerated e, HoldSpec spec, int start, Xoshiro256PlusPlus rng, long cap)
    {
        ReversibleChain chain = e.Chain;
        var tally = new HoldTally();
        tally.Begin(spec, e.Space.States[start].Terms);
        int x = start;
        while (tally.Sweeps < cap && tally.Return is null)
        {
            double u = rng.NextDouble();
            bool moved = false;
            for (int k = chain.Starts[x]; k < chain.Starts[x + 1]; k++)
            {
                u -= chain.Rate[k];
                if (u >= 0.0) continue;
                x = chain.To[k];
                moved = true;
                break;
            }
            tally.Tick(moved ? e.Space.States[x].Terms : null);
        }
        return tally;
    }

    // One production chain through the tally until its first return.
    private static HoldTally Produce(HoldSpec spec, int seed, long cap)
    {
        int[] seeds = SeedTree.Derive(seed, 2);
        var chain = new BassChain(spec.Model, spec.Terms, spec.LogTau, new Xoshiro256PlusPlus(seeds[0]), new Xoshiro256PlusPlus(seeds[1]));
        var tally = new HoldTally();
        IReadOnlyList<BassTerm> last = chain.Current.Terms;
        tally.Begin(spec, last);
        while (tally.Sweeps < cap && tally.Return is null)
        {
            IReadOnlyList<BassTerm> terms = chain.Step().Terms;
            if (ReferenceEquals(terms, last)) tally.Tick(null);
            else
            {
                last = terms;
                tally.Tick(terms);
            }
        }
        return tally;
    }

    private static (double Mean, double Se, double Cv) Moments(double[] t)
    {
        double mean = t.Average(), sd = Math.Sqrt(t.Sum(x => (x - mean) * (x - mean)) / (t.Length - 1));
        return (mean, sd / Math.Sqrt(t.Length), sd / mean);
    }

    // The hold and return tallies on the production chain against the same tally on kernel walks, on enumerated pairs: a two-well pair whose larger side can hold.
    [Fact]
    public void Gate_HoldTally()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        PathCsv validation = Validation();
        const int passages = 1000;
        const long cap = 2_000_000;
        var options = new ParallelOptions { MaxDegreeOfParallelism = HoldBook.Threads };
        string[][] rows = [.. new PathCsv("exact", ["cap"], 5, "waiting-time").Rows().Where(r => r[0] == "3" && r[4] == "basin" && r[9].Length > 0 && D(r[9]) is >= 5.0 and <= 10.0)];
        int cases = 0;
        foreach (var bySeed in rows.GroupBy(r => r[1]))
        {
            if (cases >= 2) break;
            Enumerated e = Enumerated.Of(int.Parse(bySeed.Key, CultureInfo.InvariantCulture), 3);
            double[] pi = e.Chain.Pi();
            foreach (string[] row in bySeed.OrderBy(r => D(r[9])))
            {
                int from = int.Parse(row[2], CultureInfo.InvariantCulture), to = int.Parse(row[3], CultureInfo.InvariantCulture);
                int[] target = e.Codes[to];
                double share = Enumerable.Range(0, pi.Length).Where(i => e.Codes[i].SequenceEqual(target)).Sum(i => pi[i]);
                int[] returns = [.. e.Codes[from].Except(target).Distinct()];
                if (share is < 0.55 or > 0.97 || returns.Length == 0) continue;
                HoldSpec Spec(int chain) => CheckSpec(e, e.Surface.Held, e.Space.States[from].Terms, Math.Log(e.Surface.HeldTau), target, returns, chain, "H");
                // A pilot on the kernel: the pair has to hold and return within reach.
                int[] pilot = SeedTree.Derive(9500 + cases, 40);
                var piloted = new HoldTally[pilot.Length];
                Parallel.For(0, pilot.Length, options, i => piloted[i] = Walk(e, Spec(i), from, new Xoshiro256PlusPlus(pilot[i]), cap));
                if (piloted.Any(t => t.Return is null) || piloted.Average(t => t.Sweeps) > 100_000) continue;

                var clock = Stopwatch.StartNew();
                string name = $"seed {e.Seed} {from} to {to}";
                int[] seeds = SeedTree.Derive(9600 + cases, 2 * passages);
                var walks = new HoldTally[passages];
                var chains = new HoldTally[passages];
                Parallel.For(0, passages, options, i => walks[i] = Walk(e, Spec(i), from, new Xoshiro256PlusPlus(seeds[i]), 20 * cap));
                Parallel.For(0, passages, options, i => chains[i] = Produce(Spec(i), seeds[passages + i], 20 * cap));
                foreach ((string reading, Func<HoldTally, double> of) in new (string, Func<HoldTally, double>)[]
                {
                    ("first entry", t => t.FirstEntry!.Value), ("hold", t => t.Hold!.Value), ("return", t => t.Return!.Value - t.PostFrom!.Value), ("entries", t => t.Entries), ("touches", t => t.Touches!.Value), ("hold at (200, 0.25)", t => t.Held[2]!.Value),
                })
                {
                    (double mw, double sw, double cw) = Moments([.. walks.Select(of)]);
                    (double mc, double sc, double cc) = Moments([.. chains.Select(of)]);
                    double z = mw == mc ? 0.0 : Math.Abs(mw - mc) / Math.Sqrt(sw * sw + sc * sc);
                    validation.Append($"hold tally {reading}", name, F.I(passages), F.F(z), F.F(mw), F.F(mc), F.F(sw), F.F(sc), $"{(z <= 2.5758 ? "agrees" : "OUTSIDE")}; kernel walk then production chain; cv {cw:F3} and {cc:F3}");
                    _out.WriteLine($"{name} {reading}: walk {mw:F1} ± {sw:F1}, chain {mc:F1} ± {sc:F1}; z {z:F2}.");
                    Assert.True(z <= 3.5, $"{name} {reading}: walk {mw} ± {sw}, chain {mc} ± {sc}.");
                }
                _out.WriteLine($"{name}: share {share:F3}, {clock.Elapsed.TotalSeconds:F0} s.");
                cases++;
                break;
            }
        }
        Assert.True(cases >= 1, "No enumerated pair both holds and returns within reach.");
    }

    // The long-run share of a basin on the production chain, arm H, against exact π(basin), and arm S's against τ integrated by quadrature.
    [Fact]
    public void Gate_Shares()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        PathCsv validation = Validation();
        const int chains = 16;
        const long sweeps = 3_000_000;
        const double t99 = 2.947;
        var options = new ParallelOptions { MaxDegreeOfParallelism = HoldBook.Threads };
        int[] fixtures = [.. new PathCsv("exact", ["cap"], 5, "waiting-time").Rows().Where(r => r[0] == "3").Select(r => int.Parse(r[1], CultureInfo.InvariantCulture)).Distinct().Order()];
        double[] Run(HoldSpec[] specs, int salt)
        {
            int[] seeds = SeedTree.Derive(9700 + salt, specs.Length);
            var share = new double[specs.Length];
            Parallel.For(0, specs.Length, options, i =>
            {
                int[] pair = SeedTree.Derive(seeds[i], 2);
                var chain = new BassChain(specs[i].Model, specs[i].Terms, specs[i].LogTau, new Xoshiro256PlusPlus(pair[0]), new Xoshiro256PlusPlus(pair[1]));
                var tally = new HoldTally();
                IReadOnlyList<BassTerm> last = chain.Current.Terms;
                tally.Begin(specs[i], last);
                while (tally.Sweeps < sweeps)
                {
                    IReadOnlyList<BassTerm> terms = chain.Step().Terms;
                    if (ReferenceEquals(terms, last)) tally.Tick(null);
                    else
                    {
                        last = terms;
                        tally.Tick(terms);
                    }
                }
                share[i] = tally.PostSweeps == 0 ? double.NaN : tally.PostTarget / (double)tally.PostSweeps;
            });
            return share;
        }
        int read = 0;
        foreach (int seed in fixtures)
        {
            if (read >= 4) break;
            Enumerated e = Enumerated.Of(seed, 3);
            (int[] code, double exact) = e.Basins()[0];
            if (exact is < 0.3 or > 0.97) continue;
            int start = Enumerable.Range(0, e.Codes.Length).Where(i => !e.Codes[i].SequenceEqual(code)).MaxBy(i => e.Chain.LogPi[i]);
            HoldSpec[] specs = [.. Enumerable.Range(0, chains).Select(c => CheckSpec(e, e.Surface.Held, e.Space.States[start].Terms, Math.Log(e.Surface.HeldTau), code, [1, 2], c, "H", [0.05, 0.1]))];
            var clock = Stopwatch.StartNew();
            double[] share = Run(specs, read);
            (double mean, double se, _) = Moments(share);
            bool agrees = Math.Abs(mean - exact) <= t99 * se;
            validation.Append("share, arm H", $"seed {seed}", F.I(chains), F.F(Math.Abs(mean - exact) / se), F.F(exact), F.F(mean), F.F(se), F.F(sweeps), $"{(agrees ? "agrees" : "OUTSIDE")}; exact π(basin), the chains' mean and its se, sweeps a chain; 99%");
            _out.WriteLine($"seed {seed}, arm H: exact {exact:F4}, chains {mean:F4} ± {se:F4}, {clock.Elapsed.TotalSeconds:F0} s.");
            Assert.True(Math.Abs(mean - exact) <= 1.3 * t99 * se, $"seed {seed}: exact {exact}, chains {mean} ± {se}.");
            read++;
        }
        Assert.True(read >= 3, $"{read} fixtures have a basin the check can read.");

        // Arm S: τ integrated by quadrature over π(S, τ) on one fixture.
        {
            int seed = fixtures[0];
            BassConfig config = BassSweepTests.Config();
            BassData data = BassSweepTests.Data(false, seed);
            var surface = new PathSurface(config, data);
            Assert.True(surface.Model.Target.SamplesTau);
            BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Held);
            BassTarget target = surface.Sampled!.Target;
            double centre = Math.Log(surface.HeldTau);
            const int grid = 1200;
            const double low = -25.0, step = 0.05;
            var byBasin = new Dictionary<string, double>();
            var codes = new Dictionary<string, int[]>();
            foreach (BassKernelEnumeration.State state in space.States)
            {
                BassStructureFit fit = surface.Score(state.Terms, profile: false).Fit;
                double mass = double.NegativeInfinity;
                for (int k = 0; k <= grid; k++)
                {
                    BassEvaluation at = target.Evaluate(fit, centre + low + k * step);
                    mass = LogMath.Add(mass, at.LogMarginal + at.LogPrior);
                }
                int[] code = BassWaitingTime.Code(state.Terms);
                string key = string.Join(".", code);
                codes[key] = code;
                byBasin[key] = LogMath.Add(byBasin.GetValueOrDefault(key, double.NegativeInfinity), mass);
            }
            double total = LogMath.Sum([.. byBasin.Values]);
            (string top, double lnTop) = byBasin.MaxBy(kv => kv.Value) is var best ? (best.Key, best.Value) : default;
            double exact = Math.Exp(lnTop - total);
            var fixture = new BassMixing.Fixture("grid", data, config, null);
            BassTerm[] start = space.States.Where(st => string.Join(".", BassWaitingTime.Code(st.Terms)) != top).MaxBy(st => st.LogPi)!.Terms;
            HoldSpec[] specs = [.. Enumerable.Range(0, chains).Select(c => new HoldSpec
            {
                Start = $"seed {seed}", Arm = "S", Chain = c, Model = surface.Model, Terms = start, LogTau = centre, Target = codes[top], Basin = fixture.Basin, Windows = CheckWindows, Levels = [0.05, 0.1], Gap = 20,
            })];
            var clock = Stopwatch.StartNew();
            double[] share = Run(specs, 50);
            (double mean, double se, _) = Moments(share);
            bool agrees = Math.Abs(mean - exact) <= t99 * se;
            validation.Append("share, arm S", $"seed {seed}", F.I(chains), F.F(Math.Abs(mean - exact) / se), F.F(exact), F.F(mean), F.F(se), F.F(sweeps), $"{(agrees ? "agrees" : "OUTSIDE")}; τ integrated over {grid + 1} points of ln τ, the chains' mean and its se, sweeps a chain; 99%");
            _out.WriteLine($"seed {seed}, arm S: quadrature {exact:F4}, chains {mean:F4} ± {se:F4}, {clock.Elapsed.TotalSeconds:F0} s.");
            Assert.True(Math.Abs(mean - exact) <= 1.3 * t99 * se, $"arm S: quadrature {exact}, chains {mean} ± {se}.");
        }
    }

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
        Validation().Append("starts", $"fresh design {Stamp}", F.I(count), F.F(worst), "", "", "", "", "each start and a1's exemplar at the held τ against the stored π");
        _out.WriteLine($"{count} structures: worst {worst:E2}.");
        Assert.True(worst <= 2e-6);
    }

    // ---- The flat histogram's sides and shares, before any chain ----

    [Fact]
    public void FlatHistogram_SidesAndShares()
    {
        if (!BassQualification.Enabled) return;   // heavy: arm (i) again, about ten minutes
        Quietly();
        S10 s = Load();
        PathCsv csv = Csv("flathistogram", FlatHeader, 2);
        if (csv.Has("tau", "the run's τ")) return;
        BassConfig config = BassFlatHistogramQualificationTests.ArmConfig(true, null, BassFlatHistogramQualificationTests.ArmOneSeed);
        Assert.True(BassModel.Create(config, s.Fixture.Data).Target.SamplesTau, "The run's τ is held: arm S's share has no comparison.");
        BassFit fit = Bass.Fit(config, s.Fixture.Data);
        var draws = fit.Diagnostics.FlatHistogram!.Draws!;
        double[] weights = fit.DrawWeights!;
        int chains = fit.Draws.Count;
        var basins = new Dictionary<string, (double[] W, long N, bool Both)>();
        var families = new Dictionary<int, (double[] W, long N)>();
        double[] Blank() => new double[chains];
        for (int c = 0; c < chains; c++)
            for (int d = 0; d < draws.Draws; d++)
            {
                BassTerm[] terms = [.. fit.Draws[c][d].Structure.TermIds.Select(id => fit.Terms[id])];
                int[] code = BassWaitingTime.Code(terms);
                string basin = s.Fixture.Basin(terms);
                double w = weights[c * draws.Draws + d];
                bool both = code.Contains(1) && code.Contains(2);
                int n1 = code.Count(m => m == 1), n2 = code.Count(m => m == 2), n3 = code.Count(m => m == 3);
                int family = Math.Min(n1, 7) + 8 * Math.Min(n2, 7) + 64 * Math.Min(n3, 7);
                (double[] W, long N, bool Both) held = basins.TryGetValue(basin, out var b) ? b : (Blank(), 0, both);
                held.W[c] += w;
                basins[basin] = (held.W, held.N + 1, both);
                (double[] W, long N) fam = families.TryGetValue(family, out var f) ? f : (Blank(), 0);
                fam.W[c] += w;
                families[family] = (fam.W, fam.N + 1);
            }
        // Each chain's weights sum to its share of the run's mass, as the run gives them.
        double[] byChain = [.. Enumerable.Range(0, chains).Select(c => basins.Values.Sum(b => b.W[c]))];
        double total = byChain.Sum();
        double Of(IEnumerable<double[]> set) => set.Sum(w => w.Sum()) / total;
        string PerChain(IEnumerable<double[]> set, Func<int, double> over) => string.Join(";", Enumerable.Range(0, chains).Select(c => (set.Sum(w => w[c]) / over(c)).ToString("F6", CultureInfo.InvariantCulture)));
        var notBoth = basins.Values.Where(b => !b.Both).Select(b => b.W).ToArray();
        var bothSide = basins.Values.Where(b => b.Both).Select(b => b.W).ToArray();
        double massNot = Of(notBoth), massBoth = Of(bothSide);
        double[] notByChain = [.. Enumerable.Range(0, chains).Select(c => notBoth.Sum(w => w[c]))];
        string a1 = s.Fixture.Basin(s.A1.Terms);
        var neither = families.Where(f => (f.Key & 7) == 0 && (f.Key >> 3 & 7) == 0).Select(f => f.Value.W).ToArray();
        csv.Append("tau", "the run's τ", "", "", "", "", "", "", "drawn: the model's coefficient prior samples τ, so arm S's share is the one to compare");
        csv.Append("mass", "both main effects", "", "yes", F.F(massBoth), "", F.I(basins.Values.Where(b => b.Both).Sum(b => b.N)), PerChain(bothSide, c => byChain[c]), "the return's side");
        csv.Append("mass", "not both", "", "no", F.F(massNot), "", F.I(basins.Values.Where(b => !b.Both).Sum(b => b.N)), PerChain(notBoth, c => byChain[c]), "a1's side, the primary reading");
        csv.Append("mass", "neither main effect", "", "no", F.F(Of(neither)), "", F.I(families.Where(f => (f.Key & 7) == 0 && (f.Key >> 3 & 7) == 0).Sum(f => f.Value.N)), PerChain(neither, c => byChain[c]), "a narrower side, beside it");
        csv.Append("mass", "a1's basin", a1, "no", F.F(Of([basins[a1].W])), F.F(basins[a1].W.Sum() / total / massNot), F.I(basins[a1].N), PerChain([basins[a1].W], c => byChain[c]), "the narrowest side");
        csv.Append("ratio", "ln(not both / both)", "", "", F.F(Math.Log(massNot / massBoth)), "", "", "", "the balance's primary ratio: ln(return) ≈ ln(first entry) + this");
        csv.Append("ratio", "ln(neither / both)", "", "", F.F(Math.Log(Of(neither) / massBoth)), "", "", "", "beside it");
        csv.Append("ratio", "ln(a1's basin / both)", "", "", F.F(Math.Log(Of([basins[a1].W]) / massBoth)), "", "", "", "beside it");
        csv.Append("conditional", "a1's basin given not both", a1, "no", F.F(Of([basins[a1].W])), F.F(basins[a1].W.Sum() / notBoth.Sum(w => w.Sum())), F.I(basins[a1].N), PerChain([basins[a1].W], c => notByChain[c]), "by chain, each chain's own conditional");
        int place = 0;
        foreach ((string name, (double[] W, long N, bool Both) held) in basins.OrderByDescending(b => b.Value.W.Sum()).Take(12))
        {
            place++;
            csv.Append("largest", $"largest {place}", name, held.Both ? "yes" : "no", F.F(held.W.Sum() / total), held.Both ? "" : F.F(held.W.Sum() / notBoth.Sum(w => w.Sum())), F.I(held.N), held.Both ? PerChain([held.W], c => byChain[c]) : PerChain([held.W], c => notByChain[c]), held.Both ? "by chain, the share" : "by chain, the conditional");
        }
        foreach ((int family, (double[] W, long N) held) in families.OrderByDescending(f => f.Value.W.Sum()).Take(30))
            csv.Append("family", HoldTally.FamilyName(family), "", (family & 7) > 0 && (family >> 3 & 7) > 0 ? "yes" : "no", F.F(held.W.Sum() / total), "", F.I(held.N), PerChain([held.W], c => byChain[c]), "(x1, x2, x1·x2)");
        _out.WriteLine($"mass: both {massBoth:F5}, not both {massNot:F5}; ln ratio {Math.Log(massNot / massBoth):F3}; a1 {Of([basins[a1].W]):F5}, given not both {basins[a1].W.Sum() / notBoth.Sum(w => w.Sum()):F5}; total weight {total:F6}.");
        Ledger("flat histogram");
    }

    // ---- Parts 4 and 5 ----

    private static int[][] Largest()
    {
        string[][] rows = [.. Csv("flathistogram", FlatHeader, 2).Rows().Where(r => r[0] == "largest")];
        if (rows.Length != 12) throw new InvalidOperationException("The flat histogram's sides and shares come before any chain.");
        return [.. rows.Select(r => BasinSeries.CodeOf(r[2]))];
    }

    private static HoldSpec Chain(S10 s, string label, string arm, int chain, int[][] largest)
    {
        BassPathProfileTests.Exemplar e = s.Modes[(label, "1")];
        return new HoldSpec
        {
            Start = label, Arm = arm, Chain = chain, Model = arm == "H" ? s.Surface.Held : s.Surface.Model, Terms = e.Terms, LogTau = arm == "H" ? Math.Log(s.Surface.HeldTau) : s.Surface.Score(e.Terms).BestLogTau,
            Target = BassWaitingTime.Code(s.A1.Terms), Basin = s.Fixture.Basin, Largest = largest,
        };
    }

    // Sweeps a second on both arms under four concurrent chains, for the plan's cost.
    [Fact]
    public void Rate_BothArms()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        int[][] largest = Largest();
        HoldSpec[] specs = [Chain(s, "A", "S", 900, largest), Chain(s, "chain1", "S", 901, largest), Chain(s, "A", "H", 900, largest), Chain(s, "chain1", "H", 901, largest)];
        const long sweeps = 1_000_000;
        var rate = new double[specs.Length];
        Parallel.For(0, specs.Length, new ParallelOptions { MaxDegreeOfParallelism = HoldBook.Threads }, i =>
        {
            var tally = new HoldTally();
            tally.Advance(specs[i], sweeps);
            rate[i] = tally.Sweeps / tally.Seconds;
            Validation().Append("rate", $"{specs[i].Start} arm {specs[i].Arm} {Stamp}", F.I(tally.Sweeps), F.F(rate[i]), F.F(tally.Accepted / (double)tally.Sweeps), F.F(tally.Seconds), "", "", "sweeps a second under four concurrent chains at below-normal priority, with the hold tally on; the acceptance");
        });
        double armS = (rate[0] + rate[1]) / 2, armH = (rate[2] + rate[3]) / 2;
        double hoursS25 = 48 * 25.0 * Million / armS / 3600, hoursH25 = 24 * 25.0 * Million / armH / 3600;
        double hoursS100 = 48 * 75.0 * Million / armS / 3600, hoursH100 = 24 * 75.0 * Million / armH / 3600;
        Validation().Append("plan", $"cost {Stamp}", "72", F.F(hoursS25 + hoursH25 + hoursS100 + hoursH100), F.F(hoursS25), F.F(hoursH25), F.F(hoursS100), F.F(hoursH100), "core-hours: arm S to 2.5e7; arm H to 2.5e7; arm S on to 1e8; arm H on to 1e8. The budget is 24");
        _out.WriteLine($"arm S {armS:N0} a second, arm H {armH:N0}. Plan: S to 2.5e7 {hoursS25:F1}, H to 2.5e7 {hoursH25:F1}, S to 1e8 {hoursS100:F1} more, H to 1e8 {hoursH100:F1} more; {hoursS25 + hoursH25 + hoursS100 + hoursH100:F1} of 24 core-hours.");
    }

    private static (HoldSpec[] ArmS, HoldSpec[] ArmH) Plan(S10 s)
    {
        int[][] largest = Largest();
        return ([.. Deep.SelectMany(l => Enumerable.Range(0, 8).Select(c => Chain(s, l, "S", c, largest)))],
                [.. new[] { "A", "B", "chain0" }.SelectMany(l => Enumerable.Range(0, 8).Select(c => Chain(s, l, "H", c, largest)))]);
    }

    [Fact]
    public void Part4_ChainsToAQuarter()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        (HoldSpec[] armS, HoldSpec[] armH) = Plan(s);
        var book = new HoldBook();
        if (book.Advance(armS, 25 * Million, _out)) book.Advance(armH, 25 * Million, _out);
    }

    [Fact]
    public void Part4_ChainsToAHundredMillion()
    {
        if (!BassQualification.Enabled) return;
        Quietly();
        S10 s = Load();
        (HoldSpec[] armS, HoldSpec[] armH) = Plan(s);
        var book = new HoldBook();
        // Rungs of 5 × 10⁶, so where the budget stops an arm its starts are level: a rung to 5 × 10⁷ would run four starts' chains and leave the rest at 2.5 × 10⁷.
        foreach (HoldSpec[] arm in new[] { armS, armH })
            for (long limit = 30 * Million; limit <= 100 * Million; limit += 5 * Million)
                if (!book.Advance(arm, limit, _out)) return;
    }
}
