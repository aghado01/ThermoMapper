using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>A sum over a term's class, and how many members it read and admitted.</summary>
internal readonly record struct ClassRead(double Ln, long Members, long Admitted);

/// <summary>
/// Sums over a term's class (the basin-hold brief; gap G1, reading (f)). With no window the class is whole measure: every sign and knot on the term's
/// variables. With a window, its signs and the knots within it per factor, as <see cref="BassPathProfile.Neighbours"/>. Every member is loaded on the
/// structure's design, so admission is production's.
/// </summary>
internal static class ClassSums
{
    /// <summary>Past this a class would need the brief's sample; no route's does, since no route holds a three-factor term.</summary>
    internal const long MaxMembers = 5_000_000;

    private static HingeFactor[][] Options(PathSurface surface, BassTerm t, int? window)
    {
        if (t.Categorical.Count != 0 || t.Functional.Count != 0) throw new ArgumentException("These fixtures hold continuous terms alone.", nameof(t));
        return [.. t.Continuous.Select(h =>
        {
            int n = BassMixing.Knots(surface.Model, BassPart.Continuous, h.Variable);
            if (window is not { } w) return new[] { -1, 1 }.SelectMany(sign => Enumerable.Range(0, n).Select(k => new HingeFactor(h.Variable, sign, k))).ToArray();
            int low = Math.Max(0, h.Knot - w), high = Math.Min(n - 1, h.Knot + w);
            return Enumerable.Range(low, high - low + 1).Select(k => h with { Knot = k }).ToArray();
        })];
    }

    /// <summary>ln Σ exp(<paramref name="summand"/>(π(rest + t′))) over t's class, t′ not already in <paramref name="rest"/>, at the held τ. The term's measure is one number across a class.</summary>
    internal static ClassRead Read(PathSurface surface, IReadOnlyList<BassTerm> rest, BassTerm t, int? window, Func<double, double> summand, int threads = 1)
    {
        HingeFactor[][] options = Options(surface, t, window);
        long count = options.Aggregate(1L, (a, o) => a * o.Length);
        if (count > MaxMembers) throw new InvalidOperationException($"{BassMixing.Key([t])}'s class holds {count:N0} members, past the whole sum's reach.");
        double logTermPrior = surface.LogTermPrior(rest) + surface.Held.Prior.LogPrior(t);
        double sum = double.NegativeInfinity;
        long members = 0, admitted = 0;
        object gate = new();
        void Range(long from, long to)
        {
            BassDesign design = surface.Design;
            if (design.Load(rest) is { } refusal) throw new InvalidOperationException($"The class's structure is refused ({refusal}).");
            double local = double.NegativeInfinity;
            long seen = 0, kept = 0;
            var hinges = new HingeFactor[options.Length];
            for (long index = from; index < to; index++)
            {
                long left = index;
                for (int f = options.Length - 1; f >= 0; f--)
                {
                    hinges[f] = options[f][left % options[f].Length];
                    left /= options[f].Length;
                }
                var near = new BassTerm(hinges, t.Categorical, null);
                if (rest.Contains(near)) continue;
                seen++;
                DesignProposal proposal = design.ProposeBirth(near);
                PathScore score = surface.Score(BassStructureFit.Of(proposal, proposal.Admitted ? logTermPrior : double.NaN), profile: false);
                if (!score.Admitted) continue;
                kept++;
                local = LogMath.Add(local, summand(score.Pi(PathReading.Held)));
            }
            lock (gate)
            {
                sum = LogMath.Add(sum, local);
                members += seen;
                admitted += kept;
            }
        }
        if (threads <= 1 || count < 4096) Range(0, count);
        else Parallel.ForEach(Partitioner.Create(0L, count, Math.Max(256L, count / (threads * 16L))), new ParallelOptions { MaxDegreeOfParallelism = threads }, r => Range(r.Item1, r.Item2));
        return new ClassRead(sum, members, admitted);
    }

    /// <summary>The edge's own class: ln Σ min(π(lower) + ln q_b, π(lower + t′) + ln q_d), with <paramref name="lower"/>, the route's other terms, fixed.</summary>
    internal static ClassRead Own(PathSurface surface, IReadOnlyList<BassTerm> lower, BassTerm t, int? window, int threads = 1)
    {
        double pi = surface.Score(lower, profile: false).Pi(PathReading.Held), birth = surface.LogBirth(lower, lower.Count, t), death = surface.LogDeath(lower.Count + 1);
        return Read(surface, lower, t, window, p => Math.Min(pi + birth, p + death), threads);
    }

    /// <summary>A held term's class: ln Σ π(x with t → t′) / π(x).</summary>
    internal static ClassRead Swap(PathSurface surface, IReadOnlyList<BassTerm> x, BassTerm t, int? window, int threads = 1)
    {
        double pi = surface.Score(x, profile: false).Pi(PathReading.Held);
        BassTerm[] rest = [.. x.Where(u => !u.Equals(t))];
        if (rest.Length != x.Count - 1) throw new ArgumentException("The structure does not hold the term.", nameof(t));
        return Read(surface, rest, t, window, p => p - pi, threads);
    }
}

/// <summary>One basin-changing hop, read along the route: the term added or deleted, the edge's lower state (the route's other terms), the terms the route added
/// earlier and still holds there, the basins each side, the exact flow, and the ±25 class flow and size the path profile stored.</summary>
internal sealed record BasinEdge(PathHopKind Kind, BassTerm Term, BassTerm[] Lower, BassTerm[] Earlier, int[] From, int[] To, double LnFlow, double LnFlowClass, int ClassSize);

/// <summary>A route's step: a basin-changing edge, or a slide's Σ 1/c_e, which stays inside a basin.</summary>
internal readonly record struct BasinStep(BasinEdge? Edge, double LnSlide);

/// <summary>A route read toward a target, stopped at its first state in the target's basin when one is given: π(S), its steps, and the structures of its basin sequence.</summary>
internal sealed record BasinRoute(double PiFrom, BasinStep[] Steps, BassTerm[][] Along)
{
    public BasinEdge[] Edges => [.. Steps.Where(s => s.Edge is not null).Select(s => s.Edge!)];
}

/// <summary>Builds a <see cref="BasinRoute"/> from a route's states, keeping what it has added.</summary>
internal sealed class BasinRouteBuilder(double piFrom, BassTerm[] first)
{
    private readonly List<BasinStep> _steps = [];
    private readonly List<BassTerm[]> _along = [first];
    private readonly List<BassTerm> _added = [];

    public void Slide(double lnResistance) => _steps.Add(new BasinStep(null, lnResistance));

    public void SlideTo(BassTerm[] state) => _along[^1] = state;

    public void Hop(BassTerm[] from, BassTerm[] to, double lnFlow, double lnFlowClass, int classSize)
    {
        bool adds = to.Length > from.Length;
        BassTerm term = adds ? to.Single(t => !from.Contains(t)) : from.Single(t => !to.Contains(t));
        BassTerm[] lower = adds ? from : to;
        _steps.Add(new BasinStep(new BasinEdge(adds ? PathHopKind.Add : PathHopKind.Delete, term, lower, [.. _added.Where(lower.Contains)], BassWaitingTime.Code(from), BassWaitingTime.Code(to), lnFlow, lnFlowClass, classSize), 0.0));
        _along.Add(to);
        if (adds) _added.Add(term);
        else _added.Remove(term);
    }

    public BasinRoute Build() => new(piFrom, [.. _steps], [.. _along]);
}

/// <summary>The series a route's steps read (the basin-hold brief). How each kind of edge is credited: a slide is dropped by <see cref="Own"/> and <see cref="Held"/>, since a hop inside a basin is
/// not a step of the sequence, and kept by the readings that reduce to the path profile's. An addition and a deletion take the same c_e, the sum over the upper state's own term with the lower state fixed, so
/// an edge reads alike from either end. The credit for the terms added earlier is read at the lower state, the one c_e's sum holds fixed, one term at a time. All are ln π(S)R, read at the held τ.</summary>
internal static class BasinSeries
{
    /// <summary>The route's own reading: exact flows, slides kept.</summary>
    public static double Route(BasinRoute r)
    {
        double w = double.NegativeInfinity;
        foreach (BasinStep s in r.Steps) w = LogMath.Add(w, r.PiFrom + (s.Edge is { } e ? -e.LnFlow : s.LnSlide));
        return w;
    }

    /// <summary>The route's class reading at any window: the class flow on additions, exact elsewhere, slides kept; the path profile's <c>LnWaitClass</c> at ±25.</summary>
    public static double OnAdditions(BasinRoute r, Func<BasinEdge, double> lnClass)
    {
        double w = double.NegativeInfinity;
        foreach (BasinStep s in r.Steps)
            w = LogMath.Add(w, r.PiFrom + (s.Edge is not { } e ? s.LnSlide : e.Kind == PathHopKind.Add ? -Math.Max(e.LnFlow, lnClass(e)) : -e.LnFlow));
        return w;
    }

    /// <summary>W_held at a window, as the path profile's "every addition's class credited to all later edges": the class flow on an addition, and each earlier addition's class size, ln of its members, on every later step. Its π ratios are taken as equal.</summary>
    public static double CreditedToLater(BasinRoute r, Func<BasinEdge, (double LnClass, long Size)> at)
    {
        double w = double.NegativeInfinity, credit = 0.0;
        foreach (BasinStep s in r.Steps)
        {
            if (s.Edge is not { } e)
            {
                w = LogMath.Add(w, r.PiFrom + s.LnSlide - credit);
                continue;
            }
            if (e.Kind != PathHopKind.Add)
            {
                w = LogMath.Add(w, r.PiFrom - e.LnFlow - credit);
                continue;
            }
            (double lnClass, long size) = at(e);
            w = LogMath.Add(w, r.PiFrom - Math.Max(e.LnFlow, lnClass) - credit);
            credit += Math.Log(Math.Max(1L, size));
        }
        return w;
    }

    /// <summary>W_own: over the basin sequence, c_e the sum over the edge's own class.</summary>
    public static double Own(BasinRoute r, Func<BasinEdge, double> own)
    {
        double w = double.NegativeInfinity;
        foreach (BasinEdge e in r.Edges) w = LogMath.Add(w, r.PiFrom - own(e));
        return w;
    }

    /// <summary>W_held: W_own, each c_e also multiplied, for each earlier term still held, by its class's π ratios at the lower state.</summary>
    public static double Held(BasinRoute r, Func<BasinEdge, double> own, Func<BasinEdge, BassTerm, double> swap)
    {
        double w = double.NegativeInfinity;
        foreach (BasinEdge e in r.Edges) w = LogMath.Add(w, r.PiFrom - own(e) - e.Earlier.Sum(t => swap(e, t)));
        return w;
    }

    /// <summary>The series over the route's basin sequence with lumped flows: ln c(from, to) given.</summary>
    public static double Lumped(BasinRoute r, Func<int[], int[], double> lnFlow)
    {
        double w = double.NegativeInfinity;
        foreach (BasinEdge e in r.Edges) w = LogMath.Add(w, r.PiFrom - lnFlow(e.From, e.To));
        return w;
    }

    /// <summary>A route from <see cref="BassPathProfile.Build"/>, read from one end, stopped at the first state in <paramref name="cut"/> when given.</summary>
    public static BasinRoute FromRoute(PathRoute route, bool forward, int[]? cut)
    {
        PathHop[] hops = BassWaitingTime.Directed(route, forward);
        int[] masks = forward ? route.Masks : [.. Enumerable.Reverse(route.Masks)];
        var builder = new BasinRouteBuilder(hops[0].PiFrom, route.Lattice.Terms(masks[0]));
        int at = 0;
        foreach (PathHop h in hops)
        {
            if (h.Kind == PathHopKind.Slide)
            {
                builder.Slide(-h.LnFlow);
                continue;
            }
            BassTerm[] from = route.Lattice.Terms(masks[at]), to = route.Lattice.Terms(masks[at + 1]);
            at++;
            builder.Hop(from, to, h.LnFlow, h.LnFlowClass, h.ClassSize);
            if (cut is not null && BassWaitingTime.Code(to).AsSpan().SequenceEqual(cut)) break;
        }
        return builder.Build();
    }

    /// <summary>
    /// A stored route replayed from <c>path-profile-steps.csv</c>, read from S, with each state's π checked against its row, so N = 500's routes take no lattice. A slide taken first moves S's terms that neither T holds nor the
    /// lattice's rows touch onto T's that S does not hold and the rows do not touch; a slide taken last lies past the cut.
    /// </summary>
    public static BasinRoute FromSteps(PathSurface surface, StepRow[] rows, BassTerm[] s, BassTerm[] t, int[]? cut)
    {
        var onLattice = new HashSet<BassTerm>(rows.Where(r => r.Kind != "slide").Select(r => BassWaitingTime.ParseTerm(r.Term)));
        BassTerm[] x = s;
        var builder = new BasinRouteBuilder(rows[0].PiFrom, x);
        void Check(BassTerm[] state, double pi)
        {
            double read = surface.Score(state, profile: false).Pi(PathReading.Held);
            if (Math.Abs(read - pi) > 2e-5) throw new InvalidOperationException($"A replayed state reads {read}, its row {pi}.");
        }
        Check(x, rows[0].PiFrom);
        for (int k = 0; k < rows.Length; k++)
        {
            StepRow row = rows[k];
            if (row.Kind == "slide")
            {
                builder.Slide(row.LnResistance);
                x = k == 0 ? [.. s.Where(u => t.Contains(u) || onLattice.Contains(u)), .. t.Where(u => !s.Contains(u) && !onLattice.Contains(u))] : t;
                Check(x, row.PiTo);
                builder.SlideTo(x);
                continue;
            }
            BassTerm term = BassWaitingTime.ParseTerm(row.Term);
            BassTerm[] to = row.Kind == "add" ? [.. x, term] : [.. x.Where(u => !u.Equals(term))];
            if (to.Length != x.Length + (row.Kind == "add" ? 1 : -1)) throw new InvalidOperationException($"A row {row.Kind}s {row.Term} on a state that cannot take it.");
            Check(to, row.PiTo);
            builder.Hop(x, to, row.LnFlow, row.LnFlowClass, row.ClassSize);
            x = to;
            if (cut is not null && BassWaitingTime.Code(to).AsSpan().SequenceEqual(cut)) break;
        }
        return builder.Build();
    }

    /// <summary>The sorted masks of a basin as <see cref="BassMixing.Fixture.Basin"/> wrote it: names are one above their variables.</summary>
    public static int[] CodeOf(string basin)
    {
        if (basin == "∅") return [];
        var code = new List<int>();
        foreach (string token in basin.Split(" + "))
        {
            string[] parts = token.Split(" × ");
            int times = parts.Length == 2 ? int.Parse(parts[0], CultureInfo.InvariantCulture) : 1;
            int mask = 0;
            foreach (string name in parts[^1].Split('·')) mask |= 1 << (int.Parse(name[1..], CultureInfo.InvariantCulture) - 1);
            for (int i = 0; i < times; i++) code.Add(mask);
        }
        code.Sort();
        return [.. code];
    }
}

/// <summary>A hold chain: its name, model and start, the target's basin, and what is tallied. The defaults are the brief's: holds at W of 2,000, 10⁴ and 10⁵ and c of ¼ and ½, primary (10⁵, ½), and a return to a state holding an x1 term and an x2 term.</summary>
internal sealed class HoldSpec
{
    public required string Start { get; init; }
    public required string Arm { get; init; }
    public required int Chain { get; init; }
    public required BassModel Model { get; init; }
    public required BassTerm[] Terms { get; init; }
    public required double LogTau { get; init; }
    public required int[] Target { get; init; }
    public required Func<IEnumerable<BassTerm>, string> Basin { get; init; }

    /// <summary>The flat histogram's largest basins, whose occupancy is kept.</summary>
    public int[][] Largest { get; init; } = [];

    public int[] Windows { get; init; } = [2_000, 10_000, 100_000];
    public double[] Levels { get; init; } = [0.25, 0.5];

    /// <summary>The masks a state must hold to count as a return.</summary>
    public int[] Returns { get; init; } = [1, 2];

    /// <summary>The masks whose counts name a state's family: (x1, x2, x1·x2).</summary>
    public int[] Family { get; init; } = [1, 2, 3];

    /// <summary>A visit ends, and a stretch outside the start's basin counts as a departure, at this many sweeps outside.</summary>
    public int Gap { get; init; } = PassageTally.DepartureLength;

    public long Block { get; init; } = 1_000_000;

    public string Name => $"{Start}|{Arm}|{Chain}";

    public int PrimaryWindow => Windows[^1];

    public double PrimaryLevel => Levels[^1];
}

/// <summary>
/// One chain's cumulative record for the wait to hold, as its rows store it. s_W(t) is the fraction of sweeps in [t, t + W) in the target; a chain holds at (W, c) at the first t ≥ 1 with s_W(t) ≥ c. A visit runs from an entry until
/// <see cref="HoldSpec.Gap"/> sweeps outside, so one spurious term does not end it; a touch is a visit closed before the hold. Past the hold's window, T_hold + W, the share and the first return are read, so the window that selected the chain stays out.
/// Sweeps are fed one at a time, so a kernel walk reads through the same tally as the production chain. Region counts run from a departure's recognition, its first <see cref="HoldSpec.Gap"/> sweeps left out. A row stores the sliding window's recent
/// stretches and the open visit, so a chain resumed from a row reads as if it had not stopped; its closed visits are in the visits file, and <see cref="HoldFrom"/> reads blank if the hold's start lies in a visit closed before a restart.
/// </summary>
internal sealed class HoldTally
{
    private const int Families = 512;

    internal static readonly string[] Header =
    [
        "start", "arm", "chain", "sweeps", "status", "segment", "seconds", "accepted", "structure", "log_tau", "first_exit", "departure", "out_since", "first_entry", "entry_from", "stay_sweeps", "stay_seen",
        "entries", "visits", "touches", "hold_from", "hold_w0_c0", "hold_w0_c1", "hold_w1_c0", "hold_w1_c1", "hold_w2_c0", "hold_w2_c1", "last_both_pre", "return", "pre_sweeps", "post_sweeps", "post_notboth", "post_target",
        "pre_family", "pre_largest", "region_sweeps", "region_family", "region_basins", "post_family", "post_largest", "restarts", "recent", "open_visit",
    ];

    internal static readonly string[] VisitHeader = ["start", "arm", "chain", "index", "first", "last", "in_sweeps", "entries", "from"];

    internal static readonly string[] BlockHeader = ["start", "arm", "chain", "block_end", "sweeps", "not_both", "target", "largest"];

    private sealed class Ring(int width)
    {
        public readonly int Width = width;
        public readonly bool[] In = new bool[width];
        public int Count;
    }

    private sealed class Visit(long first, string from)
    {
        public long First = first, Last = first, In, Entries;
        public string From = from;
    }

    private sealed class RegionCell(string label)
    {
        public readonly string Label = label;
        public long Count;
    }

    private HoldSpec _spec = null!;
    private int[] _code = new int[64], _startCode = [];
    private int _m, _family, _largest = -1;
    private bool _inTarget, _inStart, _both, _ringOn;
    private IReadOnlyList<BassTerm>? _current, _before;
    private Ring[]? _rings;
    private string _recent = "", _openVisit = "";
    private Visit? _visit;
    private readonly List<(long First, long Last, string From)> _closed = [];
    private readonly Dictionary<long, RegionCell> _region = [];
    private RegionCell? _cell;
    private long _blockSweeps, _blockNotBoth, _blockTarget;
    private long[] _blockLargest = [];

    public long Sweeps { get; set; }
    public int Segment { get; set; }
    public double Seconds { get; set; }
    public long Accepted { get; set; }
    public string Structure { get; set; } = "";
    public double LogTau { get; set; } = double.NaN;
    public long? FirstExit { get; set; }
    public long? Departure { get; set; }
    public long OutSince { get; set; }
    public long? FirstEntry { get; set; }
    public string EntryFrom { get; set; } = "";
    public long StaySweeps { get; set; }
    public long StaySeen { get; set; }
    public long Entries { get; set; }
    public long Visits { get; set; }
    public long? Touches { get; set; }
    public string HoldFrom { get; set; } = "";
    public long?[] Held { get; private set; } = new long?[6];
    public long? LastBothPre { get; set; }
    public long? Return { get; set; }
    public long PreSweeps { get; set; }
    public long PostSweeps { get; set; }
    public long PostNotBoth { get; set; }
    public long PostTarget { get; set; }
    public long RegionSweeps { get; set; }
    public long[] PreFamily { get; private set; } = new long[Families];
    public long[] RegionFamily { get; private set; } = new long[Families];
    public long[] PostFamily { get; private set; } = new long[Families];
    public long[] PreLargest { get; private set; } = [];
    public long[] PostLargest { get; private set; } = [];
    public int Restarts { get; set; }

    /// <summary>The primary hold's start.</summary>
    public long? Hold => Held[^1];

    /// <summary>The first sweep of the share: T_hold + W.</summary>
    public long? PostFrom => Hold is { } t ? t + _spec.PrimaryWindow : null;

    public List<string[]> PendingVisits { get; } = [];

    public List<string[]> PendingBlocks { get; } = [];

    internal static string FamilyName(int f) => $"{f & 7}/{f >> 3 & 7}/{f >> 6 & 7}";

    private static string Map(long[] family)
    {
        var cells = new List<string>();
        for (int f = 0; f < family.Length; f++)
            if (family[f] > 0) cells.Add($"{FamilyName(f)}={family[f].ToString(CultureInfo.InvariantCulture)}");
        return string.Join(";", cells);
    }

    private static void ReadMap(string cell, long[] family)
    {
        if (cell.Length == 0) return;
        foreach (string entry in cell.Split(';'))
        {
            string[] kv = entry.Split('=');
            int[] n = [.. kv[0].Split('/').Select(int.Parse)];
            family[n[0] + 8 * n[1] + 64 * n[2]] = long.Parse(kv[1], CultureInfo.InvariantCulture);
        }
    }

    private static string L(long? v) => v is { } x ? BassPathProfile.I(x) : "";

    private static string Join(long[] a) => string.Join(";", a.Select(BassPathProfile.I));

    public string Status => Hold is not null ? "held" : FirstEntry is not null ? "entered" : "open";

    /// <summary>The region's largest cells, by sweeps.</summary>
    public (string Basin, long Sweeps)[] RegionBasins(int top) => [.. _region.Values.OrderByDescending(c => c.Count).Take(top).Select(c => (c.Label, c.Count))];

    internal string[] Row(HoldSpec spec)
    {
        // The row keeps the region's largest cells, enough to rank the five largest.
        string region = string.Join(";", RegionBasins(8).Select(c => $"{c.Basin}={c.Sweeps.ToString(CultureInfo.InvariantCulture)}"));
        string recent = "";
        if (_rings is not null && _ringOn)
        {
            Ring wide = _rings[^1];
            var runs = new List<string>();
            long begin = -1;
            for (long t = Math.Max(1, Sweeps - wide.Width + 1); t <= Sweeps + 1; t++)
            {
                bool on = t <= Sweeps && wide.In[t % wide.Width];
                if (on && begin < 0) begin = t;
                if (on || begin < 0) continue;
                runs.Add($"{begin.ToString(CultureInfo.InvariantCulture)}-{(t - 1).ToString(CultureInfo.InvariantCulture)}");
                begin = -1;
            }
            recent = string.Join(";", runs);
        }
        string open = _visit is { } v ? $"{v.First.ToString(CultureInfo.InvariantCulture)};{v.Last.ToString(CultureInfo.InvariantCulture)};{v.In.ToString(CultureInfo.InvariantCulture)};{v.Entries.ToString(CultureInfo.InvariantCulture)};{v.From}" : "";
        return
        [
            spec.Start, spec.Arm, BassPathProfile.I(spec.Chain), BassPathProfile.I(Sweeps), Status, BassPathProfile.I(Segment), Seconds.ToString("F1", CultureInfo.InvariantCulture), BassPathProfile.I(Accepted), Structure, BassPathProfile.F(LogTau),
            L(FirstExit), L(Departure), BassPathProfile.I(OutSince), L(FirstEntry), EntryFrom, BassPathProfile.I(StaySweeps), BassPathProfile.I(StaySeen), BassPathProfile.I(Entries), BassPathProfile.I(Visits), L(Touches), HoldFrom,
            L(Held[0]), L(Held[1]), L(Held[2]), L(Held[3]), L(Held[4]), L(Held[5]), L(LastBothPre), L(Return), BassPathProfile.I(PreSweeps), BassPathProfile.I(PostSweeps), BassPathProfile.I(PostNotBoth), BassPathProfile.I(PostTarget),
            Map(PreFamily), Join(PreLargest), BassPathProfile.I(RegionSweeps), Map(RegionFamily), region, Map(PostFamily), Join(PostLargest), BassPathProfile.I(Restarts), recent, open,
        ];
    }

    internal static HoldTally Parse(string[] r)
    {
        static long? O(string v) => v.Length == 0 ? null : long.Parse(v, CultureInfo.InvariantCulture);
        static long N(string v) => v.Length == 0 ? 0L : long.Parse(v, CultureInfo.InvariantCulture);
        static long[] A(string v) => v.Length == 0 ? [] : [.. v.Split(';').Select(x => long.Parse(x, CultureInfo.InvariantCulture))];
        var t = new HoldTally
        {
            Sweeps = N(r[3]), Segment = (int)N(r[5]), Seconds = double.Parse(r[6], CultureInfo.InvariantCulture), Accepted = N(r[7]), Structure = r[8], LogTau = r[9].Length == 0 ? double.NaN : double.Parse(r[9], CultureInfo.InvariantCulture),
            FirstExit = O(r[10]), Departure = O(r[11]), OutSince = N(r[12]), FirstEntry = O(r[13]), EntryFrom = r[14], StaySweeps = N(r[15]), StaySeen = N(r[16]), Entries = N(r[17]), Visits = N(r[18]), Touches = O(r[19]), HoldFrom = r[20],
            Held = [O(r[21]), O(r[22]), O(r[23]), O(r[24]), O(r[25]), O(r[26])], LastBothPre = O(r[27]), Return = O(r[28]), PreSweeps = N(r[29]), PostSweeps = N(r[30]), PostNotBoth = N(r[31]), PostTarget = N(r[32]),
            PreLargest = A(r[34]), RegionSweeps = N(r[35]), PostLargest = A(r[39]), Restarts = (int)N(r[40]), _recent = r[41], _openVisit = r[42],
        };
        ReadMap(r[33], t.PreFamily);
        ReadMap(r[36], t.RegionFamily);
        ReadMap(r[38], t.PostFamily);
        if (r[37].Length > 0)
            foreach (string entry in r[37].Split(';'))
            {
                int eq = entry.LastIndexOf('=');
                t._region[Hash(BasinSeries.CodeOf(entry[..eq]))] = new RegionCell(entry[..eq]) { Count = long.Parse(entry[(eq + 1)..], CultureInfo.InvariantCulture) };
            }
        return t;
    }

    private static long Hash(ReadOnlySpan<int> code)
    {
        long h = unchecked((long)14695981039346656037UL);
        foreach (int c in code) h = unchecked((h ^ c) * 1099511628211L);
        return h;
    }

    /// <summary>Starts reading a stream at <paramref name="terms"/>: a new chain's start, or a stored row's state.</summary>
    internal void Begin(HoldSpec spec, IReadOnlyList<BassTerm> terms)
    {
        _spec = spec;
        _startCode = BassWaitingTime.Code(spec.Terms);
        if (spec.Windows.Length * spec.Levels.Length != Held.Length) throw new ArgumentException("A tally holds three windows at two levels.", nameof(spec));
        if (PreLargest.Length != spec.Largest.Length) PreLargest = new long[spec.Largest.Length];
        if (PostLargest.Length != spec.Largest.Length) PostLargest = new long[spec.Largest.Length];
        if (_blockLargest.Length != spec.Largest.Length) _blockLargest = new long[spec.Largest.Length];
        if (_rings is null && Held.Any(h => h is null))
        {
            _rings = [.. spec.Windows.Select(w => new Ring(w))];
            if (_recent.Length > 0)
            {
                _ringOn = true;
                foreach (string run in _recent.Split(';'))
                {
                    string[] ends = run.Split('-');
                    long from = long.Parse(ends[0], CultureInfo.InvariantCulture), to = long.Parse(ends[1], CultureInfo.InvariantCulture);
                    foreach (Ring ring in _rings)
                        for (long t = Math.Max(from, Sweeps - ring.Width + 1); t <= to; t++)
                        {
                            ring.In[t % ring.Width] = true;
                            ring.Count++;
                        }
                }
            }
        }
        if (_openVisit.Length > 0)
        {
            string[] v = _openVisit.Split(';');
            _visit = new Visit(long.Parse(v[0], CultureInfo.InvariantCulture), v[4]) { Last = long.Parse(v[1], CultureInfo.InvariantCulture), In = long.Parse(v[2], CultureInfo.InvariantCulture), Entries = long.Parse(v[3], CultureInfo.InvariantCulture) };
            _openVisit = "";
        }
        _recent = "";
        _current = terms;
        _before = null;
        _cell = null;
        Classify(terms);
    }

    private void Classify(IReadOnlyList<BassTerm> terms)
    {
        int m = terms.Count;
        if (m > _code.Length) _code = new int[2 * m];
        int[] code = _code, family = _spec.Family;
        int n0 = 0, n1 = 0, n2 = 0;
        for (int i = 0; i < m; i++)
        {
            int mask = BassWaitingTime.Mask(terms[i]);
            code[i] = mask;
            if (mask == family[0]) n0++;
            else if (mask == family[1]) n1++;
            else if (mask == family[2]) n2++;
        }
        Array.Sort(code, 0, m);
        _m = m;
        _family = Math.Min(n0, 7) + 8 * Math.Min(n1, 7) + 64 * Math.Min(n2, 7);
        ReadOnlySpan<int> held = code.AsSpan(0, m);
        _inTarget = held.SequenceEqual(_spec.Target);
        _inStart = held.SequenceEqual(_startCode);
        _both = true;
        foreach (int mask in _spec.Returns)
            if (held.IndexOf(mask) < 0)
            {
                _both = false;
                break;
            }
        _largest = -1;
        for (int b = 0; b < _spec.Largest.Length; b++)
            if (held.SequenceEqual(_spec.Largest[b]))
            {
                _largest = b;
                break;
            }
    }

    private string From() => _before is null ? "" : _spec.Basin(_before);

    /// <summary>One sweep: the state it ended on is the same, or <paramref name="changed"/>.</summary>
    internal void Tick(IReadOnlyList<BassTerm>? changed)
    {
        long t = ++Sweeps;
        bool wasIn = _inTarget;
        if (changed is not null)
        {
            bool wasStart = _inStart;
            _before = _current;
            _current = changed;
            Accepted++;
            Classify(changed);
            _cell = null;
            if (wasStart && !_inStart)
            {
                OutSince = t;
                FirstExit ??= t;
            }
        }
        if (Departure is null && !_inStart && t - OutSince >= _spec.Gap) Departure = OutSince;
        if (_inTarget && !wasIn)
        {
            if (FirstEntry is null)
            {
                FirstEntry = t;
                EntryFrom = From();
            }
            if (Hold is null)
            {
                Entries++;
                _visit ??= new Visit(t, From());
                _visit.Entries++;
            }
        }
        if (FirstEntry is null)
        {
            PreSweeps++;
            PreFamily[_family]++;
            if (_largest >= 0) PreLargest[_largest]++;
            if (_both) LastBothPre = t;
            if (Departure is not null)
            {
                RegionSweeps++;
                RegionFamily[_family]++;
                (_cell ??= Cell()).Count++;
            }
        }
        else if (StaySeen < PassageTally.Window && t > FirstEntry)
        {
            StaySeen++;
            if (_inTarget) StaySweeps++;
        }
        if (_visit is { } visit)
        {
            if (_inTarget)
            {
                visit.In++;
                visit.Last = t;
            }
            else if (t - visit.Last >= _spec.Gap) CloseVisit(visit);
        }
        if (_rings is not null && (_ringOn || _inTarget))
        {
            _ringOn = true;
            Slide(t);
        }
        if (Hold is not null && t >= PostFrom) Past(t);
    }

    private RegionCell Cell()
    {
        long h = Hash(_code.AsSpan(0, _m));
        if (!_region.TryGetValue(h, out RegionCell? cell)) _region[h] = cell = new RegionCell(_spec.Basin(_current!));
        return cell;
    }

    private void CloseVisit(Visit v)
    {
        Visits++;
        _closed.Add((v.First, v.Last, v.From));
        PendingVisits.Add([_spec.Start, _spec.Arm, BassPathProfile.I(_spec.Chain), BassPathProfile.I(Visits), BassPathProfile.I(v.First), BassPathProfile.I(v.Last), BassPathProfile.I(v.In), BassPathProfile.I(v.Entries), v.From]);
        _visit = null;
    }

    private void Slide(long t)
    {
        Ring[] rings = _rings!;
        for (int w = 0; w < rings.Length; w++)
        {
            Ring ring = rings[w];
            ref bool slot = ref ring.In[t % ring.Width];
            if (slot != _inTarget)
            {
                ring.Count += _inTarget ? 1 : -1;
                slot = _inTarget;
            }
            if (t < ring.Width) continue;
            for (int c = 0; c < _spec.Levels.Length; c++)
            {
                int at = w * _spec.Levels.Length + c;
                if (Held[at] is not null || ring.Count < Math.Ceiling(_spec.Levels[c] * ring.Width - 1e-9)) continue;
                Held[at] = t - ring.Width + 1;
                if (at == Held.Length - 1) Holds(Held[at]!.Value);
            }
        }
        if (Held.All(h => h is not null)) _rings = null;
    }

    private void Holds(long start)
    {
        long after = _closed.Count(c => c.Last >= start);
        Touches = Visits - after;
        HoldFrom = _closed.Where(c => c.Last >= start).Select(c => c.From).FirstOrDefault(_visit?.From ?? "");
        if (_visit is { } open) CloseVisit(open);
        _closed.Clear();
    }

    private void Past(long t)
    {
        PostSweeps++;
        PostFamily[_family]++;
        _blockSweeps++;
        if (_largest >= 0)
        {
            PostLargest[_largest]++;
            _blockLargest[_largest]++;
        }
        if (!_both)
        {
            PostNotBoth++;
            _blockNotBoth++;
        }
        if (_inTarget)
        {
            PostTarget++;
            _blockTarget++;
        }
        if (_both) Return ??= t;
        if (t % _spec.Block != 0) return;
        PendingBlocks.Add([_spec.Start, _spec.Arm, BassPathProfile.I(_spec.Chain), BassPathProfile.I(t), BassPathProfile.I(_blockSweeps), BassPathProfile.I(_blockNotBoth), BassPathProfile.I(_blockTarget), Join(_blockLargest)]);
        _blockSweeps = _blockNotBoth = _blockTarget = 0;
        Array.Clear(_blockLargest);
    }

    /// <summary>The chain's next segment, from its stored structure and τ or its start, on streams of its own name and segment.</summary>
    internal void Advance(HoldSpec spec, long limit)
    {
        var clock = Stopwatch.StartNew();
        int[] seeds = SeedTree.Derive(BassWaitingTime.SeedOf($"basin-hold|{spec.Name}|{Segment}"), 2);
        BassTerm[] from = Structure.Length == 0 ? spec.Terms : BassWaitingTime.ParseStructure(Structure);
        double logTau = double.IsNaN(LogTau) ? spec.LogTau : LogTau;
        var chain = new BassChain(spec.Model, from, logTau, new Xoshiro256PlusPlus(seeds[0]), new Xoshiro256PlusPlus(seeds[1]));
        IReadOnlyList<BassTerm> last = chain.Current.Terms;
        Begin(spec, last);
        long begin = Sweeps;
        while (Sweeps < limit)
        {
            IReadOnlyList<BassTerm> terms = chain.Step().Terms;
            if (ReferenceEquals(terms, last)) Tick(null);
            else
            {
                last = terms;
                Tick(terms);
            }
        }
        Structure = BassMixing.Key(chain.Current.Terms);
        LogTau = chain.Current.LogTau;
        if (chain.Moves.Attempts != Sweeps - begin) throw new InvalidOperationException($"{chain.Moves.Attempts} attempts in {Sweeps - begin} sweeps: a sweep is not one attempt.");
        Segment++;
        Seconds += clock.Elapsed.TotalSeconds;
    }
}

/// <summary>The hold chains' book: each chain's tally, a row at every checkpoint, its visits and blocks, and the core-hours spent.</summary>
internal sealed class HoldBook
{
    /// <summary>The brief's budget for Parts 4 and 5.</summary>
    internal const double BudgetHours = 24.0;

    internal const int Threads = 4;

    /// <summary>A chain's structure and τ are saved this often.</summary>
    internal const long Checkpoint = 10_000_000;

    private readonly PathCsv _chains = new("chains", HoldTally.Header, 4, "basin-hold");
    private readonly PathCsv _visits = new("visits", HoldTally.VisitHeader, 4, "basin-hold");
    private readonly PathCsv _blocks = new("blocks", HoldTally.BlockHeader, 4, "basin-hold");
    private readonly ConcurrentDictionary<string, HoldTally> _live = new();

    internal HoldBook()
    {
        var latest = new Dictionary<string, string[]>();
        foreach (string[] r in _chains.Rows()) latest[$"{r[0]}|{r[1]}|{r[2]}"] = r;
        foreach ((string name, string[] row) in latest)
        {
            HoldTally tally = HoldTally.Parse(row);
            tally.Restarts++;
            _live[name] = tally;
        }
    }

    public double CoreHours => _live.Values.Sum(t => t.Seconds) / 3600.0;

    /// <summary>Spent, with the segments a check may leave in flight at the rate so far, so the budget is a ceiling.</summary>
    public bool Exhausted => CoreHours + Threads * CoreHours / _live.Values.Sum(t => t.Sweeps) * Checkpoint >= BudgetHours;

    /// <summary>Every chain of <paramref name="specs"/> run to <paramref name="limit"/> sweeps, four at a time, a row at each checkpoint; false when the budget stopped it.</summary>
    public bool Advance(IReadOnlyList<HoldSpec> specs, long limit, ITestOutputHelper? log = null)
    {
        var clock = Stopwatch.StartNew();
        bool stopped = false;
        HoldSpec[] open = [.. specs.Where(s => !_live.TryGetValue(s.Name, out HoldTally? t) || t.Sweeps < limit)];
        Parallel.ForEach(Partitioner.Create(open, EnumerablePartitionerOptions.NoBuffering), new ParallelOptions { MaxDegreeOfParallelism = Threads }, spec =>
        {
            HoldTally tally = _live.GetOrAdd(spec.Name, _ => new HoldTally());
            while (tally.Sweeps < limit)
            {
                if (Exhausted)
                {
                    stopped = true;
                    return;
                }
                tally.Advance(spec, Math.Min(limit, (tally.Sweeps / Checkpoint + 1) * Checkpoint));
                foreach (string[] v in tally.PendingVisits) _visits.Append(v);
                foreach (string[] b in tally.PendingBlocks) _blocks.Append(b);
                tally.PendingVisits.Clear();
                tally.PendingBlocks.Clear();
                _chains.Append(tally.Row(spec));
            }
        });
        log?.WriteLine($"{open.Length} chains to {limit:N0}: {clock.Elapsed.TotalMinutes:F1} min; {CoreHours:F2} core-hours spent{(stopped ? "; the budget stopped it" : "")}.");
        return !stopped;
    }
}
