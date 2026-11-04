using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>What one start's hitting time reads: ln h = ln E_S τ_A; ln of the escape probability P_S(τ_A &lt; τ_S⁺), r, so that
/// π(S)R_eff = 1/r in any normalization; ln M/π(S) from the voltage's sum, and ln of the forward reduction's mean excursion, which
/// the identity E_S τ_A = R_eff · M says is the same number.</summary>
internal readonly record struct HittingReading(double LnHit, double LnEscape, double LnMassRatio, double LnExcursion)
{
    /// <summary>ln π(S)R_eff: the wait a single edge of the whole network's conductance would read.</summary>
    public double LnWaitEffective => -LnEscape;

    /// <summary>ln(R_eff · M), read from the voltage.</summary>
    public double LnHitFromVoltage => LnMassRatio - LnEscape;
}

/// <summary>
/// Mean hitting times, effective resistances and voltages of a reversible chain, by state reduction with no subtraction
/// (Grassmann, Taksar and Heyman, as <see cref="EliminatedChain"/>): each pivot is the sum of the rates left in its row, the rate
/// into the absorbed set among them, so every factor keeps K's relative accuracy however long the wait. Three stages share the work:
/// the states that are nobody's start or target go once; a target's other states go once per target; each start's fellow starts
/// go per start, leaving the censored chain on the start and the target, whose rate is the escape probability.
/// </summary>
internal sealed class HittingReduction
{
    private const double PivotFloor = 1e-290;

    private readonly ReversibleChain _chain;
    private readonly int _threads;
    private readonly Stage _free;
    private readonly double[][] _block;      // the kept states' rates among themselves, by position in _kept
    private readonly double[] _time;         // the kept states' mean sweeps per visit
    private readonly int[] _kept;
    private readonly Dictionary<int, int> _position;

    // One stage's factors: States[..Count] went in order; Upper[k] holds k's rates to States[k + 1 ..].
    internal sealed record Stage(int[] States, int Count, double[][] Upper, double[] Pivot);

    /// <summary>The first stage on the irreducible <paramref name="chain"/>: everything but <paramref name="kept"/> eliminated, ascending π.</summary>
    public HittingReduction(ReversibleChain chain, IEnumerable<int> kept, int threads = 1)
    {
        _chain = chain;
        _threads = threads;
        int n = chain.Size;
        var keeping = new bool[n];
        foreach (int k in kept) keeping[k] = true;
        int[] order = [.. Enumerable.Range(0, n).OrderBy(i => keeping[i]).ThenBy(i => chain.LogPi[i]).ThenBy(i => i)];
        int free = n - keeping.Count(b => b);
        var position = new int[n];
        for (int k = 0; k < n; k++) position[order[k]] = k;
        var rows = new double[n][];
        for (int k = 0; k < n; k++)
        {
            var row = new double[n];
            int i = order[k];
            for (int e = chain.Starts[i]; e < chain.Starts[i + 1]; e++) row[position[chain.To[e]]] = chain.Rate[e];
            rows[k] = row;
        }
        var time = new double[n];
        Array.Fill(time, 1.0);
        _free = Eliminate(order, rows, null, time, free, threads);
        _kept = order[free..];
        _position = [];
        for (int a = 0; a < _kept.Length; a++) _position.Add(_kept[a], a);
        _block = new double[_kept.Length][];
        for (int a = 0; a < _kept.Length; a++) _block[a] = rows[free + a][free..];
        _time = time[free..];
    }

    // Eliminates States[..count]; each row's rate into the absorbed set is root[i], null when none is absorbed yet.
    private static Stage Eliminate(int[] states, double[][] rows, double[]? root, double[] time, int count, int threads)
    {
        int n = states.Length;
        var upper = new double[count][];
        var pivot = new double[count];
        var options = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) };
        for (int k = 0; k < count; k++)
        {
            double[] rk = rows[k];
            double d = root is null ? 0.0 : root[k];
            for (int j = k + 1; j < n; j++) d += rk[j];
            if (!(d > PivotFloor)) throw new InvalidOperationException($"A pivot of {d:E2} at state {states[k]}: past the doubles' range, or the chain is not irreducible.");
            pivot[k] = d;
            int at = k;
            void Row(int i)
            {
                double[] ri = rows[i];
                double f = ri[at] / d;
                if (f == 0.0) return;
                int j = at + 1, width = Vector<double>.Count;
                var scale = new Vector<double>(f);
                for (; j <= n - width; j += width)
                    (new Vector<double>(ri, j) + scale * new Vector<double>(rk, j)).CopyTo(ri, j);
                for (; j < n; j++) ri[j] += f * rk[j];
                time[i] += f * time[at];
                if (root is not null) root[i] += f * root[at];
            }
            if (threads > 1 && n - k > 128) Parallel.For(k + 1, n, options, Row);
            else for (int i = k + 1; i < n; i++) Row(i);
            upper[k] = rk[(k + 1)..];
            rows[k] = [];
        }
        return new Stage(states, count, upper, pivot);
    }

    /// <summary>The second stage: <paramref name="target"/> absorbed, and all kept states but <paramref name="starts"/> eliminated.</summary>
    public Target Toward(IReadOnlyCollection<int> target, IReadOnlyCollection<int> starts)
    {
        var absorbing = new HashSet<int>(target);
        var starting = new HashSet<int>(starts);
        starting.ExceptWith(absorbing);
        int[] remaining = [.. _kept.Where(s => !absorbing.Contains(s)).OrderBy(starting.Contains)];
        int count = remaining.Length - starting.Count;
        var rows = new double[remaining.Length][];
        var root = new double[remaining.Length];
        var time = new double[remaining.Length];
        for (int a = 0; a < remaining.Length; a++)
        {
            double[] from = _block[_position[remaining[a]]];
            var row = new double[remaining.Length];
            for (int b = 0; b < remaining.Length; b++) row[b] = from[_position[remaining[b]]];
            foreach (int c in absorbing) root[a] += from[_position[c]];
            rows[a] = row;
            time[a] = _time[_position[remaining[a]]];
        }
        Stage stage = Eliminate(remaining, rows, root, time, count, _threads);
        return new Target(this, stage, remaining[count..], [.. rows[count..].Select(r => r[count..])], root[count..], time[count..]);
    }

    /// <summary>A target's factors and its starts' remaining block.</summary>
    internal sealed class Target
    {
        private readonly HittingReduction _owner;
        private readonly Stage _stage;
        private readonly int[] _starts;
        private readonly double[][] _block;
        private readonly double[] _root;
        private readonly double[] _time;

        internal Target(HittingReduction owner, Stage stage, int[] starts, double[][] block, double[] root, double[] time)
        {
            _owner = owner;
            _stage = stage;
            _starts = starts;
            _block = block;
            _root = root;
            _time = time;
        }

        /// <summary>What <paramref name="start"/>, one of the starts given, reads, and its voltage by state: v(start) = 1, 0 on the
        /// target. Null when the chain cannot reach the target.</summary>
        public (HittingReading Reading, double[] Voltage)? Read(int start)
        {
            int[] states = [.. _starts.Where(s => s != start), start];
            var index = new Dictionary<int, int>(_starts.Length);
            for (int a = 0; a < _starts.Length; a++) index[_starts[a]] = a;
            int n = states.Length;
            var rows = new double[n][];
            var root = new double[n];
            var time = new double[n];
            for (int a = 0; a < n; a++)
            {
                int i = index[states[a]];
                rows[a] = [.. states.Select(s => _block[i][index[s]])];
                root[a] = _root[i];
                time[a] = _time[i];
            }
            Stage last = Eliminate(states, rows, root, time, n - 1, 1);
            double r = root[n - 1];
            if (!(r > 0.0)) return null;

            ReversibleChain chain = _owner._chain;
            var v = new double[chain.Size];
            v[start] = 1.0;
            Back(last, v);
            Back(_stage, v);
            Back(_owner._free, v);
            double reference = chain.LogPi[start], ratio = 0.0;
            for (int x = 0; x < v.Length; x++)
                if (v[x] > 0.0) ratio += Math.Exp(chain.LogPi[x] - reference) * v[x];
            double lnEscape = Math.Log(r), lnExcursion = Math.Log(time[n - 1]);
            return (new HittingReading(lnExcursion - lnEscape, lnEscape, Math.Log(ratio), lnExcursion), v);
        }

        private static void Back(Stage stage, double[] v)
        {
            for (int k = stage.Count - 1; k >= 0; k--)
            {
                double[] u = stage.Upper[k];
                double sum = 0.0;
                for (int m = 0; m < u.Length; m++) sum += u[m] * v[stage.States[k + 1 + m]];
                v[stage.States[k]] = sum / stage.Pivot[k];
            }
        }
    }

    /// <summary>A plain solve of (I − K)h = 1 off <paramref name="target"/>, with subtraction and no pivoting, for the check: h by state,
    /// 0 on the target. It keeps six digits only while the wait stays well under e³⁰.</summary>
    public static double[] Dense(ReversibleChain chain, IReadOnlyCollection<int> target, int threads = 1)
    {
        var absorbing = new HashSet<int>(target);
        int[] states = [.. Enumerable.Range(0, chain.Size).Where(i => !absorbing.Contains(i))];
        var index = new int[chain.Size];
        Array.Fill(index, -1);
        for (int k = 0; k < states.Length; k++) index[states[k]] = k;
        int n = states.Length;
        var a = new double[n][];
        var h = new double[n];
        Array.Fill(h, 1.0);
        for (int k = 0; k < n; k++)
        {
            var row = new double[n];
            int i = states[k];
            row[k] = 1.0 - chain.Hold[i];
            for (int e = chain.Starts[i]; e < chain.Starts[i + 1]; e++)
                if (index[chain.To[e]] >= 0) row[index[chain.To[e]]] -= chain.Rate[e];
            a[k] = row;
        }
        var options = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) };
        for (int k = 0; k < n; k++)
        {
            double[] rk = a[k];
            double d = rk[k];
            int at = k;
            void Row(int i)
            {
                double[] ri = a[i];
                double f = ri[at] / d;
                if (f == 0.0) return;
                int j = at + 1, width = Vector<double>.Count;
                var scale = new Vector<double>(f);
                for (; j <= n - width; j += width)
                    (new Vector<double>(ri, j) - scale * new Vector<double>(rk, j)).CopyTo(ri, j);
                for (; j < n; j++) ri[j] -= f * rk[j];
                h[i] -= f * h[at];
            }
            if (threads > 1 && n - k > 128) Parallel.For(k + 1, n, options, Row);
            else for (int i = k + 1; i < n; i++) Row(i);
        }
        for (int k = n - 1; k >= 0; k--)
        {
            double sum = h[k];
            double[] rk = a[k];
            for (int j = k + 1; j < n; j++) sum -= rk[j] * h[j];
            h[k] = sum / rk[k];
        }
        var byState = new double[chain.Size];
        for (int k = 0; k < n; k++) byState[states[k]] = h[k];
        return byState;
    }
}

/// <summary>First passages simulated from the enumerated kernel: a hold's geometric length from one draw, then the edge left by.</summary>
internal sealed class KernelWalker(ReversibleChain chain)
{
    private readonly double[] _logStay = [.. chain.Hold.Select((_, i) =>
    {
        double leaving = 0.0;
        for (int e = chain.Starts[i]; e < chain.Starts[i + 1]; e++) leaving += chain.Rate[e];
        return leaving >= 1.0 ? double.NegativeInfinity : Math.Log(1.0 - leaving);
    })];

    public long Passage(int start, bool[] target, Xoshiro256PlusPlus rng)
    {
        long t = 0;
        int x = start;
        while (!target[x])
        {
            double stay = _logStay[x];
            t += double.IsNegativeInfinity(stay) ? 1L : Math.Max(1L, (long)Math.Ceiling(Math.Log(1.0 - rng.NextDouble()) / stay));
            double leaving = 0.0;
            for (int e = chain.Starts[x]; e < chain.Starts[x + 1]; e++) leaving += chain.Rate[e];
            double u = rng.NextDouble() * leaving;
            int to = chain.To[chain.Starts[x + 1] - 1];
            for (int e = chain.Starts[x]; e < chain.Starts[x + 1]; e++)
            {
                u -= chain.Rate[e];
                if (u <= 0.0)
                {
                    to = chain.To[e];
                    break;
                }
            }
            x = to;
        }
        return t;
    }
}

/// <summary>A route's or a route part's series read from its start: ln W = ln Σ_e exp(π(S) − ln c_e), and with the class flow on additions
/// as the path profile takes it; ln of the slide's own Σ_e 1/c_e; the lowest π, the highest M.</summary>
internal readonly record struct RouteSeries(double PiFrom, double LnWait, double LnWaitClass, double LnSlideResistance, double Low, int HighestM, int Hops, int SlideHops)
{
    public double Depth => PiFrom - Low;
}

/// <summary>The harness's shared readings (the waiting-time brief; gap G1): a basin's code, a route's series and its cut at a basin, the
/// resistance of a lattice, and a stored structure's terms.</summary>
internal static class BassWaitingTime
{
    private static readonly Regex Factor = new(@"^x(\d+)([+-])(\d+)$", RegexOptions.Compiled);

    /// <summary>A continuous structure's basin as sorted masks of its terms' variables: equal exactly when
    /// <see cref="BassMixing.Fixture.Basin"/> is, with no string built.</summary>
    internal static int[] Code(IEnumerable<BassTerm> terms)
    {
        int[] code = [.. terms.Select(Mask)];
        Array.Sort(code);
        return code;
    }

    internal static int Mask(BassTerm term)
    {
        int mask = 0;
        foreach (HingeFactor h in term.Continuous) mask |= 1 << h.Variable;
        return mask;
    }

    internal static BassTerm ParseTerm(string key)
    {
        HingeFactor[] hinges = [.. key.Split('·').Select(f =>
        {
            Match m = Factor.Match(f);
            if (!m.Success) throw new FormatException($"A term's factor reads {f}.");
            return new HingeFactor(int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture), m.Groups[2].Value == "+" ? 1 : -1, int.Parse(m.Groups[3].Value, CultureInfo.InvariantCulture));
        })];
        return new BassTerm(hinges, null, null);
    }

    /// <summary>The structure <see cref="BassMixing.Key(IEnumerable{BassTerm})"/> wrote.</summary>
    internal static BassTerm[] ParseStructure(string key) => [.. key.Split(" + ").Select(ParseTerm)];

    /// <summary>The route's hops from <paramref name="forward"/>'s end, as the path profile reads them.</summary>
    internal static PathHop[] Directed(PathRoute route, bool forward) => forward ? route.Hops : [.. Enumerable.Reverse(route.Hops).Select(h => h.Reversed)];

    /// <summary>The basin code of each state along the route, S's first, in the route's own order: a slide changes knots and signs, never a
    /// term's variables, so only the lattice's hops move it.</summary>
    internal static int[][] CodesAlong(PathRoute route)
    {
        var codes = new int[route.Hops.Length + 1][];
        int lattice = 0;
        codes[0] = Code(route.Lattice.Terms(route.Masks[0]));
        for (int i = 0; i < route.Hops.Length; i++)
            codes[i + 1] = route.Hops[i].Kind == PathHopKind.Slide ? codes[i] : Code(route.Lattice.Terms(route.Masks[++lattice]));
        return codes;
    }

    /// <summary>The hops before the route's first state in <paramref name="target"/>'s basin, read from that end; null when the start is
    /// already in it.</summary>
    internal static PathHop[]? CutAt(PathRoute route, bool forward, int[] target)
    {
        int[][] codes = CodesAlong(route);
        if (!forward) Array.Reverse(codes);
        int cut = Array.FindIndex(codes, c => c.AsSpan().SequenceEqual(target));
        if (cut < 0) throw new InvalidOperationException("The route never reaches its end's basin.");
        return cut == 0 ? null : [.. Directed(route, forward).Take(cut)];
    }

    /// <summary>The series of <paramref name="hops"/>, read from their first; the class flow applies to an addition in this direction, as in
    /// <see cref="BassPathProfile.Read"/>.</summary>
    internal static RouteSeries Series(IReadOnlyList<PathHop> hops)
    {
        double piFrom = hops[0].PiFrom, wait = double.NegativeInfinity, cls = double.NegativeInfinity, slide = double.NegativeInfinity, low = piFrom;
        int highest = hops[0].MFrom, slides = 0;
        foreach (PathHop h in hops)
        {
            double lnFlow = h.LnFlow, lnClass = h.Kind == PathHopKind.Add && !double.IsNaN(h.LnFlowClass) ? Math.Max(lnFlow, h.LnFlowClass) : lnFlow;
            wait = LogMath.Add(wait, piFrom - lnFlow);
            cls = LogMath.Add(cls, piFrom - lnClass);
            low = Math.Min(low, h.PiTo);
            highest = Math.Max(highest, h.MTo);
            if (h.Kind != PathHopKind.Slide) continue;
            slides++;
            slide = LogMath.Add(slide, -lnFlow);
        }
        return new RouteSeries(piFrom, wait, cls, slide, low, highest, hops.Count, slides);
    }

    /// <summary>
    /// ln π(S)R over the lattice between one mask and a set of masks, the route's kept slide in series: the lattice's birth-and-death edges
    /// at the production flows, by state reduction. Null when no route joins them. By Rayleigh it lies between the whole network's and
    /// the route's.
    /// </summary>
    internal static double? LatticeWait(PathSurface surface, PathLattice lattice, int from, IReadOnlyCollection<int> target, double piStart, double lnSlideResistance, int threads = 1)
    {
        double[] pi = lattice.Held;
        var ids = new int[pi.Length];
        var kept = new List<int>();
        for (int mask = 0; mask < pi.Length; mask++)
        {
            ids[mask] = double.IsNegativeInfinity(pi[mask]) ? -1 : kept.Count;
            if (ids[mask] >= 0) kept.Add(mask);
        }
        var flows = new Dictionary<(int From, int To), double>();
        foreach (int lower in kept)
            for (int i = 0; i < lattice.Dim; i++)
            {
                int upper = lower | 1 << i;
                if (upper == lower || ids[upper] < 0) continue;
                int m = lattice.Core.Length + BitOperations.PopCount((uint)lower);
                double lnFlow = Math.Min(pi[lower] + surface.LogBirth(lattice.Terms(lower), m, lattice.Bits[i]), pi[upper] + surface.LogDeath(m + 1));
                flows[(ids[lower], ids[upper])] = lnFlow - pi[lower];
                flows[(ids[upper], ids[lower])] = lnFlow - pi[upper];
            }
        ReversibleChain whole = ReversibleChain.Of([.. kept.Select(m => pi[m])], [flows]);
        int[]? home = whole.Classes().FirstOrDefault(c => c.Contains(ids[from]));
        if (home is null) return null;
        var position = new Dictionary<int, int>(home.Length);
        for (int k = 0; k < home.Length; k++) position[home[k]] = k;
        int[] ends = [.. target.Where(m => ids[m] >= 0 && position.ContainsKey(ids[m])).Select(m => position[ids[m]])];
        if (ends.Length == 0) return null;
        int start = position[ids[from]];
        double lnLattice = double.NegativeInfinity;
        if (!ends.Contains(start))
        {
            ReversibleChain chain = whole.On(home);
            var solved = new HittingReduction(chain, [start, .. ends], threads).Toward(ends, [start]).Read(start);
            if (solved is null) return null;
            lnLattice = piStart - pi[from] - solved.Value.Reading.LnEscape;   // 1/r is the wait from the lattice's own start
        }
        return LogMath.Add(lnLattice, piStart + lnSlideResistance);
    }

    /// <summary>The masks of <paramref name="lattice"/> in the basin <paramref name="target"/>.</summary>
    internal static int[] MasksIn(PathLattice lattice, int[] target)
        => [.. Enumerable.Range(0, 1 << lattice.Dim).Where(m => !double.IsNegativeInfinity(lattice.Held[m]) && Code(lattice.Terms(m)).AsSpan().SequenceEqual(target))];

    /// <summary>A key's own seed, so a chain's stream depends on its name and its segment alone.</summary>
    internal static int SeedOf(string key)
    {
        uint h = 2166136261;
        foreach (char c in key)
        {
            h ^= c;
            h *= 16777619;
        }
        return (int)(h & 0x7fffffff);
    }
}

/// <summary>One row of the path profile's steps, read from S to T.</summary>
internal sealed record StepRow(string Kind, string Term, double PiFrom, double PiTo, double LnFlow, double LnFlowClass, double PhaseLow, double LnResistance, int MFrom, int MTo);

/// <summary>A route's series read back from <c>path-profile-steps.csv</c>, so the N = 500 and N = 125 readings take no new fit: full, and cut at
/// a basin by replaying the terms the rows add and delete.</summary>
internal static class StepSeries
{
    internal static StepRow[] Read(string n, string data, string pair, string exemplar, string alignment)
    {
        var csv = new PathCsv("steps", ["n"], 7);
        static double D(string v) => v.Length == 0 ? double.NaN : double.Parse(v, CultureInfo.InvariantCulture);
        return [.. csv.Rows().Where(r => r[0] == n && r[1] == data && r[2] == pair && r[3] == exemplar && r[4] == alignment && r[5] == nameof(PathReading.Held))
                   .OrderBy(r => int.Parse(r[6], CultureInfo.InvariantCulture))
                   .Select(r => new StepRow(r[7], r[8], D(r[15]), D(r[16]), D(r[19]), D(r[20]), D(r[23]), D(r[24]), int.Parse(r[9], CultureInfo.InvariantCulture), int.Parse(r[10], CultureInfo.InvariantCulture)))];
    }

    /// <summary>The series from S (<paramref name="forward"/>) or T, whose terms <paramref name="startTerms"/> are, up to the first state in
    /// <paramref name="target"/>'s basin (the whole route when null).</summary>
    internal static RouteSeries Walk(StepRow[] rows, bool forward, IEnumerable<BassTerm> startTerms, int[]? target)
    {
        List<int> held = [.. startTerms.Select(BassWaitingTime.Mask)];
        IEnumerable<StepRow> along = forward ? rows : Enumerable.Reverse(rows);
        double piFrom = forward ? rows[0].PiFrom : rows[^1].PiTo;
        double wait = double.NegativeInfinity, cls = double.NegativeInfinity, slide = double.NegativeInfinity, low = piFrom;
        int highest = forward ? rows[0].MFrom : rows[^1].MTo, hops = 0, slides = 0;
        bool In() => target is not null && held.Order().SequenceEqual(target);
        foreach (StepRow row in along)
        {
            if (In()) break;
            hops++;
            highest = Math.Max(highest, Math.Max(row.MFrom, row.MTo));
            low = Math.Min(low, row.PhaseLow);
            if (row.Kind == "slide")
            {
                slides++;
                wait = LogMath.Add(wait, piFrom + row.LnResistance);
                cls = LogMath.Add(cls, piFrom + row.LnResistance);
                slide = LogMath.Add(slide, row.LnResistance);
                continue;
            }
            int mask = BassWaitingTime.Mask(BassWaitingTime.ParseTerm(row.Term));
            bool adds = (row.Kind == "add") == forward;
            if (adds) held.Add(mask);
            else held.Remove(mask);
            double lnClass = adds && !double.IsNaN(row.LnFlowClass) ? Math.Max(row.LnFlow, row.LnFlowClass) : row.LnFlow;
            wait = LogMath.Add(wait, piFrom - row.LnFlow);
            cls = LogMath.Add(cls, piFrom - lnClass);
        }
        return new RouteSeries(piFrom, wait, cls, slide, low, highest, hops, slides);
    }
}

/// <summary>A chain to run: its name, its model and start, b, the start's and the target's basins and the target's exact structure.</summary>
internal sealed class PassageSpec
{
    public required string Part { get; init; }
    public required string Start { get; init; }
    public required string Arm { get; init; }
    public required string Rung { get; init; }
    public required int Chain { get; init; }
    public required double B { get; init; }
    public required PathSurface Surface { get; init; }
    public required BassModel Model { get; init; }
    public required BassTerm[] Terms { get; init; }
    public required double LogTau { get; init; }
    public required BassTerm[] TargetExact { get; init; }
    public required Func<IEnumerable<BassTerm>, string> Basin { get; init; }

    /// <summary>Stop at the exact structure, not the basin: a start already in the target's basin, or the check against an exact h.</summary>
    public bool ToExact { get; init; }

    /// <summary>Sweeps kept after first entry, to read whether the chain stays.</summary>
    public int StayWindow { get; init; } = PassageTally.Window;

    /// <summary>Whether to trace the last <see cref="PassageTally.Window"/> sweeps before first entry.</summary>
    public bool Trace { get; init; }

    public string Name => $"{Part}|{Start}|{Arm}|{Rung}|{Chain}";

    public int[] StartCode => BassWaitingTime.Code(Terms);

    public int[] TargetCode => BassWaitingTime.Code(TargetExact);
}

/// <summary>
/// One chain's cumulative record, as its rows store it, and its run. π here is the chain's own density in arm H, π_b; in arm S the structure's
/// π at the held τ, as the reading has it. Counts run over the sweeps' states, so the start's own sweep 0 is not among them. "Top" is the
/// highest structure of the start's basin before first entry, its dwell read since first reached: the fraction it holds of the sweeps since
/// then reads ln(M / π(top)) while v(top) ≈ 1, since a sweep's occupation before crossing is π(x)v(x)/M; of the sweeps in the basin since then,
/// ln(mass / π(top)). A departure's stretch restarts when a chain resumes outside the start's basin, so rows run before the column existed read
/// it late: part <c>3d</c>'s departure is not read.
/// </summary>
internal sealed class PassageTally
{
    /// <summary>The brief's window: first entry's stay, and the trace before it.</summary>
    internal const int Window = 10_000;

    /// <summary>The shortest stretch outside the start's basin that counts as a departure: one spurious term is gone in a few dozen sweeps.</summary>
    internal const int DepartureLength = 1_000;

    public long Sweeps { get; set; }
    public int Segment { get; set; }
    public long? FirstEntry { get; set; }
    public long StaySweeps { get; set; }
    public long StaySeen { get; set; }
    public long? FirstExit { get; set; }
    public string ExitTo { get; set; } = "";
    public long? ExactHit { get; set; }
    public double BestPi { get; set; } = double.NegativeInfinity;
    public long BestDwell { get; set; }
    public long BestSince { get; set; }
    public string BestKey { get; set; } = "";
    public long PreSweeps { get; set; }
    public long InStartPre { get; set; }
    public double TopPi { get; set; } = double.NegativeInfinity;
    public long TopDwell { get; set; }
    public long TopSinceAll { get; set; }
    public long TopSinceBasin { get; set; }
    public string TopKey { get; set; } = "";
    public long Accepted { get; set; }
    public double Seconds { get; set; }
    public double LogTau { get; set; } = double.NaN;
    public string Structure { get; set; } = "";
    public long? Departure { get; set; }
    public string EntryFrom { get; set; } = "";
    public double EntryPi { get; set; } = double.NaN;

    internal static readonly string[] Header =
    [
        "part", "start", "arm", "rung", "chain", "sweeps", "status", "b", "segment", "first_entry", "stay_sweeps", "stay_seen", "first_exit", "exit_to", "exact_hit",
        "best_pi", "best_dwell", "best_since", "best_key", "pre_sweeps", "in_start_pre", "top_pi", "top_dwell", "top_since_all", "top_since_basin", "top_key", "accepted", "seconds", "log_tau", "structure", "departure", "entry_from", "entry_pi",
    ];

    private static string L(long? v) => v is { } x ? BassPathProfile.I(x) : "";

    public bool Crossed(PassageSpec spec) => FirstEntry is not null && StaySeen >= spec.StayWindow;

    internal string[] Row(PassageSpec spec) =>
    [
        spec.Part, spec.Start, spec.Arm, spec.Rung, BassPathProfile.I(spec.Chain), BassPathProfile.I(Sweeps), Crossed(spec) ? "crossed" : "open", BassPathProfile.F(spec.B), BassPathProfile.I(Segment),
        L(FirstEntry), L(StaySweeps), L(StaySeen), L(FirstExit), ExitTo, L(ExactHit), BassPathProfile.F(BestPi), L(BestDwell), L(BestSince), BestKey, L(PreSweeps), L(InStartPre), BassPathProfile.F(TopPi),
        L(TopDwell), L(TopSinceAll), L(TopSinceBasin), TopKey, L(Accepted), Seconds.ToString("F1", CultureInfo.InvariantCulture), BassPathProfile.F(LogTau), Structure, L(Departure), EntryFrom, BassPathProfile.F(EntryPi),
    ];

    internal static PassageTally Parse(string[] r)
    {
        static long? O(string v) => v.Length == 0 ? null : long.Parse(v, CultureInfo.InvariantCulture);
        static long N(string v) => v.Length == 0 ? 0L : long.Parse(v, CultureInfo.InvariantCulture);
        static double D(string v) => v == "-inf" ? double.NegativeInfinity : double.Parse(v, CultureInfo.InvariantCulture);
        return new PassageTally
        {
            Sweeps = N(r[5]), Segment = (int)N(r[8]), FirstEntry = O(r[9]), StaySweeps = N(r[10]), StaySeen = N(r[11]), FirstExit = O(r[12]), ExitTo = r[13], ExactHit = O(r[14]),
            BestPi = D(r[15]), BestDwell = N(r[16]), BestSince = N(r[17]), BestKey = r[18], PreSweeps = N(r[19]), InStartPre = N(r[20]), TopPi = D(r[21]), TopDwell = N(r[22]), TopSinceAll = N(r[23]),
            TopSinceBasin = N(r[24]), TopKey = r[25], Accepted = N(r[26]), Seconds = D(r[27]), LogTau = D(r[28]), Structure = r[29], Departure = O(r[30]), EntryFrom = r[31], EntryPi = r[32].Length == 0 ? double.NaN : D(r[32]),
        };
    }

    /// <summary>
    /// The chain's next segment, from its stored structure or its start, on streams of its own name and segment: a chain resumed so is
    /// Markov, so first passage stays unbiased. It stops at <paramref name="limit"/> sweeps, or once first entry's stay is read.
    /// </summary>
    internal void Advance(PassageSpec spec, long limit, Action<(long Sweep, string Kind, string Out, string In, double Pi, string Basin, long Moves)>? trace)
    {
        var clock = Stopwatch.StartNew();
        int[] seeds = SeedTree.Derive(BassWaitingTime.SeedOf($"waiting-time|{spec.Name}|{Segment}"), 2);
        BassTerm[] from = Structure.Length == 0 ? spec.Terms : BassWaitingTime.ParseStructure(Structure);
        double logTau = double.IsNaN(LogTau) ? spec.LogTau : LogTau;
        var chain = new BassChain(spec.Model, from, logTau, new Xoshiro256PlusPlus(seeds[0]), new Xoshiro256PlusPlus(seeds[1]));
        if (spec.B != 1.0) chain.Moves.SetTemperature(new Temperature(spec.B, TemperingPath.Likelihood));
        bool held = spec.Model.Target.HeldTau is not null;
        int[] startCode = spec.StartCode, targetCode = spec.TargetCode;
        var exact = new HashSet<BassTerm>(spec.TargetExact);
        var code = new int[1024];
        long begin = Sweeps;
        var ring = new Queue<(long Sweep, IReadOnlyList<BassTerm> Terms, double Pi)>();

        bool inStart = false, inTarget = false, onBest = false, onTop = false;
        long outSince = Sweeps;
        IReadOnlyList<BassTerm>? before = null;
        double PiOf(BassState state) => held ? chain.Moves.CurrentLogDensity : spec.Surface.Score(state.Fit, profile: false).Pi(PathReading.Held);
        void See(BassState state, bool first)
        {
            IReadOnlyList<BassTerm> terms = state.Terms;
            int m = terms.Count;
            for (int i = 0; i < m; i++) code[i] = BassWaitingTime.Mask(terms[i]);
            Array.Sort(code, 0, m);
            bool wasIn = inStart;
            inStart = code.AsSpan(0, m).SequenceEqual(startCode);
            if (wasIn && !inStart) outSince = Sweeps;
            bool inBasin = code.AsSpan(0, m).SequenceEqual(targetCode), isExact = inBasin && m == exact.Count && terms.All(exact.Contains);
            inTarget = spec.ToExact ? isExact : inBasin;
            double pi = PiOf(state);
            if (spec.Trace && FirstEntry is null)
            {
                ring.Enqueue((Sweeps, terms, pi));
                while (ring.Peek().Sweep < Sweeps - Window) ring.Dequeue();
            }
            if (!inStart && FirstExit is null && !first)
            {
                FirstExit = Sweeps;
                ExitTo = spec.Basin(terms);
            }
            if (isExact && ExactHit is null && !first) ExactHit = Sweeps;
            if (inTarget && FirstEntry is null && !first)
            {
                FirstEntry = Sweeps;
                EntryFrom = before is null ? "" : spec.Basin(before);
                EntryPi = pi;
                if (trace is not null) Emit(ring, trace, spec);
            }
            if (pi > BestPi + 1e-9)
            {
                (BestPi, BestDwell, BestSince, BestKey) = (pi, 0, 0, BassMixing.Key(terms));
                onBest = true;
            }
            else onBest = Math.Abs(pi - BestPi) <= 1e-9;
            if (FirstEntry is null && inStart && pi > TopPi + 1e-9)
            {
                (TopPi, TopDwell, TopSinceAll, TopSinceBasin, TopKey) = (pi, 0, 0, 0, BassMixing.Key(terms));
                onTop = true;
            }
            else onTop = FirstEntry is null && inStart && Math.Abs(pi - TopPi) <= 1e-9;
        }

        IReadOnlyList<BassTerm> last = chain.Current.Terms;
        See(chain.Current, first: true);
        while (FirstEntry is null ? Sweeps < limit : StaySeen < spec.StayWindow)
        {
            BassState state = chain.Step();
            Sweeps++;
            if (!ReferenceEquals(state.Terms, last))
            {
                before = last;
                last = state.Terms;
                Accepted++;
                See(state, first: false);
            }
            if (FirstEntry is null)
            {
                if (!inStart && Departure is null && Sweeps - outSince >= DepartureLength) Departure = outSince;
                PreSweeps++;
                if (inStart) InStartPre++;
                TopSinceAll++;
                if (inStart) TopSinceBasin++;
                if (onTop) TopDwell++;
            }
            else if (FirstEntry != Sweeps)
            {
                StaySeen++;
                if (inTarget) StaySweeps++;
            }
            BestSince++;
            if (onBest) BestDwell++;
        }
        Structure = BassMixing.Key(chain.Current.Terms);
        LogTau = chain.Current.LogTau;
        if (chain.Moves.Attempts != Sweeps - begin) throw new InvalidOperationException($"{chain.Moves.Attempts} attempts in {Sweeps - begin} sweeps: a sweep is not one attempt.");
        Segment++;
        Seconds += clock.Elapsed.TotalSeconds;
    }

    // The moves that changed the basin in the window before entry, term by term; the knot moves between them as a count.
    private static void Emit(Queue<(long Sweep, IReadOnlyList<BassTerm> Terms, double Pi)> ring, Action<(long, string, string, string, double, string, long)> trace, PassageSpec spec)
    {
        (long Sweep, IReadOnlyList<BassTerm> Terms, double Pi)[] held = [.. ring];
        long moves = 0;
        for (int i = 1; i < held.Length; i++)
        {
            moves++;
            BassTerm[] gone = [.. held[i - 1].Terms.Where(t => !held[i].Terms.Contains(t))], come = [.. held[i].Terms.Where(t => !held[i - 1].Terms.Contains(t))];
            if (gone.Length == come.Length && gone.Length == 1 && BassWaitingTime.Mask(gone[0]) == BassWaitingTime.Mask(come[0])) continue;
            trace((held[i].Sweep, come.Length > gone.Length ? "add" : gone.Length > come.Length ? "delete" : "change", BassMixing.Key(gone), BassMixing.Key(come), held[i].Pi, spec.Basin(held[i].Terms), moves));
            moves = 0;
        }
    }
}

/// <summary>The chains' book: each chain's latest row, appended as it is reached and read back on a rerun, its traces, and the core-hours spent.</summary>
internal sealed class PassageBook
{
    /// <summary>The brief's budget for chains.</summary>
    internal const double BudgetHours = 36.0;

    internal const int Threads = 4;

    private static readonly string[] TraceHeader = ["part", "start", "arm", "rung", "chain", "sweep", "kind", "out", "in", "pi", "basin", "moves_since"];

    private readonly PathCsv _chains = new("chains", PassageTally.Header, 6, "waiting-time");
    private readonly PathCsv _traces = new("traces", TraceHeader, 6, "waiting-time");
    private readonly ConcurrentDictionary<string, PassageTally> _latest = new();

    internal PassageBook()
    {
        foreach (string[] r in _chains.Rows()) _latest[string.Join("|", r.Take(5))] = PassageTally.Parse(r);
    }

    /// <summary>Core-hours the chains have taken so far.</summary>
    public double CoreHours => _latest.Values.Sum(t => t.Seconds) / 3600.0;

    public bool Exhausted => CoreHours >= BudgetHours;

    public PassageTally? Latest(PassageSpec spec) => _latest.GetValueOrDefault(spec.Name);

    /// <summary>Every chain of <paramref name="specs"/> run to <paramref name="limit"/> sweeps or its crossing, four at a time, a row each; false
    /// when the budget stopped it before starting.</summary>
    public bool Advance(IReadOnlyList<PassageSpec> specs, long limit, ITestOutputHelper? log = null)
    {
        if (Exhausted)
        {
            log?.WriteLine($"Budget: {CoreHours:F1} core-hours spent; nothing further starts.");
            return false;
        }
        PassageSpec[] open = [.. specs.Where(s => Latest(s) is not { } t || (!t.Crossed(s) && t.Sweeps < limit))];
        var clock = Stopwatch.StartNew();
        Parallel.ForEach(Partitioner.Create(open, EnumerablePartitionerOptions.NoBuffering), new ParallelOptions { MaxDegreeOfParallelism = Threads }, spec =>
        {
            PassageTally tally = Latest(spec) ?? new PassageTally();
            tally.Advance(spec, limit, spec.Trace ? e => _traces.Append(spec.Part, spec.Start, spec.Arm, spec.Rung, BassPathProfile.I(spec.Chain), BassPathProfile.I(e.Sweep), e.Kind, e.Out, e.In, BassPathProfile.F(e.Pi), e.Basin, BassPathProfile.I(e.Moves)) : null);
            _latest[spec.Name] = tally;
            _chains.Append(tally.Row(spec));
        });
        log?.WriteLine($"{open.Length} chains to {limit:N0}: {clock.Elapsed.TotalMinutes:F1} min; {CoreHours:F2} core-hours spent.");
        return true;
    }
}

/// <summary>A lattice scored once and read at any b: its masks' evaluations, each edge's proposal, and each addition's class, kept as evaluations, so no b takes a new fit.</summary>
internal sealed class LadderLattice
{
    private readonly PathSurface _surface;
    private readonly BassEvaluation?[] _evaluations;
    private readonly Dictionary<(int Lower, int Bit), (double Birth, double Death)> _proposals = [];
    private readonly Dictionary<(int Lower, int Bit), BassEvaluation[]> _classes = [];

    internal LadderLattice(PathSurface surface, IReadOnlyList<BassTerm> core, IReadOnlyList<BassTerm> bits)
    {
        _surface = surface;
        Scored = new PathLattice(surface, core, bits);
        _evaluations = new BassEvaluation?[1 << Scored.Dim];
        Parallel.For(0, _evaluations.Length, new ParallelOptions { MaxDegreeOfParallelism = PassageBook.Threads }, mask =>
        {
            PathScore s = surface.Score(Scored.Terms(mask), profile: false);
            _evaluations[mask] = s.Admitted ? surface.Held.Target.Evaluate(s.Fit) : null;
        });
    }

    public PathLattice Scored { get; }

    /// <summary>ln π_b by mask; −∞ outside F.</summary>
    public double[] Pi(Temperature at)
    {
        var pi = new double[_evaluations.Length];
        for (int m = 0; m < pi.Length; m++) pi[m] = _evaluations[m] is { } e ? at.LogDensity(e) : double.NegativeInfinity;
        return pi;
    }

    public (double Birth, double Death) Proposal(int lower, int bit)
    {
        if (_proposals.TryGetValue((lower, bit), out (double, double) held)) return held;
        int m = Scored.Core.Length + BitOperations.PopCount((uint)lower);
        return _proposals[(lower, bit)] = (_surface.LogBirth(Scored.Terms(lower), m, Scored.Bits[bit]), _surface.LogDeath(m + 1));
    }

    /// <summary>The class flow at b for adding the bit at <paramref name="lower"/>, as <see cref="BassPathProfile.ClassFlow"/>.</summary>
    public double ClassFlow(int lower, int bit, double piLower, Temperature at)
    {
        (double birth, double death) = Proposal(lower, bit);
        if (!_classes.TryGetValue((lower, bit), out BassEvaluation[]? near))
        {
            BassTerm[] x = Scored.Terms(lower);
            BassDesign design = _surface.Design;
            if (design.Load(x) is { } refusal) throw new InvalidOperationException($"A lattice state is refused ({refusal}).");
            BassTerm t = Scored.Bits[bit];
            double logTermPrior = _surface.LogTermPrior(x) + _surface.Held.Prior.LogPrior(t);
            var found = new List<BassEvaluation>();
            foreach (BassTerm neighbour in BassPathProfile.Neighbours(_surface, t))
            {
                if (x.Contains(neighbour)) continue;
                DesignProposal proposal = design.ProposeBirth(neighbour);
                PathScore score = _surface.Score(BassStructureFit.Of(proposal, proposal.Admitted ? logTermPrior : double.NaN), profile: false);
                if (score.Admitted) found.Add(_surface.Held.Target.Evaluate(score.Fit));
            }
            _classes[(lower, bit)] = near = [.. found];
        }
        double sum = double.NegativeInfinity;
        double up = piLower + birth;
        foreach (BassEvaluation e in near) sum = LogMath.Add(sum, Math.Min(up, at.LogDensity(e) + death));
        return sum;
    }
}

/// <summary>What a pair's aligned route reads at one b, from S: the position its slide takes, the whole route's series and the series cut at the target's basin, exact and by class.</summary>
internal readonly record struct LadderReading(double B, string Position, double PiFrom, double LnWaitFull, double LnWaitClassFull, double LnWait, double LnWaitClass, double Low, double LowFull, int Hops, int HighestM);

/// <summary>
/// One pair's aligned route at any b (the ladder of the waiting-time brief): the alignment and the two slides as b = 1 has them, the slide's states re-read at b, and the lattice re-searched at b, the
/// position taken as <see cref="BassPathProfile.Build"/> takes it, by the higher lowest π.
/// </summary>
internal sealed class LadderPair
{
    private readonly (BassEvaluation State, double LnQ)[] _first, _last;
    private readonly BassEvaluation _start, _end;
    private readonly LadderLattice _after, _before;
    private readonly int _afterFrom, _afterTo, _beforeFrom, _beforeTo;
    private readonly int[] _target;
    private readonly int _startM;
    private readonly bool[] _afterIn, _beforeIn;

    internal LadderPair(PathSurface surface, BassTerm[] s, BassTerm[] t)
    {
        PathRoute route = BassPathProfile.Build(surface, s, t, false, PathReading.Held, [], priced: false) ?? throw new InvalidOperationException("No aligned route at b = 1.");
        PathAlignment a = route.Aligned;
        BassTerm[] bits = [.. a.OnlyFrom, .. a.OnlyTo], fromKnots = [.. a.Matched.Select(m => m.From)], toKnots = [.. a.Matched.Select(m => m.To)];
        _first = Slide(surface, s, a.Matched);
        _last = Slide(surface, t, [.. a.Matched.Select(m => (m.To, m.From))]);
        _start = Evaluate(surface, s);
        _end = Evaluate(surface, t);
        _startM = s.Length;
        _after = new LadderLattice(surface, [.. a.Common, .. toKnots], bits);
        _before = a.Matched.Length == 0 ? _after : new LadderLattice(surface, [.. a.Common, .. fromKnots], bits);
        _afterFrom = _after.Scored.Mask(a.OnlyFrom);
        _afterTo = _after.Scored.Mask(a.OnlyTo);
        _beforeFrom = _before.Scored.Mask(a.OnlyFrom);
        _beforeTo = _before.Scored.Mask(a.OnlyTo);
        _target = BassWaitingTime.Code(t);
        _afterIn = Inside(_after);
        _beforeIn = Inside(_before);
        Position = route.Position;
    }

    /// <summary>The position b = 1's search took.</summary>
    public string Position { get; }

    private bool[] Inside(LadderLattice lattice)
    {
        var inside = new bool[1 << lattice.Scored.Dim];
        for (int m = 0; m < inside.Length; m++) inside[m] = BassWaitingTime.Code(lattice.Scored.Terms(m)).AsSpan().SequenceEqual(_target);
        return inside;
    }

    private static BassEvaluation Evaluate(PathSurface surface, IReadOnlyList<BassTerm> terms)
    {
        PathScore score = surface.Score(terms, profile: false);
        if (!score.Admitted) throw new InvalidOperationException($"A slide's state is refused ({score.Refusal}).");
        return surface.Held.Target.Evaluate(score.Fit);
    }

    private static (BassEvaluation State, double LnQ)[] Slide(PathSurface surface, BassTerm[] start, (BassTerm From, BassTerm To)[] moves)
    {
        var visited = new List<BassTerm[]>();
        PathHop[] hops = BassPathProfile.Slide(surface, start, moves, PathReading.Held, out _, visited) ?? throw new InvalidOperationException("A slide b = 1 took no longer moves.");
        return [.. hops.Select((h, i) => (Evaluate(surface, visited[i + 1]), h.LnQForward))];
    }

    // The lattice's route at b, as (π from, π to, ln c_e, ln c_e by class, into the target's basin, M to).
    private static List<(double From, double To, double Flow, double Class, bool In, int M)>? Hops(LadderLattice lattice, Temperature at, int from, int to, bool[] inside)
    {
        double[] pi = lattice.Pi(at);
        int[]? masks = PathSearch.Route(pi, lattice.Scored.Dim, from, to, out _);
        if (masks is null) return null;
        var hops = new List<(double, double, double, double, bool, int)>(masks.Length);
        for (int k = 0; k + 1 < masks.Length; k++)
        {
            int u = masks[k], v = masks[k + 1], bit = BitOperations.TrailingZeroCount(u ^ v);
            int lower = Math.Min(u, v), upper = Math.Max(u, v);
            (double birth, double death) = lattice.Proposal(lower, bit);
            double flow = Math.Min(pi[lower] + birth, pi[upper] + death);
            double cls = v > u ? Math.Max(flow, lattice.ClassFlow(lower, bit, pi[lower], at)) : flow;
            hops.Add((pi[u], pi[v], flow, cls, inside[v], lattice.Scored.Core.Length + BitOperations.PopCount((uint)v)));
        }
        return hops;
    }

    /// <summary>The route read at <paramref name="b"/>; null when no route joins the ends there.</summary>
    public LadderReading? Read(double b)
    {
        var at = new Temperature(b, TemperingPath.Likelihood);
        double piS = at.LogDensity(_start), piT = at.LogDensity(_end);
        var afterHops = Hops(_after, at, _afterFrom, _afterTo, _afterIn);
        var beforeHops = Hops(_before, at, _beforeFrom, _beforeTo, _beforeIn);
        double SlideLow(IEnumerable<(BassEvaluation State, double LnQ)> slide, double from)
        {
            double low = from;
            foreach ((BassEvaluation state, _) in slide) low = Math.Min(low, at.LogDensity(state));
            return low;
        }
        double firstLow = double.NegativeInfinity, lastLow = double.NegativeInfinity;
        if (afterHops is not null) firstLow = Math.Min(SlideLow(_first, piS), afterHops.Count == 0 ? double.PositiveInfinity : afterHops.Min(h => Math.Min(h.From, h.To)));
        if (beforeHops is not null) lastLow = Math.Min(SlideLow(_last, piT), beforeHops.Count == 0 ? double.PositiveInfinity : beforeHops.Min(h => Math.Min(h.From, h.To)));
        if (afterHops is null && beforeHops is null) return null;
        bool first = firstLow >= lastLow;

        double wait = double.NegativeInfinity, cls = double.NegativeInfinity, fullWait = double.NegativeInfinity, fullClass = double.NegativeInfinity, low = piS, lowFull = piS;
        int hops = 0, highest = _startM;
        bool reached = false;
        void Add(double flow, double flowClass, double from, double to, bool intoTarget, int m)
        {
            fullWait = LogMath.Add(fullWait, piS - flow);
            fullClass = LogMath.Add(fullClass, piS - flowClass);
            lowFull = Math.Min(lowFull, Math.Min(from, to));
            if (reached) return;
            wait = LogMath.Add(wait, piS - flow);
            cls = LogMath.Add(cls, piS - flowClass);
            low = Math.Min(low, Math.Min(from, to));
            highest = Math.Max(highest, m);
            hops++;
            reached = intoTarget;
        }
        if (first)
        {
            double from = piS;
            foreach ((BassEvaluation state, double lnQ) in _first)
            {
                double to = at.LogDensity(state), flow = Math.Min(from, to) + lnQ;
                Add(flow, flow, from, to, false, _startM);
                from = to;
            }
            foreach (var h in afterHops!) Add(h.Flow, h.Class, h.From, h.To, h.In, h.M);
        }
        else
        {
            foreach (var h in beforeHops!) Add(h.Flow, h.Class, h.From, h.To, h.In, h.M);
            // The slide from T read backwards: every state of it is in T's basin, so the cut never holds it.
            double[] along = [piT, .. _last.Select(p => at.LogDensity(p.State))];
            for (int k = _last.Length - 1; k >= 0; k--)
            {
                double flow = Math.Min(along[k + 1], along[k]) + _last[k].LnQ;
                Add(flow, flow, along[k + 1], along[k], false, _startM);
            }
        }
        return new LadderReading(b, first ? "first" : "last", piS, fullWait, fullClass, wait, cls, low, lowFull, hops, highest);
    }
}
