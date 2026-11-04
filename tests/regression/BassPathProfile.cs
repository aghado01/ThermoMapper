using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>Which τ a structure's π is read at: held at unit information (SC-0's setting, the primary reading), or at each
/// structure's own best τ, its τ density included, closer to a chain whose τ moves every sweep.</summary>
internal enum PathReading { Held, Profiled }

/// <summary>A step along a route, read from S toward T: a term added, a term deleted, or a hinge's knot stepped.</summary>
internal enum PathHopKind { Add, Delete, Slide }

/// <summary>A structure's account at both τ readings, from <see cref="BassTarget"/>: refused outside F, or its fit, its log
/// structure prior, its log marginal at the held τ, and at its best τ the log marginal and τ's log density.</summary>
internal readonly record struct PathScore(BassRefusal? Refusal, BassStructureFit Fit, double LogStructurePrior, double HeldMarginal,
                                          double ProfiledMarginal, double ProfiledTauDensity, double BestLogTau)
{
    public bool Admitted => Refusal is null;

    public double Marginal(PathReading r) => r == PathReading.Held ? HeldMarginal : ProfiledMarginal;

    public double Prior(PathReading r) => r == PathReading.Held ? LogStructurePrior : LogStructurePrior + ProfiledTauDensity;

    /// <summary>ln π at the reading; −∞ outside F.</summary>
    public double Pi(PathReading r) => Admitted ? Marginal(r) + Prior(r) : double.NegativeInfinity;
}

/// <summary>
/// One surface's scorer (the path-profile brief): its data, the fixture's model for the proposals, a model holding τ at unit
/// information, and the sampled model the profile reads where the fixture samples τ. Structures load whole on a per-thread
/// design, so admission is production's: the cap, the support, repeated terms and the angle floor.
/// </summary>
internal sealed class PathSurface
{
    private const int GridPoints = 12;
    private const double GoldenTolerance = 1e-4;

    private readonly IReadOnlyList<IRjMove<BassState>> _palette;
    private readonly ThreadLocal<BassDesign> _design;
    private readonly ConcurrentDictionary<(string Key, int M), double> _select = new();
    private readonly object _edges = new();

    internal PathSurface(BassConfig config, BassData data)
    {
        Assert.Null(config.Sampling.LocalSteps);
        Assert.Null(config.Sampling.FactorMoves);
        Model = BassModel.Create(config, data);
        Held = config.Coefficients.Scale.IsHeld ? Model : BassModel.Create(config with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) }, data);
        Sampled = Model.Target.SamplesTau ? Model : null;
        HeldTau = Held.Target.HeldTau ?? throw new InvalidOperationException("The held model holds no τ.");
        _palette = new BassStructureMoves(Model, [], Math.Log(HeldTau)).Moves;
        _design = new ThreadLocal<BassDesign>(() => Held.NewDesign());
    }

    /// <summary>The fixture's model: the proposals' and the moves' weights.</summary>
    public BassModel Model { get; }

    /// <summary>The model with τ held at unit information.</summary>
    public BassModel Held { get; }

    /// <summary>The model the profile reads; null when the fixture holds τ, and the profile with it.</summary>
    public BassModel? Sampled { get; }

    public double HeldTau { get; }

    public int Observations => Model.Target.Observations;

    /// <summary>Superset edges read against SC-5's inequality, and the least slack among them.</summary>
    public long SupersetEdges { get; private set; }

    public double SupersetLeast { get; private set; } = double.PositiveInfinity;

    internal void CountSuperset(long edges, double least)
    {
        lock (_edges)
        {
            SupersetEdges += edges;
            SupersetLeast = Math.Min(SupersetLeast, least);
        }
    }

    /// <summary>This thread's design, for a proposal or a load.</summary>
    public BassDesign Design => _design.Value!;

    /// <summary>c from SC-5's form, ½ ln(1 + 1/τ), at u = ln τ: the g-form's Occam cost per column.</summary>
    public static double Occam(double logTau) => 0.5 * Math.Log(1.0 + Math.Exp(-logTau));

    public double HeldOccam => Occam(Math.Log(HeldTau));

    public double LogTermPrior(IEnumerable<BassTerm> terms) => BassMixing.LogTermPrior(Held, terms);

    /// <summary><paramref name="terms"/> loaded whole, then scored; the profile only when <paramref name="profile"/>.</summary>
    public PathScore Score(IReadOnlyList<BassTerm> terms, bool profile = true)
    {
        BassDesign design = Design;
        if (design.Load(terms) is { } refusal) return new PathScore(refusal, default, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN);
        return Score(BassStructureFit.Of(design, LogTermPrior(terms)), profile);
    }

    /// <summary>A fit scored through the target at the held τ and, when asked, at its best τ.</summary>
    public PathScore Score(BassStructureFit fit, bool profile = true)
    {
        BassTarget held = Held.Target;
        if (held.Refusal(fit) is { } refusal) return new PathScore(refusal, fit, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN);
        BassEvaluation at = held.Evaluate(fit);
        double marginal = at.LogMarginal, prior = at.LogPrior;
        if (!profile || Sampled is null) return new PathScore(null, fit, prior, marginal, double.NaN, double.NaN, double.NaN);
        double u = BestLogTau(fit);
        BassTarget sampled = Sampled.Target;
        return new PathScore(null, fit, prior, marginal, sampled.Evaluate(fit, u).LogMarginal, sampled.LogTauDensity(u), u);
    }

    /// <summary>The log marginal at u, for SC-5's check at a state's own best τ.</summary>
    public double MarginalAt(BassStructureFit fit, double logTau) => Sampled!.Target.Evaluate(fit, logTau).LogMarginal;

    // The best u for the fit: the coarse grid's best, then golden section around it. Log marginal plus τ's log density, the
    // structure prior being constant in u.
    private double BestLogTau(BassStructureFit fit)
    {
        BassTarget target = Sampled!.Target;
        double F(double u) => target.Evaluate(fit, u).LogMarginal + target.LogTauDensity(u);
        double centre = Math.Log(HeldTau), low = centre - 8.0, high = centre + 14.0, step = (high - low) / (GridPoints - 1);
        double best = low, top = double.NegativeInfinity;
        for (int i = 0; i < GridPoints; i++)
        {
            double u = low + i * step, value = F(u);
            if (value > top) (best, top) = (u, value);
        }
        double a = best - step, b = best + step, ratio = (Math.Sqrt(5.0) - 1.0) / 2.0;
        double c = b - ratio * (b - a), d = a + ratio * (b - a), fc = F(c), fd = F(d);
        while (b - a > GoldenTolerance)
        {
            if (fc >= fd)
            {
                (b, d, fd) = (d, c, fc);
                c = b - ratio * (b - a);
                fc = F(c);
            }
            else
            {
                (a, c, fc) = (c, d, fd);
                d = a + ratio * (b - a);
                fd = F(d);
            }
        }
        return fc >= fd ? c : d;
    }

    /// <summary>ln P(the sweep's move is <paramref name="key"/>) at M terms, from the palette's own weights.</summary>
    public double LogSelect(string key, int m) => _select.GetOrAdd((key, m), k =>
    {
        var state = new BassState(new BassTerm[k.M], default, 0.0);
        double total = _palette.Sum(move => Math.Max(0.0, move.Weight(state)));
        return Math.Log(BassMixing.Move(_palette, k.Key).Weight(state) / total);
    });

    /// <summary>ln q(x → x + t): the birth's selection at x and D01's draw of t under x's counts.</summary>
    public double LogBirth(IEnumerable<BassTerm> x, int m, BassTerm t) => LogSelect("birth", m) + Model.Proposal.LogBirth(t, BassMixing.Counts(Model, x));

    /// <summary>ln q(y → y − t) at M terms: the death's selection and its pick 1/M.</summary>
    public double LogDeath(int m) => LogSelect("death", m) - Math.Log(m);

    /// <summary>ln q of a redraw at M terms that lands on a given sign and knot of a factor of <paramref name="term"/> on
    /// <paramref name="variable"/>: the change's selection, 1/M, 1/K and 1/(2n_v), D01's change with the local steps off.</summary>
    public double LogRedraw(int m, BassTerm term, BassPart part, int variable)
        => LogSelect("change", m) - Math.Log(m) - Math.Log(BassMixing.Factors(term).Count()) - Math.Log(2.0 * BassMixing.Knots(Model, part, variable));
}

/// <summary>
/// L for a pair: every structure holding <see cref="Core"/> plus a subset of <see cref="Bits"/>, scored eagerly, bit i set when
/// Bits[i] is held. SC-5's inequality is read here on every superset edge at the held τ, and a violation stops the run.
/// </summary>
internal sealed class PathLattice
{
    internal PathLattice(PathSurface surface, IReadOnlyList<BassTerm> core, IReadOnlyList<BassTerm> bits)
    {
        Core = [.. core];
        Bits = [.. bits];
        Dim = Bits.Length;
        if (Dim > BassPathProfile.MaxDim) throw new ArgumentException($"|D| = {Dim} is past the enumeration's {BassPathProfile.MaxDim}; the route needs a beam.", nameof(bits));
        int size = 1 << Dim;
        Held = new double[size];
        Profiled = new double[size];
        var marginal = new double[size];
        var refusals = new BassRefusal?[size];
        Parallel.For(0, size, new ParallelOptions { MaxDegreeOfParallelism = BassPathProfile.Threads }, mask =>
        {
            PathScore s = surface.Score(Terms(mask));
            refusals[mask] = s.Refusal;
            Held[mask] = s.Pi(PathReading.Held);
            Profiled[mask] = surface.Sampled is null ? double.NaN : s.Pi(PathReading.Profiled);
            marginal[mask] = s.Admitted ? s.HeldMarginal : double.NaN;
        });
        foreach (BassRefusal? r in refusals)
            if (r is { } reason) Refused[reason] = Refused.GetValueOrDefault(reason) + 1;

        double c = surface.HeldOccam, least = double.PositiveInfinity;
        long edges = 0;
        for (int x = 0; x < size; x++)
        {
            if (refusals[x] is not null) continue;
            for (int i = 0; i < Dim; i++)
            {
                int y = x | (1 << i);
                if (y == x || refusals[y] is not null) continue;
                double slack = marginal[y] - (marginal[x] - c);
                edges++;
                if (slack < least) least = slack;
                if (slack < -1e-7) throw new InvalidOperationException($"SC-5's inequality fails by {-slack:E3} on a superset edge at the held τ: L(M + 1, R²′) = {marginal[y]:R}, L(M, R²) − c = {marginal[x] - c:R}, adding {BassMixing.Key([Bits[i]])}.");
            }
        }
        surface.CountSuperset(edges, least);
    }

    public BassTerm[] Core { get; }

    public BassTerm[] Bits { get; }

    public int Dim { get; }

    public double[] Held { get; }

    public double[] Profiled { get; }

    /// <summary>The states each check refused, by reason.</summary>
    public Dictionary<BassRefusal, int> Refused { get; } = [];

    public double[] Pi(PathReading r) => r == PathReading.Held ? Held : Profiled;

    public BassTerm[] Terms(int mask)
    {
        var terms = new List<BassTerm>(Core.Length + Dim);
        terms.AddRange(Core);
        for (int i = 0; i < Dim; i++)
            if ((mask >> i & 1) != 0) terms.Add(Bits[i]);
        return [.. terms];
    }

    public int Mask(IEnumerable<BassTerm> held)
    {
        int mask = 0;
        foreach (BassTerm t in held)
            if (Array.IndexOf(Bits, t) is var i && i >= 0) mask |= 1 << i;
        return mask;
    }
}

/// <summary>The maximin route through a cube of states, each state's ln π given: a widest-path search, and the references
/// that check it.</summary>
internal static class PathSearch
{
    /// <summary>The route's lowest π, the highest over every path from <paramref name="start"/> to <paramref name="end"/>
    /// along the steps <paramref name="allowed"/> admits; −∞ when none connects them.</summary>
    public static double Widest(double[] pi, int dim, int start, int end, Func<int, int, bool>? allowed = null)
    {
        var best = new double[1 << dim];
        Array.Fill(best, double.NegativeInfinity);
        var queue = new PriorityQueue<int, double>();
        best[start] = pi[start];
        queue.Enqueue(start, -pi[start]);
        while (queue.TryDequeue(out int u, out double priority))
        {
            double bottleneck = -priority;
            if (bottleneck < best[u]) continue;
            if (u == end) return bottleneck;
            for (int i = 0; i < dim; i++)
            {
                int v = u ^ (1 << i);
                if (allowed is not null && !allowed(u, v)) continue;
                double through = Math.Min(bottleneck, pi[v]);
                if (through > best[v])
                {
                    best[v] = through;
                    queue.Enqueue(v, -through);
                }
            }
        }
        return double.NegativeInfinity;
    }

    /// <summary>
    /// A maximin route as its states, start and end included, or null when none exists. Among routes reaching the maximin it
    /// takes the fewest steps, then the least Σ_e exp(maximin − min(π(x), π(y))): the widest path's value does not decide the
    /// route past its lowest state, so this reads π and the step count alone, and a route's proposal costs do not choose it.
    /// </summary>
    public static int[]? Route(double[] pi, int dim, int start, int end, out double maximin, Func<int, int, bool>? allowed = null)
    {
        double value = Widest(pi, dim, start, end, allowed);
        maximin = value;
        if (double.IsNegativeInfinity(value)) return null;
        var steps = new int[1 << dim];
        var resistance = new double[1 << dim];
        var from = new int[1 << dim];
        Array.Fill(steps, int.MaxValue);
        Array.Fill(from, -1);
        var queue = new PriorityQueue<int, (int Steps, double Resistance)>();
        steps[start] = 0;
        queue.Enqueue(start, (0, 0.0));
        while (queue.TryDequeue(out int u, out (int Steps, double Resistance) at))
        {
            if (at.Steps > steps[u] || (at.Steps == steps[u] && at.Resistance > resistance[u])) continue;
            if (u == end) break;
            for (int i = 0; i < dim; i++)
            {
                int v = u ^ (1 << i);
                if (pi[v] < value || (allowed is not null && !allowed(u, v))) continue;
                double through = at.Resistance + Math.Exp(value - Math.Min(pi[u], pi[v]));
                if (at.Steps + 1 < steps[v] || (at.Steps + 1 == steps[v] && through < resistance[v]))
                {
                    steps[v] = at.Steps + 1;
                    resistance[v] = through;
                    from[v] = u;
                    queue.Enqueue(v, (at.Steps + 1, through));
                }
            }
        }
        var route = new List<int> { end };
        while (route[^1] != start) route.Add(from[route[^1]]);
        route.Reverse();
        return [.. route];
    }

    /// <summary>The widest value by fixed-point relaxation over the whole cube, sharing nothing with the search.</summary>
    public static double Relaxed(double[] pi, int dim, int start, int end)
    {
        var best = new double[1 << dim];
        Array.Fill(best, double.NegativeInfinity);
        best[start] = pi[start];
        for (bool moved = true; moved;)
        {
            moved = false;
            for (int u = 0; u < best.Length; u++)
                for (int i = 0; i < dim; i++)
                {
                    int v = u ^ (1 << i);
                    double through = Math.Min(best[u], pi[v]);
                    if (through > best[v])
                    {
                        best[v] = through;
                        moved = true;
                    }
                }
        }
        return best[end];
    }

    /// <summary>The best lowest π over every ordering of the bits that differ, each flipped once: the brute force over
    /// monotone routes. A non-monotone route can only beat it.</summary>
    public static double Orderings(double[] pi, int start, int end)
    {
        int[] differ = [.. Enumerable.Range(0, 31).Where(i => ((start ^ end) >> i & 1) != 0)];
        double best = double.NegativeInfinity;
        void Walk(int at, int left, double low)
        {
            if (low <= best) return;
            if (left == 0)
            {
                best = low;
                return;
            }
            for (int k = 0; k < differ.Length; k++)
                if ((left >> k & 1) != 0)
                {
                    int next = at ^ (1 << differ[k]);
                    Walk(next, left & ~(1 << k), Math.Min(low, pi[next]));
                }
        }
        Walk(start, (1 << differ.Length) - 1, pi[start]);
        return best;
    }

    /// <summary>
    /// The least count of additions with Δπ &lt; 0 from <paramref name="from"/> to <paramref name="target"/> over the states
    /// at or above <paramref name="floor"/>: a 0–1 search, for how far k depends on which maximin route was taken.
    /// </summary>
    public static int LeastAdditions(double[] pi, int dim, int from, int target, double floor)
    {
        var least = new int[1 << dim];
        Array.Fill(least, int.MaxValue);
        var deque = new LinkedList<int>();
        least[from] = 0;
        deque.AddFirst(from);
        while (deque.Count > 0)
        {
            int u = deque.First!.Value;
            deque.RemoveFirst();
            if (u == target) return least[u];
            for (int i = 0; i < dim; i++)
            {
                int v = u ^ (1 << i);
                if (pi[v] < floor) continue;
                int cost = v > u && pi[v] < pi[u] ? 1 : 0;   // the bit set is the addition
                if (least[u] + cost < least[v])
                {
                    least[v] = least[u] + cost;
                    if (cost == 0) deque.AddFirst(v);
                    else deque.AddLast(v);
                }
            }
        }
        return -1;
    }
}

/// <summary>One step of a route, read from S toward T at one τ reading: the term, the changes in the log marginal, the log
/// prior and ln π, c from SC-5's form, the proposal both ways with the move's selection, ln c_e, and for a term added or
/// deleted the class flow around the lower structure's added term, labelled an approximation (brief's ±25 knots).</summary>
internal sealed record PathHop(PathHopKind Kind, string Term, int MFrom, int MTo, double PiFrom, double PiTo, double DMarginal, double DPrior, double C,
                               double LnQForward, double LnQReverse, double LnFlowClass, int ClassSize)
{
    public double DPi => PiTo - PiFrom;

    /// <summary>ln c_e = min(π(x) + ln q(x → y), π(y) + ln q(y → x)).</summary>
    public double LnFlow => Math.Min(PiFrom + LnQForward, PiTo + LnQReverse);

    /// <summary>The hop read the other way.</summary>
    public PathHop Reversed => new(Kind == PathHopKind.Slide ? PathHopKind.Slide : Kind == PathHopKind.Add ? PathHopKind.Delete : PathHopKind.Add,
                                   Term, MTo, MFrom, PiTo, PiFrom, -DMarginal, -DPrior, C, LnQReverse, LnQForward, LnFlowClass, ClassSize);
}

/// <summary>How S and T line up: the terms they share, the pairs a slide takes together, and what is left for the route.</summary>
internal sealed record PathAlignment(BassTerm[] Common, (BassTerm From, BassTerm To)[] Matched, BassTerm[] OnlyFrom, BassTerm[] OnlyTo, int Dropped);

/// <summary>A route between S and T at one reading: its hops in order from S, its lattice, where the slide went, the route's
/// maximin on the lattice, and the superset route's lowest π and k on the same lattice for comparison.</summary>
internal sealed record PathRoute(string Alignment, PathReading Reading, string Position, PathAlignment Aligned, PathHop[] Hops, PathLattice Lattice,
                                 int[] Masks, double LatticeMaximin, double SupersetLow, int SupersetKFrom, int SupersetKTo, double LnResistanceBest)
{
    public double Low => Hops.Length == 0 ? double.NaN : Math.Min(Hops.Min(h => h.PiFrom), Hops[^1].PiTo);
}

/// <summary>A route read from one end: depths, k, the highest M, the waiting time's logs and their split.</summary>
internal sealed record PathDirection(double PiFrom, double PiTo, double Low, double DepthFrom, double DepthTo, int K, int KLeast, int KLattice, double SlideLow, int DeletionsBelow,
                                     int HighestM, double LnWaitExact, double LnWaitClass, double LnWaitRouteClass, double LnWaitSlide, double LnWaitBest, double WorstFall,
                                     double WorstProposal, double WorstClassFall, double WorstClassProposal, double TotalFall, double TotalProposal, double TotalProposalClass);

/// <summary>
/// The path profile between real modes (the path-profile brief; gap G1, SC-4, SC-5): the maximin route through L between two
/// structures, exact or aligned, and its readings: k, depth, ln c_e from the production proposal, the class flow, the split
/// between the fall in π and the proposal cost, and the waiting time. It reads, it does not judge.
/// </summary>
internal static class BassPathProfile
{
    /// <summary>The largest |D| enumerated. Past it the route would need the brief's beam search.</summary>
    internal const int MaxDim = 22;

    /// <summary>The threads a lattice scores on.</summary>
    internal const int Threads = 4;

    /// <summary>The class flow's window, per factor, around a term's knots (the brief's ±25).</summary>
    internal const int ClassWindow = 25;

    /// <summary>The largest |D| whose every route is priced by its flows, for the least waiting time in L.</summary>
    private const int MaxDimPriced = 14;

    /// <summary>The largest knot jump a slide tries when every single step from its state is refused.</summary>
    private const int SlideReach = 8;

    /// <summary>What a sweep attempts: one structure move (<see cref="MoveStep{TState}"/>); τ's slice step proposes none.</summary>
    internal const double MovesPerSweep = 1.0;

    /// <summary>A term's signature: its variables, signs and categorical subsets, knots left out.</summary>
    internal static string Signature(BassTerm t) => string.Join("·",
        t.Continuous.OrderBy(h => h.Variable).Select(h => $"x{h.Variable}{(h.Sign > 0 ? "+" : "-")}")
         .Concat(t.Categorical.Select(c => $"z{c.Variable}{{{string.Join(",", c.Levels)}}}"))
         .Concat(t.Functional.OrderBy(h => h.Variable).Select(h => $"t{h.Variable}{(h.Sign > 0 ? "+" : "-")}")));

    // The knot distance between two terms of one signature.
    private static int Distance(BassTerm a, BassTerm b)
    {
        static IEnumerable<(BassPart Part, HingeFactor H)> Hinges(BassTerm t)
            => t.Continuous.Select(h => (BassPart.Continuous, h)).Concat(t.Functional.Select(h => (BassPart.Functional, h)));
        var at = Hinges(b).ToDictionary(x => (x.Part, x.H.Variable), x => x.H.Knot);
        return Hinges(a).Sum(x => Math.Abs(x.H.Knot - at[(x.Part, x.H.Variable)]));
    }

    /// <summary>
    /// S and T lined up. <paramref name="exact"/>: terms match only if identical. Otherwise terms of one signature pair up by
    /// the least total knot distance, every term of the smaller side paired. With <paramref name="closed"/>, then also a pair
    /// that test holds true for, whatever its signs.
    /// </summary>
    internal static PathAlignment Align(IReadOnlyList<BassTerm> s, IReadOnlyList<BassTerm> t, bool exact, int dropped = 0, HashSet<BassTerm>? excluded = null,
                                        Func<BassTerm, BassTerm, bool>? closed = null)
    {
        BassTerm[] common = [.. s.Where(t.Contains)];
        List<BassTerm> onlyS = [.. s.Where(x => !common.Contains(x))], onlyT = [.. t.Where(x => !common.Contains(x))];
        var matched = new List<(BassTerm, BassTerm)>();
        if (!exact)
            foreach (var group in onlyS.Where(x => excluded is null || !excluded.Contains(x)).GroupBy(Signature).ToArray())
            {
                BassTerm[] left = [.. group];
                BassTerm[] right = [.. onlyT.Where(x => Signature(x) == group.Key && (excluded is null || !excluded.Contains(x)))];
                if (right.Length == 0) continue;
                bool leftSmaller = left.Length <= right.Length;
                BassTerm[] few = leftSmaller ? left : right, many = leftSmaller ? right : left;
                foreach ((int i, int j) in Assign(few, many))
                {
                    BassTerm from = leftSmaller ? few[i] : many[j], to = leftSmaller ? many[j] : few[i];
                    matched.Add((from, to));
                    onlyS.Remove(from);
                    onlyT.Remove(to);
                }
            }
        if (!exact && closed is not null)
            foreach (BassTerm from in onlyS.ToArray())
                if (onlyT.FirstOrDefault(to => closed(from, to)) is { } to)
                {
                    matched.Add((from, to));
                    onlyS.Remove(from);
                    onlyT.Remove(to);
                }
        return new PathAlignment(common, [.. matched], [.. onlyS], [.. onlyT], dropped);
    }

    /// <summary>
    /// Whether the design refuses <paramref name="a"/> and <paramref name="b"/> together though one redraw takes one to the other:
    /// one factor apart, as a hinge at its variable's lowest knot and its mirror at the highest, one linear column twice. D06 then
    /// closes every superset holding both, and L can only pass by deleting a load-bearing term first, while production's redraw
    /// makes the exchange in one step.
    /// </summary>
    internal static bool ClosedPair(PathSurface surface, BassTerm a, BassTerm b)
    {
        if (a.Equals(b) || a.Continuous.Count != b.Continuous.Count || a.Categorical.Count != 0 || b.Categorical.Count != 0 || a.Functional.Count != 0 || b.Functional.Count != 0) return false;
        int[] va = [.. a.Continuous.Select(h => h.Variable).Order()], vb = [.. b.Continuous.Select(h => h.Variable).Order()];
        if (!va.SequenceEqual(vb)) return false;
        int apart = a.Continuous.Count(h => !b.Continuous.Contains(h));
        return apart == 1 && !surface.Score([a, b], profile: false).Admitted;
    }

    // Each of `few` to a distinct one of `many`, by the least total distance: a bitmask search, the groups being small.
    private static (int, int)[] Assign(BassTerm[] few, BassTerm[] many)
    {
        var memo = new Dictionary<(int, int), (int Cost, int Pick)>();
        (int Cost, int Pick) Best(int i, int used)
        {
            if (i == few.Length) return (0, -1);
            if (memo.TryGetValue((i, used), out var held)) return held;
            (int Cost, int Pick) best = (int.MaxValue, -1);
            for (int j = 0; j < many.Length; j++)
            {
                if ((used >> j & 1) != 0) continue;
                int cost = Distance(few[i], many[j]) + Best(i + 1, used | 1 << j).Cost;
                if (cost < best.Cost) best = (cost, j);
            }
            memo[(i, used)] = best;
            return best;
        }
        var pairs = new (int, int)[few.Length];
        int mask = 0;
        for (int i = 0; i < few.Length; i++)
        {
            int j = Best(i, mask).Pick;
            pairs[i] = (i, j);
            mask |= 1 << j;
        }
        return pairs;
    }

    /// <summary>
    /// The hops of a greedy slide from <paramref name="start"/>, each matched term stepped to its target one knot at a time,
    /// always the step with the highest π; a jump of up to <see cref="SlideReach"/> where every single step is refused; and a
    /// factor whose sign differs, a closed pair's, by one redraw. Null, with the pair's term, when a pair cannot move. Each
    /// structure along it, the start's first, goes to <paramref name="visited"/> when given.
    /// </summary>
    internal static PathHop[]? Slide(PathSurface surface, IReadOnlyList<BassTerm> start, (BassTerm From, BassTerm To)[] moves, PathReading r, out BassTerm? stuck,
                                      List<BassTerm[]>? visited = null)
    {
        stuck = null;
        BassTerm[] current = [.. start];
        visited?.Add([.. current]);
        int[] slots = [.. moves.Select(m => Array.IndexOf(current, m.From))];
        bool profile = r == PathReading.Profiled;
        PathScore here = surface.Score(current, profile);
        var hops = new List<PathHop>();
        while (true)
        {
            var pending = new List<(int Pair, int Slot, BassMixing.FactorAt F, int Remaining, HingeFactor? Redraw)>();
            for (int p = 0; p < moves.Length; p++)
            {
                BassTerm now = current[slots[p]], to = moves[p].To;
                foreach (BassMixing.FactorAt f in BassMixing.Factors(now).Where(f => f.Part != BassPart.Categorical))
                {
                    HingeFactor h = f.Part == BassPart.Continuous ? now.Continuous[f.Index] : now.Functional[f.Index];
                    IReadOnlyList<HingeFactor> targets = f.Part == BassPart.Continuous ? to.Continuous : to.Functional;
                    HingeFactor target = targets.Single(x => x.Variable == h.Variable);
                    if (target.Sign != h.Sign) pending.Add((p, slots[p], f, 1, target));
                    else if (target.Knot != h.Knot) pending.Add((p, slots[p], f, target.Knot - h.Knot, null));
                }
            }
            if (pending.Count == 0) return [.. hops];

            (int Slot, BassTerm Term, PathScore Score, BassMixing.FactorAt F)? best = null;
            for (int reach = 1; reach <= SlideReach && best is null; reach++)
                foreach ((int _, int slot, BassMixing.FactorAt f, int remaining, HingeFactor? redraw) in pending)
                {
                    if (reach > Math.Abs(remaining)) continue;
                    BassTerm stepped = redraw is { } g ? BassMixing.Replaced(current[slot], f, g, null) : BassMixing.Stepped(surface.Model, current[slot], f, Math.Sign(remaining) * reach)!;
                    if (current.Contains(stepped)) continue;
                    BassTerm[] next = [.. current];
                    next[slot] = stepped;
                    PathScore score = surface.Score(next, profile);
                    if (score.Admitted && (best is null || score.Pi(r) > best.Value.Score.Pi(r)))
                        best = (slot, stepped, score, f);
                }
            if (best is not { } taken)
            {
                stuck = moves[pending[0].Pair].From;
                return null;
            }
            double lnQ = surface.LogRedraw(current.Length, current[taken.Slot], taken.F.Part, BassMixing.Variable(current[taken.Slot], taken.F));
            double dm = taken.Score.Marginal(r) - here.Marginal(r), dp = taken.Score.Prior(r) - here.Prior(r);
            hops.Add(new PathHop(PathHopKind.Slide, BassMixing.Key([taken.Term]), current.Length, current.Length, here.Pi(r), taken.Score.Pi(r), dm, dp, 0.0, lnQ, lnQ,
                                 double.NaN, 0));
            current[taken.Slot] = taken.Term;
            visited?.Add([.. current]);
            here = taken.Score;
        }
    }

    /// <summary>The class flow for adding <paramref name="t"/> at <paramref name="x"/>: Σ over the terms of its signature within
    /// ±<see cref="ClassWindow"/> knots per factor of min(π(x) + ln q_b, π(x + t′) + ln q_d). An approximation: a neighbour's
    /// structure is not the route's next state.</summary>
    internal static (double LnFlow, int Size) ClassFlow(PathSurface surface, IReadOnlyList<BassTerm> x, BassTerm t, PathScore scoreX, PathReading r)
    {
        BassDesign design = surface.Design;
        Assert.Null(design.Load(x));
        bool profile = r == PathReading.Profiled;
        double piX = scoreX.Pi(r), lnForward = surface.LogBirth(x, x.Count, t), lnReverse = surface.LogDeath(x.Count + 1);
        double logTermPrior = surface.LogTermPrior(x) + surface.Held.Prior.LogPrior(t);
        double sum = double.NegativeInfinity;
        int size = 0;
        foreach (BassTerm near in Neighbours(surface, t))
        {
            if (x.Contains(near)) continue;
            DesignProposal proposal = design.ProposeBirth(near);
            PathScore score = surface.Score(BassStructureFit.Of(proposal, proposal.Admitted ? logTermPrior : double.NaN), profile);
            size++;
            if (!score.Admitted) continue;
            sum = LogMath.Add(sum, Math.Min(piX + lnForward, score.Pi(r) + lnReverse));
        }
        return (sum, size);
    }

    // The terms of t's signature within the window, enumerated factor by factor.
    internal static IEnumerable<BassTerm> Neighbours(PathSurface surface, BassTerm t)
    {
        IEnumerable<HingeFactor[]> Spread(IReadOnlyList<HingeFactor> hinges, BassPart part)
        {
            IEnumerable<HingeFactor[]> built = [[]];
            foreach (HingeFactor h in hinges)
            {
                int n = BassMixing.Knots(surface.Model, part, h.Variable);
                int[] knots = [.. Enumerable.Range(Math.Max(0, h.Knot - ClassWindow), Math.Min(n - 1, h.Knot + ClassWindow) - Math.Max(0, h.Knot - ClassWindow) + 1)];
                built = built.SelectMany(prefix => knots.Select(k => (HingeFactor[])[.. prefix, h with { Knot = k }]));
            }
            return built;
        }
        foreach (HingeFactor[] c in Spread(t.Continuous, BassPart.Continuous))
            foreach (HingeFactor[] f in Spread(t.Functional, BassPart.Functional))
                yield return new BassTerm(c, t.Categorical, f);
    }

    /// <summary>The route's hops through the lattice, S's mask to T's, each state loaded afresh and checked against the
    /// lattice's score; SC-5's inequality checked at the lower state's best τ on the profile.</summary>
    private static List<PathHop> LatticeHops(PathSurface surface, PathLattice lattice, int[] masks, PathReading r, bool priced)
    {
        var hops = new List<PathHop>(masks.Length);
        double[] pi = lattice.Pi(r);
        for (int k = 0; k + 1 < masks.Length; k++)
        {
            int u = masks[k], v = masks[k + 1];
            BassTerm[] x = lattice.Terms(u), y = lattice.Terms(v);
            PathScore sx = surface.Score(x), sy = surface.Score(y);
            Assert.Equal(pi[u], sx.Pi(r), 6);
            Assert.Equal(pi[v], sy.Pi(r), 6);
            int bit = System.Numerics.BitOperations.TrailingZeroCount(u ^ v);
            BassTerm t = lattice.Bits[bit];
            bool add = (v & 1 << bit) != 0;
            (BassTerm[] lower, PathScore sLower, PathScore sUpper) = add ? (x, sx, sy) : (y, sy, sx);
            double logBirth = surface.LogBirth(lower, lower.Length, t), logDeath = surface.LogDeath(lower.Length + 1);
            double u0 = r == PathReading.Held ? Math.Log(surface.HeldTau) : sLower.BestLogTau, c = PathSurface.Occam(u0);
            if (r == PathReading.Profiled)
            {
                double slack = surface.MarginalAt(sUpper.Fit, u0) - (surface.MarginalAt(sLower.Fit, u0) - c);
                if (slack < -1e-7) throw new InvalidOperationException($"SC-5's inequality fails by {-slack:E3} at the lower state's best τ, adding {BassMixing.Key([t])}.");
                surface.CountSuperset(1, slack);
            }
            (double lnClass, int size) = priced ? ClassFlow(surface, lower, t, sLower, r) : (double.NaN, 0);
            hops.Add(new PathHop(add ? PathHopKind.Add : PathHopKind.Delete, BassMixing.Key([t]), x.Length, y.Length, sx.Pi(r), sy.Pi(r),
                                 sy.Marginal(r) - sx.Marginal(r), sy.Prior(r) - sx.Prior(r), c,
                                 add ? logBirth : logDeath, add ? logDeath : logBirth, lnClass, size));
        }
        return hops;
    }

    /// <summary>
    /// The route between <paramref name="s"/> and <paramref name="t"/> at one reading. Exact: the maximin route through L on the
    /// exact symmetric difference. Aligned: the slide first or last, whichever gives the higher minimum, the route through L on
    /// what remains. A pair whose slide cannot move drops out of the matching. <paramref name="lattices"/> holds the pair's
    /// lattices across its readings.
    /// </summary>
    internal static PathRoute? Build(PathSurface surface, BassTerm[] s, BassTerm[] t, bool exact, PathReading r, Dictionary<string, PathLattice> lattices,
                                     IEnumerable<BassTerm>? unmatched = null, bool priced = true, bool redraw = false)
    {
        var dropped = new HashSet<BassTerm>(unmatched ?? []);
        Func<BassTerm, BassTerm, bool>? closed = redraw ? (a, b) => ClosedPair(surface, a, b) : null;
        while (true)
        {
            PathAlignment a = Align(s, t, exact, dropped.Count / 2, dropped, closed);
            BassTerm[] bits = [.. a.OnlyFrom, .. a.OnlyTo];
            if (a.Matched.Length == 0)
            {
                PathLattice plain = Lattice(surface, a.Common, bits, lattices);
                return Finish(surface, plain, a, r, exact, "none", [], [], s, t, priced);
            }
            BassTerm[] toKnots = [.. a.Matched.Select(m => m.To)], fromKnots = [.. a.Matched.Select(m => m.From)];
            PathHop[]? first = Slide(surface, s, a.Matched, r, out BassTerm? stuck);
            PathHop[]? last = stuck is null ? Slide(surface, t, [.. a.Matched.Select(m => (m.To, m.From))], r, out stuck) : null;
            if (first is null || last is null)
            {
                BassTerm from = stuck!;
                (BassTerm From, BassTerm To) pair = a.Matched.First(m => m.From.Equals(from) || m.To.Equals(from));
                dropped.Add(pair.From);
                dropped.Add(pair.To);
                continue;
            }
            PathLattice afterSlide = Lattice(surface, [.. a.Common, .. toKnots], bits, lattices);
            PathLattice beforeSlide = Lattice(surface, [.. a.Common, .. fromKnots], bits, lattices);
            PathRoute? slideFirst = Finish(surface, afterSlide, a, r, exact, "first", first, [], [.. a.Common, .. toKnots, .. a.OnlyFrom], t, priced);
            PathRoute? slideLast = Finish(surface, beforeSlide, a, r, exact, "last", [], [.. Enumerable.Reverse(last).Select(h => h.Reversed)], s, [.. a.Common, .. fromKnots, .. a.OnlyTo], priced);
            if (slideFirst is null) return slideLast;
            if (slideLast is null) return slideFirst;
            return slideFirst.Low >= slideLast.Low ? slideFirst : slideLast;
        }
    }

    /// <summary>
    /// The aligned route with far pairs let go, an addition to the brief's three-way reading: greedily, farthest pair first, a
    /// matched pair goes to the route as an exchange whenever that raises the route's lowest π. The slide takes the pairs together
    /// whatever their distance, so it can cost more than exchanging the far ones; the search for a higher route stops at the first
    /// structure of dropped pairs no single further drop improves.
    /// </summary>
    internal static PathRoute? BuildPartial(PathSurface surface, BassTerm[] s, BassTerm[] t, PathReading r, Dictionary<string, PathLattice> lattices, bool redraw = false)
    {
        var dropped = new HashSet<BassTerm>();
        PathRoute? best = Build(surface, s, t, false, r, lattices, priced: false, redraw: redraw);
        if (best is null) return null;
        for (bool improved = true; improved;)
        {
            improved = false;
            foreach ((BassTerm from, BassTerm to) in best.Aligned.Matched.OrderByDescending(m => Distance(m.From, m.To)))
            {
                PathRoute? trial = Build(surface, s, t, false, r, lattices, [.. dropped, from, to], priced: false, redraw: redraw);
                if (trial is null || trial.Low <= best.Low + 1e-6) continue;
                best = trial;
                dropped.Add(from);
                dropped.Add(to);
                improved = true;
                break;
            }
        }
        return Build(surface, s, t, false, r, lattices, dropped, redraw: redraw)! with { Alignment = "partial" };
    }

    private static PathLattice Lattice(PathSurface surface, BassTerm[] core, BassTerm[] bits, Dictionary<string, PathLattice> lattices)
    {
        string key = $"{BassMixing.Key(core)} | {BassMixing.Key(bits)}";
        if (!lattices.TryGetValue(key, out PathLattice? lattice)) lattices[key] = lattice = new PathLattice(surface, core, bits);
        return lattice;
    }

    private static PathRoute? Finish(PathSurface surface, PathLattice lattice, PathAlignment a, PathReading r, bool exact, string position, PathHop[] before, PathHop[] after,
                                     IReadOnlyList<BassTerm> start, IReadOnlyList<BassTerm> end, bool priced)
    {
        int from = lattice.Mask(start), to = lattice.Mask(end);
        double[] pi = lattice.Pi(r);
        int[]? masks = PathSearch.Route(pi, lattice.Dim, from, to, out double maximin);
        if (masks is null) return null;
        List<PathHop> hops = [.. before, .. LatticeHops(surface, lattice, masks, r, priced), .. after];
        for (int k = 0; k + 1 < hops.Count; k++) Assert.Equal(hops[k].PiTo, hops[k + 1].PiFrom, 6);

        // The superset route on the same lattice: every addition, then every deletion, each half at its own best order.
        int union = from | to;
        int[]? up = PathSearch.Route(pi, lattice.Dim, from, union, out double lowUp, (u, v) => v > u && (v & ~union) == 0);
        int[]? down = PathSearch.Route(pi, lattice.Dim, union, to, out double lowDown, (u, v) => v < u && (u & ~union) == 0);
        double supersetLow = double.NegativeInfinity;
        int kFrom = -1, kTo = -1;
        if (up is not null && down is not null)
        {
            supersetLow = Math.Min(lowUp, lowDown);
            int[] whole = [.. up, .. down.Skip(1)];
            int lowest = Enumerable.Range(0, whole.Length).MinBy(i => pi[whole[i]]);
            kFrom = Enumerable.Range(0, lowest).Count(i => whole[i + 1] > whole[i] && pi[whole[i + 1]] < pi[whole[i]]);
            kTo = Enumerable.Range(lowest, whole.Length - 1 - lowest).Count(i => whole[i] > whole[i + 1] && pi[whole[i]] < pi[whole[i + 1]]);
        }
        double best = priced && lattice.Dim <= MaxDimPriced ? LeastResistance(surface, lattice, r, from, to) : double.NaN;
        foreach (PathHop h in hops)
            if (h.Kind == PathHopKind.Slide) best = LogMath.Add(best, -h.LnFlow);
        return new PathRoute(exact ? "exact" : "aligned", r, position, a, [.. hops], lattice, masks, maximin, supersetLow, kFrom, kTo, best);
    }

    // ln of the least Σ_e 1/c_e over every route through L, from the production proposal's flows; the slide's edges are added
    // by the caller. Whether a faster route than the maximin one lies in L.
    private static double LeastResistance(PathSurface surface, PathLattice lattice, PathReading r, int from, int to)
    {
        double[] pi = lattice.Pi(r);
        var distance = new double[1 << lattice.Dim];
        Array.Fill(distance, double.PositiveInfinity);
        var queue = new PriorityQueue<int, double>();
        distance[from] = double.NegativeInfinity;
        queue.Enqueue(from, double.NegativeInfinity);
        while (queue.TryDequeue(out int u, out double d))
        {
            if (d > distance[u]) continue;
            if (u == to) return d;
            for (int i = 0; i < lattice.Dim; i++)
            {
                int v = u ^ (1 << i);
                if (double.IsNegativeInfinity(pi[v])) continue;
                int lower = Math.Min(u, v), upper = Math.Max(u, v), m = lattice.Core.Length + System.Numerics.BitOperations.PopCount((uint)lower);
                double lnFlow = Math.Min(pi[lower] + surface.LogBirth(lattice.Terms(lower), m, lattice.Bits[i]), pi[upper] + surface.LogDeath(m + 1));
                double through = LogMath.Add(d, -lnFlow);
                if (through < distance[v])
                {
                    distance[v] = through;
                    queue.Enqueue(v, through);
                }
            }
        }
        return double.NaN;
    }

    /// <summary>The route read from S (<paramref name="forward"/>) or from T: depths and k to the route's lowest state, the
    /// waiting time Σ_e exp(π(from) − ln c_e) in attempts, exact and by class (the class flow on additions in this direction,
    /// the exact flow elsewhere), and its split at the costliest edge and over all edges.</summary>
    internal static PathDirection Read(PathRoute route, bool forward)
    {
        PathHop[] hops = forward ? route.Hops : [.. Enumerable.Reverse(route.Hops).Select(h => h.Reversed)];
        double[] along = [hops[0].PiFrom, .. hops.Select(h => h.PiTo)];
        double piFrom = along[0], piTo = along[^1], low = along.Min();
        int lowest = Array.IndexOf(along, low);
        int k = 0, below = 0, highest = hops[0].MFrom;
        for (int i = 0; i < hops.Length; i++)
        {
            highest = Math.Max(highest, hops[i].MTo);
            if (i >= lowest || hops[i].DPi >= 0.0) continue;
            if (hops[i].Kind == PathHopKind.Add) k++;
            else if (hops[i].Kind == PathHopKind.Delete) below++;
        }
        double waitExact = double.NegativeInfinity, waitClass = double.NegativeInfinity, waitRoute = double.NegativeInfinity, waitSlide = double.NegativeInfinity;
        (double Cost, double Fall) worst = (double.NegativeInfinity, 0.0), worstClass = (double.NegativeInfinity, 0.0);
        double totalFall = 0.0, totalProposal = 0.0, totalClass = 0.0;
        foreach (PathHop h in hops)
        {
            double lnFlow = h.LnFlow, side = h.PiFrom + h.LnQForward <= h.PiTo + h.LnQReverse ? h.PiFrom : h.PiTo;
            double fall = piFrom - side, cost = piFrom - lnFlow;
            double lnClass = h.Kind == PathHopKind.Add ? Math.Max(lnFlow, h.LnFlowClass) : lnFlow, costClass = piFrom - lnClass;
            waitExact = LogMath.Add(waitExact, cost);
            waitClass = LogMath.Add(waitClass, costClass);
            if (h.Kind == PathHopKind.Slide) waitSlide = LogMath.Add(waitSlide, cost);
            else waitRoute = LogMath.Add(waitRoute, costClass);
            if (cost > worst.Cost) worst = (cost, fall);
            if (costClass > worstClass.Cost) worstClass = (costClass, fall);
            totalFall += fall;
            totalProposal += cost - fall;
            totalClass += costClass - fall;
        }
        int leastK = LeastK(route, forward, k);
        PathHop[] onLattice = [.. hops.Where(h => h.Kind != PathHopKind.Slide)];
        int kLattice = 0;
        if (onLattice.Length > 0)
        {
            double[] states = [onLattice[0].PiFrom, .. onLattice.Select(h => h.PiTo)];
            int lowestState = Array.IndexOf(states, states.Min());
            kLattice = onLattice.Take(lowestState).Count(h => h.Kind == PathHopKind.Add && h.DPi < 0.0);
        }
        PathHop[] slides = [.. hops.Where(h => h.Kind == PathHopKind.Slide)];
        double slideLow = slides.Length == 0 ? double.NaN : slides.Min(h => Math.Min(h.PiFrom, h.PiTo));
        return new PathDirection(piFrom, piTo, low, piFrom - low, piTo - low, k, leastK, kLattice, slideLow, below, highest, waitExact, waitClass, waitRoute, waitSlide,
                                 piFrom + route.LnResistanceBest, worst.Fall, worst.Cost - worst.Fall, worstClass.Fall, worstClass.Cost - worstClass.Fall, totalFall, totalProposal, totalClass);
    }

    // The least k over the lattice's maximin routes, from the lattice's entry in this direction to the route's lowest
    // lattice state; the route's own k when the route's lowest state lies on the slide.
    private static int LeastK(PathRoute route, bool forward, int k)
    {
        double[] pi = route.Lattice.Pi(route.Reading);
        int[] masks = forward ? route.Masks : [.. Enumerable.Reverse(route.Masks)];
        double low = masks.Min(m => pi[m]);
        double routeLow = route.Low;
        if (routeLow < low - 1e-9) return k;
        int target = masks.First(m => pi[m] == low);
        return PathSearch.LeastAdditions(pi, route.Lattice.Dim, masks[0], target, route.LatticeMaximin);
    }

    /// <summary>The rows' number format.</summary>
    internal static string F(double v) => double.IsNaN(v) ? "" : double.IsNegativeInfinity(v) ? "-inf" : double.IsPositiveInfinity(v) ? "inf" : v.ToString("F6", CultureInfo.InvariantCulture);

    internal static string I(long v) => v.ToString(CultureInfo.InvariantCulture);
}

/// <summary>A CSV beside the brief, appended row by row and skipping rows already present, so a dead run reruns from where it
/// stopped. A row's key is its first <c>keyColumns</c> cells.</summary>
internal sealed class PathCsv
{
    private static readonly object Gate = new();
    private readonly HashSet<string> _keys = [];
    private readonly int _keyColumns;
    private readonly int _columns;

    internal PathCsv(string file, string[] header, int keyColumns, string prefix = "path-profile")
    {
        Path = System.IO.Path.Combine(Folder, $"{prefix}-{file}.csv");
        _keyColumns = keyColumns;
        _columns = header.Length;
        lock (Gate)
        {
            if (!File.Exists(Path)) File.WriteAllText(Path, string.Join(",", header) + "\n");
            foreach (string[] row in Rows()) _keys.Add(string.Join("\u0001", row.Take(keyColumns)));
        }
    }

    public string Path { get; }

    /// <summary>The brief's folder: <c>BASS_PATH_OUT</c>, or found upward from the test's folder.</summary>
    internal static string Folder
    {
        get
        {
            if (Environment.GetEnvironmentVariable("BASS_PATH_OUT") is { Length: > 0 } declared) return declared;
            for (DirectoryInfo? at = new(AppContext.BaseDirectory); at is not null; at = at.Parent)
            {
                string candidate = System.IO.Path.Combine(at.FullName, "issues", "regression", "bass");
                if (Directory.Exists(candidate)) return candidate;
            }
            throw new DirectoryNotFoundException("No issues/regression/bass above the test's folder; set BASS_PATH_OUT.");
        }
    }

    public bool Has(params string[] key) => _keys.Contains(string.Join("\u0001", key));

    public void Append(params string[] cells)
    {
        Assert.Equal(_columns, cells.Length);
        string key = string.Join("\u0001", cells.Take(_keyColumns));
        lock (Gate)
        {
            if (!_keys.Add(key)) return;
            File.AppendAllText(Path, string.Join(",", cells) + "\n");
        }
    }

    /// <summary>The rows, the header left out. Cells hold no commas.</summary>
    public IEnumerable<string[]> Rows() => File.ReadLines(Path).Skip(1).Where(l => l.Length > 0).Select(l => l.Split(','));
}
