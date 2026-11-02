using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;
using F = Maths.Regression.Tests.BassPathProfile;

namespace Maths.Regression.Tests;

/// <summary>
/// The path profile between Friedman's real modes (the path-profile brief; gap G1).
/// <list type="bullet">
/// <item>Always run: the search against brute force on random cubes, a route's own arithmetic, and the flow against the
/// enumerated kernel.</item>
/// <item>Qualification runs (<see cref="BassQualification.Enabled"/>): the modes from one flat-histogram run and the fixed
/// states; the ground truth on enumerated fixtures; Part 1's routes; Part 2's routes against N; the stored exemplars rechecked.
/// Rows go beside the brief as <c>path-profile-*.csv</c>, appended as computed and skipped on a rerun.</item>
/// </list>
/// Choices the brief leaves open:
/// <list type="bullet">
/// <item>Exemplars rank by π at the held τ. Beside the brief's two highest, each basin also gets its highest structure whose
/// terms' signs differ from the top one's, since the two highest often differ by one knot.</item>
/// <item>The maximin value does not decide the route past its lowest state: <see cref="PathSearch.Route"/> takes the least
/// Σ_e exp(−min π), reading π alone, and <c>k_least</c> is the least k over all maximin routes.</item>
/// <item>The class flow applies to additions in the direction read; a deletion has one term to delete.</item>
/// <item>The split is read at the costliest edge, whose cost is the log waiting time to within ln(edges), with Σ over all edges
/// beside it.</item>
/// <item>Beside exact and aligned, a third reading, <c>partial</c> (<see cref="BassPathProfile.BuildPartial"/>): a matched pair goes
/// to the route as an exchange when that raises the route's lowest π. Aligned stays the primary reading.</item>
/// <item>A closed pair (<see cref="BassPathProfile.ClosedPair"/>), one column twice under two signs, shuts L's supersets, which the
/// brief's rule on depth reads as a closed cap. <c>aligned+redraw</c> and <c>partial+redraw</c> take such a pair by one redraw, the
/// move production makes it by; the plain rows stay as they were read. No Part 1 pair holds one.</item>
/// <item>A slide's edges are D01's redraw at one knot's distance, so each carries 1/(2n_v): the local steps are off in S10's
/// configuration. Their share of the waiting time is its own column.</item>
/// <item>Part 2 reads fresh cases from another seed of the fixture's generator, and the generator's own stream, which
/// extends or truncates S10's. The polish holds τ at unit information. The partial reading runs on the fresh cases at the held τ.</item>
/// </list>
/// </summary>
public sealed class BassPathProfileTests
{
    // Master seeds no pilot or earlier harness uses: Part 2's fresh data and its short chains.
    private const int FreshDataSeed = 9101, ShortChainSeed = 9131;
    private const int MaxModesPerFixture = 30;

    private static readonly string[] ModeHeader =
    [
        "label", "exemplar", "term", "factor", "variable", "sign", "knot", "source", "basin", "m", "pi_held", "pi_profiled", "best_log_tau", "marginal_held",
        "structure_prior", "draw_log_tau", "draw_level", "fresh_level", "basin_share", "basin_draws", "structure_draws", "note",
    ];

    private static readonly string[] RouteHeader =
    [
        "n", "data", "pair", "exemplar", "alignment", "reading", "direction", "position", "d", "matched", "dropped", "slide_steps", "states",
        "refused_cap", "refused_support", "refused_numerical", "refused_vif", "refused_exact", "pi_from", "pi_to", "low", "depth_from", "depth_to", "k", "k_least", "k_lattice", "slide_low",
        "deletions_below", "highest_m", "ln_wait_exact", "ln_wait_class", "ln_wait_route_class", "ln_wait_slide", "ln_wait_best_in_l", "worst_fall", "worst_proposal", "worst_class_fall",
        "worst_class_proposal", "total_fall", "total_proposal", "total_proposal_class", "lattice_maximin", "superset_low", "superset_depth_from", "superset_k", "c",
        "hops", "lattice_hops", "seconds",
    ];

    private static readonly string[] StepHeader =
    [
        "n", "data", "pair", "exemplar", "alignment", "reading", "step", "kind", "term", "m_from", "m_to", "d_marginal", "d_prior", "d_pi", "c", "pi_from", "pi_to",
        "ln_q_forward", "ln_q_reverse", "ln_flow", "ln_flow_class", "class_size", "phase_hops", "phase_low", "phase_ln_resistance",
    ];

    private static readonly string[] ValidationHeader = ["check", "case", "cases", "worst", "value_a", "value_b", "value_c", "note"];

    private static readonly string[] GroundHeader =
    [
        "fixture", "cap", "seed", "mode_a", "mode_b", "m_a", "m_b", "pi_a", "pi_b", "maximin_whole", "maximin_birth_death", "route_exact", "route_aligned",
        "route_exact_lifted", "route_aligned_lifted", "d_exact", "d_aligned", "refused_cap_exact",
    ];

    private static readonly string[] CarriedHeader =
    [
        "n", "data", "label", "exemplar", "m", "nudges", "pi_carried", "pi_polished", "pi_profiled", "best_log_tau", "passes", "at_limit", "displacement", "basin_kept",
        "chain_sweeps", "chain_share_in_basin", "chain_end_basin", "chain_top_basin", "chain_top_share", "chain_rank", "chain_median_gap", "structure",
    ];

    private readonly ITestOutputHelper _out;
    public BassPathProfileTests(ITestOutputHelper output) => _out = output;

    private static PathCsv Modes() => new("modes", ModeHeader, 4);
    private static PathCsv Routes() => new("routes", RouteHeader, 7);
    private static PathCsv Steps() => new("steps", StepHeader, 7);
    private static PathCsv Validation() => new("validation", ValidationHeader, 2);

    // ---- Always run ----

    [Fact]
    public void Search_AgreesWithBruteForce_OnRandomCubes()
    {
        var rng = new Xoshiro256PlusPlus(seed: 9151);
        int nonMonotone = 0;
        for (int trial = 0; trial < 400; trial++)
        {
            int dim = 1 + trial % 7;
            var pi = new double[1 << dim];
            for (int i = 0; i < pi.Length; i++) pi[i] = -10.0 * rng.NextDouble();
            if (trial % 5 == 0) pi[BassDraws.Below(rng, pi.Length)] = double.NegativeInfinity;   // a refused state
            int start = BassDraws.Below(rng, pi.Length), end = BassDraws.Below(rng, pi.Length);
            if (double.IsNegativeInfinity(pi[start]) || double.IsNegativeInfinity(pi[end]) || start == end) continue;
            double widest = PathSearch.Widest(pi, dim, start, end), relaxed = PathSearch.Relaxed(pi, dim, start, end), orderings = PathSearch.Orderings(pi, start, end);
            Assert.Equal(relaxed, widest);
            Assert.True(orderings <= widest, "A monotone ordering is a route.");
            if (orderings < widest) nonMonotone++;
            int[]? route = PathSearch.Route(pi, dim, start, end, out double maximin);
            Assert.Equal(widest, maximin);
            if (route is null) continue;
            Assert.Equal(start, route[0]);
            Assert.Equal(end, route[^1]);
            Assert.Equal(maximin, route.Min(m => pi[m]));
            for (int k = 0; k + 1 < route.Length; k++) Assert.Equal(1, System.Numerics.BitOperations.PopCount((uint)(route[k] ^ route[k + 1])));
        }
        _out.WriteLine($"Random cubes: {nonMonotone} of 400 are beaten by a non-monotone route, as random π allows.");
    }

    [Fact]
    public void Flow_MatchesTheEnumeratedKernel()
    {
        (double worstBirth, int births, double worstRedraw, int redraws) = FlowAgainstKernel();
        _out.WriteLine($"ln c_e against π(x)K(x, y): {births} birth edges, worst {worstBirth:E2}; {redraws} redraw edges, worst {worstRedraw:E2}.");
        Assert.True(worstBirth <= 1e-9 && worstRedraw <= 1e-9);
    }

    // Each birth edge's and redraw edge's ln c_e, against the enumerated kernel's π(x)K(x, y): the selection, the draw and the
    // acceptance are the kernel's own, so they agree with no factor left over.
    private static (double WorstBirth, int Births, double WorstRedraw, int Redraws) FlowAgainstKernel()
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s);
        var surface = new PathSurface(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        double worstBirth = 0.0, worstRedraw = 0.0;
        int births = 0, redraws = 0;
        foreach (((int i, int j), double logK) in flows[BassKernelEnumeration.BirthDeath])
        {
            BassKernelEnumeration.State x = space.States[i], y = space.States[j];
            if (y.Terms.Length != x.Terms.Length + 1) continue;
            BassTerm t = y.Terms.Single(u => !x.Terms.Contains(u));
            Assert.Equal(x.LogPi, surface.Score(x.Terms).Pi(PathReading.Held), 9);
            double lnFlow = Math.Min(x.LogPi + surface.LogBirth(x.Terms, x.Terms.Length, t), y.LogPi + surface.LogDeath(y.Terms.Length));
            worstBirth = Math.Max(worstBirth, Math.Abs(lnFlow - (x.LogPi + logK)));
            births++;
        }
        foreach (((int i, int j), double logK) in flows[BassKernelEnumeration.Redraw])
        {
            BassKernelEnumeration.State x = space.States[i], y = space.States[j];
            BassTerm old = x.Terms.Single(u => !y.Terms.Contains(u));
            BassTerm neu = y.Terms.Single(u => !x.Terms.Contains(u));
            int v = old.Continuous.Zip(neu.Continuous).Where(p => p.First != p.Second).Select(p => p.First.Variable).Single();
            double lnQ = surface.LogRedraw(x.Terms.Length, old, BassPart.Continuous, v);
            double lnFlow = Math.Min(x.LogPi, y.LogPi) + lnQ;
            worstRedraw = Math.Max(worstRedraw, Math.Abs(lnFlow - (x.LogPi + logK)));
            redraws++;
        }
        return (worstBirth, births, worstRedraw, redraws);
    }

    [Fact]
    public void Route_ReadsItsOwnArithmetic()
    {
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s);
        var surface = new PathSurface(config with { Structure = config.Structure with { MaxBasis = null } }, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(BassModel.Create(config, data));
        BassKernelEnumeration.State[] top = [.. space.States.Where(s => s.Terms.Length == 3).OrderByDescending(s => s.LogPi).Take(40)];
        BassTerm[] a = top[0].Terms, b = top.First(s => !s.Terms.Intersect(a).Any()).Terms;
        foreach (bool exact in new[] { true, false })
        {
            PathRoute route = BassPathProfile.Build(surface, a, b, exact, PathReading.Held, [])!;
            PathDirection forward = BassPathProfile.Read(route, true), reverse = BassPathProfile.Read(route, false);
            Assert.Equal(surface.Score(a).Pi(PathReading.Held), forward.PiFrom, 9);
            Assert.Equal(surface.Score(b).Pi(PathReading.Held), forward.PiTo, 9);
            Assert.Equal(forward.PiTo - forward.PiFrom, route.Hops.Sum(h => h.DPi), 9);
            foreach (PathHop h in route.Hops) Assert.Equal(h.DPi, h.DMarginal + h.DPrior, 9);
            Assert.Equal(forward.Low, reverse.Low, 12);
            Assert.Equal(forward.DepthFrom, reverse.DepthTo, 12);
            // The resistances are one sum, read from either end.
            Assert.Equal(forward.LnWaitExact - forward.PiFrom, reverse.LnWaitExact - reverse.PiFrom, 9);
            if (exact) Assert.True(route.Hops.All(h => h.Kind != PathHopKind.Slide));
            Assert.True(surface.SupersetEdges > 0);
        }
    }

    // ---- The modes ----

    private static readonly Dictionary<string, string> Basins = new()
    {
        ["a1"] = "4 × x1·x2 + 2 × x3 + x4 + x5",
        ["a1x"] = "4 × x1·x2 + 2 × x3 + x4 + x5 + x6·x8",
        ["a2"] = "4 × x1·x2 + x2 + 2 × x3 + x4 + x5",
        ["A"] = "x1 + 2 × x1·x2 + x2 + 2 × x3 + x4 + x5",
        ["B"] = "x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5",
    };

    /// <summary>One stored exemplar.</summary>
    internal sealed record Exemplar(string Label, string Rank, string Source, string Basin, BassTerm[] Terms, double PiHeld, double PiProfiled, double MarginalHeld,
                                    double DrawLogTau, double DrawLevel, double FreshLevel);

    /// <summary>The exemplars stored beside the brief, by label and rank; a basin the run missed has none.</summary>
    internal static Dictionary<(string Label, string Rank), Exemplar> ReadModes()
    {
        var modes = new Dictionary<(string, string), Exemplar>();
        foreach (var structure in Modes().Rows().Where(r => r[1] != "missing").GroupBy(r => (r[0], r[1])))
        {
            string[] first = structure.First();
            BassTerm[] terms = [.. structure.GroupBy(r => r[2]).Select(term => new BassTerm(
                [.. term.Select(r => new HingeFactor(int.Parse(r[4], CultureInfo.InvariantCulture), int.Parse(r[5], CultureInfo.InvariantCulture), int.Parse(r[6], CultureInfo.InvariantCulture)))], null, null))];
            static double D(string v) => double.Parse(v, CultureInfo.InvariantCulture);
            modes[structure.Key] = new Exemplar(first[0], first[1], first[7], first[8], terms, D(first[10]), D(first[11]), D(first[13]), D(first[15]), D(first[16]), D(first[17]));
        }
        return modes;
    }

    // One structure's rows: term by term, factor by factor, with its exact variable, sign and knot.
    private static void Store(PathCsv csv, string label, string rank, string source, string basin, PathSurface surface, IReadOnlyList<BassTerm> terms,
                              double drawLogTau, double drawLevel, double share, long basinDraws, long structureDraws, string note)
    {
        PathScore score = surface.Score(terms);
        Assert.True(score.Admitted);
        double fresh = surface.Sampled!.Target.Evaluate(score.Fit, drawLogTau).LogMarginal;
        for (int t = 0; t < terms.Count; t++)
        {
            BassTerm term = terms[t];
            Assert.Empty(term.Categorical);
            Assert.Empty(term.Functional);
            for (int f = 0; f < term.Continuous.Count; f++)
            {
                HingeFactor h = term.Continuous[f];
                csv.Append(label, rank, F.I(t), F.I(f), F.I(h.Variable), F.I(h.Sign), F.I(h.Knot), source, basin, F.I(terms.Count), F.F(score.Pi(PathReading.Held)),
                           F.F(score.Pi(PathReading.Profiled)), F.F(score.BestLogTau), F.F(score.HeldMarginal), F.F(score.LogStructurePrior), F.F(drawLogTau),
                           F.F(drawLevel), F.F(fresh), F.F(share), F.I(basinDraws), F.I(structureDraws), note);
            }
        }
    }

    [Fact]
    public void Modes_FromTheChainsAndTheFlatHistogram()
    {
        if (!BassQualification.Enabled) return;   // heavy: the flat-histogram arm is about half an hour
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        var surface = new PathSurface(fixture.Config, fixture.Data);
        PathCsv csv = Modes();
        if (!csv.Has("chain0", "1", "0", "0")) StoreChains(csv, fixture, surface);
        if (!csv.Rows().Any(r => r[0] == "a1")) StoreFlatHistogram(csv, fixture, surface);
    }

    private void StoreChains(PathCsv csv, BassMixing.Fixture fixture, PathSurface surface)
    {
        var clock = Stopwatch.StartNew();
        BassMixing.FixedState[] states = BassMixing.FixedStates(surface.Model);
        foreach (string line in BassMixing.FixedTable(surface.Model, states)) _out.WriteLine(line);
        double top = states.Max(s => s.Median);
        foreach (BassMixing.FixedState s in states)
        {
            BassState state = s.Chain.Current;
            string note = $"rank {s.Rank:F3} median {s.Median:F2} gap {s.Gap:F2} {(s.Substituted ? "the phase's most frequent structure" : "the burned state")}{(s.Median == top ? " best chain" : "")}";
            Store(csv, $"chain{s.Index}", "1", "fixed state", fixture.Basin(state.Terms), surface, state.Terms, state.LogTau, surface.Model.Target.Evaluate(state).LogMarginal,
                  double.NaN, 0, 0, note);
        }
        _out.WriteLine($"The fixed states: {clock.Elapsed.TotalSeconds:F0} s.");
    }

    private sealed class Seen(BassTerm[] terms, string basin, int chain, int draw)
    {
        public BassTerm[] Terms { get; } = terms;
        public string Basin { get; } = basin;
        public int Chain { get; } = chain;
        public int Draw { get; } = draw;
        public long Draws { get; set; }
        public double Weight { get; set; }
        public double Pi { get; set; }
    }

    private void StoreFlatHistogram(PathCsv csv, BassMixing.Fixture fixture, PathSurface surface)
    {
        var clock = Stopwatch.StartNew();
        BassConfig config = BassFlatHistogramQualificationTests.ArmConfig(true, null, BassFlatHistogramQualificationTests.ArmOneSeed);
        BassFit fit = Bass.Fit(config, fixture.Data);
        _out.WriteLine($"Arm (i), master seed {BassFlatHistogramQualificationTests.ArmOneSeed}: {clock.Elapsed.TotalSeconds:F0} s.");
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        MulticanonicalDraws draws = fh.Draws!;
        double[] weights = fit.DrawWeights!;
        var seen = new Dictionary<string, Seen>();
        for (int c = 0; c < fit.Draws.Count; c++)
            for (int s = 0; s < draws.Draws; s++)
            {
                BassDraw d = fit.Draws[c][s];
                int[] ids = [.. d.Structure.TermIds.Order()];
                string key = string.Join(",", ids);
                if (!seen.TryGetValue(key, out Seen? at))
                {
                    BassTerm[] terms = [.. ids.Select(id => fit.Terms[id])];
                    seen[key] = at = new Seen(terms, fixture.Basin(terms), c, s);
                }
                at.Draws++;
                at.Weight += weights[c * draws.Draws + s];
            }
        _out.WriteLine($"Draws {draws.Chains * draws.Draws:N0}, structures {seen.Count:N0}.");
        foreach (var basin in seen.Values.GroupBy(v => v.Basin).OrderByDescending(g => g.Sum(v => v.Weight)).Take(12))
            _out.WriteLine($"{basin.Sum(v => v.Weight):P2} | {basin.Sum(v => v.Draws):N0} draws | {basin.Count():N0} structures | {basin.Key}");

        foreach ((string label, string basin) in Basins)
        {
            Seen[] members = [.. seen.Values.Where(v => v.Basin == basin)];
            if (members.Length == 0)
            {
                csv.Append(label, "missing", "", "", "", "", "", "flat histogram", basin, "", "", "", "", "", "", "", "", "", "", "0", "", "no draw in this run");
                _out.WriteLine($"{label}: no draw in this run.");
                continue;
            }
            Parallel.ForEach(members, new ParallelOptions { MaxDegreeOfParallelism = BassPathProfile.Threads },
                             m => m.Pi = surface.Score(m.Terms, profile: false).Pi(PathReading.Held));
            Seen[] ranked = [.. members.OrderByDescending(m => m.Pi)];
            double share = members.Sum(m => m.Weight);
            long basinDraws = members.Sum(m => m.Draws);
            string Signs(Seen m) => string.Join(" ", m.Terms.Select(BassPathProfile.Signature).Order(StringComparer.Ordinal));
            var chosen = new List<(string Rank, Seen M)> { ("1", ranked[0]) };
            if (ranked.Length > 1) chosen.Add(("2", ranked[1]));
            if (ranked.FirstOrDefault(m => Signs(m) != Signs(ranked[0])) is { } distinct) chosen.Add(("d", distinct));
            foreach ((string rank, Seen m) in chosen)
                Store(csv, label, rank, "flat histogram", basin, surface, m.Terms, Math.Log(fit.Draws[m.Chain][m.Draw].Tau), draws.Levels(m.Chain)[m.Draw], share, basinDraws, m.Draws,
                      $"structures in the basin {members.Length} pi_top {ranked[0].Pi:F2}");
        }
    }

    // ---- The stored exemplars rechecked (the scorer's check) ----

    [Fact]
    public void Scorer_TheStoredExemplars()
    {
        if (!BassQualification.Enabled) return;
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        var surface = new PathSurface(fixture.Config, fixture.Data);
        var modes = ReadModes();
        Assert.NotEmpty(modes);
        double worstDraw = 0.0, worstStored = 0.0;
        foreach (Exemplar e in modes.Values)
        {
            PathScore score = surface.Score(e.Terms);
            Assert.True(score.Admitted, $"{e.Label} {e.Rank} is refused: {score.Refusal}.");
            Assert.Equal(e.Basin, fixture.Basin(e.Terms));
            double fresh = surface.Sampled!.Target.Evaluate(score.Fit, e.DrawLogTau).LogMarginal;
            worstDraw = Math.Max(worstDraw, Math.Abs(fresh - e.DrawLevel));
            worstStored = Math.Max(worstStored, Math.Abs(score.Pi(PathReading.Held) - e.PiHeld));
            worstStored = Math.Max(worstStored, Math.Abs(score.Pi(PathReading.Profiled) - e.PiProfiled));
        }
        string stamp = DateTime.UtcNow.ToString("yyyy-MM-ddTHH:mm", CultureInfo.InvariantCulture);
        Validation().Append("scorer", $"exemplars {stamp}", F.I(modes.Count), F.F(Math.Max(worstDraw, worstStored)), F.F(worstDraw), F.F(worstStored), "", "fresh design against the retained draw and against the stored π");
        _out.WriteLine($"{modes.Count} exemplars: fresh against the retained draw's level, worst {worstDraw:E2}; against the stored π, worst {worstStored:E2} (rows hold six decimals).");
        Assert.True(worstDraw <= 2e-6 && worstStored <= 2e-6);
    }

    // ---- Part 1 ----

    [Fact]
    public void Part1_TheRoutes()
    {
        if (!BassQualification.Enabled) return;   // heavy: minutes a pair where |D| is large
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        var surface = new PathSurface(fixture.Config, fixture.Data);
        var modes = ReadModes();
        Exemplar a1 = modes[("a1", "1")];
        var pairs = new List<(string Pair, string Rank, Exemplar S, Exemplar T)>();
        void Add(string label, string name, string rank = "1", string targetRank = "1")
        {
            if (modes.TryGetValue((label, rank), out Exemplar? s) && modes.TryGetValue(("a1", targetRank), out Exemplar? t))
                pairs.Add(($"{name}-a1", $"{rank}/{targetRank}", s, t));
            else _out.WriteLine($"{label} exemplar {rank} or a1 exemplar {targetRank}: no draw in this run.");
        }
        Add("A", "A");
        Add("B", "B");
        Add("a2", "a2");
        Add("a1x", "a1+x6x8");
        for (int c = 0; c < 4; c++) Add($"chain{c}", $"chain{c}");
        Add("A", "A", "2", "2");
        Add("A", "A", "d", "1");
        Add("A", "A", "1", "d");
        Add("A", "A", "d", "d");

        PathCsv routes = Routes(), steps = Steps(), validation = Validation();
        foreach ((string pair, string rank, Exemplar s, Exemplar t) in pairs)
            ReadPair(surface, routes, steps, validation, F.I(surface.Observations), "S10", pair, rank, s.Terms, t.Terms, ["aligned", "partial", "exact"]);
        validation.Append("sc5", $"N {surface.Observations} rerun {DateTime.UtcNow:yyyy-MM-ddTHH:mm}", F.I(surface.SupersetEdges), F.F(surface.SupersetLeast), F.F(surface.HeldOccam), "", "",
                          "superset edges read; the least slack L(M+1) − (L(M) − c); no violation, or the run stops");
    }

    private void ReadPair(PathSurface surface, PathCsv routes, PathCsv steps, PathCsv validation, string n, string data, string pair, string rank, BassTerm[] s, BassTerm[] t, string[] alignments,
                          PathReading[]? at = null)
    {
        var lattices = new Dictionary<string, PathLattice>();
        foreach (string alignment in alignments)
            foreach (PathReading reading in at ?? [PathReading.Held, PathReading.Profiled])
            {
                string[] key = [n, data, pair, rank, alignment, reading.ToString()];
                if (routes.Has([.. key, "from S"]) && routes.Has([.. key, "from T"])) continue;
                var clock = Stopwatch.StartNew();
                PathRoute? route;
                try
                {
                    bool redraw = alignment.EndsWith("+redraw", StringComparison.Ordinal);
                    route = alignment.StartsWith("partial", StringComparison.Ordinal) ? BassPathProfile.BuildPartial(surface, s, t, reading, lattices, redraw)
                                                                                      : BassPathProfile.Build(surface, s, t, alignment == "exact", reading, lattices, redraw: redraw);
                }
                catch (ArgumentException e) when (e.Message.Contains("needs a beam", StringComparison.Ordinal))
                {
                    _out.WriteLine($"{pair} {rank} {alignment} {reading}: {e.Message}");
                    throw;
                }
                if (route is null)
                {
                    foreach (string direction in new[] { "from S", "from T" }) routes.Append(Blank(key, direction, "no route"));
                    _out.WriteLine($"{pair} {rank} {alignment} {reading}: no route through L.");
                    continue;
                }
                if (route.Hops.Length == 0)
                {
                    foreach (string direction in new[] { "from S", "from T" }) routes.Append(Blank(key, direction, "one structure"));
                    continue;
                }
                CheckSearch(validation, route, $"{n} {data} {pair} {rank} {alignment} {reading}");
                WriteSteps(steps, key, route);
                double seconds = clock.Elapsed.TotalSeconds;
                PathDirection forward = BassPathProfile.Read(route, true), reverse = BassPathProfile.Read(route, false);
                routes.Append(Cells(key, "from S", route, forward, route.SupersetKFrom, seconds));
                routes.Append(Cells(key, "from T", route, reverse, route.SupersetKTo, seconds));
                _out.WriteLine($"{pair} {rank} {alignment} {reading}: |D| {route.Lattice.Dim}, slide {route.Position}, {route.Hops.Count(h => h.Kind == PathHopKind.Slide)} knot steps, " +
                               $"low {forward.Low:F2}, from S k {forward.K} depth {forward.DepthFrom:F2} ln wait {forward.LnWaitClass:F1}; from T k {reverse.K} depth {reverse.DepthFrom:F2} ln wait {reverse.LnWaitClass:F1}; {seconds:F1} s");
            }
    }

    private static string[] Blank(string[] key, string direction, string note)
    {
        var row = new string[RouteHeader.Length];
        Array.Fill(row, "");
        key.CopyTo(row, 0);
        row[6] = direction;
        row[7] = note;
        return row;
    }

    private static string[] Cells(string[] key, string direction, PathRoute r, PathDirection d, int supersetK, double seconds)
    {
        int Refused(BassRefusal reason) => r.Lattice.Refused.GetValueOrDefault(reason);
        return
        [
            .. key, direction, r.Position, F.I(r.Lattice.Dim), F.I(r.Aligned.Matched.Length), F.I(r.Aligned.Dropped), F.I(r.Hops.Count(h => h.Kind == PathHopKind.Slide)),
            F.I(1L << r.Lattice.Dim), F.I(Refused(BassRefusal.Cap)), F.I(Refused(BassRefusal.Support)), F.I(Refused(BassRefusal.Numerical)), F.I(Refused(BassRefusal.VarianceInflation)),
            F.I(Refused(BassRefusal.ExactFit)), F.F(d.PiFrom), F.F(d.PiTo), F.F(d.Low), F.F(d.DepthFrom), F.F(d.DepthTo), F.I(d.K), F.I(d.KLeast), F.I(d.KLattice), F.F(d.SlideLow), F.I(d.DeletionsBelow), F.I(d.HighestM),
            F.F(d.LnWaitExact), F.F(d.LnWaitClass), F.F(d.LnWaitRouteClass), F.F(d.LnWaitSlide), F.F(d.LnWaitBest), F.F(d.WorstFall), F.F(d.WorstProposal), F.F(d.WorstClassFall), F.F(d.WorstClassProposal),
            F.F(d.TotalFall), F.F(d.TotalProposal), F.F(d.TotalProposalClass), F.F(r.LatticeMaximin), F.F(r.SupersetLow), F.F(d.PiFrom - r.SupersetLow), F.I(supersetK),
            F.F(r.Hops.Where(h => h.Kind != PathHopKind.Slide).Select(h => h.C).DefaultIfEmpty(double.NaN).First()), F.I(r.Hops.Length), F.I(r.Masks.Length - 1), F.F(seconds),
        ];
    }

    // The route's rows from S to T: a term added or deleted per row, a slide's knot steps as one row with their count, lowest
    // π and ln Σ exp(−ln c_e).
    private static void WriteSteps(PathCsv steps, string[] key, PathRoute route)
    {
        int row = 0;
        for (int i = 0; i < route.Hops.Length;)
        {
            PathHop h = route.Hops[i];
            if (h.Kind != PathHopKind.Slide)
            {
                steps.Append([.. key, F.I(row++), h.Kind.ToString().ToLowerInvariant(), h.Term, F.I(h.MFrom), F.I(h.MTo), F.F(h.DMarginal), F.F(h.DPrior), F.F(h.DPi), F.F(h.C),
                              F.F(h.PiFrom), F.F(h.PiTo), F.F(h.LnQForward), F.F(h.LnQReverse), F.F(h.LnFlow), F.F(h.LnFlowClass), F.I(h.ClassSize), "1", F.F(Math.Min(h.PiFrom, h.PiTo)), F.F(-h.LnFlow)]);
                i++;
                continue;
            }
            int j = i;
            double dm = 0.0, dp = 0.0, low = double.PositiveInfinity, resistance = double.NegativeInfinity;
            for (; j < route.Hops.Length && route.Hops[j].Kind == PathHopKind.Slide; j++)
            {
                PathHop k = route.Hops[j];
                dm += k.DMarginal;
                dp += k.DPrior;
                low = Math.Min(low, Math.Min(k.PiFrom, k.PiTo));
                resistance = LogMath.Add(resistance, -k.LnFlow);
            }
            PathHop end = route.Hops[j - 1];
            steps.Append([.. key, F.I(row++), "slide", "", F.I(h.MFrom), F.I(end.MTo), F.F(dm), F.F(dp), F.F(end.PiTo - h.PiFrom), "", F.F(h.PiFrom), F.F(end.PiTo),
                          F.F(h.LnQForward), F.F(h.LnQReverse), "", "", "0", F.I(j - i), F.F(low), F.F(resistance)]);
            i = j;
        }
    }

    // The search on this route's scores, against brute force over every ordering: L itself where |D| ≤ 7, else the cube of up to
    // three of S's and four of T's terms with the rest held as in S.
    private static void CheckSearch(PathCsv validation, PathRoute route, string name)
    {
        PathLattice lattice = route.Lattice;
        double[] pi = lattice.Pi(route.Reading);
        int start = route.Masks[0], end = route.Masks[^1], dim = lattice.Dim;
        string note = "L whole";
        if (dim > 7)
        {
            int[] fromS = [.. Enumerable.Range(0, dim).Where(i => (start >> i & 1) != 0).Take(3)];
            int[] fromT = [.. Enumerable.Range(0, dim).Where(i => (end >> i & 1) != 0).Take(4)];
            int[] bits = [.. fromS, .. fromT];
            var sub = new double[1 << bits.Length];
            for (int m = 0; m < sub.Length; m++)
            {
                int mask = start;
                for (int k = 0; k < bits.Length; k++)
                    if ((m >> k & 1) != 0) mask |= 1 << bits[k];
                    else mask &= ~(1 << bits[k]);
                sub[m] = pi[mask];
            }
            pi = sub;
            start = (1 << fromS.Length) - 1;
            end = ((1 << fromT.Length) - 1) << fromS.Length;
            dim = bits.Length;
            note = "a subcube of L";
        }
        double widest = PathSearch.Widest(pi, dim, start, end), relaxed = PathSearch.Relaxed(pi, dim, start, end), orderings = PathSearch.Orderings(pi, start, end);
        if (widest != relaxed) throw new InvalidOperationException($"The search reads {widest:R}, the relaxation {relaxed:R}: {name}.");
        if (orderings > widest) throw new InvalidOperationException($"An ordering reads {orderings:R} over the search's {widest:R}: {name}.");
        validation.Append("search", name, F.I(1L << dim), F.F(widest - orderings), F.F(widest), F.F(orderings), F.F(relaxed),
                          $"{note}; worst is widest less orderings, 0 unless the best route is not monotone");
    }

    // ---- Ground truth on enumerated fixtures ----

    [Fact]
    public void GroundTruth_OnEnumeratedFixtures()
    {
        if (!BassQualification.Enabled) return;
        var csv = new PathCsv("groundtruth", GroundHeader, 5);
        foreach (int seed in Enumerable.Range(5, 12)) GroundTruth(csv, "grid", 3, seed);
        foreach (int seed in Enumerable.Range(5, 4)) GroundTruth(csv, "grid", 4, seed);
        foreach (int seed in Enumerable.Range(5, 12)) GroundTruth(csv, "large additive", 4, seed);
        (double worstBirth, int births, double worstRedraw, int redraws) = FlowAgainstKernel();
        Validation().Append("flow", "continuous fixture at cap 3", F.I(births + redraws), F.F(Math.Max(worstBirth, worstRedraw)), F.F(worstBirth), F.F(worstRedraw), "",
                            "ln c_e against the enumerated kernel's ln π(x)K(x y) on birth edges and on redraw edges; the selection is in both so no factor is left over");
    }

    private void GroundTruth(PathCsv csv, string name, int cap, int seed)
    {
        var clock = Stopwatch.StartNew();
        // BassSweepTests' grid (N = 18, with interactions), or its five-value grid (N = 75) with main effects alone, F then small enough
        // to enumerate past three terms.
        bool large = name != "grid";
        BassConfig config = BassLocalStepTests.Continuous();
        config = config with { Structure = config.Structure with { MaxBasis = cap, MaxInteraction = large ? 1 : config.Structure.MaxInteraction } };
        BassData data = BassSweepTests.Data(large, seed);
        var surface = new PathSurface(config, data);
        var lifted = new PathSurface(config with { Structure = config.Structure with { MaxBasis = null } }, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(surface.Model);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        int n = space.States.Count;
        double[] pi = [.. space.States.Select(s => s.LogPi)];
        var whole = new List<int>[n];
        var births = new List<int>[n];
        for (int i = 0; i < n; i++) (whole[i], births[i]) = ([], []);
        foreach ((string group, var edges) in flows)
            foreach ((int i, int j) in edges.Keys)
            {
                whole[i].Add(j);
                if (group == BassKernelEnumeration.BirthDeath) births[i].Add(j);
            }
        // A local mode has no higher neighbour; among ties, which equal column spaces make, the first state stands for them.
        int[] modes = [.. Enumerable.Range(0, n).Where(i => whole[i].All(j => pi[j] < pi[i] - 1e-9 || (Math.Abs(pi[j] - pi[i]) <= 1e-9 && j > i))).OrderByDescending(i => pi[i]).Take(MaxModesPerFixture)];
        _out.WriteLine($"{name}, cap {cap}, seed {seed}: {n} states, {modes.Length} local modes read.");
        for (int a = 0; a < modes.Length; a++)
        {
            double[] overWhole = Maximin(whole, pi, modes[a]), overBirths = Maximin(births, pi, modes[a]);
            for (int b = a + 1; b < modes.Length; b++)
            {
                string[] key = [name, F.I(cap), F.I(seed), F.I(modes[a]), F.I(modes[b])];
                if (csv.Has(key)) continue;
                BassTerm[] s = space.States[modes[a]].Terms, t = space.States[modes[b]].Terms;
                PathRoute? exact = BassPathProfile.Build(surface, s, t, true, PathReading.Held, []);
                PathRoute? aligned = BassPathProfile.Build(surface, s, t, false, PathReading.Held, []);
                PathRoute? exactLifted = BassPathProfile.Build(lifted, s, t, true, PathReading.Held, []);
                PathRoute? alignedLifted = BassPathProfile.Build(lifted, s, t, false, PathReading.Held, []);
                static double Low(PathRoute? r) => r?.Low ?? double.NegativeInfinity;
                csv.Append([.. key, F.I(s.Length), F.I(t.Length), F.F(pi[modes[a]]), F.F(pi[modes[b]]), F.F(overWhole[modes[b]]), F.F(overBirths[modes[b]]), F.F(Low(exact)), F.F(Low(aligned)),
                            F.F(Low(exactLifted)), F.F(Low(alignedLifted)), F.I(exactLifted?.Lattice.Dim ?? -1), F.I(alignedLifted?.Lattice.Dim ?? -1), F.I(exact?.Lattice.Refused.GetValueOrDefault(BassRefusal.Cap) ?? -1)]);
            }
        }
        _out.WriteLine($"  {clock.Elapsed.TotalSeconds:F0} s.");
    }

    // The widest value from `start` to every state of the graph.
    private static double[] Maximin(List<int>[] graph, double[] pi, int start)
    {
        var best = new double[pi.Length];
        Array.Fill(best, double.NegativeInfinity);
        var queue = new PriorityQueue<int, double>();
        best[start] = pi[start];
        queue.Enqueue(start, -pi[start]);
        while (queue.TryDequeue(out int u, out double priority))
        {
            if (-priority < best[u]) continue;
            foreach (int v in graph[u])
            {
                double through = Math.Min(-priority, pi[v]);
                if (through > best[v])
                {
                    best[v] = through;
                    queue.Enqueue(v, -through);
                }
            }
        }
        return best;
    }

    // ---- Part 2 ----

    [Fact]
    public void Part2_AgainstN()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var modes = ReadModes();
        int[] sizes = [125, 250, 500, 1000, 2000];
        int[] seeds = SeedTree.Derive(FreshDataSeed, sizes.Length);
        BassMixing.Fixture source = BassMixing.Get(BassMixing.Friedman);
        var from = new PathSurface(source.Config, source.Data);
        // The best chain, by its top median log posterior: the note on its rows says which.
        string best = Modes().Rows().First(r => r[0].StartsWith("chain", StringComparison.Ordinal) && r[21].Contains("best chain", StringComparison.Ordinal))[0];
        _out.WriteLine($"The best chain's structure: {best}.");
        var carried = new PathCsv("carried", CarriedHeader, 4);
        PathCsv routes = Routes(), steps = Steps(), validation = Validation();
        string[] labels = ["a1", "A", "B", best];
        for (int k = 0; k < sizes.Length; k++)
            foreach ((string data, int seed) in new[] { ("fresh", seeds[k]), ("stream", 51) })
            {
                if (data == "stream" && sizes[k] == 500) continue;   // Part 1's own surface
                int n = sizes[k];
                var clock = Stopwatch.StartNew();
                BassData on = BassQualificationRecoveryTests.Friedman(n, seed).Data;
                BassConfig config = BassQualificationRecoveryTests.Config(53);
                var surface = new PathSurface(config, on);
                var fixture = new BassMixing.Fixture($"N = {n}", on, config, null);
                var polished = new Dictionary<string, BassTerm[]>();
                int[] chainSeeds = SeedTree.Derive(ShortChainSeed + k + (data == "fresh" ? 0 : 100), 2 * labels.Length);
                for (int l = 0; l < labels.Length; l++)
                {
                    if (!modes.TryGetValue((labels[l], "1"), out Exemplar? e)) continue;
                    if (Carry(carried, from, surface, fixture, n, data, e, chainSeeds[2 * l], chainSeeds[2 * l + 1]) is { } end) polished[labels[l]] = end;
                }
                foreach (string label in labels.Skip(1))
                    if (polished.TryGetValue(label, out BassTerm[]? s) && polished.TryGetValue("a1", out BassTerm[]? t))
                    {
                        string pair = label == best ? "best-a1" : $"{label}-a1", rank = label == best ? best : "1";
                        ReadPair(surface, routes, steps, validation, F.I(n), data, pair, rank, s, t, ["aligned"]);
                        if (data == "fresh") ReadPair(surface, routes, steps, validation, F.I(n), data, pair, rank, s, t, ["partial"], [PathReading.Held]);
                        if (!s.Any(x => t.Any(y => BassPathProfile.ClosedPair(surface, x, y)))) continue;
                        ReadPair(surface, routes, steps, validation, F.I(n), data, pair, rank, s, t, ["aligned+redraw"]);
                        if (data == "fresh") ReadPair(surface, routes, steps, validation, F.I(n), data, pair, rank, s, t, ["partial+redraw"], [PathReading.Held]);
                    }
                validation.Append("sc5", $"N {n} {data} rerun {DateTime.UtcNow:yyyy-MM-ddTHH:mm}", F.I(surface.SupersetEdges), F.F(surface.SupersetLeast), F.F(surface.HeldOccam), "", "",
                                  "superset edges read; the least slack L(M+1) − (L(M) − c); no violation, or the run stops");
                _out.WriteLine($"N = {n}, {data}: {clock.Elapsed.TotalSeconds:F0} s.");
            }
    }

    // The exemplar carried to the surface by knot quantile, polished with P3's greedy step at the held τ, then a short chain at
    // b = 1 from it. A term that would leave F steps inward until it stays; the row records the count.
    private BassTerm[]? Carry(PathCsv carried, PathSurface from, PathSurface to, BassMixing.Fixture fixture, int n, string data, Exemplar e, int moveSeed, int blockSeed)
    {
        var terms = new List<BassTerm>();
        int nudges = 0;
        foreach (BassTerm t in e.Terms)
        {
            HingeFactor[] hinges = [.. t.Continuous.Select(h =>
            {
                int was = BassMixing.Knots(from.Model, BassPart.Continuous, h.Variable), now = BassMixing.Knots(to.Model, BassPart.Continuous, h.Variable);
                return h with { Knot = (int)Math.Round((double)h.Knot / Math.Max(1, was - 1) * (now - 1)) };
            })];
            BassTerm? placed = null;
            foreach (int shift in Enumerable.Range(0, 33).Select(d => (d + 1) / 2 * (d % 2 == 0 ? 1 : -1)))
            {
                foreach (int f in Enumerable.Range(0, hinges.Length))
                {
                    int knots = BassMixing.Knots(to.Model, BassPart.Continuous, hinges[f].Variable), knot = hinges[f].Knot + shift;
                    if (knot < 0 || knot >= knots) continue;
                    HingeFactor[] moved = [.. hinges];
                    moved[f] = moved[f] with { Knot = knot };
                    var candidate = new BassTerm(moved, null, null);
                    if (terms.Contains(candidate) || !to.Score([.. terms, candidate], profile: false).Admitted) continue;
                    placed = candidate;
                    if (shift != 0) nudges++;
                    break;
                }
                if (placed is not null) break;
            }
            if (placed is null)
            {
                _out.WriteLine($"N = {n}, {data}: {e.Label}'s term {BassMixing.Key([t])} found no place.");
                return null;
            }
            terms.Add(placed);
        }
        double u = Math.Log(to.HeldTau);
        BassState start = BassMixing.Loaded(to.Held, terms, u, out BassRefusal? refusal) ?? throw new InvalidOperationException($"The carried structure left F ({refusal}).");
        double piCarried = to.Score(terms).Pi(PathReading.Held);
        BassMixingExitProfileTests.Search search = BassMixingExitProfileTests.Greedy(to.Held, start, 50);
        BassTerm[] end = [.. search.End.Terms];
        PathScore score = to.Score(end);
        string basin = fixture.Basin(end);
        bool kept = basin == e.Basin;

        const int sweeps = BassMixing.Samples;
        var chain = new BassChain(to.Model, end, score.BestLogTau, new Xoshiro256PlusPlus(moveSeed), new Xoshiro256PlusPlus(blockSeed));
        double startPosterior = BassMixing.LogPosterior(to.Model, chain.Current);
        var posterior = new double[sweeps];
        var basins = new Dictionary<string, int>();
        string last = basin;
        for (int s = 0; s < sweeps; s++)
        {
            BassState state = chain.Step();
            posterior[s] = BassMixing.LogPosterior(to.Model, state);
            last = fixture.Basin(state.Terms);
            basins[last] = basins.GetValueOrDefault(last) + 1;
        }
        KeyValuePair<string, int> top = basins.MaxBy(b => b.Value);
        double rank = posterior.Count(p => p <= startPosterior) / (double)sweeps;
        double gap = startPosterior - BassMixing.Quantile(posterior, 0.5);
        carried.Append(F.I(n), data, e.Label, e.Rank, F.I(end.Length), F.I(nudges), F.F(piCarried), F.F(score.Pi(PathReading.Held)), F.F(score.Pi(PathReading.Profiled)), F.F(score.BestLogTau),
                       F.I(search.Passes), search.AtLimit ? "1" : "0", F.I(search.Displacement), kept ? "1" : "0", F.I(sweeps), F.F(basins.GetValueOrDefault(basin) / (double)sweeps), last, top.Key,
                       F.F(top.Value / (double)sweeps), F.F(rank), F.F(gap), BassMixing.Key(end));
        return end;
    }
}
