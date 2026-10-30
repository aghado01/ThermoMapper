using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Maths.Rng;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// Exploration, not a claim (N17 Phase 1's mean-field addition): how each sampler's crossing between the two phases of a
/// mean-field Potts model's first-order transition grows with N. Runs only when <c>FLAT_HISTOGRAM_MEANFIELD_READING</c>
/// is set (1 for the full sweep; a number below 1 scales every budget down, for a calibration). Three samplers at equal
/// single-site proposals, every one on the lumped occupation chain:
/// <list type="bullet">
/// <item>the flat histogram, production on SAMC's exact limit θ* from the exact density of states (a learned run's end
/// point, so the learning's own cost is not in the reading), four chains, in two windows:
/// <list type="bullet">
/// <item>the engine's own, from the b_t law's 0.1% quantile to the b = 1 law's 99.9%, b_lo = b_t, as a run that reweighs to
/// every b in [b_t, 1] takes it;</item>
/// <item>the transition's, Berg and Neuhaus's setup at b_t: the same engine on the target whose likelihood is b_t·ℓ₁, so
/// both extensions have slope b_t, over the b_t law's 0.1% to 99.9% quantiles, the two phases and little more;</item>
/// </list></item>
/// <item>parallel tempering on the likelihood path from b = 1 to the prior, sized in its burn-in by the existing
/// self-sizing rule (<see cref="TemperingLadder.Size"/>, <see cref="TemperingLadder.Place"/>, <see cref="TemperingLadder.Sources"/>),
/// its proposals counted across rungs. Its replicas are the states, each tagged so it is followed through the swaps:
/// a replica's round trip is a state that changes phase and changes back. A rung's own sequence of ℓ₁ also flips
/// whenever states of the two phases swap past it, without any state crossing, so the rungs' count is kept beside it
/// as a detail, not as the measure;</item>
/// <item>four fixed-b chains at b_t, two started disordered and two ordered.</item>
/// </list>
/// One crossing measure: a round trip of a chain's (or a replica's) ℓ₁ from at most ℓ_dis to at least ℓ_ord and back,
/// the two peaks of the exact law at b_t, per 10⁶ production proposals.
/// <para><b>Written before the run.</b> H_exp, for tempering and the fixed-b chain: crossing needs the barrier the exact
/// law shows, so r(N) ≈ A·N^{−a}·(dip/peak)(N), exponential in N because ln(peak/dip) grows linearly in N (about 0.0011
/// nats per site at q = 3 and 0.0065 at q = 4 by the limit's Δφ). H_pow, for the flat histogram: the walk in ℓ₁ is
/// flat across the window, so crossing costs a diffusion over the order parameter, r(N) ≈ A·N^{−α} with α near 2 and in
/// [1, 3]. Each sampler's ln r is fitted by both forms, a + b·N and a − α·ln N, over the sizes with crossings. A reading
/// distinguishes them for a sampler when one form's RMS residual is under 0.3 and under half the other's, or when a
/// power law fitted to the smaller sizes predicts more than three round trips where none came (Poisson P(0) &lt; 5%). At
/// q = 3 the barrier reaches only about e^{1.4} by N = 2048, so the forms may not separate there; q = 4 to N = 2048
/// reaches about e^{12} and should. The budget is 100·N² proposals per sampler (at least 2·10⁶), set from a calibration
/// at N ≤ 256, where the engine-window flat histogram's rate fell as about N^{−2.4}, so that it sees a few round trips
/// even at the largest size.</para>
/// </summary>
public sealed class MeanFieldTunnellingReading
{
    private readonly ITestOutputHelper _out;
    public MeanFieldTunnellingReading(ITestOutputHelper output) => _out = output;

    internal const string Variable = "FLAT_HISTOGRAM_MEANFIELD_READING";
    private const int Chains = 4;

    internal static double? Scale =>
        double.TryParse(Environment.GetEnvironmentVariable(Variable), NumberStyles.Float, CultureInfo.InvariantCulture, out double s) && s > 0 ? Math.Min(s, 1.0) : null;

    [Fact]
    public void Reading_CrossingRatesAgainstN()
    {
        if (Scale is not double scale) return;   // exploration; set FLAT_HISTOGRAM_MEANFIELD_READING to run
        var table = new StringBuilder("q\tN\tb_t\tln_peak_dip\tl_dis\tl_ord\tsampler\tproposals\trate_per_1e6\tspread\tround_trips\tdetail\n");
        var rates = new Dictionary<(int Q, string Sampler), List<(int N, double Rate, long Trips, double Proposals)>>();
        foreach ((int q, int[] sizes) in new[] { (3, new[] { 64, 128, 256, 512, 1024, 2048 }), (4, new[] { 64, 128, 256, 512, 1024, 2048 }) })
            foreach (int n in sizes)
            {
                Size size = Read(q, n, scale);
                _out.WriteLine($"[meanfield] q = {q}, N = {n}: b_t {size.Bt:F4} (b_tJ {size.Bt * size.Coupling:F4} against {MeanFieldPotts.CriticalCoupling(q):F4}), " +
                               $"ln(peak/dip) {size.LogPeakToDip:F2}, ℓ_dis {size.Disordered:F1}, ℓ_ord {size.Ordered:F1}, budget {size.Budget:E2} proposals each");
                foreach (Sampler s in size.Samplers)
                {
                    _out.WriteLine($"[meanfield]   {s.Name,-14} {s.Rate,10:F3} per 10⁶ (spread {s.Spread:F3}, {s.Trips} round trips over {s.Units}) {s.Detail}");
                    table.Append($"{q}\t{n}\t{size.Bt:R}\t{size.LogPeakToDip:R}\t{size.Disordered:R}\t{size.Ordered:R}\t{s.Name}\t{s.Proposals:R}\t{s.Rate:R}\t{s.Spread:R}\t{s.Trips}\t{s.Detail}\n");
                    if (!rates.TryGetValue((q, s.Name), out var list)) rates[(q, s.Name)] = list = [];
                    list.Add((n, s.Rate, s.Trips, s.Proposals));
                }
            }
        foreach (((int q, string sampler), var list) in rates)
            _out.WriteLine($"[meanfield] fit q = {q}, {sampler}: {Fit(list)}");
        _out.WriteLine("[meanfield] table\n" + table);
    }

    internal sealed record Sampler(string Name, double Proposals, long Trips, double Rate, double Spread, string Units, string Detail);
    internal sealed record Size(double Coupling, double Bt, double LogPeakToDip, double Disordered, double Ordered, double Budget, Sampler[] Samplers);

    // Budget per sampler: 100·N² proposals, at least 2·10⁶, scaled.
    internal static double Budget(int n, double scale) => scale * Math.Max(2e6, 100.0 * n * n);

    internal static Size Read(int q, int n, double scale)
    {
        double coupling = MeanFieldPotts.CriticalCoupling(q) / 0.6;
        var model = new MeanFieldPotts(n, q, coupling);
        var exact = new MeanFieldExact(model);
        double split = exact.Split();
        double bt = exact.EqualWeight(split);
        (int dis, int dip, int ord) = exact.Phases(bt, split);
        double[] smooth = exact.SmoothedLogLaw(bt);
        double lnRatio = Math.Min(smooth[dis], smooth[ord]) - smooth[dip];
        double lDis = exact.BinCentre(dis), lOrd = exact.BinCentre(ord);
        double budget = Budget(n, scale);
        int[] seeds = SeedTree.Derive(SeedTree.Derive(q, n + 1)[n], 4);
        var samplers = new Sampler[4];
        Parallel.Invoke(
            () => samplers[0] = FlatHistogram(model, exact, bt, split, lDis, lOrd, budget, seeds[0], transition: false),
            () => samplers[1] = FlatHistogram(model, exact, bt, split, lDis, lOrd, budget, seeds[1], transition: true),
            () => samplers[2] = Tempering(model, lDis, lOrd, budget, seeds[2]),
            () => samplers[3] = FixedB(model, bt, lDis, lOrd, budget, seeds[3]));
        return new Size(coupling, bt, lnRatio, lDis, lOrd, budget, samplers);
    }

    // Production on θ* over the window: four chains from the disordered centre, a tenth of the budget's worth of burn-in
    // each, then the budget shared among them. The transition's window runs the engine on b_t·ℓ₁, in whose units the
    // extensions' slopes b_lo = 1 and 1 are b_t in ℓ₁'s.
    private static Sampler FlatHistogram(MeanFieldPotts model, MeanFieldExact exact, double bt, double split, double lDis, double lOrd,
                                         double budget, int seed, bool transition)
    {
        double scale = transition ? bt : 1.0;
        double low = exact.Quantile(bt, 0.001), high = transition ? exact.Quantile(bt, 0.999) : exact.Quantile(1.0, 0.999);
        // Nodes a disordered-phase sd apart, so the segments resolve the sharpest part of ln g.
        double[] p = exact.Law(bt);
        double mass = 0.0, mean = 0.0, second = 0.0;
        for (int j = 0; j < p.Length; j++)
            if (exact.Levels[j] <= split) { mass += p[j]; mean += p[j] * exact.Levels[j]; second += p[j] * exact.Levels[j] * exact.Levels[j]; }
        mean /= mass;
        double sdDis = Math.Sqrt(Math.Max(second / mass - mean * mean, 1e-12));
        int nodes = Math.Max(3, (int)Math.Ceiling((high - low) / sdDis) + 1);
        double[] levels = Array.ConvertAll(exact.Levels, l => scale * l);
        MulticanonicalWeights shape = MulticanonicalWeights.Flat(scale * low, scale * high, transition ? 1.0 : bt, nodes);
        MulticanonicalWeights weights = SamcLearning.FixedPoint(shape, levels, exact.LogMasses, out _);
        ISamplingTarget<int[]> inner = transition ? SamplingTarget.Separable<int[]>(model.LogPrior, x => bt * model.Level(x)) : model.Target;
        var target = new MulticanonicalTarget<int[]>(inner, weights);
        int[] seeds = SeedTree.Derive(seed, Chains);
        long steps = (long)(budget / Chains);
        var trips = new long[Chains];
        Parallel.For(0, Chains, c =>
        {
            var chain = new ReversibleJumpChain<int[]>([new OccupationMove(model)], target, model.Disordered(), new Xoshiro256PlusPlus(seeds[c]));
            for (long i = 0; i < steps / 10; i++) chain.Step();
            var counter = new Crossings(lDis, lOrd);
            for (long i = 0; i < steps; i++)
            {
                chain.Step();
                counter.Add(((MulticanonicalEvaluation)chain.CurrentEvaluation).UntemperedLogLikelihood / scale);
            }
            trips[c] = counter.RoundTrips;
        });
        return Summarize(transition ? "flat at b_t" : "flat histogram", trips, steps, Chains * (double)steps, "chains", $"K = {nodes}, window [{low:F1}, {high:F1}]");
    }

    // Parallel tempering from b = 1 to the prior, self-sized in burn-in rounds as BASS's ladder is, then frozen. Each
    // state carries a tag, kept by its local moves and swapped with it, so the replicas are followed through the ladder.
    private static Sampler Tempering(MeanFieldPotts model, double lDis, double lOrd, double budget, int seed)
    {
        double[] ladder = Enumerable.Range(0, 8).Select(i => 1.0 - i / 7.0).ToArray();
        IReadOnlyList<IRjMove<Tagged>> Moves(int _) => [new TaggedMove(model)];
        ISamplingTarget<Tagged> target = SamplingTarget.Separable<Tagged>(t => model.LogPrior(t.Counts), t => model.Level(t.Counts));
        int nextTag = 0;
        var pt = new ParallelTempering<Tagged>(Moves, target, new Tagged(nextTag++, model.Disordered()), ladder, TemperingPath.Likelihood, seed);
        for (int i = 1; i < pt.Replicas.Count; i++)
        {
            var own = new Tagged(nextTag++, model.Disordered());
            pt.Replicas[i].SetState(own, target.Evaluate(own));
        }
        // Burn-in: a quarter of the budget, in eight rounds, each resizing on its swap rejections; a new replica forks its
        // source's state under a new tag.
        double burnProposals = budget / 4.0;
        for (int round = 0; round < 8; round++)
        {
            long scans = Math.Max(100, (long)(burnProposals / 8 / pt.Replicas.Count));
            for (long s = 0; s < scans; s++) pt.Step();
            double[] rejection = pt.RoundRejection;
            int current = pt.Replicas.Count;
            int count = round == 7 ? current : TemperingLadder.Size(TemperingLadder.Barrier(rejection), current, 2, 128);
            if (count == current)
            {
                pt.AdaptLadder();
                continue;
            }
            double[] placed = TemperingLadder.Place(pt.InverseTemperatures, rejection, count);
            int[] sources = TemperingLadder.Sources(pt.InverseTemperatures, rejection, placed);
            pt.Resize(placed, sources, x => new Tagged(nextTag++, (int[])x.Counts.Clone()), Moves);
        }
        int rungs = pt.Replicas.Count;
        long production = (long)(budget / rungs);
        var replicas = new Dictionary<int, Crossings>();
        foreach (ReversibleJumpChain<Tagged> r in pt.Replicas) replicas[r.Current.Tag] = new Crossings(lDis, lOrd);
        var slots = new Crossings[rungs];
        for (int i = 0; i < rungs; i++) slots[i] = new Crossings(lDis, lOrd);
        long ladderTrips = pt.RoundTrips;
        for (long s = 0; s < production; s++)
        {
            pt.Step();
            for (int i = 0; i < rungs; i++)
            {
                ReversibleJumpChain<Tagged> rung = pt.Replicas[i];
                double level = rung.CurrentEvaluation.LogLikelihood(1.0);
                replicas[rung.Current.Tag].Add(level);
                slots[i].Add(level);
            }
        }
        long[] trips = replicas.Values.Select(c => c.RoundTrips).ToArray();
        double swap = pt.SwapAttempts > 0 ? (double)pt.SwapAccepts / pt.SwapAttempts : double.NaN;
        return Summarize("tempering", trips, production, rungs * (double)production, "replicas",
                         $"swap acceptance {swap:F2}, {pt.RoundTrips - ladderTrips} ladder round trips (hot end to cold and back) in all, " +
                         $"{slots.Sum(c => c.RoundTrips)} rung-sequence round trips (not the measure), ladder " +
                         string.Join(" ", pt.InverseTemperatures.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))));
    }

    // A mean-field state with a tag that its local moves keep and its swaps carry, so a replica can be followed.
    internal sealed record Tagged(int Tag, int[] Counts);

    private sealed class TaggedMove(MeanFieldPotts model) : IRjMove<Tagged>
    {
        private readonly OccupationMove _move = new(model);
        public string Key => _move.Key;
        public string ReverseKey => _move.ReverseKey;
        public double Weight(Tagged state) => 1.0;

        public Proposal<Tagged>? Propose(Tagged current, Xoshiro256PlusPlus rng)
        {
            Proposal<int[]> p = _move.Propose(current.Counts, rng)!.Value;
            return new Proposal<Tagged>(current with { Counts = p.Candidate }, p.LogProposalRatio, p.LogJacobian);
        }
    }

    // Four chains at b_t on the likelihood path, two from the disordered centre and two from the ordered point at s_c.
    private static Sampler FixedB(MeanFieldPotts model, double bt, double lDis, double lOrd, double budget, int seed)
    {
        int q = model.Q, n = model.Sites;
        double sc = (q - 2.0) / (q - 1.0);
        var ordered = new int[q];
        ordered[0] = (int)Math.Round(n * (1.0 + (q - 1) * sc) / q);
        for (int k = 1; k < q; k++) ordered[k] = (n - ordered[0]) / (q - 1) + (k <= (n - ordered[0]) % (q - 1) ? 1 : 0);
        int[] seeds = SeedTree.Derive(seed, Chains);
        long steps = (long)(budget / Chains);
        var trips = new long[Chains];
        Parallel.For(0, Chains, c =>
        {
            var chain = new ReversibleJumpChain<int[]>([new OccupationMove(model)], model.Target, c < 2 ? model.Disordered() : ordered,
                                                       new Xoshiro256PlusPlus(seeds[c]), new Temperature(bt, TemperingPath.Likelihood));
            var counter = new Crossings(lDis, lOrd);
            for (long i = 0; i < steps; i++)
            {
                chain.Step();
                counter.Add(chain.CurrentEvaluation.LogLikelihood(1.0));
            }
            trips[c] = counter.RoundTrips;
        });
        return Summarize("fixed b_t", trips, steps, Chains * (double)steps, "chains", "");
    }

    // The rate per 10⁶ proposals over every unit, and its spread: the sd of the units' own rates, each per 10⁶ of its
    // steps, scaled to the pooled rate's units by the unit count.
    private static Sampler Summarize(string name, long[] trips, long stepsPerUnit, double proposals, string unit, string detail)
    {
        long total = trips.Sum();
        double[] perUnit = trips.Select(t => 1e6 * t / (double)stepsPerUnit / trips.Length).ToArray();
        double mean = perUnit.Average();
        double spread = perUnit.Length > 1 ? Math.Sqrt(perUnit.Sum(r => (r - mean) * (r - mean)) / (perUnit.Length - 1)) : double.NaN;
        return new Sampler(name, proposals, total, 1e6 * total / proposals, spread, $"{trips.Length} {unit}", detail);
    }

    // Both growth forms fitted to ln r over the sizes with crossings, by least squares, with their RMS residuals; sizes
    // with none are listed with the power law's prediction there.
    internal static string Fit(List<(int N, double Rate, long Trips, double Proposals)> points)
    {
        var used = points.Where(p => p.Trips > 0).ToList();
        if (used.Count < 3) return $"too few sizes with crossings ({used.Count}) to fit; none at N = {string.Join(", ", points.Where(p => p.Trips == 0).Select(p => p.N))}";
        double[] y = used.Select(p => Math.Log(p.Rate)).ToArray();
        (double a1, double b1, double rmsExp) = LeastSquares(used.Select(p => (double)p.N).ToArray(), y);
        (double a2, double b2, double rmsPow) = LeastSquares(used.Select(p => Math.Log(p.N)).ToArray(), y);
        var censored = new StringBuilder();
        foreach (var p in points.Where(p => p.Trips == 0))
        {
            double predictedTrips = Math.Exp(a2 + b2 * Math.Log(p.N)) * p.Proposals / 1e6;
            censored.Append($"; none at N = {p.N}, where the power law predicts {predictedTrips:F1}");
        }
        return $"exponential ln r = {a1:F2} {b1:+0.00000;-0.00000}·N, RMS {rmsExp:F3}; power ln r = {a2:F2} {b2:+0.000;-0.000}·ln N, RMS {rmsPow:F3}{censored}";
    }

    private static (double A, double B, double Rms) LeastSquares(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        double sxy = 0.0, sxx = 0.0;
        for (int i = 0; i < x.Length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) * (x[i] - mx); }
        double b = sxy / sxx, a = my - b * mx;
        double ss = 0.0;
        for (int i = 0; i < x.Length; i++) ss += Math.Pow(y[i] - a - b * x[i], 2);
        return (a, b, Math.Sqrt(ss / x.Length));
    }

    /// <summary>Round trips of one sequence of ℓ₁ between two levels: from at most the low one to at least the high one
    /// and back, counted on the return; the first visit to the low level starts the count.</summary>
    internal sealed class Crossings(double low, double high)
    {
        private sbyte _last;   // −1 last at the low level, +1 at the high one, 0 neither yet
        public long RoundTrips { get; private set; }

        public void Add(double level)
        {
            if (level <= low)
            {
                if (_last == +1) RoundTrips++;
                _last = -1;
            }
            else if (level >= high && _last == -1) _last = +1;
        }
    }
}
