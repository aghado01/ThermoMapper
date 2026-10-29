using System;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's P3 (the mode-exit pricing brief): profile deaths, an optimization diagnostic of the design draft's drift-path
/// hypothesis. From each chain's fixed state (<see cref="BassMixing.FixedStates"/>), each term is removed and the remaining
/// hinges' knots are re-optimized greedily at fixed τ: in passes over the factors, each tries its knot j ± δ for δ ∈ {1, …,
/// 8}, sign kept, and takes the best step that raises the log marginal plus the log structure prior (the knot prior is
/// uniform, so the marginal); categorical subsets stay fixed. The search stops after 5 passes or a pass without gain; where
/// it stopped at the limit, it is run on to 50 passes, so the reading says whether the limit or a local optimum ended it.
/// Each term's plain death Δℓ, its profile Δℓ against the fixed state, the passes and the knots' total displacement in
/// distinct-value steps. A control beside them: the same search on the full state, no term removed, since a sampled state
/// is not its knots' optimum; each profile death is read against it too. It bounds where a drift could end, not that a
/// path exists or how long it takes. It reads, it does
/// not judge; it asserts its own soundness: nothing committed (<see cref="BassMixing.Untouched"/>), each accepted step
/// raising the objective, and the tracked objective matching a fresh load of the end state.
/// </summary>
public sealed class BassMixingExitProfileTests
{
    private const int Passes = 5, Extended = 50, Reach = 8;
    private const double Near = -10.0, Far = -20.0;
    private readonly ITestOutputHelper _out;
    public BassMixingExitProfileTests(ITestOutputHelper output) => _out = output;

    // A greedy search's end: the state, its passes, whether its last pass still gained, its knots' total displacement.
    private sealed record Search(BassState End, int Passes, bool AtLimit, int Displacement, int Steps);

    // One term's profile death: its plain death's Δℓ, the 5-pass search, and the run on where that stopped at its limit.
    private sealed record Profiled(BassTerm Term, double Plain, Search Profile, double ProfileChange, Search? Longer, double LongerChange)
    {
        // The search's end, run on where the 5 passes stopped it at the limit.
        public double Converged => Longer is null ? ProfileChange : LongerChange;
    }

    // A state's profile deaths, and the control: the same search on the full state, 5 passes and run on.
    private sealed record StateProfile(Profiled[] Terms, double Full, double FullConverged, int FullPasses, bool FullAtLimit);

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void ProfileDeaths_FromFixedStates(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassMixing.Fixture fixture = BassMixing.Get(name);
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassMixing.FixedState[] states = BassMixing.FixedStates(model);
        double setup = clock.Elapsed.TotalSeconds;
        clock.Restart();
        var profiles = new StateProfile[states.Length];
        Parallel.For(0, states.Length, c =>
        {
            BassChain chain = states[c].Chain;
            var untouched = new BassMixing.Untouched(model, chain);
            profiles[c] = Profiles(model, chain.Current);
            untouched.Check();
        });
        _out.WriteLine($"## P3, {name}: profile deaths from each chain's fixed state (states and phases {setup:F1} s, searches {clock.Elapsed.TotalSeconds:F1} s)");
        _out.WriteLine("");
        foreach (string line in BassMixing.FixedTable(model, states)) _out.WriteLine(line);
        Summary(states, profiles);
        for (int c = 0; c < states.Length; c++) Terms(fixture, model, states[c], profiles[c]);
    }

    private static StateProfile Profiles(BassModel model, BassState state)
    {
        double ell = model.Target.Evaluate(state).LogMarginal;
        Search full = Greedy(model, state, Passes), fullLonger = full.AtLimit ? Greedy(model, state, Extended) : full;
        var profiles = new Profiled[state.Count];
        for (int i = 0; i < state.Count; i++)
        {
            BassTerm[] rest = [.. state.Terms.Where((_, j) => j != i)];
            BassState reduced = BassMixing.Loaded(model, rest, state.LogTau, out BassRefusal? refusal) ?? throw new InvalidOperationException($"A death left F ({refusal}).");
            double plain = model.Target.Evaluate(reduced).LogMarginal - ell;
            Search profile = Greedy(model, reduced, Passes);
            Search? longer = profile.AtLimit ? Greedy(model, reduced, Extended) : null;
            profiles[i] = new Profiled(state.Terms[i], plain, profile, model.Target.Evaluate(profile.End).LogMarginal - ell,
                                       longer, longer is null ? double.NaN : model.Target.Evaluate(longer.End).LogMarginal - ell);
            Assert.True(profiles[i].ProfileChange >= plain - 1e-9 * Math.Abs(ell), "A greedy search never lowers its objective.");
        }
        return new StateProfile(profiles, model.Target.Evaluate(full.End).LogMarginal - ell, model.Target.Evaluate(fullLonger.End).LogMarginal - ell,
                                fullLonger.Passes, fullLonger.AtLimit);
    }

    // The greedy knot search from `start`, at most `limit` passes, through one design's proposals; the end loaded afresh.
    private static Search Greedy(BassModel model, BassState start, int limit)
    {
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(start.Terms));
        BassState current = start;
        BassTerm[] slots = [.. start.Terms];   // each term as it moves, in the start's order
        double objective = Objective(model, current);
        int passes = 0, steps = 0;
        bool gained = true;
        while (gained && passes < limit)
        {
            gained = false;
            passes++;
            for (int s = 0; s < slots.Length; s++)
                foreach (BassMixing.FactorAt f in BassMixing.Factors(slots[s]).Where(f => f.Part != BassPart.Categorical).ToArray())
                {
                    BassTerm? best = null;
                    double bestObjective = objective;
                    foreach (int delta in Enumerable.Range(1, Reach).SelectMany(d => new[] { d, -d }))
                    {
                        if (BassMixing.Stepped(model, slots[s], f, delta) is not { } moved) continue;
                        (BassMixing.Exit exit, BassState? candidate) = BassMixing.Replace(model, design, current, Index(current, slots[s]), moved, 0.0);
                        if (candidate is null) continue;
                        double next = objective + exit.LogLikelihoodChange + exit.LogPriorRatio;
                        if (next > bestObjective + 1e-9) (best, bestObjective) = (moved, next);
                    }
                    if (best is null) continue;
                    (_, BassState? taken) = BassMixing.Replace(model, design, current, Index(current, slots[s]), best, 0.0);
                    design.Accept();
                    current = taken!;
                    slots[s] = best;
                    objective = bestObjective;
                    gained = true;
                    steps++;
                }
        }
        BassState end = BassMixing.Loaded(model, current.Terms, start.LogTau, out BassRefusal? refusal) ?? throw new InvalidOperationException($"The search left F ({refusal}).");
        double fresh = Objective(model, end);
        Assert.True(Math.Abs(fresh - objective) <= 1e-7 * Math.Abs(fresh), $"The tracked objective {objective:R} against a fresh load's {fresh:R}");
        int displacement = 0;
        for (int s = 0; s < slots.Length; s++)
            displacement += Hinges(slots[s]).Zip(Hinges(start.Terms[s]), (a, b) => Math.Abs(a.Knot - b.Knot)).Sum();
        return new Search(end, passes, gained && passes == limit, displacement, steps);
    }

    private static System.Collections.Generic.IEnumerable<HingeFactor> Hinges(BassTerm t) => t.Continuous.Concat(t.Functional);

    private static int Index(BassState state, BassTerm term)
    {
        int i = Array.IndexOf(state.TermArray, term);
        Assert.True(i >= 0, "The search lost a term.");
        return i;
    }

    // The search's objective: the log marginal plus the log structure prior.
    private static double Objective(BassModel model, BassState state) => BassMixing.LogPosterior(model, state);

    // Per chain: the plain and profile deaths near the state, and the searches' passes and stops.
    private void Summary(BassMixing.FixedState[] states, StateProfile[] profiles)
    {
        _out.WriteLine("");
        _out.WriteLine($"chain | M | best plain Δℓ | best profile Δℓ | plain > {Near}, > {Far} | profile > {Near}, > {Far} | median gain (profile − plain) | largest gain | " +
                       $"passes, median (max) | stopped at {Passes} passes | run on to {Extended}: best Δℓ, profile > {Near}, > {Far}, still at the limit");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < states.Length; c++)
        {
            Profiled[] p = profiles[c].Terms;
            Profiled[] limited = [.. p.Where(x => x.Longer is not null)];
            double[] best = [.. p.Select(x => x.Longer is null ? x.ProfileChange : x.LongerChange)];
            string longer = limited.Length == 0 ? "—"
                : $"{F(best.Max())}, {best.Count(b => b > Near)}, {best.Count(b => b > Far)}, {limited.Count(x => x.Longer!.AtLimit)} of {limited.Length}";
            _out.WriteLine($"{states[c].Label} | {p.Length} | {F(p.Max(x => x.Plain))} | {F(p.Max(x => x.ProfileChange))} | {p.Count(x => x.Plain > Near)}, {p.Count(x => x.Plain > Far)} | " +
                           $"{p.Count(x => x.ProfileChange > Near)}, {p.Count(x => x.ProfileChange > Far)} | {F(BassMixing.Quantile(p.Select(x => x.ProfileChange - x.Plain), 0.5))} | " +
                           $"{F(p.Max(x => x.ProfileChange - x.Plain))} | {BassMixing.Quantile(p.Select(x => (double)x.Profile.Passes), 0.5)} ({p.Max(x => x.Profile.Passes)}) | " +
                           $"{limited.Length} | {longer}");
        }

        _out.WriteLine("");
        _out.WriteLine($"The control, chain | the full state's own profile Δℓ, {Passes} passes | run on (passes{(profiles.Any(p => p.FullAtLimit) ? ", at the limit" : "")}) | " +
                       $"best profile death less the full profile, {Passes} passes | run on | deaths within 10, 20 of the full profile, run on");
        _out.WriteLine("--- | --- | --- | --- | --- | ---");
        for (int c = 0; c < states.Length; c++)
        {
            StateProfile s = profiles[c];
            double[] against = [.. s.Terms.Select(x => x.Converged - s.FullConverged)];
            _out.WriteLine($"{states[c].Label} | {F(s.Full)} | {F(s.FullConverged)} ({s.FullPasses}{(s.FullAtLimit ? ", at the limit" : "")}) | " +
                           $"{F(s.Terms.Max(x => x.ProfileChange) - s.Full)} | {F(against.Max())} | {against.Count(a => a > Near)}, {against.Count(a => a > Far)}");
        }
    }

    // Per chain, each term: plain against profile, the search's passes, steps and displacement.
    private void Terms(BassMixing.Fixture fixture, BassModel model, BassMixing.FixedState state, StateProfile profile)
    {
        _out.WriteLine("");
        _out.WriteLine($"Chain {state.Label}, M {profile.Terms.Length}: term removed | plain Δℓ | profile Δℓ | gain | passes | steps | displacement | " +
                       $"run on to {Extended}: Δℓ, passes, displacement | run on, less the full profile ({F(profile.FullConverged)})");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | ---");
        foreach (Profiled p in profile.Terms.OrderByDescending(p => p.Converged))
            _out.WriteLine($"{fixture.Describe(model, p.Term)} | {F(p.Plain)} | {F(p.ProfileChange)} | {F(p.ProfileChange - p.Plain)} | {p.Profile.Passes} | {p.Profile.Steps} | " +
                           $"{p.Profile.Displacement} | {(p.Longer is { } l ? $"{F(p.LongerChange)}, {l.Passes}{(l.AtLimit ? " (at the limit)" : "")}, {l.Displacement}" : "—")} | " +
                           $"{F(p.Converged - profile.FullConverged)}");
    }

    private static string F(double v) => double.IsNaN(v) ? "—" : v.ToString("F2");
}
