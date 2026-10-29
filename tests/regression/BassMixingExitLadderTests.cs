using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's P4 (the mode-exit pricing brief): the structure moves along M5's dense ladder on S10's tempered Friedman fit.
/// Per chain, the coarse ladder is self-sized through its 5,000 scans of burn-in as the run sizes it; the dense ladder takes
/// its rungs above 0.40, 31 rungs from 0.40 to 0.10 at Δb = 0.01, and its rungs below 0.10 (M5's Part 2), and is burned
/// 5,000 scans with adaptation off. At each rung's state, scored under that rung's tempered target (the ladder's
/// likelihood path, the integrated form, through the target's evaluation at b), none committed: 2,000 births, every term's
/// death and 2,000 change redraws through the replica's own moves (<see cref="BassMixing.Score"/> at the replica's
/// temperature); the local knot steps δ ∈ {±1, ±2, ±3} of every hinge, sign kept; and P2's prunes. Per rung (b, M), each
/// kind's mean α, and the deaths' Δℓ, tempered (Δℓ_b) and untempered (Δℓ₁, and b·Δℓ₁). It reads, it does not judge; it
/// asserts its own soundness: nothing committed per rung (<see cref="BassMixing.Untouched"/>), each rung scored at its
/// replica's own temperature, and the counts.
/// </summary>
public sealed class BassMixingExitLadderTests
{
    private const int Births = 2_000, Changes = 2_000;
    private const double Top = 0.40, Bottom = 0.10, Step = 0.01, Negligible = 1e-3;
    private static readonly int[] KnotSteps = [1, -1, 2, -2, 3, -3];
    private readonly ITestOutputHelper _out;
    public BassMixingExitLadderTests(ITestOutputHelper output) => _out = output;

    // One rung's reading: its b and M; each kind's acceptances (a step's by |δ|); the deaths' Δℓ_b and Δℓ₁.
    private sealed record Rung(double B, int M, double[] Births, double[] Deaths, double[] Changes, double[][] Steps, double[] Prunes,
                               double[] DeathChange, double[] DeathUntempered, int NullChanges);

    [Fact]
    public void StructureMoves_AlongTheDenseLadder()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman, new TemperingSpec());
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        Assert.Equal(TemperedForm.Integrated, model.Target.Form);
        int chains = fixture.Config.Sampling.Chains;
        var clock = System.Diagnostics.Stopwatch.StartNew();
        var rungs = new Rung[chains][];
        var seconds = new double[chains];
        int[] seeds = SeedTree.Derive(127, chains);
        Parallel.For(0, chains, c =>
        {
            var watch = System.Diagnostics.Stopwatch.StartNew();
            BassLadder coarse = BassMixing.Ladder(model, c, fixture.Config.Sampling.Tempering!);
            BassMixing.Burn(coarse);
            double[] inverseTemperatures = Dense(coarse.Tempering.InverseTemperatures);
            var spec = new TemperingSpec { Levels = inverseTemperatures.Length, Adaptive = false, MaxLevels = Math.Max(64, inverseTemperatures.Length) };
            BassLadder ladder = BassMixing.Ladder(model, c, spec, inverseTemperatures);
            BassMixing.Burn(ladder);
            Assert.Equal(inverseTemperatures, ladder.Tempering.InverseTemperatures);
            int[] rungSeeds = SeedTree.Derive(seeds[c], inverseTemperatures.Length);
            rungs[c] = [.. Enumerable.Range(0, inverseTemperatures.Length).Select(i => Read(model, ladder, i, new Xoshiro256PlusPlus(rungSeeds[i])))];
            seconds[c] = watch.Elapsed.TotalSeconds;
        });
        _out.WriteLine($"## P4, Friedman tempered: structure moves along M5's dense ladder, at the end of its {BassMixing.BurnIn:N0} scans of burn-in " +
                       $"({clock.Elapsed.TotalSeconds:F1} s; per chain {string.Join(", ", seconds.Select(s => s.ToString("F1")))} s)");
        for (int c = 0; c < chains; c++) Table(c, rungs[c]);
        Summary(rungs);
    }

    // One rung's state scored under its replica's tempered target; the design and counts come back as they were.
    private static Rung Read(BassModel model, BassLadder ladder, int i, Xoshiro256PlusPlus rng)
    {
        ReversibleJumpChain<BassState> replica = ladder.Tempering.Replicas[i];
        BassStructureMoves structure = ladder.Structures[i];
        BassState state = replica.Current;
        Temperature temperature = replica.Temperature;
        double b = temperature.InverseTemperature;
        Assert.Equal(ladder.Tempering.InverseTemperatures[i], b);
        Assert.Equal(TemperingPath.Likelihood, temperature.Path);
        // The rung's own density: the replica's cached log density is its temperature's on this state.
        BassEvaluation now = model.Target.Evaluate(state);
        Assert.Equal(BassMixing.Bits(temperature.LogDensity(now)), BassMixing.Bits(replica.CurrentLogDensity));
        var untouched = new BassMixing.Untouched(model, structure, () => replica.Current);
        IReadOnlyList<IRjMove<BassState>> moves = structure.Moves;

        // Births and changes through the replica's moves.
        IRjMove<BassState> birth = BassMixing.Move(moves, "birth"), death = BassMixing.Move(moves, "death"), change = BassMixing.Move(moves, "change");
        BassTermCounts counts = BassMixing.Counts(model, state);
        var births = new double[Births];
        for (int k = 0; k < Births; k++)
        {
            if (birth.Propose(state, rng) is not { } p) continue;   // an empty draw: no move
            BassTerm term = p.Candidate.Terms[^1];
            double within = -Math.Log(state.Count + 1.0) - model.Proposal.LogBirth(term, counts);
            births[k] = BassMixing.Score(model, moves, birth, state, p, model.Prior.LogPrior(term), within, state.Terms.Contains(term), temperature).Acceptance;
        }
        var changes = new List<double>(Changes);
        int nulls = 0, unavailable = 0;
        for (int k = 0; k < Changes; k++)
        {
            if (change.Propose(state, rng) is not { } p) { unavailable++; continue; }   // M = 0: no term to change
            if (ReferenceEquals(p.Candidate, state)) { nulls++; continue; }   // a redraw of the same factor
            BassTerm added = p.Candidate.Terms[^1];
            int at = 0;
            while (at < state.Count - 1 && ReferenceEquals(p.Candidate.Terms[at], state.Terms[at])) at++;
            BassTerm old = state.Terms[at];
            changes.Add(BassMixing.Score(model, moves, change, state, p, model.Prior.LogPrior(added) - model.Prior.LogPrior(old), p.LogProposalRatio,
                                         state.Terms.Contains(added), temperature).Acceptance);
        }

        // Every term's death, as the death move proposes it.
        BassDesign design = structure.DesignAt(state);
        int m = state.Count;
        var deaths = new double[m];
        var deathChange = new double[m];
        var deathUntempered = new double[m];
        for (int j = 0; j < m; j++)
        {
            BassTerm term = state.Terms[j];
            DesignProposal proposal = design.ProposeDeath(j);
            BassTerm[] rest = [.. state.Terms.Take(j), .. state.Terms.Skip(j + 1)];
            var candidate = new BassState(rest, BassStructureFit.Of(proposal, BassMixing.LogTermPrior(model, rest)), state.LogTau);
            double within = model.Proposal.LogDeathRatio(term, counts);
            BassMixing.Scored s = BassMixing.Score(model, moves, death, state, new Proposal<BassState>(candidate, within, 0.0), -model.Prior.LogPrior(term), within,
                                                   temperature: temperature);
            deaths[j] = s.Acceptance;
            deathChange[j] = s.LogMarginalChange;
            deathUntempered[j] = model.Target.Evaluate(candidate).LogMarginal - now.LogMarginal;
        }
        if (b < 0.99 && m > 0) Assert.Contains(deathChange.Zip(deathUntempered), d => Math.Abs(d.First - d.Second) > 1e-6);

        // The local knot steps of every hinge, and P2's prunes, at the rung's b.
        var steps = new List<double>[3];
        for (int d = 0; d < 3; d++) steps[d] = [];
        var prunes = new List<double>();
        for (int j = 0; j < m; j++)
        {
            BassTerm term = state.Terms[j];
            foreach (BassMixing.FactorAt f in BassMixing.Factors(term))
            {
                foreach (int delta in f.Part == BassPart.Categorical ? [] : KnotSteps)
                    if (BassMixing.Stepped(model, term, f, delta) is { } moved)
                        steps[Math.Abs(delta) - 1].Add(BassMixing.Replace(model, design, state, j, moved, 0.0, b).Exit.Acceptance);
                if (BassMixing.Order(term, f.Part) >= 2)
                    prunes.Add(BassMixing.Replace(model, design, state, j, BassMixing.Without(term, f)!, BassMixing.LogPruneRatio(model, term, f), b).Exit.Acceptance);
            }
        }
        untouched.Check();
        Assert.Equal(Changes, changes.Count + nulls + unavailable);
        Assert.Equal(state.Count == 0 ? Changes : 0, unavailable);
        return new Rung(b, m, births, deaths, [.. changes], [.. steps.Select(s => s.ToArray())], [.. prunes], deathChange, deathUntempered, nulls);
    }

    // M5's Part 2 rungs (BassMixingBandTests): the chain's own rungs above the window, the window at Δb = 0.01, then its own below it.
    private static double[] Dense(double[] frozen)
    {
        const double slack = 1e-9;
        int count = (int)Math.Round((Top - Bottom) / Step) + 1;
        return [.. frozen.Where(b => b > Top + slack), .. Enumerable.Range(0, count).Select(k => Math.Round(Top - k * Step, 2)), .. frozen.Where(b => b < Bottom - slack)];
    }

    // One chain's rungs: each kind's mean α, and the deaths' Δℓ.
    private void Table(int chain, Rung[] rungs)
    {
        _out.WriteLine("");
        _out.WriteLine($"Chain {chain}, {rungs.Length} rungs: rung | b | M | births α | deaths α: mean (largest) | median | terms over {Negligible:G0} | changes α | " +
                       "steps α at δ = ±1, ±2, ±3 | prunes α (n) | deaths' median Δℓ_b | median b·Δℓ₁ | median Δℓ₁");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int i = 0; i < rungs.Length; i++)
        {
            Rung r = rungs[i];
            _out.WriteLine($"{i} | {r.B:G4} | {r.M} | {G(r.Births.Average())} | {G(Mean(r.Deaths))} ({G(Max(r.Deaths))}) | {G(BassMixing.Quantile(r.Deaths, 0.5))} | " +
                           $"{r.Deaths.Count(a => a > Negligible)} | {G(Mean(r.Changes))} | " +
                           $"{string.Join(", ", r.Steps.Select(s => F3(Mean(s))))} | {G(Mean(r.Prunes))} ({r.Prunes.Length}) | " +
                           $"{F(BassMixing.Quantile(r.DeathChange, 0.5))} | {F(BassMixing.Quantile(r.DeathUntempered.Select(d => r.B * d), 0.5))} | {F(BassMixing.Quantile(r.DeathUntempered, 0.5))}");
        }
    }

    // Per chain: the coldest b whose deaths are not negligible, by their mean, by any one term and by the median term; and the
    // band's rungs, [0.10, 0.40], each kind's proposals pooled over them.
    private void Summary(Rung[][] rungs)
    {
        _out.WriteLine("");
        _out.WriteLine($"chain | rungs | coldest b with deaths' mean α > {Negligible:G0} | with a death's α > {Negligible:G0} | with the median term's α > {Negligible:G0} | " +
                       $"band [{Bottom}, {Top}], pooled: M range | births α | deaths α | term deaths over {Negligible:G0} | changes α | steps α at δ = ±1, ±2, ±3 | prunes α (n) | " +
                       $"rungs whose median term's death α ≤ {Negligible:G0}");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < rungs.Length; c++)
        {
            Rung[] r = rungs[c];
            Rung? mean = r.FirstOrDefault(x => Mean(x.Deaths) > Negligible), any = r.FirstOrDefault(x => Max(x.Deaths) > Negligible);
            Rung? median = r.FirstOrDefault(x => BassMixing.Quantile(x.Deaths, 0.5) > Negligible);
            Rung[] band = [.. r.Where(x => x.B >= Bottom - 1e-9 && x.B <= Top + 1e-9)];
            double[] deaths = [.. band.SelectMany(x => x.Deaths)], prunes = [.. band.SelectMany(x => x.Prunes)];
            string Coldest(Rung? x) => x is null ? "none" : x.B.ToString("G3");
            _out.WriteLine($"{c} | {r.Length} | {Coldest(mean)} | {Coldest(any)} | {Coldest(median)} | {band.Min(x => x.M)}–{band.Max(x => x.M)} | " +
                           $"{G(band.SelectMany(x => x.Births).Average())} | {G(Mean(deaths))} | {deaths.Count(a => a > Negligible)} of {deaths.Length} | " +
                           $"{G(band.SelectMany(x => x.Changes).Average())} | {string.Join(", ", Enumerable.Range(0, 3).Select(d => F3(band.SelectMany(x => x.Steps[d]).Average())))} | " +
                           $"{G(Mean(prunes))} ({prunes.Length}) | {band.Count(x => !(BassMixing.Quantile(x.Deaths, 0.5) > Negligible))} of {band.Length}");
        }
    }

    private static double Mean(double[] values) => values.Length == 0 ? double.NaN : values.Average();
    private static double Max(double[] values) => values.Length == 0 ? double.NaN : values.Max();
    private static string F(double v) => double.IsNaN(v) ? "—" : v.ToString("F2");
    private static string F3(double v) => double.IsNaN(v) ? "—" : v.ToString("F3");
    private static string G(double v) => double.IsNaN(v) ? "—" : v.ToString("G3");
}
