using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's M1 (the local-move measurements brief): where a sweep's cost goes. On <see cref="BassQualificationCostTests"/>'
/// surface and grid (n ∈ {200, 1,000, 5,000} × M_max ∈ {5, 20, 50} × p ∈ {2, 10}), one chain is burned 500 sweeps as the
/// run burns it, and each part of a sweep is timed at that state, one thread: the median over nine batches of a batch's
/// mean per call, after a warm-up batch, over a pool of the proposals the chain would make there. The parts' batches run
/// round-robin with the sweep's, so the machine's load falls on every part alike.
/// <list type="bullet">
///   <item><b>A birth:</b> its term draw, its within-move ratio, the new column's evaluation (the design's own steps:
///   support, centring), its cross products Zᵀz with the columns, ỹ and the residual, the bordered Cholesky (the factor's
///   bordered proposal and its quadratic form), the rest of the design's proposal, and the candidate's bookkeeping, the
///   whole birth move less those.</item>
///   <item><b>A death's and a change's proposal,</b> each the whole move; <b>the target's evaluation</b> of a candidate;
///   <b>the commit</b> of a birth, death and change, forced on admitted proposals and restored after (its explicit residual
///   included, the factor's strict rebuild amortized over its cadence); <b>τ's block</b> with its count of evaluations;
///   <b>one composition.</b></item>
///   <item><b>The sweep:</b> the chain's step and a composition, timed whole in batches of 100, each batch on a copy of the
///   chain restored at the state with its own streams, so every batch starts there; each part weighed by its frequency per
///   sweep over those batches, and the rest (the scan's bookkeeping, the state's allocation, the refused candidates'
///   evaluation) is the sweep less the parts.</item>
/// </list>
/// Tempered: on S10's tempered Friedman, one ladder burned as the run burns it, a scan against its pairs' swap ratios
/// and its hand-overs. It reads, it does not judge; it asserts that the timed proposals committed nothing
/// (<see cref="BassMixing.Untouched"/>) and that every time is finite. The test-side replicas of the design's steps are
/// compiled fully optimized, as the engine's hot paths are by then; tiered compilation timed them unoptimized at first.
/// Run it alone, in Release.
/// </summary>
public sealed class BassMixingCostTests
{
    private const int Repeats = 9, Pool = 256, BurnIn = 500, SweepBatch = 100;
    private readonly ITestOutputHelper _out;
    public BassMixingCostTests(ITestOutputHelper output) => _out = output;

    // One grid point's parts, per call in µs, and their frequencies per sweep.
    private sealed record Parts(int N, int P, int Cap, int M)
    {
        public double Draw, Ratio, Column, Cross, Cholesky, Design, BirthMove, Death, Change, Target, CommitBirth, CommitDeath,
                      CommitChange, Rebuild, Tau, TauEvaluations, Compose, Sweep;
        public double EmptyShare, SweepEndM;
        public double BirthsPerSweep, NonEmptyPerSweep, DeathsPerSweep, ChangesPerSweep, AdmittedPerSweep, BirthAccepts, DeathAccepts, ChangeAccepts;

        // At the cap the design refuses a birth before its column, and the chain never proposes one there: the birth's
        // column, Zᵀz and Cholesky stand for one below the cap, and its residuals are not read.
        public bool AtCap => M == Cap;
        public double Candidate => AtCap ? 0.0 : (BirthMove - Draw) / (1.0 - EmptyShare) - Ratio - Design;
        public double DesignOther => AtCap ? 0.0 : Design - Column - Cross - Cholesky;

        // The sweep's parts, µs per sweep, in the table's order.
        public (string Name, double Micro)[] PerSweep()
        {
            (string, double)[] parts =
            [
                ("birth: draw", BirthsPerSweep * Draw),
                ("birth: ratio", NonEmptyPerSweep * Ratio),
                ("birth: column", NonEmptyPerSweep * Column),
                ("birth: Zᵀz", NonEmptyPerSweep * Cross),
                ("birth: Cholesky", NonEmptyPerSweep * Cholesky),
                ("birth: design other", NonEmptyPerSweep * DesignOther),
                ("birth: candidate", NonEmptyPerSweep * Candidate),
                ("death proposal", DeathsPerSweep * Death),
                ("change proposal", ChangesPerSweep * Change),
                ("target", AdmittedPerSweep * Target),
                ("commit", Weigh(BirthAccepts, CommitBirth) + Weigh(DeathAccepts, CommitDeath) + Weigh(ChangeAccepts, CommitChange)),
                ("τ block", Tau),
                ("composition", Compose),
            ];
            return [.. parts, ("rest", Sweep - parts.Sum(x => x.Item2))];
        }

        // A part's cost per sweep at its frequency; a part the chain never takes here (a birth at the cap) costs nothing.
        private static double Weigh(double frequency, double micro) => frequency == 0.0 ? 0.0 : frequency * micro;
    }

    [Fact]
    public void SweepParts_OverTheCostGrid()
    {
        if (!BassQualification.Enabled) return;   // heavy
        Measure(200, 2, 5);   // warm-up: JIT and caches
        var points = new List<Parts>();
        foreach (int n in new[] { 200, 1000, 5000 })
            foreach (int p in new[] { 2, 10 })
                foreach (int cap in new[] { 5, 20, 50 })
                    points.Add(Measure(n, p, cap));

        _out.WriteLine("## M1: each part per call, µs (median of nine batch means, one thread)");
        _out.WriteLine("n | p | M_max | M | draw | ratio | column | Zᵀz | Cholesky | design other | candidate | death | change | target | " +
                       "commit b/d/c | rebuild | τ (evals) | compose | sweep");
        _out.WriteLine(string.Join(" | ", Enumerable.Repeat("---", 19)));
        foreach (Parts x in points)
            _out.WriteLine($"{x.N} | {x.P} | {x.Cap} | {x.M} | {x.Draw:F2} | {x.Ratio:F2} | {x.Column:F2} | {x.Cross:F2} | {x.Cholesky:F2} | {(x.AtCap ? "—" : x.DesignOther.ToString("F2"))} | " +
                           $"{(x.AtCap ? "—" : x.Candidate.ToString("F2"))} | {x.Death:F2} | {x.Change:F2} | {x.Target:F2} | {x.CommitBirth:F1}/{x.CommitDeath:F1}/{x.CommitChange:F1} | " +
                           $"{x.Rebuild:F1} | {x.Tau:F2} ({x.TauEvaluations:F1}) | {x.Compose:F1} | {x.Sweep:F1}");

        _out.WriteLine("");
        _out.WriteLine("## M1: each part's share of a sweep, weighed by its frequency per sweep there");
        (string Name, double Micro)[] names = points[0].PerSweep();
        _out.WriteLine($"n | p | M_max | M | sweep µs | {string.Join(" | ", names.Select(x => x.Name))} | structure move | moves/sweep b/d/c | accepts/sweep | M at the batches' ends");
        _out.WriteLine(string.Join(" | ", Enumerable.Repeat("---", names.Length + 9)));
        foreach (Parts x in points)
        {
            (string Name, double Micro)[] parts = x.PerSweep();
            double structure = parts.Where(q => q.Name.StartsWith("birth") || q.Name is "death proposal" or "change proposal" or "target" or "commit").Sum(q => q.Micro);
            _out.WriteLine($"{x.N} | {x.P} | {x.Cap} | {x.M} | {x.Sweep:F1} | {string.Join(" | ", parts.Select(q => $"{100 * q.Micro / x.Sweep:F1}%"))} | " +
                           $"{100 * structure / x.Sweep:F1}% | {x.BirthsPerSweep:F2}/{x.DeathsPerSweep:F2}/{x.ChangesPerSweep:F2} | " +
                           $"{x.BirthAccepts + x.DeathAccepts + x.ChangeAccepts:F3} | {x.SweepEndM:F1}");
        }
    }

    // One grid point: the chain burned, every part timed at its state, then the sweep on from it.
    private Parts Measure(int n, int p, int cap)
    {
        BassModel model = BassModel.Create(Config(cap), Data(n, p));
        BassChain chain = BassMixing.Chain(model, 0);
        for (int s = 0; s < BurnIn; s++) chain.Step();
        chain.Structure.ResetTallies();
        BassState state = chain.Current;
        int m = state.Count;
        var parts = new Parts(n, p, cap, m);
        var untouched = new BassMixing.Untouched(model, chain);
        BassDesignSnapshot snapshot = untouched.Before;
        BassStructureMoves structure = chain.Structure;
        BassDesign design = structure.DesignAt(state);
        BassTermCounts counts = BassMixing.Counts(model, state);
        IReadOnlyList<IRjMove<BassState>> moves = structure.Moves;
        IRjMove<BassState> birth = BassMixing.Move(moves, "birth"), death = BassMixing.Move(moves, "death"), change = BassMixing.Move(moves, "change");
        var rng = new Xoshiro256PlusPlus(seed: 113);

        // The pool: births as the chain would draw them here, each one's column, border, cross products and b.
        var terms = new List<BassTerm>();
        int draws = 0;
        while (terms.Count < Pool)
        {
            draws++;
            if (model.Proposal.DrawBirth(rng, counts) is { } t) terms.Add(t);
        }
        parts.EmptyShare = 1.0 - (double)Pool / draws;
        double[][] current = snapshot.Columns!.Cases;
        double[] residual = design.Residual.ToArray();
        double[] centred = [.. model.Data.Response.Values.ToArray().Select(v => v - design.ResponseMean)];
        // As the design: a repeated term is refused before its column, an unsupported one after it.
        double[]?[] columns = [.. terms.Select(t => state.Terms.Contains(t) ? null : Column(model, t))];
        var borders = new double[Pool][];
        var diagonals = new double[Pool];
        var bs = new double[Pool][];
        for (int k = 0; k < Pool; k++)
        {
            if (columns[k] is not { } z) continue;
            borders[k] = [.. current.Select(a => PairwiseSum.Dot(z, a))];
            diagonals[k] = PairwiseSum.Dot(z, z);
            bs[k] = [.. snapshot.Columns.CrossResponses, PairwiseSum.Dot(z, centred)];
        }

        int i = 0;
        double sink = 0.0;
        // The design's factor, or at its capacity the factor its growth would make, a strict factorization of the same Gram.
        var factor = (BorderedCholesky)design.Factor;
        if (m >= factor.Capacity)
        {
            var gram = new double[m, m];
            for (int r = 0; r < m; r++)
                for (int c = 0; c < m; c++) gram[r, c] = factor.Entry(r, c);
            factor = new BorderedCholesky(m + 1, factor.RebuildInterval);
            factor.Decompose(gram);
        }
        // The target's evaluation of admitted candidates, each a fresh evaluation as the chain makes one.
        var candidates = new List<BassState>();
        for (int k = 0; k < 20 * Pool && candidates.Count < Pool; k++)
            foreach (IRjMove<BassState> move in moves)
                if (move.Propose(state, rng) is { } q && model.Target.Refusal(q.Candidate.Fit) is null) candidates.Add(q.Candidate);

        // The commits, forced on proposals from the state, the design restored before each and after every batch. A commit
        // whose accept fell on the factor's rebuild cadence (a rebuild leaves its drift recorded, a strict accept clears
        // it) has the rebuild taken out afterwards; every commit carries the rebuild's amortized share instead.
        BassTerm[] supported = [.. terms.Where((t, k) => columns[k] is not null)];
        Func<int, double> Commit(Func<int, DesignProposal> propose, long[] fired)
        {
            int k = 0;
            Func<int, double> each = EachBatch(() =>
            {
                design.Restore(snapshot);
                if (!propose(k++).Admitted) return -1;
                long t0 = Stopwatch.GetTimestamp();
                design.Accept();
                long ticks = Stopwatch.GetTimestamp() - t0;
                var held = (BorderedCholesky)design.Factor;
                if (held.UpdatesSinceRebuild == 0 && held.LastDrift is not null) fired[0]++;
                fired[1]++;
                return ticks;
            });
            return round =>
            {
                double mean = each(round);
                design.Restore(snapshot);
                return mean;
            };
        }
        long[] firedBirth = new long[2], firedDeath = new long[2], firedChange = new long[2];

        // τ's block on a copy of the chain restored at the state, its target's evaluations counted.
        BassChainSnapshot fixedChain = chain.Snapshot();
        var counting = new Counting(model.Target);
        var tau = new BassTauBlock(counting);
        BassChain moved = BassChain.Restore(model, fixedChain);
        long updates = 0;
        var table = new BassTermTable();

        // The sweep, a step and a composition: each batch on a copy of the chain restored at the state with its own streams,
        // so every batch starts there; the moves' frequencies over the timed batches.
        int[] seeds = SeedTree.Derive(127, 2 * (Repeats + 1));
        long sweeps = 0, births = 0, empties = 0, deaths = 0, changes = 0, refused = 0, bornAccepted = 0, diedAccepted = 0, changedAccepted = 0, endM = 0;
        double Sweeps(int round)
        {
            (ulong a0, ulong a1, ulong a2, ulong a3) = new Xoshiro256PlusPlus(seeds[2 * (round + 1)]).SaveState();
            BassChainSnapshot reseeded = fixedChain with
            {
                Chain = fixedChain.Chain with { Rng0 = a0, Rng1 = a1, Rng2 = a2, Rng3 = a3 },
                BlockRng = new Xoshiro256PlusPlus(seeds[2 * (round + 1) + 1]).SaveState(),
            };
            BassChain run = BassChain.Restore(model, reseeded);
            IReadOnlyList<MoveStat> before = run.Moves.MoveStats();
            long emptyBefore = run.Structure.EmptyBirths, refusedBefore = run.Structure.Refusals.Total;
            long t0 = Stopwatch.GetTimestamp();
            for (int k = 0; k < SweepBatch; k++)
            {
                run.Step();
                run.Compose(rng, table);
            }
            double micro = (Stopwatch.GetTimestamp() - t0) * 1e6 / Stopwatch.Frequency / SweepBatch;
            if (round < 0) return micro;   // the warm-up round
            IReadOnlyList<MoveStat> after = run.Moves.MoveStats();
            sweeps += SweepBatch;
            births += Delta(after, before, "birth", false);
            deaths += Delta(after, before, "death", false);
            changes += Delta(after, before, "change", false);
            bornAccepted += Delta(after, before, "birth", true);
            diedAccepted += Delta(after, before, "death", true);
            changedAccepted += Delta(after, before, "change", true);
            empties += run.Structure.EmptyBirths - emptyBefore;
            refused += run.Structure.Refusals.Total - refusedBefore;
            endM += run.Current.Count;
            return micro;
        }

        // Every part's batches round-robin with the sweep's, a warm-up round and then nine, so the machine's load falls on
        // every part alike.
        var timers = new List<(string Name, Func<int, double> Batch)>
        {
            ("draw", Batch(() => model.Proposal.DrawBirth(rng, counts))),
            ("ratio", Batch(() => model.Proposal.LogBirthRatio(terms[i++ % Pool], counts))),
            ("column", Batch(() =>
            {
                BassTerm t = terms[i++ % Pool];
                if (!state.Terms.Contains(t)) Column(model, t);
            })),
            ("cross", Batch(() =>
            {
                if (columns[i++ % Pool] is { } z) sink += CrossProducts(z, current, centred, residual);
            })),
            ("cholesky", Batch(() =>
            {
                int k = i++ % Pool;
                if (columns[k] is not null) sink += BorderedProposal(factor, borders[k], diagonals[k], bs[k]);
            })),
            ("design", Batch(() => design.ProposeBirth(terms[i++ % Pool]))),
            ("birth move", Batch(() => birth.Propose(state, rng))),
            ("target", Batch(() => sink += Temperature.Posterior.LogDensity(model.Target.Evaluate(candidates[i++ % candidates.Count])))),
            ("compose", Batch(() => structure.Compose(state, rng, table))),
            ("tau", Batch(() => { tau.Update(moved.Moves, rng); updates++; })),
            ("sweep", Sweeps),
        };
        if (m > 0)
        {
            timers.Add(("death", Batch(() => death.Propose(state, rng))));
            timers.Add(("change", Batch(() => change.Propose(state, rng))));
            timers.Add(("rebuild", EachBatch(() =>
            {
                BorderedCholesky copy = snapshot.Factor.Copy();
                long t0 = Stopwatch.GetTimestamp();
                copy.Rebuild();
                return Stopwatch.GetTimestamp() - t0;
            })));
            timers.Add(("commit death", Commit(k => design.ProposeDeath(k % m), firedDeath)));
            if (supported.Length > 0) timers.Add(("commit change", Commit(k => design.ProposeChange(k % m, supported[k % supported.Length]), firedChange)));
        }
        // At the cap a birth commits from one term below it: a death committed first, untimed.
        if (supported.Length > 0)
            timers.Add(("commit birth", Commit(k =>
            {
                if (m == model.MaxBasis)
                {
                    design.ProposeDeath(k % m);
                    design.Accept();
                }
                return design.ProposeBirth(supported[k % supported.Length]);
            }, firedBirth)));
        Dictionary<string, double> median = Interleave(timers);
        untouched.Check();

        double Read(string name) => median.TryGetValue(name, out double v) ? v : double.NaN;
        parts.Draw = Read("draw");
        parts.Ratio = Read("ratio");
        parts.Column = Read("column");
        parts.Cross = Read("cross");
        parts.Cholesky = Read("cholesky");
        parts.Design = Read("design");
        parts.BirthMove = Read("birth move");
        parts.Death = m > 0 ? Read("death") : 0.0;
        parts.Change = m > 0 ? Read("change") : 0.0;
        parts.Target = Read("target");
        parts.Compose = Read("compose");
        parts.Tau = Read("tau");
        parts.TauEvaluations = (double)counting.Evaluations / updates;
        parts.Rebuild = m > 0 ? Read("rebuild") : 0.0;
        double amortized = snapshot.Factor.RebuildInterval > 0 ? parts.Rebuild / snapshot.Factor.RebuildInterval : 0.0;
        double Committed(string name, long[] fired) => Read(name) - (fired[1] > 0 ? (double)fired[0] / fired[1] : 0.0) * parts.Rebuild + amortized;
        parts.CommitBirth = Committed("commit birth", firedBirth);
        parts.CommitDeath = Committed("commit death", firedDeath);
        parts.CommitChange = Committed("commit change", firedChange);
        parts.Sweep = Read("sweep");
        parts.SweepEndM = (double)endM / Repeats;
        parts.NonEmptyPerSweep = (double)births / sweeps;
        parts.BirthsPerSweep = (double)(births + empties) / sweeps;
        parts.DeathsPerSweep = (double)deaths / sweeps;
        parts.ChangesPerSweep = (double)changes / sweeps;
        parts.AdmittedPerSweep = (double)(births + deaths + changes - refused) / sweeps;
        parts.BirthAccepts = (double)bornAccepted / sweeps;
        parts.DeathAccepts = (double)diedAccepted / sweeps;
        parts.ChangeAccepts = (double)changedAccepted / sweeps;
        Assert.Same(state, chain.Current);
        Assert.True(double.IsFinite(sink));
        foreach ((string name, double micro) in parts.PerSweep()) Assert.True(double.IsFinite(micro), $"{name} is {micro}");
        return parts;
    }

    // A move's attempts or acceptances between two readings of the chain's statistics.
    private static long Delta(IReadOnlyList<MoveStat> after, IReadOnlyList<MoveStat> before, string key, bool accepted)
    {
        MoveStat a = after.Single(s => s.Key == key), b = before.Single(s => s.Key == key);
        return accepted ? a.Accepted - b.Accepted : a.Attempts - b.Attempts;
    }

    // The design's cross products for a birth's column z: Zᵀz with every current column, zᵀz, zᵀỹ and zᵀr. The test-side
    // replicas are compiled fully optimized, as the engine's hot paths are by the time they are timed.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double CrossProducts(double[] z, double[][] current, double[] centred, double[] residual)
    {
        double sum = 0.0;
        foreach (double[] a in current) sum += PairwiseSum.Dot(z, a);
        return sum + PairwiseSum.Dot(z, z) + PairwiseSum.Dot(z, centred) + PairwiseSum.Dot(z, residual);
    }

    // The factor's bordered proposal and, from a positive definite one, its quadratic form, as the design takes them.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double BorderedProposal(BorderedCholesky factor, double[] border, double diagonal, double[] b)
        => factor.ProposeBirth(border, diagonal).PositiveDefinite ? factor.ProposedInverseQuadraticForm(b) : 0.0;

    // The new column as a scalar design evaluates it: the case column, its support by count or mass, then centred; null
    // when it falls short of the support rule (BassDesign's private evaluation, step for step).
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double[]? Column(BassModel model, BassTerm term)
    {
        BassBasis basis = model.Basis;
        basis.Validate(term);
        var a = new double[basis.Inputs.Cases.Cases];
        basis.CaseColumn(term, basis.Inputs.Cases, a);
        BassAdmission rules = model.Admission;
        if (rules.Support == SupportKind.ColumnMass)
        {
            if (!(PairwiseSum.Of(a) >= rules.SupportMinimum)) return null;
        }
        else
        {
            int nonzero = 0;
            foreach (double v in a)
                if (v != 0.0) nonzero++;
            if (nonzero < rules.SupportMinimum) return null;
        }
        double mean = PairwiseSum.Of(a) / a.Length;
        for (int k = 0; k < a.Length; k++) a[k] -= mean;
        return a;
    }

    [Fact]
    public void TemperedScan_SwapsAndHandOvers()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman, new TemperingSpec());
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        BassLadder ladder = BassMixing.Ladder(model, 0, fixture.Config.Sampling.Tempering!);
        BassMixing.Burn(ladder);
        ParallelTempering<BassState> tempering = ladder.Tempering;
        int levels = tempering.Replicas.Count, pairs = levels - 1;

        // The untempered sweep (moves and τ, no composition) on the same surface, for scale.
        BassMixing.Fixture plain = BassMixing.Get(BassMixing.Friedman);
        BassChain chain = BassMixing.Chain(BassModel.Create(plain.Config, plain.Data), 0);
        for (int s = 0; s < BassMixing.BurnIn; s++) chain.Step();

        // A scan, and the swaps it accepts; a pair's swap ratio, as the ladder records it every scan for every pair; a
        // hand-over, the giver's design copied for its replica's state and adopted by a scratch copy of the giver's moves;
        // each round-robin with the others.
        long swaps = 0, scans = 0;
        double sink = 0.0;
        int i = 0;
        int[] givers = [0, levels / 2, levels - 1];
        var timers = new List<(string Name, Func<int, double> Batch)>
        {
            ("scan", round =>
            {
                long before = tempering.SwapAccepts;
                double micro = Batch(() => ladder.Step(), batch: 20)(round);
                if (round >= 0)
                {
                    swaps += tempering.SwapAccepts - before;
                    scans += 20;
                }
                return micro;
            }),
            ("ratio", Batch(() =>
            {
                var a = tempering.Replicas[i % pairs];
                var b = tempering.Replicas[i++ % pairs + 1];
                sink += a.Temperature.LogDensity(b.CurrentEvaluation) + b.Temperature.LogDensity(a.CurrentEvaluation) - a.CurrentLogDensity - b.CurrentLogDensity;
            })),
            ("sweep", Batch(() => chain.Step())),
        };
        foreach (int j in givers)
            timers.Add(($"hand-over {j}", round =>
            {
                BassState state = tempering.Replicas[j].Current;
                BassStructureMoves giver = ladder.Structures[j];
                var receiver = new BassStructureMoves(model, state, giver.Snapshot(state));
                return Batch(() => receiver.Adopt(state, giver.DesignFor(state)))(round);
            }));
        Dictionary<string, double> median = Interleave(timers);
        double scan = median["scan"], ratio = median["ratio"], sweep = median["sweep"], accepted = (double)swaps / scans;
        double handOver = givers.Average(j => median[$"hand-over {j}"]);
        string[] giverRows = [.. givers.Select(j => $"rung {j} (b {tempering.InverseTemperatures[j]:G3}, M {tempering.Replicas[j].Current.Count}) {median[$"hand-over {j}"]:F1} µs")];
        double replicas = scan - pairs * ratio - 2 * accepted * handOver;

        Assert.True(double.IsFinite(sink));
        _out.WriteLine($"## M1, tempered Friedman: chain 0's self-sized ladder, {levels} rungs, {pairs} pairs");
        _out.WriteLine("");
        _out.WriteLine("part | µs | per scan, µs | share of a scan");
        _out.WriteLine("--- | --- | --- | ---");
        _out.WriteLine($"a scan | {scan:F0} | {scan:F0} | 100%");
        _out.WriteLine($"a pair's swap ratio (every pair, every scan) | {ratio:F2} | {pairs * ratio:F1} | {100 * pairs * ratio / scan:F2}%");
        _out.WriteLine($"a hand-over (two per accepted swap at most; {accepted:F2} accepted per scan) | {handOver:F1} | {2 * accepted * handOver:F1} | {100 * 2 * accepted * handOver / scan:F2}%");
        _out.WriteLine($"the replicas' sweeps (the scan less the above) | {replicas / levels:F1} per replica | {replicas:F0} | {100 * replicas / scan:F1}%");
        _out.WriteLine($"an untempered sweep at b = 1 on the same surface (moves and τ) | {sweep:F1} | — | —");
        _out.WriteLine($"hand-over by giver: {string.Join("; ", giverRows)}");
    }

    // A batch of `batch` calls, timed whole: its mean per call in µs.
    private static Func<int, double> Batch(Action call, int batch = 200)
        => [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
        {
            long t0 = Stopwatch.GetTimestamp();
            for (int k = 0; k < batch; k++) call();
            return (Stopwatch.GetTimestamp() - t0) * 1e6 / Stopwatch.Frequency / batch;
        };

    // A batch of calls that each time their own segment in ticks, so their setup is left out, or skip themselves with −1:
    // the mean per timed call in µs, NaN when none was timed.
    private static Func<int, double> EachBatch(Func<long> call, int batch = 100)
        => [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
        {
            long ticks = 0, counted = 0;
            for (int k = 0; k < batch; k++)
                if (call() is var t && t >= 0)
                {
                    ticks += t;
                    counted++;
                }
            return counted > 0 ? ticks * 1e6 / Stopwatch.Frequency / counted : double.NaN;
        };

    // Each timer's batch run round-robin with the others' (a warm-up round, numbered −1, then nine), so load on the machine
    // falls on every part alike: each timer's median over its nine batch means, NaN batches left out.
    private static Dictionary<string, double> Interleave(IReadOnlyList<(string Name, Func<int, double> Batch)> timers)
    {
        var means = timers.ToDictionary(t => t.Name, _ => new double[Repeats]);
        for (int r = -1; r < Repeats; r++)
            foreach ((string name, Func<int, double> batch) in timers)
            {
                double mean = batch(r);
                if (r >= 0) means[name][r] = mean;
            }
        return means.ToDictionary(kv => kv.Key, kv => BassMixing.Quantile(kv.Value.Where(double.IsFinite), 0.5));
    }

    // τ's block's target, its evaluations counted.
    private sealed class Counting(BassTarget target) : ISamplingTarget<BassState>
    {
        public long Evaluations { get; private set; }

        public ITargetEvaluation Evaluate(BassState state)
        {
            Evaluations++;
            return target.Evaluate(state);
        }
    }

    // BassQualificationCostTests' config and surface.
    private static BassConfig Config(int cap) => new()
    {
        Structure = new BassStructureSpec { MaxBasis = cap },
        Sampling = new BassSampling { Chains = 2, MasterSeed = 71, BurnIn = BurnIn, Samples = BurnIn },
    };

    private static BassData Data(int n, int p)
    {
        var rng = new Xoshiro256PlusPlus(seed: 73);
        var x = new double[p][];
        for (int v = 0; v < p; v++) x[v] = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int v = 0; v < p; v++) x[v][i] = rng.NextDouble();
            y[i] = 10 * Math.Sin(Math.PI * x[0][i] * x[1][i]) + 10 * x[1][i] * x[1][i] + 3 * Math.Sin(6 * Math.PI * x[0][i])
                   + 0.05 * BassQualification.Gaussian(rng);
        }
        return new BassData([.. Enumerable.Range(0, p).Select(v => (BassPredictor)new ContinuousPredictor($"x{v + 1}", x[v]))], new ScalarResponse(y));
    }
}
