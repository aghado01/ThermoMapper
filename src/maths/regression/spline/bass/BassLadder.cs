using System;
using System.Collections.Generic;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One tempered chain (D05, D10; stage 4's S8): a <see cref="ParallelTempering{TState}"/> ladder over replicas of D04's
/// sweep on the path and form <see cref="TemperingSpec"/> declares, the b = 1 replica the sampler.
/// <list type="bullet">
/// <item><b>Replicas.</b> Each runs its own <see cref="BassStructureMoves"/> (its design, counts and tallies) and, in the g
/// mode, its own <see cref="BassTauBlock"/> after its move, at its rung: a move accepts on ℓ_b(new) − ℓ_b(old) plus the
/// untempered prior and proposal ratios, and τ's slice targets p(τ)·p_b.</item>
/// <item><b>Swaps.</b> After every sweep, DEO's adjacent pairs, each accepted on the likelihood path on ℓ_b alone, since
/// the prior cancels (D05). A state a swap hands over brings its design: the receiving replica adopts an exact copy of
/// the giver's, so a structure one replica admitted by incremental updates is neither refactored nor re-admitted by
/// another, the reference rung's replica included, whose design a checkpoint saves. A state from the reference rung,
/// whose sampler loaded it whole, is loaded whole, by the same admission.</item>
/// <item><b>The ladder.</b> Geometric from 1 down to the hot end, then b = 0 with the exact reference rung, whose replica
/// draws from <see cref="BassPriorSampler"/> every scan in place of its sweep. <see cref="Adapt"/> re-spaces the
/// interior rungs during burn-in; <see cref="EndBurnIn"/> freezes the ladder, restarts the tallies and, for an evidence,
/// starts each rung's trace (D10).</item>
/// <item><b>Streams.</b> One seed fans out, through the shared ladder's <see cref="SeedTree"/>, into each replica's move
/// stream, the swap stream, the reference rung's stream and each replica's τ stream.</item>
/// </list>
/// Step and adapt it through its own <see cref="Step"/> and <see cref="Adapt"/>: the shared ladder's, called directly, would
/// skip the hand-over and the freeze. The model is shared read-only. Not thread-safe: one per chain.
/// </summary>
public sealed class BassLadder
{
    private readonly BassStructureMoves[] _structures;

    /// <summary>
    /// A ladder of <paramref name="model"/> declared by <paramref name="spec"/>, every replica started at
    /// <paramref name="start"/> (which must lie in F) with u = <paramref name="logTau"/>, the reference rung on an exact
    /// draw; its streams derive from <paramref name="seed"/>. The model's target must be in the spec's form.
    /// </summary>
    public BassLadder(BassModel model, TemperingSpec spec, IReadOnlyList<BassTerm> start, double logTau, int seed)
        : this(model, spec, Rungs(spec), start, logTau, seed)
    {
    }

    // A ladder on the declared rungs `betas` (cold first, b = 0 last with the reference rung), for tests that place them.
    internal BassLadder(BassModel model, TemperingSpec spec, double[] betas, IReadOnlyList<BassTerm> start, double logTau, int seed)
    {
        Check(model, spec);
        ArgumentNullException.ThrowIfNull(betas);
        if (betas.Length != spec.Levels) throw new ArgumentException($"The ladder declares {spec.Levels} levels, not {betas.Length}.", nameof(betas));
        Model = model;
        Spec = spec;
        var cold = new BassStructureMoves(model, start, logTau);
        BassMovesSnapshot copy = cold.Snapshot(cold.Start);
        _structures = new BassStructureMoves[betas.Length];
        _structures[0] = cold;
        // Every replica starts on the cold one's state and a copy of its design, so the start is one structure throughout.
        IReadOnlyList<IRjMove<BassState>> ReplicaMoves(int i) => (_structures[i] ??= new BassStructureMoves(model, cold.Start, copy)).Moves;
        Reference = spec.ReferenceRung ? new BassPriorSampler(model) : null;
        Tempering = new ParallelTempering<BassState>(ReplicaMoves, model.Target, cold.Start, betas, spec.Path, seed, spec.Schedule, Reference,
                                                     Blocks(model));
    }

    private BassLadder(BassModel model, TemperingSpec spec, ParallelTempering<BassState> tempering, BassStructureMoves[] structures,
                       BassPriorSampler? reference, BassLadderMarks? burn)
    {
        Model = model;
        Spec = spec;
        Tempering = tempering;
        _structures = structures;
        Reference = reference;
        BurnMarks = burn;
    }

    /// <summary>The model.</summary>
    public BassModel Model { get; }

    /// <summary>The ladder's declaration.</summary>
    public TemperingSpec Spec { get; }

    /// <summary>The shared ladder: the replicas, their rungs, the swaps and their statistics.</summary>
    public ParallelTempering<BassState> Tempering { get; }

    /// <summary>Each replica's structure moves, cold first.</summary>
    public IReadOnlyList<BassStructureMoves> Structures => _structures;

    /// <summary>The reference rung's exact sampler of the prior on F; null without the rung.</summary>
    public BassPriorSampler? Reference { get; }

    /// <summary>The b = 1 replica, whose draws are the posterior.</summary>
    public ReversibleJumpChain<BassState> Moves => Tempering.ColdChain;

    /// <summary>The b = 1 replica's structure moves, design and tallies.</summary>
    public BassStructureMoves Structure => _structures[0];

    /// <summary>The b = 1 replica's state.</summary>
    public BassState Current => Tempering.Cold;

    /// <summary>Whether burn-in has ended: the ladder is frozen.</summary>
    public bool Frozen => BurnMarks is not null;

    // The ladder's counters when burn-in ended, from which the readout counts the sampling phase; null before.
    internal BassLadderMarks? BurnMarks { get; private set; }

    /// <summary>One sweep of every replica, then the scan's swaps, each handed-over state with its design; the b = 1
    /// state after them.</summary>
    public BassState Step()
    {
        Tempering.Step();
        HandOver();
        return Current;
    }

    /// <summary>
    /// Re-spaces the interior rungs from the swap rejections since the last call (<see cref="ParallelTempering{TState}.AdaptLadder"/>)
    /// and returns their barrier estimate Λ̂. Burn-in only: a frozen ladder refuses.
    /// </summary>
    public double Adapt()
    {
        if (Frozen) throw new InvalidOperationException("The ladder is frozen after burn-in (D10).");
        return Tempering.AdaptLadder();
    }

    /// <summary>
    /// Ends burn-in: the ladder is frozen, every replica's tallies restart, the counters are marked so the readout covers
    /// the sampling phase, and with <paramref name="traceRungs"/> each rung's trace starts, for the evidence.
    /// </summary>
    public void EndBurnIn(bool traceRungs)
    {
        if (Frozen) throw new InvalidOperationException("Burn-in has ended already.");
        foreach (BassStructureMoves s in _structures) s.ResetTallies();
        IReadOnlyList<SwapPairStat> stats = Tempering.PairStats;
        var rejection = new double[stats.Count];
        var attempts = new long[stats.Count];
        var accepts = new long[stats.Count];
        for (int i = 0; i < stats.Count; i++)
        {
            rejection[i] = Tempering.Scans > 0 ? stats[i].MeanRejection * Tempering.Scans : 0.0;
            attempts[i] = stats[i].Attempts;
            accepts[i] = stats[i].Accepts;
        }
        BurnMarks = new BassLadderMarks(Tempering.Scans, Tempering.RoundTrips, Tempering.AnnealedRestarts, rejection, attempts, accepts,
                                        Reference?.Attempts ?? 0);
        if (traceRungs) Tempering.TraceRungs();
    }

    /// <summary>The retained draw at the b = 1 state (<see cref="BassStructureMoves.Compose"/>).</summary>
    public BassDraw Compose(Xoshiro256PlusPlus rng, BassTermTable terms) => Structure.Compose(Current, rng, terms);

    /// <summary>
    /// The ladder's restart state, the ladder unchanged: the shared ladder's (each replica's state, rung, streams and
    /// counters, the swap and reference streams, the swap statistics, the adaptation's state and the traces), each
    /// replica's moves' (its design as it is, its tallies), the reference sampler's account and the burn-in's marks.
    /// </summary>
    internal BassLadderSnapshot Snapshot()
    {
        LadderSnapshot<BassState> ladder = Tempering.Snapshot();
        var structures = new BassMovesSnapshot[_structures.Length];
        for (int i = 0; i < structures.Length; i++) structures[i] = _structures[i].Snapshot(ladder.Replicas[i].Current);
        return new BassLadderSnapshot(ladder, structures, Reference?.SaveCounts(), BurnMarks);
    }

    /// <summary>A ladder of <paramref name="model"/> declared by <paramref name="spec"/> that continues bit for bit where
    /// <paramref name="snapshot"/>'s left off.</summary>
    internal static BassLadder Restore(BassModel model, TemperingSpec spec, BassLadderSnapshot snapshot)
    {
        Check(model, spec);
        ArgumentNullException.ThrowIfNull(snapshot);
        LadderSnapshot<BassState> ladder = snapshot.Tempering;
        int levels = ladder.Replicas.Count;
        if (snapshot.Structures.Length != levels) throw new ArgumentException("One moves' snapshot per replica is required.", nameof(snapshot));
        if (ladder.Path != spec.Path || ladder.Schedule != spec.Schedule || levels != spec.Levels)
            throw new ArgumentException("The snapshot is another ladder's.", nameof(snapshot));
        if ((snapshot.Reference is null) == spec.ReferenceRung)
            throw new ArgumentException("The ladder ran " + (spec.ReferenceRung ? "without" : "with") + " a reference rung.", nameof(snapshot));
        var structures = new BassStructureMoves[levels];
        for (int i = 0; i < levels; i++) structures[i] = new BassStructureMoves(model, ladder.Replicas[i].Current, snapshot.Structures[i]);
        BassPriorSampler? reference = null;
        if (snapshot.Reference is { } counts)
        {
            reference = new BassPriorSampler(model);
            reference.RestoreCounts(counts);
        }
        ParallelTempering<BassState> tempering = ParallelTempering<BassState>.Restore(i => structures[i].Moves, model.Target, ladder, reference,
                                                                                      Blocks(model));
        return new BassLadder(model, spec, tempering, structures, reference, snapshot.Burn);
    }

    // The declared ladder: geometric from 1 to the hot end, then b = 0 with the reference rung.
    private static double[] Rungs(TemperingSpec spec)
    {
        CheckSpec(spec);
        return spec.ReferenceRung
            ? ParallelTempering<BassState>.ReferenceLadder(spec.Levels, spec.HotEnd)
            : ParallelTempering<BassState>.GeometricLadder(spec.Levels, spec.HotEnd);
    }

    private static void Check(BassModel model, TemperingSpec spec)
    {
        ArgumentNullException.ThrowIfNull(model);
        CheckSpec(spec);
        if (spec.Path == TemperingPath.Likelihood && model.Target.Form != spec.Form)
            throw new ArgumentException($"The model's target is in the {model.Target.Form} form; the ladder declares the {spec.Form}.", nameof(spec));
    }

    // What the shared ladder would refuse in its own parameters' names, refused in the spec's (BassConfig.Validate's rules).
    private static void CheckSpec(TemperingSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        if (spec.Levels < (spec.ReferenceRung ? 3 : 2))
            throw new ArgumentException($"The ladder needs at least {(spec.ReferenceRung ? 3 : 2)} levels, not {spec.Levels}.", nameof(spec));
        if (!(spec.HotEnd > 0.0 && spec.HotEnd < 1.0)) throw new ArgumentException($"The hot end {spec.HotEnd} lies outside (0, 1).", nameof(spec));
        if (spec.ReferenceRung && spec.Path != TemperingPath.Likelihood)
            throw new ArgumentException("The reference rung needs the likelihood path.", nameof(spec));
    }

    // Each replica's τ block in the g mode; none under fixed g.
    private static Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>>>? Blocks(BassModel model)
        => model.Target.SamplesTau ? _ => [new BassTauBlock(model.Target)] : null;

    // Each replica a swap handed a state it does not hold adopts the giver's design for it, the reference rung's too,
    // which takes no moves but saves its design in a checkpoint. The givers' copies are all taken before any replica
    // adopts one, since a giver is usually a receiver too. A giver that hands over its pending candidate commits it, so
    // while replicas share one structure's array (their start) it can stop holding its own state; another holder of that
    // structure then hands it back in a further round. A state no replica holds came from the reference rung, and its
    // replica loads it whole when it next needs it, by the admission the reference sampler's own load passed.
    private void HandOver()
    {
        IReadOnlyList<ReversibleJumpChain<BassState>> replicas = Tempering.Replicas;
        for (int round = 0; round < replicas.Count; round++)
        {
            List<(int Receiver, BassState State, BassDesignSnapshot Design)>? handed = null;
            for (int i = 0; i < replicas.Count; i++)
            {
                BassState state = replicas[i].Current;
                if (_structures[i].Holds(state)) continue;
                for (int j = 0; j < _structures.Length; j++)
                    if (j != i && _structures[j].Holds(state))
                    {
                        (handed ??= []).Add((i, state, _structures[j].DesignFor(state)));
                        break;
                    }
            }
            if (handed is null) return;
            foreach ((int receiver, BassState state, BassDesignSnapshot design) in handed) _structures[receiver].Adopt(state, design);
        }
    }
}

/// <summary>A <see cref="BassLadder"/>'s counters when burn-in ended: scans, round trips, annealed restarts, each pair's
/// summed rejection, swap attempts and acceptances, and the reference sampler's attempts.</summary>
internal sealed record BassLadderMarks(
    long Scans, long RoundTrips, long AnnealedRestarts, double[] Rejection, long[] PairAttempts, long[] PairAccepts, long ReferenceAttempts);

/// <summary>A <see cref="BassLadder"/>'s restart state: the shared ladder's, each replica's moves', the reference
/// sampler's account (null without the rung) and the burn-in's marks (null before its end).</summary>
internal sealed record BassLadderSnapshot(
    LadderSnapshot<BassState> Tempering, BassMovesSnapshot[] Structures, BassPriorSamplerCounts? Reference, BassLadderMarks? Burn);
