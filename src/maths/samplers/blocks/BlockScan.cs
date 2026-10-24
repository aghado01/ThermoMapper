using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Samplers.Blocks;

/// <summary>
/// One block of a <see cref="BlockScan{TState}"/>: an update of some of a state's variables that leaves the target
/// invariant, either a draw from their full conditional (Gibbs) or a Metropolis step against it. The model owns its
/// blocks: their sufficient statistics, factors and scratch stay with it, and the primitives here
/// (<see cref="GaussianConditional"/>, <see cref="ConjugateUpdates"/>) do the drawing.
/// </summary>
public interface IBlock<in TState>
{
    /// <summary>The block's name in <see cref="BlockScan{TState}.BlockStats"/>.</summary>
    string Name { get; }

    /// <summary>The scan runs the block on every <c>Every</c>-th step, counting from 1; the default runs it every
    /// step.</summary>
    int Every => 1;

    /// <summary>Update <paramref name="state"/> in place on the chain's one stream. Returns whether the update was
    /// accepted; a conditional draw always is.</summary>
    bool Update(TState state, Xoshiro256PlusPlus rng);
}

/// <summary>
/// A block scan's restart state: a copy of its state, its stream's words, its step count and each block's counters.
/// The blocks are the caller's to rebuild, in the same order. Mutable block state, such as an adaptive step, is the
/// engine's to save beside it.
/// </summary>
public sealed record BlockScanSnapshot<TState>(
    TState State, ulong Rng0, ulong Rng1, ulong Rng2, ulong Rng3, long Steps, long[] BlockAttempts, long[] BlockAccepts);

/// <summary>
/// A systematic-scan chain over a model's state: each step runs the blocks in order on one stream, each block when its
/// cadence (<see cref="IBlock{TState}.Every"/>) falls due. It is Gibbs sampling when every block draws from its
/// conditional and Metropolis-within-Gibbs when some step instead; a collapsed block, one that integrates some
/// variables out before drawing, is just another block. As an <see cref="IChain{TDraw}"/> the ensemble shell drives
/// it like any other kernel: <see cref="Accepted"/> and <see cref="Attempts"/> sum the blocks' counters, and
/// <see cref="BlockStats"/> keeps them apart. The state is mutable and owned by the scan; <c>copy</c> takes the
/// independent copies that snapshots need.
/// </summary>
public sealed class BlockScan<TState> : IChain<TState> where TState : class
{
    private readonly IReadOnlyList<IBlock<TState>> _blocks;
    private readonly Func<TState, TState> _copy;
    private readonly Xoshiro256PlusPlus _rng;
    private readonly long[] _attempts;
    private readonly long[] _accepts;

    public BlockScan(IReadOnlyList<IBlock<TState>> blocks, TState state, Xoshiro256PlusPlus rng, Func<TState, TState> copy)
    {
        ArgumentNullException.ThrowIfNull(blocks);
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        ArgumentNullException.ThrowIfNull(copy);
        if (blocks.Count == 0) throw new ArgumentException("A scan needs at least one block.", nameof(blocks));
        foreach (IBlock<TState> block in blocks)
        {
            ArgumentNullException.ThrowIfNull(block, nameof(blocks));
            if (block.Every < 1) throw new ArgumentException($"Block '{block.Name}' has a cadence below 1.", nameof(blocks));
        }
        _blocks = blocks;
        State = state;
        _rng = rng;
        _copy = copy;
        _attempts = new long[blocks.Count];
        _accepts = new long[blocks.Count];
    }

    /// <summary>The current state, updated in place by each step.</summary>
    public TState State { get; }

    /// <summary>Steps taken so far.</summary>
    public long Steps { get; private set; }

    /// <inheritdoc/>
    public long Accepted { get; private set; }

    /// <inheritdoc/>
    public long Attempts { get; private set; }

    /// <summary>Run each block whose cadence falls due on this step, in order, and return the state.</summary>
    public TState Step()
    {
        long step = Steps + 1;
        for (int i = 0; i < _blocks.Count; i++)
        {
            IBlock<TState> block = _blocks[i];
            if (step % block.Every != 0) continue;
            _attempts[i]++;
            Attempts++;
            if (block.Update(State, _rng))
            {
                _accepts[i]++;
                Accepted++;
            }
        }
        Steps = step;
        return State;
    }

    /// <summary>Each block's cumulative attempts and acceptances, in scan order.</summary>
    public IReadOnlyList<MoveStat> BlockStats()
    {
        var stats = new MoveStat[_blocks.Count];
        for (int i = 0; i < _blocks.Count; i++) stats[i] = new MoveStat(_blocks[i].Name, _attempts[i], _accepts[i]);
        return stats;
    }

    /// <summary>This scan's restart state, with a copy of its state; the scan itself is unchanged.</summary>
    public BlockScanSnapshot<TState> Snapshot()
    {
        (ulong s0, ulong s1, ulong s2, ulong s3) = _rng.SaveState();
        return new BlockScanSnapshot<TState>(_copy(State), s0, s1, s2, s3, Steps,
                                             (long[])_attempts.Clone(), (long[])_accepts.Clone());
    }

    /// <summary>
    /// A scan that continues exactly where <paramref name="snapshot"/>'s left off, over <paramref name="blocks"/>, which
    /// must be the ones it ran with, in the same order. It starts from a copy of the snapshot's state, so the snapshot
    /// stays usable.
    /// </summary>
    public static BlockScan<TState> Restore(IReadOnlyList<IBlock<TState>> blocks, BlockScanSnapshot<TState> snapshot,
                                            Func<TState, TState> copy)
    {
        ArgumentNullException.ThrowIfNull(snapshot);
        ArgumentNullException.ThrowIfNull(copy);
        var scan = new BlockScan<TState>(blocks, copy(snapshot.State),
            Xoshiro256PlusPlus.Restore(snapshot.Rng0, snapshot.Rng1, snapshot.Rng2, snapshot.Rng3), copy);
        if (snapshot.BlockAttempts.Length != blocks.Count || snapshot.BlockAccepts.Length != blocks.Count)
            throw new ArgumentException("The snapshot was taken over a different block set.", nameof(snapshot));
        scan.Steps = snapshot.Steps;
        Array.Copy(snapshot.BlockAttempts, scan._attempts, blocks.Count);
        Array.Copy(snapshot.BlockAccepts, scan._accepts, blocks.Count);
        foreach (long a in scan._attempts) scan.Attempts += a;
        foreach (long a in scan._accepts) scan.Accepted += a;
        return scan;
    }
}
