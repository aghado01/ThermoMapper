using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// One chain's engine-supplied restart state: what the sampler needs to continue exactly, such as its current state,
/// every RNG stream's words, adaptation state, replica assignments and swap-schedule state, and the engine's own
/// readouts and traces. Opaque to the shell, which only hands it back to <see cref="IEnsembleModel{TDraw}.RestoreChain"/>.
/// Restart state does not pool; it is never merged.
/// </summary>
public interface IChainState
{
}

/// <summary>A chain run that can save its restart state (<see cref="IChainState"/>).</summary>
public interface IResumableChainRun<TDraw> : IChainRun<TDraw>
{
    /// <summary>A snapshot of this chain's restart state, independent of the live chain.</summary>
    IChainState SaveState();
}

/// <summary>What a checkpoint belongs to: the master seed, the chain count, and the model's fingerprint of its
/// configuration and data.</summary>
public sealed record RunIdentity(int MasterSeed, int Chains, string Fingerprint);

/// <summary>
/// One chain's part of an <see cref="EnsembleCheckpoint"/>, as separate components: the engine's restart state, the
/// shell's readout statistics for this chain (streaming moments of every functional), its diagnostic trace (the
/// monitored functionals' ordered draws), and its counters at the end of burn-in.
/// </summary>
public sealed record ChainCheckpoint(
    IChainState Engine, RunningMoments[] Moments, ChainTrace Trace, long BurnAccepted, long BurnAttempts);

/// <summary>
/// The restart state of an ensemble run after some sampling rounds: its identity, the chain seeds, the draws per
/// chain so far, and each chain's components. It owns readouts, traces and restart state side by side, never as one
/// merged structure, and it holds copies, so the live run can continue without changing it. A checkpoint is taken
/// after burn-in, at a round boundary.
/// </summary>
public sealed record EnsembleCheckpoint(
    RunIdentity Identity, int[] ChainSeeds, int DrawsPerChain, IReadOnlyList<ChainCheckpoint> Chains);

/// <summary>
/// A stable fingerprint of a run's configuration and data: FNV-1a over the exact bits of the values added, written
/// as hex. The same inputs give the same fingerprint on every platform and process.
/// </summary>
public sealed class RunFingerprint
{
    private ulong _hash = 14695981039346656037UL;

    public RunFingerprint Add(ulong bits)
    {
        for (int i = 0; i < 8; i++)
        {
            _hash ^= (bits >> (8 * i)) & 0xFF;
            _hash *= 1099511628211UL;
        }
        return this;
    }

    public RunFingerprint Add(double value) => Add((ulong)BitConverter.DoubleToInt64Bits(value));

    public RunFingerprint Add(int value) => Add((ulong)(long)value);

    public RunFingerprint Add(ReadOnlySpan<double> values)
    {
        Add(values.Length);
        foreach (double v in values) Add(v);
        return this;
    }

    public RunFingerprint Add(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        Add(text.Length);
        foreach (byte b in Encoding.UTF8.GetBytes(text)) Add((ulong)b);
        return this;
    }

    public override string ToString() => _hash.ToString("x16", CultureInfo.InvariantCulture);
}
