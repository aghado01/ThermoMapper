using System;
using System.Threading;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// What a BASS run needs that is not the model (the contract shapes): one bound on parallelism, a cancellation token, a
/// checkpoint request and progress. Live objects that ride the call, neither serialized nor fingerprinted; none of them
/// changes a result.
/// </summary>
public sealed record BassRunOptions
{
    /// <summary>
    /// The most chains that run at once; −1, the default, leaves it to the machine. It is the run's one bound on
    /// parallelism over chains, replicas, particles and workers (the posterior-ownership rules), so tempering and SMC
    /// read it too: what the chains leave of it can go to each one's rungs (<see cref="RungWorkers"/>), as what SMC's
    /// replicates leave goes to each one's particles (<see cref="Split"/>). Results do not depend on it.
    /// </summary>
    public int MaxDegreeOfParallelism { get; init; } = -1;

    /// <summary>
    /// How many of each tempered chain's rungs step at once (<see cref="ParallelTempering{TState}.MaxDegreeOfParallelism"/>),
    /// within what the chains leave of the bound (<see cref="Split"/>): 1, the default, steps them in order; −1 takes all the
    /// chains leave; n at most n. The draws are the same bit for bit at every setting. Each scan fans its rungs out and waits
    /// for the last, so it pays only where a rung's step is long against that wait and cores are free: on a shared, loaded
    /// machine, with steps of tens of microseconds, it was slower than in order.
    /// </summary>
    public int RungWorkers { get; init; } = 1;

    /// <summary>Stops the run between sweeps with an <see cref="OperationCanceledException"/>.</summary>
    public CancellationToken Cancellation { get; init; }

    /// <summary>Whether the fit carries the run's restart state (<see cref="BassFit.Checkpoint"/>), which
    /// <see cref="Bass.Resume"/> continues.</summary>
    public bool Checkpoint { get; init; }

    /// <summary>Receives the sweeps done so far over every chain, each time they pass a multiple of
    /// <see cref="ProgressInterval"/> and once at the end. It is called from the chains' threads.</summary>
    public IProgress<BassProgress>? Progress { get; init; }

    /// <summary>The sweeps, over every chain, between progress reports.</summary>
    public int ProgressInterval { get; init; } = 1000;

    /// <summary>
    /// The bound P split between <paramref name="units"/> chains (or SMC's replicates) and the work within each: min(units, P)
    /// units at once, and ⌊P / that⌋ workers, at least 1, left for each one's rungs (or particles).
    /// Unbounded, P is the machine's processor count. A pure function of the bound and the count, never of the machine's load.
    /// </summary>
    public (int Units, int Within) Split(int units)
    {
        if (units < 1) throw new ArgumentOutOfRangeException(nameof(units), "At least one unit runs.");
        if (MaxDegreeOfParallelism is 0 or < -1)
            throw new InvalidOperationException("MaxDegreeOfParallelism is −1, no bound, or at least 1.");
        int bound = MaxDegreeOfParallelism == -1 ? Environment.ProcessorCount : MaxDegreeOfParallelism;
        int outer = Math.Min(units, bound);
        return (outer, Math.Max(1, bound / outer));
    }

    /// <summary>The rungs each of <paramref name="chains"/> tempered chains steps at once: <see cref="RungWorkers"/>, within what
    /// the chains leave of the bound (<see cref="Split"/>).</summary>
    public int RungWorkersFor(int chains)
    {
        if (RungWorkers is 0 or < -1) throw new InvalidOperationException("RungWorkers is −1, what the chains leave, or at least 1.");
        if (RungWorkers == 1) return 1;
        int within = Split(chains).Within;
        return RungWorkers == -1 ? within : Math.Min(RungWorkers, within);
    }
}

/// <summary>A run's progress: the sweeps done over every chain, and the sweeps planned, which an early R̂ stop leaves
/// unfinished.</summary>
public readonly record struct BassProgress(long Sweeps, long PlannedSweeps);
