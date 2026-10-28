using System;
using System.Threading;

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
    /// read it too when they join.
    /// </summary>
    public int MaxDegreeOfParallelism { get; init; } = -1;

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
}

/// <summary>A run's progress: the sweeps done over every chain, and the sweeps planned, which an early R̂ stop leaves
/// unfinished.</summary>
public readonly record struct BassProgress(long Sweeps, long PlannedSweeps);
