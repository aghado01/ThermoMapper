using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Runtime.Scheduling;

namespace Clustering.Graphical.SPC.Runtime.Execution;

/// <summary>
/// Envelope for a batch of SPC runs executed by the kernel.
/// </summary>
public sealed record SpcBatchResult
{
    /// <summary>
    /// Runs that were executed as part of the batch (dense view of
    /// freshly-computed results, in completion order — not aligned
    /// to the task list). Tasks skipped because their checkpoint already
    /// existed are omitted from this list. For task-list-aligned access
    /// across both executed and resumed tasks, use
    /// <see cref="AlignedRuns"/>.
    /// </summary>
    public required IReadOnlyList<SpcRunResult> Runs { get; init; }

    /// <summary>
    /// Task-list-aligned view: <c>AlignedRuns[i]</c> is the result for
    /// <c>tasks[i]</c>, regardless of whether it was freshly executed
    /// this run or hydrated from a prior checkpoint via
    /// <see cref="IFrameSink.TryLoad"/>. <see langword="null"/> only when
    /// the sink reported the task as cached but could not reconstruct
    /// the result (non-persistent sink, missing tier-2 sidecar, etc.).
    /// Sweep strategies consume this so partial-grid + resume produces
    /// the same in-memory shape whether or not the run was a cold start.
    /// </summary>
    public required IReadOnlyList<SpcRunResult?> AlignedRuns { get; init; }

    /// <summary>
    /// Number of tasks requested for this batch.
    /// </summary>
    public required int RequestedTaskCount { get; init; }

    /// <summary>
    /// Number of tasks actually executed by the kernel.
    /// </summary>
    public required int ExecutedTaskCount { get; init; }

    /// <summary>
    /// Number of tasks skipped because their checkpoint already existed.
    /// </summary>
    public required int SkippedTaskCount { get; init; }

    /// <summary>
    /// Total wall-clock time taken by the batch execution.
    /// </summary>
    public required TimeSpan Elapsed { get; init; }

    /// <summary>
    /// The resolved worker budget the batch ran under — the requested policy
    /// plus the resolved worker count and ceiling/grading provenance.
    /// </summary>
    public required WorkerBudget WorkerBudget { get; init; }
}
