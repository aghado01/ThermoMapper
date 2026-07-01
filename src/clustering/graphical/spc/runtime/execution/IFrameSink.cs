using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Execution;

/// <summary>
/// Sink for completed SPC task results. Decouples the executor's task
/// orchestration from the question of "what should happen with the
/// result once it's computed" — write to disk, accumulate in memory,
/// push to a message bus, or no-op.
/// </summary>
/// <remarks>
/// <para><b>Composition.</b> The executor's parallel dispatch loop
/// consults the sink for two things: (a) whether the task has already
/// been completed (resumable runs), and (b) what to do with each fresh
/// result. The sink doesn't own the executor's in-memory result array
/// for batch returns — that's tracked independently.</para>
///
/// <para><b>Threading.</b> Implementations must be thread-safe. The
/// executor calls <see cref="Accept"/> from multiple parallel workers
/// concurrently; sinks that maintain state across calls need their own
/// synchronization (concurrent collections, locks, etc.).</para>
/// </remarks>
public interface IFrameSink
{
    /// <summary>
    /// Allow the sink to short-circuit a task whose result has already
    /// been captured (e.g. a checkpoint file is already present on
    /// disk). Returns <see langword="false"/> for sinks that don't
    /// model "already completed" — those always re-execute.
    /// </summary>
    bool TaskAlreadyCompleted(SpcTaskSpec task);

    /// <summary>
    /// Called with a completed task's result. The sink decides what to
    /// do with it. Thread-safe (called concurrently from the executor's
    /// parallel loop).
    /// </summary>
    void Accept(SpcTaskSpec task, SpcRunResult result);

    /// <summary>
    /// Hydrate a previously-persisted result for <paramref name="task"/>.
    /// Returns <see langword="null"/> when the sink does not persist or
    /// cannot reconstruct the result (e.g. a tier-2 sidecar is missing).
    /// Used by <see cref="SpcExecutor.RunBatch"/> to populate
    /// task-aligned results for tasks that
    /// <see cref="TaskAlreadyCompleted"/> reported as cached.
    /// </summary>
    /// <remarks>
    /// Default implementation returns <see langword="null"/> — non-persistent
    /// sinks (<see cref="Sinks.NullFrameSink"/>,
    /// <see cref="Sinks.InMemoryFrameSink"/>) get the right behavior for
    /// free; persistent sinks override.
    /// </remarks>
    SpcRunResult? TryLoad(SpcTaskSpec task, CsrGraph graph) => null;
}
