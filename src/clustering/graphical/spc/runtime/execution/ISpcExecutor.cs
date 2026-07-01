using System.Collections.Generic;
using System.Threading;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Execution;

public interface ISpcExecutor
{
    /// <summary>
    /// Builds a flat schedule of independent SPC tasks from the requested
    /// temperatures, replicas, cycles, and output directory.
    /// </summary>
    List<SpcTaskSpec> BuildTaskList(
        IReadOnlyList<double> temperatures,
        int numReplicas,
        int q,
        AccumulationSpec accumulation,
        RunBudget budget,
        string checkpointDirectory,
        int? baseSeed = null);

    /// <summary>
    /// Execute the provided task list. Persistence is controlled by
    /// <paramref name="executionOptions"/>.
    /// </summary>
    void Run(
        CsrGraph graph,
        IReadOnlyList<SpcTaskSpec> tasks,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null);

    /// <summary>
    /// Execute the provided task list and return a batch result envelope.
    /// Persistence is controlled by <paramref name="executionOptions"/>.
    /// </summary>
    SpcBatchResult RunBatch(
        CsrGraph graph,
        IReadOnlyList<SpcTaskSpec> tasks,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null);

    /// <summary>
    /// Execute the provided task list in memory and return the results.
    /// This is a pure execution path; it does not persist checkpoint
    /// artifacts or sidecar files.
    /// </summary>
    IReadOnlyList<SpcRunResult> ExecuteAll(
        CsrGraph graph,
        IReadOnlyList<SpcTaskSpec> tasks,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null);

    /// <summary>
    /// Execute the provided Potts run specs in memory and return the results.
    /// This is a pure execution path; it does not persist checkpoint artifacts
    /// or sidecar files.
    /// </summary>
    IReadOnlyList<SpcRunResult> ExecuteAll(
        IReadOnlyList<SwRunSpec> specs,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null);

    /// <summary>
    /// Execute the provided Potts run specs and return a batch result envelope.
    /// This is a pure execution path; it does not persist checkpoint artifacts
    /// unless <see cref="SpcExecutionOptions.FrameSink"/> is explicitly
    /// handled by a future implementation.
    /// </summary>
    SpcBatchResult ExecuteBatch(
        IReadOnlyList<SwRunSpec> specs,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null);

    /// <summary>
    /// Execute a single task spec and return the unified SPC result.
    /// </summary>
    SpcRunResult ExecuteTask(CsrGraph graph, SpcTaskSpec task, SpcExecutionOptions? executionOptions = null);

    /// <summary>
    /// Execute a single Potts run spec and return the unified SPC result.
    /// </summary>
    SpcRunResult ExecuteTask(SwRunSpec spec, SpcExecutionOptions? executionOptions = null);
}
