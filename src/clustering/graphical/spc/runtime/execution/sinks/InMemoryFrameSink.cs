using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;

namespace Clustering.Graphical.SPC.Runtime.Execution.Sinks;

/// <summary>
/// Frame sink that accumulates completed results into an in-memory
/// collection. Useful for callers (REPL, ad-hoc analysis harnesses)
/// that want to inspect every produced frame without going through
/// disk serialization.
/// </summary>
/// <remarks>
/// Results are collected in completion order — the executor's parallel
/// dispatch means this is not the task-list order. Callers that need
/// task-list-aligned results should use the executor's <c>RunBatch</c>
/// path (which builds an indexed array independent of the sink) and
/// pass <see cref="NullFrameSink"/> if no further persistence is needed.
/// </remarks>
public sealed class InMemoryFrameSink : IFrameSink
{
    private readonly ConcurrentBag<SpcRunResult> _results = new();

    public bool TaskAlreadyCompleted(SpcTaskSpec task) => false;

    public void Accept(SpcTaskSpec task, SpcRunResult result)
        => _results.Add(result);

    /// <summary>Snapshot of all results captured so far.</summary>
    public IReadOnlyList<SpcRunResult> Results => _results.ToArray();

    /// <summary>Number of results captured so far.</summary>
    public int Count => _results.Count;
}
