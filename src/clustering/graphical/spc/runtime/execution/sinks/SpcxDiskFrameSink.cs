using System;
using System.IO;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Execution.Sinks;

/// <summary>
/// Frame sink that writes each run's <see cref="Accumulator"/> to disk as one SPCX file via
/// <see cref="AccumulatorSerializer"/> — the per-edge arrays fold inline, so the old paired SPCE
/// sidecar (and its key-matching) is retired. The default sink for the executor; preserves the
/// pre-IFrameSink persistence behavior.
/// </summary>
/// <remarks>
/// <para><b>Resumable runs.</b> <see cref="TaskAlreadyCompleted"/> returns true when the checkpoint
/// file already exists at <see cref="SpcTaskSpec.CheckpointPath"/>, so re-running the same task list
/// against the same output directory skips already-completed tasks rather than redoing them.</para>
///
/// <para><b>Atomic writes.</b> The serializer writes through Archivory's
/// <see cref="Archivory.BinarySerializerBase{T}"/> temp-rename pattern; partial writes never leave a
/// half-written file at the canonical path.</para>
/// </remarks>
public sealed class SpcxDiskFrameSink : IFrameSink
{
    public static SpcxDiskFrameSink Instance { get; } = new();

    public bool TaskAlreadyCompleted(SpcTaskSpec task)
    {
        if (task is null)
            throw new ArgumentNullException(nameof(task));
        if (task.CheckpointPath is null)
            throw new InvalidOperationException(
                $"{nameof(SpcxDiskFrameSink)} requires {nameof(SpcTaskSpec.CheckpointPath)} to be non-null.");

        return File.Exists(task.CheckpointPath);
    }

    public void Accept(SpcTaskSpec task, SpcRunResult result)
    {
        if (task is null)
            throw new ArgumentNullException(nameof(task));
        if (result is null)
            throw new ArgumentNullException(nameof(result));
        if (task.CheckpointPath is null)
            throw new InvalidOperationException(
                $"{nameof(SpcxDiskFrameSink)} requires {nameof(SpcTaskSpec.CheckpointPath)} to be non-null.");

        string? dir = Path.GetDirectoryName(task.CheckpointPath);
        if (!string.IsNullOrEmpty(dir))
            Directory.CreateDirectory(dir);

        AccumulatorSerializer.Instance.WriteToFile(result.Accumulator, task.CheckpointPath);
    }

    /// <inheritdoc />
    public SpcRunResult? TryLoad(SpcTaskSpec task, CsrGraph graph)
    {
        if (task is null)
            throw new ArgumentNullException(nameof(task));
        // CsrGraph is a struct — no null check possible or needed.
        if (task.CheckpointPath is null || !File.Exists(task.CheckpointPath))
            return null;

        var accumulator = AccumulatorSerializer.Instance.ReadFromFile(task.CheckpointPath);

        return new SpcRunResult
        {
            Graph       = graph,
            Accumulator = accumulator,
        };
    }
}
