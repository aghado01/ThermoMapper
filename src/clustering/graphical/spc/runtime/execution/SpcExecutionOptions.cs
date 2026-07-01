using Clustering.Graphical.SPC.Runtime.Execution.Sinks;
using Clustering.Graphical.SPC.Runtime.Scheduling;

namespace Clustering.Graphical.SPC.Runtime.Execution;

/// <summary>
/// Execution-time options for the SPC kernel.
/// </summary>
public sealed record SpcExecutionOptions
{
    /// <summary>
    /// Worker-budget policy for the executor's flat-task pool — how many
    /// workers to use and how (reserved cores, grading). Defaults to the
    /// auto policy (reserve 2 cores, grade ≥4 tasks/worker).
    /// </summary>
    public WorkerBudgetPolicy Parallelism { get; init; } = new();

    /// <summary>
    /// Sink that receives each completed task's result. Defaults to
    /// <see cref="SpcxDiskFrameSink"/> (writes canonical SPCX/SPCE
    /// files to <see cref="SpcTaskSpec.CheckpointPath"/>). Set
    /// <see cref="NullFrameSink.Instance"/> to suppress persistence;
    /// set <see cref="InMemoryFrameSink"/> to additionally accumulate
    /// completed results in memory.
    /// </summary>
    /// <remarks>
    /// Replaces the older <c>PersistArtifacts</c> boolean. The sink
    /// owns both the "skip already-completed task" decision and the
    /// per-result side-effect, so the executor's parallel loop stays
    /// uniform regardless of persistence mode.
    /// </remarks>
    public IFrameSink? FrameSink { get; init; }

    /// <summary>
    /// Optional root seed for reproducible task-level RNG derivation.
    /// When present, the executor will derive per-task seeds from this
    /// base seed unless the task already supplies an explicit seed.
    /// </summary>
    public int? BaseSeed { get; init; }

    /// <summary>
    /// Quantization factor for temperature-based seed derivation.
    /// Higher values preserve more precision when the same temperature
    /// occurs across different tasks.
    /// </summary>
    public int TemperatureQuantizationFactor { get; init; } = 1000;
}
