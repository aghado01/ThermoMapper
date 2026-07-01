namespace Clustering.Graphical.SPC.Runtime.Execution.Sinks;

/// <summary>
/// No-op frame sink. Use when the caller only wants the in-memory
/// batch result and does not want any persistence or post-processing
/// side-effects. Replaces the old <c>PersistArtifacts = false</c>
/// shortcut on <see cref="SpcExecutionOptions"/>.
/// </summary>
public sealed class NullFrameSink : IFrameSink
{
    public static NullFrameSink Instance { get; } = new();

    public bool TaskAlreadyCompleted(SpcTaskSpec task) => false;

    public void Accept(SpcTaskSpec task, SpcRunResult result)
    {
        // Intentionally empty.
    }
}
