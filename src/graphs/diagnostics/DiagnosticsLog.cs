using System.Collections.Generic;
using System.Linq;
using System.Text.Json.Serialization;

using Graphs.Observables;

namespace Graphs.Diagnostics;

/// <summary>
/// Severity of a <see cref="DiagnosticMessage"/>. <c>Fatal</c> messages
/// never appear on a successful build — they escalate to
/// <c>GraphPathologyException</c> at the gate that produced
/// them. <c>Info</c> and <c>Warning</c> ride on
/// <c>GraphBuildResult.Diagnostics</c> for the CLI to render.
/// </summary>
public enum DiagnosticLevel { Info, Warning, Fatal }

/// <summary>
/// Single entry in the engine's forensic trail. Stage names match the
/// 5-stage pipeline ("Topology", "Filter", "Repair", "Refinement",
/// "Scaling") or a cross-stage tag ("Bandwidth", "DeltaCheck",
/// "AutoPick").
/// </summary>
public sealed record DiagnosticMessage(
    DiagnosticLevel Level,
    string          Stage,
    string          Text);

/// <summary>
/// Accumulated diagnostic messages from a <c>GraphCompiler.Build</c>
/// call, in emission order. The full log is persisted into
/// <c>GraphConstructionManifest</c> so a researcher can answer "what
/// did the engine actually do, and why?" without re-running anything.
/// </summary>
/// <remarks>
/// <para>The log is on by default — the engine always emits Info
/// entries for stage decisions ("auto-picked MutualKnnFilter:
/// skewness=4.2"), bandwidth resolution ("MAD bandwidth=0.347"),
/// repair stats ("added 3 MST bridge edges"). CLI consumers can
/// filter by <see cref="DiagnosticLevel"/> if they want a quieter
/// stream (e.g., <c>--quiet</c> drops Info).</para>
///
/// <para>Warning entries surface engine-internal heuristics — "median
/// weight 0.04 looks low, consider switching kernel" — that don't
/// rise to Fatal but warrant attention. CLI typically renders these
/// in yellow.</para>
///
/// <para>Fatal entries are not appended to a returned log — they
/// always throw <c>GraphPathologyException</c> at the gate
/// that detected them, so by definition a successful build's log
/// has no Fatal entries.</para>
/// </remarks>
public sealed class DiagnosticsLog
{
    [JsonInclude]
    public List<DiagnosticMessage> Messages { get; set; } = new();

    [JsonIgnore]
    public IEnumerable<DiagnosticMessage> Warnings =>
        Messages.Where(m => m.Level == DiagnosticLevel.Warning);

    [JsonIgnore]
    public IEnumerable<DiagnosticMessage> Infos =>
        Messages.Where(m => m.Level == DiagnosticLevel.Info);

    public void Info(string stage, string text)
        => Messages.Add(new DiagnosticMessage(DiagnosticLevel.Info, stage, text));

    public void Warning(string stage, string text)
        => Messages.Add(new DiagnosticMessage(DiagnosticLevel.Warning, stage, text));
}
