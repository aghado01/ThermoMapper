using System.IO;
using System.Text.Json;
using System.Text.Json.Serialization;
using Graphs.Diagnostics;
using Graphs.Observables;

namespace UserRepl.Commands;

/// <summary>
/// Conventional JSON persistence for <see cref="GraphHealthReport"/>.
/// Written to <c>&lt;runDir&gt;/graph_health.json</c> alongside
/// <c>manifest.json</c> on every SPC run; refreshed by
/// <see cref="ExtractCommand"/> and the standalone
/// <see cref="GraphHealthCommand"/>. The file is a snapshot of the
/// reconstructible-from-config state of the graph — keeping it on
/// disk means the user can inspect health verdicts after the fact
/// without re-running the sampler.
/// </summary>
public static class GraphHealthFile
{
    public const string FileName = "graph_health.json";

    public static string PathFor(string runDirectory) => Path.Combine(runDirectory, FileName);

    public static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented        = true,
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        // Diagnostic fields can legitimately carry NaN / ±Infinity
        // (NeighborhoodScale ratios on degenerate graphs, MstBridge
        // skewness with a single bridge edge, AlgebraicConnectivity
        // fallback paths). Allow the named-literal form rather than
        // failing the whole write.
        NumberHandling       = JsonNumberHandling.AllowNamedFloatingPointLiterals,
        TypeInfoResolver     = UserReplJsonContext.Default,
    };

    public static string WriteTo(string runDirectory, GraphHealthReport report)
    {
        Directory.CreateDirectory(runDirectory);
        string path = PathFor(runDirectory);
        File.WriteAllText(path, JsonSerializer.Serialize(report, JsonOptions));
        return path;
    }
}
