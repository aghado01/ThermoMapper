using System;
using System.IO;
using System.Text.Json;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;

namespace UserRepl.Commands;

/// <summary>
/// Persisted metadata for a constructed graph artifact. Lives alongside
/// <c>graph.bin</c> and the dataset fingerprint so repeated SPC sweeps can
/// reuse the same graph without rebuilding it.
/// </summary>
public sealed record GraphConstructionManifest(
    string              SchemaVersion,
    DateTime            CreatedUtc,
    string              DatasetFingerprint,
    GraphCompilerConfig Config,
    MetricProperties?   Metric,
    DiagnosticsLog      Diagnostics)
{
    public const string CurrentSchemaVersion = "1.1";
    public const string FileName = "graph.manifest.json";

    public static string PathFor(string graphDirectory) =>
        Path.Combine(graphDirectory, FileName);

    public void WriteTo(string graphDirectory)
    {
        if (string.IsNullOrWhiteSpace(graphDirectory))
            throw new ArgumentException("Graph directory must be provided.", nameof(graphDirectory));

        Directory.CreateDirectory(graphDirectory);
        string path = PathFor(graphDirectory);
        File.WriteAllText(path, JsonSerializer.Serialize(this, JsonOptions));
    }

    public static GraphConstructionManifest ReadFrom(string graphDirectory)
    {
        if (string.IsNullOrWhiteSpace(graphDirectory))
            throw new ArgumentException("Graph directory must be provided.", nameof(graphDirectory));

        string path = PathFor(graphDirectory);
        if (!File.Exists(path))
            throw new FileNotFoundException($"No graph manifest at {path}.", path);

        string json = File.ReadAllText(path);
        return JsonSerializer.Deserialize<GraphConstructionManifest>(json, JsonOptions)
            ?? throw new InvalidDataException($"Graph manifest at {path} deserialized to null.");
    }

    public static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = true,
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        TypeInfoResolver = UserReplJsonContext.Default,
    };
}
