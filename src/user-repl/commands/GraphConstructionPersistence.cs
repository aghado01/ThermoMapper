using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;
using Graphs.Primitives;

namespace UserRepl.Commands;

internal static class GraphConstructionPersistence
{
    public const string GraphFileName = "graph.bin";
    public const string DatasetFingerprintFileName = "dataset.fingerprint";

    public static string GetGraphPath(string datasetRoot) => Path.Combine(datasetRoot, GraphFileName);
    public static string GetDatasetFingerprintPath(string datasetRoot) => Path.Combine(datasetRoot, DatasetFingerprintFileName);
    public static string GetManifestPath(string datasetRoot) => GraphConstructionManifest.PathFor(datasetRoot);

    public static string ComputeDatasetFingerprint(SpcUserDataset dataset, IDistanceMetric? metric)
    {
        if (dataset is null) throw new ArgumentNullException(nameof(dataset));

        using var hasher = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        AppendUtf8(hasher, dataset.Metadata.TryGetValue("Source", out var rawSource)
            ? rawSource?.ToString() ?? string.Empty
            : "unknown");
        AppendUtf8(hasher, dataset.Metadata.TryGetValue("Generator", out var rawGenerator)
            ? rawGenerator?.ToString() ?? string.Empty
            : string.Empty);
        AppendUtf8(hasher, metric is null
            ? "DefaultEuclidean"
            : metric.GetType().FullName ?? metric.ToString() ?? string.Empty);
        AppendUtf8(hasher, metric?.Properties.BandwidthStrategy.ToString() ?? string.Empty);

        AppendUtf8(hasher, dataset.Features.Length.ToString(CultureInfo.InvariantCulture));
        if (dataset.Features.Length > 0)
            AppendUtf8(hasher, dataset.Features[0].Length.ToString(CultureInfo.InvariantCulture));

        for (int i = 0; i < dataset.Features.Length; i++)
        {
            double[] row = dataset.Features[i];
            for (int j = 0; j < row.Length; j++)
                hasher.AppendData(BitConverter.GetBytes(row[j]));
        }

        return Convert.ToHexString(hasher.GetHashAndReset()).ToLowerInvariant();
    }

    public static string ChooseDatasetRoot(string baseDirectory, SpcUserDataset dataset, string datasetFingerprint, string configFingerprint)
    {
        if (string.IsNullOrWhiteSpace(baseDirectory))
            throw new ArgumentException("Base directory must be provided.", nameof(baseDirectory));
        if (dataset is null) throw new ArgumentNullException(nameof(dataset));
        if (string.IsNullOrWhiteSpace(datasetFingerprint))
            throw new ArgumentException("Dataset fingerprint must be provided.", nameof(datasetFingerprint));
        if (string.IsNullOrWhiteSpace(configFingerprint))
            throw new ArgumentException("Config fingerprint must be provided.", nameof(configFingerprint));

        string label = dataset.Metadata.TryGetValue("Generator", out var generator) && generator is string generatorName
            ? generatorName
            : dataset.Metadata.TryGetValue("CsvPath", out var csvPath) && csvPath is string path
                ? Path.GetFileNameWithoutExtension(path)
                : "dataset";

        string sanitized = SanitizePathSegment(label);
        string rootName = sanitized.Length > 0
            ? $"{sanitized}_{datasetFingerprint[..8]}_{configFingerprint[..8]}"
            : $"{datasetFingerprint[..8]}_{configFingerprint[..8]}";

        return Path.Combine(Path.GetFullPath(baseDirectory), rootName);
    }

    public static string ComputeConfigFingerprint(GraphCompilerConfig config)
    {
        if (config is null) throw new ArgumentNullException(nameof(config));

        using var hasher = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        AppendUtf8(hasher, config.Topology.Kind.ToString());
        AppendUtf8(hasher, (config.Topology.K ?? 0).ToString(CultureInfo.InvariantCulture));
        AppendUtf8(hasher, (config.Topology.Epsilon ?? 0.0).ToString(CultureInfo.InvariantCulture));
        AppendUtf8(hasher, config.Filter.Kind.ToString());
        AppendUtf8(hasher, config.Filter.MutualBandwidthSource?.ToString() ?? string.Empty);
        AppendUtf8(hasher, config.Repair.Kind.ToString());

        // The projection (kind + kernel + LMP + bandwidth-override) is the only
        // discriminated union in the config, and it already describes itself once
        // via [JsonPolymorphic]. Fingerprint it from that canonical serialized
        // form rather than re-listing the kernel variants in a switch here — the
        // last of the per-output-target re-descriptions of this union to go.
        AppendUtf8(hasher, JsonSerializer.Serialize(config.Projection, UserReplJsonContext.Default.IEdgeProjection));

        return Convert.ToHexString(hasher.GetHashAndReset()).ToLowerInvariant();
    }

    public static GraphConstructionManifest MaterializeManifest(
        GraphCompilerConfig config,
        GraphBuildResult buildResult,
        string datasetFingerprint)
    {
        if (config is null) throw new ArgumentNullException(nameof(config));
        if (buildResult is null) throw new ArgumentNullException(nameof(buildResult));
        if (string.IsNullOrWhiteSpace(datasetFingerprint))
            throw new ArgumentException("Dataset fingerprint must be provided.", nameof(datasetFingerprint));

        return new GraphConstructionManifest(
            SchemaVersion: GraphConstructionManifest.CurrentSchemaVersion,
            CreatedUtc:    DateTime.UtcNow,
            DatasetFingerprint: datasetFingerprint,
            Config:        config,
            Metric:        buildResult.Metric,
            Diagnostics:   buildResult.Diagnostics ?? new DiagnosticsLog());
    }

    public static bool TryLoadGraphArtifact(
        string datasetRoot,
        out CsrGraph graph,
        out GraphConstructionManifest? manifest)
    {
        graph = default;
        manifest = null;

        string graphPath = GetGraphPath(datasetRoot);
        string manifestPath = GetManifestPath(datasetRoot);
        string fingerprintPath = GetDatasetFingerprintPath(datasetRoot);

        if (!File.Exists(graphPath) || !File.Exists(manifestPath) || !File.Exists(fingerprintPath))
            return false;

        try
        {
            manifest = GraphConstructionManifest.ReadFrom(datasetRoot);
            string storedFingerprint = File.ReadAllText(fingerprintPath).Trim();
            if (!string.Equals(storedFingerprint, manifest.DatasetFingerprint, StringComparison.OrdinalIgnoreCase))
                return false;

            using var stream = File.Open(graphPath, FileMode.Open, FileAccess.Read, FileShare.Read);
            using var reader = new BinaryReader(stream);
            graph = CsrGraph.FromBinary(reader);
            return true;
        }
        catch
        {
            graph = default;
            manifest = null;
            return false;
        }
    }

    public static void WriteGraphArtifact(
        string datasetRoot,
        GraphBuildResult buildResult,
        GraphConstructionManifest manifest)
    {
        if (string.IsNullOrWhiteSpace(datasetRoot))
            throw new ArgumentException("Dataset root must be provided.", nameof(datasetRoot));
        if (buildResult is null) throw new ArgumentNullException(nameof(buildResult));
        if (manifest is null) throw new ArgumentNullException(nameof(manifest));

        Directory.CreateDirectory(datasetRoot);
        File.WriteAllText(GetDatasetFingerprintPath(datasetRoot), manifest.DatasetFingerprint);

        string graphPath = GetGraphPath(datasetRoot);
        using var stream = File.Open(graphPath, FileMode.Create, FileAccess.Write, FileShare.None);
        using var writer = new BinaryWriter(stream);
        buildResult.Graph.WriteTo(writer);

        manifest.WriteTo(datasetRoot);
    }

    private static void AppendUtf8(IncrementalHash hasher, string text)
    {
        byte[] bytes = Encoding.UTF8.GetBytes(text);
        hasher.AppendData(bytes);
    }

    private static string SanitizePathSegment(string value)
    {
        if (string.IsNullOrWhiteSpace(value))
            return "dataset";

        char[] invalid = Path.GetInvalidFileNameChars();
        var builder = new StringBuilder(value.Length);
        bool lastWasSeparator = false;

        foreach (char c in value)
        {
            bool isInvalid = Array.IndexOf(invalid, c) >= 0;
            bool isSeparator = isInvalid || char.IsWhiteSpace(c) || c == '.';
            if (isSeparator)
            {
                if (!lastWasSeparator)
                {
                    builder.Append('-');
                    lastWasSeparator = true;
                }
                continue;
            }

            builder.Append(c);
            lastWasSeparator = false;
        }

        string sanitized = builder.ToString().Trim('-');
        return string.IsNullOrWhiteSpace(sanitized) ? "dataset" : sanitized;
    }

}
