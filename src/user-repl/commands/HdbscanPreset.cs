using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;

namespace UserRepl.Commands;

/// <summary>
/// JSON preset for the <c>userrepl hdbscan</c> subcommand. Sister of
/// <see cref="SpcPreset"/> — same rules: all-nullable, applied before
/// CLI args parse, snake_case JSON. Use to save HDBSCAN configurations
/// (minPts, metric, etc.) for repeated comparison against SPC runs.
/// </summary>
public sealed record HdbscanPreset
{
    // ── Dataset ──────────────────────────────────────────────────────────
    public string?                              Dataset            { get; init; }
    public Dictionary<string, object?>?         Params             { get; init; }
    public string?                              DatasetFile        { get; init; }
    public string?                              LabelColumn        { get; init; }
    public bool?                                HasHeader          { get; init; }
    public string?                              Delimiter          { get; init; }
    public int?                                 Seed               { get; init; }

    // ── Output ───────────────────────────────────────────────────────────
    public string?                              BaseDir            { get; init; }
    public string?                              RunName            { get; init; }
    public bool?                                NoGuid             { get; init; }

    // ── HDBSCAN ──────────────────────────────────────────────────────────
    public int?                                 MinPts             { get; init; }
    public int?                                 MinClusterSize     { get; init; }
    public bool?                                AllowSingleCluster { get; init; }
    public string?                              DistanceMetric     { get; init; }   // e.g. "minkowski:p=1.5"

    internal static HdbscanPreset Load(string path)
    {
        if (!File.Exists(path))
            throw new FileNotFoundException($"Preset file not found: {path}", path);

        string json = File.ReadAllText(path);
        return JsonSerializer.Deserialize<HdbscanPreset>(json, JsonOptions)
            ?? throw new InvalidDataException($"Preset at {path} deserialized to null.");
    }

    internal void ApplyTo(HdbscanCommand.Options options)
    {
        // Dataset
        if (Dataset is not null)            options.DatasetKind        = HdbscanCommand.ParseDatasetKind(Dataset);
        if (Params is not null)
            foreach (var kvp in Params)     options.GeneratorParameters[kvp.Key] = kvp.Value;
        if (DatasetFile is not null)        options.DatasetFile        = DatasetFile;
        if (LabelColumn is not null)        options.LabelColumn        = LabelColumn;
        if (HasHeader.HasValue)             options.HasHeader          = HasHeader.Value;
        if (Delimiter is not null)          options.Delimiter          = HdbscanCommand.ParseDelimiter(Delimiter);
        if (Seed.HasValue)                  options.Seed               = Seed.Value;

        // Output
        if (BaseDir is not null)            options.BaseDirectory      = BaseDir;
        if (RunName is not null)            options.RunName            = RunName;
        if (NoGuid.HasValue)                options.NoGuid             = NoGuid.Value;

        // HDBSCAN
        if (MinPts.HasValue)                options.MinPts             = MinPts.Value;
        if (MinClusterSize.HasValue)        options.MinClusterSize     = MinClusterSize.Value;
        if (AllowSingleCluster.HasValue)    options.AllowSingleCluster = AllowSingleCluster.Value;
        if (DistanceMetric is not null)     options.DistanceMetricSpec = DistanceMetric;
    }

    private static readonly JsonSerializerOptions JsonOptions = new()
    {
        PropertyNamingPolicy        = JsonNamingPolicy.SnakeCaseLower,
        ReadCommentHandling         = JsonCommentHandling.Skip,
        AllowTrailingCommas         = true,
        PropertyNameCaseInsensitive = true,
        TypeInfoResolver            = UserReplJsonContext.Default,
    };
}
