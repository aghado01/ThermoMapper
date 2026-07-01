using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;

namespace UserRepl.Commands;

/// <summary>
/// JSON preset for the <c>userrepl spc</c> subcommand. All fields are
/// nullable — present fields override the defaults; absent fields leave
/// the defaults in place. Loaded via <c>--config &lt;path&gt;</c> before
/// other CLI args are processed, so any flag passed alongside
/// <c>--config</c> overrides the corresponding preset entry.
/// </summary>
/// <remarks>
/// <para><b>Field-name convention.</b> JSON properties are snake_case
/// (matching <c>manifest.json</c>'s convention). Enum-like fields
/// (<c>topology</c>, <c>filter</c>, <c>kernel</c>, <c>schedule</c>,
/// <c>analyzer</c>, <c>partition_strategy</c>) are stored as the same
/// strings the CLI accepts and parsed through <see cref="SpcCommand"/>'s
/// shared parsers, so the validation rules don't drift between argv and
/// JSON.</para>
///
/// <para><b>Minimal preset.</b> Every field is optional, so a preset
/// can be as small as <c>{ "dataset": "BlattHierarchy" }</c> — useful
/// for "give me the standard run on this dataset" presets.</para>
/// </remarks>
public sealed record SpcPreset
{
    // ── Dataset ──────────────────────────────────────────────────────────
    public string?                              Dataset           { get; init; }
    public Dictionary<string, object?>?         Params            { get; init; }
    public string?                              DatasetFile       { get; init; }
    public string?                              LabelColumn       { get; init; }
    public bool?                                HasHeader         { get; init; }
    public string?                              Delimiter         { get; init; }
    public int?                                 Seed              { get; init; }

    // ── Output ───────────────────────────────────────────────────────────
    public string?                              BaseDir           { get; init; }
    public string?                              RunName           { get; init; }
    public bool?                                NoGuid            { get; init; }
    public string?                              CheckpointDir     { get; init; }

    // ── Graph topology ───────────────────────────────────────────────────
    public string?                              Topology         { get; init; }
    public string?                              Filter           { get; init; }
    public int?                                 K                 { get; init; }
    public double?                              Epsilon           { get; init; }
    public string?                              DistanceMetric    { get; init; }
    public bool?                                EnsureConnected   { get; init; }
    public bool?                                Lmp               { get; init; }

    // ── Coupling kernel ──────────────────────────────────────────────────
    public string?                              Kernel            { get; init; }
    public double?                              Bandwidth         { get; init; }
    public string?                              Mixture           { get; init; }
    public string?                              MixtureBandwidth  { get; init; }

    // ── Schedule ─────────────────────────────────────────────────────────
    public string?                              Schedule          { get; init; }
    public string?                              Temperatures      { get; init; }
    public int?                                 Replicas          { get; init; }
    public string?                              Sweep             { get; init; }
    public string?                              Equilibrium       { get; init; }
    public string?                              Accumulation      { get; init; }
    public int?                                 Q                 { get; init; }

    // ── Kernel bandwidth + landmark channel ──────────────────────────────
    public string?                              BandwidthStrategy { get; init; }
    public string?                              Susceptibility    { get; init; }

    // ── Profile analyzer + cut ───────────────────────────────────────────
    public string?                              Analyzer             { get; init; }
    public string?                              PartitionStrategy    { get; init; }
    public string?                              HierarchicalStrategy { get; init; }
    public double?                              Theta                { get; init; }
    public bool?                                PeripheralCapture    { get; init; }

    internal static SpcPreset Load(string path)
    {
        if (!File.Exists(path))
            throw new FileNotFoundException($"Preset file not found: {path}", path);

        string json = File.ReadAllText(path);
        var preset = JsonSerializer.Deserialize<SpcPreset>(json, JsonOptions)
            ?? throw new InvalidDataException($"Preset at {path} deserialized to null.");
        return preset;
    }

    /// <summary>
    /// Apply this preset to the given options instance. Only non-null
    /// preset fields are written; null fields leave the options
    /// untouched. Designed to run BEFORE the main argv loop so any CLI
    /// flag passed alongside <c>--config</c> overrides the preset's
    /// value.
    /// </summary>
    internal void ApplyTo(SpcCommand.Options options)
    {
        // Dataset
        if (Dataset is not null)          options.DatasetKind          = SpcCommand.ParseDatasetKind(Dataset);
        if (Params is not null)
            foreach (var kvp in Params)   options.GeneratorParameters[kvp.Key] = kvp.Value;
        if (DatasetFile is not null)      options.DatasetFile          = DatasetFile;
        if (LabelColumn is not null)      options.LabelColumn          = LabelColumn;
        if (HasHeader.HasValue)           options.HasHeader            = HasHeader.Value;
        if (Delimiter is not null)        options.Delimiter            = SpcCommand.ParseDelimiter(Delimiter);
        if (Seed.HasValue)                options.Seed                 = Seed.Value;

        // Output
        if (BaseDir is not null)          options.BaseDirectory        = BaseDir;
        if (RunName is not null)          options.RunName              = RunName;
        if (NoGuid.HasValue)              options.NoGuid               = NoGuid.Value;
        if (CheckpointDir is not null)    options.CheckpointDir        = CheckpointDir;

        // Graph topology
        if (Topology is not null)         options.TopologyKind         = SpcCommand.ParseTopologyKind(Topology);
        if (Filter is not null)           options.FilterKind           = SpcCommand.ParseFilterKind(Filter);
        if (K.HasValue)                   options.K                    = K.Value;
        if (Epsilon.HasValue)             options.Epsilon              = Epsilon.Value;
        if (DistanceMetric is not null)   options.DistanceMetricSpec   = DistanceMetric;
        if (EnsureConnected.HasValue)     options.EnsureConnected      = EnsureConnected.Value;
        if (Lmp.HasValue)                 options.ApplyLmp             = Lmp.Value;

        // Coupling kernel
        if (Kernel is not null)           options.KernelType           = SpcCommand.ParseKernelType(Kernel);
        if (Bandwidth.HasValue)           options.Bandwidth            = Bandwidth.Value;
        if (Mixture is not null)          options.MixtureSpec          = Mixture;
        if (MixtureBandwidth is not null) options.MixtureBandwidthSpec = MixtureBandwidth;

        // Schedule
        if (Schedule is not null)         options.Schedule             = SpcCommand.ParseScheduleMode(Schedule);
        if (Temperatures is not null)     options.TemperaturesSpec     = Temperatures;
        if (Replicas.HasValue)            options.Replicas             = Replicas.Value;
        if (Sweep is not null)            options.SweepBudget          = SpcCommand.ParseRunBudget(Sweep, "sweep");
        if (Equilibrium is not null)      options.EquilibriumBudget    = SpcCommand.ParseRunBudget(Equilibrium, "equilibrium");
        if (Accumulation is not null)     options.Accumulation         = SpcCommand.ParseAccumulationSpec(Accumulation);
        if (Q.HasValue)                   options.Q                    = Q.Value;

        // Kernel bandwidth + landmark channel
        if (BandwidthStrategy is not null)    options.BandwidthStrategy        = SpcCommand.ParseBandwidthStrategy(BandwidthStrategy);
        if (Susceptibility is not null)       options.SusceptibilityKind       = SpcCommand.ParseSusceptibilityKind(Susceptibility);

        // Profile analyzer + cut
        if (Analyzer is not null)             options.AnalyzerKind             = SpcCommand.ParseAnalyzerKind(Analyzer);
        if (PartitionStrategy is not null)    options.PartitionStrategyKind    = SpcCommand.ParsePartitionStrategyKind(PartitionStrategy);
        if (HierarchicalStrategy is not null) options.HierarchicalStrategyKind = SpcCommand.ParseHierarchicalStrategyKind(HierarchicalStrategy);
        if (Theta.HasValue)                   options.Theta                    = Theta.Value;
        if (PeripheralCapture.HasValue)       options.PeripheralCapture        = PeripheralCapture.Value;
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
