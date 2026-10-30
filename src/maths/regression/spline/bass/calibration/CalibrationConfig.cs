using System;
using System.Collections.Generic;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// A calibration's configuration (D09): the parameters, the measurement-error model, the covariance formulation and the
/// run. Every field is a value, so the record serializes to JSON whole and is its own manifest, in
/// <see cref="BassConfig"/>'s style. The observations and the emulator ride the call.
/// </summary>
public sealed record CalibrationConfig
{
    /// <summary>The calibration parameters θ, each naming an emulator input.</summary>
    public required IReadOnlyList<CalibrationParameter> Parameters { get; init; }

    /// <summary>The measurement-error model (C-X1).</summary>
    public required MeasurementErrorModel Measurement { get; init; }

    /// <summary>The covariance formulation, declared (C1): functional or multivariate.</summary>
    public required CovarianceFormulation Covariance { get; init; }

    /// <summary>The most categorical level combinations whose joint conditional is enumerated whole (D09: 1,024, well
    /// above Francom et al.'s 162); above it, single-site updates take over.</summary>
    public int EnumerationCap { get; init; } = 1024;

    /// <summary>The conditional sampler's settings.</summary>
    public CalibrationSampling Sampling { get; init; } = new();

    /// <summary>The nested reference sampler's settings: the outer draws and each inner run's ensemble.</summary>
    public ReferenceSampling Reference { get; init; } = new();

    /// <summary>The cut SMC's settings (stage 8, <see cref="CutSmcSampler"/>); null, the default, leaves the cut posterior
    /// to the nested reference sampler (C3).</summary>
    public CutSmcSpec? Smc { get; init; }

    /// <summary>The canonical JSON options, <see cref="BassConfig.JsonOptions"/>' settings: enums by name, compact, an
    /// unknown field refused, a missing one defaulted (a required one refused), a variant's <c>Kind</c> anywhere.</summary>
    public static JsonSerializerOptions JsonOptions { get; } = CreateJsonOptions();

    private static JsonSerializerOptions CreateJsonOptions()
    {
        var options = new JsonSerializerOptions
        {
            Converters = { new JsonStringEnumConverter() },
            UnmappedMemberHandling = JsonUnmappedMemberHandling.Disallow,
            AllowOutOfOrderMetadataProperties = true,
        };
        options.MakeReadOnly(populateMissingResolver: true);
        return options;
    }

    /// <summary>This config as canonical JSON.</summary>
    public string ToJson() => JsonSerializer.Serialize(this, JsonOptions);

    /// <summary>A config read from <paramref name="json"/>, checked as <see cref="Validate"/> checks a built one.</summary>
    public static CalibrationConfig FromJson(string json)
    {
        CalibrationConfig config = JsonSerializer.Deserialize<CalibrationConfig>(json, JsonOptions)
                                   ?? throw new ArgumentException("The JSON holds no config.", nameof(json));
        config.Validate();
        return config;
    }

    /// <summary>
    /// Checks every field that needs no data, naming it in the reason; the checks against the observations and the
    /// emulator (the outputs' count, the groups' coverage, the shapes, positive definiteness) are the built pieces', at
    /// construction (<see cref="MeasurementModel"/>, <see cref="CalibrationCovariance.Create"/>).
    /// </summary>
    public void Validate()
    {
        try
        {
            _ = new ParameterLayout(Parameters ?? throw new ArgumentException("is missing"));
        }
        catch (ArgumentException e)
        {
            throw Bad("Parameters", e.Message);
        }
        try
        {
            if (Measurement is not GaussianMeasurement gaussian) throw new ArgumentException("is missing or of no known kind");
            MeasurementModel.Validate(gaussian);
        }
        catch (ArgumentException e)
        {
            throw Bad("Measurement", e.Message);
        }
        switch (Covariance)
        {
            case FunctionalFormulation f:
                if (f.Truncation is not (TruncationTerm.ResidualSet or TruncationTerm.Diagonal)) throw Bad("Covariance.Truncation", "is not a known term");
                if (f.Discrepancy is { } d)
                {
                    if (d.Basis is null || d.Basis.Length == 0) throw Bad("Covariance.Discrepancy.Basis", "is missing");
                    if (d.Scale is not { } s) throw Bad("Covariance.Discrepancy.Scale", "has no default: declare s_v (1 reproduces R)");
                    if (!(s > 0.0 && double.IsFinite(s))) throw Bad("Covariance.Discrepancy.Scale", "must be positive and finite");
                }
                break;
            case MultivariateFormulation v:
                if (v.Truncation is not (TruncationTerm.ResidualSet or TruncationTerm.Diagonal)) throw Bad("Covariance.Truncation", "is not a known term");
                if (v.Correlation is null || v.Correlation.Length == 0) throw Bad("Covariance.Correlation", "is missing");
                break;
            default:
                throw Bad("Covariance", "is missing or of no known kind");
        }
        if (EnumerationCap < 1) throw Bad("EnumerationCap", "must be at least 1");
        (Sampling ?? throw Bad("Sampling", "is missing")).Validate();
        (Reference ?? throw Bad("Reference", "is missing")).Validate();
        Smc?.Validate();
    }

    private static ArgumentException Bad(string field, string reason) => new($"CalibrationConfig.{field}: {reason}");
}
