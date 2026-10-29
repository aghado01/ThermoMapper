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
    /// emulator (the outputs' count, the groups' coverage, the shapes) are the model's, at construction.
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
    }

    private static ArgumentException Bad(string field, string reason) => new($"CalibrationConfig.{field}: {reason}");
}
