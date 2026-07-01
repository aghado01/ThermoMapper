using System.Text.Json.Serialization;

namespace Viz
{
    /// <summary>
    /// Describes one interactive control for a generator parameter.
    /// The viewer reads this to construct the regen UI dynamically —
    /// no hardcoded HTML per generator.
    /// </summary>
    public sealed class ParamSpec
    {
        /// <summary>Camel-case key matching the property name in generator_params JSON.</summary>
        [JsonPropertyName("key")] public string Key { get; init; } = "";
        /// <summary>Human-readable label shown in the UI.</summary>
        [JsonPropertyName("label")] public string Label { get; init; } = "";
        /// <summary>
        /// Control type: "int" | "float" | "bool" | "enum" | "vec3".
        /// int   → number input
        /// float → range slider + value span
        /// bool  → checkbox
        /// enum  → select dropdown (values from EnumValues)
        /// vec3  → 3 range sliders (labels from VecLabels, default a/b/c)
        /// </summary>
        [JsonPropertyName("type")] public string Type { get; init; } = "float";
        [JsonPropertyName("min")] public double? Min { get; init; }
        [JsonPropertyName("max")] public double? Max { get; init; }
        [JsonPropertyName("step")] public double? Step { get; init; }
        [JsonPropertyName("enum_values")] public string[]? EnumValues { get; init; }
        /// <summary>Component labels for vec3 params. Defaults to ["a","b","c"].</summary>
        [JsonPropertyName("vec_labels")] public string[]? VecLabels { get; init; }
        /// <summary>
        /// Optional display labels for float slider positions (index == slider integer position).
        /// When present, the viewer shows display_values[round(pos)] instead of the raw number.
        /// </summary>
        [JsonPropertyName("display_values")] public string[]? DisplayValues { get; init; }
        /// <summary>Optional hint text shown below the control (e.g. "0 = auto-estimate").</summary>
        [JsonPropertyName("description")] public string? Description { get; init; }
    }

    /// <summary>Named group of parameters shown as a collapsible section header.</summary>
    public sealed class ParamSection
    {
        [JsonPropertyName("label")] public string Label { get; init; } = "";
        [JsonPropertyName("params")] public ParamSpec[] Params { get; init; } = [];
    }

    /// <summary>
    /// Full schema for one generator's regen panel.
    /// Emitted as SCENE.generator_param_schema in the scene JSON.
    /// </summary>
    public sealed class GeneratorParamSchema
    {
        [JsonPropertyName("sections")] public ParamSection[] Sections { get; init; } = [];
    }
}
