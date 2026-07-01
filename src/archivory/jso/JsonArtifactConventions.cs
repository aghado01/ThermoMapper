using System.Text.Json;
using System.Text.Json.Serialization;
using System.Text.Json.Serialization.Metadata;

namespace Archivory.Jso;

/// <summary>
/// The single source of truth for JSON artifact conventions across the
/// repository — naming policy, named-float handling, indentation, and the
/// type-info resolver. Replaces the per-assembly option bundles that drifted
/// (e.g. <c>UserReplJson.*Options</c>, <c>RunManifest.JsonOptions</c>).
/// </summary>
/// <remarks>
/// <para>Schema-agnostic by construction: it knows nothing about any artifact
/// type. The caller supplies its own source-generated
/// <see cref="JsonSerializerContext"/> via the <c>typeInfoResolver</c> argument,
/// so a new artifact type is registered in <i>its own</i> domain context — never
/// here. When no resolver is supplied a reflection-based
/// <see cref="DefaultJsonTypeInfoResolver"/> is used (ergonomic for ad-hoc
/// artifacts; pass a context for AOT-safe paths).</para>
///
/// <para>The flags express every profile currently in the tree:
/// the canonical artifact profile is the all-defaults call; a resolver-only
/// bundle is <c>Create(ctx, indented: false, snakeCase: false,
/// allowNamedFloatingPointLiterals: false)</c>; and so on.</para>
/// </remarks>
public static class JsonArtifactConventions
{
    public static JsonSerializerOptions Create(
        IJsonTypeInfoResolver? typeInfoResolver = null,
        bool indented = true,
        bool snakeCase = true,
        bool allowNamedFloatingPointLiterals = true)
    {
        var options = new JsonSerializerOptions
        {
            WriteIndented = indented,
            TypeInfoResolver = typeInfoResolver ?? new DefaultJsonTypeInfoResolver(),
        };

        if (snakeCase)
            options.PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower;
        if (allowNamedFloatingPointLiterals)
            options.NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals;

        return options;
    }
}
