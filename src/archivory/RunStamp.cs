using System;
using System.Globalization;

namespace Archivory;

/// <summary>
/// Canonical invocation timestamp for artifact directories: UTC, ISO-8601 basic,
/// millisecond precision, dot-suffixed — e.g. <c>20260606T164703.741Z</c>.
/// Lexicographically sortable, unambiguous (<c>T</c> and trailing <c>Z</c> mark
/// UTC), and filesystem-clean. One stamp shape for every artifact writer so every
/// run tree is named the same way.
/// </summary>
/// <remarks>
/// Millisecond precision is deliberate: "rarely collides" is not a uniqueness
/// guarantee. The stamp identifies an <em>invocation</em> and lives once, as the
/// child of the family folder (<c>{base}/{family}/{stamp}/</c>); role-named
/// sub-scopes never restate it.
/// </remarks>
public static class RunStamp
{
    /// <summary>The .NET format string: ISO-8601 basic, millisecond, UTC.</summary>
    public const string Pattern = "yyyyMMdd'T'HHmmss'.'fff'Z'";

    /// <summary>Stamp for the current UTC instant.</summary>
    public static string Now() => For(DateTime.UtcNow);

    /// <summary>Format a specific instant (converted to UTC) as a run stamp.</summary>
    public static string For(DateTime instant)
        => instant.ToUniversalTime().ToString(Pattern, CultureInfo.InvariantCulture);
}
