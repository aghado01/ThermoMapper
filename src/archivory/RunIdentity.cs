using System;
using System.Text;

namespace Archivory;

/// <summary>
/// Resolved identity for a run's root scope — the family folder name plus the
/// provenance of how it was chosen, recorded requested-vs-resolved like any
/// auto-resolved configuration value. Either an explicit request, or (when none
/// is given) the hierarchy-aware caller stub: the outermost caller's semantic
/// name. Pairs with <see cref="ArtifactScope.Root"/>, which takes the family.
/// </summary>
/// <param name="Family">The path-safe family folder name actually used.</param>
/// <param name="Source">How it was chosen: <c>explicit</c> or <c>auto:caller={stub}</c>.</param>
/// <param name="Requested">The raw requested name, or <see langword="null"/> when none was given.</param>
public sealed record RunIdentity(string Family, string Source, string? Requested)
{
    /// <summary>
    /// Resolve a run's family. An explicit <paramref name="requested"/> name wins; otherwise
    /// the <paramref name="callerStub"/> names the family. The result is sanitized to a single
    /// path-safe <c>snake_case</c> segment. There is deliberately <b>no subject inference</b> —
    /// rooting a run on a dataset is the caller's job, via an explicit request.
    /// </summary>
    public static RunIdentity Resolve(string? requested, string callerStub)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(callerStub);
        bool isExplicit = !string.IsNullOrWhiteSpace(requested);
        string family = Sanitize(isExplicit ? requested! : callerStub);
        string source = isExplicit ? "explicit" : $"auto:caller={Sanitize(callerStub)}";
        return new RunIdentity(family, source, isExplicit ? requested : null);
    }

    // Collapse whitespace and path-invalid characters to a single '_' (snake_case,
    // never doubled), so a family is always exactly one safe path segment.
    private static string Sanitize(string value)
    {
        char[] invalid = System.IO.Path.GetInvalidFileNameChars();
        var sb = new StringBuilder(value.Length);
        bool lastWasSep = false;
        foreach (char c in value)
        {
            bool isSep = char.IsWhiteSpace(c) || Array.IndexOf(invalid, c) >= 0;
            if (isSep)
            {
                if (!lastWasSep) { sb.Append('_'); lastWasSep = true; }
                continue;
            }
            sb.Append(c);
            lastWasSep = false;
        }
        string sanitized = sb.ToString().Trim('_');
        return sanitized.Length == 0 ? "run" : sanitized;
    }
}
