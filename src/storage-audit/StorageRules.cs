using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace StorageAudit;

/// <summary>What the audit proposes for an item.</summary>
public enum Category
{
    /// <summary>A build output or cache that a build or run recreates; plan deletes it.</summary>
    Regenerable,
    /// <summary>Worth keeping but not here; plan moves it to the archive.</summary>
    Archivable,
    /// <summary>Stays where it is.</summary>
    Keep,
    /// <summary>Unclear; the user decides.</summary>
    Ask,
}

/// <summary>
/// One line of the rules file: a glob under a named root, the category it assigns, and why.
/// A path matched by a deeper rule is carved out of the item a shallower rule made; among rules
/// matching the same path, the first in the file wins.
/// </summary>
public sealed record StorageRule
{
    public required string Glob { get; init; }

    public required Category Category { get; init; }

    public required string Reason { get; init; }

    /// <summary>A short label for the report (build, cache, run, snapshot, chip, temp, …).</summary>
    public string Kind { get; init; } = "other";

    /// <summary><c>repo</c>, or the name of one of the rules file's temp roots.</summary>
    public string Root { get; init; } = StorageRules.RepoRootName;

    /// <summary>
    /// Among this rule's matches that share a parent folder, the newest this many (by their newest
    /// file) are kept rather than given the rule's category.
    /// </summary>
    public int? KeepNewest { get; init; }

    /// <summary>Lengthens the rules' recent window for this rule's items (a session's temp may be resumed days later); it never shortens it.</summary>
    public double? RecentWindowHours { get; init; }
}

/// <summary>
/// A folder outside the repository that the repository's sessions write temp data to. Its path
/// expands environment variables (<c>%TEMP%</c>) and <c>{repoSlug}</c>, the repository's full path
/// with every character other than a letter or digit turned into <c>-</c>.
/// </summary>
public sealed record TempRoot
{
    public required string Name { get; init; }

    public required string Path { get; init; }
}

/// <summary>A root the profiler walks: the repository, or a temp root resolved to its full path.</summary>
public sealed record StorageRoot(string Name, string FullPath, bool IsRepo);

/// <summary>
/// The editable rules file (JSON with comments): the paths and globs that map to categories, the
/// temp roots the audit may reach outside the repository, the folders it never touches, and the
/// window within which a written path counts as active.
/// </summary>
public sealed record StorageRules
{
    public const string RepoRootName = "repo";

    public static readonly JsonSerializerOptions JsonOptions = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        ReadCommentHandling = JsonCommentHandling.Skip,
        AllowTrailingCommas = true,
        WriteIndented = true,
        Converters = { new JsonStringEnumConverter(JsonNamingPolicy.CamelCase) },
    };

    /// <summary>Hours within which a written path counts as active and is left alone (at least 1).</summary>
    public double RecentWindowHours { get; init; } = 24;

    /// <summary>
    /// Repository-relative folders the audit never acts on, nor on anything holding them
    /// (<c>issues/</c>, <c>docs/</c>, …). <c>.git/</c> is always protected besides.
    /// </summary>
    public IReadOnlyList<string> Protect { get; init; } = [];

    public IReadOnlyList<TempRoot> TempRoots { get; init; } = [];

    public IReadOnlyList<StorageRule> Rules { get; init; } = [];

    public static StorageRules Load(string path) => Parse(File.ReadAllText(path));

    public static StorageRules Parse(string json)
    {
        StorageRules rules = JsonSerializer.Deserialize<StorageRules>(json, JsonOptions)
            ?? throw new InvalidDataException("The rules file is empty.");
        rules.Validate();
        return rules;
    }

    public static string Sha256Of(string path)
        => Convert.ToHexStringLower(SHA256.HashData(File.ReadAllBytes(path)));

    /// <summary>The repository and every temp root, with paths expanded; missing temp roots are kept.</summary>
    public IReadOnlyList<StorageRoot> ResolveRoots(string repoRoot)
    {
        string repo = System.IO.Path.GetFullPath(repoRoot).TrimEnd('\\', '/');
        var roots = new List<StorageRoot> { new(RepoRootName, repo, IsRepo: true) };
        foreach (TempRoot t in TempRoots)
        {
            string expanded = Environment.ExpandEnvironmentVariables(t.Path).Replace("{repoSlug}", RepoSlug(repo));
            roots.Add(new StorageRoot(t.Name, System.IO.Path.GetFullPath(expanded).TrimEnd('\\', '/'), IsRepo: false));
        }
        return roots;
    }

    /// <summary>The repository's full path with each non-alphanumeric character turned into <c>-</c>.</summary>
    public static string RepoSlug(string repoRoot)
    {
        var sb = new StringBuilder(repoRoot.Length);
        foreach (char c in repoRoot) sb.Append(char.IsAsciiLetterOrDigit(c) ? c : '-');
        return sb.ToString();
    }

    /// <summary>Normalizes a protect entry to a relative folder prefix ending in <c>/</c>.</summary>
    public static string NormalizePrefix(string folder)
        => folder.Replace('\\', '/').Trim('/') + "/";

    private void Validate()
    {
        if (RecentWindowHours < 1)
            throw new InvalidDataException("recentWindowHours must be at least 1.");

        var names = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { RepoRootName };
        foreach (TempRoot t in TempRoots)
        {
            if (string.IsNullOrWhiteSpace(t.Name) || string.IsNullOrWhiteSpace(t.Path))
                throw new InvalidDataException("Each temp root needs a name and a path.");
            if (!names.Add(t.Name))
                throw new InvalidDataException($"The temp root name '{t.Name}' is used twice or is reserved.");
        }

        foreach (string p in Protect)
            if (string.IsNullOrWhiteSpace(p) || p.Contains("..") || System.IO.Path.IsPathRooted(p))
                throw new InvalidDataException($"The protect entry '{p}' must be a folder relative to the repository.");

        foreach (StorageRule r in Rules)
        {
            if (!names.Contains(r.Root))
                throw new InvalidDataException($"The rule '{r.Glob}' names the unknown root '{r.Root}'.");
            if (r.Glob.Contains(".."))
                throw new InvalidDataException($"The rule '{r.Glob}' climbs out of its root.");
            if (string.IsNullOrWhiteSpace(r.Reason))
                throw new InvalidDataException($"The rule '{r.Glob}' needs a reason.");
            if (r.KeepNewest is < 1)
                throw new InvalidDataException($"The rule '{r.Glob}' keeps fewer than one; omit keepNewest instead.");
            if (r.RecentWindowHours is { } w && w < RecentWindowHours)
                throw new InvalidDataException($"The rule '{r.Glob}' shortens the recent window; a rule may only lengthen it.");
            try { _ = new PathGlob(r.Glob); }
            catch (ArgumentException e) { throw new InvalidDataException(e.Message); }
        }

        if (Rules.Select(r => (r.Root.ToLowerInvariant(), r.Glob.ToLowerInvariant())).Distinct().Count() != Rules.Count)
            throw new InvalidDataException("A glob appears twice under the same root.");
    }
}
