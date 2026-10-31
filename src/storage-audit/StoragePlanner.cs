using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace StorageAudit;

public enum PlanActionKind
{
    /// <summary>Remove a regenerable item.</summary>
    Delete,
    /// <summary>Copy an archivable item to the archive, verify it, then remove the source.</summary>
    Archive,
}

/// <summary>One proposed action. The user reviews the plan file and sets <see cref="Enabled"/> to false to skip one.</summary>
public sealed record PlanAction
{
    public required int Id { get; init; }

    public bool Enabled { get; init; } = true;

    public required PlanActionKind Action { get; init; }

    public required string Root { get; init; }

    public required string Path { get; init; }

    public required Category Category { get; init; }

    public required string Kind { get; init; }

    public required string Reason { get; init; }

    public required long Bytes { get; init; }

    public required int Files { get; init; }

    public required DateTime NewestUtc { get; init; }

    public required string Drive { get; init; }
}

/// <summary>An item the plan leaves alone, with why: a refusal, or an <see cref="Category.Ask"/> for the user.</summary>
public sealed record HeldItem(string Root, string Path, Category Category, long Bytes, string Why);

/// <summary>The plan file: the actions <c>apply</c> may take, what it holds back, and git advice it never runs.</summary>
public sealed record StoragePlan
{
    public required DateTime CreatedUtc { get; init; }

    public required string RepoRoot { get; init; }

    public string? RulesPath { get; init; }

    public string? RulesSha256 { get; init; }

    public required double RecentWindowHours { get; init; }

    public required IReadOnlyList<PlanAction> Actions { get; init; }

    public IReadOnlyList<HeldItem> Held { get; init; } = [];

    public IReadOnlyList<string> GitAdvice { get; init; } = [];

    public void Write(string path) => File.WriteAllText(path, JsonSerializer.Serialize(this, StorageRules.JsonOptions));

    public static StoragePlan Read(string path)
        => JsonSerializer.Deserialize<StoragePlan>(File.ReadAllText(path), StorageRules.JsonOptions)
           ?? throw new InvalidDataException($"The plan file {path} is empty.");
}

/// <summary>
/// The refusals <c>plan</c> and <c>apply</c> share, so a plan never proposes what <c>apply</c> would
/// refuse. Each returns why the item must be left alone, or null.
/// </summary>
public static class Refusals
{
    public static string? Check(StorageItem item, PlanActionKind action, IReadOnlyList<string> protect, string? ownBinaryDirectory, StorageProfile profile)
    {
        Category wanted = action == PlanActionKind.Delete ? Category.Regenerable : Category.Archivable;
        if (item.Category != wanted)
            return $"the rules now call it {item.Category.ToString().ToLowerInvariant()}, not {wanted.ToString().ToLowerInvariant()}";
        if (item.Active)
            return $"written within the last {item.WindowHours:0.#} hours (newest {item.NewestUtc.ToLocalTime():yyyy-MM-dd HH:mm})";
        if (item.TrackedFiles > 0)
            return $"holds {item.TrackedFiles} tracked file(s)";
        if (item.UnignoredFiles > 0)
            return $"holds {item.UnignoredFiles} file(s) git does not ignore (notes or work in progress)";
        if (item.Links > 0)
            return $"holds {item.Links} link(s) (junction or symbolic link); remove them by hand";

        bool isRepo = string.Equals(item.Root, StorageRules.RepoRootName, StringComparison.OrdinalIgnoreCase);
        if (isRepo)
        {
            foreach (string prefix in protect.Select(StorageRules.NormalizePrefix).Append(".git/"))
            {
                string path = item.IsDirectory ? item.Path : item.Path + "/";
                if (path.StartsWith(prefix, StringComparison.OrdinalIgnoreCase)
                    || prefix.StartsWith(path, StringComparison.OrdinalIgnoreCase))
                    return $"protected: it is or holds {prefix}";
            }
        }

        if (ownBinaryDirectory is not null)
        {
            StorageRoot root = profile.Roots.First(r => string.Equals(r.Name, item.Root, StringComparison.OrdinalIgnoreCase));
            string full = System.IO.Path.GetFullPath(System.IO.Path.Combine(root.FullPath, item.Path)).TrimEnd('\\', '/');
            string own = System.IO.Path.GetFullPath(ownBinaryDirectory).TrimEnd('\\', '/');
            if (own.Equals(full, StringComparison.OrdinalIgnoreCase)
                || own.StartsWith(full + System.IO.Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase))
                return "holds the running storage-audit binary";
        }
        return null;
    }
}

/// <summary>Turns a profile into a plan: a delete per regenerable item, an archive per archivable one, unless refused.</summary>
public static class StoragePlanner
{
    public static StoragePlan Build(StorageProfile profile, StorageRules rules, string? ownBinaryDirectory)
    {
        var actions = new List<PlanAction>();
        var held = new List<HeldItem>();
        foreach (StorageItem item in profile.Items.OrderByDescending(i => i.Bytes))
        {
            if (item.Category is Category.Keep || item.Files == 0) continue;
            if (item.Category is Category.Ask)
            {
                held.Add(new HeldItem(item.Root, item.Path, item.Category, item.Bytes, $"ask: {item.Reason}"));
                continue;
            }

            PlanActionKind kind = item.Category == Category.Regenerable ? PlanActionKind.Delete : PlanActionKind.Archive;
            string? refusal = Refusals.Check(item, kind, rules.Protect, ownBinaryDirectory, profile);
            if (refusal is not null)
            {
                held.Add(new HeldItem(item.Root, item.Path, item.Category, item.Bytes, refusal));
                continue;
            }

            actions.Add(new PlanAction
            {
                Id = actions.Count + 1,
                Action = kind,
                Root = item.Root,
                Path = item.Path,
                Category = item.Category,
                Kind = item.Kind,
                Reason = item.Reason,
                Bytes = item.Bytes,
                Files = item.Files,
                NewestUtc = item.NewestUtc,
                Drive = item.Drive,
            });
        }

        return new StoragePlan
        {
            CreatedUtc = profile.CreatedUtc,
            RepoRoot = profile.RepoRoot,
            RulesPath = profile.RulesPath,
            RulesSha256 = profile.RulesSha256,
            RecentWindowHours = profile.RecentWindowHours,
            Actions = actions,
            Held = held,
            GitAdvice = GitAdvice(profile),
        };
    }

    /// <summary>Git maintenance to propose; the tool never runs it.</summary>
    public static List<string> GitAdvice(StorageProfile profile)
    {
        var advice = new List<string>();
        IReadOnlyDictionary<string, string> c = profile.Git.CountObjects;
        long Kib(string key) => c.TryGetValue(key, out string? v) && long.TryParse(v, out long n) ? n : 0;

        long loose = Kib("count"), looseKib = Kib("size"), packs = Kib("packs"), garbageKib = Kib("size-garbage");
        if (loose > 0 || packs > 1)
            advice.Add($"git gc: packs {loose} loose object(s) ({looseKib / 1024.0:0.0} MiB loose) and merges {packs} pack(s); "
                     + "frees up to roughly the loose size. Run it yourself when no session is committing; it rewrites no history.");
        if (Kib("prune-packable") > 0 || garbageKib > 0)
            advice.Add($"git prune-packed / gc: {Kib("prune-packable")} prune-packable object(s), {garbageKib / 1024.0:0.0} MiB garbage.");

        List<Worktree> prunable = profile.Git.Worktrees.Where(w => w.Prunable).ToList();
        if (prunable.Count > 0)
            advice.Add($"git worktree prune: {prunable.Count} registered worktree(s) whose folder is gone ({string.Join(", ", prunable.Select(w => w.Path))}).");
        if (profile.OrphanedWorktreeFolders.Count > 0)
            advice.Add($"Orphaned worktree folder(s), not registered with git: {string.Join(", ", profile.OrphanedWorktreeFolders)}. Check them and remove by hand.");
        if (profile.Git.ReflogEntries > 0)
            advice.Add($"Reflog: {profile.Git.ReflogEntries} HEAD entries. Expiring them frees little and drops recovery points; leave it unless asked.");
        return advice;
    }
}
