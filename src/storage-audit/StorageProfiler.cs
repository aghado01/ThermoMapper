using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json.Serialization;

namespace StorageAudit;

/// <summary>
/// A path a rule classified: a folder (its path ends in <c>/</c>) or a file, under one root. It owns
/// every file and folder beneath it except those a deeper rule carved out (<see cref="Nested"/>), and
/// its size, file count and newest write cover what it owns.
/// </summary>
public sealed class StorageItem
{
    public required string Root { get; init; }

    public required string Path { get; init; }

    public required bool IsDirectory { get; init; }

    public required Category Category { get; set; }

    public required string Kind { get; init; }

    public required string Reason { get; set; }

    public required string Drive { get; init; }

    public long Bytes { get; set; }

    public int Files { get; set; }

    /// <summary>The newest write to a file or folder the item owns; a folder's time moves when entries come or go.</summary>
    public DateTime NewestUtc { get; set; }

    /// <summary>The recent window for this item: the rules', or its rule's longer one.</summary>
    public double WindowHours { get; set; }

    /// <summary>Written within <see cref="WindowHours"/>, so left alone.</summary>
    public bool Active { get; set; }

    public int TrackedFiles { get; set; }

    /// <summary>Repository files git does not ignore: notes or work in progress, never acted on.</summary>
    public int UnignoredFiles { get; set; }

    /// <summary>Junctions and symbolic links, which the walk does not follow.</summary>
    public int Links { get; set; }

    public List<string> Nested { get; } = [];

    [JsonIgnore] internal int RuleIndex { get; init; }

    [JsonIgnore] internal List<string> OwnedFiles { get; } = [];

    [JsonIgnore] internal List<string> OwnedDirectories { get; } = [];

    internal void Touch(DateTime writeUtc)
    {
        if (writeUtc > NewestUtc) NewestUtc = writeUtc;
    }
}

/// <summary>Where files no rule claims are counted.</summary>
public enum RollupBucket
{
    /// <summary>Tracked by git: source, briefs, docs. Kept.</summary>
    Tracked,
    /// <summary>Untracked and not ignored: notes and work in progress. Never acted on.</summary>
    Untracked,
    /// <summary>Ignored by git but matched by no rule: a candidate for a new rule.</summary>
    IgnoredNoRule,
    /// <summary>Under a temp root but matched by no rule.</summary>
    TempNoRule,
}

/// <summary>Files no rule claims, summed per root, top-level entry and bucket. Report-only: nothing acts on a rollup.</summary>
public sealed class Rollup
{
    public required string Root { get; init; }

    public required string Top { get; init; }

    public required RollupBucket Bucket { get; init; }

    public long Bytes { get; set; }

    public int Files { get; set; }

    public DateTime NewestUtc { get; set; }
}

/// <summary>Everything under one top-level entry of a root (<c>./</c> collects the root's own files).</summary>
public sealed class DirectoryTotal
{
    public required string Root { get; init; }

    public required string Name { get; init; }

    public long Bytes { get; set; }

    public int Files { get; set; }

    public DateTime NewestUtc { get; set; }
}

public sealed record DriveSpace(string Drive, long FreeBytes, long TotalBytes);

/// <summary>The survey: every classified item, the rollups, the totals per top-level entry, git and the drives.</summary>
public sealed class StorageProfile
{
    public required DateTime CreatedUtc { get; init; }

    public required string RepoRoot { get; init; }

    public string? RulesPath { get; init; }

    public string? RulesSha256 { get; init; }

    public required double RecentWindowHours { get; init; }

    public required IReadOnlyList<StorageRoot> Roots { get; init; }

    public required IReadOnlyList<string> MissingRoots { get; init; }

    public required IReadOnlyList<StorageItem> Items { get; init; }

    public required IReadOnlyList<Rollup> Rollups { get; init; }

    public required IReadOnlyList<DirectoryTotal> TopDirectories { get; init; }

    public required IReadOnlyList<DriveSpace> Drives { get; init; }

    public required GitStatistics Git { get; init; }

    /// <summary>Folders under <c>.claude/worktrees/</c> that git no longer lists as worktrees.</summary>
    public required IReadOnlyList<string> OrphanedWorktreeFolders { get; init; }

    public required IReadOnlyList<string> Unreadable { get; init; }

    public StorageItem? Find(string root, string path)
        => Items.FirstOrDefault(i => string.Equals(i.Root, root, StringComparison.OrdinalIgnoreCase)
                                  && string.Equals(i.Path, path, StringComparison.OrdinalIgnoreCase));
}

/// <summary>
/// Walks the repository and the rules' temp roots once, matching every folder and file against the
/// rules: a match makes an item, and the walk continues beneath it so deeper rules can carve their
/// own items out. Files no rule claims go to the <see cref="Rollup"/>s. Links are counted, not followed.
/// </summary>
public static class StorageProfiler
{
    private static readonly EnumerationOptions Options = new()
    {
        RecurseSubdirectories = false,
        IgnoreInaccessible = true,
        AttributesToSkip = 0,
        ReturnSpecialDirectories = false,
    };

    public static StorageProfile Run(string repoRoot, StorageRules rules, DateTime nowUtc, string? rulesPath = null)
    {
        IReadOnlyList<StorageRoot> roots = rules.ResolveRoots(repoRoot);
        StorageRoot repo = roots[0];
        if (!Directory.Exists(repo.FullPath))
            throw new DirectoryNotFoundException($"The repository {repo.FullPath} does not exist.");

        // A temp root walks without git's protections, so it must not reach into the repository or another root.
        foreach (StorageRoot temp in roots.Skip(1))
            foreach (StorageRoot other in roots.Where(r => r.Name != temp.Name))
                if (Overlaps(temp.FullPath, other.FullPath))
                    throw new InvalidDataException(
                        $"The temp root '{temp.Name}' ({temp.FullPath}) overlaps the root '{other.Name}' ({other.FullPath}); roots must be disjoint.");

        var walker = new Walker(rules, RepoGit.Read(repo.FullPath));
        var missing = new List<string>();
        foreach (StorageRoot root in roots)
        {
            if (Directory.Exists(root.FullPath)) walker.Walk(root, new DirectoryInfo(root.FullPath), "", null, null);
            else missing.Add(root.Name);
        }

        walker.KeepNewest();
        foreach (StorageItem item in walker.Items)
        {
            item.WindowHours = Math.Max(rules.RecentWindowHours, rules.Rules[item.RuleIndex].RecentWindowHours ?? 0);
            item.Active = item.NewestUtc > nowUtc.AddHours(-item.WindowHours);
        }

        GitStatistics git = RepoGit.ReadStatistics(repo.FullPath);
        return new StorageProfile
        {
            CreatedUtc = nowUtc,
            RepoRoot = repo.FullPath,
            RulesPath = rulesPath,
            RulesSha256 = rulesPath is not null && File.Exists(rulesPath) ? StorageRules.Sha256Of(rulesPath) : null,
            RecentWindowHours = rules.RecentWindowHours,
            Roots = roots,
            MissingRoots = missing,
            Items = walker.Items,
            Rollups = walker.Rollups.Values.OrderByDescending(r => r.Bytes).ToList(),
            TopDirectories = walker.Tops.Values.OrderByDescending(t => t.Bytes).ToList(),
            Drives = ReadDrives(roots),
            Git = git,
            OrphanedWorktreeFolders = OrphanedWorktrees(repo.FullPath, git),
            Unreadable = walker.Unreadable,
        };
    }

    private static bool Overlaps(string a, string b)
        => a.Equals(b, StringComparison.OrdinalIgnoreCase)
        || a.StartsWith(b + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase)
        || b.StartsWith(a + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase);

    private static List<DriveSpace> ReadDrives(IEnumerable<StorageRoot> roots)
    {
        var drives = new List<DriveSpace>();
        foreach (string drive in roots.Select(r => Path.GetPathRoot(r.FullPath)!).Distinct(StringComparer.OrdinalIgnoreCase))
        {
            try
            {
                var info = new DriveInfo(drive);
                drives.Add(new DriveSpace(drive, info.AvailableFreeSpace, info.TotalSize));
            }
            catch (Exception e) when (e is IOException or ArgumentException or UnauthorizedAccessException) { }
        }
        return drives;
    }

    private static List<string> OrphanedWorktrees(string repoRoot, GitStatistics git)
    {
        string folder = Path.Combine(repoRoot, ".claude", "worktrees");
        if (!Directory.Exists(folder)) return [];
        var registered = new HashSet<string>(
            git.Worktrees.Select(w => Path.GetFullPath(w.Path).TrimEnd('\\', '/')), StringComparer.OrdinalIgnoreCase);
        return Directory.EnumerateDirectories(folder)
            .Where(d => !registered.Contains(Path.GetFullPath(d).TrimEnd('\\', '/')))
            .Select(d => Path.GetRelativePath(repoRoot, d).Replace('\\', '/') + "/")
            .ToList();
    }

    private sealed class Walker(StorageRules rules, RepoGit git)
    {
        private readonly List<(int Index, StorageRule Rule, PathGlob Glob)> _rules =
            rules.Rules.Select((r, i) => (i, r, new PathGlob(r.Glob))).ToList();

        public List<StorageItem> Items { get; } = [];

        public Dictionary<(string, string, RollupBucket), Rollup> Rollups { get; } = new();

        public Dictionary<(string, string), DirectoryTotal> Tops { get; } = new();

        public List<string> Unreadable { get; } = [];

        public void Walk(StorageRoot root, DirectoryInfo directory, string relative, StorageItem? owner, string? top)
        {
            List<FileSystemInfo> entries;
            try { entries = directory.EnumerateFileSystemInfos("*", Options).ToList(); }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                Unreadable.Add($"{root.Name}:{relative}");
                return;
            }

            foreach (FileSystemInfo entry in entries)
            {
                bool isDirectory = entry is DirectoryInfo;
                string childRelative = relative + entry.Name + (isDirectory ? "/" : "");
                // The repository totals per top-level entry; a temp root is one total.
                string childTop = !root.IsRepo ? "(all)" : top ?? (isDirectory ? childRelative : "./");
                StorageItem? matched = Match(root, childRelative, isDirectory, owner);
                StorageItem? childOwner = matched ?? owner;

                if (entry.Attributes.HasFlag(FileAttributes.ReparsePoint))
                {
                    if (childOwner is not null) childOwner.Links++;
                    continue;
                }

                if (entry is DirectoryInfo sub)
                {
                    if (childOwner is not null)
                    {
                        childOwner.OwnedDirectories.Add(childRelative);
                        childOwner.Touch(sub.LastWriteTimeUtc);
                    }
                    Walk(root, sub, childRelative, childOwner, childTop);
                    continue;
                }

                var file = (FileInfo)entry;
                long length = file.Length;
                DateTime written = file.LastWriteTimeUtc;
                Total(root.Name, childTop, length, written);

                bool tracked = root.IsRepo && git.IsTracked(childRelative);
                bool ignored = !root.IsRepo || git.IsIgnored(childRelative);
                if (childOwner is not null)
                {
                    childOwner.OwnedFiles.Add(childRelative);
                    childOwner.Bytes += length;
                    childOwner.Files++;
                    childOwner.Touch(written);
                    if (tracked) childOwner.TrackedFiles++;
                    else if (!ignored) childOwner.UnignoredFiles++;
                    continue;
                }

                RollupBucket bucket = !root.IsRepo ? RollupBucket.TempNoRule
                    : tracked ? RollupBucket.Tracked
                    : ignored ? RollupBucket.IgnoredNoRule
                    : RollupBucket.Untracked;
                var key = (root.Name, childTop, bucket);
                if (!Rollups.TryGetValue(key, out Rollup? rollup))
                    Rollups[key] = rollup = new Rollup { Root = root.Name, Top = childTop, Bucket = bucket };
                rollup.Bytes += length;
                rollup.Files++;
                if (written > rollup.NewestUtc) rollup.NewestUtc = written;
            }
        }

        /// <summary>Hands the rule's category to all but the newest few of each keepNewest group.</summary>
        public void KeepNewest()
        {
            var groups = Items
                .Where(i => _rules[i.RuleIndex].Rule.KeepNewest is not null)
                .GroupBy(i => (i.RuleIndex, Parent(i.Path)));
            foreach (var group in groups)
            {
                int keep = _rules[group.Key.RuleIndex].Rule.KeepNewest!.Value;
                foreach (StorageItem item in group
                    .OrderByDescending(i => i.NewestUtc)
                    .ThenByDescending(i => i.Path, StringComparer.Ordinal)
                    .Take(keep))
                {
                    item.Category = Category.Keep;
                    item.Reason += $" (one of the newest {keep}, kept)";
                }
            }
        }

        private StorageItem? Match(StorageRoot root, string relative, bool isDirectory, StorageItem? owner)
        {
            foreach ((int index, StorageRule rule, PathGlob glob) in _rules)
            {
                if (!string.Equals(rule.Root, root.Name, StringComparison.OrdinalIgnoreCase)) continue;
                if (!glob.IsMatch(relative, isDirectory)) continue;

                var item = new StorageItem
                {
                    Root = root.Name,
                    Path = relative,
                    IsDirectory = isDirectory,
                    Category = rule.Category,
                    Kind = rule.Kind,
                    Reason = rule.Reason,
                    Drive = Path.GetPathRoot(root.FullPath)!,
                    RuleIndex = index,
                };
                Items.Add(item);
                owner?.Nested.Add(relative);
                return item;
            }
            return null;
        }

        private void Total(string root, string top, long length, DateTime written)
        {
            if (!Tops.TryGetValue((root, top), out DirectoryTotal? total))
                Tops[(root, top)] = total = new DirectoryTotal { Root = root, Name = top };
            total.Bytes += length;
            total.Files++;
            if (written > total.NewestUtc) total.NewestUtc = written;
        }

        private static string Parent(string path)
        {
            int slash = path.TrimEnd('/').LastIndexOf('/');
            return slash < 0 ? "" : path[..(slash + 1)];
        }
    }
}
