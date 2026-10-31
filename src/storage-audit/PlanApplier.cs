using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text.Json;

namespace StorageAudit;

public sealed record ApplyOptions
{
    public required string PlanPath { get; init; }

    public required string RulesPath { get; init; }

    public required string LedgerPath { get; init; }

    /// <summary>Act; without it, <c>apply</c> is a dry run that reports what it would do.</summary>
    public bool Execute { get; init; }

    /// <summary>Where archivable items go, outside every root; required to execute an archive.</summary>
    public string? ArchiveTo { get; init; }

    public string? OwnBinaryDirectory { get; init; }

    /// <summary>The repository the plan must be for; a plan naming another root is refused whole.</summary>
    public string? ExpectedRepoRoot { get; init; }

    public DateTime NowUtc { get; init; } = DateTime.UtcNow;

    /// <summary>A test seam: runs after an action is vetted and before it acts.</summary>
    internal Action<PlanAction>? BeforeAct { get; init; }
}

/// <summary>What happened to one plan action; each is also one line of the ledger.</summary>
public sealed record ApplyResult(
    DateTime TimeUtc,
    string Mode,
    string Plan,
    int Id,
    PlanActionKind Action,
    string Root,
    string Path,
    string Outcome,
    string Detail,
    long Bytes,
    int Files,
    string? Destination);

public sealed record ManifestFile(string Path, long Bytes, string Sha256);

public sealed record ManifestItem(string Root, string Path, string Destination, long Bytes, int Files, IReadOnlyList<ManifestFile> FileList);

/// <summary>The archive's record of what it holds and where each item came from, written before any source is removed.</summary>
public sealed record ArchiveManifest(DateTime CreatedUtc, string RepoRoot, string ArchiveTo, IReadOnlyList<ManifestItem> Items);

/// <summary>Checks an archived copy against its manifest: every file present, with its size and SHA-256.</summary>
public static class ArchiveVerifier
{
    public static List<string> Verify(IReadOnlyList<ManifestFile> files, string destinationBase)
    {
        var problems = new List<string>();
        foreach (ManifestFile f in files)
        {
            var copy = new FileInfo(Path.Combine(destinationBase, f.Path));
            if (!copy.Exists) problems.Add($"missing {f.Path}");
            else if (copy.Length != f.Bytes) problems.Add($"size of {f.Path}: {copy.Length} != {f.Bytes}");
            else if (!string.Equals(HashFile(copy.FullName), f.Sha256, StringComparison.Ordinal)) problems.Add($"hash of {f.Path}");
        }
        return problems;
    }

    public static string HashFile(string path)
    {
        using FileStream stream = File.OpenRead(path);
        return Convert.ToHexStringLower(SHA256.HashData(stream));
    }

    /// <summary>Copies a file, hashing the source as it streams, and returns its length and hash.</summary>
    internal static (long Bytes, string Sha256) CopyHashed(string source, string destination)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        using (FileStream input = File.OpenRead(source))
        using (FileStream output = new(destination, FileMode.CreateNew, FileAccess.Write))
        {
            byte[] buffer = new byte[1 << 20];
            int read;
            while ((read = input.Read(buffer, 0, buffer.Length)) > 0)
            {
                hash.AppendData(buffer, 0, read);
                output.Write(buffer, 0, read);
            }
        }
        File.SetLastWriteTimeUtc(destination, File.GetLastWriteTimeUtc(source));
        return (new FileInfo(destination).Length, Convert.ToHexStringLower(hash.GetHashAndReset()));
    }
}

/// <summary>
/// Executes a reviewed plan file and nothing else. It profiles afresh with the current rules and
/// re-checks every action before acting: the path stays inside its root, the rules still give the
/// action's category, nothing was written within the recent window, nothing tracked or un-ignored
/// is held, and the item has not grown since the plan. A dry run unless told to execute; every
/// action, refusal and dry-run line goes to the ledger.
/// </summary>
public static class PlanApplier
{
    private static readonly JsonSerializerOptions LedgerOptions = new(StorageRules.JsonOptions) { WriteIndented = false };

    public static IReadOnlyList<ApplyResult> Apply(ApplyOptions options, TextWriter log)
    {
        string planPath = Path.GetFullPath(options.PlanPath);
        StoragePlan plan = StoragePlan.Read(planPath);
        if (options.ExpectedRepoRoot is { } expected
            && !Path.GetFullPath(expected).TrimEnd('\\', '/').Equals(Path.GetFullPath(plan.RepoRoot).TrimEnd('\\', '/'), StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException($"The plan is for {plan.RepoRoot}, not {expected}.");
        StorageRules rules = StorageRules.Load(options.RulesPath);
        StorageProfile fresh = StorageProfiler.Run(plan.RepoRoot, rules, options.NowUtc, options.RulesPath);
        if (plan.RulesSha256 is not null && fresh.RulesSha256 is not null && plan.RulesSha256 != fresh.RulesSha256)
            log.WriteLine("note: the rules file changed since the plan; every action is re-checked against the current rules.");

        string mode = options.Execute ? "execute" : "dry-run";
        string? archiveTo = options.ArchiveTo is null ? null : Path.GetFullPath(options.ArchiveTo).TrimEnd('\\', '/');
        string manifestPath = archiveTo is null ? "" : Path.Combine(archiveTo,
            $"storage-audit-manifest-{options.NowUtc.ToLocalTime():yyyyMMdd_HHmmss}.json");
        var manifestItems = new List<ManifestItem>();
        var results = new List<ApplyResult>();

        foreach (PlanAction action in plan.Actions)
        {
            ApplyResult Result(string outcome, string detail, long bytes = 0, int files = 0, string? destination = null)
                => new(DateTime.UtcNow, mode, planPath, action.Id, action.Action, action.Root, action.Path,
                       outcome, detail, bytes, files, destination);

            ApplyResult result;
            if (!action.Enabled)
                result = Result("skipped", "disabled in the plan");
            else if (Vet(action, fresh, rules, options, archiveTo) is { } refusal)
                result = Result("refused", refusal.Why, refusal.Item?.Bytes ?? 0, refusal.Item?.Files ?? 0);
            else
            {
                options.BeforeAct?.Invoke(action);
                StorageItem item = fresh.Find(action.Root, action.Path)!;
                StorageRoot root = RootOf(fresh, action.Root)!;
                string? destination = action.Action == PlanActionKind.Archive && archiveTo is not null
                    ? Destination(archiveTo, root, item) : null;

                if (!options.Execute)
                    result = Result(action.Action == PlanActionKind.Delete ? "would-delete" : "would-archive",
                        action.Action == PlanActionKind.Archive && archiveTo is null ? "executing needs --archive-to" : "",
                        item.Bytes, item.Files, destination);
                else if (action.Action == PlanActionKind.Delete)
                {
                    Deletion d = DeleteOwned(root, item, options.NowUtc);
                    result = Result(d.Complete ? "deleted" : "partial", d.Describe(), d.Freed, item.Files);
                }
                else
                    result = Archive(root, item, destination!, plan, archiveTo!, manifestPath, manifestItems, options.NowUtc, Result);
            }

            results.Add(result);
            AppendLedger(options.LedgerPath, result);
            log.WriteLine($"{result.Outcome,-13} #{action.Id,-3} {Format.Size(result.Bytes),10}  {action.Root}:{action.Path}"
                        + (result.Detail.Length > 0 ? $"  ({result.Detail})" : "")
                        + (result.Destination is not null ? $"  -> {result.Destination}" : ""));
        }
        return results;
    }

    private sealed record Refusal(string Why, StorageItem? Item);

    private static Refusal? Vet(PlanAction action, StorageProfile fresh, StorageRules rules, ApplyOptions options, string? archiveTo)
    {
        if (string.IsNullOrWhiteSpace(action.Path) || Path.IsPathRooted(action.Path)
            || action.Path.Split('/', '\\').Contains(".."))
            return new Refusal("the path escapes its root", null);

        StorageRoot? root = RootOf(fresh, action.Root);
        if (root is null)
            return new Refusal($"the root '{action.Root}' is not in the rules; outside the repository only listed temp roots are reachable", null);
        string full = Path.GetFullPath(Path.Combine(root.FullPath, action.Path)).TrimEnd('\\', '/');
        if (!full.StartsWith(root.FullPath + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase))
            return new Refusal("the path escapes its root", null);

        StorageItem? item = fresh.Find(action.Root, action.Path);
        if (item is null)
            return new Refusal(File.Exists(full) || Directory.Exists(full)
                ? "no rule classifies this path now" : "the path no longer exists", null);
        if (item.Files == 0)
            return new Refusal("nothing left to act on", item);
        if (Refusals.Check(item, action.Action, rules.Protect, options.OwnBinaryDirectory, fresh) is { } why)
            return new Refusal(why, item);
        if (item.Bytes > action.Bytes + Math.Max(action.Bytes / 10, 1L << 20))
            return new Refusal($"grew from {Format.Size(action.Bytes)} to {Format.Size(item.Bytes)} since the plan; plan again", item);

        if (action.Action != PlanActionKind.Archive) return null;
        if (archiveTo is null)
            return options.Execute ? new Refusal("archiving needs --archive-to", item) : null;
        foreach (StorageRoot r in fresh.Roots)
            if (archiveTo.Equals(r.FullPath, StringComparison.OrdinalIgnoreCase)
                || archiveTo.StartsWith(r.FullPath + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase))
                return new Refusal($"the archive is inside the root '{r.Name}'; choose a path outside it", item);
        string destination = Destination(archiveTo, root, item);
        if (File.Exists(destination) || Directory.Exists(destination))
            return new Refusal($"the archive already holds {destination}", item);
        if (options.Execute && Path.GetPathRoot(archiveTo) is { } drive && Directory.Exists(drive))
        {
            long free = new DriveInfo(drive).AvailableFreeSpace;
            if (free < item.Bytes) return new Refusal($"the archive drive has {Format.Size(free)} free, the item needs {Format.Size(item.Bytes)}", item);
        }
        return null;
    }

    private static ApplyResult Archive(StorageRoot root, StorageItem item, string destination, StoragePlan plan,
        string archiveTo, string manifestPath, List<ManifestItem> manifestItems, DateTime nowUtc,
        Func<string, string, long, int, string?, ApplyResult> result)
    {
        string destinationBase = item.IsDirectory ? destination : Path.GetDirectoryName(destination)!;
        var files = new List<ManifestFile>();
        try
        {
            foreach (string relative in item.OwnedFiles)
            {
                string inside = item.IsDirectory ? relative[item.Path.Length..] : Path.GetFileName(relative);
                string target = Path.Combine(destinationBase, inside);
                Directory.CreateDirectory(Path.GetDirectoryName(target)!);
                (long bytes, string sha) = ArchiveVerifier.CopyHashed(Path.Combine(root.FullPath, relative), target);
                files.Add(new ManifestFile(inside.Replace('\\', '/'), bytes, sha));
            }
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            return result("failed", $"copy stopped ({e.Message}); the source is kept, the partial copy left at the destination",
                0, files.Count, destination);
        }

        List<string> problems = ArchiveVerifier.Verify(files, destinationBase);
        if (files.Count != item.OwnedFiles.Count)
            problems.Add($"copied {files.Count} of {item.OwnedFiles.Count} files");
        if (item.IsDirectory && Directory.EnumerateFiles(destination, "*", SearchOption.AllDirectories).Count() != files.Count)
            problems.Add("the destination holds a different number of files than were copied");
        if (problems.Count > 0)
            return result("failed", $"verification failed ({string.Join("; ", problems.Take(3))}); the source is kept",
                0, files.Count, destination);

        manifestItems.Add(new ManifestItem(item.Root, item.Path, destination, files.Sum(f => f.Bytes), files.Count, files));
        File.WriteAllText(manifestPath, JsonSerializer.Serialize(
            new ArchiveManifest(DateTime.UtcNow, plan.RepoRoot, archiveTo, manifestItems), StorageRules.JsonOptions));

        Deletion d = DeleteOwned(root, item, nowUtc);
        return result(d.Complete ? "archived" : "partial",
            d.Complete ? "verified" : $"verified and copied, but {d.Describe()}", d.Freed, files.Count, destination);
    }

    private sealed record Deletion(long Freed, int Failed, int Recent)
    {
        public bool Complete => Failed == 0 && Recent == 0;

        public string Describe() => string.Join("; ", new[]
        {
            Failed > 0 ? $"{Failed} file(s) could not be deleted (in use?)" : "",
            Recent > 0 ? $"{Recent} file(s) written since the profile were kept" : "",
        }.Where(s => s.Length > 0));
    }

    /// <summary>
    /// Deletes the files the item owns, then its folders left empty, deepest first; carved-out items
    /// stay. Each file's write time is read again first, so one written since the profile is kept.
    /// </summary>
    private static Deletion DeleteOwned(StorageRoot root, StorageItem item, DateTime nowUtc)
    {
        DateTime cutoff = nowUtc.AddHours(-item.WindowHours);
        long freed = 0;
        int failed = 0, recent = 0;
        foreach (string relative in item.OwnedFiles)
        {
            var file = new FileInfo(Path.Combine(root.FullPath, relative));
            try
            {
                if (!file.Exists) continue;
                if (file.LastWriteTimeUtc > cutoff) { recent++; continue; }
                long length = file.Length;
                if (file.IsReadOnly) file.IsReadOnly = false;
                file.Delete();
                freed += length;
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException) { failed++; }
        }

        foreach (string relative in item.OwnedDirectories.OrderByDescending(d => d.Length))
        {
            string directory = Path.Combine(root.FullPath, relative);
            try
            {
                if (Directory.Exists(directory) && !Directory.EnumerateFileSystemEntries(directory).Any())
                    Directory.Delete(directory, recursive: false);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException) { }
        }
        return new Deletion(freed, failed, recent);
    }

    private static string Destination(string archiveTo, StorageRoot root, StorageItem item)
    {
        string label = root.IsRepo ? Path.GetFileName(root.FullPath) : root.Name;
        return Path.GetFullPath(Path.Combine(archiveTo, label, item.Path)).TrimEnd('\\', '/');
    }

    private static StorageRoot? RootOf(StorageProfile profile, string name)
        => profile.Roots.FirstOrDefault(r => string.Equals(r.Name, name, StringComparison.OrdinalIgnoreCase));

    private static void AppendLedger(string ledgerPath, ApplyResult result)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(ledgerPath))!);
        string line = JsonSerializer.Serialize(result, LedgerOptions);
        File.AppendAllText(ledgerPath, line + "\n");
    }
}
