using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;

namespace StorageAudit;

internal static class Format
{
    public static string Size(long bytes) => bytes switch
    {
        >= 1L << 30 => (bytes / (double)(1L << 30)).ToString("0.00", CultureInfo.InvariantCulture) + " GB",
        >= 1L << 20 => (bytes / (double)(1L << 20)).ToString("0.0", CultureInfo.InvariantCulture) + " MB",
        >= 1L << 10 => (bytes / (double)(1L << 10)).ToString("0.0", CultureInfo.InvariantCulture) + " KB",
        _ => bytes + " B",
    };

    public static string When(DateTime utc) => utc == default ? "—" : utc.ToLocalTime().ToString("yyyy-MM-dd HH:mm", CultureInfo.InvariantCulture);

    public static string AgeDays(DateTime utc, DateTime nowUtc) => utc == default ? "—" : ((nowUtc - utc).TotalDays).ToString("0.0", CultureInfo.InvariantCulture);

    public static string Name(Category c) => c.ToString().ToLowerInvariant();

    public static string Cell(string text) => text.Replace("|", "\\|");
}

/// <summary>Renders the profile and the plan as markdown for review, and the profile as JSON.</summary>
public static class ReportWriter
{
    private const int LargestItems = 40;

    public static void WriteProfile(StorageProfile profile, string directory)
    {
        Directory.CreateDirectory(directory);
        File.WriteAllText(Path.Combine(directory, "profile.json"), JsonSerializer.Serialize(profile, StorageRules.JsonOptions));
        File.WriteAllText(Path.Combine(directory, "profile.md"), ProfileMarkdown(profile));
    }

    public static void WritePlan(StoragePlan plan, string directory)
    {
        Directory.CreateDirectory(directory);
        plan.Write(Path.Combine(directory, "plan.json"));
        File.WriteAllText(Path.Combine(directory, "plan.md"), PlanMarkdown(plan, Path.Combine(directory, "plan.json")));
    }

    public static string ProfileMarkdown(StorageProfile p)
    {
        var sb = new StringBuilder();
        sb.Append($"# Storage profile: {p.RepoRoot}\n\n");
        sb.Append($"Created {Format.When(p.CreatedUtc)}. Rules: `{p.RulesPath ?? "(in memory)"}`. ");
        sb.Append($"Anything written within {p.RecentWindowHours:0.#} hours (longer where a rule says) counts as active.\n\n");

        sb.Append("## Drives\n\n| Drive | Free | Total |\n|---|---:|---:|\n");
        foreach (DriveSpace d in p.Drives) sb.Append($"| {d.Drive} | {Format.Size(d.FreeBytes)} | {Format.Size(d.TotalBytes)} |\n");
        foreach (StorageRoot r in p.Roots) sb.Append($"\nRoot `{r.Name}`: `{r.FullPath}`{(p.MissingRoots.Contains(r.Name) ? " (missing)" : "")}");
        sb.Append("\n\n");

        sb.Append("## By category\n\n| Category | Size | Items | Active | By drive |\n|---|---:|---:|---:|---|\n");
        foreach (Category c in Enum.GetValues<Category>())
        {
            List<StorageItem> items = p.Items.Where(i => i.Category == c).ToList();
            string drives = string.Join(", ", items.GroupBy(i => i.Drive).Select(g => $"{g.Key} {Format.Size(g.Sum(i => i.Bytes))}"));
            sb.Append($"| {Format.Name(c)} | {Format.Size(items.Sum(i => i.Bytes))} | {items.Count} | {Format.Size(items.Where(i => i.Active).Sum(i => i.Bytes))} | {drives} |\n");
        }
        foreach (var bucket in p.Rollups.GroupBy(r => r.Bucket))
            sb.Append($"| *{BucketName(bucket.Key)}* | {Format.Size(bucket.Sum(r => r.Bytes))} | — | — | |\n");
        sb.Append('\n');

        sb.Append("## By top-level entry\n\n| Root | Entry | Size | Files | Newest |\n|---|---|---:|---:|---|\n");
        foreach (DirectoryTotal t in p.TopDirectories)
            sb.Append($"| {t.Root} | `{Format.Cell(t.Name)}` | {Format.Size(t.Bytes)} | {t.Files} | {Format.When(t.NewestUtc)} |\n");
        sb.Append('\n');

        sb.Append($"## Largest items\n\nThe {LargestItems} largest classified items. Age is days since the newest write.\n\n");
        sb.Append("| Size | Files | Newest | Age | Category | Kind | Path | Notes |\n|---:|---:|---|---:|---|---|---|---|\n");
        foreach (StorageItem i in p.Items.OrderByDescending(i => i.Bytes).Take(LargestItems))
            sb.Append($"| {Format.Size(i.Bytes)} | {i.Files} | {Format.When(i.NewestUtc)} | {Format.AgeDays(i.NewestUtc, p.CreatedUtc)} | "
                    + $"{Format.Name(i.Category)} | {i.Kind} | `{i.Root}:{Format.Cell(i.Path)}` | {Format.Cell(Notes(i))} |\n");
        sb.Append('\n');

        sb.Append("## Files no rule claims\n\nReport-only; nothing acts on these. An ignored bulk here is a candidate for a new rule.\n\n");
        sb.Append("| Root | Entry | Bucket | Size | Files | Newest |\n|---|---|---|---:|---:|---|\n");
        foreach (Rollup r in p.Rollups.Take(30))
            sb.Append($"| {r.Root} | `{Format.Cell(r.Top)}` | {BucketName(r.Bucket)} | {Format.Size(r.Bytes)} | {r.Files} | {Format.When(r.NewestUtc)} |\n");
        sb.Append('\n');

        sb.Append("## Git\n\n| count-objects | value |\n|---|---|\n");
        foreach (var kv in p.Git.CountObjects) sb.Append($"| {kv.Key} | {kv.Value} |\n");
        sb.Append("\nSizes are in KiB.\n\n| Worktree | Branch | Prunable | Locked |\n|---|---|---|---|\n");
        foreach (Worktree w in p.Git.Worktrees) sb.Append($"| `{w.Path}` | {w.Branch ?? "—"} | {w.Prunable} | {w.Locked} |\n");
        sb.Append($"\nReflog: {p.Git.ReflogEntries} HEAD entries. Orphaned worktree folders: "
                + (p.OrphanedWorktreeFolders.Count == 0 ? "none" : string.Join(", ", p.OrphanedWorktreeFolders.Select(f => $"`{f}`"))) + ".\n");

        if (p.Unreadable.Count > 0)
            sb.Append($"\n## Unreadable\n\n{string.Join("\n", p.Unreadable.Select(u => $"- `{u}`"))}\n");
        return sb.ToString();
    }

    public static string PlanMarkdown(StoragePlan plan, string planJsonPath)
    {
        var sb = new StringBuilder();
        sb.Append($"# Storage plan: {plan.RepoRoot}\n\n");
        sb.Append($"Created {Format.When(plan.CreatedUtc)}. Review it; to skip an action, set its `\"enabled\": false` in `{planJsonPath}`. Then:\n\n");
        sb.Append($"```powershell\nscripts/storage-audit.ps1 apply '{planJsonPath}'            # dry run\n");
        sb.Append($"scripts/storage-audit.ps1 apply '{planJsonPath}' -Execute   # add -ArchiveTo <path> for archives\n```\n\n");

        sb.Append("## Frees, by drive\n\n| Drive | Delete | Archive |\n|---|---:|---:|\n");
        foreach (var g in plan.Actions.GroupBy(a => a.Drive))
            sb.Append($"| {g.Key} | {Format.Size(g.Where(a => a.Action == PlanActionKind.Delete).Sum(a => a.Bytes))} | "
                    + $"{Format.Size(g.Where(a => a.Action == PlanActionKind.Archive).Sum(a => a.Bytes))} |\n");
        sb.Append('\n');

        sb.Append("## Actions\n\n| # | Action | Size | Files | Newest | Path | Reason |\n|---:|---|---:|---:|---|---|---|\n");
        foreach (PlanAction a in plan.Actions)
            sb.Append($"| {a.Id} | {a.Action.ToString().ToLowerInvariant()} | {Format.Size(a.Bytes)} | {a.Files} | {Format.When(a.NewestUtc)} | "
                    + $"`{a.Root}:{Format.Cell(a.Path)}` | {Format.Cell(a.Reason)} |\n");
        sb.Append('\n');

        sb.Append("## Held back\n\n| Path | Category | Size | Why |\n|---|---|---:|---|\n");
        foreach (HeldItem h in plan.Held)
            sb.Append($"| `{h.Root}:{Format.Cell(h.Path)}` | {Format.Name(h.Category)} | {Format.Size(h.Bytes)} | {Format.Cell(h.Why)} |\n");
        sb.Append('\n');

        sb.Append("## Git (advice only; the tool never runs these)\n\n");
        sb.Append(plan.GitAdvice.Count == 0 ? "Nothing to propose.\n" : string.Join("\n", plan.GitAdvice.Select(a => $"- {a}")) + "\n");
        return sb.ToString();
    }

    private static string Notes(StorageItem i)
    {
        var notes = new List<string>();
        if (i.Active) notes.Add("active");
        if (i.TrackedFiles > 0) notes.Add($"{i.TrackedFiles} tracked");
        if (i.UnignoredFiles > 0) notes.Add($"{i.UnignoredFiles} not ignored");
        if (i.Links > 0) notes.Add($"{i.Links} link(s)");
        if (i.Nested.Count > 0) notes.Add($"carves out {i.Nested.Count}");
        notes.Add(i.Reason);
        return string.Join("; ", notes);
    }

    private static string BucketName(RollupBucket b) => b switch
    {
        RollupBucket.Tracked => "tracked (keep)",
        RollupBucket.Untracked => "untracked, not ignored (never acted on)",
        RollupBucket.IgnoredNoRule => "ignored, no rule",
        _ => "temp, no rule",
    };
}
