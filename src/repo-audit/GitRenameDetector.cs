using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace RepoAudit;

public sealed record GitRename(string OldAbsPath, string NewAbsPath);

public static class GitRenameDetector
{
    public static IReadOnlyList<GitRename> TryGetRenames(string repoRoot)
    {
        if (!IsGitRepo(repoRoot)) return Array.Empty<GitRename>();

        // old-relative → new-relative; later entries win (more recent diff wins)
        var seen = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

        // Last commit
        RunGitDiff(repoRoot, "diff --name-status HEAD~1 HEAD", seen);
        // Staged changes
        RunGitDiff(repoRoot, "diff --cached --name-status", seen);
        // Working-tree changes (with rename detection)
        RunGitDiff(repoRoot, "diff --find-renames --name-status HEAD", seen);

        return seen
            .Select(kv => new GitRename(
                Path.GetFullPath(Path.Combine(repoRoot, kv.Key)),
                Path.GetFullPath(Path.Combine(repoRoot, kv.Value))))
            .ToList()
            .AsReadOnly();
    }

    public static IReadOnlyList<Violation> DetectRenameGaps(
        FileOwnershipResult ownership,
        IReadOnlyList<GitRename> renames,
        string repoRoot)
    {
        if (renames.Count == 0 || ownership.BrokenIncludes.Count == 0)
            return Array.Empty<Violation>();

        var oldToNew = renames.ToDictionary(r => r.OldAbsPath, r => r.NewAbsPath,
            StringComparer.OrdinalIgnoreCase);
        var orphanedSet = new HashSet<string>(ownership.OrphanedFiles,
            StringComparer.OrdinalIgnoreCase);

        var violations = new List<Violation>();

        foreach (BrokenInclude broken in ownership.BrokenIncludes)
        {
            if (!oldToNew.TryGetValue(broken.MissingAbsPath, out string? newAbsPath)) continue;

            string oldRel = Path.GetRelativePath(repoRoot, broken.MissingAbsPath).Replace('\\', '/');
            string newRel = Path.GetRelativePath(repoRoot, newAbsPath).Replace('\\', '/');

            string hint = orphanedSet.Contains(newAbsPath)
                ? "new path is unclaimed — update the glob or <Compile Include>"
                : "new path is already claimed by another project";

            violations.Add(new Violation(
                ViolationKind.RenameGap,
                broken.AssemblyName,
                $"Broken include '{oldRel}' — git shows rename to '{newRel}' ({hint})",
                FilePath: broken.MissingAbsPath));
        }

        return violations.AsReadOnly();
    }

    private static bool IsGitRepo(string repoRoot)
    {
        try
        {
            using Process? p = Start(repoRoot, "rev-parse --git-dir");
            p?.WaitForExit(3000);
            return p?.ExitCode == 0;
        }
        catch { return false; }
    }

    private static void RunGitDiff(string repoRoot, string args,
        Dictionary<string, string> seen)
    {
        try
        {
            using Process? p = Start(repoRoot, args);
            if (p == null) return;
            string output = p.StandardOutput.ReadToEnd();
            p.WaitForExit(5000);
            if (p.ExitCode != 0) return;

            foreach (string line in output.Split('\n', StringSplitOptions.RemoveEmptyEntries))
            {
                // Format: R<similarity>\t<old>\t<new>
                if (line.Length < 2 || (line[0] != 'R' && line[0] != 'r')) continue;
                string[] parts = line.Split('\t');
                if (parts.Length < 3) continue;
                seen[parts[1].Trim()] = parts[2].Trim('\r', '\n', ' ');
            }
        }
        catch { /* git unavailable or repo error — skip silently */ }
    }

    private static Process? Start(string workDir, string args)
    {
        var psi = new ProcessStartInfo("git", args)
        {
            WorkingDirectory       = workDir,
            RedirectStandardOutput = true,
            RedirectStandardError  = true,
            UseShellExecute        = false
        };
        return Process.Start(psi);
    }
}
