using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace StorageAudit;

/// <summary>A registered worktree, from <c>git worktree list --porcelain</c>.</summary>
public sealed record Worktree(string Path, string? Branch, bool Prunable, bool Locked);

/// <summary>Git's object store and worktrees, as the profile reports them.</summary>
public sealed record GitStatistics(
    IReadOnlyDictionary<string, string> CountObjects,
    IReadOnlyList<Worktree> Worktrees,
    int ReflogEntries);

/// <summary>
/// What git knows about the repository's files: which are tracked and which are ignored. The
/// audit acts only on ignored, untracked files, so anything git does not ignore (the user's notes,
/// work in progress) stays out of reach. Git is required; a failed query throws rather than
/// leaving the audit to guess.
/// </summary>
public sealed class RepoGit
{
    private readonly HashSet<string> _trackedFiles;
    private readonly HashSet<string> _trackedDirectories;
    private readonly HashSet<string> _ignored;

    private RepoGit(HashSet<string> trackedFiles, HashSet<string> trackedDirectories, HashSet<string> ignored)
    {
        _trackedFiles = trackedFiles;
        _trackedDirectories = trackedDirectories;
        _ignored = ignored;
    }

    public static RepoGit Read(string repoRoot)
    {
        var tracked = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        var trackedDirs = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (string path in Split0(Run(repoRoot, "ls-files", "-z")))
        {
            tracked.Add(path);
            for (int slash = path.IndexOf('/'); slash >= 0; slash = path.IndexOf('/', slash + 1))
                trackedDirs.Add(path[..(slash + 1)]);
        }

        // Traditional mode collapses a wholly ignored folder to one "!! dir/" entry.
        var ignored = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        string[] status = Split0(Run(repoRoot, "status", "--porcelain=v1", "-z", "--ignored=traditional", "--untracked-files=normal"));
        for (int i = 0; i < status.Length; i++)
        {
            string entry = status[i];
            if (entry.Length < 4) continue;
            if (entry[0] is 'R' or 'C' || entry[1] is 'R' or 'C') i++;   // a rename or copy carries its source next
            if (entry.StartsWith("!! ", StringComparison.Ordinal)) ignored.Add(entry[3..]);
        }
        return new RepoGit(tracked, trackedDirs, ignored);
    }

    public bool IsTracked(string relativeFile) => _trackedFiles.Contains(relativeFile);

    /// <summary>Whether a folder (relative, ending in <c>/</c>) holds any tracked file.</summary>
    public bool HoldsTracked(string relativeDirectory) => _trackedDirectories.Contains(relativeDirectory);

    /// <summary>Whether git ignores a relative file, itself or through an ignored folder above it.</summary>
    public bool IsIgnored(string relativeFile)
    {
        if (_ignored.Contains(relativeFile)) return true;
        for (int slash = relativeFile.IndexOf('/'); slash >= 0; slash = relativeFile.IndexOf('/', slash + 1))
            if (_ignored.Contains(relativeFile[..(slash + 1)])) return true;
        return false;
    }

    public static GitStatistics ReadStatistics(string repoRoot)
    {
        var counts = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (string line in Lines(Run(repoRoot, "count-objects", "-v")))
        {
            int colon = line.IndexOf(':');
            if (colon > 0) counts[line[..colon].Trim()] = line[(colon + 1)..].Trim();
        }

        var worktrees = new List<Worktree>();
        string? path = null, branch = null;
        bool prunable = false, locked = false;
        // Blocks are separated by blank lines, so this split keeps them.
        string[] porcelain = Run(repoRoot, "worktree", "list", "--porcelain")
            .Split('\n').Select(l => l.TrimEnd('\r')).ToArray();
        foreach (string line in porcelain.Append(""))
        {
            if (line.Length == 0)
            {
                if (path is not null) worktrees.Add(new Worktree(path, branch, prunable, locked));
                path = branch = null;
                prunable = locked = false;
            }
            else if (line.StartsWith("worktree ", StringComparison.Ordinal)) path = line[9..];
            else if (line.StartsWith("branch ", StringComparison.Ordinal)) branch = line[7..];
            else if (line.StartsWith("prunable", StringComparison.Ordinal)) prunable = true;
            else if (line.StartsWith("locked", StringComparison.Ordinal)) locked = true;
        }

        // HEAD's reflog is one line per entry; reading it saves a git call. A worktree's .git is a
        // file, and an empty repository has no log yet: both count as zero.
        string headLog = Path.Combine(repoRoot, ".git", "logs", "HEAD");
        int reflog = File.Exists(headLog) ? File.ReadLines(headLog).Count(l => l.Length > 0) : 0;

        return new GitStatistics(counts, worktrees, reflog);
    }

    private static string Run(string repoRoot, params string[] args)
    {
        var psi = new ProcessStartInfo("git")
        {
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            UseShellExecute = false,
            CreateNoWindow = true,
        };
        psi.ArgumentList.Add("-C");
        psi.ArgumentList.Add(repoRoot);
        psi.ArgumentList.Add("-c");
        psi.ArgumentList.Add("core.quotepath=off");
        foreach (string a in args) psi.ArgumentList.Add(a);

        using Process process = Process.Start(psi)
            ?? throw new InvalidOperationException("git could not be started.");
        var stderr = process.StandardError.ReadToEndAsync();
        string stdout = process.StandardOutput.ReadToEnd();
        process.WaitForExit();
        if (process.ExitCode != 0)
            throw new InvalidOperationException(
                $"git {string.Join(' ', args)} failed in {repoRoot}: {stderr.Result.Trim()}");
        return stdout;
    }

    private static string[] Split0(string output)
        => output.Split('\0', StringSplitOptions.RemoveEmptyEntries);

    private static string[] Lines(string output)
        => output.Split('\n').Select(l => l.TrimEnd('\r')).Where(l => l.Length > 0).ToArray();
}
