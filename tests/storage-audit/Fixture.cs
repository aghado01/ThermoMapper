using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using StorageAudit;
using Xunit;

namespace StorageAudit.Tests;

/// <summary>
/// A throwaway git repository under the system temp folder, with a scratch folder beside it for the
/// rules file, plan, ledger, archive and temp roots, so nothing a test writes lands in this repository.
/// It starts with a committed README and a .gitignore that ignores <c>artifacts/</c>, <c>lake/</c> and
/// <c>issues/scratch/</c>.
/// </summary>
internal sealed class Fixture : IDisposable
{
    private static readonly EnumerationOptions NoLinks = new()
    {
        RecurseSubdirectories = true,
        AttributesToSkip = FileAttributes.ReparsePoint,
    };

    // Git calls are slow on Windows, so one repository is initialized per test run and copied per test.
    private static readonly Lazy<string> Template = new(CreateTemplate);

    public Fixture()
    {
        Base = Path.Combine(Path.GetTempPath(), "storage-audit-tests", Guid.NewGuid().ToString("N")[..12]);
        Repo = Path.Combine(Base, "repo");
        Scratch = Path.Combine(Base, "scratch");
        Directory.CreateDirectory(Scratch);
        CopyTree(Template.Value, Repo);
    }

    private static string CreateTemplate()
    {
        string root = Path.Combine(Path.GetTempPath(), "storage-audit-tests", "template-" + Guid.NewGuid().ToString("N")[..12]);
        Directory.CreateDirectory(root);
        AppDomain.CurrentDomain.ProcessExit += (_, _) => Delete(root);

        RunGit(root, "init", "-q");
        RunGit(root, "config", "user.email", "fixture@example.invalid");
        RunGit(root, "config", "user.name", "fixture");
        File.WriteAllText(Path.Combine(root, "README.md"), "fixture\n");
        File.WriteAllText(Path.Combine(root, ".gitignore"), "artifacts/\nlake/\nissues/scratch/\n");
        RunGit(root, "add", "--", "README.md", ".gitignore");
        RunGit(root, "commit", "-q", "-m", "fixture");
        return root;
    }

    private static void CopyTree(string source, string destination)
    {
        Directory.CreateDirectory(destination);
        foreach (string d in Directory.EnumerateDirectories(source, "*", NoLinks))
            Directory.CreateDirectory(Path.Combine(destination, Path.GetRelativePath(source, d)));
        foreach (string f in Directory.EnumerateFiles(source, "*", NoLinks))
            File.Copy(f, Path.Combine(destination, Path.GetRelativePath(source, f)));
    }

    public string Base { get; }

    public string Repo { get; }

    public string Scratch { get; }

    public DateTime Now { get; } = DateTime.UtcNow;

    public string RulesPath => Path.Combine(Scratch, "rules.json");

    public string PlanPath => Path.Combine(Scratch, "plan.json");

    public string LedgerPath => Path.Combine(Scratch, "ledger.jsonl");

    public string Full(string relative) => Path.Combine(Repo, relative.Replace('/', Path.DirectorySeparatorChar));

    /// <summary>Writes a file of the given size under the repository (or an absolute path) and returns its full path.</summary>
    public string Write(string relative, int bytes = 1024, byte fill = (byte)'x')
    {
        string full = Path.IsPathRooted(relative) ? relative : Full(relative);
        Directory.CreateDirectory(Path.GetDirectoryName(full)!);
        File.WriteAllBytes(full, Enumerable.Repeat(fill, bytes).ToArray());
        return full;
    }

    public void Commit(params string[] relatives)
    {
        Git(["add", "-f", "--", .. relatives]);
        Git("commit", "-q", "-m", "fixture");
    }

    /// <summary>Backdates every file and folder (git's own excepted) under a root, deepest first.</summary>
    public void Age(double days, string? root = null)
    {
        root ??= Repo;
        DateTime then = DateTime.UtcNow.AddDays(-days);
        string git = Path.Combine(Repo, ".git");
        foreach (string f in Directory.EnumerateFiles(root, "*", NoLinks).Where(p => !p.StartsWith(git, StringComparison.OrdinalIgnoreCase)))
            File.SetLastWriteTimeUtc(f, then);
        foreach (string d in Directory.EnumerateDirectories(root, "*", NoLinks)
                     .Where(p => !p.StartsWith(git, StringComparison.OrdinalIgnoreCase))
                     .OrderByDescending(p => p.Length))
            Directory.SetLastWriteTimeUtc(d, then);
        Directory.SetLastWriteTimeUtc(root, then);
    }

    public static StorageRule Rule(string glob, Category category, int? keepNewest = null, string root = StorageRules.RepoRootName)
        => new() { Glob = glob, Category = category, Reason = $"test rule {glob}", Kind = "test", KeepNewest = keepNewest, Root = root };

    public static StorageRules Rules(params StorageRule[] rules)
        => new() { Rules = rules, Protect = ["issues/"] };

    public StorageRules SaveRules(StorageRules rules)
    {
        File.WriteAllText(RulesPath, System.Text.Json.JsonSerializer.Serialize(rules, StorageRules.JsonOptions));
        return StorageRules.Load(RulesPath);
    }

    public StorageProfile Profile(StorageRules rules) => StorageProfiler.Run(Repo, rules, Now, RulesPath);

    /// <summary>Saves the rules, profiles, plans and writes the plan file.</summary>
    public StoragePlan Plan(StorageRules rules, string? ownBinaryDirectory = null)
    {
        StorageRules saved = SaveRules(rules);
        StoragePlan plan = StoragePlanner.Build(Profile(saved), saved, ownBinaryDirectory);
        plan.Write(PlanPath);
        return plan;
    }

    public IReadOnlyList<ApplyResult> Apply(bool execute = false, string? archiveTo = null, string? expectedRepo = null,
        Action<PlanAction>? beforeAct = null)
        => PlanApplier.Apply(new ApplyOptions
        {
            PlanPath = PlanPath,
            RulesPath = RulesPath,
            LedgerPath = LedgerPath,
            Execute = execute,
            ArchiveTo = archiveTo,
            ExpectedRepoRoot = expectedRepo,
            BeforeAct = beforeAct,
        }, TextWriter.Null);

    /// <summary>Asserts a single result was refused and returns why.</summary>
    public static string Refused(IReadOnlyList<ApplyResult> results)
    {
        ApplyResult r = Assert.Single(results);
        Assert.Equal("refused", r.Outcome);
        return r.Detail;
    }

    public string Git(params string[] args) => RunGit(Repo, args);

    private static string RunGit(string repo, params string[] args)
    {
        var psi = new ProcessStartInfo("git")
        {
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            UseShellExecute = false,
            CreateNoWindow = true,
        };
        psi.ArgumentList.Add("-C");
        psi.ArgumentList.Add(repo);
        foreach (string a in args) psi.ArgumentList.Add(a);
        using Process p = Process.Start(psi)!;
        var err = p.StandardError.ReadToEndAsync();
        string output = p.StandardOutput.ReadToEnd();
        p.WaitForExit();
        if (p.ExitCode != 0) throw new InvalidOperationException($"git {string.Join(' ', args)}: {err.Result}");
        return output;
    }

    /// <summary>Creates a directory junction (no privilege needed), or returns false where mklink is unavailable.</summary>
    public static bool TryJunction(string link, string target)
    {
        if (!OperatingSystem.IsWindows()) return false;
        var psi = new ProcessStartInfo("cmd.exe") { UseShellExecute = false, CreateNoWindow = true, RedirectStandardOutput = true };
        psi.ArgumentList.Add("/c");
        psi.ArgumentList.Add("mklink");
        psi.ArgumentList.Add("/J");
        psi.ArgumentList.Add(link);
        psi.ArgumentList.Add(target);
        using Process p = Process.Start(psi)!;
        p.StandardOutput.ReadToEnd();
        p.WaitForExit();
        return p.ExitCode == 0 && Directory.Exists(link);
    }

    public void Dispose() => Delete(Base);

    /// <summary>Removes the links first (never following them), clears git's read-only objects, then deletes.</summary>
    private static void Delete(string root)
    {
        try
        {
            RemoveLinks(root);
            foreach (string f in Directory.EnumerateFiles(root, "*", NoLinks))
                File.SetAttributes(f, FileAttributes.Normal);
            Directory.Delete(root, recursive: true);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException) { }
    }

    private static void RemoveLinks(string directory)
    {
        foreach (string sub in Directory.EnumerateDirectories(directory, "*", new EnumerationOptions { AttributesToSkip = 0 }))
        {
            if (File.GetAttributes(sub).HasFlag(FileAttributes.ReparsePoint)) Directory.Delete(sub, recursive: false);
            else RemoveLinks(sub);
        }
    }
}
