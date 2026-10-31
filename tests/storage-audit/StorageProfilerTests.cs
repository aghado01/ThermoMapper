using System;
using System.IO;
using System.Linq;
using StorageAudit;
using Xunit;
using static StorageAudit.Tests.Fixture;

namespace StorageAudit.Tests;

public class StorageProfilerTests
{
    [Fact]
    public void ADeeperRuleCarvesItsItemOutOfAShallowerOne()
    {
        using var f = new Fixture();
        f.Write("artifacts/loose.txt", 100);
        f.Write("artifacts/bin/a.dll", 1000);
        f.Write("artifacts/bin/sub/b.dll", 2000);
        f.Age(10);

        StorageProfile p = f.Profile(Rules(Rule("artifacts/", Category.Ask), Rule("artifacts/bin/", Category.Regenerable)));

        StorageItem parent = p.Find("repo", "artifacts/")!;
        StorageItem bin = p.Find("repo", "artifacts/bin/")!;
        Assert.Equal(100, parent.Bytes);
        Assert.Equal(1, parent.Files);
        Assert.Equal(["artifacts/bin/"], parent.Nested);
        Assert.Equal(3000, bin.Bytes);
        Assert.Equal(2, bin.Files);
        Assert.Equal(Category.Regenerable, bin.Category);
    }

    [Fact]
    public void AmongRulesMatchingOnePath_TheFirstListedWins()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");

        StorageProfile p = f.Profile(Rules(Rule("artifacts/bin/", Category.Keep), Rule("artifacts/*/", Category.Regenerable)));

        Assert.Equal(Category.Keep, p.Find("repo", "artifacts/bin/")!.Category);
        Assert.Single(p.Items);
    }

    [Fact]
    public void KeepNewestKeepsTheNewestOfEachParentsMatches()
    {
        using var f = new Fixture();
        foreach (string run in new[] { "r1", "r2", "r3" }) f.Write($"artifacts/runs/{run}/log.txt");
        f.Write("artifacts/other/r9/log.txt");
        f.Age(10);
        File.SetLastWriteTimeUtc(f.Full("artifacts/runs/r2/log.txt"), f.Now.AddDays(-2));

        StorageProfile p = f.Profile(Rules(Rule("artifacts/*/*/", Category.Regenerable, keepNewest: 1)));

        Assert.Equal(Category.Keep, p.Find("repo", "artifacts/runs/r2/")!.Category);
        Assert.Equal(Category.Regenerable, p.Find("repo", "artifacts/runs/r1/")!.Category);
        Assert.Equal(Category.Regenerable, p.Find("repo", "artifacts/runs/r3/")!.Category);
        Assert.Equal(Category.Keep, p.Find("repo", "artifacts/other/r9/")!.Category);
        Assert.Contains("kept", p.Find("repo", "artifacts/runs/r2/")!.Reason);
    }

    [Fact]
    public void FileGlobsMakeFileItems()
    {
        using var f = new Fixture();
        f.Write("artifacts/table.tsv", 300);
        f.Write("artifacts/notes.md", 50);

        StorageProfile p = f.Profile(Rules(Rule("artifacts/*.tsv", Category.Archivable)));

        StorageItem item = Assert.Single(p.Items);
        Assert.Equal("artifacts/table.tsv", item.Path);
        Assert.False(item.IsDirectory);
        Assert.Equal(300, item.Bytes);
    }

    [Fact]
    public void FilesNoRuleClaimsAreRolledUpByTrackedUntrackedAndIgnored()
    {
        using var f = new Fixture();
        f.Write("src/a.cs", 400);
        f.Commit("src/a.cs");
        f.Write("issues/note.md", 70);
        f.Write("lake/cache.bin", 900);

        StorageProfile p = f.Profile(Rules());

        Assert.Empty(p.Items);
        Assert.Equal(400, p.Rollups.Single(r => r.Top == "src/").Bytes);
        Assert.Equal(RollupBucket.Tracked, p.Rollups.Single(r => r.Top == "src/").Bucket);
        Assert.Equal(RollupBucket.Untracked, p.Rollups.Single(r => r.Top == "issues/").Bucket);
        Assert.Equal(RollupBucket.IgnoredNoRule, p.Rollups.Single(r => r.Top == "lake/").Bucket);
        Assert.Equal(900, p.TopDirectories.Single(t => t.Name == "lake/").Bytes);
    }

    [Fact]
    public void ItemsCountTheirTrackedAndUnignoredFiles()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/forced.dll");
        f.Commit("artifacts/bin/forced.dll");
        f.Write("artifacts/bin/plain.dll");
        f.Write("notes/draft.md");

        StorageProfile p = f.Profile(Rules(Rule("artifacts/bin/", Category.Regenerable), Rule("notes/", Category.Regenerable)));

        Assert.Equal(1, p.Find("repo", "artifacts/bin/")!.TrackedFiles);
        Assert.Equal(1, p.Find("repo", "notes/")!.UnignoredFiles);
    }

    [Fact]
    public void AnItemWrittenWithinTheWindowIsActive_FolderWritesIncluded()
    {
        using var f = new Fixture();
        f.Write("artifacts/old/a.txt");
        f.Write("artifacts/fresh/a.txt");
        f.Write("artifacts/emptied/a.txt");
        f.Age(10);
        File.SetLastWriteTimeUtc(f.Full("artifacts/fresh/a.txt"), f.Now.AddHours(-1));
        Directory.SetLastWriteTimeUtc(f.Full("artifacts/emptied"), f.Now.AddHours(-1));

        StorageProfile p = f.Profile(Rules(Rule("artifacts/*/", Category.Regenerable)));

        Assert.False(p.Find("repo", "artifacts/old/")!.Active);
        Assert.True(p.Find("repo", "artifacts/fresh/")!.Active);
        Assert.True(p.Find("repo", "artifacts/emptied/")!.Active);
    }

    [Fact]
    public void ARuleCanLengthenTheWindowForItsItems()
    {
        using var f = new Fixture();
        f.Write("artifacts/a/x.txt");
        f.Write("artifacts/b/x.txt");
        f.Age(3);

        StorageProfile p = f.Profile(Rules(
            Rule("artifacts/a/", Category.Regenerable),
            Rule("artifacts/b/", Category.Regenerable) with { RecentWindowHours = 168 }));

        Assert.False(p.Find("repo", "artifacts/a/")!.Active);
        Assert.True(p.Find("repo", "artifacts/b/")!.Active);
        Assert.Equal(168, p.Find("repo", "artifacts/b/")!.WindowHours);
    }

    [Fact]
    public void TempRootsAreWalkedUnderTheirOwnName()
    {
        using var f = new Fixture();
        string temp = Path.Combine(f.Scratch, "temp");
        f.Write(Path.Combine(temp, "session-1", "scratch.txt"), 500);
        StorageRules rules = Rules(Rule("*/", Category.Regenerable, root: "t")) with
        {
            TempRoots = [new TempRoot { Name = "t", Path = temp }],
        };

        StorageProfile p = f.Profile(rules);

        StorageItem item = p.Find("t", "session-1/")!;
        Assert.Equal(500, item.Bytes);
        Assert.Equal(Path.GetPathRoot(temp), item.Drive);
        Assert.Empty(p.MissingRoots);
        Assert.Equal("(all)", p.TopDirectories.Single(t => t.Root == "t").Name);
    }

    [Fact]
    public void ATempRootOverlappingTheRepositoryIsRejected()
    {
        using var f = new Fixture();
        foreach (string overlap in new[] { f.Repo, f.Base, Path.Combine(f.Repo, "artifacts") })
        {
            StorageRules rules = Rules(Rule("*/", Category.Regenerable, root: "t")) with
            {
                TempRoots = [new TempRoot { Name = "t", Path = overlap }],
            };
            Assert.Throws<InvalidDataException>(() => f.Profile(rules));
        }
    }

    [Fact]
    public void LinksAreCountedNotFollowed()
    {
        using var f = new Fixture();
        string target = Path.Combine(f.Scratch, "elsewhere");
        f.Write(Path.Combine(target, "big.bin"), 5000);
        f.Write("artifacts/cache/own.bin", 10);
        if (!OperatingSystem.IsWindows()) return;   // junctions are Windows links
        Assert.True(TryJunction(f.Full("artifacts/cache/link"), target), "mklink /J failed");

        StorageProfile p = f.Profile(Rules(Rule("artifacts/cache/", Category.Regenerable)));

        StorageItem item = p.Find("repo", "artifacts/cache/")!;
        Assert.Equal(10, item.Bytes);
        Assert.Equal(1, item.Links);
    }

    [Fact]
    public void GitStatisticsAndWorktreesAreReported()
    {
        using var f = new Fixture();

        StorageProfile p = f.Profile(Rules());

        Assert.True(p.Git.CountObjects.ContainsKey("count"));
        Assert.Single(p.Git.Worktrees);
        Assert.True(p.Git.ReflogEntries >= 1);
        Assert.Empty(p.OrphanedWorktreeFolders);
    }

    [Fact]
    public void AFolderUnderClaudeWorktreesThatGitDoesNotListIsOrphaned()
    {
        using var f = new Fixture();
        f.Write(".claude/worktrees/stale/file.txt");

        StorageProfile p = f.Profile(Rules());

        Assert.Equal([".claude/worktrees/stale/"], p.OrphanedWorktreeFolders);
    }
}
