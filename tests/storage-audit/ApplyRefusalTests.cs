using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using StorageAudit;
using Xunit;
using static StorageAudit.Tests.Fixture;

namespace StorageAudit.Tests;

/// <summary>Everything that changed since the plan was reviewed is re-checked, and refused.</summary>
public class ApplyRefusalTests
{
    [Fact]
    public void RefusesATrackedFileAddedSinceThePlan()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        f.Commit("artifacts/bin/a.dll");

        Assert.Contains("tracked", Refused(f.Apply(execute: true)));
        Assert.True(File.Exists(f.Full("artifacts/bin/a.dll")));
    }

    [Fact]
    public void RefusesAPathWrittenWithinTheWindowSinceThePlan()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        File.SetLastWriteTimeUtc(f.Full("artifacts/bin/a.dll"), DateTime.UtcNow);

        Assert.Contains("written within", Refused(f.Apply(execute: true)));
        Assert.True(File.Exists(f.Full("artifacts/bin/a.dll")));
    }

    [Fact]
    public void KeepsAFileWrittenBetweenTheFreshProfileAndTheDelete()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Write("artifacts/bin/b.dll");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));

        ApplyResult r = Assert.Single(f.Apply(execute: true,
            beforeAct: _ => File.SetLastWriteTimeUtc(f.Full("artifacts/bin/b.dll"), DateTime.UtcNow)));

        Assert.Equal("partial", r.Outcome);
        Assert.Contains("written since the profile", r.Detail);
        Assert.False(File.Exists(f.Full("artifacts/bin/a.dll")));
        Assert.True(File.Exists(f.Full("artifacts/bin/b.dll")));
    }

    [Fact]
    public void RefusesAPlanForAnotherRepositoryWhole()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));

        Assert.Throws<InvalidDataException>(() => f.Apply(execute: true, expectedRepo: f.Scratch));
        Assert.True(File.Exists(f.Full("artifacts/bin/a.dll")));
        Assert.False(File.Exists(f.LedgerPath));
    }

    [Fact]
    public void RefusesWhenTheRulesNoLongerGiveTheActionsCategory()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        f.SaveRules(Rules(Rule("artifacts/bin/", Category.Keep)));

        Assert.Contains("now call it keep", Refused(f.Apply(execute: true)));
    }

    [Fact]
    public void RefusesAnItemThatGrewSinceThePlan()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll", 1000);
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        f.Write("artifacts/bin/b.dll", 3 << 20);
        f.Age(10);

        Assert.Contains("grew from", Refused(f.Apply(execute: true)));
        Assert.True(File.Exists(f.Full("artifacts/bin/a.dll")));
    }

    [Theory]
    [InlineData("../outside/")]
    [InlineData("artifacts/../../outside/")]
    [InlineData(@"C:\Windows\")]
    [InlineData("")]
    public void RefusesAPathThatEscapesItsRoot(string path)
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        StoragePlan plan = f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        (plan with { Actions = [plan.Actions[0] with { Path = path }] }).Write(f.PlanPath);

        Assert.Contains("escapes its root", Refused(f.Apply(execute: true)));
    }

    [Fact]
    public void RefusesARootTheRulesDoNotList_SoNothingOutsideTheRepoIsReachableUnlisted()
    {
        using var f = new Fixture();
        string temp = Path.Combine(f.Scratch, "temp");
        f.Write(Path.Combine(temp, "s1", "x.txt"));
        f.Age(10, temp);
        StorageRules rules = Rules(Rule("*/", Category.Regenerable, root: "t")) with { TempRoots = [new TempRoot { Name = "t", Path = temp }] };
        f.Plan(rules);
        f.SaveRules(Rules());

        Assert.Contains("not in the rules", Refused(f.Apply(execute: true)));
        Assert.True(File.Exists(Path.Combine(temp, "s1", "x.txt")));
    }

    [Fact]
    public void ExecutesUnderAListedTempRoot()
    {
        using var f = new Fixture();
        string temp = Path.Combine(f.Scratch, "temp");
        f.Write(Path.Combine(temp, "old-session", "x.txt"), 400);
        f.Write(Path.Combine(temp, "live-session", "x.txt"), 400);
        f.Age(10, temp);
        File.SetLastWriteTimeUtc(Path.Combine(temp, "live-session", "x.txt"), DateTime.UtcNow);
        StorageRules rules = Rules(Rule("*/", Category.Regenerable, root: "t")) with { TempRoots = [new TempRoot { Name = "t", Path = temp }] };
        StoragePlan plan = f.Plan(rules);

        PlanAction action = Assert.Single(plan.Actions);
        Assert.Equal("old-session/", action.Path);
        Assert.Equal("deleted", Assert.Single(f.Apply(execute: true)).Outcome);
        Assert.False(Directory.Exists(Path.Combine(temp, "old-session")));
        Assert.True(Directory.Exists(Path.Combine(temp, "live-session")));
    }

    [Fact]
    public void RefusesAPathThatNoLongerExists()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        Directory.Delete(f.Full("artifacts/bin"), recursive: true);

        Assert.Contains("no longer exists", Refused(f.Apply(execute: true)));
    }

    [Fact]
    public void RefusesAPathHandEditedIntoThePlanThatNoRuleClaims()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Write("src/Program.cs");
        f.Commit("src/Program.cs");
        f.Age(10);
        StoragePlan plan = f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        (plan with { Actions = [plan.Actions[0] with { Path = "src/" }] }).Write(f.PlanPath);

        Assert.Contains("no rule classifies", Refused(f.Apply(execute: true)));
        Assert.True(File.Exists(f.Full("src/Program.cs")));
    }

    [Fact]
    public void RefusesAnItemHoldingALink()
    {
        using var f = new Fixture();
        string target = Path.Combine(f.Scratch, "elsewhere");
        f.Write(Path.Combine(target, "precious.bin"));
        f.Write("artifacts/cache/a.bin");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/cache/", Category.Regenerable)));
        if (!OperatingSystem.IsWindows()) return;   // junctions are Windows links
        Assert.True(TryJunction(f.Full("artifacts/cache/link"), target), "mklink /J failed");
        Directory.SetLastWriteTimeUtc(f.Full("artifacts/cache"), DateTime.UtcNow.AddDays(-10));

        Assert.Contains("link", Refused(f.Apply(execute: true)));
        Assert.True(File.Exists(Path.Combine(target, "precious.bin")));
    }
}
