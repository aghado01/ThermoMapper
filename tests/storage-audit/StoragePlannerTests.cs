using System;
using System.IO;
using System.Linq;
using StorageAudit;
using Xunit;
using static StorageAudit.Tests.Fixture;

namespace StorageAudit.Tests;

public class StoragePlannerTests
{
    [Fact]
    public void RegenerablesAreDeleted_ArchivablesArchived_KeepAndAskLeftAlone()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll", 3000);
        f.Write("artifacts/spc/run.json", 2000);
        f.Write("artifacts/kept/x.txt", 1000);
        f.Write("artifacts/odd/x.txt", 500);
        f.Age(10);

        StoragePlan plan = f.Plan(Rules(
            Rule("artifacts/bin/", Category.Regenerable),
            Rule("artifacts/spc/", Category.Archivable),
            Rule("artifacts/kept/", Category.Keep),
            Rule("artifacts/odd/", Category.Ask)));

        Assert.Equal(2, plan.Actions.Count);
        Assert.Equal((PlanActionKind.Delete, "artifacts/bin/", 3000L), (plan.Actions[0].Action, plan.Actions[0].Path, plan.Actions[0].Bytes));
        Assert.Equal((PlanActionKind.Archive, "artifacts/spc/"), (plan.Actions[1].Action, plan.Actions[1].Path));
        Assert.All(plan.Actions, a => Assert.True(a.Enabled));
        HeldItem ask = Assert.Single(plan.Held);
        Assert.StartsWith("ask:", ask.Why);
        Assert.DoesNotContain(plan.Held, h => h.Path == "artifacts/kept/");
    }

    [Fact]
    public void ThePlanFileRoundTrips()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);

        StoragePlan plan = f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        StoragePlan read = StoragePlan.Read(f.PlanPath);

        Assert.Equal(plan.Actions.Single(), read.Actions.Single());
        Assert.Equal(plan.RepoRoot, read.RepoRoot);
        Assert.Equal(plan.RulesSha256, read.RulesSha256);
        Assert.Contains("\"action\": \"delete\"", File.ReadAllText(f.PlanPath));
    }

    [Fact]
    public void ActiveTrackedUnignoredAndProtectedItemsAreHeldWithTheirReason()
    {
        using var f = new Fixture();
        f.Write("artifacts/fresh/a.txt");
        f.Write("artifacts/forced/a.dll");
        f.Commit("artifacts/forced/a.dll");
        f.Write("notes/draft.md");
        f.Write("issues/scratch/cache.bin");
        f.Age(10);
        File.SetLastWriteTimeUtc(f.Full("artifacts/fresh/a.txt"), f.Now.AddHours(-2));

        StoragePlan plan = f.Plan(Rules(
            Rule("artifacts/fresh/", Category.Regenerable),
            Rule("artifacts/forced/", Category.Regenerable),
            Rule("notes/", Category.Regenerable),
            Rule("issues/scratch/", Category.Regenerable)));

        Assert.Empty(plan.Actions);
        Assert.Contains("written within", Why(plan, "artifacts/fresh/"));
        Assert.Contains("tracked", Why(plan, "artifacts/forced/"));
        Assert.Contains("does not ignore", Why(plan, "notes/"));
        Assert.Contains("protected", Why(plan, "issues/scratch/"));
    }

    [Fact]
    public void AnItemHoldingTheRunningBinaryIsHeld()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/StorageAudit/StorageAudit.dll");
        f.Age(10);

        StoragePlan plan = f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)),
            ownBinaryDirectory: f.Full("artifacts/bin/StorageAudit"));

        Assert.Empty(plan.Actions);
        Assert.Contains("running storage-audit binary", Why(plan, "artifacts/bin/"));
    }

    [Fact]
    public void AnItemWhoseContentIsAllCarvedOutProposesNothing()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);

        StoragePlan plan = f.Plan(Rules(Rule("artifacts/", Category.Regenerable), Rule("artifacts/bin/", Category.Keep)));

        Assert.Empty(plan.Actions);
        Assert.Empty(plan.Held);
    }

    [Fact]
    public void GitAdviceProposesGcForLooseObjects_AndNeverAnAction()
    {
        using var f = new Fixture();

        StoragePlan plan = f.Plan(Rules());

        Assert.Contains(plan.GitAdvice, a => a.StartsWith("git gc", StringComparison.Ordinal));
        Assert.Empty(plan.Actions);
    }

    private static string Why(StoragePlan plan, string path) => plan.Held.Single(h => h.Path == path).Why;
}
