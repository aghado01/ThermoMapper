using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using StorageAudit;
using Xunit;
using static StorageAudit.Tests.Fixture;

namespace StorageAudit.Tests;

public class PlanApplierTests
{
    [Fact]
    public void ADryRunChangesNothing_AndLedgersWhatItWouldDo()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll", 2000);
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));

        ApplyResult r = Assert.Single(f.Apply());

        Assert.Equal("would-delete", r.Outcome);
        Assert.Equal(2000, r.Bytes);
        Assert.True(File.Exists(f.Full("artifacts/bin/a.dll")));
        string line = Assert.Single(File.ReadAllLines(f.LedgerPath));
        Assert.Equal("dry-run", JsonDocument.Parse(line).RootElement.GetProperty("mode").GetString());
    }

    [Fact]
    public void ExecuteDeletesTheItem_KeepsWhatADeeperRuleCarvedOut_AndPrunesEmptyFolders()
    {
        using var f = new Fixture();
        f.Write("artifacts/runs/old/deep/log.txt", 300);
        f.Write("artifacts/runs/old/top.txt", 200);
        f.Write("artifacts/runs/keep/table.tsv", 100);
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/runs/", Category.Regenerable), Rule("artifacts/runs/keep/", Category.Keep)));

        ApplyResult r = Assert.Single(f.Apply(execute: true));

        Assert.Equal("deleted", r.Outcome);
        Assert.Equal(500, r.Bytes);
        Assert.False(Directory.Exists(f.Full("artifacts/runs/old")));
        Assert.True(File.Exists(f.Full("artifacts/runs/keep/table.tsv")));
        Assert.Contains("\"outcome\":\"deleted\"", File.ReadAllText(f.LedgerPath));
    }

    [Fact]
    public void ADisabledActionIsSkipped()
    {
        using var f = new Fixture();
        f.Write("artifacts/bin/a.dll");
        f.Age(10);
        StoragePlan plan = f.Plan(Rules(Rule("artifacts/bin/", Category.Regenerable)));
        (plan with { Actions = [plan.Actions[0] with { Enabled = false }] }).Write(f.PlanPath);

        Assert.Equal("skipped", Assert.Single(f.Apply(execute: true)).Outcome);
        Assert.True(File.Exists(f.Full("artifacts/bin/a.dll")));
    }
}
