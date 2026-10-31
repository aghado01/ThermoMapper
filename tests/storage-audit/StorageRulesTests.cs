using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using StorageAudit;
using Xunit;

namespace StorageAudit.Tests;

public class StorageRulesTests
{
    private static string Json(string rules, string extra = "")
        => "{ " + extra + " \"rules\": [ " + rules + " ] }";

    [Fact]
    public void ParsesCommentsTrailingCommasAndCategoryNames()
    {
        StorageRules rules = StorageRules.Parse("""
            // a comment
            {
              "recentWindowHours": 48,
              "protect": ["issues/"],
              "tempRoots": [ { "name": "t", "path": "C:\\tmp" }, ],
              "rules": [
                { "glob": "artifacts/bin/", "category": "regenerable", "reason": "build", "kind": "build" },
                { "glob": "*/", "root": "t", "category": "archivable", "reason": "temp", "keepNewest": 2 },
              ]
            }
            """);

        Assert.Equal(48, rules.RecentWindowHours);
        Assert.Equal(Category.Regenerable, rules.Rules[0].Category);
        Assert.Equal("t", rules.Rules[1].Root);
        Assert.Equal(2, rules.Rules[1].KeepNewest);
        Assert.Equal(StorageRules.RepoRootName, rules.Rules[0].Root);
    }

    [Theory]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r", "root": "nowhere" }""", "")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r", "keepNewest": 0 }""", "")]
    [InlineData("""{ "glob": "../x/", "category": "keep", "reason": "r" }""", "")]
    [InlineData("""{ "glob": "/", "category": "keep", "reason": "r" }""", "")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": " " }""", "")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r" }, { "glob": "X/", "category": "ask", "reason": "r" }""", "")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r" }""", "\"recentWindowHours\": 0,")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r", "recentWindowHours": 2 }""", "")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r" }""", "\"protect\": [\"C:\\\\Windows\"],")]
    [InlineData("""{ "glob": "x/", "category": "keep", "reason": "r" }""", "\"tempRoots\": [ { \"name\": \"repo\", \"path\": \"C:\\\\t\" } ],")]
    public void RejectsMalformedRules(string rules, string extra)
        => Assert.Throws<InvalidDataException>(() => StorageRules.Parse(Json(rules, extra)));

    [Fact]
    public void RequiresGlobCategoryAndReason()
        => Assert.ThrowsAny<JsonException>(() => StorageRules.Parse(Json("""{ "glob": "x/", "category": "keep" }""")));

    [Fact]
    public void RepoSlugTurnsEveryNonAlphanumericIntoADash()
        => Assert.Equal("D--aghado01-ThermoMapper", StorageRules.RepoSlug(@"D:\aghado01\ThermoMapper"));

    [Fact]
    public void TempRootPathsExpandEnvironmentVariablesAndTheRepoSlug()
    {
        Environment.SetEnvironmentVariable("STORAGE_AUDIT_TEST_ROOT", @"C:\scratch");
        StorageRules rules = new() { TempRoots = [new TempRoot { Name = "t", Path = @"%STORAGE_AUDIT_TEST_ROOT%\claude\{repoSlug}" }] };

        StorageRoot temp = rules.ResolveRoots(@"D:\work\Repo").Single(r => r.Name == "t");

        Assert.Equal(@"C:\scratch\claude\D--work-Repo", temp.FullPath);
        Assert.False(temp.IsRepo);
    }

    [Fact]
    public void TheShippedRulesFileParses()
    {
        // The test binary sits under the repository's artifacts root; the rules file is above it.
        string? path = null;
        for (string? dir = AppContext.BaseDirectory.TrimEnd('\\', '/'); dir is not null && path is null; dir = Path.GetDirectoryName(dir))
        {
            string candidate = Path.Combine(dir, "presets", "storage-audit", "rules.jsonc");
            if (File.Exists(candidate)) path = candidate;
        }
        Assert.True(path is not null, "presets/storage-audit/rules.jsonc not found above the test binary");

        StorageRules rules = StorageRules.Load(path!);
        Assert.Contains(rules.Rules, r => r.Glob == "artifacts/bin/" && r.Category == Category.Regenerable);
        Assert.Contains(rules.Rules, r => r.Glob == "artifacts/chip-*/" && r.Category == Category.Ask);
        Assert.Contains("issues/", rules.Protect);
    }
}
