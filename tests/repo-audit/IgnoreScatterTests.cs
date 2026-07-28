using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using RepoAudit;
using Xunit;

namespace RepoAudit.Tests;

// Regression coverage for Stage 4 (Gather-Scatter) of the gitignore pipeline.
//
// The stage keys compiled regex pairs by a signature over the node's effective pattern
// sets, so sibling directories that inherit the same rules compile once and share. The
// cache is an optimisation only: it decides whether a regex pair is recompiled, never
// whether a node lands in the compiled lookup. An earlier revision stamped the result
// dictionary inside the cache-miss branch, which would have silently dropped every
// cache-hit node — and with it, all filtering for that directory. These tests pin both
// halves of the contract so the two can't drift apart again.
public class IgnoreScatterTests
{
    private static readonly string Root = Path.Combine(Path.GetTempPath(), "ignore-scatter-fixture");

    // ── GatherScatter (internal; reachable via InternalsVisibleTo) ──────────

    [Fact]
    public void IdenticalSignatures_StampEveryNode_AndShareOneRegexPair()
    {
        // Two sibling directories with byte-identical effective pattern sets: the second
        // is a cache hit, and must still appear in the lookup with the shared regexes.
        GitIgnoreNode a = Node("a/", positives: ["*.log", "tmp/"], exceptions: ["keep.log"]);
        GitIgnoreNode b = Node("b/", positives: ["*.log", "tmp/"], exceptions: ["keep.log"]);

        Dictionary<string, CompiledIgnoreState> compiled =
            GitIgnoreCompiler.GatherScatter([a, b], GlobSemantics.Ignore);

        Assert.Equal(2, compiled.Count);
        Assert.True(compiled.ContainsKey("a/"));
        Assert.True(compiled.ContainsKey("b/"));

        // Scatter: one compilation, two references — not two equivalent instances.
        Assert.NotNull(compiled["a/"].Positives);
        Assert.Same(compiled["a/"].Positives,  compiled["b/"].Positives);
        Assert.Same(compiled["a/"].Exceptions, compiled["b/"].Exceptions);
    }

    [Fact]
    public void SignatureIsOrderInsensitive_SoReorderedPatternsStillShare()
    {
        // The signature sorts before joining, so inheritance order can't fragment the cache.
        GitIgnoreNode a = Node("a/", positives: ["*.log", "tmp/"], exceptions: []);
        GitIgnoreNode b = Node("b/", positives: ["tmp/", "*.log"], exceptions: []);

        Dictionary<string, CompiledIgnoreState> compiled =
            GitIgnoreCompiler.GatherScatter([a, b], GlobSemantics.Ignore);

        Assert.Same(compiled["a/"].Positives, compiled["b/"].Positives);
    }

    [Fact]
    public void DifferentSignatures_CompileIndependently()
    {
        // The converse guard: sharing must be driven by the signature, not applied blindly.
        GitIgnoreNode a = Node("a/", positives: ["*.log"], exceptions: []);
        GitIgnoreNode b = Node("b/", positives: ["*.tmp"], exceptions: []);

        Dictionary<string, CompiledIgnoreState> compiled =
            GitIgnoreCompiler.GatherScatter([a, b], GlobSemantics.Ignore);

        Assert.Equal(2, compiled.Count);
        Assert.NotSame(compiled["a/"].Positives, compiled["b/"].Positives);

        Regex positives = compiled["a/"].Positives!;
        Assert.Matches(positives, "a/build.log");
        Assert.DoesNotMatch(positives, "a/build.tmp");
    }

    // ── End-to-end through the public Filter entry point ────────────────────

    [Fact]
    public void SiblingsInheritingOneRootPattern_AreBothFiltered()
    {
        // Root injects '*.log'; the Walk stage propagates it verbatim to both children, so
        // all three nodes share a signature. Every one of them must still filter its files —
        // a node missing from the compiled lookup is dropped from the result entirely.
        CrawlerResult crawl = Crawl(
            ("",   []),
            ("a/", ["build.log", "main.cs"]),
            ("b/", ["build.log", "util.cs"]));

        CrawlerResult filtered = GitIgnoreCompiler.Filter(
            crawl, sentinelNames: [], extraPatterns: ["*.log"]);

        Assert.Equal(["main.cs"], FileNames(filtered, "a/"));
        Assert.Equal(["util.cs"], FileNames(filtered, "b/"));
        Assert.Equal(2, filtered.FileCount);
    }

    // ── Helpers ─────────────────────────────────────────────────────────────

    private static GitIgnoreNode Node(string nodePath, string[] positives, string[] exceptions) =>
        new(new CrawlNode(nodePath, Path.Combine(Root, nodePath.TrimEnd('/')), 1, []))
        {
            EffectivePositives  = positives,
            EffectiveExceptions = exceptions,
        };

    private static CrawlerResult Crawl(params (string NodePath, string[] Files)[] dirs)
    {
        var graph = new Dictionary<string, CrawlNode>(StringComparer.Ordinal);

        foreach ((string nodePath, string[] files) in dirs)
        {
            string trimmed = nodePath.TrimEnd('/');
            string abs     = trimmed.Length == 0 ? Root : Path.Combine(Root, trimmed);
            int    depth   = trimmed.Length == 0 ? 0 : trimmed.Split('/').Length;

            graph[nodePath] = new CrawlNode(nodePath, abs, depth,
                [.. files.Select(f => new CrawlFile(Path.Combine(abs, f), 0L))]);
        }

        return new CrawlerResult(Root, graph, graph.Count, graph.Values.Sum(n => n.Files.Count), []);
    }

    private static string[] FileNames(CrawlerResult crawl, string nodePath) =>
        [.. crawl.Graph[nodePath].Files.Select(f => Path.GetFileName(f.AbsolutePath)).OrderBy(n => n, StringComparer.Ordinal)];
}
