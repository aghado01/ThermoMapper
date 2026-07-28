using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace RepoAudit;

// ── Public surface ─────────────────────────────────────────────────────────

/// <summary>
/// Controls how compiled glob patterns are interpreted at filter time.
/// </summary>
public enum GlobSemantics
{
    /// <summary>Matching paths are excluded — standard gitignore behaviour.</summary>
    Ignore,
    /// <summary>Matching paths are kept; all non-matching paths are excluded.</summary>
    Include
}

/// <summary>Immutable compiled ignore state for a single directory node.</summary>
public sealed record CompiledIgnoreState(
    string        NodePath,
    Regex?        Positives,
    Regex?        Exceptions,
    GlobSemantics Semantics);

/// <summary>
/// Gitignore-style filter compiler.
///
/// Runs a five-stage pipeline — Normalize → Coalesce → Walk → Reduce → Gather-Scatter —
/// against the nodes of a <see cref="CrawlerResult"/>, then returns a new
/// <see cref="CrawlerResult"/> with ignored files and pruned directory branches removed.
///
/// Pattern compilation (<c>TranslateGlob</c>, <c>CompileGlobs</c>, <c>GlobSubsumes</c>)
/// is fully delegated to <see cref="GlobCompiler"/>; this class owns only the
/// sentinel scan, inheritance walk, subsumption reduce, and prune passes.
/// </summary>
public static class GitIgnoreCompiler
{
    /// <summary>Sentinel file names scanned in each directory by default.</summary>
    public static readonly string[] DefaultSentinels = [".gitignore", ".snapignore"];

    /// <summary>Root-level glob patterns applied by default.</summary>
    public static readonly string[] DefaultPatterns  = [".git/", "node_modules/"];

    private static readonly Regex SepNorm = new(@"[/\\]+", RegexOptions.Compiled);

    // ── Entry point ───────────────────────────────────────────────────────

    /// <summary>
    /// Compiles ignore rules and returns a filtered <see cref="CrawlerResult"/>.
    /// </summary>
    /// <param name="crawl">Crawler output to filter.</param>
    /// <param name="sentinelNames">
    /// Names of ignore-sentinel files to discover in each directory.
    /// <see langword="null"/> → use <see cref="DefaultSentinels"/>.
    /// Empty array → skip sentinel discovery entirely.
    /// </param>
    /// <param name="extraPatterns">
    /// Additional root-level glob patterns injected before sentinel patterns.
    /// <see langword="null"/> → use <see cref="DefaultPatterns"/>.
    /// Empty array → suppress default root patterns.
    /// </param>
    /// <param name="semantics">
    /// <see cref="GlobSemantics.Ignore"/>: matching paths are excluded (default).<br/>
    /// <see cref="GlobSemantics.Include"/>: matching paths are kept; all others excluded.
    /// Both modes run the full five-stage pipeline.
    /// </param>
    public static CrawlerResult Filter(
        CrawlerResult crawl,
        string[]?     sentinelNames = null,
        string[]?     extraPatterns = null,
        GlobSemantics semantics     = GlobSemantics.Ignore)
    {
        sentinelNames ??= DefaultSentinels;
        extraPatterns ??= DefaultPatterns;

        // Build mutable working nodes in BFS order (root first, then children by depth)
        List<GitIgnoreNode> nodes = [.. crawl.Graph.Values
            .OrderBy(n => n.NodeDepth)
            .ThenBy(n => n.NodePath, StringComparer.Ordinal)
            .Select(n => new GitIgnoreNode(n))];

        var lookup = nodes.ToDictionary(n => n.NodePath, StringComparer.Ordinal);

        // Inject extra patterns as a virtual IgnoreFiles entry on the root node,
        // positioned before any sentinel patterns so sentinels can override them.
        if (extraPatterns.Length > 0 && lookup.TryGetValue("", out GitIgnoreNode? rootNode))
            rootNode.IgnoreFiles.Insert(0, new IgnoreFileEntry("<patterns>", extraPatterns));

        // Sentinel scan — populates IgnoreFiles and strips sentinels from Files
        if (sentinelNames.Length > 0)
            ScanSentinels(nodes, sentinelNames);

        // Stage 0: Normalize — separator collapse, degenerate rejection
        Normalize(nodes);

        // Stage 1: Coalesce — per-node source merge, annihilation, anchor-prefix
        Coalesce(nodes);

        // Stage 2: Walk — BFS inheritance with depth-annotated dictionaries
        Walk(nodes, lookup);

        // Stage 3: Reduce — depth precedence, subsumption heuristic
        Reduce(nodes);

        // Stage 4: Gather-Scatter — signature-keyed regex compilation + scatter
        Dictionary<string, CompiledIgnoreState> compiledLookup = GatherScatter(nodes, semantics);

        // Post-compile: prune ignored directory branches.
        // Skipped for Include mode — file-targeted patterns (e.g. *.cs) would
        // prune parent directories before their children can be evaluated.
        if (semantics == GlobSemantics.Ignore)
            Prune(nodes, compiledLookup);

        return BuildResult(crawl, nodes, compiledLookup);
    }

    /// <summary>
    /// Returns <see langword="true"/> when <paramref name="relativePath"/> should be
    /// excluded from results according to <paramref name="state"/>.
    /// </summary>
    public static bool TestPath(string relativePath, CompiledIgnoreState state)
    {
        if (state.Semantics == GlobSemantics.Include)
        {
            // match = KEEP   →   non-match = excluded
            if (state.Positives is null)                              return true;   // no keep patterns → exclude all
            if (!state.Positives.IsMatch(relativePath))               return true;   // unmatched → excluded
            if (state.Exceptions?.IsMatch(relativePath) is true)      return true;   // exception undoes keep
            return false;
        }
        else // Ignore
        {
            // match = EXCLUDED   →   non-match = kept
            if (state.Positives is null)                              return false;  // no ignore rules → keep all
            if (!state.Positives.IsMatch(relativePath))               return false;  // unmatched → kept
            if (state.Exceptions?.IsMatch(relativePath) is true)      return false;  // exception rescues
            return true;
        }
    }

    // ── Sentinel scan ─────────────────────────────────────────────────────

    private static void ScanSentinels(List<GitIgnoreNode> nodes, string[] sentinelNames)
    {
        var sentinelSet = new HashSet<string>(sentinelNames, StringComparer.OrdinalIgnoreCase);

        foreach (GitIgnoreNode node in nodes)
        {
            var remaining = new List<CrawlFile>(node.Files.Count);

            foreach (CrawlFile f in node.Files)
            {
                string fname = Path.GetFileName(f.AbsolutePath);
                if (!sentinelSet.Contains(fname))
                {
                    remaining.Add(f);
                    continue;
                }

                // Sentinel discovered — read patterns; consumed from Files regardless of outcome
                try
                {
                    string[] globs = File.ReadAllLines(f.AbsolutePath)
                        .Select(l => l.TrimEnd())
                        .Where(l => l.Length > 0 && l[0] != '#')
                        .ToArray();
                    node.IgnoreFiles.Add(new IgnoreFileEntry(fname, globs));
                }
                catch { /* non-fatal — sentinel is still consumed */ }
            }

            node.Files = remaining;
        }
    }

    // ── Stage 0: Normalize ────────────────────────────────────────────────

    private static void Normalize(List<GitIgnoreNode> nodes)
    {
        foreach (GitIgnoreNode node in nodes)
        foreach (IgnoreFileEntry entry in node.IgnoreFiles)
        {
            entry.Globs = [.. entry.Globs
                .Select(NormalizeGlob)
                .OfType<string>()];
        }
    }

    /// <returns>Normalised pattern, or <see langword="null"/> for degenerates/comments.</returns>
    private static string? NormalizeGlob(string raw)
    {
        string p = raw.Trim();
        if (string.IsNullOrWhiteSpace(p)) return null;
        if (p[0] == '#') return null;                         // discard comments

        string prefix = "";
        if (p[0] == '!')
        {
            prefix = "!";
            p = p[1..];
        }

        p = SepNorm.Replace(p, "/");                          // collapse /\ runs to /

        if (string.IsNullOrWhiteSpace(p) || p == "/")
        {
            Console.Error.WriteLine($"GitIgnoreCompiler: discarding degenerate pattern '{raw}'");
            return null;
        }

        return prefix + p;
    }

    // ── Stage 1: Coalesce ─────────────────────────────────────────────────

    private static void Coalesce(List<GitIgnoreNode> nodes)
    {
        foreach (GitIgnoreNode node in nodes)
        {
            var positives  = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            var exceptions = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

            foreach (IgnoreFileEntry entry in node.IgnoreFiles)
            foreach (string glob in entry.Globs)
            {
                string t = glob.Trim();
                if (string.IsNullOrWhiteSpace(t) || t[0] == '#') continue;

                if (t[0] == '!') exceptions.Add(t[1..]);
                else             positives.Add(t);
            }

            // Exact-match annihilation: a pattern and its own negation cancel out
            var annihilated = exceptions.Where(positives.Contains).ToList();
            foreach (string a in annihilated) { positives.Remove(a); exceptions.Remove(a); }

            // Anchor-prefix: anchored patterns declared in a non-root node are
            // prepended with the node's path so they apply only within that subtree
            if (!string.IsNullOrEmpty(node.NodePath))
            {
                node.LocalPositives  = [.. positives .Select(g => AnchorPrefix(g, node.NodePath))];
                node.LocalExceptions = [.. exceptions.Select(g => AnchorPrefix(g, node.NodePath))];
            }
            else
            {
                node.LocalPositives  = [.. positives];
                node.LocalExceptions = [.. exceptions];
            }
        }
    }

    private static string AnchorPrefix(string glob, string nodePath)
    {
        if (string.IsNullOrEmpty(glob)) return glob;
        bool hadLeadingSlash = glob[0] == '/';
        if (hadLeadingSlash) glob = glob[1..];
        bool isAnchored = hadLeadingSlash || glob.TrimEnd('/').Contains('/');
        return isAnchored ? nodePath + glob : glob;
    }

    // ── Stage 2: Walk ─────────────────────────────────────────────────────

    private static void Walk(List<GitIgnoreNode> nodes, Dictionary<string, GitIgnoreNode> lookup)
    {
        // nodes is already in BFS (depth-ascending) order
        foreach (GitIgnoreNode node in nodes)
        {
            string? parentPath = GetParentPath(node.NodePath);

            if (parentPath is not null && lookup.TryGetValue(parentPath, out GitIgnoreNode? parent))
            {
                // Deep-clone parent's active sets so each node owns its own copy
                node.ActiveIgnores    = new Dictionary<string, int>(parent.ActiveIgnores,    StringComparer.OrdinalIgnoreCase);
                node.ActiveExceptions = new Dictionary<string, int>(parent.ActiveExceptions, StringComparer.OrdinalIgnoreCase);
            }

            int depth = node.NodeDepth;

            foreach (string glob in node.LocalPositives)
                node.ActiveIgnores[glob] = depth;

            foreach (string glob in node.LocalExceptions)
            {
                node.ActiveIgnores.Remove(glob);           // cross-depth annihilation
                node.ActiveExceptions[glob] = depth;
            }
        }
    }

    // ── Stage 3: Reduce ───────────────────────────────────────────────────

    private static void Reduce(List<GitIgnoreNode> nodes)
    {
        foreach (GitIgnoreNode node in nodes)
        {
            var survivingExceptions = new List<string>();

            foreach ((string exGlob, int exDepth) in node.ActiveExceptions)
            {
                // An exception is dominated (and therefore dropped) when a deeper
                // positive pattern subsumes it — the more-specific rule wins
                bool dominated = node.ActiveIgnores.Any(ig =>
                    ig.Value > exDepth && GlobCompiler.GlobSubsumes(ig.Key, exGlob));

                if (!dominated) survivingExceptions.Add(exGlob);
            }

            node.EffectivePositives  = [.. node.ActiveIgnores.Keys];
            node.EffectiveExceptions = [.. survivingExceptions];
        }
    }

    // ── Stage 4: Gather-Scatter ───────────────────────────────────────────

    // internal rather than private so the scatter contract — every node stamped,
    // identical signatures sharing one regex pair — is directly testable
    internal static Dictionary<string, CompiledIgnoreState> GatherScatter(
        List<GitIgnoreNode> nodes,
        GlobSemantics semantics)
    {
        // Signature-keyed regex cache: nodes sharing identical effective pattern sets
        // reuse the same compiled Regex objects (scatter)
        var cache  = new Dictionary<string, (Regex? Pos, Regex? Ex)>(StringComparer.Ordinal);
        var result = new Dictionary<string, CompiledIgnoreState>(StringComparer.Ordinal);

        foreach (GitIgnoreNode node in nodes)
        {
            string sig = BuildSignature(node.EffectivePositives, node.EffectiveExceptions);

            if (!cache.TryGetValue(sig, out (Regex? Pos, Regex? Ex) compiled))
            {
                compiled = (
                        GlobCompiler.CompileGlobs(node.EffectivePositives),
                        GlobCompiler.CompileGlobs(node.EffectiveExceptions));
                cache[sig] = compiled;
            }

            // Every node is stamped, cache hit or miss — the cache decides whether the
            // regex pair is recompiled, never whether the node appears in the result
            result[node.NodePath] = new CompiledIgnoreState(
                node.NodePath, compiled.Pos, compiled.Ex, semantics);
        }

        return result;
    }

    private static string BuildSignature(string[] pos, string[] ex) =>
        string.Join("|", pos.OrderBy(s => s, StringComparer.Ordinal))
        + "||" +
        string.Join("|", ex.OrderBy(s => s, StringComparer.Ordinal));

    // ── Post-compile: Prune ───────────────────────────────────────────────

    private static void Prune(
        List<GitIgnoreNode>                     nodes,
        Dictionary<string, CompiledIgnoreState> compiledLookup)
    {
        var pruned = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

        foreach (GitIgnoreNode node in nodes)   // already depth-sorted
        {
            if (node.NodeDepth == 0) continue;

            // Propagate: if any ancestor was pruned, so is this node
            if (HasPrunedAncestor(node.NodePath, pruned))
            {
                pruned.Add(node.NodePath);
                continue;
            }

            // Test the directory path against its direct parent's compiled state
            string? parentPath = GetParentPath(node.NodePath);
            if (parentPath is not null
                && compiledLookup.TryGetValue(parentPath, out CompiledIgnoreState? parentState)
                && TestPath(node.NodePath, parentState))
            {
                pruned.Add(node.NodePath);
            }
        }

        foreach (string p in pruned)
            compiledLookup.Remove(p);
    }

    private static bool HasPrunedAncestor(string nodePath, HashSet<string> pruned)
    {
        string? anc = GetParentPath(nodePath);
        while (anc is not null)
        {
            if (pruned.Contains(anc)) return true;
            anc = GetParentPath(anc);
        }
        return false;
    }

    // ── Result assembly ───────────────────────────────────────────────────

    private static CrawlerResult BuildResult(
        CrawlerResult                           crawl,
        List<GitIgnoreNode>                     nodes,
        Dictionary<string, CompiledIgnoreState> compiledLookup)
    {
        string rootPath = crawl.RootPath.TrimEnd('/').TrimEnd('\\');
        var newGraph = new Dictionary<string, CrawlNode>(StringComparer.Ordinal);

        foreach (GitIgnoreNode node in nodes)
        {
            if (!compiledLookup.TryGetValue(node.NodePath, out CompiledIgnoreState? state))
                continue;   // directory was pruned

            // Filter individual files; node.Files already has sentinels removed
            var keepFiles = node.Files
                .Where(f =>
                {
                    string rel = Path.GetRelativePath(rootPath, f.AbsolutePath).Replace('\\', '/');
                    return !TestPath(rel, state);
                })
                .ToList();

            newGraph[node.NodePath] = new CrawlNode(
                node.NodePath, node.AbsolutePath, node.NodeDepth, keepFiles);
        }

        int dirCount  = newGraph.Count;
        int fileCount = newGraph.Values.Sum(n => n.Files.Count);

        return new CrawlerResult(crawl.RootPath, newGraph, dirCount, fileCount, crawl.Skipped);
    }

    // ── Helpers ───────────────────────────────────────────────────────────

    /// <returns>
    /// Parent node path, or <c>""</c> for a first-level segment, or
    /// <see langword="null"/> for the root node itself.
    /// </returns>
    private static string? GetParentPath(string nodePath)
    {
        if (string.IsNullOrEmpty(nodePath)) return null;
        string trimmed = nodePath.TrimEnd('/');
        int lastSlash  = trimmed.LastIndexOf('/');
        return lastSlash < 0 ? "" : trimmed[..(lastSlash + 1)];
    }
}

// ── Internal working types ─────────────────────────────────────────────────

internal sealed class IgnoreFileEntry
{
    public string   Source { get; }
    public string[] Globs  { get; set; }    // mutable: Stage 0 replaces the array in-place

    public IgnoreFileEntry(string source, string[] globs)
    {
        Source = source;
        Globs  = globs;
    }
}

internal sealed class GitIgnoreNode
{
    // ── Identity (immutable) ──────────────────────────────────────────────
    public string NodePath     { get; }
    public string AbsolutePath { get; }
    public int    NodeDepth    { get; }

    // ── Sentinel scan output ──────────────────────────────────────────────
    public List<CrawlFile>       Files      { get; set; }   // sentinels removed during scan
    public List<IgnoreFileEntry> IgnoreFiles { get; } = [];

    // ── Stage 1 (Coalesce) ────────────────────────────────────────────────
    public List<string> LocalPositives  { get; set; } = [];
    public List<string> LocalExceptions { get; set; } = [];

    // ── Stage 2 (Walk) — depth-annotated ─────────────────────────────────
    public Dictionary<string, int> ActiveIgnores    { get; set; } = new(StringComparer.OrdinalIgnoreCase);
    public Dictionary<string, int> ActiveExceptions { get; set; } = new(StringComparer.OrdinalIgnoreCase);

    // ── Stage 3 (Reduce) ──────────────────────────────────────────────────
    public string[] EffectivePositives  { get; set; } = [];
    public string[] EffectiveExceptions { get; set; } = [];

    public GitIgnoreNode(CrawlNode node)
    {
        NodePath     = node.NodePath;
        AbsolutePath = node.AbsolutePath;
        NodeDepth    = node.NodeDepth;
        Files        = new List<CrawlFile>(node.Files);
    }
}
