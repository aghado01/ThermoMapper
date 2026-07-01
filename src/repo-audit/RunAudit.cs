using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace RepoAudit
{
    internal static class Program
    {
        public static int Main(string[] args)
        {
    string repoRoot = Directory.GetCurrentDirectory();
    string? impactNs = null;
    string? attemptDir = null;
    bool validate = false;
    bool useGit = true;
    bool semantic = true;
    bool link = false;

    for (int i = 0; i < args.Length; i++)
    {
        switch (args[i])
        {
            case "--impact" when i + 1 < args.Length:
                impactNs = args[++i];
                break;
            case "--attempt-dir" when i + 1 < args.Length:
                attemptDir = Path.GetFullPath(args[++i]);
                break;
            case "--validate":
                validate = true;
                break;
            case "--no-git":
                useGit = false;
                break;
            case "--semantic":
                semantic = true;
                break;
            case "--no-semantic":
                semantic = false;
                break;
            case "--link":
                link = true;
                break;
            default:
                if (!args[i].StartsWith("--"))
                    repoRoot = Path.GetFullPath(args[i]);
                else
                {
                    Console.Error.WriteLine($"Unknown option: {args[i]}");
                    return 2;
                }
                break;
        }
    }

    string srcRoot = Path.Combine(repoRoot, "src");
    string projectsRoot = Path.Combine(repoRoot, "projects");
    string testsRoot = Path.Combine(repoRoot, "tests");

    if (!Directory.Exists(srcRoot))
    {
        Console.Error.WriteLine($"[fatal] src not found under {repoRoot}");
        return 1;
    }

    string artifactsDir = attemptDir ?? Path.Combine(
        repoRoot, "artifacts", "repo-audit",
        DateTime.Now.ToString("yyyyMMdd_HHmmss"));
    Directory.CreateDirectory(artifactsDir);

    Console.WriteLine("repo-audit");
    Console.WriteLine($"  repo:    {repoRoot}");
    Console.WriteLine($"  attempt: {Path.GetRelativePath(repoRoot, artifactsDir).Replace('\\', '/')}");
    Console.WriteLine();

    // [crawl] — src/ is production source; tests/ is a first-class source root
    // (shared fixtures + per-project test/smoke sources). Both feed file-ownership.
    var crawls = new List<CrawlerResult> { new FileSystemCrawler(srcRoot).Invoke() };
    if (Directory.Exists(testsRoot))
        crawls.Add(new FileSystemCrawler(testsRoot).Invoke());

    int crawlFiles   = crawls.Sum(c => c.FileCount);
    int crawlDirs    = crawls.Sum(c => c.DirectoryCount);
    int crawlSkipped = crawls.Sum(c => c.Skipped.Count);
    Console.WriteLine($"[crawl]    OK    {crawlFiles} files in {crawlDirs} dirs" +
        (crawlSkipped > 0 ? $" ({crawlSkipped} skipped)" : ""));

    // [projects] — csproj live under projects/ (production) and tests/ (test & smoke harnesses)
    var csprojList = new List<string>();
    if (Directory.Exists(projectsRoot))
        csprojList.AddRange(Directory.GetFiles(projectsRoot, "*.csproj", SearchOption.AllDirectories));
    if (Directory.Exists(testsRoot))
        csprojList.AddRange(Directory.GetFiles(testsRoot, "*.csproj", SearchOption.AllDirectories));
    string[] csprojFiles = csprojList.ToArray();
    Console.WriteLine($"[projects] OK    {csprojFiles.Length} csproj files");

    // [syntax]
    var allGlobs = csprojFiles.Select(CsprojParser.Parse).ToList();
    var ownership = FileOwnership.Build(allGlobs, crawls);
    var globsByPath = allGlobs.ToDictionary(g => g.CsprojPath, StringComparer.OrdinalIgnoreCase);
    var analyses = ownership.Projects
        .Where(p => globsByPath.ContainsKey(p.CsprojPath))
        .Select(p => ProjectAnalyzer.Analyze(p, globsByPath[p.CsprojPath]))
        .ToList();

    var xref = CrossReference.Build(analyses);
    if (useGit)
    {
        var renames = GitRenameDetector.TryGetRenames(repoRoot);
        var renameGaps = GitRenameDetector.DetectRenameGaps(ownership, renames, repoRoot);
        if (renameGaps.Count > 0)
            xref = xref with { Violations = xref.Violations.Concat(renameGaps).ToList().AsReadOnly() };
    }
    var driftViolations = CrossReference.DetectNamespaceDrift(analyses, srcRoot);
    if (driftViolations.Count > 0)
        xref = xref with { Violations = xref.Violations.Concat(driftViolations).ToList().AsReadOnly() };
    var orphanViolations = CrossReference.DetectPartialOrphans(analyses);
    if (orphanViolations.Count > 0)
        xref = xref with { Violations = xref.Violations.Concat(orphanViolations).ToList().AsReadOnly() };
    var solutionGapViolations = CrossReference.DetectSolutionGaps(repoRoot, csprojFiles);
    if (solutionGapViolations.Count > 0)
        xref = xref with { Violations = xref.Violations.Concat(solutionGapViolations).ToList().AsReadOnly() };

    Console.WriteLine($"[syntax]   OK    {analyses.Count} analyses, " +
        $"{ownership.OrphanedFiles.Count} orphaned, " +
        $"{ownership.BrokenIncludes.Count} broken includes, " +
        $"{xref.Violations.Count} violation{(xref.Violations.Count == 1 ? "" : "s")}");

    if (impactNs != null)
    {
        Console.WriteLine();
        RunImpact(impactNs, analyses, xref, repoRoot);
        return 0;
    }

    // [semantic]
    IReadOnlyList<ProjectDiagnostics>? semanticResults = null;
    string? semanticFailurePath = null;
    bool semanticCrashed = false;

    if (semantic)
    {
        // Semantic in-process compilation is a fidelity check over the PRODUCTION source
        // graph. Test/smoke projects under tests/ pull in xunit + the test SDK and lean on
        // transitive project references the in-process compiler does not resolve, so compiling
        // them here yields only false CS0246/CS8805 noise — their real build is authoritatively
        // validated by `dotnet test` / CI. They stay fully covered by the syntax/ownership/
        // cross-reference analysis above; only the in-process compile skips them.
        string testsPrefix = Path.GetFullPath(testsRoot).TrimEnd('/', '\\') + Path.DirectorySeparatorChar;
        var productionAnalyses = analyses
            .Where(a => !Path.GetFullPath(a.CsprojPath).StartsWith(testsPrefix, StringComparison.OrdinalIgnoreCase))
            .ToList();
        var sorted = CrossReference.TopologicalSort(productionAnalyses);
        Console.WriteLine($"[semantic] running {sorted.Count} projects (tests/ excluded — validated by dotnet test)...");
        try
        {
            // Only surface per-project lines for non-clean compiles
            void Progress(string line)
            {
                if (!line.TrimEnd().EndsWith("clean", StringComparison.Ordinal))
                    Console.WriteLine(line);
            }

            semanticResults = InProcessCompiler.Compile(sorted, progress: Progress);
            int totalErrors   = semanticResults.Sum(r => r.Errors.Count);
            int totalWarnings = semanticResults.Sum(r => r.Warnings.Count);
            int cleanProjects = semanticResults.Count(r => r.Errors.Count == 0 && r.Warnings.Count == 0);

            Console.WriteLine($"[semantic] OK    {cleanProjects}/{sorted.Count} clean, " +
                $"{totalErrors} error{(totalErrors == 1 ? "" : "s")}, " +
                $"{totalWarnings} warning{(totalWarnings == 1 ? "" : "s")}");

            ArtifactWriter.WriteSemantic(artifactsDir, repoRoot, semanticResults);
        }
        catch (Exception ex)
        {
            semanticCrashed = true;
            semanticFailurePath = ArtifactWriter.WriteSemanticFailure(artifactsDir, repoRoot, ex);
            Console.WriteLine($"[semantic] FAIL  {ex.GetType().Name}: {ex.Message}");
            Console.WriteLine($"                 trace: {Path.GetRelativePath(repoRoot, semanticFailurePath).Replace('\\', '/')}");
        }
    }
    else
    {
        Console.WriteLine("[semantic] SKIP  --no-semantic");
    }

    // [refine] — opportunistic: drop StaleProjectReference violations that
    // semantic analysis confirms are actually used (var-inferred return types
    // and other transitive uses that don't surface as syntactic `using` directives).
    if (semanticResults != null)
    {
        var usedByConsumer = semanticResults.ToDictionary(
            r => r.AssemblyName,
            r => r.UsedProjectReferences,
            StringComparer.Ordinal);

        int before = xref.Violations.Count;
        var filtered = xref.Violations.Where(v =>
        {
            if (v.Kind != ViolationKind.StaleProjectReference) return true;
            if (v.TargetAssembly is null) return true;
            return !(usedByConsumer.TryGetValue(v.Project, out var used)
                     && used.Contains(v.TargetAssembly));
        }).ToList().AsReadOnly();

        int suppressed = before - filtered.Count;
        if (suppressed > 0)
        {
            xref = xref with { Violations = filtered };
            Console.WriteLine($"[refine]   OK    {suppressed} stale-ref violation{(suppressed == 1 ? "" : "s")} suppressed by semantic-confirmed usage");
        }
    }

    // [write]
    var written = new List<string>();
    try
    {
        var (jsonPath, healthPath) = ArtifactWriter.Write(artifactsDir, repoRoot, analyses, ownership, xref);
        string violationsPath = ArtifactWriter.WriteViolations(
            artifactsDir, repoRoot, analyses, xref, semanticResults);
        written.Add(jsonPath);
        written.Add(healthPath);
        written.Add(violationsPath);
        if (semantic)
            written.Add(Path.Combine(artifactsDir, "semantic-analysis.json"));
        Console.WriteLine("[write]    OK");
    }
    catch (Exception ex)
    {
        Console.Error.WriteLine($"[write]    FAIL  {ex.GetType().Name}: {ex.Message}");
        return 1;
    }

    foreach (string p in written)
        Console.WriteLine($"  {Path.GetRelativePath(repoRoot, p).Replace('\\', '/')}");

    // [link] (optional, summary only)
    if (link)
    {
        var linkage = SymbolLinker.BuildLinkage(analyses);
        var unlinked = SymbolLinker.FindEffectivelyUnlinkedProjects(linkage);
        Console.WriteLine($"[link]     OK    {linkage.Count} in graph, {unlinked.Count} effectively unlinked");
    }

    // Final summary
    int semErrs  = semanticResults?.Sum(r => r.Errors.Count)   ?? 0;
    int semWarns = semanticResults?.Sum(r => r.Warnings.Count) ?? 0;
    Console.WriteLine();
    Console.WriteLine($"result: {xref.Violations.Count} cross-ref violation{(xref.Violations.Count == 1 ? "" : "s")}, " +
        $"{semErrs} error{(semErrs == 1 ? "" : "s")}, " +
        $"{semWarns} warning{(semWarns == 1 ? "" : "s")}");

    if (validate)
    {
        if (semanticCrashed)
        {
            Console.Error.WriteLine("[validate] FAIL  semantic analysis crashed");
            return 1;
        }
        if (semErrs > 0)
        {
            Console.Error.WriteLine($"[validate] FAIL  {semErrs} compilation error{(semErrs == 1 ? "" : "s")}");
            return 1;
        }
        if (xref.Violations.Count > 0)
        {
            Console.Error.WriteLine($"[validate] FAIL  {xref.Violations.Count} violation{(xref.Violations.Count == 1 ? "" : "s")}");
            return 1;
        }
        Console.WriteLine("[validate] OK");
    }

    return 0;
}

// ── --impact helper ────────────────────────────────────────────────────────────
private static void RunImpact(string ns, IReadOnlyList<ProjectAnalysis> analyses,
               CrossReferenceIndex xref, string repoRoot)
{
    Console.WriteLine($"Impact analysis for: {ns}");

    var declarers = analyses
        .Where(a => a.DeclaredNamespaces.Contains(ns))
        .OrderBy(a => a.AssemblyName)
        .ToList();

    if (declarers.Count == 0)
    {
        Console.WriteLine($"Namespace {ns} is not declared in any project.");
        return;
    }

    Console.WriteLine("Declared in:");
    foreach (var a in declarers)
    {
        var files = a.FileAnalyses
            .Where(f => f.DeclaredNamespaces.Contains(ns))
            .Select(f => Path.GetRelativePath(repoRoot, f.FilePath).Replace('\\', '/'))
            .OrderBy(f => f)
            .ToList();
        Console.WriteLine($"  {a.AssemblyName} ({files.Count} file{(files.Count != 1 ? "s" : "")})");
        foreach (string f in files) Console.WriteLine($"    {f}");
    }

    if (!xref.ImportersOf.TryGetValue(ns, out var importerNames) || importerNames.Count == 0)
    {
        Console.WriteLine("Imported by: none");
        return;
    }

    var importerLookup = analyses.ToDictionary(a => a.AssemblyName, StringComparer.Ordinal);
    Console.WriteLine("Imported by:");
    foreach (string name in importerNames.OrderBy(n => n))
    {
        if (!importerLookup.TryGetValue(name, out var a)) continue;
        var files = a.FileAnalyses
            .Where(f => f.ImportedNamespaces.Contains(ns))
            .Select(f => Path.GetRelativePath(repoRoot, f.FilePath).Replace('\\', '/'))
            .OrderBy(f => f)
            .ToList();
        Console.WriteLine($"  {name} ({files.Count} file{(files.Count != 1 ? "s" : "")})");
        foreach (string f in files) Console.WriteLine($"    {f}");
    }
}
}
}
