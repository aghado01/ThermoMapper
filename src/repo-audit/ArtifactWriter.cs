using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;

namespace RepoAudit;

public static class ArtifactWriter
{
    public static (string jsonPath, string healthPath) Write(
        string artifactsDir,
        string repoRoot,
        IReadOnlyList<ProjectAnalysis> analyses,
        FileOwnershipResult ownership,
        CrossReferenceIndex xref)
    {
        Directory.CreateDirectory(artifactsDir);
        string jsonPath   = Path.Combine(artifactsDir, "project-structure.json");
        string healthPath = Path.Combine(artifactsDir, "project-health.md");

        WriteJson(jsonPath, repoRoot, analyses, ownership, xref);
        WriteHealth(healthPath, repoRoot, analyses, ownership, xref);

        return (jsonPath, healthPath);
    }

    // ── JSON ─────────────────────────────────────────────────────────────

    private static void WriteJson(
        string path, string repoRoot,
        IReadOnlyList<ProjectAnalysis> analyses,
        FileOwnershipResult ownership,
        CrossReferenceIndex xref)
    {
        int totalTypes   = analyses.Sum(a => a.Types.Count);
        int partialCount = analyses.Sum(a => a.Types.Count(t => t.IsPartial));

        var doc = new JsonObject
        {
            ["generatedAt"]  = DateTime.UtcNow.ToString("O"),
            ["summary"] = new JsonObject
            {
                ["projectCount"]       = analyses.Count,
                ["typeCount"]          = totalTypes,
                ["partialFamilyCount"] = partialCount,
                ["orphanedFileCount"]  = ownership.OrphanedFiles.Count,
                ["brokenIncludeCount"] = ownership.BrokenIncludes.Count,
                ["violationCount"]     = xref.Violations.Count
            },
            ["projects"]      = BuildProjectsNode(analyses, repoRoot),
            ["orphanedFiles"] = new JsonArray(ownership.OrphanedFiles.OrderBy(f => f)
                                    .Select(f => (JsonNode)Rel(f, repoRoot)).ToArray()),
            ["brokenIncludes"] = new JsonArray(ownership.BrokenIncludes
                                    .Select(b => (JsonNode)new JsonObject
                                    {
                                        ["assembly"]    = b.AssemblyName,
                                        ["missingPath"] = Rel(b.MissingAbsPath, repoRoot)
                                    }).ToArray())
        };

        File.WriteAllText(path, doc.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
    }

    private static JsonObject BuildProjectsNode(IReadOnlyList<ProjectAnalysis> analyses, string repoRoot)
    {
        var root = new JsonObject();
        foreach (ProjectAnalysis a in analyses.OrderBy(a => a.AssemblyName))
        {
            var typesByNs = new JsonObject();
            foreach (IGrouping<string, TypeEntry> nsGroup in a.Types
                .GroupBy(t => t.ContainingNamespace)
                .OrderBy(g => g.Key))
            {
                typesByNs[nsGroup.Key] = new JsonArray(nsGroup.OrderBy(t => t.Name).Select(t =>
                    (JsonNode)new JsonObject
                    {
                        ["name"]         = t.Name,
                        ["kind"]         = t.Kind.ToString(),
                        ["modifiers"]    = new JsonArray(t.Modifiers
                                               .Where(m => m is not "public" and not "internal")
                                               .Select(m => (JsonNode)JsonValue.Create(m)!).ToArray()),
                        ["isPartial"]    = t.IsPartial,
                        ["files"]        = new JsonArray(t.ContributingFiles
                                               .Select(f => (JsonNode)Rel(f, repoRoot)).ToArray()),
                        ["baseTypes"]    = new JsonArray(t.BaseTypes
                                               .Select(b => (JsonNode)JsonValue.Create(b)!).ToArray()),
                        ["genericArity"] = t.GenericArity
                    }).ToArray());
            }

            root[a.AssemblyName] = new JsonObject
            {
                ["csprojPath"]         = Rel(a.CsprojPath, repoRoot),
                ["declaredNamespaces"] = JsonArr(a.DeclaredNamespaces.OrderBy(n => n)),
                ["importedNamespaces"] = JsonArr(a.ImportedNamespaces.OrderBy(n => n)),
                ["projectReferences"]  = JsonArr(a.ProjectReferences),
                ["packageReferences"]  = new JsonArray(a.PackageReferences.OrderBy(p => p.Name)
                                             .Select(p => (JsonNode)new JsonObject
                                             {
                                                 ["name"]    = p.Name,
                                                 ["version"] = p.Version
                                             }).ToArray()),
                ["compilerSettings"]   = new JsonObject
                {
                    ["langVersion"]      = a.CompilerSettings.LangVersion,
                    ["nullable"]         = a.CompilerSettings.Nullable,
                    ["allowUnsafeBlocks"] = a.CompilerSettings.AllowUnsafeBlocks
                },
                ["namespaces"]         = typesByNs
            };
        }
        return root;
    }

    // ── Structured repair hints (consumed by WriteViolations) ─────────────

    private static JsonObject? BuildHint(
        Violation v,
        Dictionary<string, ProjectAnalysis> byName,
        string repoRoot) => v.Kind switch
    {
        ViolationKind.GhostDependency        => GhostHint(v, byName),
        ViolationKind.StaleProjectReference  => StaleHint(v, byName),
        ViolationKind.RenameGap              => RenameHint(v),
        ViolationKind.PackageVersionConflict => PkgConflictHint(v),
        _                                    => null
    };

    // Detail format: "Imports 'ns' (declared by 'asm'[; also declared by ...]) with no <ProjectReference>"
    // The chosen declaring assembly is carried structurally in Violation.TargetAssembly — only
    // the namespace is recovered from Detail, where the leading "Imports '...'" is unambiguous.
    private static JsonObject? GhostHint(
        Violation v,
        Dictionary<string, ProjectAnalysis> byName)
    {
        string? ns           = ExtractQuoted(v.Detail, "Imports '");
        string? declaringAsm = v.TargetAssembly;
        if (ns is null || declaringAsm is null) return null;

        // More than one assembly declares this namespace (none referenced) — the suggested
        // ref names the natural owner, but the pick is less certain. See PickPrimaryDeclarer.
        bool ambiguous = v.Detail.Contains("also declared by", StringComparison.Ordinal);

        var hint = new JsonObject
        {
            ["affectedNamespace"] = ns,
            ["declaringAssembly"] = declaringAsm,
            ["confidence"]        = ambiguous ? "medium" : "high"
        };

        if (byName.TryGetValue(v.Project, out ProjectAnalysis? importing)
            && byName.TryGetValue(declaringAsm, out ProjectAnalysis? declaring))
        {
            string importingDir = Path.GetDirectoryName(importing.CsprojPath)!;
            hint["suggestedProjectReference"] =
                Path.GetRelativePath(importingDir, declaring.CsprojPath).Replace('\\', '/');
        }

        return hint;
    }

    // Detail format: "<ProjectReference> to 'asm' but none of [ns1, ns2] are imported"
    private static JsonObject? StaleHint(
        Violation v,
        Dictionary<string, ProjectAnalysis> byName)
    {
        string? staleAsm = ExtractQuoted(v.Detail, "<ProjectReference> to '");
        if (staleAsm is null) return null;

        var hint = new JsonObject
        {
            ["staleAssembly"] = staleAsm,
            ["confidence"]    = "high"
        };

        if (byName.TryGetValue(staleAsm, out ProjectAnalysis? stale))
            hint["unusedNamespaces"] = new JsonArray(
                stale.DeclaredNamespaces.OrderBy(n => n)
                    .Select(n => (JsonNode)JsonValue.Create(n)!).ToArray());

        return hint;
    }

    // Detail format: "Broken include 'old' — git shows rename to 'new' (hint text)"
    private static JsonObject? RenameHint(Violation v)
    {
        string? oldPath = ExtractQuoted(v.Detail, "Broken include '");
        string? newPath = ExtractQuoted(v.Detail, "rename to '");
        if (oldPath is null || newPath is null) return null;

        return new JsonObject
        {
            ["brokenPath"]       = oldPath,
            ["suggestedNewPath"] = newPath,
            ["newPathStatus"]    = v.Detail.Contains("unclaimed", StringComparison.Ordinal)
                                       ? "orphaned" : "claimed",
            ["confidence"]       = "high"
        };
    }

    // Detail format: "Version mismatch [v1 vs v2]: Asm1@v1, Asm2@v2"
    private static JsonObject? PkgConflictHint(Violation v)
    {
        int start = v.Detail.IndexOf('[');
        int end   = v.Detail.IndexOf(']');
        if (start < 0 || end <= start) return null;

        string[] versions = v.Detail[(start + 1)..end]
            .Split(" vs ", StringSplitOptions.RemoveEmptyEntries)
            .Select(s => s.Trim())
            .ToArray();

        return new JsonObject
        {
            ["versions"]   = new JsonArray(versions.Select(ver => (JsonNode)JsonValue.Create(ver)!).ToArray()),
            ["confidence"] = "medium"
        };
    }

    // ── Collated violations (file-keyed + assembly-level bucket) ──────────

    public static string WriteViolations(
        string artifactsDir,
        string repoRoot,
        IReadOnlyList<ProjectAnalysis> analyses,
        CrossReferenceIndex xref,
        IReadOnlyList<ProjectDiagnostics>? semanticResults)
    {
        Directory.CreateDirectory(artifactsDir);
        string path = Path.Combine(artifactsDir, "violations.json");

        var byName = analyses.ToDictionary(a => a.AssemblyName, StringComparer.Ordinal);

        // Bucket syntax violations: file-pinned vs assembly-level.
        var syntaxByFile = new SortedDictionary<string, List<Violation>>(StringComparer.Ordinal);
        var syntaxAssemblyLevel = new List<Violation>();
        foreach (Violation v in xref.Violations)
        {
            if (v.FilePath is null) { syntaxAssemblyLevel.Add(v); continue; }
            string rel = Rel(v.FilePath, repoRoot);
            if (!syntaxByFile.TryGetValue(rel, out List<Violation>? list))
                syntaxByFile[rel] = list = [];
            list.Add(v);
        }

        // Bucket semantic diagnostics: file-pinned vs project-level (null FilePath).
        var semanticByFile = new SortedDictionary<string, List<(ProjectDiagnostics Pd, CompilerDiagnostic D)>>(StringComparer.Ordinal);
        var semanticAssemblyLevel = new List<(ProjectDiagnostics Pd, CompilerDiagnostic D)>();
        int semErrorTotal = 0, semWarnTotal = 0;
        if (semanticResults is not null)
        {
            foreach (ProjectDiagnostics pd in semanticResults)
            {
                semErrorTotal += pd.Errors.Count;
                semWarnTotal  += pd.Warnings.Count;
                foreach (CompilerDiagnostic d in pd.Errors.Concat(pd.Warnings))
                {
                    if (d.FilePath is null) { semanticAssemblyLevel.Add((pd, d)); continue; }
                    string rel = Rel(d.FilePath, repoRoot);
                    if (!semanticByFile.TryGetValue(rel, out var list))
                        semanticByFile[rel] = list = [];
                    list.Add((pd, d));
                }
            }
        }

        // byFile: union of file keys; per file, emit syntax[] + rolled-up semantic[]
        var allFileKeys = new SortedSet<string>(StringComparer.Ordinal);
        foreach (string k in syntaxByFile.Keys)   allFileKeys.Add(k);
        foreach (string k in semanticByFile.Keys) allFileKeys.Add(k);

        var byFile = new JsonObject();
        foreach (string file in allFileKeys)
        {
            var syntaxArr = new JsonArray();
            if (syntaxByFile.TryGetValue(file, out var sList))
                foreach (Violation v in sList) syntaxArr.Add(BuildViolationNode(v, byName, repoRoot));

            var semanticArr = new JsonArray();
            if (semanticByFile.TryGetValue(file, out var dList))
                foreach (JsonObject entry in RollupSemanticByFile(dList))
                    semanticArr.Add(entry);

            byFile[file] = new JsonObject
            {
                ["syntax"]   = syntaxArr,
                ["semantic"] = semanticArr
            };
        }

        // assemblyLevel: dict keyed by ViolationKind for syntax; synthetic
        // "SemanticDiagnostic" key for any project-level semantic findings.
        var assemblyLevel = new JsonObject();
        foreach (var kindGroup in syntaxAssemblyLevel
            .GroupBy(v => v.Kind)
            .OrderBy(g => g.Key.ToString(), StringComparer.Ordinal))
        {
            var arr = new JsonArray();
            foreach (Violation v in kindGroup.OrderBy(v => v.Project, StringComparer.Ordinal))
                arr.Add(BuildViolationNode(v, byName, repoRoot));
            assemblyLevel[kindGroup.Key.ToString()] = arr;
        }
        if (semanticAssemblyLevel.Count > 0)
        {
            var arr = new JsonArray();
            foreach (JsonObject entry in RollupSemanticAssemblyLevel(semanticAssemblyLevel))
                arr.Add(entry);
            assemblyLevel["SemanticDiagnostic"] = arr;
        }

        var doc = new JsonObject
        {
            ["generatedAt"] = DateTime.UtcNow.ToString("O"),
            ["summary"] = new JsonObject
            {
                ["syntaxViolationCount"] = xref.Violations.Count,
                ["semanticErrorCount"]   = semErrorTotal,
                ["semanticWarningCount"] = semWarnTotal,
                ["filesAffected"]        = allFileKeys.Count,
                ["semanticIncluded"]     = semanticResults is not null
            },
            ["byFile"]        = byFile,
            ["assemblyLevel"] = assemblyLevel
        };

        File.WriteAllText(path, doc.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
        return path;
    }

    // Per-file semantic rollup: group by (id, severity, project) → inside each,
    // list distinct messages → each message carries its {line, column} locations.
    private static IEnumerable<JsonObject> RollupSemanticByFile(
        List<(ProjectDiagnostics Pd, CompilerDiagnostic D)> diags)
    {
        var groups = diags
            .GroupBy(t => (t.D.Id, t.D.Severity, t.Pd.AssemblyName))
            .OrderBy(g => SeverityRank(g.Key.Severity))
            .ThenBy(g => g.Key.Id,           StringComparer.Ordinal)
            .ThenBy(g => g.Key.AssemblyName, StringComparer.Ordinal);

        foreach (var group in groups)
        {
            var messages = new JsonArray();
            foreach (var msgGroup in group
                .GroupBy(t => t.D.Message, StringComparer.Ordinal)
                .OrderBy(g => g.Key, StringComparer.Ordinal))
            {
                var locations = new JsonArray();
                foreach (var item in msgGroup.OrderBy(t => t.D.Line).ThenBy(t => t.D.Column))
                    locations.Add(new JsonObject
                    {
                        ["line"]   = item.D.Line,
                        ["column"] = item.D.Column
                    });
                messages.Add(new JsonObject
                {
                    ["text"]      = msgGroup.Key,
                    ["locations"] = locations
                });
            }
            yield return new JsonObject
            {
                ["source"]   = "semantic",
                ["id"]       = group.Key.Id,
                ["severity"] = group.Key.Severity,
                ["project"]  = group.Key.AssemblyName,
                ["count"]    = group.Count(),
                ["messages"] = messages
            };
        }
    }

    // Assembly-level semantic (no FilePath): no locations, group by full key.
    private static IEnumerable<JsonObject> RollupSemanticAssemblyLevel(
        List<(ProjectDiagnostics Pd, CompilerDiagnostic D)> diags)
    {
        var groups = diags
            .GroupBy(t => (t.D.Id, t.D.Severity, t.Pd.AssemblyName, t.D.Message))
            .OrderBy(g => SeverityRank(g.Key.Severity))
            .ThenBy(g => g.Key.Id,           StringComparer.Ordinal)
            .ThenBy(g => g.Key.AssemblyName, StringComparer.Ordinal);

        foreach (var group in groups)
        {
            yield return new JsonObject
            {
                ["source"]   = "semantic",
                ["id"]       = group.Key.Id,
                ["severity"] = group.Key.Severity,
                ["project"]  = group.Key.AssemblyName,
                ["count"]    = group.Count(),
                ["message"]  = group.Key.Message
            };
        }
    }

    private static int SeverityRank(string severity) =>
        severity.Equals("Error",   StringComparison.OrdinalIgnoreCase) ? 0 :
        severity.Equals("Warning", StringComparison.OrdinalIgnoreCase) ? 1 :
        2;

    private static JsonObject BuildViolationNode(
        Violation v,
        Dictionary<string, ProjectAnalysis> byName,
        string repoRoot)
    {
        var node = new JsonObject
        {
            ["source"]  = "syntax",
            ["kind"]    = v.Kind.ToString(),
            ["project"] = v.Project,
            ["detail"]  = v.Detail
        };
        JsonObject? hint = BuildHint(v, byName, repoRoot);
        if (hint is not null) node["hint"] = hint;
        return node;
    }

    // ── Markdown health report ────────────────────────────────────────────

    private static void WriteHealth(
        string path, string repoRoot,
        IReadOnlyList<ProjectAnalysis> analyses,
        FileOwnershipResult ownership,
        CrossReferenceIndex xref)
    {
        var sb = new StringBuilder();
        var analysisByName = analyses.ToDictionary(a => a.AssemblyName, StringComparer.Ordinal);

        int totalTypes   = analyses.Sum(a => a.Types.Count);
        int partialCount = analyses.Sum(a => a.Types.Count(t => t.IsPartial));

        sb.AppendLine("# Project Health Report");
        sb.AppendLine($"Generated: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
        sb.AppendLine();

        // Summary
        sb.AppendLine("## Summary");
        sb.AppendLine();
        sb.AppendLine($"| Metric | Value |");
        sb.AppendLine($"|--------|-------|");
        sb.AppendLine($"| Projects | {analyses.Count} |");
        sb.AppendLine($"| Types declared | {totalTypes} |");
        sb.AppendLine($"| Partial families | {partialCount} |");
        sb.AppendLine($"| Orphaned files | {ownership.OrphanedFiles.Count} |");
        sb.AppendLine($"| Broken includes | {ownership.BrokenIncludes.Count} |");
        sb.AppendLine($"| Violations | {xref.Violations.Count} |");
        var artifactDir = Path.GetDirectoryName(path) ?? Path.Combine(repoRoot, "artifacts");
        string semanticStatus;
        string semanticFile = Path.Combine(artifactDir, "semantic-analysis.json");
        if (File.Exists(semanticFile))
        {
            semanticStatus = "errors/warnings";
        }
        else
        {
            semanticStatus = "disabled";
        }
        sb.AppendLine($"| Semantic analysis | {semanticStatus} |");
        sb.AppendLine();

        // Violations
        if (xref.Violations.Count > 0)
        {
            sb.AppendLine("## Violations");
            sb.AppendLine();

            var ghosts   = xref.Violations.Where(v => v.Kind == ViolationKind.GhostDependency).ToList();
            var stales   = xref.Violations.Where(v => v.Kind == ViolationKind.StaleProjectReference).ToList();
            var renames  = xref.Violations.Where(v => v.Kind == ViolationKind.RenameGap).ToList();
            var pkgConfs = xref.Violations.Where(v => v.Kind == ViolationKind.PackageVersionConflict).ToList();
            var others   = xref.Violations.Where(v => v.Kind is not ViolationKind.GhostDependency
                                                                and not ViolationKind.StaleProjectReference
                                                                and not ViolationKind.RenameGap
                                                                and not ViolationKind.PackageVersionConflict).ToList();

            if (ghosts.Count > 0)
            {
                sb.AppendLine($"### Ghost Dependencies ({ghosts.Count})");
                sb.AppendLine();
                sb.AppendLine("These projects import namespaces from an assembly with no `<ProjectReference>` declared.");
                sb.AppendLine();
                foreach (Violation v in ghosts.OrderBy(v => v.Project))
                {
                    sb.AppendLine($"**`{v.Project}`** — {v.Detail}");

                    string? declaringAssembly = v.TargetAssembly;
                    if (declaringAssembly != null && analysisByName.TryGetValue(v.Project, out ProjectAnalysis? importing)
                        && analysisByName.TryGetValue(declaringAssembly, out ProjectAnalysis? declaring))
                    {
                        string importingDir = Path.GetDirectoryName(importing.CsprojPath)!;
                        string relRef = Path.GetRelativePath(importingDir, declaring.CsprojPath).Replace('\\', '/');
                        sb.AppendLine($"```xml");
                        sb.AppendLine($"<ProjectReference Include=\"{relRef}\" />");
                        sb.AppendLine($"```");
                    }
                    sb.AppendLine();
                }
            }

            if (stales.Count > 0)
            {
                sb.AppendLine($"### Stale Project References ({stales.Count})");
                sb.AppendLine();
                sb.AppendLine("These projects declare a `<ProjectReference>` but import none of its namespaces.");
                sb.AppendLine();
                foreach (Violation v in stales.OrderBy(v => v.Project))
                {
                    sb.AppendLine($"**`{v.Project}`** — {v.Detail}");
                    sb.AppendLine();
                }
            }

            if (pkgConfs.Count > 0)
            {
                sb.AppendLine($"### Package Version Conflicts ({pkgConfs.Count})");
                sb.AppendLine();
                sb.AppendLine("The same NuGet package is referenced at different versions across projects.");
                sb.AppendLine();
                foreach (Violation v in pkgConfs.OrderBy(v => v.Project))
                    sb.AppendLine($"**`{v.Project}`** — {v.Detail}");
                sb.AppendLine();
            }

            if (renames.Count > 0)
            {
                sb.AppendLine($"### Rename Gaps ({renames.Count})");
                sb.AppendLine();
                sb.AppendLine("These broken includes correspond to files that git shows were recently renamed.");
                sb.AppendLine();
                foreach (Violation v in renames.OrderBy(v => v.Project))
                    sb.AppendLine($"**`{v.Project}`** — {v.Detail}");
                sb.AppendLine();
            }

            if (others.Count > 0)
            {
                sb.AppendLine($"### Other Violations ({others.Count})");
                sb.AppendLine();
                foreach (Violation v in others.OrderBy(v => v.Kind.ToString()).ThenBy(v => v.Project))
                    sb.AppendLine($"- `[{v.Kind}]` **{v.Project}**: {v.Detail}");
                sb.AppendLine();
            }
        }
        else
        {
            sb.AppendLine("## Violations");
            sb.AppendLine();
            sb.AppendLine("None. All project references match imported namespaces.");
            sb.AppendLine();
        }

        // Orphaned files
        if (ownership.OrphanedFiles.Count > 0)
        {
            sb.AppendLine($"## Orphaned Files ({ownership.OrphanedFiles.Count})");
            sb.AppendLine();
            sb.AppendLine("These `.cs` files exist under `src/` but are not claimed by any `<Compile Include>` pattern.");
            sb.AppendLine();
            foreach (string f in ownership.OrphanedFiles.OrderBy(f => f))
                sb.AppendLine($"- `{Rel(f, repoRoot)}`");
            sb.AppendLine();
        }

        // Broken includes
        if (ownership.BrokenIncludes.Count > 0)
        {
            sb.AppendLine($"## Broken Includes ({ownership.BrokenIncludes.Count})");
            sb.AppendLine();
            foreach (BrokenInclude b in ownership.BrokenIncludes)
                sb.AppendLine($"- `{Rel(b.MissingAbsPath, repoRoot)}` (`{b.AssemblyName}`) — does not exist on disk");
            sb.AppendLine();
        }

        // Partial families
        var partialFamilies = analyses
            .SelectMany(a => a.Types.Where(t => t.IsPartial)
                .Select(t => (Assembly: a.AssemblyName, Type: t)))
            .OrderBy(x => x.Assembly).ThenBy(x => x.Type.Name)
            .ToList();

        if (partialFamilies.Count > 0)
        {
            sb.AppendLine($"## Partial Families ({partialFamilies.Count})");
            sb.AppendLine();
            foreach (var (assembly, t) in partialFamilies)
            {
                sb.AppendLine($"**`{t.ContainingNamespace}.{t.Name}`** ({assembly}) — {t.ContributingFiles.Count} files");
                foreach (string f in t.ContributingFiles)
                    sb.AppendLine($"  - `{Rel(f, repoRoot)}`");
                sb.AppendLine();
            }
        }

        // Hotspot namespaces — widest blast radius for structural change
        var hotNs = xref.ImportersOf
            .Where(kv => kv.Value.Count >= 2)
            .OrderByDescending(kv => kv.Value.Count)
            .ThenBy(kv => kv.Key, StringComparer.Ordinal)
            .Take(10)
            .ToList();

        if (hotNs.Count > 0)
        {
            sb.AppendLine("## Hotspot Namespaces");
            sb.AppendLine();
            sb.AppendLine("Namespaces imported by multiple assemblies — changes here have the widest blast radius.");
            sb.AppendLine();
            sb.AppendLine("| Namespace | Importers | Imported By |");
            sb.AppendLine("|-----------|-----------|-------------|");
            foreach (var kv in hotNs)
            {
                string importers = string.Join(", ", kv.Value.OrderBy(n => n).Select(n => $"`{n}`"));
                sb.AppendLine($"| `{kv.Key}` | {kv.Value.Count} | {importers} |");
            }
            sb.AppendLine();
        }

        // Project overview table
        sb.AppendLine("## Project Overview");
        sb.AppendLine();
        sb.AppendLine("| Project | Namespaces | Types | Partial | References |");
        sb.AppendLine("|---------|-----------|-------|---------|-----------|");
        foreach (ProjectAnalysis a in analyses.OrderBy(a => a.AssemblyName))
        {
            int partial = a.Types.Count(t => t.IsPartial);
            sb.AppendLine($"| `{a.AssemblyName}` | {a.DeclaredNamespaces.Count} | {a.Types.Count} | {partial} | {a.ProjectReferences.Count} |");
        }
        sb.AppendLine();

        File.WriteAllText(path, sb.ToString());
    }

    // ── Semantic diagnostics ──────────────────────────────────────────────

    public static string WriteSemantic(
        string artifactsDir,
        string repoRoot,
        IReadOnlyList<ProjectDiagnostics> results)
    {
        Directory.CreateDirectory(artifactsDir);
        string path = Path.Combine(artifactsDir, "semantic-analysis.json");

        int totalErrors   = results.Sum(r => r.Errors.Count);
        int totalWarnings = results.Sum(r => r.Warnings.Count);

        var projectsNode = new JsonObject();
        foreach (ProjectDiagnostics pd in results)
        {
            projectsNode[pd.AssemblyName] = new JsonObject
            {
                ["errors"]   = DiagArray(pd.Errors,   repoRoot),
                ["warnings"] = DiagArray(pd.Warnings, repoRoot)
            };
        }

        var doc = new JsonObject
        {
            ["generatedAt"] = DateTime.UtcNow.ToString("O"),
            ["summary"] = new JsonObject
            {
                ["projectsWithErrors"]   = results.Count(r => r.Errors.Count > 0),
                ["projectsWithWarnings"] = results.Count(r => r.Warnings.Count > 0),
                ["totalErrors"]          = totalErrors,
                ["totalWarnings"]        = totalWarnings
            },
            ["projects"] = projectsNode
        };

        doc["status"] = "succeeded";
        File.WriteAllText(path, doc.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
        return path;
    }

    public static string WriteSemanticFailure(string artifactsDir, string repoRoot, Exception ex)
    {
        Directory.CreateDirectory(artifactsDir);
        string path = Path.Combine(artifactsDir, "semantic-analysis.json");

        var doc = new JsonObject
        {
            ["generatedAt"] = DateTime.UtcNow.ToString("O"),
            ["status"]      = "failed",
            ["error"]       = ex.ToString()
        };

        File.WriteAllText(path, doc.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
        return path;
    }

    private static JsonArray DiagArray(IReadOnlyList<CompilerDiagnostic> diags, string repoRoot)
    {
        return new JsonArray(diags.Select(d => (JsonNode)new JsonObject
        {
            ["id"]       = d.Id,
            ["message"]  = d.Message,
            ["file"]     = d.FilePath != null ? Rel(d.FilePath, repoRoot) : null,
            ["line"]     = d.Line,
            ["column"]   = d.Column
        }).ToArray());
    }

    // ── Helpers ───────────────────────────────────────────────────────────

    private static string Rel(string absolutePath, string repoRoot) =>
        Path.GetRelativePath(repoRoot, absolutePath).Replace('\\', '/');

    private static JsonArray JsonArr(IEnumerable<string> items) =>
        new(items.Select(s => (JsonNode)JsonValue.Create(s)!).ToArray());

    private static string EscapeMd(string s) => s.Replace("|", "\\|");

    private static string? ExtractQuoted(string text, string prefix)
    {
        int start = text.IndexOf(prefix, StringComparison.Ordinal);
        if (start < 0) return null;
        start += prefix.Length;
        int end = text.IndexOf('\'', start);
        return end > start ? text[start..end] : null;
    }
}
