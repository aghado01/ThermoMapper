using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace RepoAudit;

public enum ViolationKind
{
    StaleProjectReference,
    GhostDependency,
    NamespacePathDrift,
    PartialFamilyOrphan,
    RenameGap,
    PackageVersionConflict,
    SolutionGap
}

public sealed record Violation(
    ViolationKind Kind,
    string Project,
    string Detail,
    string? FilePath = null,
    string? TargetAssembly = null);

public sealed record CrossReferenceIndex(
    IReadOnlyDictionary<string, IReadOnlySet<string>> DeclaredByAssembly,
    IReadOnlyDictionary<string, IReadOnlyList<string>> ImportersOf,
    IReadOnlyList<Violation> Violations);

public static class CrossReference
{
    public static CrossReferenceIndex Build(IReadOnlyList<ProjectAnalysis> projects)
    {
        var declaredByAssembly = projects.ToDictionary(
            p => p.AssemblyName,
            p => p.DeclaredNamespaces,
            System.StringComparer.Ordinal);

        var importersOf = new Dictionary<string, List<string>>(System.StringComparer.Ordinal);
        foreach (ProjectAnalysis project in projects)
            foreach (string ns in project.ImportedNamespaces)
            {
                if (!importersOf.TryGetValue(ns, out List<string>? list))
                    importersOf[ns] = list = [];
                list.Add(project.AssemblyName);
            }

        var violations = new List<Violation>();

        // Inverse index: namespace -> assemblies that declare it. A namespace can be
        // declared by more than one assembly in this repo (shared abstraction files are
        // compiled into multiple projects — e.g. the Graphs.Distance abstractions live
        // in both Graphs.Distance and Graphs.Primitives), so ghost detection must
        // consider ALL declarers rather than assume a 1:1 namespace↔assembly mapping.
        var declarersOf = new Dictionary<string, List<string>>(System.StringComparer.Ordinal);
        foreach (ProjectAnalysis declarer in projects)
            foreach (string ns in declarer.DeclaredNamespaces)
            {
                if (!declarersOf.TryGetValue(ns, out List<string>? list))
                    declarersOf[ns] = list = [];
                list.Add(declarer.AssemblyName);
            }

        foreach (ProjectAnalysis project in projects)
        {
            // Stale reference: <ProjectReference> declared but none of its namespaces are imported
            foreach (string refAssembly in project.ProjectReferences)
            {
                if (!declaredByAssembly.TryGetValue(refAssembly, out IReadOnlySet<string>? refNs)) continue;
                if (project.ImportedNamespaces.Count == 0) continue;

                if (!project.ImportedNamespaces.Any(ns => refNs.Contains(ns)))
                    violations.Add(new Violation(
                        ViolationKind.StaleProjectReference,
                        project.AssemblyName,
                        $"<ProjectReference> to '{refAssembly}' but none of [{string.Join(", ", refNs)}] are imported",
                        TargetAssembly: refAssembly));
            }

            // Ghost dependency: imports a namespace whose declaring assembly is reachable
            // through neither an in-assembly declaration nor a direct <ProjectReference>.
            // A namespace satisfied in-assembly (shared files compiled directly in) or by
            // a directly-referenced project that also declares it is NOT a ghost — only an
            // import with no available satisfier is.
            foreach (string usedNs in project.ImportedNamespaces)
            {
                if (!declarersOf.TryGetValue(usedNs, out List<string>? declarers)) continue;  // external/framework ns
                if (project.DeclaredNamespaces.Contains(usedNs)) continue;                     // satisfied in-assembly
                if (declarers.Any(a => project.ProjectReferences.Contains(a, System.StringComparer.Ordinal)))
                    continue;                                                                  // satisfied by a direct ref

                // A namespace can have several unreferenced declarers (shared source compiled
                // into multiple assemblies). Name the natural owner deterministically, and when
                // there is more than one candidate list the rest so the report can never quietly
                // point at the wrong project. The chosen owner also travels in TargetAssembly so
                // downstream consumers read it structurally instead of re-parsing Detail.
                string primary = PickPrimaryDeclarer(usedNs, declarers);
                string detail;
                if (declarers.Count == 1)
                {
                    detail = $"Imports '{usedNs}' (declared by '{primary}') with no <ProjectReference>";
                }
                else
                {
                    string others = string.Join(", ", declarers
                        .Where(a => a != primary)
                        .OrderBy(a => a, System.StringComparer.Ordinal)
                        .Select(a => $"'{a}'"));
                    detail = $"Imports '{usedNs}' (declared by '{primary}'; also declared by {others}) " +
                             "with no <ProjectReference>";
                }

                violations.Add(new Violation(
                    ViolationKind.GhostDependency,
                    project.AssemblyName,
                    detail,
                    TargetAssembly: primary));
            }
        }

        // Package version conflicts: same package, different versions across projects
        var pkgUsages = new Dictionary<string, List<(string Version, string Assembly)>>(
            System.StringComparer.OrdinalIgnoreCase);

        foreach (ProjectAnalysis project in projects)
            foreach (PackageRef pkg in project.PackageReferences)
            {
                if (!pkgUsages.TryGetValue(pkg.Name, out var list))
                    pkgUsages[pkg.Name] = list = [];
                list.Add((pkg.Version, project.AssemblyName));
            }

        foreach ((string pkgName, var usages) in pkgUsages)
        {
            var versions = usages
                .Select(u => u.Version)
                .Distinct(System.StringComparer.OrdinalIgnoreCase)
                .OrderBy(v => v)
                .ToList();
            if (versions.Count <= 1) continue;

            string detail = string.Join(", ",
                usages.OrderBy(u => u.Assembly).Select(u => $"{u.Assembly}@{u.Version}"));
            violations.Add(new Violation(
                ViolationKind.PackageVersionConflict,
                pkgName,
                $"Version mismatch [{string.Join(" vs ", versions)}]: {detail}"));
        }

        return new CrossReferenceIndex(
            declaredByAssembly.ToDictionary(
                kvp => kvp.Key,
                kvp => (IReadOnlySet<string>)kvp.Value),
            importersOf.ToDictionary(
                kvp => kvp.Key,
                kvp => (IReadOnlyList<string>)kvp.Value.AsReadOnly()),
            violations.AsReadOnly());
    }

    // Returns projects sorted so every dependency appears before its dependents.
    // Cycles are broken by alphabetical order — the first encountered node wins.
    public static IReadOnlyList<ProjectAnalysis> TopologicalSort(
        IReadOnlyList<ProjectAnalysis> projects)
    {
        var byName  = projects.ToDictionary(p => p.AssemblyName, System.StringComparer.Ordinal);
        var result  = new List<ProjectAnalysis>(projects.Count);
        var visited  = new HashSet<string>(System.StringComparer.Ordinal);
        var visiting = new HashSet<string>(System.StringComparer.Ordinal);

        void Visit(ProjectAnalysis p)
        {
            if (visited.Contains(p.AssemblyName))  return;
            if (visiting.Contains(p.AssemblyName)) { result.Add(p); visited.Add(p.AssemblyName); return; }

            visiting.Add(p.AssemblyName);
            foreach (string dep in p.ProjectReferences.OrderBy(d => d, System.StringComparer.Ordinal))
                if (byName.TryGetValue(dep, out ProjectAnalysis? depProject))
                    Visit(depProject);
            visiting.Remove(p.AssemblyName);
            visited.Add(p.AssemblyName);
            result.Add(p);
        }

        foreach (ProjectAnalysis p in projects.OrderBy(p => p.AssemblyName, System.StringComparer.Ordinal))
            Visit(p);

        return result.AsReadOnly();
    }

    // ── NamespacePathDrift ────────────────────────────────────────────────
    // Flags files whose directory path (relative to srcRoot) contains a
    // sub-folder component that does not appear — even loosely — in any of
    // the file's declared namespaces.  Only fires when a file is more than
    // one level deep from srcRoot so flat projects are never touched.

    public static IReadOnlyList<Violation> DetectNamespaceDrift(
        IReadOnlyList<ProjectAnalysis> projects,
        string srcRoot)
    {
        var violations = new List<Violation>();

        // Folder names that are structural (no namespace meaning)
        var structural = new System.Collections.Generic.HashSet<string>(System.StringComparer.OrdinalIgnoreCase)
            { "src", "lib", "source", "code", "tests", "test", "unit", "integration", "bin", "obj", "debug", "release" };

        foreach (ProjectAnalysis project in projects)
        foreach (FileAnalysis file in project.FileAnalyses)
        {
            if (file.DeclaredNamespaces.Count == 0) continue;

            string rel = Path.GetRelativePath(srcRoot, file.FilePath).Replace('\\', '/');
            if (rel.StartsWith("..")) continue;            // outside srcRoot — skip

            string[] parts = rel.Split('/');
            if (parts.Length <= 2) continue;              // ≤1 dir level — nothing to check

            // Skip the outermost component (project root folder under src)
            // and the filename; check only intermediate directory components
            IEnumerable<string> dirComponents = parts[1..^1]
                .Where(p => p.Length > 0 && !structural.Contains(p));

            foreach (string dir in dirComponents)
            {
                string normDir = NormalizeName(dir);
                bool matched  = file.DeclaredNamespaces.Any(ns =>
                    ns.Split('.').Any(seg => NamesOverlap(NormalizeName(seg), normDir)));

                if (!matched)
                {
                    violations.Add(new Violation(
                        ViolationKind.NamespacePathDrift,
                        project.AssemblyName,
                        $"'{rel}' is in folder '{dir}' but no declared namespace reflects this " +
                        $"(declared: {string.Join(", ", file.DeclaredNamespaces)})",
                        FilePath: file.FilePath));
                    break;  // one violation per file is enough
                }
            }
        }

        return violations.AsReadOnly();
    }

    // ── PartialFamilyOrphan ───────────────────────────────────────────────
    // Flags when the same fully-qualified partial type name is claimed by
    // more than one assembly — the partial family is split across a project
    // boundary, which will either fail to compile or silently diverge.

    public static IReadOnlyList<Violation> DetectPartialOrphans(
        IReadOnlyList<ProjectAnalysis> projects)
    {
        var violations = new List<Violation>();

        // Group partial types by fully-qualified type name across all projects
        var byFqtn = new Dictionary<string, List<string>>(System.StringComparer.Ordinal);

        foreach (ProjectAnalysis project in projects)
        foreach (TypeEntry type in project.Types.Where(t => t.IsPartial))
        {
            string fqtn = string.IsNullOrEmpty(type.ContainingNamespace)
                ? type.Name
                : $"{type.ContainingNamespace}.{type.Name}";

            if (!byFqtn.TryGetValue(fqtn, out List<string>? assemblyList))
                byFqtn[fqtn] = assemblyList = [];

            assemblyList.Add(project.AssemblyName);
        }

        foreach ((string fqtn, List<string> assemblies) in byFqtn)
        {
            var distinct = assemblies.Distinct(System.StringComparer.Ordinal).OrderBy(a => a).ToList();
            if (distinct.Count <= 1) continue;

            string asmList = string.Join(", ", distinct.Select(a => $"'{a}'"));
            violations.Add(new Violation(
                ViolationKind.PartialFamilyOrphan,
                distinct[0],
                $"Partial type '{fqtn}' is split across assemblies: {asmList}"));
        }

        return violations.AsReadOnly();
    }

    // ── SolutionGap ───────────────────────────────────────────────────────
    // Flags any projects/ or tests/ *.csproj that exists on disk but is not listed
    // in a .sln at repoRoot. If no .sln is present, the check is skipped.

    private static readonly System.Text.RegularExpressions.Regex SolutionProjectLine =
        new(@"Project\(""\{[A-F0-9\-]+\}""\)\s*=\s*""[^""]*"",\s*""(?<path>[^""]+.csproj)""",
            System.Text.RegularExpressions.RegexOptions.IgnoreCase |
            System.Text.RegularExpressions.RegexOptions.Compiled);

    public static IReadOnlyList<Violation> DetectSolutionGaps(
        string repoRoot,
        IReadOnlyList<string> csprojFiles)
    {
        var violations = new List<Violation>();

        string[] slnFiles = Directory.Exists(repoRoot)
            ? Directory.GetFiles(repoRoot, "*.sln", SearchOption.TopDirectoryOnly)
            : [];
        if (slnFiles.Length == 0) return violations.AsReadOnly();

        var listed = new System.Collections.Generic.HashSet<string>(System.StringComparer.OrdinalIgnoreCase);
        foreach (string slnPath in slnFiles)
        {
            foreach (string line in File.ReadLines(slnPath))
            {
                var m = SolutionProjectLine.Match(line);
                if (!m.Success) continue;
                string rel = m.Groups["path"].Value.Replace('/', '\\');
                listed.Add(Path.GetFullPath(Path.Combine(repoRoot, rel)));
            }
        }

        // csproj are first-class under both projects/ (production) and tests/ (test & smoke
        // harnesses); a project missing from the .sln under either root is a gap.
        string projectsRoot = Path.Combine(repoRoot, "projects");
        string testsRoot    = Path.Combine(repoRoot, "tests");
        foreach (string csproj in csprojFiles)
        {
            string full = Path.GetFullPath(csproj);
            bool inScope =
                full.StartsWith(projectsRoot, System.StringComparison.OrdinalIgnoreCase) ||
                full.StartsWith(testsRoot,    System.StringComparison.OrdinalIgnoreCase);
            if (!inScope) continue;
            if (listed.Contains(full)) continue;

            string assemblyName = Path.GetFileNameWithoutExtension(csproj);
            string relForReport = Path.GetRelativePath(repoRoot, full).Replace('\\', '/');
            string slnList = string.Join(", ",
                slnFiles.Select(s => Path.GetFileName(s)));
            violations.Add(new Violation(
                ViolationKind.SolutionGap,
                assemblyName,
                $"'{relForReport}' exists on disk but is not listed in {slnList}",
                FilePath: full));
        }

        return violations.AsReadOnly();
    }

    // ── Helpers ───────────────────────────────────────────────────────────

    // Picks which declaring assembly to name when a namespace is declared by more than
    // one unreferenced project. Deterministic and intent-revealing — prefer the assembly
    // that owns the namespace by name, so the suggested <ProjectReference> points at the
    // natural home rather than an arbitrary co-declarer:
    //   1. exact:    assembly name == namespace            (ns 'Graphs.Distance' -> asm 'Graphs.Distance')
    //   2. prefix:   assembly is a dotted ancestor of ns   (ns 'Foo.Bar.Baz' -> asm 'Foo.Bar' over 'Zoo')
    //   3. fallback: ordinal-first
    // Every tier resolves ties by ordinal name order, so the result is stable across runs
    // and platforms regardless of project enumeration order.
    internal static string PickPrimaryDeclarer(string ns, IReadOnlyList<string> declarers)
    {
        var ordered = declarers.OrderBy(a => a, System.StringComparer.Ordinal).ToList();

        string? exact = ordered.FirstOrDefault(a => string.Equals(a, ns, System.StringComparison.Ordinal));
        if (exact != null) return exact;

        string? prefix = ordered
            .Where(a => ns.StartsWith(a + ".", System.StringComparison.Ordinal))
            .OrderByDescending(a => a.Length)
            .FirstOrDefault();
        if (prefix != null) return prefix;

        return ordered[0];
    }

    // Lowercase + keep only alphanumerics — normalises hyphens, dots, underscores away
    private static string NormalizeName(string s) =>
        new string(s.ToLowerInvariant().Where(char.IsLetterOrDigit).ToArray());

    // True when either normalised string is a substring of the other
    private static bool NamesOverlap(string a, string b) =>
        a.Contains(b, System.StringComparison.Ordinal) ||
        b.Contains(a, System.StringComparison.Ordinal);
}
