using System.Collections.Generic;
using System.Linq;

namespace RepoAudit;

// Per-dependency slot: which types from a referenced project are
// actually referenced by name in type positions in the dependent.
public sealed record LinkedDependency(
    string ReferencedAssembly,
    IReadOnlyList<string> UsedTypes,
    bool IsTypeEffective);   // true when UsedTypes.Count > 0

// Linking footprint for one project — one entry per ProjectReference.
public sealed record ProjectLinkage(
    string Assembly,
    IReadOnlyList<LinkedDependency> Dependencies);

public static class SymbolLinker
{
    // Build per-project linking report.
    //
    // For each <ProjectReference> in a project, we intersect the set of type
    // names referenced in type-position syntax nodes across all of that
    // project's source files with the set of type names declared by the
    // referenced assembly.  The result tells you which types are actually
    // pulling in each dependency — or that none are (candidate dead reference).
    //
    // Note: this is a syntactic approximation.  A non-empty UsedTypes list is
    // strong evidence of real usage.  An empty list is a useful signal but can
    // have false positives when types are accessed via fully-qualified names
    // that weren't resolved to their rightmost segment, or via reflection.
    public static IReadOnlyList<ProjectLinkage> BuildLinkage(
        IReadOnlyList<ProjectAnalysis> projects)
    {
        // Assembly name → set of simple type names it declares
        var typesByAssembly = projects.ToDictionary(
            p => p.AssemblyName,
            p => (IReadOnlySet<string>) new HashSet<string>(
                p.Types.Select(t => t.Name), System.StringComparer.Ordinal),
            System.StringComparer.Ordinal);

        var result = new List<ProjectLinkage>();

        foreach (ProjectAnalysis project in projects)
        {
            if (project.ProjectReferences.Count == 0) continue;

            // Union all type-position names referenced anywhere in this project
            var allRefNames = new HashSet<string>(System.StringComparer.Ordinal);
            foreach (FileAnalysis fa in project.FileAnalyses)
                foreach (string name in fa.ReferencedTypeNames)
                    allRefNames.Add(name);

            var deps = new List<LinkedDependency>(project.ProjectReferences.Count);
            foreach (string refAssembly in project.ProjectReferences)
            {
                if (!typesByAssembly.TryGetValue(refAssembly, out IReadOnlySet<string>? aTypes))
                    continue;

                var used = allRefNames
                    .Where(n => aTypes.Contains(n))
                    .OrderBy(n => n, System.StringComparer.Ordinal)
                    .ToList();

                deps.Add(new LinkedDependency(refAssembly, used.AsReadOnly(), used.Count > 0));
            }

            if (deps.Count > 0)
                result.Add(new ProjectLinkage(project.AssemblyName, deps.AsReadOnly()));
        }

        return result.AsReadOnly();
    }

    // Returns projects whose every ProjectReference has IsTypeEffective == false.
    // These are candidates for having all their references cleaned up or reviewed.
    public static IReadOnlyList<string> FindEffectivelyUnlinkedProjects(
        IReadOnlyList<ProjectLinkage> linkage) =>
        linkage
            .Where(pl => pl.Dependencies.Count > 0 && pl.Dependencies.All(d => !d.IsTypeEffective))
            .Select(pl => pl.Assembly)
            .OrderBy(n => n, System.StringComparer.Ordinal)
            .ToList()
            .AsReadOnly();
}
