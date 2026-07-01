using System;
using System.Collections.Generic;
using System.Linq;
using RepoAudit;
using Xunit;

namespace RepoAudit.Tests;

// Regression coverage for the ghost-dependency detector and its declarer selection.
//
// Two facts make "is this a ghost?" subtler than a 1:1 namespace↔assembly lookup, and
// both have bitten this analyzer before:
//   * a namespace can be declared by more than one assembly (shared abstraction files
//     are compiled into multiple projects in this repo), and
//   * project references flow transitively at build time.
// These tests pin the resulting rules so a future refactor can't silently regress them.
public class CrossReferenceGhostTests
{
    // ── PickPrimaryDeclarer (internal; reachable via InternalsVisibleTo) ────

    [Fact]
    public void PickPrimaryDeclarer_PrefersExactNamespaceOwner()
    {
        // 'Graphs.Distance' is declared by both its owning assembly and a co-declarer
        // that compiles the shared abstraction in. The exact name match must win.
        string pick = CrossReference.PickPrimaryDeclarer(
            "Graphs.Distance",
            new[] { "Graphs.Primitives", "Graphs.Distance" });
        Assert.Equal("Graphs.Distance", pick);
    }

    [Fact]
    public void PickPrimaryDeclarer_PrefersLongestDottedPrefix()
    {
        string pick = CrossReference.PickPrimaryDeclarer(
            "Foo.Bar.Baz",
            new[] { "Zoo", "Foo", "Foo.Bar" });
        Assert.Equal("Foo.Bar", pick);
    }

    [Fact]
    public void PickPrimaryDeclarer_FallsBackToOrdinalFirst()
    {
        // No exact or dotted-prefix match -> deterministic ordinal-first.
        string pick = CrossReference.PickPrimaryDeclarer(
            "Xyz.Thing",
            new[] { "Beta", "Alpha" });
        Assert.Equal("Alpha", pick);
    }

    [Fact]
    public void PickPrimaryDeclarer_IsInputOrderIndependent()
    {
        string a = CrossReference.PickPrimaryDeclarer("N.S", new[] { "B", "A", "N.S" });
        string b = CrossReference.PickPrimaryDeclarer("N.S", new[] { "N.S", "A", "B" });
        Assert.Equal(a, b);
        Assert.Equal("N.S", a);
    }

    // ── GhostDependency detection via CrossReference.Build ──────────────────

    [Fact]
    public void GenuineGhost_Fires_AndCarriesDeclaringAssemblyStructurally()
    {
        // A imports namespace 'B' (declared only by assembly B) with no reference to it.
        var projects = new[]
        {
            Proj("A", declared: new[] { "A" }, imported: new[] { "B" }, projectRefs: Array.Empty<string>()),
            Proj("B", declared: new[] { "B" }, imported: Array.Empty<string>(), projectRefs: Array.Empty<string>()),
        };

        Violation ghost = Ghosts(projects).Single();
        Assert.Equal("A", ghost.Project);
        Assert.Equal("B", ghost.TargetAssembly);   // structural, not re-parsed from Detail
        Assert.Contains("Imports 'B'", ghost.Detail);
    }

    [Fact]
    public void InAssemblyDeclaration_IsNotAGhost()
    {
        // A imports 'Shared' and also declares it (shared file compiled into A) -> in-assembly.
        var projects = new[]
        {
            Proj("A",      declared: new[] { "A", "Shared" }, imported: new[] { "Shared" }, projectRefs: Array.Empty<string>()),
            Proj("Shared", declared: new[] { "Shared" },      imported: Array.Empty<string>(), projectRefs: Array.Empty<string>()),
        };
        Assert.Empty(Ghosts(projects));
    }

    [Fact]
    public void DirectReferenceThatDeclaresNamespace_IsNotAGhost()
    {
        // A imports 'Shared'; A references Lib; Lib declares 'Shared'. Satisfied.
        // (The Graphs.Proximity -> Graphs.Primitives pattern.)
        var projects = new[]
        {
            Proj("A",   declared: new[] { "A" },      imported: new[] { "Shared" }, projectRefs: new[] { "Lib" }),
            Proj("Lib", declared: new[] { "Shared" }, imported: Array.Empty<string>(), projectRefs: Array.Empty<string>()),
        };
        Assert.Empty(Ghosts(projects));
    }

    [Fact]
    public void ExternalNamespace_DeclaredByNoProject_IsNotAGhost()
    {
        // Framework / NuGet namespaces are not project ghosts (that's a PackageReference concern).
        var projects = new[]
        {
            Proj("A", declared: new[] { "A" },
                 imported: new[] { "System.Text", "Newtonsoft.Json" },
                 projectRefs: Array.Empty<string>()),
        };
        Assert.Empty(Ghosts(projects));
    }

    [Fact]
    public void MultipleUnreferencedDeclarers_NamesOwner_AndListsAlternates()
    {
        // 'Graphs.Distance' is declared by both assemblies; A references neither. The report
        // must name the owner deterministically and surface the alternate, not pick silently.
        var projects = new[]
        {
            Proj("A",                 declared: new[] { "A" },               imported: new[] { "Graphs.Distance" }, projectRefs: Array.Empty<string>()),
            Proj("Graphs.Distance",   declared: new[] { "Graphs.Distance" }, imported: Array.Empty<string>(),       projectRefs: Array.Empty<string>()),
            Proj("Graphs.Primitives", declared: new[] { "Graphs.Distance" }, imported: Array.Empty<string>(),       projectRefs: Array.Empty<string>()),
        };

        Violation ghost = Ghosts(projects).Single();
        Assert.Equal("Graphs.Distance", ghost.TargetAssembly);                 // owner, deterministically
        Assert.Contains("also declared by 'Graphs.Primitives'", ghost.Detail); // alternate surfaced
    }

    [Fact]
    public void TransitivelyReachableNamespace_StillFires()
    {
        // A -> C -> B; B declares 'B'; A imports 'B' but has no DIRECT reference to B.
        // A real build compiles (transitive flow), but the lint intentionally flags reliance
        // on the transitive graph — declare the dependency explicitly. This is by design.
        var projects = new[]
        {
            Proj("A", declared: new[] { "A" }, imported: new[] { "B" }, projectRefs: new[] { "C" }),
            Proj("C", declared: new[] { "C" }, imported: Array.Empty<string>(), projectRefs: new[] { "B" }),
            Proj("B", declared: new[] { "B" }, imported: Array.Empty<string>(), projectRefs: Array.Empty<string>()),
        };

        Violation ghost = Ghosts(projects).Single(g => g.Project == "A");
        Assert.Equal("B", ghost.TargetAssembly);
    }

    // ── Helpers ─────────────────────────────────────────────────────────────

    private static IReadOnlyList<Violation> Ghosts(IReadOnlyList<ProjectAnalysis> projects) =>
        CrossReference.Build(projects).Violations
            .Where(v => v.Kind == ViolationKind.GhostDependency)
            .ToList();

    private static ProjectAnalysis Proj(
        string name, string[] declared, string[] imported, string[] projectRefs) =>
        new(
            CsprojPath:         $"/repo/projects/{name}/{name}.csproj",
            AssemblyName:       name,
            DeclaredNamespaces: new HashSet<string>(declared, StringComparer.Ordinal),
            ImportedNamespaces: new HashSet<string>(imported, StringComparer.Ordinal),
            ProjectReferences:  projectRefs,
            Types:              Array.Empty<TypeEntry>(),
            FileAnalyses:       Array.Empty<FileAnalysis>(),
            PackageReferences:  Array.Empty<PackageRef>(),
            CompilerSettings:   new CsprojCompilerSettings(null, null, null));
}
