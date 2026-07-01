using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Loader;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;

namespace RepoAudit;

public sealed record CompilerDiagnostic(
    string Id,
    string Severity,
    string Message,
    string? FilePath,
    int Line,
    int Column);

public sealed record ProjectDiagnostics(
    string AssemblyName,
    IReadOnlyList<CompilerDiagnostic> Errors,
    IReadOnlyList<CompilerDiagnostic> Warnings,
    IReadOnlySet<string> UsedProjectReferences);

public static class InProcessCompiler
{
    private const string TargetFramework = "net10.0";

    internal sealed record GeneratorBundle(
        IReadOnlyList<ISourceGenerator> Generators,
        AnalyzerAssemblyLoadContext? LoadContext);

    // sortedProjects must be in topological (dependency-first) order.
    // progress receives a formatted status line after each project compiles.
    public static IReadOnlyList<ProjectDiagnostics> Compile(
        IReadOnlyList<ProjectAnalysis> sortedProjects,
        string? nugetCacheOverride = null,
        Action<string>? progress = null)
    {
        string nugetCache = nugetCacheOverride ?? NuGetResolver.GetDefaultCacheDir();
        var frameworkRefs = FindFrameworkReferences();

        // Load source generators from the framework's analyzers/ dir
        // (System.Text.Json.SourceGeneration, Microsoft.Interop, etc.).
        // Failing to honor these surfaces as false-positive CS0117/CS0534
        // errors against generator-emitted members.
        var sourceGeneratorBundle = LoadFrameworkSourceGenerators();

        var compilations = new Dictionary<string, CSharpCompilation>(StringComparer.Ordinal);
        var results      = new List<ProjectDiagnostics>();
        var byName       = sortedProjects.ToDictionary(p => p.AssemblyName, StringComparer.Ordinal);

        foreach (ProjectAnalysis project in sortedProjects)
        {
            var parseOptions = BuildParseOptions(project.CompilerSettings);

            var syntaxTrees = project.FileAnalyses
                .Select(fa => TryParse(fa.FilePath, parseOptions))
                .OfType<SyntaxTree>()
                .ToList();

            var refs = new List<MetadataReference>(frameworkRefs);

            foreach (PackageRef pkg in project.PackageReferences)
                foreach (string dll in NuGetResolver.Resolve(pkg, nugetCache, TargetFramework))
                {
                    var mr = TryLoadMetadata(dll);
                    if (mr != null) refs.Add(mr);
                }

            // Add the full transitive closure of project references, not just the
            // direct edges. A CompilationReference exposes only the types declared in
            // that assembly — never its own references — so a type used from a project
            // reached only transitively would raise a false CS0012. A real `dotnet
            // build` flows transitive <ProjectReference>s into the compiler's reference
            // set; mirror that. Topological order guarantees every transitive dep has
            // already been compiled and cached by the time we get here.
            foreach (string dep in TransitiveProjectRefs(project, byName))
                if (compilations.TryGetValue(dep, out CSharpCompilation? depComp))
                    refs.Add(depComp.ToMetadataReference());

            var compilation = CSharpCompilation.Create(
                project.AssemblyName,
                syntaxTrees,
                refs,
                new CSharpCompilationOptions(
                    OutputKind.DynamicallyLinkedLibrary,
                    nullableContextOptions: ParseNullable(project.CompilerSettings.Nullable),
                    allowUnsafe: project.CompilerSettings.AllowUnsafeBlocks ?? false));

            // Run source generators so [JsonSerializable]-emitted members,
            // [LibraryImport]-emitted stubs, etc. are part of the semantic
            // model. Skip when no generators were discoverable — falls
            // back to ungenerated compilation (the historical behavior).
            Compilation effective = compilation;
            if (sourceGeneratorBundle.Generators.Count > 0)
            {
                IDisposable? reflectionScope = null;
                if (sourceGeneratorBundle.LoadContext != null)
                    reflectionScope = sourceGeneratorBundle.LoadContext.EnterContextualReflection();

                try
                {
                    var driver = CSharpGeneratorDriver.Create(
                        generators:      sourceGeneratorBundle.Generators,
                        additionalTexts: null,
                        parseOptions:    parseOptions,
                        optionsProvider: null);
                    try
                    {
                        driver.RunGeneratorsAndUpdateCompilation(
                            compilation, out Compilation augmented, out _);
                        effective = augmented;
                    }
                    catch
                    {
                        // A misbehaving generator shouldn't crash the audit —
                        // fall back to the un-augmented compilation.
                        effective = compilation;
                    }
                }
                finally
                {
                    reflectionScope?.Dispose();
                }
            }

            // The cached compilation is the post-generator one so downstream
            // projects can see generated types from their references.
            compilations[project.AssemblyName] =
                effective as CSharpCompilation ?? compilation;

            var diags = effective.GetDiagnostics();
            var errors = diags.Where(d => d.Severity == DiagnosticSeverity.Error)
                               .Select(ToRecord).ToList().AsReadOnly();
            var warnings = diags.Where(d => d.Severity == DiagnosticSeverity.Warning
                                            && d.Id != "CS1701")
                               .Select(ToRecord).ToList().AsReadOnly();

            // Harvest assembly names of project refs the compilation actually uses.
            // GetUsedAssemblyReferences returns only the references whose types/members
            // appeared in semantic resolution — perfect for detecting truly-stale
            // <ProjectReference> declarations that pass syntactic checks via var/transitive.
            var usedProjectRefs = new HashSet<string>(StringComparer.Ordinal);
            try
            {
                foreach (MetadataReference r in effective.GetUsedAssemblyReferences())
                    if (r is CompilationReference cr && cr.Compilation.AssemblyName is { Length: > 0 } name)
                        usedProjectRefs.Add(name);
            }
            catch
            {
                // GetUsedAssemblyReferences can fail on partial / broken compilations.
                // Leave the set empty — the downstream filter is opportunistic.
            }

            var result = new ProjectDiagnostics(project.AssemblyName, errors, warnings, usedProjectRefs);
            results.Add(result);

            if (progress != null)
            {
                string status = errors.Count > 0
                    ? $"{errors.Count} error{(errors.Count == 1 ? "" : "s")}"
                      + (warnings.Count > 0 ? $"  {warnings.Count} warning{(warnings.Count == 1 ? "" : "s")}" : "")
                    : warnings.Count > 0
                        ? $"{warnings.Count} warning{(warnings.Count == 1 ? "" : "s")}"
                        : "clean";
                progress($"  {project.AssemblyName,-34} {status}");
            }
        }

        return results.AsReadOnly();
    }

    // Walks the project-reference graph to its full transitive closure.
    // Cycle-safe via the visited set; emission order is irrelevant — callers
    // look each name up in the compilation cache. Mirrors MSBuild's default
    // transitive ProjectReference flow.
    private static IEnumerable<string> TransitiveProjectRefs(
        ProjectAnalysis project,
        IReadOnlyDictionary<string, ProjectAnalysis> byName)
    {
        var seen  = new HashSet<string>(StringComparer.Ordinal);
        var stack = new Stack<string>(project.ProjectReferences);
        while (stack.Count > 0)
        {
            string name = stack.Pop();
            if (!seen.Add(name)) continue;
            yield return name;
            if (byName.TryGetValue(name, out ProjectAnalysis? dep))
                foreach (string next in dep.ProjectReferences)
                    stack.Push(next);
        }
    }

    // ── Source generator discovery ────────────────────────────────────────

    /// <summary>
    /// Load source generators that ship alongside the framework (the
    /// <c>analyzers/dotnet/cs/</c> sibling of the ref-assembly dir). Covers
    /// <c>System.Text.Json.SourceGeneration</c> (the <c>[JsonSerializable]</c>
    /// generator), <c>Microsoft.Interop.SourceGeneration</c>
    /// (<c>[LibraryImport]</c>), and any future framework generators. NuGet-
    /// package-shipped generators aren't loaded here yet — extend
    /// <see cref="NuGetResolver"/> if a project pulls one in.
    /// </summary>
    /// <remarks>
    /// Returns an empty list on any failure. The audit falls back to
    /// pre-generator semantic analysis in that case (the historical
    /// behavior); generator-emitted members may appear as false-positive
    /// CS0117/CS0534 errors when this happens.
    /// </remarks>
    internal static GeneratorBundle LoadFrameworkSourceGenerators()
    {
        string? analyzerDir = FindAnalyzerDir();
        if (analyzerDir == null) return new(Array.Empty<ISourceGenerator>(), null);

        var generators = new List<ISourceGenerator>();
        var loadContext = new AnalyzerAssemblyLoadContext(analyzerDir);
        using (loadContext.EnterContextualReflection())
        {
            foreach (string dll in Directory.EnumerateFiles(analyzerDir, "*.dll"))
            {
                Assembly asm;
                try { asm = loadContext.LoadFromAssemblyPath(dll); }
                catch { continue; }

                Type[] types;
                try { types = asm.GetTypes(); }
                catch (ReflectionTypeLoadException ex)
                { types = ex.Types.Where(t => t != null).ToArray()!; }
                catch { continue; }

                foreach (Type t in types)
                {
                    try
                    {
                        if (t.IsAbstract || t.IsInterface) continue;
                        if (t.GetConstructor(Type.EmptyTypes) == null) continue;

                        object? instance = null;
                        try
                        {
                            if (typeof(ISourceGenerator).IsAssignableFrom(t))
                                instance = Activator.CreateInstance(t);
                            else if (typeof(IIncrementalGenerator).IsAssignableFrom(t))
                            {
                                var ig = (IIncrementalGenerator)Activator.CreateInstance(t)!;
                                instance = ig.AsSourceGenerator();
                            }
                        }
                        catch { continue; }

                        if (instance is ISourceGenerator sg)
                            generators.Add(sg);
                    }
                    catch { continue; }
                }
            }
        }

        return new(generators, loadContext);
    }

    internal sealed class AnalyzerAssemblyLoadContext : AssemblyLoadContext
    {
        private readonly string _analyzerDir;

        public AnalyzerAssemblyLoadContext(string analyzerDir)
            : base(isCollectible: false)
        {
            _analyzerDir = analyzerDir;
        }

        protected override Assembly? Load(AssemblyName assemblyName)
        {
            string candidate = Path.Combine(_analyzerDir, assemblyName.Name + ".dll");
            if (File.Exists(candidate))
            {
                try
                {
                    return LoadFromAssemblyPath(candidate);
                }
                catch
                {
                    // If an analyzer assembly has a transitive runtime dependency that is
                    // not in the analyzer folder, allow the default load context to resolve it.
                }
            }

            try
            {
                return Default.LoadFromAssemblyName(assemblyName);
            }
            catch
            {
                return null;
            }
        }
    }

    private static string? FindAnalyzerDir()
    {
        string? packsDir = FindPacksDir();
        if (packsDir == null) return null;

        string refPackDir = Path.Combine(packsDir, "Microsoft.NETCore.App.Ref");
        if (!Directory.Exists(refPackDir)) return null;

        string versionPrefix = TargetFramework.Replace("net", "");
        string? best = Directory.EnumerateDirectories(refPackDir)
            .Where(d => Path.GetFileName(d).StartsWith(versionPrefix))
            .OrderByDescending(d => Path.GetFileName(d))
            .FirstOrDefault();
        if (best == null) return null;

        string analyzerDir = Path.Combine(best, "analyzers", "dotnet", "cs");
        return Directory.Exists(analyzerDir) ? analyzerDir : null;
    }

    // ── Framework resolution ──────────────────────────────────────────────

    public static IReadOnlyList<MetadataReference> FindFrameworkReferences()
    {
        string? refDir = FindRefAssemblyDir();
        string searchDir = refDir ?? Path.GetDirectoryName(typeof(object).Assembly.Location)!;

        return Directory.GetFiles(searchDir, "*.dll")
            .Where(f => !f.EndsWith(".resources.dll", StringComparison.OrdinalIgnoreCase))
            .Select(TryLoadMetadata)
            .OfType<MetadataReference>()
            .ToList()
            .AsReadOnly();
    }

    private static string? FindRefAssemblyDir()
    {
        string? packsDir = FindPacksDir();
        if (packsDir == null) return null;

        string refPackDir = Path.Combine(packsDir, "Microsoft.NETCore.App.Ref");
        if (!Directory.Exists(refPackDir)) return null;

        // e.g. net10.0 → look for 10.0.x directories
        string versionPrefix = TargetFramework.Replace("net", "");
        string? best = Directory.EnumerateDirectories(refPackDir)
            .Where(d => Path.GetFileName(d).StartsWith(versionPrefix))
            .OrderByDescending(d => Path.GetFileName(d))
            .FirstOrDefault();
        if (best == null) return null;

        string refDir = Path.Combine(best, "ref", TargetFramework);
        return Directory.Exists(refDir) ? refDir : null;
    }

    private static string? FindPacksDir()
    {
        string? dotnetRoot = TryGetDotnetRootFromEnvironment();
        if (dotnetRoot != null)
        {
            string packs = Path.Combine(dotnetRoot, "packs");
            if (Directory.Exists(packs)) return packs;
        }

        // Navigate up from typeof(object) location until we find a sibling packs/ dir
        string? dir = Path.GetDirectoryName(typeof(object).Assembly.Location);
        for (int i = 0; i < 5 && dir != null; i++)
        {
            dir = Path.GetDirectoryName(dir);
            if (dir == null) break;
            string packs = Path.Combine(dir, "packs");
            if (Directory.Exists(packs)) return packs;
        }
        return null;
    }

    private static string? TryGetDotnetRootFromEnvironment()
    {
        string? dotnetRoot = Environment.GetEnvironmentVariable("DOTNET_ROOT");
        if (!string.IsNullOrWhiteSpace(dotnetRoot) && Directory.Exists(dotnetRoot))
            return dotnetRoot;

        string? installDir = Environment.GetEnvironmentVariable("DOTNET_INSTALL_DIR");
        if (!string.IsNullOrWhiteSpace(installDir) && Directory.Exists(installDir))
            return installDir;

        return null;
    }

    // ── Parse / compile options ───────────────────────────────────────────

    private static CSharpParseOptions BuildParseOptions(CsprojCompilerSettings settings)
    {
        LanguageVersion langVersion = settings.LangVersion?.ToLowerInvariant() switch
        {
            "preview"  => LanguageVersion.Preview,
            "latest" or null => LanguageVersion.Latest,
            "13" or "13.0"   => LanguageVersion.CSharp13,
            "12" or "12.0"   => LanguageVersion.CSharp12,
            "11" or "11.0"   => LanguageVersion.CSharp11,
            "10" or "10.0"   => LanguageVersion.CSharp10,
            "9"  or "9.0"    => LanguageVersion.CSharp9,
            "8"  or "8.0"    => LanguageVersion.CSharp8,
            _                => LanguageVersion.Latest
        };
        return new CSharpParseOptions(langVersion);
    }

    private static NullableContextOptions ParseNullable(string? value) =>
        value?.ToLowerInvariant() switch
        {
            "enable"      => NullableContextOptions.Enable,
            "warnings"    => NullableContextOptions.Warnings,
            "annotations" => NullableContextOptions.Annotations,
            "disable"     => NullableContextOptions.Disable,
            // Default matches this repo's Directory.Build.props
            _             => NullableContextOptions.Enable
        };

    // ── Helpers ───────────────────────────────────────────────────────────

    private static SyntaxTree? TryParse(string filePath, CSharpParseOptions options)
    {
        try
        {
            return CSharpSyntaxTree.ParseText(File.ReadAllText(filePath), options, filePath);
        }
        catch { return null; }
    }

    private static MetadataReference? TryLoadMetadata(string dllPath)
    {
        // CreateFromFile memory-maps the DLL and holds the OS file lock for the
        // lifetime of the process — blocks any concurrent build/clean of the
        // referenced assemblies (framework refs, NuGet cache, etc). Read the
        // bytes into memory first so the file handle releases immediately.
        try
        {
            byte[] image = File.ReadAllBytes(dllPath);
            return MetadataReference.CreateFromImage(image, filePath: dllPath);
        }
        catch { return null; }
    }

    private static CompilerDiagnostic ToRecord(Diagnostic d)
    {
        FileLinePositionSpan span = d.Location.GetLineSpan();
        return new CompilerDiagnostic(
            d.Id,
            d.Severity.ToString(),
            d.GetMessage(),
            span.IsValid ? span.Path : null,
            span.IsValid ? span.StartLinePosition.Line + 1 : 0,
            span.IsValid ? span.StartLinePosition.Character + 1 : 0);
    }
}
