using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Linq;

namespace RepoAudit;

public sealed record PackageRef(string Name, string Version);

public sealed record CsprojCompilerSettings(
    string? LangVersion,
    string? Nullable,
    bool? AllowUnsafeBlocks);

public sealed record CsprojGlobs(
    string CsprojPath,
    string ProjectDir,
    IReadOnlyList<string> IncludePatterns,
    IReadOnlyList<string> ExcludePatterns,
    IReadOnlyList<string> ProjectReferences,
    IReadOnlyList<PackageRef> PackageReferences,
    CsprojCompilerSettings CompilerSettings);

public static class CsprojParser
{
    private static readonly string[] SdkImplicitIncludes = ["**/*.cs"];
    private static readonly string[] SdkImplicitExcludes = ["obj/**", "bin/**"];

    public static CsprojGlobs Parse(string csprojPath)
    {
        string absPath = Path.GetFullPath(csprojPath);
        string projectDir = Path.GetDirectoryName(absPath)!.Replace('\\', '/');
        if (!projectDir.EndsWith('/')) projectDir += '/';

        XDocument doc = XDocument.Load(absPath);
        XNamespace ns = doc.Root?.Name.Namespace ?? XNamespace.None;

        bool enableDefaultCompile = true;
        XElement? defaultCompileEl = doc.Descendants(ns + "EnableDefaultCompileItems").FirstOrDefault();
        if (defaultCompileEl != null &&
            defaultCompileEl.Value.Equals("false", StringComparison.OrdinalIgnoreCase))
            enableDefaultCompile = false;

        var includes = new List<string>();
        var excludes = new List<string>();

        if (enableDefaultCompile)
        {
            foreach (string p in SdkImplicitIncludes)
                includes.Add(ResolvePattern(p, projectDir));
            foreach (string p in SdkImplicitExcludes)
                excludes.Add(ResolvePattern(p, projectDir));
        }

        foreach (XElement item in doc.Descendants(ns + "Compile"))
        {
            // Include/Exclude/Remove are MSBuild item-list attributes: each may carry
            // multiple ';'-separated patterns. Split so every pattern is honored
            // (a single-pattern attribute splits to a one-element list — unchanged).
            foreach (string inc in SplitItemList((string?)item.Attribute("Include")))
                includes.Add(ResolvePattern(inc, projectDir));

            foreach (string exc in SplitItemList((string?)item.Attribute("Exclude")))
                excludes.Add(ResolvePattern(exc, projectDir));

            foreach (string rem in SplitItemList((string?)item.Attribute("Remove")))
                excludes.Add(ResolvePattern(rem, projectDir));
        }

        // Fold in <Compile> items injected by the nearest Directory.Build.props and
        // Directory.Build.targets (MSBuild imports the nearest of each, walking up).
        // These carry shared compile-links the bare .csproj never names — e.g. the
        // test-harness HarnessArtifacts.cs that tests/Directory.Build.targets links into
        // every project under tests/. Without this, such files read as false orphans and
        // their owning projects miss source they actually compile.
        foreach (string dbFile in FindNearestDirectoryBuildFiles(projectDir))
        {
            XDocument dbDoc;
            try { dbDoc = XDocument.Load(dbFile); }
            catch { continue; }

            XNamespace dbNs = dbDoc.Root?.Name.Namespace ?? XNamespace.None;
            string dbDir = Path.GetDirectoryName(dbFile)!.Replace('\\', '/').TrimEnd('/') + "/";

            foreach (XElement item in dbDoc.Descendants(dbNs + "Compile"))
            foreach (string inc in SplitItemList((string?)item.Attribute("Include")))
            {
                string? resolved = ResolveThisFileDirMacro(inc, dbDir);
                if (resolved is null) continue;                    // unresolved $(...) macro
                // Concrete includes here are conditionally injected (Condition="Exists(...)");
                // only claim them when present so we never fabricate a broken-include.
                if (IndexOfWildcard(resolved) < 0 && !File.Exists(resolved)) continue;
                includes.Add(ResolvePattern(resolved, projectDir));
            }
        }

        var projectRefs = new List<string>();
        foreach (XElement item in doc.Descendants(ns + "ProjectReference"))
        {
            string? inc = (string?)item.Attribute("Include");
            if (inc == null) continue;
            projectRefs.Add(Path.GetFileNameWithoutExtension(inc.Replace('\\', '/')));
        }

        var packageRefs = new List<PackageRef>();
        foreach (XElement item in doc.Descendants(ns + "PackageReference"))
        {
            string? name    = (string?)item.Attribute("Include");
            string? version = (string?)item.Attribute("Version")
                           ?? (string?)item.Element(ns + "Version");
            if (name != null)
                packageRefs.Add(new PackageRef(name, version ?? ""));
        }

        string? langVersion = (string?)doc.Descendants(ns + "LangVersion").FirstOrDefault();
        string? nullable    = (string?)doc.Descendants(ns + "Nullable").FirstOrDefault();
        bool? allowUnsafe = null;
        string? unsafeValue = (string?)doc.Descendants(ns + "AllowUnsafeBlocks").FirstOrDefault();
        if (!string.IsNullOrWhiteSpace(unsafeValue) && bool.TryParse(unsafeValue, out bool parsedUnsafe))
            allowUnsafe = parsedUnsafe;
        var compilerSettings = new CsprojCompilerSettings(langVersion, nullable, allowUnsafe);

        return new CsprojGlobs(absPath, projectDir,
            includes.AsReadOnly(), excludes.AsReadOnly(),
            projectRefs.AsReadOnly(), packageRefs.AsReadOnly(),
            compilerSettings);
    }

    // Resolves a csproj-relative pattern to an absolute-base + glob-suffix form.
    // Literal prefix is resolved via Path.GetFullPath; wildcard suffix is preserved as-is.
    private static string ResolvePattern(string pattern, string projectDir)
    {
        int wildcardIdx = IndexOfWildcard(pattern);

        if (wildcardIdx < 0)
            return Path.GetFullPath(Path.Combine(projectDir, pattern)).Replace('\\', '/');

        int lastSlash = pattern[..wildcardIdx].LastIndexOf('/');
        string literalPart  = lastSlash >= 0 ? pattern[..(lastSlash + 1)] : "";
        string globSuffix   = lastSlash >= 0 ? pattern[(lastSlash + 1)..] : pattern;

        string resolvedBase = literalPart.Length > 0
            ? Path.GetFullPath(Path.Combine(projectDir, literalPart)).Replace('\\', '/').TrimEnd('/') + "/"
            : projectDir;

        return resolvedBase + globSuffix;
    }

    private static int IndexOfWildcard(string pattern)
    {
        for (int i = 0; i < pattern.Length; i++)
            if (pattern[i] == '*' || pattern[i] == '?') return i;
        return -1;
    }

    // Walks up from projectDir returning the nearest Directory.Build.props and the nearest
    // Directory.Build.targets (each independently — MSBuild imports only the closest of
    // each, not the whole ancestor chain).
    private static IEnumerable<string> FindNearestDirectoryBuildFiles(string projectDir)
    {
        bool foundProps = false, foundTargets = false;
        DirectoryInfo? dir = new(projectDir.TrimEnd('/'));
        while (dir is not null && (!foundProps || !foundTargets))
        {
            if (!foundProps)
            {
                string p = Path.Combine(dir.FullName, "Directory.Build.props");
                if (File.Exists(p)) { foundProps = true; yield return p; }
            }
            if (!foundTargets)
            {
                string t = Path.Combine(dir.FullName, "Directory.Build.targets");
                if (File.Exists(t)) { foundTargets = true; yield return t; }
            }
            dir = dir.Parent;
        }
    }

    // Substitutes $(MSBuildThisFileDirectory) with the directory of the importing
    // Directory.Build.* file. Returns null if any other unresolved $(...) macro remains,
    // so a half-resolved path is never emitted into ownership matching.
    private static string? ResolveThisFileDirMacro(string pattern, string thisFileDir)
    {
        string resolved = System.Text.RegularExpressions.Regex
            .Replace(pattern, @"\$\(MSBuildThisFileDirectory\)", thisFileDir,
                System.Text.RegularExpressions.RegexOptions.IgnoreCase)
            .Replace('\\', '/');
        return resolved.Contains("$(") ? null : resolved;
    }

    // Splits a ';'-separated MSBuild item-list attribute into normalized
    // (forward-slash) patterns. Empty/whitespace entries are dropped.
    private static IEnumerable<string> SplitItemList(string? attribute)
    {
        if (string.IsNullOrWhiteSpace(attribute))
            yield break;

        foreach (string part in attribute.Split(
            ';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
            yield return part.Replace('\\', '/');
    }
}
