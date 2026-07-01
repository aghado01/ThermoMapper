using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Xml.Linq;

namespace RepoAudit;

public static class NuGetResolver
{
    // TFM compatibility chain for net10.0 — checked in order, first hit wins
    private static readonly string[] TfmChain =
    [
        "net10.0", "net9.0", "net8.0", "net7.0", "net6.0", "net5.0",
        "netcoreapp3.1", "netcoreapp3.0", "netcoreapp2.1", "netcoreapp2.0",
        "netstandard2.1", "netstandard2.0"
    ];

    public static string GetDefaultCacheDir()
    {
        if (TryResolveCacheDir(Environment.GetEnvironmentVariable("NUGET_PACKAGES"), out var dir))
            return dir;

        if (TryResolveCacheDir(TryGetCacheDirFromDotnetCliHome(), out dir))
            return dir;

        if (TryResolveCacheDir(TryGetCacheDirFromDotnetRoot(), out dir))
            return dir;

        if (TryResolveCacheDir(Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), ".nuget", "packages"), out dir))
            return dir;

        if (TryResolveCacheDir(TryGetCacheDirFromProcessPath(), out dir))
            return dir;

        return Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), ".nuget", "packages");
    }

    private static bool TryResolveCacheDir(string? candidate, out string resolved)
    {
        if (!string.IsNullOrWhiteSpace(candidate) && Directory.Exists(candidate))
        {
            resolved = Path.GetFullPath(candidate);
            return true;
        }

        resolved = string.Empty;
        return false;
    }

    private static string? TryGetCacheDirFromDotnetCliHome()
    {
        string? cliHome = Environment.GetEnvironmentVariable("DOTNET_CLI_HOME");
        if (string.IsNullOrWhiteSpace(cliHome))
            return null;
        return Path.Combine(cliHome, ".nuget", "packages");
    }

    private static string? TryGetCacheDirFromDotnetRoot()
    {
        string? dotnetRoot = TryGetDotnetRootFromEnvironment();
        if (string.IsNullOrWhiteSpace(dotnetRoot))
            return null;

        var candidate = Path.GetFullPath(Path.Combine(dotnetRoot, "..", "cli-home", ".nuget", "packages"));
        return Directory.Exists(candidate) ? candidate : null;
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

    private static string? TryGetCacheDirFromProcessPath()
    {
        try
        {
            string? mainModule = Process.GetCurrentProcess().MainModule?.FileName;
            var seedPaths = new List<string?>
            {
                mainModule is null ? null : Path.GetDirectoryName(mainModule),
                AppContext.BaseDirectory,
                Environment.CurrentDirectory,
            };

            foreach (string? seed in seedPaths)
            {
                if (string.IsNullOrWhiteSpace(seed))
                    continue;

                var dir = new DirectoryInfo(seed);
                while (dir != null)
                {
                    var candidate1 = Path.Combine(dir.FullName, ".nuget", "packages");
                    if (Directory.Exists(candidate1))
                        return candidate1;

                    var candidate2 = Path.Combine(dir.FullName, "cli-home", ".nuget", "packages");
                    if (Directory.Exists(candidate2))
                        return candidate2;

                    if (dir.Parent != null)
                    {
                        var candidate3 = Path.Combine(dir.Parent.FullName, "dotnet", "cli-home", ".nuget", "packages");
                        if (Directory.Exists(candidate3))
                            return candidate3;
                    }

                    dir = dir.Parent;
                }
            }
        }
        catch
        {
            // best-effort only
        }

        return null;
    }

    // Returns all non-resource DLL paths for the package and its transitive dependencies
    // under the best matching TFM.
    public static IReadOnlyList<string> Resolve(
        PackageRef pkg, string cacheDir, string targetFramework = "net10.0")
    {
        var results = new List<string>();
        ResolvePackage(pkg, cacheDir, targetFramework, results, new HashSet<string>(StringComparer.OrdinalIgnoreCase));
        return results.AsReadOnly();
    }

    private static void ResolvePackage(
        PackageRef pkg,
        string cacheDir,
        string targetFramework,
        List<string> resolved,
        HashSet<string> visited)
    {
        string pkgKey = $"{pkg.Name.ToLowerInvariant()}|{pkg.Version}";
        if (!visited.Add(pkgKey))
            return;

        string pkgDir = Path.Combine(cacheDir, pkg.Name.ToLowerInvariant());
        if (!Directory.Exists(pkgDir))
            return;

        string? versionDir = FindVersionDir(pkgDir, pkg.Version);
        if (versionDir == null)
            return;

        string libDir = Path.Combine(versionDir, "lib");
        if (Directory.Exists(libDir))
        {
            string? tfmDir = FindTfmDir(libDir, targetFramework);
            if (tfmDir != null)
            {
                foreach (string dll in Directory.GetFiles(tfmDir, "*.dll")
                    .Where(f => !f.EndsWith(".resources.dll", StringComparison.OrdinalIgnoreCase)))
                {
                    if (!resolved.Contains(dll, StringComparer.OrdinalIgnoreCase))
                        resolved.Add(dll);
                }
            }
        }

        foreach (PackageRef dep in FindPackageDependencies(versionDir, targetFramework))
            ResolvePackage(dep, cacheDir, targetFramework, resolved, visited);
    }

    private static IReadOnlyList<PackageRef> FindPackageDependencies(string versionDir, string targetFramework)
    {
        string? nuspecPath = Directory.GetFiles(versionDir, "*.nuspec").FirstOrDefault();
        if (nuspecPath is null)
            return Array.Empty<PackageRef>();

        try
        {
            var doc = XDocument.Load(nuspecPath);
            XNamespace ns = doc.Root?.Name.Namespace ?? XNamespace.None;
            var metadata = doc.Root?.Element(ns + "metadata");
            if (metadata is null)
                return Array.Empty<PackageRef>();

            var dependencyGroups = metadata.Element(ns + "dependencies")?.Elements(ns + "group").ToArray();
            IEnumerable<XElement> dependencyElements;

            if (dependencyGroups != null && dependencyGroups.Length > 0)
            {
                XElement? group = FindCompatibleDependencyGroup(dependencyGroups, targetFramework);
                if (group == null)
                    group = dependencyGroups.FirstOrDefault(g => g.Attribute("targetFramework") == null);

                dependencyElements = group?.Elements(ns + "dependency") ?? Enumerable.Empty<XElement>();
            }
            else
            {
                dependencyElements = metadata.Element(ns + "dependencies")?.Elements(ns + "dependency") ?? Enumerable.Empty<XElement>();
            }

            return dependencyElements
                .Select(e => new PackageRef(
                    (string?)e.Attribute("id") ?? string.Empty,
                    NormalizeVersionSpecifier((string?)e.Attribute("version") ?? "*")))
                .Where(pkg => !string.IsNullOrWhiteSpace(pkg.Name))
                .ToArray();
        }
        catch
        {
            return Array.Empty<PackageRef>();
        }
    }

    private static XElement? FindCompatibleDependencyGroup(
        XElement[] groups,
        string targetFramework)
    {
        string normalizedTarget = NormalizeTargetFramework(targetFramework);

        XElement? exactMatch = groups.FirstOrDefault(g =>
            string.Equals(NormalizeTargetFramework((string?)g.Attribute("targetFramework")),
                          normalizedTarget,
                          StringComparison.OrdinalIgnoreCase));
        if (exactMatch != null)
            return exactMatch;

        foreach (string tfm in TfmChain)
        {
            XElement? fallback = groups.FirstOrDefault(g =>
                string.Equals(NormalizeTargetFramework((string?)g.Attribute("targetFramework")),
                              NormalizeTargetFramework(tfm),
                              StringComparison.OrdinalIgnoreCase));
            if (fallback != null)
                return fallback;
        }

        return null;
    }

    private static string NormalizeTargetFramework(string? tfm)
        => tfm?.TrimStart('.')?.Trim() ?? string.Empty;

    private static string NormalizeVersionSpecifier(string version)
    {
        string trimmed = version.Trim();
        if (trimmed.StartsWith("[") || trimmed.StartsWith("(") || trimmed.StartsWith("*"))
        {
            trimmed = trimmed.Trim('[', ']', '(', ')');
            int commaIndex = trimmed.IndexOf(',');
            if (commaIndex >= 0)
                trimmed = trimmed.Substring(0, commaIndex).Trim();
        }
        return string.IsNullOrEmpty(trimmed) ? "*" : trimmed;
    }

    // ── Internals ─────────────────────────────────────────────────────────

    private static string? FindVersionDir(string pkgDir, string versionGlob)
    {
        if (!versionGlob.Contains('*') && !versionGlob.Contains('?'))
        {
            string exact = Path.Combine(pkgDir, versionGlob);
            return Directory.Exists(exact) ? exact : null;
        }

        return Directory.EnumerateDirectories(pkgDir)
            .Where(d => VersionMatchesGlob(Path.GetFileName(d), versionGlob))
            .OrderByDescending(d => Path.GetFileName(d), VersionComparer.Instance)
            .FirstOrDefault();
    }

    private static string? FindTfmDir(string libDir, string targetFramework)
    {
        var candidates = new List<string> { targetFramework };
        foreach (string tfm in TfmChain)
            if (!candidates.Contains(tfm, StringComparer.OrdinalIgnoreCase))
                candidates.Add(tfm);

        foreach (string tfm in candidates)
        {
            string dir = Path.Combine(libDir, tfm);
            if (Directory.Exists(dir)) return dir;
        }
        return null;
    }

    private static bool VersionMatchesGlob(string version, string glob)
    {
        string[] gParts = glob.Split('.');
        string[] vParts = version.Split('.');
        for (int i = 0; i < gParts.Length; i++)
        {
            if (gParts[i] == "*") return true;
            if (i >= vParts.Length) return false;
            if (!gParts[i].Equals(vParts[i], StringComparison.OrdinalIgnoreCase)) return false;
        }
        return true;
    }

    private sealed class VersionComparer : IComparer<string>
    {
        public static readonly VersionComparer Instance = new();

        public int Compare(string? x, string? y)
        {
            if (x == y) return 0;
            int[] px = Parse(x), py = Parse(y);
            int len = Math.Max(px.Length, py.Length);
            for (int i = 0; i < len; i++)
            {
                int va = i < px.Length ? px[i] : 0;
                int vb = i < py.Length ? py[i] : 0;
                if (va != vb) return va.CompareTo(vb);
            }
            return 0;
        }

        private static int[] Parse(string? s) =>
            s?.Split('.').Select(p => int.TryParse(p, out int n) ? n : 0).ToArray()
            ?? Array.Empty<int>();
    }
}
