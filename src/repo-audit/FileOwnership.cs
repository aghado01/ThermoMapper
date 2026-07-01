using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace RepoAudit;

public sealed record ProjectOwnership(
    string CsprojPath,
    string AssemblyName,
    IReadOnlySet<string> OwnedFiles);

public sealed record BrokenInclude(string AssemblyName, string MissingAbsPath);

public sealed record FileOwnershipResult(
    IReadOnlyList<ProjectOwnership> Projects,
    IReadOnlyList<string> OrphanedFiles,
    IReadOnlyList<BrokenInclude> BrokenIncludes);

public static class FileOwnership
{
    public static FileOwnershipResult Build(IEnumerable<CsprojGlobs> allGlobs, IReadOnlyList<CrawlerResult> crawlerResults)
    {
        var allCsFiles = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (CrawlerResult crawlerResult in crawlerResults)
            foreach (CrawlNode node in crawlerResult.Graph.Values)
                foreach (CrawlFile file in node.Files)
                    if (file.AbsolutePath.EndsWith(".cs", StringComparison.OrdinalIgnoreCase))
                        allCsFiles.Add(file.AbsolutePath);

        var projects = new List<ProjectOwnership>();
        var claimed  = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        var brokenIncludes = new List<BrokenInclude>();

        foreach (CsprojGlobs globs in allGlobs)
        {
            string assemblyName = Path.GetFileNameWithoutExtension(globs.CsprojPath);
            var owned = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

            foreach (string pattern in globs.IncludePatterns)
            {
                int wildcardIdx = IndexOfWildcard(pattern);

                if (wildcardIdx < 0)
                {
                    if (allCsFiles.Contains(pattern))
                    {
                        if (!IsExcluded(pattern, globs.ExcludePatterns)) owned.Add(pattern);
                    }
                    else if (pattern.EndsWith(".cs", StringComparison.OrdinalIgnoreCase))
                    {
                        brokenIncludes.Add(new BrokenInclude(assemblyName, pattern));
                    }
                    continue;
                }

                int lastSlash = pattern[..wildcardIdx].LastIndexOf('/');
                string literalBase = lastSlash >= 0 ? pattern[..(lastSlash + 1)] : "";
                string globSuffix  = lastSlash >= 0 ? pattern[(lastSlash + 1)..] : pattern;
                bool anchored = true; // suffix is matched relative to literalBase, so anchor at its start: '*' stays non-recursive, '**' recursive (MSBuild semantics)

                Regex globRegex = new(GlobCompiler.TranslateGlob(globSuffix, anchored),
                    RegexOptions.Compiled | RegexOptions.IgnoreCase);

                foreach (string file in allCsFiles)
                {
                    if (literalBase.Length > 0 &&
                        !file.StartsWith(literalBase, StringComparison.OrdinalIgnoreCase)) continue;

                    string relative = literalBase.Length > 0 ? file[literalBase.Length..] : file;
                    if (!globRegex.IsMatch(relative)) continue;
                    if (IsExcluded(file, globs.ExcludePatterns)) continue;

                    owned.Add(file);
                }
            }

            projects.Add(new ProjectOwnership(globs.CsprojPath, assemblyName, owned));
            foreach (string f in owned) claimed.Add(f);
        }

        var orphaned = new List<string>();
        foreach (string f in allCsFiles)
            if (!claimed.Contains(f)) orphaned.Add(f);

        return new FileOwnershipResult(projects, orphaned, brokenIncludes);
    }

    private static bool IsExcluded(string absoluteFile, IReadOnlyList<string> excludePatterns)
    {
        foreach (string pattern in excludePatterns)
        {
            int wildcardIdx = IndexOfWildcard(pattern);

            if (wildcardIdx < 0)
            {
                if (absoluteFile.Equals(pattern, StringComparison.OrdinalIgnoreCase)) return true;
                continue;
            }

            int lastSlash = pattern[..wildcardIdx].LastIndexOf('/');
            string literalBase = lastSlash >= 0 ? pattern[..(lastSlash + 1)] : "";
            string globSuffix  = lastSlash >= 0 ? pattern[(lastSlash + 1)..] : pattern;

            if (literalBase.Length > 0 &&
                !absoluteFile.StartsWith(literalBase, StringComparison.OrdinalIgnoreCase)) continue;

            string relative = literalBase.Length > 0 ? absoluteFile[literalBase.Length..] : absoluteFile;
            bool anchored = true; // suffix is matched relative to literalBase, so anchor at its start: '*' stays non-recursive, '**' recursive (MSBuild semantics)
            if (Regex.IsMatch(relative, GlobCompiler.TranslateGlob(globSuffix, anchored),
                RegexOptions.IgnoreCase)) return true;
        }
        return false;
    }

    private static int IndexOfWildcard(string pattern)
    {
        for (int i = 0; i < pattern.Length; i++)
            if (pattern[i] == '*' || pattern[i] == '?') return i;
        return -1;
    }
}
