using System;
using System.Text;
using System.Text.RegularExpressions;

namespace StorageAudit;

/// <summary>
/// A rule's path pattern, matched against root-relative paths with forward slashes and anchored
/// at the root: <c>*</c> and <c>?</c> stay within one segment, <c>**</c> spans segments, and a
/// trailing <c>/</c> makes the pattern match directories only (without one it matches files only).
/// Matching ignores case, as the Windows file system does. The translation follows repo-audit's
/// <c>GlobCompiler</c>, without its unanchored gitignore forms.
/// </summary>
public sealed class PathGlob
{
    private readonly Regex _regex;

    public PathGlob(string pattern)
    {
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("A glob cannot be empty.", nameof(pattern));

        Pattern = pattern;
        string glob = pattern.Replace('\\', '/').TrimStart('/');
        DirectoryOnly = glob.EndsWith('/');
        glob = glob.TrimEnd('/');
        if (glob.Length == 0)
            throw new ArgumentException($"The glob '{pattern}' names the root itself.", nameof(pattern));

        _regex = new Regex("^" + Translate(glob) + "$",
            RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
    }

    public string Pattern { get; }

    public bool DirectoryOnly { get; }

    /// <summary>Whether a root-relative path (no trailing slash) of the given kind matches.</summary>
    public bool IsMatch(string relativePath, bool isDirectory)
        => isDirectory == DirectoryOnly && _regex.IsMatch(relativePath.TrimEnd('/'));

    private static string Translate(string glob)
    {
        var sb = new StringBuilder(glob.Length * 2);
        int i = 0;
        while (i < glob.Length)
        {
            char c = glob[i];
            if (c == '*' && i + 1 < glob.Length && glob[i + 1] == '*')
            {
                i += 2;
                if (i < glob.Length && glob[i] == '/') { i++; sb.Append("(?:.+/)?"); }
                else sb.Append(".*");
                continue;
            }
            sb.Append(c switch
            {
                '*' => "[^/]*",
                '?' => "[^/]",
                _ => Regex.Escape(c.ToString()),
            });
            i++;
        }
        return sb.ToString();
    }
}
