using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;

namespace RepoAudit;

public static class GlobCompiler
{
    public static Regex? CompileGlobs(IEnumerable<string> globs)
    {
        var fragments = new List<string>();
        foreach (string g in globs)
        {
            bool anchored = g.TrimEnd('/').Contains('/');
            fragments.Add(TranslateGlob(g, anchored));
        }
        if (fragments.Count == 0) return null;
        return new Regex("(" + string.Join("|", fragments) + ")",
            RegexOptions.Compiled | RegexOptions.IgnoreCase);
    }

    public static string TranslateGlob(string glob, bool anchored)
    {
        bool dirOnly = glob.EndsWith('/');
        if (dirOnly) glob = glob.TrimEnd('/');

        if (glob.StartsWith('/'))
        {
            glob = glob.TrimStart('/');
            anchored = true;
        }

        char[] chars = glob.ToCharArray();
        int len = chars.Length;
        var sb = new StringBuilder(len * 2);
        int i = 0;

        while (i < len)
        {
            char c = chars[i];
            switch (c)
            {
                case '*':
                    if (i + 1 < len && chars[i + 1] == '*')
                    {
                        i += 2;
                        if (i < len && chars[i] == '/') { i++; sb.Append("(.+/)?"); }
                        else sb.Append(".*");
                        continue;
                    }
                    sb.Append("[^/]*");
                    i++;
                    break;
                case '?':
                    sb.Append("[^/]");
                    i++;
                    break;
                case '[':
                    sb.Append('[');
                    i++;
                    if (i < len && chars[i] == '!') { sb.Append('^'); i++; }
                    while (i < len && chars[i] != ']')
                    {
                        if (chars[i] == '\\' && i + 1 < len) { sb.Append('\\'); i++; sb.Append(chars[i]); }
                        else sb.Append(chars[i]);
                        i++;
                    }
                    if (i < len) { sb.Append(']'); i++; }
                    break;
                case '\\':
                    i++;
                    if (i < len) { sb.Append(Regex.Escape(chars[i].ToString())); i++; }
                    break;
                case '.': sb.Append(@"."); i++; break;
                case '+': sb.Append(@"\+"); i++; break;
                case '^': sb.Append(@"\^"); i++; break;
                case '$': sb.Append(@"\$"); i++; break;
                case '{': sb.Append(@"\{"); i++; break;
                case '}': sb.Append(@"\}"); i++; break;
                case '(': sb.Append(@"\("); i++; break;
                case ')': sb.Append(@"\)"); i++; break;
                case '|': sb.Append(@"\|"); i++; break;
                default:  sb.Append(c); i++; break;
            }
        }

        string pattern = sb.ToString();
        pattern = anchored ? "^" + pattern : "(^|/)" + pattern;
        return dirOnly ? pattern + "/" : pattern + "$";
    }

    public static bool GlobSubsumes(string broad, string narrow)
    {
        if (broad == narrow) return true;
        try
        {
            bool anchored = broad.TrimEnd('/').Contains('/');
            var regex = new Regex(TranslateGlob(broad, anchored), RegexOptions.IgnoreCase);
            return regex.IsMatch(narrow.TrimEnd('/'));
        }
        catch { return false; }
    }
}
