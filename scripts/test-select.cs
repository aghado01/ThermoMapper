// The name graph behind scripts/test-select.ps1, compiled by it with Add-Type. Comments and strings are blanked first, so
// that names and braces are read from code only.
#nullable enable
using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;

public sealed class SourceFacts
{
    /// <summary>The names other files can reach this file by: its top-level types and delegates, and ".Name" for each
    /// extension method it declares.</summary>
    public HashSet<string> Declared { get; } = new(StringComparer.Ordinal);
    /// <summary>The names it reads: identifiers that are not a member access or a call, and ".Name" for each member
    /// access or call, which is how an extension method is read.</summary>
    public HashSet<string> Named { get; } = new(StringComparer.Ordinal);
    /// <summary>"root|Name" for each name read after a dotted chain that starts at an identifier root: a qualified type,
    /// as in Maths.Distributions.Normal.SampleStandard, when root is a namespace's first segment.</summary>
    public HashSet<string> Qualified { get; } = new(StringComparer.Ordinal);
    /// <summary>The first segment of each namespace the file declares.</summary>
    public HashSet<string> NamespaceRoots { get; } = new(StringComparer.Ordinal);
    /// <summary>The types that hold a [Fact] or a [Theory] directly, namespace-qualified.</summary>
    public List<string> TestClasses { get; } = new();
}

public static class SourceGraph
{
    private const RegexOptions O = RegexOptions.Compiled | RegexOptions.CultureInvariant;
    private static readonly Regex Blank = new(
        "\"{3,}[\\s\\S]*?\"{3,}|(?:\\$@|@\\$?)\"(?:\"\"|[^\"])*\"|\"(?:\\\\.|[^\"\\\\\\n])*\"|'(?:\\\\.|[^'\\\\\\n])+'|//[^\\n]*|/\\*[\\s\\S]*?\\*/", O);
    private static readonly Regex Type = new(@"\b(?:class|struct|interface|enum|record)\s+(?:(?:struct|class)\s+)?([A-Za-z_]\w*)", O);
    private static readonly Regex Delegate = new(@"\bdelegate\s+[\w<>\[\]?,.\s]+?\s([A-Za-z_]\w*)\s*(?:<[^<>]*>)?\s*\(", O);
    private static readonly Regex Extension = new(@"\bstatic\b[^;{}()=]*?\b([A-Za-z_]\w*)\s*(?:<[^<>()]*>)?\s*\(\s*this\b", O);
    private static readonly Regex Token = new(@"[A-Za-z_]\w*", O);
    private static readonly Regex UsingTarget = new(@"^\s*(?:global\s+)?using\s+(?:static\s+[\w.]*\.|\w+\s*=\s*[\w.]*\.)(\w+)", O | RegexOptions.Multiline);
    private static readonly Regex Test = new(@"\[(?:Fact|Theory)\b", O);
    private static readonly Regex Namespace = new(@"\bnamespace\s+([\w.]+)", O);
    private static readonly Regex Directive = new(@"^\s*(?:global\s+)?using\s+[\w.]+\s*;|\bnamespace\s+[\w.]+", O | RegexOptions.Multiline);

    private sealed record Span(string Name, int Start, int End, bool HasBody);

    public static SourceFacts Analyze(string text)
    {
        string code = Blank.Replace(text, m => m.Value[0] == '$' || m.Value.StartsWith("@$", StringComparison.Ordinal)
                                               || (m.Index > 0 && text[m.Index - 1] == '$')
                                               ? Holes(m.Value) : new string(' ', m.Length));
        var facts = new SourceFacts();
        List<Span> spans = TypeSpans(code);

        foreach (Span s in spans)
            if (!spans.Exists(o => o != s && o.HasBody && o.Start < s.Start && s.Start < o.End))
                facts.Declared.Add(s.Name);
        foreach (Match m in Delegate.Matches(code)) facts.Declared.Add(m.Groups[1].Value);
        foreach (Match m in Extension.Matches(code)) facts.Declared.Add("." + m.Groups[1].Value);

        foreach (Match m in Namespace.Matches(code)) facts.NamespaceRoots.Add(m.Groups[1].Value.Split('.')[0]);
        foreach (Match m in UsingTarget.Matches(code)) facts.Named.Add(m.Groups[1].Value);

        // Namespace declarations and plain using directives name namespaces, not types: their segments are left out.
        string body = Directive.Replace(code, m => new string(' ', m.Length));
        foreach (Match m in Token.Matches(body))
        {
            int before = m.Index - 1;
            while (before >= 0 && char.IsWhiteSpace(body[before])) before--;
            if (before >= 0 && (char.IsLetterOrDigit(body[before]) || body[before] == '_') && before == m.Index - 1) continue;
            if (before >= 0 && body[before] == '.')
            {
                facts.Named.Add("." + m.Value);
                if (ChainRoot(body, before) is { } root) facts.Qualified.Add(root + "|" + m.Value);
                continue;
            }
            if (IsCall(body, m) && !AfterNew(body, m.Index)) continue;
            facts.Named.Add(m.Value);
        }

        if (Test.IsMatch(code))
        {
            Match ns = Namespace.Match(code);
            string prefix = ns.Success ? ns.Groups[1].Value + "." : "";
            var seen = new HashSet<string>(StringComparer.Ordinal);
            foreach (Match a in Test.Matches(code))
            {
                // The innermost type around the attribute, named as the runner names it: Outer+Inner when nested.
                var around = spans.FindAll(s => s.HasBody && s.Start < a.Index && a.Index < s.End);
                if (around.Count == 0) continue;
                around.Sort((x, y) => x.Start.CompareTo(y.Start));
                string name = string.Join("+", around.ConvertAll(s => s.Name));
                if (seen.Add(name)) facts.TestClasses.Add(prefix + name);
            }
        }
        return facts;
    }

    // An interpolated string's holes kept as code, its text and its hole braces blanked, the length unchanged.
    private static string Holes(string s)
    {
        var chars = s.ToCharArray();
        int depth = 0;
        for (int i = 0; i < chars.Length; i++)
        {
            char c = chars[i];
            if (depth == 0 && (c == '{' || c == '}') && i + 1 < chars.Length && chars[i + 1] == c) { chars[i] = chars[i + 1] = ' '; i++; }
            else if (c == '{') { if (depth++ == 0) chars[i] = ' '; }
            else if (c == '}' && depth > 0) { if (--depth == 0) chars[i] = ' '; }
            else if (depth == 0) chars[i] = ' ';
        }
        return new string(chars);
    }

    // A type's span runs from its keyword to its closing brace, or to its semicolon when it has no body.
    private static List<Span> TypeSpans(string code)
    {
        var close = new Dictionary<int, int>();
        var stack = new Stack<int>();
        for (int i = 0; i < code.Length; i++)
        {
            if (code[i] == '{') stack.Push(i);
            else if (code[i] == '}' && stack.Count > 0) close[stack.Pop()] = i;
        }
        var spans = new List<Span>();
        foreach (Match m in Type.Matches(code))
        {
            int depth = 0, i = m.Index + m.Length;
            for (; i < code.Length; i++)
            {
                char c = code[i];
                if (c == '(' || c == '[') depth++;
                else if (c == ')' || c == ']') depth--;
                else if (depth == 0 && (c == '{' || c == ';')) break;
            }
            if (i < code.Length && code[i] == '{' && close.TryGetValue(i, out int end)) spans.Add(new Span(m.Groups[1].Value, m.Index, end, true));
            else spans.Add(new Span(m.Groups[1].Value, m.Index, i, false));
        }
        return spans;
    }

    // Whether the identifier is called: followed by '(' or by type arguments and then '('.
    private static bool IsCall(string code, Match m)
    {
        int i = m.Index + m.Length;
        while (i < code.Length && char.IsWhiteSpace(code[i])) i++;
        if (i < code.Length && code[i] == '<')
        {
            int depth = 0;
            for (; i < code.Length; i++)
            {
                char c = code[i];
                if (c == '<') depth++;
                else if (c == '>') { depth--; if (depth == 0) { i++; break; } }
                else if (c == ';' || c == '{' || c == '}' || c == '(' || c == ')' || c == '=') return false;
            }
            while (i < code.Length && char.IsWhiteSpace(code[i])) i++;
        }
        return i < code.Length && code[i] == '(';
    }

    // The first identifier of the dotted chain whose last dot is at <paramref name="dot"/>, or null when the chain starts
    // at anything else (a call's result, an indexer, a literal).
    private static string? ChainRoot(string code, int dot)
    {
        string? root = null;
        while (true)
        {
            int end = dot - 1;
            while (end >= 0 && char.IsWhiteSpace(code[end])) end--;
            int start = end;
            while (start >= 0 && (char.IsLetterOrDigit(code[start]) || code[start] == '_')) start--;
            if (start == end) return root is null ? null : root;
            root = code.Substring(start + 1, end - start);
            int prior = start;
            while (prior >= 0 && char.IsWhiteSpace(code[prior])) prior--;
            if (prior < 0 || code[prior] != '.') return char.IsLetter(root[0]) || root[0] == '_' ? root : null;
            dot = prior;
        }
    }

    private static bool AfterNew(string code, int index)
    {
        int i = index - 1;
        while (i >= 0 && char.IsWhiteSpace(code[i])) i--;
        return i >= 2 && code[i] == 'w' && code[i - 1] == 'e' && code[i - 2] == 'n' && (i == 2 || !char.IsLetterOrDigit(code[i - 3]));
    }

    /// <summary>The files a change reaches: each changed file, then every file that names a name an affected file
    /// declares, to the fixed point. <paramref name="why"/> takes, for each affected file, the name it was reached by
    /// (or "changed"), and for each name, ":" + name → the file that declared it.</summary>
    public static HashSet<string> Affected(IDictionary<string, SourceFacts> files, IEnumerable<string> changed,
                                           IEnumerable<string> extraNames, IDictionary<string, string> why)
    {
        var roots = new HashSet<string>(StringComparer.Ordinal);
        foreach (SourceFacts f in files.Values) roots.UnionWith(f.NamespaceRoots);
        var index = new Dictionary<string, List<string>>(StringComparer.Ordinal);
        void Index(string n, string file)
        {
            if (!index.TryGetValue(n, out List<string>? list)) index[n] = list = new List<string>();
            if (list.Count == 0 || list[^1] != file) list.Add(file);
        }
        foreach (KeyValuePair<string, SourceFacts> f in files)
        {
            foreach (string n in f.Value.Named) Index(n, f.Key);
            foreach (string q in f.Value.Qualified)
            {
                int bar = q.IndexOf('|');
                if (roots.Contains(q[..bar]) && !f.Value.Named.Contains(q[(bar + 1)..])) Index(q[(bar + 1)..], f.Key);
            }
        }
        var affected = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        var names = new HashSet<string>(StringComparer.Ordinal);
        var queue = new Queue<string>();
        void Reach(string file, string declarer)
        {
            foreach (string d in files[file].Declared)
                if (names.Add(d)) { queue.Enqueue(d); why[":" + d] = declarer; }
        }
        foreach (string c in changed)
            if (files.ContainsKey(c) && affected.Add(c)) { why[c] = "changed"; Reach(c, c); }
        foreach (string e in extraNames)
            if (names.Add(e)) { queue.Enqueue(e); why[":" + e] = "(deleted)"; }
        while (queue.Count > 0)
        {
            string n = queue.Dequeue();
            if (!index.TryGetValue(n, out List<string>? readers)) continue;
            foreach (string file in readers)
                if (affected.Add(file)) { why[file] = n; Reach(file, file); }
        }
        return affected;
    }
}
