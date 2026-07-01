using System;
using System.IO;
using System.Text;

namespace TestHarness.Runner;

internal static class PathUtil
{
    public static string FindRepositoryRoot()
    {
        DirectoryInfo? current = new DirectoryInfo(AppContext.BaseDirectory);
        while (current is not null)
        {
            if (File.Exists(Path.Combine(current.FullName, "Directory.Build.props")) &&
                File.Exists(Path.Combine(current.FullName, "changelog.md")))
                return current.FullName;
            current = current.Parent;
        }

        throw new DirectoryNotFoundException(
            "Could not locate repository root — expected Directory.Build.props and changelog.md on the path up from the binary.");
    }

    public static string BuildRunRoot(string repoRoot, DateTime stamp) =>
        Path.Combine(repoRoot, "artifacts", "test-runs", "TestHarness.Runner",
            stamp.ToString("yyyyMMdd_HHmmssfff"));

    public static string BuildSuiteDirectory(string repoRoot, string suiteName, string runStamp) =>
        Path.Combine(repoRoot, "artifacts", "test-runs", Sanitize(suiteName), runStamp);

    public static string SuiteNameFromFact(string fullyQualifiedName)
    {
        ReadOnlySpan<char> span = fullyQualifiedName.AsSpan();
        int last = span.LastIndexOf('.');
        if (last <= 0) return Sanitize(fullyQualifiedName);
        int prev = span[..last].LastIndexOf('.');
        return prev < 0
            ? span[..last].ToString()
            : span[(prev + 1)..last].ToString();
    }

    public static string Sanitize(string value)
    {
        char[] invalid = Path.GetInvalidFileNameChars();
        var sb = new StringBuilder(value.Length);
        bool lastWasSep = false;
        foreach (char c in value)
        {
            bool isSep = Array.IndexOf(invalid, c) >= 0 || char.IsWhiteSpace(c) || c == '.';
            if (isSep)
            {
                if (!lastWasSep) { sb.Append('-'); lastWasSep = true; }
                continue;
            }
            sb.Append(c);
            lastWasSep = false;
        }
        string result = sb.ToString().Trim('-');
        return string.IsNullOrWhiteSpace(result) ? "artifact" : result;
    }

    public static void EnsureDirectory(string path) => Directory.CreateDirectory(path);
}
