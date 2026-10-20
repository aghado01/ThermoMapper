using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace TestHarness.Runner;

internal static class FactRunner
{
    private const string SharedRunDirectoryVar = "PWSHSPC_SHARED_RUN_DIRECTORY";
    private const string SharedManifestPathVar  = "PWSHSPC_SHARED_MANIFEST_PATH";

    public static async Task<int> RunAsync(
        string project,
        string factName,
        string configuration,
        IReadOnlyList<string> commonArgs,
        bool noRestore,
        string stdoutLogPath,
        string stderrLogPath,
        string runSettingsPath,
        string sharedRunDirectory,
        string manifestPath,
        CancellationToken ct)
    {
        var psi = new ProcessStartInfo("dotnet")
        {
            UseShellExecute = false,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
        };

        WriteCaseRunSettings(runSettingsPath, factName);

        foreach (string a in new[]
            { "test", project, "--no-build", "--nologo", "--verbosity", "minimal",
              "--configuration", configuration,
              "--settings", runSettingsPath })
            psi.ArgumentList.Add(a);

        if (noRestore) psi.ArgumentList.Add("--no-restore");
        foreach (string a in commonArgs) psi.ArgumentList.Add(a);

        // Per-worker env vars injected on child only; parent process never carries these.
        psi.EnvironmentVariables[SharedRunDirectoryVar] = sharedRunDirectory;
        psi.EnvironmentVariables[SharedManifestPathVar]  = manifestPath;

        PathUtil.EnsureDirectory(Path.GetDirectoryName(stdoutLogPath)!);

        using var stdoutWriter = new StreamWriter(stdoutLogPath, append: false, Encoding.UTF8);
        using var stderrWriter = new StreamWriter(stderrLogPath, append: false, Encoding.UTF8);

        using var process = new Process { StartInfo = psi, EnableRaisingEvents = true };

        // Each stream's DataReceived fires sequentially from its own background thread — no lock needed.
        process.OutputDataReceived += (_, e) => { if (e.Data is not null) stdoutWriter.WriteLine(e.Data); };
        process.ErrorDataReceived  += (_, e) => { if (e.Data is not null) stderrWriter.WriteLine(e.Data); };

        process.Start();
        process.BeginOutputReadLine();
        process.BeginErrorReadLine();

        try
        {
            await process.WaitForExitAsync(ct);
        }
        catch (OperationCanceledException)
        {
            process.Kill(entireProcessTree: true);
            // Wait for the process to fully exit so all pending DataReceived events drain before the writers close.
            await process.WaitForExitAsync(CancellationToken.None);
            return -1;
        }

        return process.ExitCode;
    }

    /// <summary>
    /// Selects one case by its exact display name through a runsettings file. A theory display
    /// name carries its arguments (<c>Method(shape: "sine")</c>), so the case is matched on
    /// <c>DisplayName</c>, not <c>FullyQualifiedName</c>, and the filter travels as XML: on the
    /// command line it would pass through MSBuild property parsing and the VSTest task's
    /// re-quoting, which break on quotes, semicolons and backslash-quote pairs.
    /// <c>TreatNoTestsAsError</c> makes a filter that matches nothing fail instead of passing.
    /// </summary>
    private static void WriteCaseRunSettings(string path, string displayName)
    {
        PathUtil.EnsureDirectory(Path.GetDirectoryName(path)!);
        CaseRunSettings(displayName).Save(path);
    }

    /// <summary>The runsettings document that selects the case with this display name.</summary>
    internal static XElement CaseRunSettings(string displayName) =>
        new("RunSettings",
            new XElement("RunConfiguration",
                new XElement("TestCaseFilter", "DisplayName=" + EscapeFilterValue(displayName)),
                new XElement("TreatNoTestsAsError", "true")));

    /// <summary>Backslash-escapes the characters the VSTest filter grammar treats as syntax.</summary>
    private static string EscapeFilterValue(string value)
    {
        var sb = new StringBuilder(value.Length);
        foreach (char c in value)
        {
            if (c is '\\' or '(' or ')' or '&' or '|' or '=' or '!' or '~') sb.Append('\\');
            sb.Append(c);
        }
        return sb.ToString();
    }
}
