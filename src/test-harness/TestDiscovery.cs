using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace TestHarness.Runner;

internal static class TestDiscovery
{
    private const string AvailableTestsSentinel = "The following Tests are available:";

    public static IReadOnlyList<string> DiscoverFacts(
        string project,
        string filter,
        string configuration,
        bool noRestore,
        bool noBuild)
    {
        var argList = new List<string>
        {
            "test", project,
            "--list-tests",
            "--nologo",
            "--verbosity", "minimal",
            "--configuration", configuration,
        };

        if (noBuild)    argList.Add("--no-build");
        if (noRestore)  argList.Add("--no-restore");

        if (!string.IsNullOrWhiteSpace(filter))
        {
            argList.Add("--filter");
            argList.Add(filter);
        }

        var psi = new ProcessStartInfo("dotnet")
        {
            UseShellExecute = false,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
        };
        foreach (string a in argList) psi.ArgumentList.Add(a);

        var outputLines = new List<string>();
        using var process = new Process { StartInfo = psi };
        process.OutputDataReceived += (_, e) => { if (e.Data is not null) lock (outputLines) outputLines.Add(e.Data); };
        process.ErrorDataReceived  += (_, e) => { };

        process.Start();
        process.BeginOutputReadLine();
        process.BeginErrorReadLine();
        process.WaitForExit();

        if (process.ExitCode != 0)
            throw new InvalidOperationException(
                $"Test discovery exited with code {process.ExitCode}. " +
                $"Check that --project exists and the test binary is built.");

        return ParseFacts(outputLines);
    }

    private static IReadOnlyList<string> ParseFacts(List<string> outputLines)
    {
        bool capture = false;
        var facts = new List<string>();

        foreach (string line in outputLines)
        {
            if (!capture)
            {
                if (line.Contains(AvailableTestsSentinel, StringComparison.Ordinal))
                    capture = true;
                continue;
            }

            string trimmed = line.Trim();
            if (trimmed.Length > 0)
                facts.Add(trimmed);
        }

        return facts.Distinct().OrderBy(f => f, StringComparer.Ordinal).ToList();
    }
}
