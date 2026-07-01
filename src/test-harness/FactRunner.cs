using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

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

        foreach (string a in new[]
            { "test", project, "--no-build", "--nologo", "--verbosity", "minimal",
              "--configuration", configuration,
              "--filter", $"FullyQualifiedName={factName}" })
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
}
