using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace Maths.Oracle.Tests;

/// <summary>
/// Bridge to the opt-in portable R toolchain (<c>$PORTABLE_ROOT/rlang</c>). Resolves Rscript,
/// runs an oracle script with the <c>r/</c> project as the working dir (so a renv .Rprofile would
/// auto-activate), and parses the JSON it emits. Parity tests gate on <see cref="IsAvailable"/> and
/// skip when R is absent, so the normal <c>dotnet test</c> run is unaffected.
/// </summary>
internal static class ROracle
{
    private static readonly string Rlang =
        Path.Combine(Environment.GetEnvironmentVariable("PORTABLE_ROOT") ?? string.Empty, "rlang");
    private static readonly string? RHome = FindRHome(Rlang);
    private static readonly string? Rscript = RHome is null ? null : Path.Combine(RHome, "bin", "Rscript.exe");

    public static bool IsAvailable => Rscript is not null && File.Exists(Rscript);

    /// <summary>
    /// Runs <paramref name="oracleRelPath"/> (relative to r/) as <c>Rscript &lt;oracle&gt; &lt;fixture&gt;
    /// &lt;out.json&gt; [tail…]</c> and returns the parsed root element. The output path is managed here.
    /// </summary>
    public static JsonElement Run(string oracleRelPath, string fixturePath, params string[] tailArgs)
    {
        string outPath = Path.Combine(Path.GetTempPath(), $"oracle_{Guid.NewGuid():N}.json");
        var psi = new ProcessStartInfo
        {
            FileName = Rscript!,
            WorkingDirectory = RDir(),
            RedirectStandardError = true,
            RedirectStandardOutput = true,
            UseShellExecute = false,
        };
        psi.ArgumentList.Add(oracleRelPath);
        psi.ArgumentList.Add(fixturePath);
        psi.ArgumentList.Add(outPath);
        foreach (string a in tailArgs) psi.ArgumentList.Add(a);

        // Mirror env-Rlang.ps1 so the oracle finds the base library without shell pre-activation.
        psi.Environment["R_HOME"] = RHome!;
        psi.Environment["R_LIBS"] = Path.Combine(Rlang, "library");
        psi.Environment["RENV_PATHS_ROOT"] = Path.Combine(Rlang, "renv");
        psi.Environment["RENV_PATHS_CACHE"] = Path.Combine(Rlang, "renv", "cache");

        using var proc = Process.Start(psi)!;
        string stderr = proc.StandardError.ReadToEnd();
        proc.WaitForExit();
        if (proc.ExitCode != 0 || !File.Exists(outPath))
            throw new InvalidOperationException(
                $"R oracle '{oracleRelPath}' failed (exit {proc.ExitCode}):\n{stderr}");

        string json = File.ReadAllText(outPath);
        File.Delete(outPath);
        using var doc = JsonDocument.Parse(json);
        return doc.RootElement.Clone();
    }

    private static string? FindRHome(string rlang) =>
        Directory.Exists(rlang) ? Directory.GetDirectories(rlang, "R-*").FirstOrDefault() : null;

    // Walk up from the test bin dir to the repo's r/ (identified by its oracles/ subdir).
    private static string RDir()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir is not null)
        {
            if (Directory.Exists(Path.Combine(dir.FullName, "r", "oracles")))
                return Path.Combine(dir.FullName, "r");
            dir = dir.Parent;
        }
        throw new DirectoryNotFoundException("Could not locate r/oracles above " + AppContext.BaseDirectory);
    }
}
