using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text.Json;

namespace Maths.Oracle.Tests;

/// <summary>
/// Bridge to the opt-in R toolchain. Resolves Rscript from the live/user environment or PDenv,
/// runs an oracle script with the <c>r/</c> project as the working dir (so a renv .Rprofile would
/// auto-activate), and parses the JSON it emits. Parity tests gate on <see cref="IsAvailable"/> and
/// skip when R is absent, so the normal <c>dotnet test</c> run is unaffected.
/// </summary>
internal static class ROracle
{
    private const int WindowsAccessViolation = unchecked((int)0xC0000005);
    private const uint SEM_FAILCRITICALERRORS = 0x0001;
    private const uint SEM_NOGPFAULTERRORBOX = 0x0002;
    private const uint SEM_NOOPENFILEERRORBOX = 0x8000;

    private sealed record RToolchain(string RHome, string Rscript, string? RLibs, string? RenvRoot, string? RenvCache);

    private static readonly RToolchain? Toolchain = FindRToolchain();

    public static bool IsAvailable => Toolchain is not null && File.Exists(Toolchain.Rscript);

    /// <summary>
    /// Runs <paramref name="oracleRelPath"/> (relative to r/) as <c>Rscript &lt;oracle&gt; &lt;fixture&gt;
    /// &lt;out.json&gt; [tail…]</c> and returns the parsed root element. The output path is managed here.
    /// </summary>
    public static JsonElement Run(string oracleRelPath, string fixturePath, params string[] tailArgs)
    {
        string outPath = Path.Combine(Path.GetTempPath(), $"oracle_{Guid.NewGuid():N}.json");
        var psi = new ProcessStartInfo
        {
            FileName = Toolchain!.Rscript,
            WorkingDirectory = RDir(),
            RedirectStandardError = true,
            RedirectStandardOutput = true,
            UseShellExecute = false,
        };
        psi.ArgumentList.Add(oracleRelPath);
        psi.ArgumentList.Add(fixturePath);
        psi.ArgumentList.Add(outPath);
        foreach (string a in tailArgs) psi.ArgumentList.Add(a);

        // Mirror r/scripts/r-session.ps1 so the oracle finds the base library without shell pre-activation.
        psi.Environment["R_HOME"] = Toolchain.RHome;
        if (Toolchain.RLibs is not null) psi.Environment["R_LIBS"] = Toolchain.RLibs;
        if (Toolchain.RenvRoot is not null) psi.Environment["RENV_PATHS_ROOT"] = Toolchain.RenvRoot;
        if (Toolchain.RenvCache is not null) psi.Environment["RENV_PATHS_CACHE"] = Toolchain.RenvCache;
        psi.Environment["RENV_CONFIG_SANDBOX_ENABLED"] =
            GetEnv("RENV_CONFIG_SANDBOX_ENABLED") ?? "FALSE";

        uint? previousErrorMode = null;
        if (OperatingSystem.IsWindows())
            previousErrorMode = SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);

        Process? started = null;
        try { started = Process.Start(psi); }
        finally
        {
            if (previousErrorMode is not null)
                SetErrorMode(previousErrorMode.Value);
        }

        using var proc = started!;

        var stdoutTask = proc.StandardOutput.ReadToEndAsync();
        var stderrTask = proc.StandardError.ReadToEndAsync();
        if (!proc.WaitForExit(120_000))
        {
            try { proc.Kill(entireProcessTree: true); }
            catch (InvalidOperationException) { }
            throw new TimeoutException($"R oracle '{oracleRelPath}' timed out after 120 seconds.");
        }

        string stdout = stdoutTask.GetAwaiter().GetResult();
        string stderr = stderrTask.GetAwaiter().GetResult();
        if (!File.Exists(outPath))
            throw new InvalidOperationException(
                $"R oracle '{oracleRelPath}' failed (exit {proc.ExitCode}):\n{stderr}\n{stdout}");

        string json = File.ReadAllText(outPath);
        File.Delete(outPath);
        using var doc = JsonDocument.Parse(json);
        if (proc.ExitCode != 0 && !IsKnownOracleExitAfterJson(oracleRelPath, proc.ExitCode))
            throw new InvalidOperationException(
                $"R oracle '{oracleRelPath}' failed (exit {proc.ExitCode}):\n{stderr}\n{stdout}");

        return doc.RootElement.Clone();
    }

    private static bool IsKnownOracleExitAfterJson(string oracleRelPath, int exitCode)
    {
        // TDAstats / Ripser and Riemann's Grassmann median currently write complete JSON and then
        // exit Rscript with 0xC0000005 on Windows. The JSON contract is complete, so tolerate only
        // these known oracle scripts.
        string normalized = oracleRelPath.Replace('\\', '/');
        return OperatingSystem.IsWindows() &&
               exitCode == WindowsAccessViolation &&
               (normalized.EndsWith("oracles/tda_oracle.R", StringComparison.OrdinalIgnoreCase) ||
                normalized.EndsWith("oracles/mom_oracle.R", StringComparison.OrdinalIgnoreCase));
    }

    private static RToolchain? FindRToolchain()
    {
        foreach (string rHome in CandidateRHomes().Distinct(StringComparer.OrdinalIgnoreCase))
        {
            string? rscript = FindRscript(rHome);
            if (rscript is null) continue;

            string rlang = Directory.GetParent(rHome)?.FullName ?? string.Empty;
            string? rLibs = GetEnv("R_LIBS") ?? ExistingDir(Path.Combine(rlang, "library"));
            string? renvRoot = GetEnv("RENV_PATHS_ROOT") ?? ExistingDir(Path.Combine(rlang, "renv"));
            string? renvCache = GetEnv("RENV_PATHS_CACHE") ??
                                (renvRoot is null ? null : Path.Combine(renvRoot, "cache"));
            return new RToolchain(rHome, rscript, rLibs, renvRoot, renvCache);
        }

        return null;
    }

    private static IEnumerable<string> CandidateRHomes()
    {
        string? rHome = GetEnv("R_HOME");
        if (rHome is not null) yield return rHome;

        string? portableRoot = GetEnv("PORTABLE_ROOT");
        if (portableRoot is not null)
        {
            foreach (string home in RlangHomes(Path.Combine(portableRoot, "rlang"))) yield return home;
        }

        foreach (string root in PdenvRoots())
        {
            foreach (string home in RlangHomes(Path.Combine(root, "rlang"))) yield return home;
        }

        foreach (string path in Paths())
        {
            string candidate = Path.Combine(path, "Rscript.exe");
            if (!File.Exists(candidate)) continue;

            string? home = RHomeFromRscript(candidate);
            if (home is not null) yield return home;
        }
    }

    private static IEnumerable<string> RlangHomes(string rlang)
    {
        if (!Directory.Exists(rlang)) yield break;

        foreach (string dir in Directory.GetDirectories(rlang, "R-*").OrderByDescending(x => x))
            yield return dir;
    }

    private static IEnumerable<string> PdenvRoots()
    {
        foreach (string profile in ProfileRoots())
            yield return Path.Combine(profile, "PDenv");

        foreach (string name in new[] { "CLAUDE_CODE_SHELL", "CLAUDE_CODE_GIT_BASH_PATH" })
        {
            string? value = GetEnv(name);
            if (value is null) continue;

            string normalized = value.Replace('/', Path.DirectorySeparatorChar);
            string marker = Path.DirectorySeparatorChar + "PDenv" + Path.DirectorySeparatorChar;
            int index = normalized.IndexOf(marker, StringComparison.OrdinalIgnoreCase);
            if (index < 0) continue;

            yield return normalized[..(index + marker.Length - 1)];
        }
    }

    private static IEnumerable<string> ProfileRoots()
    {
        foreach (string name in new[] { "USERPROFILE", "HOME" })
        {
            string? value = GetEnv(name);
            if (value is not null) yield return value;
        }

        string profile = Environment.GetFolderPath(Environment.SpecialFolder.UserProfile);
        if (!string.IsNullOrWhiteSpace(profile)) yield return profile;

        foreach (string name in new[] { "CLAUDE_CONFIG_DIR", "CODEX_HOME" })
        {
            string? value = GetEnv(name);
            string? parent = value is null ? null : Directory.GetParent(value)?.FullName;
            if (parent is not null) yield return parent;
        }
    }

    private static IEnumerable<string> Paths()
    {
        string? process = Environment.GetEnvironmentVariable("PATH");
        string? user = Environment.GetEnvironmentVariable("PATH", EnvironmentVariableTarget.User);
        string joined = string.Join(
            Path.PathSeparator.ToString(),
            new[] { process, user }
                .Where(x => !string.IsNullOrWhiteSpace(x))
                .Select(x => Environment.ExpandEnvironmentVariables(x!)));
        foreach (string path in joined.Split(Path.PathSeparator))
            if (!string.IsNullOrWhiteSpace(path)) yield return path;
    }

    private static string? FindRscript(string rHome)
    {
        string[] rels = { Path.Combine("bin", "x64", "Rscript.exe"), Path.Combine("bin", "Rscript.exe") };
        foreach (string rel in rels)
        {
            string candidate = Path.Combine(rHome, rel);
            if (File.Exists(candidate)) return candidate;
        }

        return null;
    }

    private static string? RHomeFromRscript(string rscript)
    {
        DirectoryInfo? dir = Directory.GetParent(rscript);
        if (dir is null) return null;
        if (dir.Name.Equals("x64", StringComparison.OrdinalIgnoreCase) &&
            dir.Parent?.Name.Equals("bin", StringComparison.OrdinalIgnoreCase) == true)
            return dir.Parent.Parent?.FullName;
        if (dir.Name.Equals("bin", StringComparison.OrdinalIgnoreCase))
            return dir.Parent?.FullName;
        return null;
    }

    private static string? GetEnv(string name)
    {
        string? value = Environment.GetEnvironmentVariable(name);
        if (!string.IsNullOrWhiteSpace(value)) return Environment.ExpandEnvironmentVariables(value);

        value = Environment.GetEnvironmentVariable(name, EnvironmentVariableTarget.User);
        return string.IsNullOrWhiteSpace(value) ? null : Environment.ExpandEnvironmentVariables(value);
    }

    private static string? ExistingDir(string path) => Directory.Exists(path) ? path : null;

    [DllImport("kernel32.dll")]
    private static extern uint SetErrorMode(uint uMode);

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
