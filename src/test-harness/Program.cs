using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Threading;
using System.Threading.Tasks;

namespace TestHarness.Runner;

// This class needs to be renamed
internal static class Program
{
    private static readonly JsonSerializerOptions JsonOpts = new() { WriteIndented = true };

    public static async Task<int> Main(string[] args)
    {
        HarnessOptions opts;
        try
        {
            opts = HarnessOptions.Parse(args);
        }
        catch (ArgumentException ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 2;
        }

        string repoRoot = PathUtil.FindRepositoryRoot();

        string projectPath = Path.IsPathRooted(opts.Project)
            ? opts.Project
            : Path.GetFullPath(Path.Combine(repoRoot, opts.Project));

        // Build
        if (!opts.NoBuild)
        {
            Console.WriteLine($"Building {Path.GetFileName(projectPath)}...");
            int buildExit = DotnetSync("build", projectPath, "--nologo", "--verbosity", "minimal", "--configuration", opts.Configuration);
            if (buildExit != 0)
            {
                Console.Error.WriteLine("Build failed.");
                return buildExit;
            }
        }

        // Discover
        IReadOnlyList<string> facts;
        try
        {
            facts = TestDiscovery.DiscoverFacts(projectPath, opts.Filter, opts.Configuration, opts.NoRestore, noBuild: true);
        }
        catch (InvalidOperationException ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }

        if (facts.Count == 0)
        {
            Console.Error.WriteLine("No tests matched the given filter.");
            return 1;
        }

        if (opts.ListOnly)
        {
            Console.WriteLine($"Project\t{projectPath}");
            Console.WriteLine($"Filter\t{opts.Filter}");
            Console.WriteLine($"Workers\t{opts.MaxWorkers}");
            Console.WriteLine($"Facts\t{facts.Count}");
            foreach (string f in facts) Console.WriteLine($"  {f}");
            return 0;
        }

        // Plan
        string runStamp    = DateTime.Now.ToString("yyyyMMdd_HHmmssfff");
        Guid runId         = Guid.NewGuid();
        string runDirectory = opts.RunDirectory ?? PathUtil.BuildRunRoot(repoRoot, DateTime.Now);
        PathUtil.EnsureDirectory(runDirectory);

        var suiteGroups        = facts.GroupBy(PathUtil.SuiteNameFromFact).ToDictionary(g => g.Key, g => (IReadOnlyList<string>)g.ToList());
        var suiteDirs          = new Dictionary<string, string>();
        var suiteManifestPaths = new Dictionary<string, string>();

        foreach (var (suiteName, suiteTests) in suiteGroups)
        {
            string suiteDir      = PathUtil.BuildSuiteDirectory(repoRoot, suiteName, runStamp);
            string manifestPath  = Path.Combine(suiteDir, PathUtil.Sanitize(suiteName) + ".manifest.json");
            PathUtil.EnsureDirectory(suiteDir);
            WriteSuiteManifest(manifestPath, suiteName, projectPath, opts.Configuration, suiteTests, runId, runStamp, runDirectory, DateTime.Now, status: "Planned");
            suiteDirs[suiteName]          = suiteDir;
            suiteManifestPaths[suiteName] = manifestPath;
        }

        // Run
        var results  = new ConcurrentBag<FactResult>();
        int ordinal  = 0;
        using var cts = new CancellationTokenSource();

        Console.CancelKeyPress += (_, e) =>
        {
            e.Cancel = true;
            Console.WriteLine("\nCancelling — waiting for in-flight facts to stop...");
            cts.Cancel();
        };

        DateTime startedAt = DateTime.Now;
        int workers = Math.Min(opts.MaxWorkers, facts.Count);
        Console.WriteLine($"Running {facts.Count} fact(s) across {suiteGroups.Count} suite(s) with {workers} worker(s).");

        try
        {
            await Parallel.ForEachAsync(
                facts,
                new ParallelOptions { MaxDegreeOfParallelism = workers, CancellationToken = cts.Token },
                async (fact, ct) =>
                {
                    string suiteName    = PathUtil.SuiteNameFromFact(fact);
                    string suiteDir     = suiteDirs[suiteName];
                    string manifestPath = suiteManifestPaths[suiteName];

                    int    n        = Interlocked.Increment(ref ordinal);
                    string stem     = PathUtil.Sanitize(fact);
                    string stdout   = Path.Combine(runDirectory, $"{n:D3}__{stem}.stdout.log");
                    string stderr   = Path.Combine(runDirectory, $"{n:D3}__{stem}.stderr.log");

                    Console.WriteLine($"  START  [{n:D3}] {fact}");
                    DateTime factStart = DateTime.Now;

                    int exitCode;
                    try
                    {
                        exitCode = await FactRunner.RunAsync(
                            projectPath, fact, opts.Configuration, opts.CommonArgs,
                            opts.NoRestore, stdout, stderr, suiteDir, manifestPath, ct);
                    }
                    catch (OperationCanceledException)
                    {
                        exitCode = -1;
                    }

                    double elapsed = (DateTime.Now - factStart).TotalSeconds;
                    string status  = exitCode == 0 ? "PASS" : exitCode == -1 ? "CANCEL" : "FAIL";
                    Console.WriteLine($"  {status,-6} [{n:D3}] {fact} ({elapsed:F1}s)");

                    results.Add(new FactResult(fact, suiteName, exitCode, elapsed, stdout, stderr));
                });
        }
        catch (OperationCanceledException)
        {
            // Expected on Ctrl-C; fall through to finalize and summary.
        }

        // Finalize suite manifests
        var orderedResults = results.OrderBy(r => r.FactName, StringComparer.Ordinal).ToList();
        DateTime completedAt = DateTime.Now;

        foreach (var (suiteName, suiteTests) in suiteGroups)
        {
            var suiteResults = orderedResults.Where(r => r.SuiteName == suiteName).ToList();
            bool anyFailed   = suiteResults.Any(r => r.ExitCode != 0);
            WriteSuiteManifest(suiteManifestPaths[suiteName], suiteName, projectPath, opts.Configuration,
                suiteTests, runId, runStamp, runDirectory, startedAt, anyFailed ? "Failed" : "Completed", suiteResults);
        }

        // Summary
        int passed   = orderedResults.Count(r => r.ExitCode == 0);
        int failed   = orderedResults.Count(r => r.ExitCode > 0);
        int canceled = orderedResults.Count(r => r.ExitCode == -1);
        double wallSec = (completedAt - startedAt).TotalSeconds;

        var summary = new
        {
            RunId            = runId,
            RunStamp         = runStamp,
            RunDirectory     = runDirectory,
            Project          = projectPath,
            Configuration    = opts.Configuration,
            Workers          = workers,
            RequestedAtLocal = startedAt,
            CompletedAtLocal = completedAt,
            WallSeconds      = Math.Round(wallSec, 2),
            Total            = orderedResults.Count,
            Passed           = passed,
            Failed           = failed,
            Canceled         = canceled,
            Results          = orderedResults.Select(r => new
            {
                r.FactName,
                r.SuiteName,
                r.ExitCode,
                ElapsedSeconds = Math.Round(r.ElapsedSeconds, 2),
                r.StdoutLogPath,
                r.StderrLogPath,
            }),
        };

        string summaryPath = Path.Combine(runDirectory, "summary.json");
        File.WriteAllText(summaryPath, JsonSerializer.Serialize(summary, JsonOpts));

        Console.WriteLine();
        Console.WriteLine($"Results  {passed} passed  {failed} failed  {canceled} canceled  ({wallSec:F1}s)");
        Console.WriteLine($"RunRoot  {runDirectory}");
        Console.WriteLine($"Summary  {summaryPath}");

        return failed > 0 ? 1 : canceled > 0 ? 2 : 0;
    }

    private static void WriteSuiteManifest(
        string path,
        string suiteName,
        string project,
        string configuration,
        IReadOnlyList<string> tests,
        Guid runId,
        string runStamp,
        string runDirectory,
        DateTime requestedAt,
        string status,
        IReadOnlyList<FactResult>? results = null)
    {
        object manifest = results is null
            ? (object)new
            {
                RunId            = runId,
                RunStamp         = runStamp,
                RunDirectory     = runDirectory,
                SuiteName        = suiteName,
                Project          = project,
                Configuration    = configuration,
                RequestedAtLocal = requestedAt,
                Tests            = tests,
                Status           = status,
            }
            : new
            {
                RunId            = runId,
                RunStamp         = runStamp,
                RunDirectory     = runDirectory,
                SuiteName        = suiteName,
                Project          = project,
                Configuration    = configuration,
                RequestedAtLocal = requestedAt,
                CompletedAtLocal = DateTime.Now,
                Tests            = tests,
                Status           = status,
                Results          = results.Select(r => new
                {
                    r.FactName,
                    r.ExitCode,
                    ElapsedSeconds = Math.Round(r.ElapsedSeconds, 2),
                }),
            };

        File.WriteAllText(path, JsonSerializer.Serialize(manifest, JsonOpts));
    }

    private static int DotnetSync(params string[] args)
    {
        var psi = new ProcessStartInfo("dotnet") { UseShellExecute = false };
        foreach (string a in args) psi.ArgumentList.Add(a);
        using var p = Process.Start(psi)!;
        p.WaitForExit();
        return p.ExitCode;
    }
}

internal sealed record FactResult(
    string FactName,
    string SuiteName,
    int    ExitCode,
    double ElapsedSeconds,
    string StdoutLogPath,
    string StderrLogPath);
