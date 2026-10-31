using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace StorageAudit;

/// <summary>
/// The command line: <c>profile</c> surveys the repository and its temp roots, <c>plan</c> writes the
/// proposed actions for review, and <c>apply</c> executes a reviewed plan (a dry run unless
/// <c>--execute</c>). Translates arguments into the engine's calls; the engine holds no CLI state.
/// </summary>
internal static class Program
{
    private const string Usage = """
        usage:
          storage-audit profile [--repo <dir>] [--rules <file>] [--out <dir>]
          storage-audit plan    [--repo <dir>] [--rules <file>] [--out <dir>]
          storage-audit apply   <plan.json> [--repo <dir>] [--execute] [--archive-to <dir>] [--rules <file>] [--ledger <file>]
        """;

    public static int Main(string[] args)
    {
        if (args.Length == 0 || args[0] is "-h" or "--help" or "help")
        {
            Console.WriteLine(Usage);
            return args.Length == 0 ? 2 : 0;
        }

        string command = args[0];
        string repo = Directory.GetCurrentDirectory();
        string? rules = null, output = null, archiveTo = null, ledger = null, planPath = null;
        bool execute = false, repoGiven = false;

        for (int i = 1; i < args.Length; i++)
        {
            switch (args[i])
            {
                case "--repo" when i + 1 < args.Length: repo = Path.GetFullPath(args[++i]); repoGiven = true; break;
                case "--rules" when i + 1 < args.Length: rules = Path.GetFullPath(args[++i]); break;
                case "--out" when i + 1 < args.Length: output = Path.GetFullPath(args[++i]); break;
                case "--archive-to" when i + 1 < args.Length: archiveTo = Path.GetFullPath(args[++i]); break;
                case "--ledger" when i + 1 < args.Length: ledger = Path.GetFullPath(args[++i]); break;
                case "--execute": execute = true; break;
                default:
                    if (!args[i].StartsWith("--", StringComparison.Ordinal) && planPath is null)
                        planPath = Path.GetFullPath(args[i]);
                    else
                    {
                        Console.Error.WriteLine($"Unknown option: {args[i]}\n{Usage}");
                        return 2;
                    }
                    break;
            }
        }

        try
        {
            return command switch
            {
                "profile" => Survey(repo, rules, output, withPlan: false),
                "plan" => Survey(repo, rules, output, withPlan: true),
                "apply" when planPath is not null => Apply(planPath, repoGiven ? repo : null, rules, ledger, archiveTo, execute),
                _ => UsageError(),
            };
        }
        catch (Exception e) when (e is IOException or InvalidDataException or InvalidOperationException
                                     or System.Text.Json.JsonException or UnauthorizedAccessException)
        {
            Console.Error.WriteLine($"[fatal] {e.Message}");
            return 1;
        }
    }

    private static int UsageError()
    {
        Console.Error.WriteLine(Usage);
        return 2;
    }

    private static string DefaultRules(string repo) => Path.Combine(repo, "presets", "storage-audit", "rules.jsonc");

    private static int Survey(string repo, string? rulesPath, string? output, bool withPlan)
    {
        rulesPath ??= DefaultRules(repo);
        StorageRules rules = StorageRules.Load(rulesPath);
        output ??= Path.Combine(repo, "artifacts", "storage-audit", "runs", DateTime.Now.ToString("yyyyMMdd_HHmmss"));

        StorageProfile profile = StorageProfiler.Run(repo, rules, DateTime.UtcNow, rulesPath);
        ReportWriter.WriteProfile(profile, output);

        Console.WriteLine($"storage-audit {(withPlan ? "plan" : "profile")}");
        Console.WriteLine($"  repo:  {profile.RepoRoot}");
        Console.WriteLine($"  rules: {rulesPath}");
        foreach (Category c in Enum.GetValues<Category>())
        {
            List<StorageItem> items = profile.Items.Where(i => i.Category == c).ToList();
            Console.WriteLine($"  {Format.Name(c),-12} {Format.Size(items.Sum(i => i.Bytes)),10}  ({items.Count} items)");
        }
        foreach (DriveSpace d in profile.Drives)
            Console.WriteLine($"  drive {d.Drive} free {Format.Size(d.FreeBytes)} of {Format.Size(d.TotalBytes)}");

        if (withPlan)
        {
            StoragePlan plan = StoragePlanner.Build(profile, rules, AppContext.BaseDirectory);
            ReportWriter.WritePlan(plan, output);
            foreach (var g in plan.Actions.GroupBy(a => a.Drive))
                Console.WriteLine($"  plan frees on {g.Key}: {Format.Size(g.Sum(a => a.Bytes))} ({g.Count()} actions)");
            Console.WriteLine($"  held back: {plan.Held.Count}");
            Console.WriteLine($"  review: {Path.Combine(output, "plan.md")}");
        }
        Console.WriteLine($"  report: {Path.Combine(output, "profile.md")}");
        return 0;
    }

    private static int Apply(string planPath, string? repo, string? rulesPath, string? ledger, string? archiveTo, bool execute)
    {
        StoragePlan plan = StoragePlan.Read(planPath);
        rulesPath ??= DefaultRules(repo ?? plan.RepoRoot);
        ledger ??= Path.Combine(plan.RepoRoot, "artifacts", "storage-audit", "ledger.jsonl");

        Console.WriteLine($"storage-audit apply ({(execute ? "EXECUTE" : "dry run; add --execute to act")})");
        Console.WriteLine($"  plan:   {planPath}");
        Console.WriteLine($"  rules:  {rulesPath}");
        Console.WriteLine($"  ledger: {ledger}");
        if (archiveTo is not null) Console.WriteLine($"  archive to: {archiveTo}");
        Console.WriteLine();

        IReadOnlyList<ApplyResult> results = PlanApplier.Apply(new ApplyOptions
        {
            PlanPath = planPath,
            RulesPath = rulesPath,
            LedgerPath = ledger,
            Execute = execute,
            ArchiveTo = archiveTo,
            OwnBinaryDirectory = AppContext.BaseDirectory,
            ExpectedRepoRoot = repo,
        }, Console.Out);

        Console.WriteLine();
        foreach (var g in results.GroupBy(r => r.Outcome))
            Console.WriteLine($"  {g.Key,-13} {g.Count(),4}  {Format.Size(g.Sum(r => r.Bytes))}");
        return results.Any(r => r.Outcome is "refused" or "failed" or "partial") ? 3 : 0;
    }
}
