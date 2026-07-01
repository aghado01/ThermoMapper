using System;
using System.Collections.Generic;

namespace TestHarness.Runner;

internal sealed record HarnessOptions(
    string Project,
    string Filter,
    string Configuration,
    int MaxWorkers,
    IReadOnlyList<string> CommonArgs,
    bool NoRestore,
    bool NoBuild,
    bool ListOnly,
    string? RunDirectory)
{
    public static HarnessOptions Parse(string[] args)
    {
        string? project = null;
        string? fixture = null;
        string? filter = null;
        // Release by default: this is a parallel performance runner, and unoptimized Debug
        // builds push heavy compute facts (e.g. dense eigensolves) past their fixture timeouts.
        string configuration = "Release";
        int maxWorkers = Environment.ProcessorCount;
        var commonArgs = new List<string>();
        bool noRestore = false;
        bool noBuild = false;
        bool listOnly = false;
        string? runDirectory = null;

        for (int i = 0; i < args.Length; i++)
        {
            string arg = args[i];
            bool hasNext = i + 1 < args.Length;

            switch (arg)
            {
                case "--project" when hasNext:      project = args[++i]; break;
                case "--fixture" when hasNext:      fixture = args[++i]; break;
                case "--filter"  when hasNext:      filter  = args[++i]; break;
                case "--configuration" when hasNext: configuration = args[++i]; break;
                case "--max-workers"   when hasNext: maxWorkers = ParseInt("--max-workers", args[++i]); break;
                case "--common-arg"    when hasNext: commonArgs.Add(args[++i]); break;
                case "--run-directory" when hasNext: runDirectory = args[++i]; break;
                case "--no-restore": noRestore = true; break;
                case "--no-build":   noBuild   = true; break;
                case "--list-only":  listOnly  = true; break;
                default:
                    throw new ArgumentException(
                        $"Unknown or incomplete argument: '{arg}'\n" + Usage());
            }
        }

        if (string.IsNullOrWhiteSpace(project))
            throw new ArgumentException("--project is required.\n" + Usage());

        if (fixture is not null && filter is not null)
            throw new ArgumentException("--fixture and --filter are mutually exclusive.");

        string resolvedFilter = fixture is not null
            ? $"FullyQualifiedName~{fixture}"
            : filter ?? string.Empty;

        int clampedWorkers = Math.Max(1, Math.Min(maxWorkers, Environment.ProcessorCount));

        return new HarnessOptions(
            Project: project,
            Filter: resolvedFilter,
            Configuration: configuration,
            MaxWorkers: clampedWorkers,
            CommonArgs: commonArgs,
            NoRestore: noRestore,
            NoBuild: noBuild,
            ListOnly: listOnly,
            RunDirectory: runDirectory);
    }

    private static int ParseInt(string argName, string value)
    {
        if (!int.TryParse(value, out int n))
            throw new ArgumentException($"{argName} must be an integer, got '{value}'.");
        return n;
    }

    private static string Usage() =>
        "Usage: TestHarness.Runner --project <csproj> " +
        "[--fixture <ClassName> | --filter <expression>] " +
        "[--configuration Debug|Release] [--max-workers N] " +
        "[--common-arg <arg>]... [--no-restore] [--no-build] [--list-only] " +
        "[--run-directory <path>]";
}
