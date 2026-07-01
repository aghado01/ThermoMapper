using System;
using System.Linq;
using UserRepl.Commands;

namespace UserRepl;

/// <summary>
/// Top-level argv router for the UserRepl CLI. Dispatches based on the
/// first positional argument to one of the per-algorithm or per-action
/// subcommands.
/// </summary>
/// <remarks>
/// <para>Subcommands are kept self-contained — each owns its own argv
/// parsing, defaults, and help text. The router's job is only to peel
/// the leading verb off and forward the rest. This keeps adding a new
/// algorithm cheap (one new command class + one switch arm) and avoids
/// cross-coupling between unrelated CLI surfaces.</para>
/// </remarks>
public static class SubcommandRouter
{
    public static int Run(string[] args)
    {
        if (args.Length == 0 || IsHelpFlag(args[0]))
        {
            PrintRootHelp();
            return args.Length == 0 ? 1 : 0;
        }

        string subcommand = args[0];
        string[] rest = args.Skip(1).ToArray();

        return subcommand.ToLowerInvariant() switch
        {
            "spc"          => SpcCommand.Run(rest),
            "hdbscan"      => HdbscanCommand.Run(rest),
            "extract"      => ExtractCommand.Run(rest),
            "graph-health" => GraphHealthCommand.Run(rest),
            _ => UnknownSubcommand(subcommand),
        };
    }

    private static bool IsHelpFlag(string s)
        => s is "--help" or "-h" or "help";

    private static int UnknownSubcommand(string s)
    {
        Console.Error.WriteLine($"Unknown subcommand '{s}'.");
        PrintRootHelp();
        return 1;
    }

    private static void PrintRootHelp()
    {
        Console.WriteLine("Usage: userrepl <subcommand> [options]");
        Console.WriteLine();
        Console.WriteLine("Subcommands:");
        Console.WriteLine("  spc       Run an SPC clustering session (adaptive or fixed-grid sweep)");
        Console.WriteLine("  hdbscan   Run HDBSCAN on the same dataset shapes the SPC subcommand accepts");
        Console.WriteLine("  extract   Rebuild CSVs from an existing checkpoint directory without re-running the sampler");
        Console.WriteLine("  graph-health  Re-evaluate and refresh <run-dir>/graph_health.json from manifest.json");
        Console.WriteLine();
        Console.WriteLine("Use 'userrepl <subcommand> --help' for subcommand-specific options.");
    }
}
