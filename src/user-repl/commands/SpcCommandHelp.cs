using System;
using System.Collections.Generic;
using System.Linq;

namespace UserRepl.Commands;

internal static class SpcCommandHelp
{
    public static void PrintHelp(IReadOnlyList<string> datasetKinds)
    {
        Console.WriteLine("Usage: userrepl spc --dataset <generator> [options]");
        Console.WriteLine("       userrepl spc --dataset-file <path> [options]");
        Console.WriteLine("       userrepl spc --config <preset.json> [options]");
        Console.WriteLine();
        Console.WriteLine("Presets:");
        Console.WriteLine("  --config <path>               Load a JSON preset. Pre-applied before other flags");
        Console.WriteLine("                                so any CLI flag passed alongside --config overrides");
        Console.WriteLine("                                the corresponding preset entry. Multiple --config");
        Console.WriteLine("                                files stack in argv order (later wins).");
        Console.WriteLine();
        Console.WriteLine("Available synthetic generators:");
        foreach (var kind in datasetKinds) Console.WriteLine($"  {kind}");
        Console.WriteLine();
        Console.WriteLine("Dataset:");
        Console.WriteLine("  --dataset <name>              Synthetic generator (see list above)");
        Console.WriteLine("  --param <name>=<value>        Per-generator parameter; may repeat");
        Console.WriteLine("  --list-generators             Print generators and exit");
        Console.WriteLine("  --generator-schema <name>     Print a generator's parameter schema and exit");
        Console.WriteLine("  --dataset-file <path>         CSV input instead of synthetic data");
        Console.WriteLine("  --label-column <name|idx>     Label column (default: last)");
        Console.WriteLine("  --delimiter <char|tab>        CSV separator (default: ,)");
        Console.WriteLine("  --no-header                   Treat CSV as having no header row");
        Console.WriteLine("  --seed <int>                  RNG seed (default: 42)");
        Console.WriteLine();
        Console.WriteLine("Output:");
        Console.WriteLine("  --base-dir <path>             Run-root base directory (default: artifacts/spc-user)");
        Console.WriteLine("  --run-name <name>             Subdirectory name (default: dataset kind)");
        Console.WriteLine("  --no-guid                     Don't append a GUID to the run folder");
        Console.WriteLine("  --checkpoint-dir <path>       Persist per-task .spcx/.spce here for resume.");
        Console.WriteLine("                                Defaults to <runDir>/checkpoints. Pointing at an");
        Console.WriteLine("                                existing dir resumes any completed work.");
        Console.WriteLine();
        Console.WriteLine("Graph topology:");
        Console.WriteLine("  --topology <knn|epsilon>            Graph topology kind (default: knn)");
        Console.WriteLine("  --filter <or_rule|mutualknn>        Graph filter kind (default: or_rule)");
        Console.WriteLine("  --k <int>                             KNN k (default: 10)");
        Console.WriteLine("  --epsilon <double>                    Epsilon-ball radius (epsilon mode only)");
        Console.WriteLine("  --distance-metric <spec>              euclidean | manhattan | minkowski:p=N");
        Console.WriteLine("                                        | poincare | cosine | hamming.");
        Console.WriteLine("                                        Default: euclidean (inline, no metric obj).");
        Console.WriteLine("  --ensure-connected                    MST-repair the graph (Blatt: typically on)");
        Console.WriteLine("  --lmp                                 Apply Local Mutual Proximity scaling");
        Console.WriteLine();
        Console.WriteLine("Coupling kernel:");
        Console.WriteLine("  --kernel <gaussian|cauchy|laplacian|linear>   Single-kernel coupling (default: gaussian)");
        Console.WriteLine("  --bandwidth <double>                          Kernel bandwidth (0=auto)");
        Console.WriteLine("  --mixture <spec>                              Kernel mixture weights, e.g.");
        Console.WriteLine("                                                'gauss=0.5,cauchy=0.3,laplace=0.2'.");
        Console.WriteLine("                                                When set, overrides --kernel.");
        Console.WriteLine("  --mixture-bandwidth <spec>                    Per-component bandwidths, same format.");
        Console.WriteLine("                                                Omit to auto-estimate per component.");
        Console.WriteLine();
        Console.WriteLine("Schedule:");
        Console.WriteLine("  --schedule <fixed-grid>              Default: fixed-grid (adaptive is parked)");
        Console.WriteLine("  --temperatures <spec>                Required for fixed-grid. Formats:");
        Console.WriteLine("                                         linspace:Tmin,Tmax,N");
        Console.WriteLine("                                         logspace:Tmin,Tmax,N");
        Console.WriteLine("                                         0.01,0.05,0.1,0.5  (explicit list)");
        Console.WriteLine("  --replicas <int>                     Independent replicas per T (default: 1)");
        Console.WriteLine("  --sweep burnin=N,cycles=N            Sweep-probe budget (default: burnin=200,cycles=1000)");
        Console.WriteLine("  --equilibrium burnin=N,cycles=N      Equilibrium-pass budget (default: burnin=1000,cycles=5000)");
        Console.WriteLine("  --accumulation <spec>                Per-edge currencies to collect (default: none)");
        Console.WriteLine("                                         none | affinities | alignments | comembership");
        Console.WriteLine("                                         or comma-joined, e.g. affinities,alignments,comembership");
        Console.WriteLine("  --q <int>                            Potts q (number of colors; default: 20)");
        Console.WriteLine();
        Console.WriteLine("Profile analyzer + cut:");
        Console.WriteLine("  --analyzer <chi-peak>                             Default: chi-peak");
        Console.WriteLine("  --susceptibility <fk-cluster|fk-reduced|magnetization>  Peak-driver χ; default: fk-cluster");
        Console.WriteLine("                                                    (all three always emitted as profile channels)");
        Console.WriteLine("  --partition-strategy <co-membership|spin-agreement|bond-frequency>  Default: co-membership");
        Console.WriteLine("  --peripheral-capture                              Union each node with its max-G neighbor");
        Console.WriteLine("                                                    post-threshold (Domany1999 step 2).");
        Console.WriteLine("                                                    Off by default (strict BWD1995 parity).");
        Console.WriteLine("  --hierarchical-strategy <none|blatt>              Opt-in cross-T partition tree.");
        Console.WriteLine("                                                    'blatt' detects pseudo-transitions on");
        Console.WriteLine("                                                    χ_m, cuts at each stable phase, and");
        Console.WriteLine("                                                    writes partition_hierarchy.json.");
        Console.WriteLine("  --theta <double>                                  Cut threshold (default: 0.5)");
        Console.WriteLine();
        Console.WriteLine("T-stack resolver (post-sweep; resolves the whole stack, not one cut; off by default):");
        Console.WriteLine("  --resolver <none|thermal-eom|hierarchy|lineage>   thermal-eom: merge tree from per-edge G(T)");
        Console.WriteLine("                                                    hierarchy:   Blatt/Domany dendrogram-across-T");
        Console.WriteLine("                                                    lineage:     select lineages by persistence");
        Console.WriteLine("                                                    thermal-eom/hierarchy need --accumulation");
        Console.WriteLine("                                                    comembership,cluster-size-landscape;");
        Console.WriteLine("                                                    lineage needs only comembership.");
        Console.WriteLine("  --min-cluster-size <int>                          Selection floor (default: 1)");
        Console.WriteLine("  --periphery <none|ascend>                         Complete abstains by modal ascent (default: none)");
        Console.WriteLine();
        Console.WriteLine("  --help, -h                    Show this help text");
    }

    public static void PrintGeneratorList()
    {
        Console.WriteLine("Available synthetic generators:");
        foreach (var g in SpcUserSession.ListAvailableSyntheticGenerators())
        {
            Console.WriteLine($"- {g.GeneratorName}");
            if (!string.IsNullOrWhiteSpace(g.Description))
                Console.WriteLine($"    {g.Description}");
        }
    }

    public static void PrintGeneratorSchema(string generatorName)
    {
        var g = SpcUserSession.ListAvailableSyntheticGenerators()
            .FirstOrDefault(d => string.Equals(d.GeneratorName, generatorName, StringComparison.OrdinalIgnoreCase));
        if (g is null)
            throw new ArgumentException($"Unknown generator '{generatorName}'.", nameof(generatorName));

        Console.WriteLine($"Generator: {g.GeneratorName}");
        if (!string.IsNullOrWhiteSpace(g.Description)) Console.WriteLine(g.Description);
        Console.WriteLine("Parameters:");
        foreach (var p in g.Parameters)
        {
            string def = p.DefaultValue is null ? "" : $" [default={p.DefaultValue}]";
            string opt = p.IsOptional ? " (optional)" : "";
            string desc = string.IsNullOrWhiteSpace(p.Description) ? "" : $" - {p.Description}";
            Console.WriteLine($"  {p.Name}: {p.TypeName}{opt}{def}{desc}");
        }
    }
}
