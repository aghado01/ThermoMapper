using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Repo.TestHarness;

namespace Spc.BlattAnalyze;

/// <summary>
/// Offline one-off analysis: reads an SpcScheduler checkpoint directory
/// (a flat set of <c>.spcx</c> files), groups by temperature, aggregates
/// any replicas at each T, and emits the canonical SPC diagnostics — FK
/// susceptibility, energy mean, specific heat, and a rough fragmentation
/// proxy (count of unique spin colors) — as CSV.
/// </summary>
/// <remarks>
/// <para><b>Aggregation.</b> Multiple replicas at the same temperature
/// are independent MC chains; their accumulators add directly. Total
/// cycles used as the denominator is the sum across replicas, so chains
/// of different lengths combine correctly.</para>
///
/// <para><b>Inferred-k caveat.</b> The CSV column <c>unique_spins</c> is
/// the number of distinct Potts color values present in the final
/// configuration — a rough fragmentation indicator, not a real cluster
/// count. The proper cluster count comes from a cut policy on the FK
/// bond-cluster decomposition (not yet implemented). Color collisions
/// (two distinct bond clusters drawing the same Potts color) inflate
/// this number; the FK accumulator does not have that bias because it
/// reads UF roots directly.</para>
///
/// <para><b>Formulas.</b>
/// <code>
///   χ_FK(T)  = (Σ_replicas RunningSumSqClusterSizes) / (Σ cycles · N)
///   ⟨E⟩(T)   = (Σ_replicas RunningSumEnergy) / (Σ cycles)
///   ⟨E²⟩(T)  = (Σ_replicas RunningSumEnergySq) / (Σ cycles)
///   C(T)     = (⟨E²⟩ - ⟨E⟩²) / T²
/// </code>
/// </para>
/// </remarks>
internal static class Program
{
    private static int Main(string[] args)
    {
        if (args.Length == 0)
        {
            Console.Error.WriteLine("usage: Spc.BlattAnalyze <checkpoint-directory> [output.csv]");
            return 2;
        }

        string dir = args[0];
        if (!Directory.Exists(dir))
        {
            Console.Error.WriteLine($"directory not found: {dir}");
            return 1;
        }

        var files = Directory.GetFiles(dir, "*.spcx");
        if (files.Length == 0)
        {
            Console.Error.WriteLine($"no .spcx files in {dir}");
            return 1;
        }

        var checkpoints = new List<(Accumulator Main, Accumulator? Edges)>(files.Length);
        foreach (var f in files)
        {
            try
            {
                var main = AccumulatorSerializer.Instance.ReadFromFile(f);

                // Per-edge arrays now travel inside the .spcx accumulator (the .spce
                // sidecar is retired); Edges is the same accumulator when it carries them.
                Accumulator? edges = main.BondFormedCount is not null ? main : null;

                checkpoints.Add((main, edges));
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine($"skip {Path.GetFileName(f)}: {ex.Message}");
            }
        }
        if (checkpoints.Count == 0) return 1;

        int n = checkpoints[0].Main.Spins.Length;
        int sidecarCount = checkpoints.Count(c => c.Edges is not null);
        var byTemperature = checkpoints
            .GroupBy(c => c.Main.Temperature)
            .OrderByDescending(g => g.Key)   // top-down per SPC convention
            .ToList();

        Console.Error.WriteLine(
            $"# {checkpoints.Count} checkpoint(s) ({sidecarCount} with .spce edge sidecars), " +
            $"N = {n}, {byTemperature.Count} unique temperature(s)");

        ArtifactRun run;
        string csvPath;
        if (args.Length > 1)
        {
            csvPath = Path.GetFullPath(args[1]);
            string? outputDirectory = Path.GetDirectoryName(csvPath);
            run = HarnessArtifacts.Attach(
                runKind: "spc",
                suiteName: "Spc.BlattAnalyze",
                runName: "Main",
                runDirectory: string.IsNullOrEmpty(outputDirectory) ? Directory.GetCurrentDirectory() : outputDirectory,
                metadata: new Dictionary<string, object?>
                {
                    ["InputDirectory"] = Path.GetFullPath(dir),
                    ["OutputMode"] = "user-specified-csv",
                });
        }
        else
        {
            run = HarnessArtifacts.Create(
                runKind: "spc",
                suiteName: "Spc.BlattAnalyze",
                runName: "Main",
                metadata: new Dictionary<string, object?>
                {
                    ["InputDirectory"] = Path.GetFullPath(dir),
                    ["OutputMode"] = "default",
                });
            csvPath = SpcOutputPathHelper.GetAnalysisCsvPath(run.RunDirectory);
        }

        using TextWriter writer = csvPath is null
            ? Console.Out
            : new StreamWriter(csvPath);

        writer.WriteLine(
            "T,replicas,total_cycles," +
            "chi_fk,chi_excl,energy_mean,specific_heat,unique_spins," +
            "mean_omega,max_omega,mean_spin_agreement,max_spin_agreement");

        foreach (var group in byTemperature)
        {
            double t = group.Key;
            int replicas = group.Count();

            long   totalCycles        = group.Sum(c => (long)c.Main.DrawCount);
            double sumSqClusters      = group.Sum(c => c.Main.RunningSumSqClusterSizes);
            double sumSqClustersExcl  = group.Sum(c => c.Main.RunningSumSqClusterSizesExcl);
            double sumEnergy          = group.Sum(c => c.Main.RunningSumEnergy);
            double sumEnergySq        = group.Sum(c => c.Main.RunningSumEnergySq);

            double chi          = sumSqClusters     / (totalCycles * (double)n);
            double chiExcl      = sumSqClustersExcl / (totalCycles * (double)n);
            double energyMean   = sumEnergy         /  totalCycles;
            double energySqMean = sumEnergySq       /  totalCycles;
            double specificHeat = (energySqMean - energyMean * energyMean) / (t * t);

            // Rough fragmentation proxy across replicas: mean unique-spin count.
            double uniqueSpins = group.Average(c => (double)c.Main.Spins.Distinct().Count());

            // Per-edge observables (only when .spce sidecars are present).
            // Symmetric CSR stores each undirected edge twice but only the
            // j > i slot is ever incremented; so denominator for "real edges"
            // is csrLength / 2.
            string edgeStats = ",,,";
            var withEdges = group.Where(c => c.Edges is not null).ToList();
            if (withEdges.Count > 0)
            {
                int csrLength = withEdges[0].Edges!.BondFormedCount.Length;
                var bondSum = new long[csrLength];
                var spinSum = new long[csrLength];
                foreach (var c in withEdges)
                {
                    var e = c.Edges!;
                    for (int idx = 0; idx < csrLength; idx++)
                    {
                        bondSum[idx] += e.BondFormedCount[idx];
                        spinSum[idx] += e.SpinAgreementCount[idx];
                    }
                }

                double cyclesD = totalCycles;
                int realEdges  = csrLength / 2;
                double sumOmega = 0, sumSpin = 0, maxOmega = 0, maxSpin = 0;
                for (int idx = 0; idx < csrLength; idx++)
                {
                    double om = bondSum[idx] / cyclesD;
                    double sp = spinSum[idx] / cyclesD;
                    sumOmega += om;
                    sumSpin  += sp;
                    if (om > maxOmega) maxOmega = om;
                    if (sp > maxSpin)  maxSpin  = sp;
                }
                double meanOmega = sumOmega / realEdges;
                double meanSpin  = sumSpin  / realEdges;

                edgeStats = string.Format(CultureInfo.InvariantCulture,
                    "{0:F6},{1:F6},{2:F6},{3:F6}",
                    meanOmega, maxOmega, meanSpin, maxSpin);
            }

            writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                "{0:F5},{1},{2},{3:F6},{4:F6},{5:F4},{6:F4},{7:F2},{8}",
                t, replicas, totalCycles, chi, chiExcl, energyMean, specificHeat, uniqueSpins,
                edgeStats));
        }

        writer.Flush();

        string summaryPath = run.WriteJson(
            "summary.json",
            new
            {
                InputDirectory = Path.GetFullPath(dir),
                CheckpointCount = checkpoints.Count,
                SidecarCount = sidecarCount,
                NodeCount = n,
                UniqueTemperatureCount = byTemperature.Count,
                CsvPath = csvPath,
            });

        Console.WriteLine($"RunRoot\t{run.RunDirectory}");
        Console.WriteLine($"Manifest\t{run.ManifestPath}");
        Console.WriteLine($"Csv\t{csvPath}");
        Console.WriteLine($"Summary\t{summaryPath}");
        return 0;
    }
}
