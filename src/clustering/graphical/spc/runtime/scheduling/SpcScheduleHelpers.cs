using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Models.Potts;
using Graphs.Observables;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Reusable helpers for SPC schedule shapes and schedule-envelope derivation.
/// </summary>
public static class SpcScheduleHelpers
{
    public static double[] LogSpaceGrid(double lo, double hi, int steps)
    {
        if (steps <= 1) return new[] { 0.5 * (lo + hi) };
        if (lo <= 0.0 || hi <= 0.0)
            throw new ArgumentException("LogSpaceGrid requires strictly positive endpoints.");
        if (hi <= lo) (lo, hi) = (hi, lo);

        double logLo = Math.Log(lo);
        double logHi = Math.Log(hi);
        var grid = new double[steps];
        for (int i = 0; i < steps; i++)
            grid[i] = Math.Exp(logLo + i * (logHi - logLo) / (steps - 1));
        return grid;
    }

    public static double Estimate(CsrGraph graph, int q)
    {
        if (graph.NodeCount == 0 || graph.Targets.Length == 0 || q < 2)
            return 1.0;

        EdgeWeightSummary weights = EdgeWeights.Summary(graph);
        double meanWeightedDegree = 2.0 * weights.MeanWeight * weights.EdgeCount / graph.NodeCount;
        double logQ = Math.Log(q);
        if (logQ <= 0.0 || weights.MedianWeight <= 0.0)
            return 1.0;

        return weights.MedianWeight * meanWeightedDegree / logQ;
    }

    public static (double TempMin, double TempMax) EstimateBracket(
        CsrGraph graph,
        int q,
        double coldOvershoot = 0.05,
        double hotOvershoot = 5.0)
    {
        if (coldOvershoot <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(coldOvershoot));
        if (hotOvershoot <= coldOvershoot)
            throw new ArgumentException("hotOvershoot must exceed coldOvershoot.");

        double tc = Estimate(graph, q);
        return (tc * coldOvershoot, tc * hotOvershoot);
    }

    /// <summary>
    /// The temperature bracket a graph declares through its <see cref="ResolvedConstruction"/>; the
    /// one function the REPL's <c>auto</c> and <see cref="AutoGridFixedSweep"/> both call.
    /// Couplings in the classic units (<see cref="CouplingResolution.HasBwdUnits"/>) take BWD's
    /// q-only T_ps(q) bracket. Any other graph takes the data heuristic of <see cref="Estimate"/>
    /// divided by the gauge: the heuristic goes as J² while a gauge g scales J by g, so dividing
    /// keeps the bracket equivariant, <c>J → gJ</c> moving it to <c>T → gT</c>, which is all a gauge is.
    /// <paramref name="construction"/> is null for a graph with no declaration (a Mapper patch);
    /// it takes the heuristic as it stands. <c>Basis</c> names the estimator for the run manifest.
    /// </summary>
    public static (double Lo, double Hi, string Basis) AutoBracket(
        CsrGraph graph,
        int q,
        ResolvedConstruction? construction,
        double coldOvershoot = 0.05,
        double hotOvershoot = 5.0)
    {
        if (construction?.Coupling is { HasBwdUnits: true })
        {
            (double tpsLo, double tpsHi) = BwdPottsCriticalEstimate.Bracket(q);
            return (tpsLo, tpsHi, $"bwd-t_ps(q={q})");
        }

        double gauge = construction?.Coupling is { GaugeFactor: > 0.0 } c ? c.GaugeFactor : 1.0;
        (double lo, double hi) = EstimateBracket(graph, q, coldOvershoot, hotOvershoot);
        string basis = gauge == 1.0
            ? "edge-statistics"
            : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"edge-statistics/gauge={gauge:G6}");
        return (lo / gauge, hi / gauge, basis);
    }

    public static IReadOnlyList<SchedulePartitionRollup> BuildPartitionScheduleRollups(
        CsrGraph graph,
        IEnumerable<Accumulator> observables,
        int[][] groundTruthLevels,
        double theta = 0.5,
        IReadOnlyList<string>? levelNames = null)
    {
        ValidateGraph(graph, nameof(graph));
        if (observables is null) throw new ArgumentNullException(nameof(observables));
        if (groundTruthLevels is null) throw new ArgumentNullException(nameof(groundTruthLevels));

        levelNames ??= new[] { "coarse", "medium", "fine" };

        var grouped = observables
            .GroupBy(o => o.Temperature)
            .OrderByDescending(g => g.Key);

        var rollups = new List<SchedulePartitionRollup>();
        foreach (var group in grouped)
        {
            var replicas = group.ToList();
            foreach (var e in replicas)
            {
                if (e.SpinAgreementCount is null)
                    throw new InvalidOperationException(
                        $"Partition rollups require Tier 1 edge observables, but the T={group.Key} " +
                        "group contains a Standard-tier accumulator (SpinAgreementCount is null).");
            }
            int csrLength = replicas[0].SpinAgreementCount!.Length;
            var pooledSpin = new int[csrLength];
            int pooledCycles = 0;

            foreach (var e in replicas)
            {
                pooledCycles += e.DrawCount;
                for (int i = 0; i < csrLength; i++)
                    pooledSpin[i] += e.SpinAgreementCount![i];
            }

            // Spin-only diagnostic: threshold the pooled spin-agreement rate directly via the
            // shared close (the IPartitionStrategy contract requires an Affinities this rollup
            // doesn't mint; AffinityThreshold.Connect is exactly what ThresholdSpinAgreement calls).
            double invCycles = pooledCycles > 0 ? 1.0 / pooledCycles : 0.0;
            var pooledSpinRate = new double[csrLength];
            for (int i = 0; i < csrLength; i++)
                pooledSpinRate[i] = pooledSpin[i] * invCycles;
            var partition = AffinityThreshold.Connect(graph, pooledSpinRate, theta);
            var purities = new double[groundTruthLevels.Length];
            for (int li = 0; li < groundTruthLevels.Length; li++)
                purities[li] = Purity.Compute(partition.Labels, groundTruthLevels[li]);

            rollups.Add(new SchedulePartitionRollup(
                group.Key,
                replicas.Count,
                partition.Count,
                pooledCycles,
                purities,
                levelNames.Take(groundTruthLevels.Length).ToArray()));
        }

        return rollups;
    }

    private static void ValidateGraph(CsrGraph graph, string paramName)
    {
        if (graph.Targets is null)
            throw new ArgumentException("CSR graph Targets buffer must be initialized.", paramName);
        if (graph.Weights is null)
            throw new ArgumentException("CSR graph Weights buffer must be initialized.", paramName);
        if (graph.RowPointers is null)
            throw new ArgumentException("CSR graph RowPointers buffer must be initialized.", paramName);
        if (graph.NodeCount < 0)
            throw new ArgumentException("CSR graph NodeCount must be non-negative.", paramName);
        if (graph.RowPointers.Length != graph.NodeCount + 1)
            throw new ArgumentException($"CSR graph RowPointers length must equal NodeCount + 1 ({graph.NodeCount + 1}).", paramName);
        if (graph.Targets.Length != graph.Weights.Length)
            throw new ArgumentException("CSR graph Targets and Weights must have the same length.", paramName);
        if (graph.RowPointers.Length > 0 && graph.RowPointers[^1] != graph.Targets.Length)
            throw new ArgumentException("CSR graph RowPointers must end at the number of target entries.", paramName);
        for (int i = 1; i < graph.RowPointers.Length; i++)
        {
            if (graph.RowPointers[i] < graph.RowPointers[i - 1])
                throw new ArgumentException("CSR graph RowPointers must be non-decreasing.", paramName);
        }
    }
}
