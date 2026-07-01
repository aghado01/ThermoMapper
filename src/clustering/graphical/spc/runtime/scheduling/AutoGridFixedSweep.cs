using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

using Graphs.Models.Potts;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Builds a <see cref="FixedGridSweepConfig"/> whose temperature grid is auto-bracketed
/// from the graph via <see cref="SpcScheduleHelpers"/> — the "no grid needed" convenience
/// the parked adaptive scheduler offered, minus the signal-driven refinement. The grid is a
/// plain log-spaced span over the estimated bracket; everything downstream is an ordinary
/// fixed-grid run. Defaults are tuned for small induced subgraphs (e.g. Mapper patches).
/// </summary>
/// <remarks>
/// This is the bridge for consumers that previously leaned on the adaptive strategy purely
/// to avoid supplying a temperature grid. It deliberately does <b>not</b> refine — it picks a
/// fixed grid once and runs it. Callers that want to <i>see</i> the curves should supply their
/// own grid via <see cref="FixedGridSweepConfig"/> directly.
/// </remarks>
public static class AutoGridFixedSweep
{
    /// <summary>
    /// Build a fixed-grid config with a log-spaced grid auto-bracketed from
    /// <paramref name="graph"/>. <paramref name="sampler"/> null uses the default
    /// <see cref="PottsModelConfig"/>.
    /// </summary>
    public static FixedGridSweepConfig BuildConfig(
        CsrGraph graph,
        PottsModelConfig? sampler = null,
        int gridSteps = 12,
        double coldOvershoot = 0.05,
        double hotOvershoot = 5.0,
        int burnInCycles = 200,
        int measureCycles = 400,
        int equilibriumBurnIn = 1000,
        int equilibriumCycles = 5000,
        int replicas = 1,
        AccumulationSpec accumulation = default,
        int? baseSeed = null,
        string? checkpointDirectory = null)
    {
        sampler ??= new PottsModelConfig();

        (double lo, double hi) = SpcScheduleHelpers.EstimateBracket(graph, sampler.Q, coldOvershoot, hotOvershoot);
        double[] grid = SpcScheduleHelpers.LogSpaceGrid(lo, hi, gridSteps);

        return new FixedGridSweepConfig
        {
            Temperatures        = grid,
            Replicas            = replicas,
            SweepBudget         = new RunBudget(burnInCycles, measureCycles),
            Sampler             = sampler,
            EquilibriumBudget   = new RunBudget(equilibriumBurnIn, equilibriumCycles),
            Accumulation        = accumulation,
            BaseSeed            = baseSeed,
            CheckpointDirectory = checkpointDirectory,
        };
    }
}
