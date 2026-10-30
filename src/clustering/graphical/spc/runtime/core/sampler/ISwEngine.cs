using System;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Type-erased contract over <see cref="PottsModel{TConfig}"/>. The public
/// <see cref="SwRunner"/> façade holds an <see cref="ISwEngine"/> so that
/// the susceptibility-config choice (made once at construction) does not leak
/// into the public API surface. All members forward to the underlying generic
/// runner via virtual dispatch; this dispatch is one indirection per public
/// call, not per SWCycle, so it is paid at the orchestration boundary, not in
/// the hot path.
/// </summary>
internal interface ISwEngine
{
    void Draw();
    void Run(int cycleCount);
    void BurnIn(int cycles);

    int DrawCount { get; }
    double Temperature { get; }
    int Q { get; }
    int N { get; }

    ReadOnlySpan<int> Spins { get; }
    ReadOnlySpan<long> ClusterSizeHistogram { get; }

    double RunningSumSqClusterSizes { get; }
    double RunningSumSqClusterSizesExcl { get; }
    double RunningSumEnergy { get; }
    double RunningSumEnergySq { get; }
    double RunningSumMag { get; }
    double RunningSumMagSq { get; }

    // The per-draw record, one entry per measured draw (see Accumulator for what each describes).
    ReadOnlySpan<double> EnergyPerDraw { get; }
    ReadOnlySpan<double> MagPerDraw { get; }
    ReadOnlySpan<long> SumSqClusterSizesPerDraw { get; }
    ReadOnlySpan<int> MaxClusterSizePerDraw { get; }

    /// <summary>
    /// Snapshot the full sampler state — scalar sufficient-statistics, the per-draw record + resume
    /// state, plus each optional array (the per-edge counts such as <see cref="Accumulator.BondFormedCount"/>,
    /// the landscapes, the burn-in trace, the exact histogram) when the run tracks it,
    /// <see langword="null"/> otherwise.
    /// </summary>
    Accumulator GetCheckpoint();
    void Restore(Accumulator result);
}
