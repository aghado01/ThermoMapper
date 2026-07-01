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
    ReadOnlySpan<int> ClusterSizeHistogram { get; }

    double RunningSumSqClusterSizes { get; }
    double RunningSumSqClusterSizesExcl { get; }
    double RunningSumEnergy { get; }
    double RunningSumEnergySq { get; }
    double RunningSumMag { get; }
    double RunningSumMagSq { get; }

    /// <summary>
    /// Snapshot the full sampler state — scalar sufficient-statistics + resume state, plus the
    /// per-edge arrays (<see cref="Accumulator.BondFormedCount"/>/<see cref="Accumulator.SpinAgreementCount"/>)
    /// when this specialization tracks them, <see langword="null"/> otherwise.
    /// </summary>
    Accumulator GetCheckpoint();
    void Restore(Accumulator result);
}
