using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Execution;

/// <summary>
/// Bundle of equilibrium inputs consumed by post-orchestrator analyses
/// (<see cref="Partitions.IPartitionStrategy"/>) — the structural input plus the run's
/// <see cref="Accumulator"/>. (Graph signals no longer read this bundle: they take the
/// model-agnostic per-edge currency directly — see <c>Graphs.Observables.IGraphSignal</c>.)
/// Scalar moments are always populated; the per-edge arrays
/// (<see cref="Accumulator.BondFormedCount"/>/<see cref="Accumulator.SpinAgreementCount"/>) are
/// present only when the run tracked edge observables, <see langword="null"/> otherwise — strategies
/// that need them branch on null, not on a tier flag.
/// </summary>
/// <remarks>
/// <b>What this is not.</b> Not the orchestrator's return type (that additionally carries diagnostics,
/// per-T frames, and sweep-strategy state). This is the minimum bundle a downstream analysis needs to
/// consume the equilibrium output of one run.
/// </remarks>
public sealed record SpcRunResult
{
    public required CsrGraph    Graph       { get; init; }
    public required Accumulator Accumulator { get; init; }
}
