using System;
using Maths.Information;

namespace Graphs.Observables;

/// <summary>
/// Shannon entropy (nats) of a cluster-size distribution — a degree-0 reduction over a size
/// histogram (count of clusters in each size bin). The graph-intrinsic, model-agnostic kernel for
/// "how concentrated is the partition into a few large clusters vs many small ones."
/// </summary>
/// <remarks>
/// <b>Binning-agnostic</b> — takes any histogram, so the one definition serves both feeds: an SW
/// ensemble's pooled <c>ClusterSizeHistogram</c> (averaged over draws, assembled in profiling) and a
/// future single-<c>Assignment</c> structural size count. The histogram is the additive,
/// commutative sufficient-statistic; this entropy is the nonlinear reduction applied once to the
/// pooled distribution (never per-draw — see <see cref="AffinityEntropy"/>).
/// </remarks>
public static class ClusterSizeEntropy
{
    /// <summary>Entropy of the size histogram in nats.</summary>
    public static double EntropyNats(ReadOnlySpan<int> histogram)
        => Shannon.EntropyNats(histogram);
}
