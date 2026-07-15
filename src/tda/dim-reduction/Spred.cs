using System.Threading;
using Maths.Geometry.DimReduction;

namespace TDA.DimReduction;

/// <summary>
/// SPRED — shape-preserving linear dimensionality reduction (Yu &amp; You, arXiv:2106.02096). Anneals a
/// k×d orthonormal projection that minimizes the persistent-homology <see cref="PersistenceObjective"/>,
/// i.e. that best preserves the barcode of the ambient cloud in the low-dimensional embedding.
///
/// <para>This is the consumer that wires the two halves the split keeps apart: the
/// <see cref="SubspaceAnnealer"/> engine (Maths.Geometry) and the PH objective (built from TDA.Ph +
/// Graphs). It lives here, above both, because it depends on each.</para>
/// </summary>
public static class Spred
{
    /// <summary>
    /// Reduce <paramref name="data"/> to <paramref name="targetDim"/> dimensions, preserving persistent
    /// homology per <paramref name="objective"/>.
    /// </summary>
    /// <param name="data">Row-major ambient samples.</param>
    /// <param name="targetDim">Embedding dimension k (e.g. 2 or 3).</param>
    /// <param name="objective">The PH objective recipe (graph construction, matched dimensions, …).</param>
    /// <param name="maxIters">Simulated-annealing steps.</param>
    /// <param name="seed">RNG seed for a reproducible anneal; null draws OS entropy.</param>
    /// <param name="cancellationToken">Cancellation observed around objective setup and between annealing steps.</param>
    /// <returns>The best k×d orthonormal projection found.</returns>
    public static double[][] Compute(
        double[][] data, int targetDim, PersistenceObjectiveConfig objective,
        int maxIters = 1000, int? seed = null, CancellationToken cancellationToken = default)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var ph = new PersistenceObjective(data, objective);
        cancellationToken.ThrowIfCancellationRequested();
        return SubspaceAnnealer.Compute(data, targetDim, ph.Evaluate, maxIters, seed, cancellationToken);
    }
}
