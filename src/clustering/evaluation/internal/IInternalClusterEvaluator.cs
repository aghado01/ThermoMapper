namespace Clustering.Evaluation.Internal;

/// <summary>
/// Internal cluster validation index. Measures partition quality using
/// only the data and the partition itself — no ground truth required.
/// </summary>
/// <remarks>
/// <para><b>Input contract.</b> <c>data</c> is row-major
/// (one observation per row, equal-length rows); <c>labels</c>
/// is parallel to <c>data</c>, one cluster label per
/// observation. Labels need not be densely-numbered; implementations
/// densify internally.</para>
///
/// <para><b>Family.</b> Internal indices (Silhouette, Davies-Bouldin,
/// Calinski-Harabasz) live here. Indices that compare two label arrays
/// (Purity, NMI, ARI) live in
/// <see cref="External.IExternalClusterEvaluator"/>; domain-specific
/// indices that consume additional structure (e.g. graph + edge
/// observables) live in their owning domain assembly
/// (<c>Clustering.Graphical.SPC.Evaluators</c>, etc.).</para>
///
/// <para><b>Score direction.</b> Higher-or-lower-is-better is
/// index-specific; each concrete implementation documents its
/// direction.</para>
/// </remarks>
public interface IInternalClusterEvaluator
{
    /// <param name="data">Row-major data matrix; one observation per row.</param>
    /// <param name="labels">Per-observation cluster label.</param>
    /// <returns>A scalar quality score; direction is index-specific.</returns>
    double Evaluate(double[][] data, int[] labels);
}
