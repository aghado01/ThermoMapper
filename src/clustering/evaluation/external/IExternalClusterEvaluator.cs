namespace Clustering.Evaluation.External;

/// <summary>
/// External cluster validation index. Measures agreement between a
/// predicted partition and a reference labeling (typically a ground
/// truth, but any reference partition works).
/// </summary>
/// <remarks>
/// <para><b>Input contract.</b> Both label arrays must have equal
/// length; labels need not be densely-numbered. Implementations
/// densify internally as needed.</para>
///
/// <para><b>Family.</b> External indices (Purity, Normalized Mutual
/// Information, Adjusted Rand Index, V-measure, Fowlkes-Mallows) live
/// here. Indices that measure quality without a reference live in
/// <see cref="Internal.IInternalClusterEvaluator"/>.</para>
///
/// <para><b>Symmetry.</b> External indices are not in general
/// symmetric in their two arguments (Purity definitely isn't; NMI and
/// ARI are). Concrete implementations document which argument plays
/// which role.</para>
/// </remarks>
public interface IExternalClusterEvaluator
{
    /// <summary>
    /// Stable identifier used as the dictionary key when the score
    /// lands in a result envelope (e.g.
    /// <c>SpcSessionResult.EvaluatorScores</c>) and as the CSV column
    /// header in downstream exports. Defaults to the runtime type name;
    /// concrete implementations override with a short canonical form
    /// (e.g. <c>"NMI"</c>, <c>"ARI"</c>) when the type name is verbose.
    /// </summary>
    string Name => GetType().Name;

    /// <param name="predictedLabels">Per-observation predicted cluster label.</param>
    /// <param name="referenceLabels">Per-observation reference / ground-truth label.</param>
    /// <returns>A scalar agreement score; direction is index-specific.</returns>
    double Evaluate(int[] predictedLabels, int[] referenceLabels);
}
