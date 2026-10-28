using Graphs.Primitives;

namespace Clustering.Evaluation.External;

/// <summary>
/// External cluster-validation index over a <b>graph</b>: agreement between a
/// predicted partition and a reference labeling, counted on the graph's edges
/// rather than over all point pairs — the graph-restricted sibling of
/// <see cref="IExternalClusterEvaluator"/>, as
/// <see cref="Internal.IGraphPartitionEvaluator"/> is of the point-based internal
/// indices. On a sparse coupling graph it reads where two partitions draw their
/// boundaries; on the complete graph it reduces to the all-pairs index.
/// </summary>
/// <remarks>
/// <para><b>Input contract.</b> Both label arrays have length
/// <see cref="CsrGraph.NodeCount"/>; labels need not be densely numbered.
/// Edges with an endpoint whose predicted label is
/// <see cref="Clustering.Primitives.Assignment.Unassigned"/> are not scored,
/// matching the assigned-subset rule of the all-pairs indices.</para>
/// </remarks>
public interface IGraphExternalClusterEvaluator
{
    /// <summary>Short stable identifier (CSV column key); overrides the type name.</summary>
    string Name => GetType().Name;

    /// <param name="graph">Symmetric CSR graph (evaluators walk the <c>j &gt; i</c> half).</param>
    /// <param name="predictedLabels">Per-node predicted cluster label.</param>
    /// <param name="referenceLabels">Per-node reference label.</param>
    /// <returns>A scalar agreement score; direction is index-specific.</returns>
    double Evaluate(CsrGraph graph, int[] predictedLabels, int[] referenceLabels);
}
