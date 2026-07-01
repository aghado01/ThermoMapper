using Graphs.Primitives;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Internal cluster-validation index over a <b>weighted graph</b> partition —
/// the graph-based sibling of <see cref="IInternalClusterEvaluator"/> (which is
/// point/distance-based). Scores a labeling using only the graph, a per-edge
/// weight field, and the labels — no ground truth, no domain types. The edge
/// weight is whatever the producer commits to its currency (SW bond frequency,
/// PKWang closed-form affinity, …), so a single implementation scores every
/// inference method's output identically.
/// </summary>
public interface IGraphPartitionEvaluator
{
    /// <summary>Short stable identifier (CSV column key); overrides the type name.</summary>
    string Name => GetType().Name;

    /// <param name="graph">Symmetric CSR graph (evaluators walk the <c>j &gt; i</c> half).</param>
    /// <param name="edgeWeight">Per-CSR-slot edge weight, parallel to <see cref="CsrGraph.Targets"/>.</param>
    /// <param name="labels">Per-node cluster label, length <see cref="CsrGraph.NodeCount"/>.</param>
    /// <param name="clusterCount">Number of distinct dense labels in <paramref name="labels"/>.</param>
    /// <returns>A scalar quality score; direction is index-specific.</returns>
    double Evaluate(CsrGraph graph, double[] edgeWeight, int[] labels, int clusterCount);
}
