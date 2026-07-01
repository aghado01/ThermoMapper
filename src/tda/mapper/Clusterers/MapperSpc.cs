using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Primitives;
using TDA.Mapper;

using Graphs.Models.Potts;

namespace TDA.Mapper.Clusterers;

/// <summary>
/// Graph-input MAPPER clusterer that runs an SPC fixed-grid temperature sweep on each
/// preimage's induced subgraph and produces cluster labels via an SPC partitioning strategy
/// keyed on equilibrium spin agreement.
/// </summary>
/// <remarks>
/// <para>This is the clusterer that demonstrates the SPC-MAPPER thesis: the same nerve
/// construction that <see cref="ConnectedComponentsClusterer"/> produces with pure topology
/// can be produced with thermodynamic clustering that simultaneously surfaces a hierarchy of
/// substructures across temperature. Drop this in alongside the connected-components clusterer
/// to cross-check at the methodological level.</para>
///
/// <para><b>Grid.</b> Each patch needs a temperature grid; rather than burden the caller with
/// one per patch, the grid is auto-bracketed from the patch subgraph via
/// <see cref="AutoGridFixedSweep"/> (a log-spaced span over the estimated bracket). This is the
/// fixed-grid replacement for the parked adaptive scheduler's coarse→dense pipeline: it picks a
/// grid and runs it, without signal-driven refinement.</para>
///
/// <para><b>Index translation.</b> The <see cref="IGraphClusterer"/> contract expects cluster
/// labels in <em>preimage-index order</em>: <c>labels[i]</c> corresponds to
/// <c>preimageIndices[i]</c>. The induced subgraph is densely renumbered by
/// <see cref="CsrGraph.InducedSubgraph"/>, so the cut policy returns labels in subgraph-local
/// order — this clusterer translates back using the <c>oldToNew</c> output of the subgraph
/// extraction.</para>
/// </remarks>
public sealed class MapperSpcClusterer : IGraphClusterer
{
    private readonly PottsModelConfig _sampler;
    private readonly int _gridSteps;
    private readonly double _cutTheta;

    /// <summary>
    /// Construct a MAPPER SPC clusterer.
    /// </summary>
    /// <param name="sampler">
    /// Potts sampler configuration (number of colors, etc.). <c>null</c> uses defaults.
    /// </param>
    /// <param name="cutTheta">
    /// Friends-of-friends threshold for the spin-agreement partition strategy.
    /// Edges with average spin agreement > <paramref name="cutTheta"/> become
    /// "friends" and union their endpoints. Default 0.5 per Blatt 1996; the
    /// paper reports the cut is robust across <c>0.2 &lt;= θ &lt;= 0.9</c>.
    /// </param>
    /// <param name="gridSteps">
    /// Number of log-spaced temperatures in the per-patch auto-bracketed grid.
    /// </param>
    public MapperSpcClusterer(
        PottsModelConfig? sampler = null,
        double cutTheta = 0.5,
        int gridSteps = 12)
    {
        if (cutTheta < 0.0 || cutTheta > 1.0)
            throw new ArgumentOutOfRangeException(nameof(cutTheta), "θ must be in [0, 1].");
        if (gridSteps < 2)
            throw new ArgumentOutOfRangeException(nameof(gridSteps), "gridSteps must be at least 2.");
        _sampler = sampler ?? new PottsModelConfig();
        _gridSteps = gridSteps;
        _cutTheta = cutTheta;
    }

    public string Name => $"SPC fixed-grid (Blatt FoF θ={_cutTheta:F2})";

    /// <summary>
    /// Last <see cref="SweepSummary"/> from the most recent
    /// <see cref="ClusterInduced"/> call. Useful for introspection and for the
    /// "did we find T_c with confidence" question. Thread-affine — if you
    /// run MAPPER patches in parallel this will be racy; pull status
    /// information from a wrapper instead in that case.
    /// </summary>
    public SweepSummary? LastScheduleSummary { get; private set; }

    public ClusterResult ClusterInduced(CsrGraph graph, IReadOnlyList<int> preimageIndices)
    {
        if (preimageIndices is null) throw new ArgumentNullException(nameof(preimageIndices));
        int k = preimageIndices.Count;
        if (k == 0) return new ClusterResult(Array.Empty<int>(), 0);
        if (k == 1) return new ClusterResult(new[] { 0 }, 1);

        // Build the induced-subgraph mask from the preimage.
        var mask = new bool[graph.NodeCount];
        foreach (int idx in preimageIndices)
        {
            if (idx < 0 || idx >= graph.NodeCount)
                throw new ArgumentOutOfRangeException(nameof(preimageIndices),
                    $"preimage index {idx} out of range [0, {graph.NodeCount}).");
            mask[idx] = true;
        }

        var subgraph = graph.InducedSubgraph(mask, out int[] newToOld, out int[] oldToNew);

        if (subgraph.NodeCount < 2)
            return SingletonsPerPoint(k);

        // Run an SPC fixed-grid sweep on an auto-bracketed grid for this patch.
        var cfg = AutoGridFixedSweep.BuildConfig(subgraph, _sampler, gridSteps: _gridSteps);
        var result = new FixedGridSweepStrategy(cfg).Run(subgraph);
        LastScheduleSummary = result.Summary;

        if (result.ChosenAffinities.G.Length == 0 || result.ChosenAlignments is null)
            return SingletonsPerPoint(k);

        // Apply the SPC partition strategy backed by threshold spin agreement.
        var partition = new ThresholdSpinAgreement { Theta = _cutTheta }
            .Apply(result.Graph, result.ChosenAffinities, result.ChosenAlignments);

        // Translate subgraph-local labels back to preimage-index order.
        // IGraphClusterer contract: labels[i] is the cluster for preimageIndices[i].
        var preimageLabels = new int[k];
        for (int i = 0; i < k; i++)
        {
            int origIdx = preimageIndices[i];
            int subIdx  = oldToNew[origIdx];
            preimageLabels[i] = partition.Labels[subIdx];
        }
        return new ClusterResult(preimageLabels, partition.Count);
    }

    /// <summary>
    /// Trivial-cluster fallback: assign every point its own cluster.
    /// Returned when the subgraph is degenerate or the sweep produced
    /// no observables; the conservative read is "no claimed structure."
    /// </summary>
    private static ClusterResult SingletonsPerPoint(int k)
    {
        var labels = new int[k];
        for (int i = 0; i < k; i++) labels[i] = i;
        return new ClusterResult(labels, k);
    }
}
