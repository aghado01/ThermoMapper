using System;
using Clustering.Dendrograms;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// R1 completion stage: the periphery policies. Hand-computed oracles pin
/// the ascent semantics (uphill adoption, chaining through unassigned
/// terrain, local-max abstain, valley respect) and the capture semantics
/// (edge-greedy chains, mutual-best ORBITS abstain) — including the runtime
/// realization of the capture ≁ ascent separation the lemma stack predicted.
/// </summary>
public sealed class PeripheryPolicyTests
{
    private static CsrGraph PathGraph(int n)
    {
        var edges = new Edge[n - 1];
        for (int i = 0; i < n - 1; i++) edges[i] = new Edge(i, i + 1, 1.0);
        return CsrGraph.FromEdges(edges, n);
    }

    private static Assignment Partial(params int[] labels) => new()
    {
        Labels = labels,
        Count = 1 + Math.Max(-1, System.Linq.Enumerable.Max(labels)),
    };

    [Fact]
    public void Ascend_AdoptsUphill_ChainingThroughUnassignedTerrain()
    {
        // 0(A) — 1 — 2 with L = [3,2,1]: node 1 steps uphill to 0 (assigned);
        // node 2 steps to 1 (unassigned) then to 0 — chains resolve against
        // ORIGINAL labels, so both adopt A.
        var graph = PathGraph(3);
        var resolved = PeripheryPolicies.Ascend(
            Partial(0, Assignment.Unassigned, Assignment.Unassigned), graph, new[] { 3.0, 2.0, 1.0 });

        Assert.Equal(new[] { 0, 0, 0 }, resolved.Labels);
    }

    [Fact]
    public void Ascend_LocalMaxAbstains()
    {
        // L = [1,5,1]: node 1 is a mode of its own — honest abstain; node 2
        // ascends into that unassigned mode and abstains too.
        var graph = PathGraph(3);
        var resolved = PeripheryPolicies.Ascend(
            Partial(0, Assignment.Unassigned, Assignment.Unassigned), graph, new[] { 1.0, 5.0, 1.0 });

        Assert.Equal(new[] { 0, Assignment.Unassigned, Assignment.Unassigned }, resolved.Labels);
    }

    [Fact]
    public void Ascend_RespectsTheValley()
    {
        // 0(A) — 1 — 2 — 3 — 4(B), L = [5,2,1,3,4]: the valley floor at node 2
        // ascends toward B's side; node 1 toward A — neither crosses the floor.
        var graph = PathGraph(5);
        var resolved = PeripheryPolicies.Ascend(
            Partial(0, Assignment.Unassigned, Assignment.Unassigned, Assignment.Unassigned, 1),
            graph,
            new[] { 5.0, 2.0, 1.0, 3.0, 4.0 });

        Assert.Equal(new[] { 0, 0, 1, 1, 1 }, resolved.Labels);
    }

    [Fact]
    public void Capture_OrbitsAbstain_WhereAscentResolves()
    {
        // 0(A) — 1 — 2 — 3: the (1,2) edge is the strongest field on both
        // sides, so capture chains 1→2→1 — a mutual-best ORBIT — and 3 chains
        // into it: all three abstain. Ascent on a monotone landscape resolves
        // all three. The lemma-stack separation (capture ≁ ascent), live.
        var graph = PathGraph(4);
        var partial = Partial(0, Assignment.Unassigned, Assignment.Unassigned, Assignment.Unassigned);

        var edgeField = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            edgeField[edge.Slot] = (edge.Source, edge.Target) switch
            {
                (0, 1) or (1, 0) => 0.1,
                (1, 2) or (2, 1) => 0.9,
                _                => 0.8,
            };

        var captured = PeripheryPolicies.Capture(partial, graph, edgeField);
        Assert.Equal(
            new[] { 0, Assignment.Unassigned, Assignment.Unassigned, Assignment.Unassigned },
            captured.Labels);

        var ascended = PeripheryPolicies.Ascend(partial, graph, new[] { 9.0, 3.0, 2.0, 1.0 });
        Assert.Equal(new[] { 0, 0, 0, 0 }, ascended.Labels);
    }

    [Fact]
    public void Capture_ChainsToTheAssignedSide()
    {
        // Field rises toward the assigned core: 1's best edge points at 0.
        var graph = PathGraph(3);
        var edgeField = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            edgeField[edge.Slot] = (edge.Source, edge.Target) is (0, 1) or (1, 0) ? 0.9 : 0.4;

        var resolved = PeripheryPolicies.Capture(
            Partial(0, Assignment.Unassigned, Assignment.Unassigned), graph, edgeField);

        Assert.Equal(new[] { 0, 0, 0 }, resolved.Labels);
    }

    [Fact]
    public void SelectByExcessOfMass_MinClusterSizeFloorsEligibility()
    {
        // FourLeafTree masses [6,4,4]: with minClusterSize=3 the size-2 pairs
        // are ineligible (their mass passes through, not blocks); the root is
        // the only size-3+ node — ineligible by default, selectable with
        // allowRoot.
        var tree = new Dendrogram(
            new[]
            {
                new DendrogramNode(0, 1, 1.0, 2),
                new DendrogramNode(2, 3, 2.0, 2),
                new DendrogramNode(4, 5, 4.0, 4),
            },
            LeafCount: 4,
            CostAxis: "temperature");
        var mass = new[] { 6.0, 4.0, 4.0 };

        Assert.Equal(
            new[] { false, false, false },
            LandscapeWalk.SelectByExcessOfMass(tree, mass, minClusterSize: 3));
        Assert.Equal(
            new[] { false, false, true },
            LandscapeWalk.SelectByExcessOfMass(tree, mass, allowRoot: true, minClusterSize: 3));
    }
}
