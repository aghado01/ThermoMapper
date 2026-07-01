#nullable enable
using System;
using Graphs.Primitives;
using TDA.Primitives;
using Xunit;

namespace TDA.Mapper.Tests;

/// <summary>
/// Direct structural tests for <see cref="NerveDiff"/> — node matching,
/// edge diff, and component event classification.
/// </summary>
public sealed class NerveDiffTests
{
    // Frame 0: two isolated nerve nodes (0→{0,1}, 1→{2,3}).
    // Frame 1: same nodes connected by one edge, same member sets.
    // Expected diff: 2 node matches, 1 born edge, 0 died edges, 1 Merge event.

    [Fact]
    public void PathGraph_NerveDiff_NodeMatchesEdgesAndMergeEvent()
    {
        var frame0 = new NerveFiltrationFrame(
            ParameterValue: 0.0,
            Nerve: CsrGraph.FromEdges(Array.Empty<Edge>(), 2),
            NodeMemberIndices: new[] { new[] { 0, 1 }, new[] { 2, 3 } },
            FrameIndex: 0);

        var frame1 = new NerveFiltrationFrame(
            ParameterValue: 1.0,
            Nerve: CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0) }, 2),
            NodeMemberIndices: new[] { new[] { 0, 1 }, new[] { 2, 3 } },
            FrameIndex: 1);

        var filtration = new NerveFiltration(new[] { frame0, frame1 }, "T");
        var diffs = filtration.ComputeDiffs();

        Assert.Single(diffs);
        var diff = diffs[0];

        // Frame parameters.
        Assert.Equal(0.0, diff.ParameterFrom);
        Assert.Equal(1.0, diff.ParameterTo);

        // Node matching: each from-node has an identical to-node.
        Assert.Equal(2, diff.NodeMatches.Count);

        // Edge diff: edge (0,1) is new; no edges died.
        Assert.Equal((0, 1), Assert.Single(diff.BornEdges));
        Assert.Empty(diff.DiedEdges);

        // Component events: the two isolated CCs merge into one.
        var evt = Assert.Single(diff.ComponentEvents);
        Assert.Equal(ComponentEventKind.Merge, evt.Kind);
        Assert.Equal(2, evt.CcsFrom.Count);
        Assert.Single(evt.CcsTo);
    }

    // Frame 0: one connected pair (0→{0,1}, 1→{2,3}, edge 0-1).
    // Frame 1: the pair splits into two isolated nodes with the same member sets.
    // Expected diff: 2 node matches, 0 born edges, 1 died edge, 1 Split event.

    [Fact]
    public void DisconnectingGraph_NerveDiff_SplitEvent()
    {
        var frame0 = new NerveFiltrationFrame(
            ParameterValue: 0.0,
            Nerve: CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0) }, 2),
            NodeMemberIndices: new[] { new[] { 0, 1 }, new[] { 2, 3 } },
            FrameIndex: 0);

        var frame1 = new NerveFiltrationFrame(
            ParameterValue: 1.0,
            Nerve: CsrGraph.FromEdges(Array.Empty<Edge>(), 2),
            NodeMemberIndices: new[] { new[] { 0, 1 }, new[] { 2, 3 } },
            FrameIndex: 1);

        var filtration = new NerveFiltration(new[] { frame0, frame1 }, "T");
        var diff = filtration.ComputeDiffs()[0];

        Assert.Equal(2, diff.NodeMatches.Count);
        Assert.Empty(diff.BornEdges);
        Assert.Equal((0, 1), Assert.Single(diff.DiedEdges));

        var evt = Assert.Single(diff.ComponentEvents);
        Assert.Equal(ComponentEventKind.Split, evt.Kind);
        Assert.Single(evt.CcsFrom);
        Assert.Equal(2, evt.CcsTo.Count);
    }
}
