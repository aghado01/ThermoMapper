#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;

namespace TDA.Primitives;

// ── NodeMatch ─────────────────────────────────────────────────────────────────

/// <summary>
/// How one nerve node in the earlier frame maps to a nerve node in the later frame,
/// measured by shared member-point count.
/// </summary>
public sealed record NodeMatch(
    int NodeFrom,      // index in the earlier frame's nerve
    int NodeTo,        // index in the later frame's nerve
    int OverlapCount); // number of shared member-point indices

// ── ComponentEventKind ────────────────────────────────────────────────────────

public enum ComponentEventKind
{
    Birth,        // CC appears with no predecessor in the earlier frame
    Death,        // CC disappears with no successor in the later frame
    Continuation, // 1:1 mapping (CC continues, possibly reshaped)
    Merge,        // multiple earlier CCs → one later CC (elder rule applies for H0)
    Split,        // one earlier CC → multiple later CCs
}

// ── ComponentEvent ────────────────────────────────────────────────────────────

/// <summary>
/// A connected-component topology event between two consecutive nerve frames.
/// CC indices are frame-local (valid within their respective frame).
/// </summary>
public sealed record ComponentEvent(
    ComponentEventKind Kind,
    IReadOnlyList<int> CcsFrom,  // earlier-frame CC indices (empty for Birth)
    IReadOnlyList<int> CcsTo);   // later-frame CC indices (empty for Death)

// ── NerveDiff ─────────────────────────────────────────────────────────────────

/// <summary>
/// Structural diff between two consecutive frames in a <see cref="NerveFiltration"/>.
/// <para>
/// Node and edge indices are frame-local: indices in <see cref="BornNodes"/>,
/// <see cref="DiedNodes"/>, and <see cref="NodeMatches"/> refer to the nerve
/// node numbering within their respective frames.  Edge pairs <c>(U, V)</c>
/// with <c>U &lt; V</c> use the same convention.
/// </para>
/// <para>
/// <see cref="ComponentEvents"/> is the derived topology view: birth, death,
/// continuation, merge, and split events at the connected-component level,
/// computed from <see cref="NodeMatches"/> and the CC structure of each frame.
/// The H0 persistence barcode is a function of this event sequence.
/// </para>
/// </summary>
public sealed class NerveDiff
{
    public double ParameterFrom { get; }
    public double ParameterTo { get; }
    public int FrameIndexFrom { get; }
    public int FrameIndexTo { get; }

    // Node-level structural changes.
    public IReadOnlyList<NodeMatch> NodeMatches { get; }
    public IReadOnlyList<int> BornNodes { get; }    // later-frame nodes with no predecessor
    public IReadOnlyList<int> DiedNodes { get; }    // earlier-frame nodes with no successor

    // Edge-level structural changes (U < V, frame-local node indices).
    public IReadOnlyList<(int U, int V)> BornEdges { get; }
    public IReadOnlyList<(int U, int V)> DiedEdges { get; }

    // Connected-component topology events derived from node matching + CC structure.
    public IReadOnlyList<ComponentEvent> ComponentEvents { get; }

    private NerveDiff(
        double paramFrom, double paramTo, int frameFrom, int frameTo,
        IReadOnlyList<NodeMatch> nodeMatches,
        IReadOnlyList<int> bornNodes, IReadOnlyList<int> diedNodes,
        IReadOnlyList<(int, int)> bornEdges, IReadOnlyList<(int, int)> diedEdges,
        IReadOnlyList<ComponentEvent> componentEvents)
    {
        ParameterFrom = paramFrom; ParameterTo = paramTo;
        FrameIndexFrom = frameFrom; FrameIndexTo = frameTo;
        NodeMatches = nodeMatches;
        BornNodes = bornNodes; DiedNodes = diedNodes;
        BornEdges = bornEdges; DiedEdges = diedEdges;
        ComponentEvents = componentEvents;
    }

    // ── Factory ───────────────────────────────────────────────────────────────

    internal static NerveDiff Compute(NerveFiltrationFrame from, NerveFiltrationFrame to)
    {
        var nodeMatches = MatchNodes(from, to);

        int nFrom = from.Nerve.NodeCount;
        int nTo = to.Nerve.NodeCount;

        var matchedFromSet = new HashSet<int>(nodeMatches.Select(m => m.NodeFrom));
        var matchedToSet = new HashSet<int>(nodeMatches.Select(m => m.NodeTo));

        var bornNodes = new List<int>();
        for (int j = 0; j < nTo; j++)
            if (!matchedToSet.Contains(j)) bornNodes.Add(j);

        var diedNodes = new List<int>();
        for (int i = 0; i < nFrom; i++)
            if (!matchedFromSet.Contains(i)) diedNodes.Add(i);

        var (bornEdges, diedEdges) = DiffEdges(from, to, nodeMatches);
        var componentEvents = ComputeComponentEvents(from, to, nodeMatches);

        return new NerveDiff(
            from.ParameterValue, to.ParameterValue,
            from.FrameIndex, to.FrameIndex,
            nodeMatches, bornNodes, diedNodes,
            bornEdges, diedEdges, componentEvents);
    }

    // ── Node matching ─────────────────────────────────────────────────────────

    private static IReadOnlyList<NodeMatch> MatchNodes(
        NerveFiltrationFrame from, NerveFiltrationFrame to)
    {
        int nFrom = from.Nerve.NodeCount;
        int nTo = to.Nerve.NodeCount;
        if (nFrom == 0 || nTo == 0) return Array.Empty<NodeMatch>();

        // Pre-build to-node member sets for O(1) lookup during overlap counting.
        var toSets = new HashSet<int>[nTo];
        for (int j = 0; j < nTo; j++)
            toSets[j] = new HashSet<int>(to.NodeMemberIndices[j]);

        var matches = new List<NodeMatch>(nFrom);
        for (int i = 0; i < nFrom; i++)
        {
            int bestJ = -1, bestOverlap = 0;
            foreach (int pt in from.NodeMemberIndices[i])
            {
                for (int j = 0; j < nTo; j++)
                {
                    if (toSets[j].Contains(pt))
                    {
                        // Count full overlap lazily when we find a candidate.
                        int overlap = CountIntersection(from.NodeMemberIndices[i], toSets[j]);
                        if (overlap > bestOverlap) { bestOverlap = overlap; bestJ = j; }
                        break; // a point lives in at most one Mapper cluster
                    }
                }
            }
            if (bestJ >= 0)
                matches.Add(new NodeMatch(i, bestJ, bestOverlap));
        }
        return matches;
    }

    // ── Edge diff ─────────────────────────────────────────────────────────────

    private static (List<(int, int)> Born, List<(int, int)> Died) DiffEdges(
        NerveFiltrationFrame from, NerveFiltrationFrame to,
        IReadOnlyList<NodeMatch> nodeMatches)
    {
        // Build forward and reverse node maps.
        var matchFrom = new Dictionary<int, int>(); // from-node → to-node
        var matchTo = new Dictionary<int, int>();   // to-node → from-node (first match wins for merges)
        foreach (var m in nodeMatches)
        {
            matchFrom[m.NodeFrom] = m.NodeTo;
            matchTo.TryAdd(m.NodeTo, m.NodeFrom);
        }

        var fromEdges = CollectEdgeSet(from.Nerve);
        var toEdges = CollectEdgeSet(to.Nerve);

        var born = new List<(int, int)>();
        foreach (var (u, v) in toEdges)
        {
            bool continued = matchTo.TryGetValue(u, out int uFrom)
                          && matchTo.TryGetValue(v, out int vFrom)
                          && fromEdges.Contains((Math.Min(uFrom, vFrom), Math.Max(uFrom, vFrom)));
            if (!continued) born.Add((u, v));
        }

        var died = new List<(int, int)>();
        foreach (var (u, v) in fromEdges)
        {
            bool continued = matchFrom.TryGetValue(u, out int uTo)
                          && matchFrom.TryGetValue(v, out int vTo)
                          && toEdges.Contains((Math.Min(uTo, vTo), Math.Max(uTo, vTo)));
            if (!continued) died.Add((u, v));
        }

        return (born, died);
    }

    private static HashSet<(int, int)> CollectEdgeSet(CsrGraph nerve)
    {
        var set = new HashSet<(int, int)>();
        for (int u = 0; u < nerve.NodeCount; u++)
        {
            int start = nerve.RowPointers[u], end = nerve.RowPointers[u + 1];
            for (int e = start; e < end; e++)
            {
                int v = nerve.Targets[e];
                if (u < v) set.Add((u, v));
            }
        }
        return set;
    }

    // ── Component events ──────────────────────────────────────────────────────

    private static IReadOnlyList<ComponentEvent> ComputeComponentEvents(
        NerveFiltrationFrame from, NerveFiltrationFrame to,
        IReadOnlyList<NodeMatch> nodeMatches)
    {
        int nFromCcs = CountCcs(from.Nerve, out int[] ccsFrom);
        int nToCcs = CountCcs(to.Nerve, out int[] ccsTo);

        // Build CC-level connection maps.
        // fromCcToToCcCount[ccFrom][ccTo] = number of matched nodes bridging that CC pair.
        var fromCcToToCcCount = new Dictionary<int, Dictionary<int, int>>();
        var toCcToFromCcs = new Dictionary<int, HashSet<int>>();
        for (int cc = 0; cc < nFromCcs; cc++) fromCcToToCcCount[cc] = new Dictionary<int, int>();
        for (int cc = 0; cc < nToCcs; cc++) toCcToFromCcs[cc] = new HashSet<int>();

        foreach (var m in nodeMatches)
        {
            int cfrom = ccsFrom[m.NodeFrom];
            int cto = ccsTo[m.NodeTo];
            fromCcToToCcCount[cfrom].TryGetValue(cto, out int cur);
            fromCcToToCcCount[cfrom][cto] = cur + 1;
            toCcToFromCcs[cto].Add(cfrom);
        }

        var events = new List<ComponentEvent>();
        var handledFromCcs = new HashSet<int>();
        var handledToCcs = new HashSet<int>();

        // Merges: multiple from-CCs → one to-CC.
        for (int ccTo = 0; ccTo < nToCcs; ccTo++)
        {
            var sources = toCcToFromCcs[ccTo];
            if (sources.Count <= 1) continue;

            var ccsFrom2 = sources.OrderBy(x => x).ToArray();
            events.Add(new ComponentEvent(ComponentEventKind.Merge, ccsFrom2, new[] { ccTo }));
            foreach (int s in ccsFrom2) handledFromCcs.Add(s);
            handledToCcs.Add(ccTo);
        }

        // Splits: one from-CC → multiple to-CCs.
        for (int ccFrom = 0; ccFrom < nFromCcs; ccFrom++)
        {
            if (handledFromCcs.Contains(ccFrom)) continue;
            var targets = fromCcToToCcCount[ccFrom];
            if (targets.Count <= 1) continue;

            // Order targets by overlap count descending so CcsTo[0] has the highest overlap.
            var orderedTargets = targets.OrderByDescending(kv => kv.Value).Select(kv => kv.Key).ToArray();
            events.Add(new ComponentEvent(ComponentEventKind.Split, new[] { ccFrom }, orderedTargets));
            handledFromCcs.Add(ccFrom);
            foreach (int t in orderedTargets) handledToCcs.Add(t);
        }

        // Deaths: from-CCs with no matching to-nodes.
        for (int ccFrom = 0; ccFrom < nFromCcs; ccFrom++)
        {
            if (handledFromCcs.Contains(ccFrom)) continue;
            if (fromCcToToCcCount[ccFrom].Count == 0)
            {
                events.Add(new ComponentEvent(ComponentEventKind.Death, new[] { ccFrom }, Array.Empty<int>()));
                handledFromCcs.Add(ccFrom);
            }
        }

        // Continuations: 1:1 remaining.
        for (int ccFrom = 0; ccFrom < nFromCcs; ccFrom++)
        {
            if (handledFromCcs.Contains(ccFrom)) continue;
            var targets = fromCcToToCcCount[ccFrom];
            if (targets.Count != 1) continue;
            int ccTo = targets.Keys.First();
            if (handledToCcs.Contains(ccTo)) continue;

            events.Add(new ComponentEvent(ComponentEventKind.Continuation, new[] { ccFrom }, new[] { ccTo }));
            handledFromCcs.Add(ccFrom);
            handledToCcs.Add(ccTo);
        }

        // Births: to-CCs with no predecessors.
        for (int ccTo = 0; ccTo < nToCcs; ccTo++)
        {
            if (!handledToCcs.Contains(ccTo))
                events.Add(new ComponentEvent(ComponentEventKind.Birth, Array.Empty<int>(), new[] { ccTo }));
        }

        return events;
    }

    // ── Shared utilities ──────────────────────────────────────────────────────

    internal static int CountCcs(CsrGraph nerve, out int[] labels)
    {
        int n = nerve.NodeCount;
        if (n == 0) { labels = Array.Empty<int>(); return 0; }

        var uf = new UnionFind(n);
        for (int u = 0; u < n; u++)
        {
            int start = nerve.RowPointers[u], end = nerve.RowPointers[u + 1];
            for (int e = start; e < end; e++)
            {
                int v = nerve.Targets[e];
                if (u < v) uf.Union(u, v);
            }
        }

        var rootToIdx = new Dictionary<int, int>();
        labels = new int[n];
        int nextCc = 0;
        for (int i = 0; i < n; i++)
        {
            int root = uf.Find(i);
            if (!rootToIdx.TryGetValue(root, out int idx))
                rootToIdx[root] = idx = nextCc++;
            labels[i] = idx;
        }
        return nextCc;
    }

    private static int CountIntersection(int[] a, HashSet<int> bSet)
    {
        int count = 0;
        foreach (int x in a)
            if (bSet.Contains(x)) count++;
        return count;
    }
}
