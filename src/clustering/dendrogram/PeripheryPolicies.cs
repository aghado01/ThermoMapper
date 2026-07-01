using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// The completion stage of the resolution layer: policies that resolve a
/// partial <see cref="Assignment"/>'s Unassigned points AFTER a selector
/// (cut / EOM / peak-select) has decided which clusters exist. Inputs stay
/// heterogeneous by design — ascent consumes a per-node landscape slice
/// (ordinal: only the field's order matters, so it is gauge-free), capture
/// consumes a per-edge field — and only the output unifies. Selector
/// verdicts are final: policies fill Unassigned slots, never relabel, and
/// every walk resolves against the ORIGINAL labels, so the result is
/// independent of node order.
/// </summary>
public static class PeripheryPolicies
{
    /// <summary>
    /// Height-greedy modal ascent (the quick-shift family): each unassigned
    /// node follows strictly-uphill steps in the landscape slice — ties
    /// broken by lower index, the same total order as the ascent lemmas —
    /// until it reaches an assigned node (adopts its label) or a local
    /// maximum (stays abstained: being a mode of one's own is an honest
    /// answer). Cannot cross a valley by construction.
    /// </summary>
    public static Assignment Ascend(Assignment partial, CsrGraph graph, double[] landscapeSlice)
    {
        ArgumentNullException.ThrowIfNull(partial);
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(landscapeSlice);
        int n = graph.NodeCount;
        if (partial.Labels.Length != n)
            throw new InvalidOperationException(
                $"Assignment length ({partial.Labels.Length}) != graph node count ({n}).");
        if (landscapeSlice.Length != n)
            throw new InvalidOperationException(
                $"Landscape slice length ({landscapeSlice.Length}) != graph node count ({n}).");

        int[] original = partial.Labels;
        var resolved = (int[])original.Clone();
        int[] rowPtr = graph.RowPointers;
        int[] targets = graph.Targets;

        for (int i = 0; i < n; i++)
        {
            if (original[i] != Assignment.Unassigned) continue;

            int cur = i;
            // Strict ascent terminates in ≤ n steps (the key strictly increases).
            for (int step = 0; step <= n; step++)
            {
                int best = cur;
                for (int s = rowPtr[cur]; s < rowPtr[cur + 1]; s++)
                {
                    int j = targets[s];
                    // key(j) > key(best): higher landscape, or equal and lower index.
                    if (landscapeSlice[j] > landscapeSlice[best]
                        || (landscapeSlice[j] == landscapeSlice[best] && j < best))
                    {
                        best = j;
                    }
                }
                if (best == cur) break;                      // local max — honest abstain
                if (original[best] != Assignment.Unassigned) // adopt the selector's verdict
                {
                    resolved[i] = original[best];
                    break;
                }
                cur = best;
            }
        }

        return new Assignment { Labels = resolved, Count = partial.Count };
    }

    /// <summary>
    /// Edge-greedy capture (Domany's step 2 recast as a completion): each
    /// unassigned node follows max-field edges until it reaches an assigned
    /// node or revisits one (an orbit — e.g. duplicate pairs that are each
    /// other's best neighbor — stays abstained). No potential is followed:
    /// capture can chain across a landscape valley, which is exactly its
    /// A/B against <see cref="Ascend"/>.
    /// </summary>
    public static Assignment Capture(Assignment partial, CsrGraph graph, double[] edgeField)
    {
        ArgumentNullException.ThrowIfNull(partial);
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(edgeField);
        int n = graph.NodeCount;
        if (partial.Labels.Length != n)
            throw new InvalidOperationException(
                $"Assignment length ({partial.Labels.Length}) != graph node count ({n}).");
        if (edgeField.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"Edge field length ({edgeField.Length}) != CSR slot count ({graph.Targets.Length}).");

        // The per-edge field lives on the canonical j>i slots; mirror it so
        // every directed slot carries its undirected edge's value.
        var slotField = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            double value = edgeField[edge.Slot];
            slotField[edge.Slot] = value;
            for (int s = graph.RowPointers[edge.Target]; s < graph.RowPointers[edge.Target + 1]; s++)
                if (graph.Targets[s] == edge.Source) { slotField[s] = value; break; }
        }

        int[] original = partial.Labels;
        var resolved = (int[])original.Clone();
        int[] rowPtr = graph.RowPointers;
        int[] targets = graph.Targets;
        var visited = new HashSet<int>();

        for (int i = 0; i < n; i++)
        {
            if (original[i] != Assignment.Unassigned) continue;

            visited.Clear();
            int cur = i;
            while (visited.Add(cur))
            {
                int bestTarget = -1;
                double bestValue = double.NegativeInfinity;
                for (int s = rowPtr[cur]; s < rowPtr[cur + 1]; s++)
                {
                    int j = targets[s];
                    double value = slotField[s];
                    if (value > bestValue || (value == bestValue && j < bestTarget))
                    {
                        bestValue = value;
                        bestTarget = j;
                    }
                }
                if (bestTarget < 0) break;                         // isolated node
                if (original[bestTarget] != Assignment.Unassigned) // adopt
                {
                    resolved[i] = original[bestTarget];
                    break;
                }
                cur = bestTarget;                                  // chain on; orbits abstain
            }
        }

        return new Assignment { Labels = resolved, Count = partial.Count };
    }
}
