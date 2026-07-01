#nullable enable
using System;
using System.Collections.Generic;

namespace Graphs.Primitives.Traversal;

/// <summary>
/// Full-graph bounded single-source shortest path over a <see cref="CsrGraph"/>.
/// Designed for parallel execution by injecting thread-local scratch space.
/// </summary>
/// <remarks>
/// <para>For edge-local refinement of an existing neighbor graph, see
/// the <c>Graphs.Pipeline.Refinement.PathNeighborRefiner</c>
/// implementation in <c>Graphs.Pipeline.Refinement</c>. It uses a
/// target-masked search with early exit, not this helper.</para>
/// </remarks>
public static class Dijkstra
{
    /// <summary>
    /// Computes bounded shortest-path distances from a source across the full CSR graph.
    /// </summary>
    /// <param name="graph">The global topology and spatial weights.</param>
    /// <param name="sourceNode">The index of the starting node.</param>
    /// <param name="distances">Thread-local scratch span (length >= NodeCount). Will be filled with geodesic distances.</param>
    /// <param name="hops">Thread-local scratch span (length >= NodeCount). Tracks topological depth.</param>
    /// <param name="scratchQueue">Thread-local priority queue. Must be empty when passed.</param>
    /// <param name="maxDistance">Abort traversing paths longer than this metric distance.</param>
    /// <param name="maxHops">Abort traversing paths deeper than this topological hop-count.</param>
    public static void ComputeBoundedDistances(
        CsrGraph graph,
        int sourceNode,
        Span<double> distances,
        Span<int> hops,
        PriorityQueue<int, double> scratchQueue,
        double maxDistance = double.PositiveInfinity,
        int maxHops = int.MaxValue)
    {
        if (sourceNode < 0 || sourceNode >= graph.NodeCount)
            throw new ArgumentOutOfRangeException(nameof(sourceNode));

        // 1. Initialize scratch spaces
        distances.Fill(double.PositiveInfinity);
        hops.Fill(int.MaxValue);

        distances[sourceNode] = 0.0;
        hops[sourceNode] = 0;

        scratchQueue.Clear();
        scratchQueue.Enqueue(sourceNode, 0.0);

        ReadOnlySpan<int> rowPtrs = graph.RowPointers.AsSpan();
        ReadOnlySpan<int> targets = graph.Targets.AsSpan();
        ReadOnlySpan<double> weights = graph.Weights.AsSpan();

        // 2. Greedy traversal
        while (scratchQueue.TryDequeue(out int u, out double d))
        {
            // Standard stale-entry check for binary heaps without Decrease-Key support
            if (d > distances[u])
                continue;

            int currentHops = hops[u];

            // Enforce bounds: If we hit the hop limit or distance limit,
            // we do not explore outgoing edges from this node.
            if (currentHops >= maxHops || d >= maxDistance)
                continue;

            int start = rowPtrs[u];
            int end = rowPtrs[u + 1];

            for (int e = start; e < end; e++)
            {
                int v = targets[e];
                double weight = weights[e];

                // Safety guard against degenerate negative weights
                if (weight < 0.0) weight = 0.0;

                double altDistance = d + weight;
                int altHops = currentHops + 1;

                // Relaxation step
                if (altDistance < distances[v])
                {
                    distances[v] = altDistance;
                    hops[v] = altHops;

                    // Note: .NET's PriorityQueue allows duplicate elements.
                    // The stale-entry check at the top of the loop handles the ghost entries.
                    scratchQueue.Enqueue(v, altDistance);
                }
            }
        }

        // Ensure queue is clean for the next parallel thread that uses it
        scratchQueue.Clear();
    }

    /// <summary>
    /// Computes bounded shortest-path distances from a source across the full CSR graph,
    /// with optional early termination once all masked targets are settled.
    /// </summary>
    /// <param name="graph">The global topology and spatial weights.</param>
    /// <param name="sourceNode">The index of the starting node.</param>
    /// <param name="distances">Thread-local scratch span (length >= NodeCount). Will be filled with geodesic distances.</param>
    /// <param name="hops">Thread-local scratch span (length >= NodeCount). Tracks topological depth.</param>
    /// <param name="scratchQueue">Thread-local priority queue. Must be empty when passed.</param>
    /// <param name="targetMask">Early-exit mask. Traversal is not restricted by the mask; it only stops once all
    /// marked targets have been settled.</param>
    /// <param name="maxDistance">Abort traversing paths longer than this metric distance.</param>
    /// <param name="maxHops">Abort traversing paths deeper than this topological hop-count.</param>
    public static void ComputeBoundedDistances(
        CsrGraph graph,
        int sourceNode,
        Span<double> distances,
        Span<int> hops,
        PriorityQueue<int, double> scratchQueue,
        ReadOnlySpan<bool> targetMask,
        double maxDistance = double.PositiveInfinity,
        int maxHops = int.MaxValue)
    {
        if (sourceNode < 0 || sourceNode >= graph.NodeCount)
            throw new ArgumentOutOfRangeException(nameof(sourceNode));

        int remainingTargets = 0;
        for (int i = 0; i < targetMask.Length; i++)
        {
            if (targetMask[i])
                remainingTargets++;
        }

        distances.Fill(double.PositiveInfinity);
        hops.Fill(int.MaxValue);

        distances[sourceNode] = 0.0;
        hops[sourceNode] = 0;

        scratchQueue.Clear();
        scratchQueue.Enqueue(sourceNode, 0.0);

        ReadOnlySpan<int> rowPtrs = graph.RowPointers.AsSpan();
        ReadOnlySpan<int> targets = graph.Targets.AsSpan();
        ReadOnlySpan<double> weights = graph.Weights.AsSpan();

        while (scratchQueue.TryDequeue(out int u, out double d))
        {
            if (d > distances[u])
                continue;

            int currentHops = hops[u];

            if (currentHops >= maxHops || d >= maxDistance)
                continue;

            if (remainingTargets > 0 && u < targetMask.Length && targetMask[u])
            {
                remainingTargets--;
                if (remainingTargets == 0)
                    break;
            }

            int start = rowPtrs[u];
            int end = rowPtrs[u + 1];

            for (int e = start; e < end; e++)
            {
                int v = targets[e];
                double weight = weights[e];

                if (weight < 0.0) weight = 0.0;

                double altDistance = d + weight;
                int altHops = currentHops + 1;

                if (altDistance < distances[v])
                {
                    distances[v] = altDistance;
                    hops[v] = altHops;
                    scratchQueue.Enqueue(v, altDistance);
                }
            }
        }

        scratchQueue.Clear();
    }
}
