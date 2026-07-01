#nullable enable
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Threading.Tasks;
using Graphs.Primitives;
using Graphs.Primitives.Traversal;

namespace Graphs.Pipeline.Refinement;

/// <summary>
/// Stage 4 — path-neighbor refinement.
/// Recomputes edge distances over the repaired directed topology using
/// bounded single-source shortest-path search. The adjacency topology is
/// preserved; only the per-edge distance values are updated.
/// </summary>
/// <remarks>
/// <para>Uses an inline target-masked SSSP (early exit once all declared
/// edges are refined). The shared <see cref="Graphs.Primitives.Traversal.Dijkstra"/>
/// helper computes full-graph distances for callers that need every reachable
/// node.</para>
/// </remarks>
public sealed class PathNeighborRefiner : IMetricRefiner
{
    private readonly double? _maxDistance;

    public PathNeighborRefiner(double? maxDistance = null)
    {
        if (maxDistance.HasValue && maxDistance.Value <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(maxDistance), "MaxDistance must be positive when specified.");

        _maxDistance = maxDistance;
    }

    public NeighborSelection Refine(NeighborSelection input, int n)
    {
        if (input.AllNeighbors is null)
            throw new ArgumentNullException(nameof(input));
        if (n < 0)
            throw new ArgumentOutOfRangeException(nameof(n));
        if (input.AllNeighbors.Length != n)
            throw new ArgumentException("NeighborSelection row count must match node count.", nameof(n));
        if (n == 0)
            return input;

        CsrGraph graph = BuildWeightedGraph(input, n);
        var refinedRows = new Neighbor[n][];
        var nearestSample = new double[n];

        Parallel.For(0, n,
            () => new ThreadScratch(n),
            (i, _, scratch) =>
            {
                var row = input.AllNeighbors[i];
                int rowLength = row.Length;
                if (rowLength == 0)
                {
                    refinedRows[i] = Array.Empty<Neighbor>();
                    nearestSample[i] = double.PositiveInfinity;
                    return scratch;
                }

                var refinedRow = new Neighbor[rowLength];
                row.CopyTo(refinedRow, 0);
                refinedRows[i] = refinedRow;

                scratch.Reset(n, _maxDistance);
                foreach (var neighbor in refinedRow)
                    scratch.TargetMask[neighbor.Index] = true;

                Dijkstra.ComputeBoundedDistances(
                    graph,
                    i,
                    scratch.Distances.AsSpan(0, n),
                    scratch.Hops.AsSpan(0, n),
                    scratch.Queue,
                    scratch.TargetMask.AsSpan(0, n),
                    scratch.MaxDistance);

                for (int index = 0; index < rowLength; index++)
                {
                    int target = refinedRow[index].Index;
                    double original = refinedRow[index].Distance;
                    double refinedDistance = Math.Min(original, scratch.Distances[target]);
                    refinedRow[index].Distance = refinedDistance;
                }

                double best = double.PositiveInfinity;
                for (int j = 0; j < rowLength; j++)
                {
                    double d = refinedRow[j].Distance;
                    if (d < best) best = d;
                }
                nearestSample[i] = best;
                return scratch;
            },
            scratch => scratch.Dispose());

        return new NeighborSelection(refinedRows, nearestSample, input.KthNeighborDistances);
    }

    private static CsrGraph BuildWeightedGraph(NeighborSelection selection, int n)
    {
        int edgeCount = 0;
        for (int i = 0; i < n; i++)
            edgeCount += selection.AllNeighbors[i].Length;

        var rowPointers = new int[n + 1];
        var targets = new int[edgeCount];
        var weights = new double[edgeCount];

        int pos = 0;
        for (int i = 0; i < n; i++)
        {
            rowPointers[i] = pos;
            foreach (var neighbor in selection.AllNeighbors[i])
            {
                targets[pos] = neighbor.Index;
                weights[pos] = neighbor.Distance < 0.0 ? 0.0 : neighbor.Distance;
                pos++;
            }
        }
        rowPointers[n] = pos;

        return new CsrGraph
        {
            NodeCount = n,
            RowPointers = rowPointers,
            Targets = targets,
            Weights = weights,
        };
    }

    private sealed class ThreadScratch : IDisposable
    {
        public readonly double[] Distances;
        public readonly int[] Hops;
        public readonly bool[] TargetMask;
        public readonly PriorityQueue<int, double> Queue;
        public double MaxDistance;

        public ThreadScratch(int n)
        {
            Distances = ArrayPool<double>.Shared.Rent(n);
            Hops = ArrayPool<int>.Shared.Rent(n);
            TargetMask = ArrayPool<bool>.Shared.Rent(n);
            Queue = new PriorityQueue<int, double>();
            MaxDistance = double.PositiveInfinity;
        }

        public void Reset(int n, double? maxDistance)
        {
            MaxDistance = maxDistance ?? double.PositiveInfinity;
            Array.Fill(Distances, double.PositiveInfinity, 0, n);
            Array.Fill(Hops, int.MaxValue, 0, n);
            Array.Fill(TargetMask, false, 0, n);
            Queue.Clear();
        }

        public void Dispose()
        {
            ArrayPool<double>.Shared.Return(Distances);
            ArrayPool<int>.Shared.Return(Hops);
            ArrayPool<bool>.Shared.Return(TargetMask);
        }
    }
}
