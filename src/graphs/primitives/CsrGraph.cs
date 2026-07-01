using System;
using System.IO;

namespace Graphs.Primitives
{
    /// <summary>
    /// Weighted undirected edge between two data points. <c>J</c> is the coupling
    /// strength (kernel-transformed distance), not the raw distance. Lives here
    /// because <see cref="CsrGraph.FromEdges"/> is its primary consumer; any
    /// builder that produces a graph passes an array of these to be densified
    /// into CSR form.
    /// </summary>
    public struct Edge
    {
        public int Source;
        public int Target;
        public double J;

        public Edge(int source, int target, double coupling)
        {
            Source = source;
            Target = target;
            J = coupling;
        }
    }

    /// <summary>
    /// Symmetric CSR adjacency graph. Each undirected edge is stored in both
    /// endpoint rows. The Swendsen-Wang inner loop skips target &lt;= source to
    /// process each edge exactly once. Also used as the shared adjacency
    /// structure for TDA passes (LocalTangent, BFS orientation propagation).
    /// </summary>
    public struct CsrGraph
    {
        public int[] Targets;
        public double[] Weights;
        public int[] RowPointers;
        public int NodeCount;

        /// <summary>
        /// Build a symmetric CSR graph from an undirected Edge[] where each
        /// edge appears once (source &lt; target is not required, but duplicates
        /// must not exist).
        /// </summary>
        public static CsrGraph FromEdges(Edge[] edges, int nodeCount)
        {
            int[] degree = new int[nodeCount];
            for (int e = 0; e < edges.Length; e++)
            {
                degree[edges[e].Source]++;
                degree[edges[e].Target]++;
            }

            int[] rowPointers = new int[nodeCount + 1];
            for (int node = 0; node < nodeCount; node++)
                rowPointers[node + 1] = rowPointers[node] + degree[node];

            int totalEntries = rowPointers[nodeCount];
            int[] targets = new int[totalEntries];
            double[] weights = new double[totalEntries];

            int[] cursor = new int[nodeCount];
            Array.Copy(rowPointers, cursor, nodeCount);

            for (int e = 0; e < edges.Length; e++)
            {
                int src = edges[e].Source;
                int tgt = edges[e].Target;
                double w = edges[e].J;

                targets[cursor[src]] = tgt;
                weights[cursor[src]] = w;
                cursor[src]++;

                targets[cursor[tgt]] = src;
                weights[cursor[tgt]] = w;
                cursor[tgt]++;
            }

            return new CsrGraph
            {
                Targets = targets,
                Weights = weights,
                RowPointers = rowPointers,
                NodeCount = nodeCount
            };
        }

        /// <summary>
        /// Unweighted degree of a node — the count of CSR row entries for it.
        /// Replaces the inline <c>RowPointers[node+1] - RowPointers[node]</c>
        /// idiom scattered across diagnostics and the LMP post-filter.
        /// </summary>
        public int Degree(int node) => RowPointers[node + 1] - RowPointers[node];

        /// <summary>
        /// Sum of edge weights incident to a node. The Laplacian-diagonal
        /// quantity D[v,v]; consumed by <c>AlgebraicConnectivity</c>'s
        /// Laplacian construction and by any other diagnostic that wants
        /// volume-style summaries.
        /// </summary>
        public double WeightedDegree(int node)
        {
            double sum = 0.0;
            int start = RowPointers[node];
            int end   = RowPointers[node + 1];
            for (int e = start; e < end; e++) sum += Weights[e];
            return sum;
        }

        /// <summary>
        /// CSR slot of the directed edge <c>(row → target)</c>; throws if absent
        /// (an asymmetric CSR, which this type does not produce).
        /// </summary>
        public int FindSlot(int row, int target)
        {
            int rowEnd = RowPointers[row + 1];
            for (int e = RowPointers[row]; e < rowEnd; e++)
                if (Targets[e] == target) return e;
            throw new InvalidOperationException(
                $"No CSR slot for directed edge ({row}→{target}); CSR is not symmetric.");
        }

        /// <summary>
        /// For each directed CSR slot, the slot of its reverse direction — the
        /// mirror map of a symmetric CSR. Lets a consumer reconcile directed
        /// per-edge fields (e.g. a per-site field whose two endpoints disagree)
        /// without re-walking the graph each time. Symmetric CSR required.
        /// </summary>
        public int[] BuildReverseSlotMap()
        {
            var mirror = new int[Targets.Length];
            for (int i = 0; i < NodeCount; i++)
            {
                int rowEnd = RowPointers[i + 1];
                for (int e = RowPointers[i]; e < rowEnd; e++)
                    mirror[e] = FindSlot(Targets[e], i);
            }
            return mirror;
        }

        /// <summary>
        /// Persist this CSR graph in a compact binary form.
        /// </summary>
        public void WriteTo(BinaryWriter writer)
        {
            if (writer is null) throw new ArgumentNullException(nameof(writer));
            writer.Write(NodeCount);
            writer.Write(RowPointers.Length);
            writer.Write(Targets.Length);
            writer.Write(Weights.Length);

            for (int i = 0; i < RowPointers.Length; i++)
                writer.Write(RowPointers[i]);

            for (int i = 0; i < Targets.Length; i++)
                writer.Write(Targets[i]);

            for (int i = 0; i < Weights.Length; i++)
                writer.Write(Weights[i]);
        }

        /// <summary>
        /// Rehydrate a CSR graph from the binary format written by <see cref="WriteTo"/>.
        /// </summary>
        public static CsrGraph FromBinary(BinaryReader reader)
        {
            if (reader is null) throw new ArgumentNullException(nameof(reader));
            int nodeCount = reader.ReadInt32();
            int rowPointersLength = reader.ReadInt32();
            int targetsLength = reader.ReadInt32();
            int weightsLength = reader.ReadInt32();

            if (rowPointersLength != nodeCount + 1)
                throw new InvalidDataException(
                    $"Corrupt CSR graph: expected RowPointers length {nodeCount + 1}, got {rowPointersLength}.");
            if (targetsLength != weightsLength)
                throw new InvalidDataException(
                    $"Corrupt CSR graph: Targets length ({targetsLength}) must equal Weights length ({weightsLength}).");

            int[] rowPointers = new int[rowPointersLength];
            for (int i = 0; i < rowPointersLength; i++)
                rowPointers[i] = reader.ReadInt32();

            int[] targets = new int[targetsLength];
            for (int i = 0; i < targetsLength; i++)
                targets[i] = reader.ReadInt32();

            double[] weights = new double[weightsLength];
            for (int i = 0; i < weightsLength; i++)
                weights[i] = reader.ReadDouble();

            return new CsrGraph
            {
                NodeCount = nodeCount,
                RowPointers = rowPointers,
                Targets = targets,
                Weights = weights,
            };
        }

        /// <summary>
        /// Build a new <see cref="CsrGraph"/> on the subset of nodes selected by
        /// <paramref name="nodeMask"/>. Only edges whose <em>both</em> endpoints
        /// are masked are retained; edge weights are copied verbatim. Output
        /// node indices are densely renumbered in <c>[0, subN)</c>; use
        /// <paramref name="newToOld"/> to translate back, or
        /// <paramref name="oldToNew"/> for the reverse direction
        /// (<c>-1</c> for nodes that were filtered out).
        /// </summary>
        /// <remarks>
        /// O(N + E) time and space. The MAPPER-SPC patch driver consumes this:
        /// given a Mapper cover bin's preimage, restrict the full proximity
        /// graph to that bin's nodes and hand the result to the adaptive SPC
        /// scheduler. Output preserves the input's CSR symmetry — if both
        /// (i,j) and (j,i) survive the mask, both entries appear in the result.
        /// </remarks>
        public CsrGraph InducedSubgraph(
            bool[] nodeMask,
            out int[] newToOld,
            out int[] oldToNew)
        {
            if (nodeMask is null)
                throw new ArgumentNullException(nameof(nodeMask));
            if (nodeMask.Length != NodeCount)
                throw new ArgumentException(
                    $"nodeMask.Length ({nodeMask.Length}) must equal NodeCount ({NodeCount}).",
                    nameof(nodeMask));

            // Pass 1: build dense index maps. New node ids are assigned in the
            // order the mask is scanned, so the output indexing is stable and
            // deterministic for a given mask.
            oldToNew = new int[NodeCount];
            int subN = 0;
            for (int v = 0; v < NodeCount; v++)
                oldToNew[v] = nodeMask[v] ? subN++ : -1;

            newToOld = new int[subN];
            for (int v = 0; v < NodeCount; v++)
                if (nodeMask[v]) newToOld[oldToNew[v]] = v;

            if (subN == 0)
            {
                return new CsrGraph
                {
                    Targets     = Array.Empty<int>(),
                    Weights     = Array.Empty<double>(),
                    RowPointers = new int[1],
                    NodeCount   = 0,
                };
            }

            // Pass 2: count induced degree per new node (i.e., count surviving
            // edges per row before we know where each row starts).
            int[] subDegree = new int[subN];
            for (int u = 0; u < NodeCount; u++)
            {
                if (!nodeMask[u]) continue;
                int newU = oldToNew[u];
                int rowEnd = RowPointers[u + 1];
                for (int e = RowPointers[u]; e < rowEnd; e++)
                    if (nodeMask[Targets[e]]) subDegree[newU]++;
            }

            // Pass 3: prefix-sum degrees into RowPointers.
            int[] subRowPointers = new int[subN + 1];
            for (int i = 0; i < subN; i++)
                subRowPointers[i + 1] = subRowPointers[i] + subDegree[i];

            int totalEntries = subRowPointers[subN];
            int[] subTargets  = new int[totalEntries];
            double[] subWeights = new double[totalEntries];

            // Pass 4: fill targets + weights using a per-row write cursor.
            int[] cursor = new int[subN];
            for (int u = 0; u < NodeCount; u++)
            {
                if (!nodeMask[u]) continue;
                int newU = oldToNew[u];
                int rowEnd = RowPointers[u + 1];
                for (int e = RowPointers[u]; e < rowEnd; e++)
                {
                    int v = Targets[e];
                    if (!nodeMask[v]) continue;
                    int pos = subRowPointers[newU] + cursor[newU]++;
                    subTargets[pos] = oldToNew[v];
                    subWeights[pos] = Weights[e];
                }
            }

            return new CsrGraph
            {
                Targets     = subTargets,
                Weights     = subWeights,
                RowPointers = subRowPointers,
                NodeCount   = subN,
            };
        }
    }
}
