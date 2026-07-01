using System.Runtime.CompilerServices;

namespace Graphs.Primitives;

/// <summary>One undirected edge of a symmetric CSR graph: its endpoints and the
/// CSR slot of its <c>j &gt; i</c> half (the slot every per-edge field is keyed by).</summary>
public readonly struct UndirectedEdge
{
    public readonly int Source;
    public readonly int Target;
    public readonly int Slot;

    public UndirectedEdge(int source, int target, int slot)
    {
        Source = source;
        Target = target;
        Slot = slot;
    }
}

/// <summary>
/// Zero-allocation walk over each undirected edge of a symmetric CSR graph
/// exactly once (the <c>j &gt; i</c> upper-triangular half) — the single home of
/// the <c>for i; for e in row; if Targets[e] &lt;= i continue</c> idiom hand-rolled
/// across the codebase. A <see langword="ref"/> <see langword="struct"/> so it
/// cannot escape to the heap; <c>MoveNext</c>/<c>Current</c> are inlined, so a
/// <c>foreach</c> over it lowers to the same nested loop with no alloc and no
/// delegate.
/// </summary>
/// <remarks>
/// Adoption note: safe wherever the walk is the outer structure. **Hot fused
/// loops (SW <c>Draw</c>, PKWang <c>BuildHcum</c>/<c>Solve</c>) keep their
/// hand-rolled loops until a benchmark confirms this lowers identically** — the
/// loop is intricate and perf-critical, and "zero-cost" is asserted here, not yet
/// measured.
/// </remarks>
public readonly ref struct UndirectedEdgeWalk
{
    private readonly int[] _rowPtr;
    private readonly int[] _targets;
    private readonly int _nodeCount;

    internal UndirectedEdgeWalk(int[] rowPtr, int[] targets, int nodeCount)
    {
        _rowPtr = rowPtr;
        _targets = targets;
        _nodeCount = nodeCount;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public Enumerator GetEnumerator() => new Enumerator(_rowPtr, _targets, _nodeCount);

    public ref struct Enumerator
    {
        private readonly int[] _rowPtr;
        private readonly int[] _targets;
        private readonly int _nodeCount;
        private int _i;
        private int _e;
        private int _rowEnd;
        private UndirectedEdge _current;

        internal Enumerator(int[] rowPtr, int[] targets, int nodeCount)
        {
            _rowPtr = rowPtr;
            _targets = targets;
            _nodeCount = nodeCount;
            _i = -1;
            _e = 0;
            _rowEnd = 0;
            _current = default;
        }

        public readonly UndirectedEdge Current => _current;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool MoveNext()
        {
            while (true)
            {
                if (_e >= _rowEnd)
                {
                    _i++;
                    if (_i >= _nodeCount) return false;
                    _e = _rowPtr[_i];
                    _rowEnd = _rowPtr[_i + 1];
                    continue;
                }

                int slot = _e++;
                int j = _targets[slot];
                if (j <= _i) continue;
                _current = new UndirectedEdge(_i, j, slot);
                return true;
            }
        }
    }
}

/// <summary>Extension entry point: <c>foreach (var edge in graph.UndirectedEdges())</c>.</summary>
public static class CsrGraphEdgeExtensions
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static UndirectedEdgeWalk UndirectedEdges(this CsrGraph graph)
        => new UndirectedEdgeWalk(graph.RowPointers, graph.Targets, graph.NodeCount);
}
