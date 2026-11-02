#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace TDA.Ph;

/// <summary>One graph of a sequence over stable vertex ids: the vertices present and its undirected edges,
/// each edge's endpoints among them. Self-loops are dropped and duplicates collapse.</summary>
public sealed record GraphFrame(IReadOnlyList<int> Vertices, IReadOnlyList<(int i, int j)> Edges);

/// <summary>
/// How consecutive frames connect in the zigzag. The two agree on every frame (Carlsson and de Silva's Diamond
/// Principle) and differ only in the transients between frames, one degree apart.
/// </summary>
public enum FrameTransition
{
    /// <summary><c>K_t ⊂ K_t ∪ K_{t+1} ⊃ K_{t+1}</c>: additions first. A transient is a class only the union holds: a
    /// cycle, or under <see cref="FrameComplex.Flag"/> also a void.</summary>
    Union,
    /// <summary><c>K_t ⊃ K_t ∩ K_{t+1} ⊂ K_{t+1}</c>: deletions first. A transient is a class only the intersection
    /// holds: a component, or under <see cref="FrameComplex.Flag"/> also a cycle.</summary>
    Intersection,
}

/// <summary>Which complex each frame contributes.</summary>
public enum FrameComplex
{
    /// <summary>The graph itself. A cycle lasts until an edge leaves it. Read by <see cref="GraphZigzag"/> or
    /// <see cref="FastZigzag"/>.</summary>
    Graph,
    /// <summary>The graph's flag complex through triangles. A cycle also dies when triangles fill it, and a hole opens
    /// when they leave. Read by <see cref="FastZigzag"/> with <c>maxDimension: 1</c>: H₀ and H₁ are exact, while H₂ has
    /// nothing above to fill it. Under <see cref="FrameTransition.Union"/> the between-frames complex is the union of
    /// the two flag complexes, so it holds no triangle that neither frame holds.</summary>
    Flag,
}

/// <summary>What starts or ends a class between two frames. It does not depend on the
/// <see cref="FrameTransition"/>.</summary>
public enum FrameEvent
{
    /// <summary>The bar is at the sequence's first or last frame, so nothing observed starts or ends it there.</summary>
    Unobserved,
    /// <summary>Cells arrive. A birth: a vertex enters, or an edge closes a cycle. A death: a merge, or triangles fill
    /// a cycle.</summary>
    Arrival,
    /// <summary>Cells leave. A birth: an edge leaves and splits a component, or triangles leave and open a hole. A
    /// death: an edge leaves and breaks a cycle, or a vertex leaves.</summary>
    Departure,
}

/// <summary>A zigzag bar read on the frames: its class exists at frames <paramref name="FirstFrame"/> through
/// <paramref name="LastFrame"/>, inclusive. <paramref name="Start"/> and <paramref name="End"/> say what began and ended
/// it, so a cycle filled and a cycle broken read apart. <paramref name="StepBar"/> is the bar as the engine reported it.</summary>
public readonly record struct FrameBar(int Dimension, int FirstFrame, int LastFrame, FrameEvent Start, FrameEvent End, Bar StepBar);

/// <summary>A zigzag bar present between frames <paramref name="AfterFrame"/> and <paramref name="AfterFrame"/> + 1
/// and at no frame. <paramref name="StepBar"/> is the bar as the engine reported it.</summary>
public readonly record struct TransientBar(int Dimension, int AfterFrame, Bar StepBar);

/// <summary>
/// A frame sequence as a <see cref="ZigzagFiltration"/>, one cell per step from empty back to empty, with the
/// level of each frame so a barcode reads back on them. A cell's step order inside one batch is an artifact of
/// the refinement, so a bar that reaches neither a frame nor a between-frames complex is dropped by the views.
/// </summary>
public sealed class GraphFrameFiltration
{
    readonly int[] _frameLevels;
    readonly int[] _betweenLevels;
    readonly int[][] _cellVertices;

    /// <summary>The filtration to hand to an engine: <see cref="GraphZigzag"/>, <see cref="FastZigzag"/>, an oracle.</summary>
    public ZigzagFiltration Filtration { get; }

    public FrameTransition Transition { get; }

    public FrameComplex Complex { get; }

    public int FrameCount => _frameLevels.Length;

    internal GraphFrameFiltration(ZigzagFiltration filtration, FrameTransition transition, FrameComplex complex, int[] frameLevels, int[] betweenLevels, int[][] cellVertices)
    {
        Filtration = filtration;
        Transition = transition;
        Complex = complex;
        _frameLevels = frameLevels;
        _betweenLevels = betweenLevels;
        _cellVertices = cellVertices;
    }

    /// <summary>The vertices of a cell id in <see cref="Filtration"/>, so a bar's <see cref="Bar.Cycle"/> reads in frame
    /// terms: a vertex's own id, an edge's two, a triangle's three. That cycle belongs to the step complex where the bar
    /// begins, which can still hold cells that leave before the bar's first frame.</summary>
    public IReadOnlyList<int> VerticesOf(int cellId) => _cellVertices[cellId];

    /// <summary>The bars that exist at some frame, read on the frames. A view over the step-level barcode: it
    /// drops the transients (<see cref="TransientBars"/>) and the step-order artifacts.</summary>
    public IReadOnlyList<FrameBar> FrameBars(Barcode barcode)
    {
        ArgumentNullException.ThrowIfNull(barcode);

        var bars = new List<FrameBar>();
        foreach (Bar bar in barcode.Bars)
        {
            var (first, last) = LevelsOf(bar);
            int a = LowerBound(_frameLevels, first);
            int b = LowerBound(_frameLevels, last + 1) - 1;
            if (a > b) continue;

            FrameEvent start = a == 0 ? FrameEvent.Unobserved : EventOf(Filtration[first - 1]);
            FrameEvent end = b == FrameCount - 1 ? FrameEvent.Unobserved : EventOf(Filtration[last]);
            bars.Add(new FrameBar(bar.Dimension, a, b, start, end, bar));
        }

        return bars;
    }

    static FrameEvent EventOf(ZigzagStep step) =>
        step.Direction == ZigzagDirection.Add ? FrameEvent.Arrival : FrameEvent.Departure;

    /// <summary>The bars present between two frames and at neither. A class that reaches a between-frames complex
    /// without a frame covers exactly one, since frames separate them.</summary>
    public IReadOnlyList<TransientBar> TransientBars(Barcode barcode)
    {
        ArgumentNullException.ThrowIfNull(barcode);

        var bars = new List<TransientBar>();
        foreach (Bar bar in barcode.Bars)
        {
            var (first, last) = LevelsOf(bar);
            if (LowerBound(_frameLevels, first) <= LowerBound(_frameLevels, last + 1) - 1)
                continue;

            int t = LowerBound(_betweenLevels, first);
            if (t < _betweenLevels.Length && _betweenLevels[t] <= last)
                bars.Add(new TransientBar(bar.Dimension, t, bar));
        }

        return bars;
    }

    // Level ℓ is the complex after ℓ steps. A bar reports its creating step (level − 1) and its last level,
    // by ZigzagBarcodeNaive's convention, so it covers levels [Birth + 1, Death].
    static (int First, int Last) LevelsOf(in Bar bar) => ((int)bar.Birth + 1, (int)bar.Death);

    static int LowerBound(int[] sorted, int value)
    {
        int index = Array.BinarySearch(sorted, value);
        if (index < 0) return ~index;
        while (index > 0 && sorted[index - 1] == value) index--;
        return index;
    }
}

/// <summary>
/// The input boundary of the zigzag engines for a sequence of graphs: frames over stable vertex ids become a
/// <see cref="ZigzagFiltration"/>, read by the engines unchanged, the core staying free of <c>Graphs.*</c>. A frame
/// sequence is whatever a caller's sweep produces: graphs at successive parameter values, a sliding window, a
/// thresholded field along any path. Nested frames reduce to ordinary persistence; frames that lose edges are what
/// needs the zigzag.
/// </summary>
public static class GraphFrameZigzag
{
    /// <summary>Builds the filtration for <paramref name="frames"/>. Vertex <c>v</c> is cell <c>v</c>; an edge or a
    /// triangle keeps one cell id across re-entries. Faces enter before cofaces and leave after them, each batch in id
    /// order.</summary>
    /// <param name="n">Vertex id bound: every id is in <c>[0, n)</c>.</param>
    /// <param name="frames">The sequence.</param>
    /// <param name="transition">How consecutive frames connect.</param>
    /// <param name="complex">Which complex each frame contributes.</param>
    public static GraphFrameFiltration Build(
        int n,
        IReadOnlyList<GraphFrame> frames,
        FrameTransition transition = FrameTransition.Union,
        FrameComplex complex = FrameComplex.Graph)
    {
        ArgumentNullException.ThrowIfNull(frames);
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

        var filtration = new ZigzagFiltration();
        var frameLevels = new int[frames.Count];
        var betweenLevels = new int[Math.Max(0, frames.Count - 1)];
        var ids = new CellIds(n);
        if (frames.Count == 0)
            return new GraphFrameFiltration(filtration, transition, complex, frameLevels, betweenLevels, ids.CellVertices());

        var cells = FrameCells.Empty;

        for (int t = 0; t < frames.Count; t++)
        {
            FrameCells next = Normalize(frames[t], n, t, complex);

            if (t == 0)
            {
                Enter(filtration, next, ids);
            }
            else if (transition == FrameTransition.Union)
            {
                Enter(filtration, next.Except(cells), ids);
                betweenLevels[t - 1] = filtration.Count;
                Leave(filtration, cells.Except(next), ids);
            }
            else
            {
                Leave(filtration, cells.Except(next), ids);
                betweenLevels[t - 1] = filtration.Count;
                Enter(filtration, next.Except(cells), ids);
            }

            frameLevels[t] = filtration.Count;
            cells = next;
        }

        Leave(filtration, cells, ids);
        return new GraphFrameFiltration(filtration, transition, complex, frameLevels, betweenLevels, ids.CellVertices());
    }

    /// <summary><see cref="Build"/> for frames that all hold every vertex <c>0..n-1</c>.</summary>
    public static GraphFrameFiltration FromEdgeSets(
        int n,
        IReadOnlyList<IReadOnlyList<(int i, int j)>> edgeSets,
        FrameTransition transition = FrameTransition.Union,
        FrameComplex complex = FrameComplex.Graph)
    {
        ArgumentNullException.ThrowIfNull(edgeSets);
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

        var all = new int[n];
        for (int v = 0; v < n; v++) all[v] = v;

        var frames = new GraphFrame[edgeSets.Count];
        for (int t = 0; t < frames.Length; t++)
            frames[t] = new GraphFrame(all, edgeSets[t]);

        return Build(n, frames, transition, complex);
    }

    /// <summary><see cref="Build"/> for a sequence of graphs on one vertex set; weights are not read.</summary>
    public static GraphFrameFiltration FromGraphs(
        IReadOnlyList<CsrGraph> graphs,
        FrameTransition transition = FrameTransition.Union,
        FrameComplex complex = FrameComplex.Graph)
    {
        ArgumentNullException.ThrowIfNull(graphs);
        if (graphs.Count == 0) return FromEdgeSets(0, Array.Empty<IReadOnlyList<(int i, int j)>>(), transition, complex);

        int n = graphs[0].NodeCount;
        var edgeSets = new IReadOnlyList<(int i, int j)>[graphs.Count];
        for (int t = 0; t < edgeSets.Length; t++)
        {
            CsrGraph g = graphs[t];
            if (g.NodeCount != n)
                throw new ArgumentException($"Frame {t} has {g.NodeCount} nodes against {n}.", nameof(graphs));

            var edges = new List<(int i, int j)>();
            for (int u = 0; u < n; u++)
                for (int e = g.RowPointers[u]; e < g.RowPointers[u + 1]; e++)
                    if (g.Targets[e] > u)
                        edges.Add((u, g.Targets[e]));
            edgeSets[t] = edges;
        }

        return FromEdgeSets(n, edgeSets, transition, complex);
    }

    /// <summary>
    /// <see cref="Build"/> for one substrate graph and a per-edge field per frame: frame <c>t</c> holds the
    /// substrate's edges with <c>fields[t][slot] ≥ thresholds[t]</c>, reading the <c>j &gt; i</c> slot that every
    /// per-edge currency is keyed by. For a sublevel admission, negate the field and the threshold.
    /// </summary>
    /// <param name="substrate">The graph every frame is a subgraph of.</param>
    /// <param name="fields">Per frame, a field parallel to <see cref="CsrGraph.Targets"/>.</param>
    /// <param name="thresholds">Per frame, the admission bound. A constant, a warped path and a non-monotone ridge
    /// are all the caller's choice; nesting is not assumed.</param>
    /// <param name="transition">How consecutive frames connect.</param>
    /// <param name="complex">Which complex each frame contributes.</param>
    public static GraphFrameFiltration FromEdgeField(
        CsrGraph substrate,
        IReadOnlyList<double[]> fields,
        IReadOnlyList<double> thresholds,
        FrameTransition transition = FrameTransition.Union,
        FrameComplex complex = FrameComplex.Graph)
    {
        ArgumentNullException.ThrowIfNull(fields);
        ArgumentNullException.ThrowIfNull(thresholds);
        if (fields.Count != thresholds.Count)
            throw new ArgumentException($"{fields.Count} fields against {thresholds.Count} thresholds.", nameof(thresholds));

        int n = substrate.NodeCount;
        var edgeSets = new IReadOnlyList<(int i, int j)>[fields.Count];
        for (int t = 0; t < edgeSets.Length; t++)
        {
            double[] field = fields[t];
            if (field.Length != substrate.Targets.Length)
                throw new ArgumentException($"Field {t} has {field.Length} slots against {substrate.Targets.Length}.", nameof(fields));

            var edges = new List<(int i, int j)>();
            for (int u = 0; u < n; u++)
                for (int e = substrate.RowPointers[u]; e < substrate.RowPointers[u + 1]; e++)
                    if (substrate.Targets[e] > u && field[e] >= thresholds[t])
                        edges.Add((u, substrate.Targets[e]));
            edgeSets[t] = edges;
        }

        return FromEdgeSets(n, edgeSets, transition, complex);
    }

    // One frame's cells by dimension, ordered so each batch enters and leaves in a fixed order.
    sealed record FrameCells(SortedSet<int> Vertices, SortedSet<(int Lo, int Hi)> Edges, SortedSet<(int A, int B, int C)> Triangles)
    {
        public static FrameCells Empty => new(new(), new(), new());

        public FrameCells Except(FrameCells other) => new(
            Minus(Vertices, other.Vertices),
            Minus(Edges, other.Edges),
            Minus(Triangles, other.Triangles));

        static SortedSet<T> Minus<T>(SortedSet<T> keep, SortedSet<T> minus)
        {
            var result = new SortedSet<T>(keep);
            result.ExceptWith(minus);
            return result;
        }
    }

    // Vertex v is cell v; edges and triangles take the next free id at first appearance and keep it.
    sealed class CellIds(int vertexCount)
    {
        readonly Dictionary<(int Lo, int Hi), int> _edges = new();
        readonly Dictionary<(int A, int B, int C), int> _triangles = new();

        int Next => vertexCount + _edges.Count + _triangles.Count;

        public int Edge((int Lo, int Hi) edge)
        {
            if (!_edges.TryGetValue(edge, out int id))
                _edges[edge] = id = Next;
            return id;
        }

        public int Triangle((int A, int B, int C) triangle)
        {
            if (!_triangles.TryGetValue(triangle, out int id))
                _triangles[triangle] = id = Next;
            return id;
        }

        public int[][] CellVertices()
        {
            var vertices = new int[Next][];
            for (int v = 0; v < vertexCount; v++)
                vertices[v] = new[] { v };
            foreach (var (edge, id) in _edges)
                vertices[id] = new[] { edge.Lo, edge.Hi };
            foreach (var (triangle, id) in _triangles)
                vertices[id] = new[] { triangle.A, triangle.B, triangle.C };
            return vertices;
        }
    }

    static FrameCells Normalize(GraphFrame frame, int n, int t, FrameComplex complex)
    {
        ArgumentNullException.ThrowIfNull(frame);

        var vertices = new SortedSet<int>();
        foreach (int v in frame.Vertices)
        {
            if ((uint)v >= (uint)n)
                throw new ArgumentOutOfRangeException(nameof(frame), $"Frame {t}: vertex {v} out of range for n={n}.");
            vertices.Add(v);
        }

        var edges = new SortedSet<(int Lo, int Hi)>();
        foreach (var (i, j) in frame.Edges)
        {
            if (i == j) continue;
            if (!vertices.Contains(i) || !vertices.Contains(j))
                throw new ArgumentException($"Frame {t}: edge ({i},{j}) has an endpoint the frame does not hold.", nameof(frame));
            edges.Add(i < j ? (i, j) : (j, i));
        }

        var triangles = new SortedSet<(int A, int B, int C)>();
        if (complex == FrameComplex.Flag && edges.Count > 0)
        {
            var src = new int[edges.Count];
            var dst = new int[edges.Count];
            int e = 0;
            foreach (var (lo, hi) in edges)
            {
                src[e] = lo;
                dst[e++] = hi;
            }

            int[] flat = FlagComplex.FromEdges(src, dst);
            for (int f = 0; f < flat.Length; f += 3)
                triangles.Add((flat[f], flat[f + 1], flat[f + 2]));
        }

        return new FrameCells(vertices, edges, triangles);
    }

    static void Enter(ZigzagFiltration filtration, FrameCells cells, CellIds ids)
    {
        foreach (int v in cells.Vertices)
            filtration.Add(v, Array.Empty<int>());

        foreach (var edge in cells.Edges)
            filtration.Add(ids.Edge(edge), new[] { edge.Lo, edge.Hi });

        foreach (var (a, b, c) in cells.Triangles)
            filtration.Add(ids.Triangle((a, b, c)), new[] { ids.Edge((a, b)), ids.Edge((a, c)), ids.Edge((b, c)) });
    }

    static void Leave(ZigzagFiltration filtration, FrameCells cells, CellIds ids)
    {
        foreach (var triangle in cells.Triangles)
            filtration.Delete(ids.Triangle(triangle));

        foreach (var edge in cells.Edges)
            filtration.Delete(ids.Edge(edge));

        foreach (int v in cells.Vertices)
            filtration.Delete(v);
    }
}
