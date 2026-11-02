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
    /// <summary><c>K_t ⊂ K_t ∪ K_{t+1} ⊃ K_{t+1}</c>: additions first. A transient is a cycle only the union holds.</summary>
    Union,
    /// <summary><c>K_t ⊃ K_t ∩ K_{t+1} ⊂ K_{t+1}</c>: deletions first. A transient is a component only the intersection holds.</summary>
    Intersection,
}

/// <summary>A zigzag bar read on the frames: its class exists at frames <paramref name="FirstFrame"/> through
/// <paramref name="LastFrame"/>, inclusive. <paramref name="StepBar"/> is the bar as the engine reported it.</summary>
public readonly record struct FrameBar(int Dimension, int FirstFrame, int LastFrame, Bar StepBar);

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

    /// <summary>The filtration to hand to an engine: <see cref="GraphZigzag"/>, <see cref="FastZigzag"/>, an oracle.</summary>
    public ZigzagFiltration Filtration { get; }

    public FrameTransition Transition { get; }

    public int FrameCount => _frameLevels.Length;

    internal GraphFrameFiltration(ZigzagFiltration filtration, FrameTransition transition, int[] frameLevels, int[] betweenLevels)
    {
        Filtration = filtration;
        Transition = transition;
        _frameLevels = frameLevels;
        _betweenLevels = betweenLevels;
    }

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
            if (a <= b)
                bars.Add(new FrameBar(bar.Dimension, a, b, bar));
        }

        return bars;
    }

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
    /// <summary>Builds the filtration for <paramref name="frames"/>. Vertex <c>v</c> is cell <c>v</c>; an edge keeps
    /// one cell id across re-entries. Faces enter before cofaces and leave after them, each batch in id order.</summary>
    /// <param name="n">Vertex id bound: every id is in <c>[0, n)</c>.</param>
    /// <param name="frames">The sequence.</param>
    /// <param name="transition">How consecutive frames connect.</param>
    public static GraphFrameFiltration Build(int n, IReadOnlyList<GraphFrame> frames, FrameTransition transition = FrameTransition.Union)
    {
        ArgumentNullException.ThrowIfNull(frames);
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

        var filtration = new ZigzagFiltration();
        var frameLevels = new int[frames.Count];
        var betweenLevels = new int[Math.Max(0, frames.Count - 1)];
        if (frames.Count == 0)
            return new GraphFrameFiltration(filtration, transition, frameLevels, betweenLevels);

        var edgeCell = new Dictionary<(int Lo, int Hi), int>();
        var vertices = new SortedSet<int>();
        var edges = new SortedSet<(int Lo, int Hi)>();

        for (int t = 0; t < frames.Count; t++)
        {
            var (nextVertices, nextEdges) = Normalize(frames[t], n, t);

            if (t == 0)
            {
                Enter(filtration, nextVertices, nextEdges, edgeCell, n);
            }
            else if (transition == FrameTransition.Union)
            {
                Enter(filtration, Except(nextVertices, vertices), Except(nextEdges, edges), edgeCell, n);
                betweenLevels[t - 1] = filtration.Count;
                Leave(filtration, Except(vertices, nextVertices), Except(edges, nextEdges), edgeCell);
            }
            else
            {
                Leave(filtration, Except(vertices, nextVertices), Except(edges, nextEdges), edgeCell);
                betweenLevels[t - 1] = filtration.Count;
                Enter(filtration, Except(nextVertices, vertices), Except(nextEdges, edges), edgeCell, n);
            }

            frameLevels[t] = filtration.Count;
            vertices = nextVertices;
            edges = nextEdges;
        }

        Leave(filtration, vertices, edges, edgeCell);
        return new GraphFrameFiltration(filtration, transition, frameLevels, betweenLevels);
    }

    /// <summary><see cref="Build"/> for frames that all hold every vertex <c>0..n-1</c>.</summary>
    public static GraphFrameFiltration FromEdgeSets(
        int n,
        IReadOnlyList<IReadOnlyList<(int i, int j)>> edgeSets,
        FrameTransition transition = FrameTransition.Union)
    {
        ArgumentNullException.ThrowIfNull(edgeSets);
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

        var all = new int[n];
        for (int v = 0; v < n; v++) all[v] = v;

        var frames = new GraphFrame[edgeSets.Count];
        for (int t = 0; t < frames.Length; t++)
            frames[t] = new GraphFrame(all, edgeSets[t]);

        return Build(n, frames, transition);
    }

    /// <summary><see cref="Build"/> for a sequence of graphs on one vertex set; weights are not read.</summary>
    public static GraphFrameFiltration FromGraphs(IReadOnlyList<CsrGraph> graphs, FrameTransition transition = FrameTransition.Union)
    {
        ArgumentNullException.ThrowIfNull(graphs);
        if (graphs.Count == 0) return FromEdgeSets(0, Array.Empty<IReadOnlyList<(int i, int j)>>(), transition);

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

        return FromEdgeSets(n, edgeSets, transition);
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
    public static GraphFrameFiltration FromEdgeField(
        CsrGraph substrate,
        IReadOnlyList<double[]> fields,
        IReadOnlyList<double> thresholds,
        FrameTransition transition = FrameTransition.Union)
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

        return FromEdgeSets(n, edgeSets, transition);
    }

    static (SortedSet<int> Vertices, SortedSet<(int Lo, int Hi)> Edges) Normalize(GraphFrame frame, int n, int t)
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

        return (vertices, edges);
    }

    static void Enter(
        ZigzagFiltration filtration,
        IEnumerable<int> vertices,
        IEnumerable<(int Lo, int Hi)> edges,
        Dictionary<(int Lo, int Hi), int> edgeCell,
        int n)
    {
        foreach (int v in vertices)
            filtration.Add(v, Array.Empty<int>());

        foreach (var edge in edges)
        {
            if (!edgeCell.TryGetValue(edge, out int cell))
            {
                cell = n + edgeCell.Count;
                edgeCell[edge] = cell;
            }
            filtration.Add(cell, new[] { edge.Lo, edge.Hi });
        }
    }

    static void Leave(
        ZigzagFiltration filtration,
        IEnumerable<int> vertices,
        IEnumerable<(int Lo, int Hi)> edges,
        Dictionary<(int Lo, int Hi), int> edgeCell)
    {
        foreach (var edge in edges)
            filtration.Delete(edgeCell[edge]);

        foreach (int v in vertices)
            filtration.Delete(v);
    }

    static List<T> Except<T>(SortedSet<T> keep, SortedSet<T> minus)
    {
        var result = new List<T>();
        foreach (T x in keep)
            if (!minus.Contains(x))
                result.Add(x);
        return result;
    }
}
