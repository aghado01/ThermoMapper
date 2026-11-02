#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Maths.Topology;

namespace TDA.Ph;

/// <summary>
/// P0 of the Conditioned-Persistence synthesis — the undirected, single-parameter (δ)
/// content core. Unions a backbone prior <c>K₀</c> (present from filtration 0) with
/// similarity/content edges carried at their metric distance, then reads the result
/// through the existing Rips → Z₂ engine. A similarity chord that closes a path along
/// the backbone is an H₁ generator born at that chord's distance — a <em>return relative
/// to the prior</em>. SIFTS is the <c>τ≡0</c> degenerate (reading-order backbone, zero prior).
/// <para>The union convention (<see cref="BuildGraph"/>) is the only new code here; everything
/// downstream is reuse — <see cref="RipsFiltration.GraphRips"/> and
/// <see cref="PersistentInvolutedHomology"/>. Two readings of one filtration:
/// <see cref="ComputeBarcode"/> is backbone-at-base <c>H(K_δ)</c>, which equals the relative one when
/// <c>K₀</c> is connected, spanning and acyclic; <see cref="ComputeRelativeBarcode"/> is
/// backbone-quotiented <c>H(K_δ, K₀)</c> for any backbone. A non-monotone sweep reads through
/// <see cref="BuildSweep"/>. Directedness, the gauge <c>q</c> and the sheaf / <c>λ_q</c> reading are
/// P2–P4 and deliberately absent.</para>
/// </summary>
public static class ConditionedFiltration
{
    /// <summary>
    /// Union the backbone prior (at <paramref name="baseWeight"/> = ε₀, present from
    /// filtration 0 so H₀ is anchored) with the similarity/content edges (at their distance,
    /// the δ scale) into one symmetric distance-weighted <see cref="CsrGraph"/>. The single
    /// new idea of P0 — the rest of the pipeline is the built engine.
    /// </summary>
    /// <param name="n">Vertex count.</param>
    /// <param name="backbone">Prior <c>K₀</c> edges <c>(i, j)</c>; admitted at ε₀ = <paramref name="baseWeight"/>.</param>
    /// <param name="similarity">Content edges <c>(i, j, d)</c>; admitted at distance <c>d</c>.</param>
    /// <param name="baseWeight">ε₀, the backbone birth value (default 0). A pair present in
    /// both sets takes the smaller weight — the backbone's ε₀ wins, so a coincident edge is
    /// anchored, not re-born at its distance.</param>
    public static CsrGraph BuildGraph(
        int n,
        IReadOnlyList<(int i, int j)> backbone,
        IReadOnlyList<(int i, int j, double d)> similarity,
        double baseWeight = 0.0)
    {
        ArgumentNullException.ThrowIfNull(backbone);
        ArgumentNullException.ThrowIfNull(similarity);
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
        if (baseWeight < 0.0)
            throw new ArgumentOutOfRangeException(nameof(baseWeight),
                "ε₀ must be non-negative — Rips filtration values are distances.");

        // Accumulate as (lo, hi) → min weight so duplicates within and across the two edge
        // sets collapse deterministically and the earliest birth wins (ε₀ ≤ d).
        var weights = new Dictionary<(int Lo, int Hi), double>(backbone.Count + similarity.Count);

        foreach (var (i, j) in backbone)
            Accumulate(weights, i, j, baseWeight, n);

        foreach (var (i, j, d) in similarity)
        {
            if (d < 0.0)
                throw new ArgumentException("Similarity edge distances must be non-negative.", nameof(similarity));
            Accumulate(weights, i, j, d, n);
        }

        var edges = new Edge[weights.Count];
        int e = 0;
        foreach (var ((lo, hi), w) in weights)
            edges[e++] = new Edge(lo, hi, w);

        return CsrGraph.FromEdges(edges, n);
    }

    static void Accumulate(Dictionary<(int Lo, int Hi), double> weights, int i, int j, double w, int n)
    {
        if (i == j) return; // no self-loops in a filtration skeleton
        if ((uint)i >= (uint)n || (uint)j >= (uint)n)
            throw new ArgumentOutOfRangeException(nameof(i), $"Edge ({i},{j}) out of range for n={n}.");

        (int Lo, int Hi) key = i < j ? (i, j) : (j, i);
        weights[key] = weights.TryGetValue(key, out double existing) ? Math.Min(existing, w) : w;
    }

    /// <summary>
    /// Convenience: <see cref="BuildGraph"/> → <see cref="RipsFiltration.GraphRips"/> (raw
    /// distance filtration values) → involuted persistence → <see cref="TDA.Ph.Barcode"/>. H₁ bars
    /// carry representative cycles, so <see cref="BarCycleEdges.H1Loops"/> recovers which chord
    /// closed each return. Reads the built engine unchanged — no reduction code lives here.
    /// </summary>
    public static Barcode ComputeBarcode(
        int n,
        IReadOnlyList<(int i, int j)> backbone,
        IReadOnlyList<(int i, int j, double d)> similarity,
        int maxDimension = 2,
        double baseWeight = 0.0)
    {
        CsrGraph graph = BuildGraph(n, backbone, similarity, baseWeight);
        SimplicialFiltration filtration = RipsFiltration.GraphRips(
            graph, FiltrationWeights.RawDistance, maxDimension, label: "conditioned");
        return PersistentInvolutedHomology.Compute(filtration, representatives: true);
    }

    /// <summary>
    /// The backbone-quotiented reading: the union's filtration relative to <c>K₀</c>, the flag complex of the
    /// backbone edges on their endpoints. <see cref="ComputeBarcode"/> reads <c>H(K_δ)</c>, which also reports
    /// <c>K₀</c>'s own cycles; this reads <c>H(K_δ, K₀)</c>, the returns the prior does not already hold. A
    /// content edge is outside <c>K₀</c> whatever its value, so a return the prior predicted (residual 0) stays.
    /// <c>K₀</c> is present throughout the pair, so ε₀ does not enter.
    /// </summary>
    /// <param name="n">Vertex count.</param>
    /// <param name="backbone">Prior <c>K₀</c> edges <c>(i, j)</c>. Empty gives <c>K₀ = ∅</c>, the absolute complex.</param>
    /// <param name="similarity">Content edges <c>(i, j, d)</c>; a pair also in the backbone belongs to <c>K₀</c>.</param>
    /// <param name="maxDimension">Skeleton dimension, as in <see cref="RipsFiltration.GraphRips"/>.</param>
    public static RelativeFiltration BuildRelative(
        int n,
        IReadOnlyList<(int i, int j)> backbone,
        IReadOnlyList<(int i, int j, double d)> similarity,
        int maxDimension = 2)
    {
        CsrGraph graph = BuildGraph(n, backbone, similarity);
        SimplicialFiltration filtration = RipsFiltration.GraphRips(
            graph, FiltrationWeights.RawDistance, maxDimension, label: "conditioned");

        var edges = new HashSet<(int Lo, int Hi)>(backbone.Count);
        var vertices = new HashSet<int>();
        foreach (var (i, j) in backbone)
        {
            if (i == j) continue;
            edges.Add(i < j ? (i, j) : (j, i));
            vertices.Add(i);
            vertices.Add(j);
        }

        return new RelativeFiltration(filtration, s => InFlagComplex(s, vertices, edges));
    }

    static bool InFlagComplex(Simplex s, HashSet<int> vertices, HashSet<(int Lo, int Hi)> edges)
    {
        int[] v = s.Vertices;
        if (v.Length == 1) return vertices.Contains(v[0]);

        for (int a = 0; a < v.Length; a++)
            for (int b = a + 1; b < v.Length; b++)
                if (!edges.Contains(v[a] < v[b] ? (v[a], v[b]) : (v[b], v[a])))
                    return false;

        return true;
    }

    /// <summary>
    /// <see cref="BuildRelative"/> → <see cref="RelativeFiltration.ComputeBarcode"/>, with representatives,
    /// reporting the degrees below the skeleton's top, which <see cref="RipsFiltration.GraphRips"/> caps at
    /// triangles. Those are exact; the top degree of a truncated skeleton is not, since nothing above it can fill
    /// it. For the top degree, read <see cref="BuildRelative"/> directly.
    /// </summary>
    public static Barcode ComputeRelativeBarcode(
        int n,
        IReadOnlyList<(int i, int j)> backbone,
        IReadOnlyList<(int i, int j, double d)> similarity,
        int maxDimension = 2)
    {
        if (maxDimension < 1)
            throw new ArgumentOutOfRangeException(nameof(maxDimension), "The skeleton needs at least edges.");

        return BuildRelative(n, backbone, similarity, maxDimension)
            .ComputeBarcode(Math.Min(maxDimension, 2) - 1, representatives: true);
    }

    /// <summary>
    /// P1c · a non-monotone sweep of the conditioned graph, as frames for the zigzag engines: frame <c>t</c> holds
    /// the backbone and the content edges of <c>content[t]</c> within slack <c>slack[t]</c>. Edges enter and leave
    /// between frames, as a monotone δ cannot show. With a connected spanning backbone H₀ is one bar, so a sweep's
    /// returns read in H₁ (<c>maxDimension ≥ 1</c>), as unfilled graph cycles.
    /// </summary>
    /// <param name="n">Vertex count; every vertex is present in every frame.</param>
    /// <param name="backbone">Prior <c>K₀</c> edges, present in every frame.</param>
    /// <param name="content">Per frame, the content edges <c>(i, j, d)</c>.</param>
    /// <param name="slack">Per frame, the admission bound: an edge enters at <c>d ≤ slack[t]</c>.</param>
    /// <param name="transition">How consecutive frames connect.</param>
    public static GraphFrameFiltration BuildSweep(
        int n,
        IReadOnlyList<(int i, int j)> backbone,
        IReadOnlyList<IReadOnlyList<(int i, int j, double d)>> content,
        IReadOnlyList<double> slack,
        FrameTransition transition = FrameTransition.Union)
    {
        ArgumentNullException.ThrowIfNull(backbone);
        ArgumentNullException.ThrowIfNull(content);
        ArgumentNullException.ThrowIfNull(slack);
        if (content.Count != slack.Count)
            throw new ArgumentException($"{content.Count} content frames against {slack.Count} slack values.", nameof(slack));

        var edgeSets = new IReadOnlyList<(int i, int j)>[content.Count];
        for (int t = 0; t < edgeSets.Length; t++)
        {
            var edges = new List<(int i, int j)>(backbone);
            foreach (var (i, j, d) in content[t])
                if (d <= slack[t])
                    edges.Add((i, j));
            edgeSets[t] = edges;
        }

        return GraphFrameZigzag.FromEdgeSets(n, edgeSets, transition);
    }
}
