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
/// <see cref="PersistentInvolutedHomology"/>. Directedness, the <c>Δ</c> reach axis, the gauge
/// <c>q</c>, the sheaf / <c>λ_q</c> reading, and zigzag slices are P1–P4 and deliberately absent.</para>
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
}
