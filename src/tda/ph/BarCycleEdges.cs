#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Normalized undirected vertex pair (Lo ≤ Hi) from an H1 cycle representative edge.
/// </summary>
public readonly record struct UndirectedEdge(int Lo, int Hi);

/// <summary>
/// Maps involuted <see cref="Bar.Cycle"/> simplex indices to LMP-consumable edge sets.
/// </summary>
public static class BarCycleEdges
{
    /// <summary>
    /// Extract 1-simplex edges from <paramref name="bar"/>.Cycle as normalized vertex pairs.
    /// Higher-dimensional simplices in the chain are skipped (H1 cycles should be edge-only).
    /// </summary>
    public static IReadOnlyList<UndirectedEdge> GetEdgePairs(in Bar bar, SimplicialFiltration filtration)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        if (bar.Cycle is null || bar.Cycle.Length == 0)
            return Array.Empty<UndirectedEdge>();

        var edges = new List<UndirectedEdge>(bar.Cycle.Length);
        foreach (int idx in bar.Cycle)
        {
            Simplex s = filtration.Simplices[idx];
            if (s.Dimension != 1 || s.Vertices.Length != 2)
                continue;

            edges.Add(new UndirectedEdge(s.Vertices[0], s.Vertices[1]));
        }

        return edges;
    }

    /// <summary>All H1 bars in <paramref name="barcode"/> with non-empty edge cycles.</summary>
    public static IEnumerable<(Bar Bar, IReadOnlyList<UndirectedEdge> Edges)> H1Loops(
        Barcode barcode,
        SimplicialFiltration filtration)
    {
        ArgumentNullException.ThrowIfNull(barcode);
        ArgumentNullException.ThrowIfNull(filtration);

        foreach (Bar bar in barcode.Bars)
        {
            if (bar.Dimension != 1)
                continue;

            IReadOnlyList<UndirectedEdge> edges = GetEdgePairs(bar, filtration);
            if (edges.Count > 0)
                yield return (bar, edges);
        }
    }
}
