using System;

namespace Clustering.Primitives;

/// <summary>
/// Provenance coordinates a landscape's topology is only defined relative to:
/// which sink built the field, on which graph, under which gauge.
/// </summary>
public sealed record LandscapeProvenance(string Sink, string GraphId, string? GaugeNote = null);

/// <summary>
/// Per-node scalar field sampled over a declared axis — the height the
/// resolution layer walks; the third resolution carrier beside
/// <see cref="Assignment"/> and <see cref="Groups"/>. Grid-major:
/// <see cref="ValuesByGridPoint"/>[g] is the per-node slice at
/// <see cref="Grid"/>[g] — the shape producers mint (one column per sweep
/// temperature) and every consumer reads (selector walks consume columns in
/// axis order; periphery policies consume a single column). A one-point grid
/// is the slice degenerate.
/// </summary>
/// <remarks>
/// Walk strategies enforce the axis-alignment law
/// (<c>Dendrogram.CostAxis == Axis</c>) and consume the field CARDINALLY —
/// the gauge rides in <see cref="Provenance"/> and masses are sink-relative
/// by design. Ascent-style periphery policies consume only the field's
/// order and are gauge-free.
/// </remarks>
public sealed record Landscape
{
    public required string Axis { get; init; }

    /// <summary>Strictly ascending axis samples.</summary>
    public required double[] Grid { get; init; }

    /// <summary>One per-node column per grid point: <c>[gridIndex][node]</c>.</summary>
    public required double[][] ValuesByGridPoint { get; init; }

    public LandscapeProvenance? Provenance { get; init; }

    public int GridCount => Grid.Length;

    public int NodeCount => ValuesByGridPoint.Length > 0 ? ValuesByGridPoint[0].Length : 0;

    /// <summary>
    /// Shape-validating factory: strictly ascending grid, one column per grid
    /// point, equal column lengths.
    /// </summary>
    public static Landscape Create(
        string axis, double[] grid, double[][] valuesByGridPoint, LandscapeProvenance? provenance = null)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(axis);
        ArgumentNullException.ThrowIfNull(grid);
        ArgumentNullException.ThrowIfNull(valuesByGridPoint);
        if (grid.Length == 0)
            throw new ArgumentException("Grid must contain at least one point.", nameof(grid));
        if (valuesByGridPoint.Length != grid.Length)
            throw new ArgumentException(
                $"One column per grid point: {valuesByGridPoint.Length} columns vs {grid.Length} grid points.",
                nameof(valuesByGridPoint));
        for (int g = 1; g < grid.Length; g++)
            if (grid[g] <= grid[g - 1])
                throw new ArgumentException("Grid must be strictly ascending.", nameof(grid));
        int nodes = valuesByGridPoint[0].Length;
        for (int g = 1; g < valuesByGridPoint.Length; g++)
            if (valuesByGridPoint[g].Length != nodes)
                throw new ArgumentException("All columns must have the same node count.", nameof(valuesByGridPoint));

        return new Landscape
        {
            Axis = axis,
            Grid = grid,
            ValuesByGridPoint = valuesByGridPoint,
            Provenance = provenance,
        };
    }
}
