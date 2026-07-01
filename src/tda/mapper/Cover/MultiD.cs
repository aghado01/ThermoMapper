// ============================================================================
// TDA.Mapper — Multi.cs
// ============================================================================
// Multi-D filter and cover implementations for product-cover MAPPER.
//
// CompositeFilter : IMultiFilter
//   Composes N existing 1-D filters into a single multi-D lens. The Dimension-th
//   coordinate of each output filter vector is the value of the Dimension-th 1-D
//   filter at that point. Example: composing PCA1 + PoincareRadial gives a 2-D
//   lens with PCA1 on axis 0 and Poincaré radial on axis 1.
//
// ProductCover : IMultiCover
//   Constructs a multi-D cover as the Cartesian product of per-dimension 1-D
//   covers (one ICover per filter dimension). Each multi-D bin is the
//   intersection of one 1-D bin per dimension. Bin counts can grow rapidly
//   (Π_d |bins_d|), so be deliberate about Dimension and per-dim NumIntervals.
//
// Reference: standard product-cover construction in the MAPPER literature
// (Singh, Mémoli, Carlsson 2007; Carrière et al. 2017).
// ============================================================================

#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Mapper;

namespace TDA.Mapper.Cover;

// ── CompositeFilter ─────────────────────────────────────────────────────────

/// <summary>
/// Composes N existing 1-D filters into a multi-D filter. Each output filter
/// vector has length equal to the number of component filters.
/// </summary>
public sealed class CompositeFilter : IMultiFilter
{
    public IReadOnlyList<IFilter> ComponentFilters { get; }
    public int Dimension { get; }

    public string Name => $"Composite[{string.Join(" × ", ComponentFilters.Select(f => f.Name))}]";

    public CompositeFilter(params IFilter[] filters)
    {
        ArgumentNullException.ThrowIfNull(filters);
        if (filters.Length < 2)
            throw new ArgumentException("CompositeFilter requires at least 2 component filters.", nameof(filters));
        foreach (var f in filters)
            ArgumentNullException.ThrowIfNull(f);

        ComponentFilters = filters;
        Dimension = filters.Length;
    }

    public double[][] Apply(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);

        // Apply each filter once over the full dataset → one double[N] per dim.
        var perDim = new double[Dimension][];
        for (int d = 0; d < Dimension; d++)
            perDim[d] = ComponentFilters[d].Apply(data);

        // Pack into row-jagged: result[i][d] = filter d's value at point i.
        var result = new double[data.Length][];
        for (int i = 0; i < data.Length; i++)
        {
            result[i] = new double[Dimension];
            for (int d = 0; d < Dimension; d++)
                result[i][d] = perDim[d][i];
        }
        return result;
    }
}

// ── ProductCover ────────────────────────────────────────────────────────────

/// <summary>
/// Cartesian product of per-dimension 1-D covers. Each multi-D bin is the
/// intersection of one 1-D bin per dimension.
///
/// Bin count grows as Π_d |bins_d|; with 12 intervals per dimension and 2-D,
/// expect up to 144 product bins (often fewer in practice due to empty
/// intersections). 3-D with 8 intervals yields up to 512. Be deliberate.
/// </summary>
public sealed class ProductCover : IMultiCover
{
    public IReadOnlyList<ICover> PerDimensionCovers { get; }
    public int Dimension { get; }

    public string Name => $"Product[{string.Join(" × ", PerDimensionCovers.Select(c => c.Name))}]";

    public ProductCover(params ICover[] perDimensionCovers)
    {
        ArgumentNullException.ThrowIfNull(perDimensionCovers);
        if (perDimensionCovers.Length < 2)
            throw new ArgumentException("ProductCover requires at least 2 per-dimension covers.", nameof(perDimensionCovers));
        foreach (var c in perDimensionCovers)
            ArgumentNullException.ThrowIfNull(c);

        PerDimensionCovers = perDimensionCovers;
        Dimension = perDimensionCovers.Length;
    }

    public CoverResult Generate(double[][] filterValues)
    {
        ArgumentNullException.ThrowIfNull(filterValues);
        if (filterValues.Length == 0)
            return new CoverResult(Array.Empty<CoverBin>(), 0.0, 0.0);

        int n = filterValues.Length;

        // 1. Slice the multi-D filter values into per-dimension scalar columns.
        var perDimValues = new double[Dimension][];
        for (int d = 0; d < Dimension; d++)
        {
            perDimValues[d] = new double[n];
            for (int p = 0; p < n; p++)
                perDimValues[d][p] = filterValues[p].Length > d ? filterValues[p][d] : 0.0;
        }

        // 2. Run each 1-D cover on its column → per-dim bin lists.
        var perDimBins = new IReadOnlyList<CoverBin>[Dimension];
        double globalMin = double.PositiveInfinity, globalMax = double.NegativeInfinity;
        for (int d = 0; d < Dimension; d++)
        {
            var result = PerDimensionCovers[d].Generate(perDimValues[d]);
            perDimBins[d] = result.Bins;
            if (result.FilterMin < globalMin) globalMin = result.FilterMin;
            if (result.FilterMax > globalMax) globalMax = result.FilterMax;
        }

        // 3. Build per-dim membership sets for fast intersection.
        //    membershipMask[d][p] = true if point p is in the currently-selected bin of dim d.
        //    We rebuild masks per (d, current bin) on demand instead of pre-allocating
        //    Π_d |bins_d| × N booleans.
        //
        //    Strategy: iterate odometer-style over the bin tuple (i_0, i_1, ..., i_{D-1}).
        //    For each tuple, compute the intersection of the corresponding 1-D bins'
        //    point sets using the smallest bin as the seed and HashSet lookups for the rest.

        var productBins = new List<CoverBin>();
        var currentIndices = new int[Dimension];
        var dimSets = new HashSet<int>[Dimension];
        int productBinId = 0;

        // Pre-size dimSets to avoid per-iteration allocation.
        // We refresh each set at every tuple step.
        for (int d = 0; d < Dimension; d++) dimSets[d] = new HashSet<int>();

        bool done = false;
        while (!done)
        {
            // Identify the smallest bin among the current tuple — seed for intersection.
            int smallestDim = 0;
            int smallestCount = perDimBins[0][currentIndices[0]].PointIndices.Count;
            for (int d = 1; d < Dimension; d++)
            {
                int c = perDimBins[d][currentIndices[d]].PointIndices.Count;
                if (c < smallestCount)
                {
                    smallestDim = d;
                    smallestCount = c;
                }
            }

            // Build hash sets for non-seed dims (only needed if Dimension > 1).
            for (int d = 0; d < Dimension; d++)
            {
                if (d == smallestDim) continue;
                dimSets[d].Clear();
                foreach (int p in perDimBins[d][currentIndices[d]].PointIndices)
                    dimSets[d].Add(p);
            }

            var seedBin = perDimBins[smallestDim][currentIndices[smallestDim]];
            var intersected = new List<int>(seedBin.PointIndices.Count);

            foreach (int p in seedBin.PointIndices)
            {
                bool inAll = true;
                for (int d = 0; d < Dimension; d++)
                {
                    if (d == smallestDim) continue;
                    if (!dimSets[d].Contains(p)) { inAll = false; break; }
                }
                if (inAll) intersected.Add(p);
            }

            if (intersected.Count > 0)
            {
                // Lower/Upper use first-dim bounds as a representative summary; the actual
                // multi-D extent of the bin is implicit in the intersection geometry.
                // Downstream viz can recover per-dim bounds from member indices if needed.
                productBins.Add(new CoverBin(
                    BinId: productBinId,
                    Lower: perDimBins[0][currentIndices[0]].Lower,
                    Upper: perDimBins[0][currentIndices[0]].Upper,
                    PointIndices: intersected));
                productBinId++;
            }

            // Advance the odometer.
            int dimAdvance = 0;
            while (dimAdvance < Dimension)
            {
                currentIndices[dimAdvance]++;
                if (currentIndices[dimAdvance] < perDimBins[dimAdvance].Count) break;
                currentIndices[dimAdvance] = 0;
                dimAdvance++;
            }
            if (dimAdvance == Dimension) done = true;
        }

        return new CoverResult(productBins, globalMin, globalMax);
    }
}
