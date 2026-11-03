#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Standard persistent homology (Z/2Z column reduction) — explicit correctness path.
/// <para>With representatives, every bar has a <see cref="Bar.Cycle"/>: the reduced column of what kills it, or, on an
/// essential bar, the chain <c>V_i</c> of its birth cell, whose boundary is that column, here zero. <c>maxCycleDimension</c>
/// is the highest degree that gets one: the bars above it stay, without a cycle, so a path that reads none of them pays for none.</para>
/// </summary>
public static class PersistentHomology
{
    public static Barcode Compute(SimplicialFiltration filtration, int maxDimension = int.MaxValue, bool representatives = false, int maxCycleDimension = int.MaxValue) =>
        Compute((IFiltration)filtration, maxDimension, representatives, maxCycleDimension);

    public static Barcode Compute(IFiltration filtration, int maxDimension = int.MaxValue, bool representatives = false, int maxCycleDimension = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        return Reduce(filtration, maxDimension, representatives ? new ReductionChains(filtration.Count) : null, maxCycleDimension);
    }

    /// <summary>
    /// <see cref="Compute(IFiltration, int, bool, int)"/> with representatives, also handing back <c>V</c>, to read any column's chain.
    /// </summary>
    internal static Barcode ComputeWithChains(IFiltration filtration, int maxDimension, out ReductionChains chains)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        chains = new ReductionChains(filtration.Count);
        return Reduce(filtration, maxDimension, chains, int.MaxValue);
    }

    static Barcode Reduce(IFiltration filtration, int maxDimension, ReductionChains? chains, int maxCycleDimension)
    {
        int n = filtration.Count;
        var bars = new List<Bar>();
        var paired = new bool[n];
        var pivotCol = new Dictionary<int, int>(n);
        var reducedCols = new SortedSet<int>[n];
        List<int>? took = chains is null ? null : new List<int>();

        for (int j = 0; j < n; j++)
        {
            var col = new SortedSet<int>(filtration.GetBoundaryIndices(j));
            took?.Clear();

            while (col.Count > 0 && pivotCol.TryGetValue(col.Max, out int i))
            {
                foreach (int r in reducedCols[i])
                {
                    if (!col.Remove(r)) col.Add(r);
                }

                took?.Add(i);
            }

            if (chains is not null && filtration.GetDimension(j) <= maxCycleDimension)
                chains.Record(j, took!);

            reducedCols[j] = col;

            if (col.Count > 0)
            {
                int pivot = col.Max;
                pivotCol[pivot] = j;
                paired[pivot] = true;
                paired[j] = true;

                int dim = filtration.GetDimension(pivot);
                if (dim <= maxDimension)
                    bars.Add(new Bar(
                        filtration.GetBirth(pivot),
                        filtration.GetBirth(j),
                        dim, pivot,
                        // R-matrix column of the negative cell IS the cycle rep (D·R = 0); the
                        // pivot is already its max element. §2.2 of 2412.02591.
                        Cycle: chains is null || dim > maxCycleDimension ? null : ColumnToArray(col)));
            }
        }

        for (int i = 0; i < n; i++)
            if (!paired[i])
            {
                int dim = filtration.GetDimension(i);
                if (dim <= maxDimension)
                    bars.Add(new Bar(
                        filtration.GetBirth(i),
                        double.PositiveInfinity, dim, i,
                        // Unpaired, so R_i = 0 and ∂V_i = 0: V_i is a cycle whose latest cell is i.
                        Cycle: dim > maxCycleDimension ? null : chains?.ChainOf(i)));
            }

        return new Barcode(bars, filtration.Label);
    }

    static int[] ColumnToArray(SortedSet<int> col)
    {
        var a = new int[col.Count];
        col.CopyTo(a);
        return a;
    }
}
