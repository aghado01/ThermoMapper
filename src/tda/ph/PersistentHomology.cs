#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Standard persistent homology (Z/2Z column reduction) — explicit correctness path.
/// </summary>
public static class PersistentHomology
{
    public static Barcode Compute(SimplicialFiltration filtration, int maxDimension = int.MaxValue, bool representatives = false) =>
        Compute((IFiltration)filtration, maxDimension, representatives);

    public static Barcode Compute(IFiltration filtration, int maxDimension = int.MaxValue, bool representatives = false)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        int n = filtration.Count;
        var bars = new List<Bar>();
        var paired = new bool[n];
        var pivotCol = new Dictionary<int, int>(n);
        var reducedCols = new SortedSet<int>[n];

        for (int j = 0; j < n; j++)
        {
            var col = new SortedSet<int>(filtration.GetBoundaryIndices(j));

            while (col.Count > 0 && pivotCol.TryGetValue(col.Max, out int i))
                foreach (int r in reducedCols[i])
                {
                    if (!col.Remove(r)) col.Add(r);
                }

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
                        Cycle: representatives ? ColumnToArray(col) : null));
            }
        }

        for (int i = 0; i < n; i++)
            if (!paired[i])
            {
                int dim = filtration.GetDimension(i);
                if (dim <= maxDimension)
                    bars.Add(new Bar(
                        filtration.GetBirth(i),
                        double.PositiveInfinity, dim, i));
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
