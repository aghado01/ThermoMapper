#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Čufar–Virk involuted persistent homology: cohomology death simplices → assisted boundary
/// reduction for fast barcodes with cycle representatives (Ripser <c>:involuted</c>).
/// </summary>
public static class PersistentInvolutedHomology
{
    public static Barcode Compute(
        SimplicialFiltration filtration,
        int maxDimension = int.MaxValue,
        double cutoff = 0.0,
        bool representatives = false)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        int n = filtration.Simplices.Count;
        var bars = new List<Bar>();
        (List<Bar> h0Bars, HashSet<int> toReduce, HashSet<int> toSkip) =
            PersistenceClearing.ComputeH0(filtration, maxDimension, cutoff);
        bars.AddRange(h0Bars);

        if (maxDimension <= 0)
            return new Barcode(bars, filtration.Label);

        int maxSimplexDim = 0;
        for (int i = 0; i < n; i++)
            if (filtration.Simplices[i].Dimension > maxSimplexDim)
                maxSimplexDim = filtration.Simplices[i].Dimension;

        int dimCap = Math.Min(maxDimension, maxSimplexDim);
        var cohomology = new FiltrationCohomologyReducer(filtration, toReduce, toSkip);

        for (int dim = 1; dim <= dimCap; dim++)
        {
            var deaths = new List<int>();
            var infiniteBirths = new List<int>();
            cohomology.RunDimension(dim, maxDimension, cutoff, representatives, null, deaths, infiniteBirths);

            bars.AddRange(ReduceBoundaryColumns(
                filtration, deaths, maxDimension, cutoff, representatives));

            foreach (int birthIdx in infiniteBirths)
            {
                int dimBirth = filtration.Simplices[birthIdx].Dimension;
                if (dimBirth > maxDimension)
                    continue;

                int[]? cycle = null;
                if (representatives)
                {
                    // H1 only — infinite H2+ voids keep birth-simplex-only reps (surface reconstruction is separate).
                    cycle = dimBirth == 1
                        ? CycleReconstruction.ReconstructH1Cycle(filtration, birthIdx)
                        : new[] { birthIdx };
                }

                bars.Add(new Bar(
                    filtration.Simplices[birthIdx].FiltrationValue,
                    double.PositiveInfinity,
                    dimBirth,
                    birthIdx,
                    Cycle: cycle));
            }

            if (!cohomology.AdvanceDimension(dim, dimCap))
                break;
        }

        return new Barcode(bars, filtration.Label);
    }

    static List<Bar> ReduceBoundaryColumns(
        SimplicialFiltration filtration,
        List<int> columns,
        int maxDimension,
        double cutoff,
        bool representatives)
    {
        var bars = new List<Bar>();
        if (columns.Count == 0)
            return bars;

        int n = filtration.Simplices.Count;
        var pivotCol = new Dictionary<int, int>();
        var reducedCols = new SortedSet<int>[n];
        var sorted = columns.OrderBy(j => j).ToList();

        foreach (int j in sorted)
        {
            var col = new SortedSet<int>(filtration.GetBoundaryIndices(j));
            while (col.Count > 0 && pivotCol.TryGetValue(col.Max, out int i))
                foreach (int r in reducedCols[i])
                {
                    if (!col.Remove(r)) col.Add(r);
                }

            reducedCols[j] = col;

            if (col.Count == 0)
                continue;

            int pivot = col.Max;
            pivotCol[pivot] = j;
            int dim = filtration.Simplices[pivot].Dimension;
            double birth = filtration.Simplices[pivot].FiltrationValue;
            double death = filtration.Simplices[j].FiltrationValue;
            if (dim <= maxDimension && PersistenceClearing.PassesCutoff(birth, death, cutoff))
                bars.Add(new Bar(
                    birth,
                    death,
                    dim,
                    pivot,
                    Cycle: FiniteCycleRep(representatives, col, pivot)));
        }

        return bars;
    }

    /// <summary>Homology cycle rep: reduced column at death plus birth pivot (Ripser explicit).</summary>
    static int[]? FiniteCycleRep(bool representatives, SortedSet<int> col, int pivot)
    {
        if (!representatives)
            return null;

        var cycle = new SortedSet<int>(col) { pivot };
        var indices = new int[cycle.Count];
        cycle.CopyTo(indices);
        return indices;
    }
}
