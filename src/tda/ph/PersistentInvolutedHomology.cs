#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Čufar–Virk involuted persistent homology: cohomology death simplices → assisted boundary
/// reduction for fast barcodes with cycle representatives (Ripser <c>:involuted</c>).
/// <para>With representatives, every bar has a <see cref="Bar.Cycle"/>. An essential H₁ bar's is a shortest path closed by its
/// birth edge. An essential bar higher up has the chain <c>V_i</c> of its birth cell, the same as
/// <see cref="PersistentHomology"/>'s, and it costs a reduction of that cell's boundary through the degree's death columns.
/// <c>maxCycleDimension</c> is the highest degree that gets one: the bars above it stay, without a cycle, so a path that reads
/// none of them pays for none.</para>
/// </summary>
public static class PersistentInvolutedHomology
{
    public static Barcode Compute(
        SimplicialFiltration filtration,
        int maxDimension = int.MaxValue,
        double cutoff = 0.0,
        bool representatives = false,
        int maxCycleDimension = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        int n = filtration.Simplices.Count;
        var bars = new List<Bar>();
        (List<Bar> h0Bars, HashSet<int> toReduce) = PersistenceClearing.ComputeH0(filtration, maxDimension, cutoff, cycles: representatives && maxCycleDimension >= 0);
        bars.AddRange(h0Bars);

        if (maxDimension <= 0)
            return new Barcode(bars, filtration.Label);

        int maxSimplexDim = 0;
        for (int i = 0; i < n; i++)
            if (filtration.Simplices[i].Dimension > maxSimplexDim)
                maxSimplexDim = filtration.Simplices[i].Dimension;

        int dimCap = Math.Min(maxDimension, maxSimplexDim);
        var cohomology = new FiltrationCohomologyReducer(filtration, toReduce);
        ReducedDeaths? deathsOfThisDegree = null;

        for (int dim = 1; dim <= dimCap; dim++)
        {
            var deaths = new List<int>();
            var infiniteBirths = new List<int>();
            cohomology.RunDimension(dim, maxDimension, cutoff, representatives, null, deaths, infiniteBirths);

            var deathsOfNextDegree = new ReducedDeaths(n, keepChains: representatives && dim + 1 <= maxCycleDimension);
            bars.AddRange(ReduceBoundaryColumns(
                filtration, deaths, maxDimension, cutoff, representatives && dim <= maxCycleDimension, deathsOfNextDegree));

            foreach (int birthIdx in infiniteBirths)
            {
                int dimBirth = filtration.Simplices[birthIdx].Dimension;
                if (dimBirth > maxDimension)
                    continue;

                int[]? cycle = null;
                if (representatives && dimBirth <= maxCycleDimension)
                {
                    cycle = dimBirth == 1
                        ? CycleReconstruction.ReconstructH1Cycle(filtration, birthIdx)
                        : deathsOfThisDegree!.EssentialCycle(filtration, birthIdx);
                }

                bars.Add(new Bar(
                    filtration.Simplices[birthIdx].FiltrationValue,
                    double.PositiveInfinity,
                    dimBirth,
                    birthIdx,
                    Cycle: cycle));
            }

            deathsOfThisDegree = deathsOfNextDegree;

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
        bool representatives,
        ReducedDeaths reduced)
    {
        var bars = new List<Bar>();
        if (columns.Count == 0)
            return bars;

        Dictionary<int, int> pivotCol = reduced.PivotCol;
        SortedSet<int>[] reducedCols = reduced.Columns;
        List<int>? took = reduced.Chains is null ? null : new List<int>();
        var sorted = columns.OrderBy(j => j).ToList();

        foreach (int j in sorted)
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

            reduced.Chains?.Record(j, took!);
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

    /// <summary>
    /// One degree's death columns, reduced. They are all that degree's non-zero columns in a full reduction, so a cell whose
    /// column reduces to zero takes the same columns through them as it would there.
    /// </summary>
    sealed class ReducedDeaths
    {
        public readonly Dictionary<int, int> PivotCol = new();
        public readonly SortedSet<int>[] Columns;
        public readonly ReductionChains? Chains;

        public ReducedDeaths(int n, bool keepChains)
        {
            Columns = new SortedSet<int>[n];
            Chains = keepChains ? new ReductionChains(n) : null;
        }

        /// <summary>The cycle of a cell of this degree that kills nothing: <c>V</c> of its column, as in a full reduction.</summary>
        public int[] EssentialCycle(SimplicialFiltration filtration, int birth)
        {
            var col = new SortedSet<int>(filtration.GetBoundaryIndices(birth));
            var took = new List<int>();

            while (col.Count > 0)
            {
                int i = PivotCol[col.Max];
                foreach (int r in Columns[i])
                {
                    if (!col.Remove(r)) col.Add(r);
                }

                took.Add(i);
            }

            Chains!.Record(birth, took);
            return Chains.ChainOf(birth);
        }
    }
}
