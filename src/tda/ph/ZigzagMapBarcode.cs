#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// General simplicial-zigzag oracle: the zigzag barcode of a sequence of complexes connected by
/// ARBITRARY cell maps (not only single-cell inclusions). Needed to validate strong-collapse
/// core-assembly (§4 of 1809.10945), whose induced maps are retractions. Reuses the homology /
/// generalized-rank / interval-decomposition machinery of <see cref="ZigzagBarcodeNaive"/> — only the
/// induced map differs (push a cycle through the cell map, then express it in the target basis).
/// </summary>
public static class ZigzagMapBarcode
{
    /// <summary>One zigzag step: its arrow direction and the induced cell map on the source complex.</summary>
    public readonly struct Step
    {
        /// <summary>true: complexes[i] → complexes[i+1]; false: complexes[i+1] → complexes[i].</summary>
        public readonly bool IsForward;
        /// <summary>CellMap[sourceCell] = the target cell it maps to, or -1 if it maps to 0 (degenerate).</summary>
        public readonly int[] CellMap;
        public Step(bool isForward, int[] cellMap) { IsForward = isForward; CellMap = cellMap; }
    }

    /// <summary>
    /// Barcode of a complex sequence. <paramref name="complexes"/> has <c>steps.Count + 1</c> entries
    /// (cell-id sets over a shared universe of <paramref name="N"/> cells with <paramref name="dimOf"/>
    /// / <paramref name="bndOf"/>); <paramref name="steps"/>[i] connects complexes[i] and complexes[i+1].
    /// </summary>
    public static Barcode Compute(
        IReadOnlyList<HashSet<int>> complexes,
        int[] dimOf,
        int[][] bndOf,
        int N,
        IReadOnlyList<Step> steps,
        int maxDimension = int.MaxValue)
    {
        int numSteps = steps.Count;
        if (numSteps == 0 || complexes.Count == 0) return new Barcode(Array.Empty<Bar>(), "ZigzagMap");

        var isForward = new bool[numSteps];
        for (int i = 0; i < numSteps; i++) isForward[i] = steps[i].IsForward;

        var bars = new List<Bar>();
        for (int p = 0; p <= maxDimension; p++)
        {
            var V = new List<List<bool[]>>();
            for (int i = 0; i <= numSteps; i++)
                V.Add(ZigzagBarcodeNaive.ComputeHomologyBasis(complexes[i], p, dimOf, bndOf, N));

            var M = new List<bool[,]>();
            for (int i = 0; i < numSteps; i++)
            {
                bool fwd = steps[i].IsForward;
                var srcV = fwd ? V[i] : V[i + 1];
                var tgtV = fwd ? V[i + 1] : V[i];
                var tgtK = fwd ? complexes[i + 1] : complexes[i];
                var mapped = MapCycles(srcV, steps[i].CellMap, N);
                M.Add(ZigzagBarcodeNaive.ComputeInducedMapMatrix(mapped, tgtV, tgtK, p + 1, dimOf, bndOf, N));
            }

            bars.AddRange(ZigzagBarcodeNaive.DecomposeDimension(numSteps, V, M, isForward, p));
        }

        return new Barcode(bars, "ZigzagMap");
    }

    /// <summary>Push each cycle through the cell map (Z/2): image[CellMap[c]] ^= 1 for each c in the cycle.</summary>
    static List<bool[]> MapCycles(List<bool[]> cycles, int[] cellMap, int N)
    {
        var result = new List<bool[]>(cycles.Count);
        foreach (var z in cycles)
        {
            var img = new bool[N];
            for (int c = 0; c < N; c++)
                if (z[c] && cellMap[c] >= 0) img[cellMap[c]] ^= true;
            result.Add(img);
        }
        return result;
    }

    /// <summary>
    /// Convenience: treat a single-cell <see cref="ZigzagFiltration"/> as complexes with INCLUSION
    /// maps (identity cell maps). Must equal <see cref="ZigzagBarcodeNaive"/> exactly — the isolation
    /// check that this general oracle's machinery is correct.
    /// </summary>
    public static Barcode ComputeFromZigzag(ZigzagFiltration f, int maxDimension = int.MaxValue)
    {
        int numSteps = f.Count;
        if (numSteps == 0) return new Barcode(Array.Empty<Bar>(), "ZigzagMap");

        int maxCellId = -1;
        foreach (var s in f) if (s.GlobalCellId > maxCellId) maxCellId = s.GlobalCellId;
        int N = maxCellId + 1;

        var dimOf = new int[N];
        var bndOf = new int[N][];
        for (int i = 0; i < N; i++) bndOf[i] = Array.Empty<int>();

        var complexes = new List<HashSet<int>> { new HashSet<int>() };
        var current = new HashSet<int>();
        for (int i = 0; i < numSteps; i++)
        {
            var s = f[i];
            if (s.Direction == ZigzagDirection.Add)
            {
                current.Add(s.GlobalCellId);
                var bnd = s.BoundaryAtAdd ?? Array.Empty<int>();
                dimOf[s.GlobalCellId] = bnd.Length > 0 ? dimOf[bnd[0]] + 1 : 0;
                bndOf[s.GlobalCellId] = bnd;
            }
            else current.Remove(s.GlobalCellId);
            complexes.Add(new HashSet<int>(current));
        }

        var identity = new int[N];
        for (int c = 0; c < N; c++) identity[c] = c;

        var steps = new List<Step>(numSteps);
        for (int i = 0; i < numSteps; i++)
            steps.Add(new Step(f[i].Direction == ZigzagDirection.Add, identity));

        return Compute(complexes, dimOf, bndOf, N, steps, maxDimension);
    }
}
