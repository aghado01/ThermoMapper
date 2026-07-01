#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace TDA.Ph;

/// <summary>
/// Position-indexed Z2 RU decomposition with vineyard transposition — Cohen-Steiner–Edelsbrunner–Morozov
/// "Vines and Vineyards" (SoCG 2006), ported from GUDHI's <c>ru_vine_swap.h</c> (Hannah Schreiber, MIT).
/// The substrate for <see cref="ReflectionZigzag"/>'s backward arrows (the Maria–Oudot zigzag engine).
///
/// <para>Columns are sorted-set Z2 vectors over positions; <c>R</c> is the reduced boundary and <c>U</c> the
/// GUDHI "mirror" (= V⁻¹). Oracle-grade (naive ops, recompute pairing from R rather than maintain it) — this
/// is the independent <i>check</i>, not the fast engine. Kept <c>internal</c>: one consumer today, do not
/// promote to a shared primitive until a second appears.</para>
///
/// <para><b>Invariant validated in isolation</b> (RuVineyardTests): a <see cref="VineSwap"/> of two adjacent,
/// face/coface-free cells leaves R reduced to exactly the pairing a from-scratch reduction of the reordered
/// filtration would produce. The transpose/barcode bookkeeping GUDHI maintains incrementally is skipped here —
/// the pairing is recomputed from R's pivots, which is equivalent and simpler.</para>
/// </summary>
internal sealed class RuVineyard
{
    readonly List<SortedSet<int>> _r;   // reduced boundary columns (row positions)
    readonly List<SortedSet<int>> _u;   // the V⁻¹ "mirror"
    readonly int[] _dim;

    public int Count => _r.Count;

    /// <param name="cells">Filtration order: cell at position p has the given dimension and a boundary listed
    /// as the positions (&lt; p) of its faces.</param>
    public RuVineyard(IReadOnlyList<(int Dim, int[] Boundary)> cells)
    {
        int n = cells.Count;
        _r = new List<SortedSet<int>>(n);
        _u = new List<SortedSet<int>>(n);
        _dim = new int[n];
        for (int p = 0; p < n; p++)
        {
            _r.Add(new SortedSet<int>(cells[p].Boundary));
            _u.Add(new SortedSet<int> { p });
            _dim[p] = cells[p].Dim;
        }
        // Standard reduction, using the same AddTo convention the swaps use (so U is maintained consistently).
        var pivotToCol = new Dictionary<int, int>();
        for (int j = 0; j < n; j++)
        {
            while (_r[j].Count > 0)
            {
                int low = _r[j].Max;
                if (pivotToCol.TryGetValue(low, out int j2)) AddTo(j2, j);
                else { pivotToCol[low] = j; break; }
            }
        }
    }

    /// <summary>Finite pairs (birthPos, deathPos, dim) and infinite ones (birthPos, -1, dim), read from R.</summary>
    public List<(int Birth, int Death, int Dim)> Pairs()
    {
        var res = new List<(int, int, int)>();
        var paired = new HashSet<int>();
        for (int q = 0; q < _r.Count; q++)
            if (_r[q].Count > 0) { int p = _r[q].Max; res.Add((p, q, _dim[p])); paired.Add(p); }
        for (int p = 0; p < _r.Count; p++)
            if (_r[p].Count == 0 && !paired.Contains(p)) res.Add((p, -1, _dim[p]));
        return res;
    }

    // ── primitives (ru_vine_swap.h) ───────────────────────────────────────────────────────────────────
    bool IsZeroColumn(int p) => _r[p].Count == 0;
    bool IsZeroEntryR(int col, int row) => !_r[col].Contains(row);
    bool IsZeroEntryU(int col, int row) => !_u[col].Contains(row);
    void ZeroEntryU(int col, int row) => _u[col].Remove(row);

    static void Xor(SortedSet<int> target, SortedSet<int> source)
    {
        foreach (int x in source) if (!target.Add(x)) target.Remove(x);
    }

    // add_to(source, target): R[target] ^= R[source];  U[source] ^= U[target]  (reversed on U — the V⁻¹ detail)
    void AddTo(int source, int target) { Xor(_r[target], _r[source]); Xor(_u[source], _u[target]); }

    void SwapAtIndex(int i)
    {
        (_r[i], _r[i + 1]) = (_r[i + 1], _r[i]);
        (_u[i], _u[i + 1]) = (_u[i + 1], _u[i]);
        (_dim[i], _dim[i + 1]) = (_dim[i + 1], _dim[i]);   // dimension follows the cell to its new position
        SwapRow(_r, i, i + 1);
        SwapRow(_u, i, i + 1);
    }

    static void SwapRow(List<SortedSet<int>> m, int a, int b)
    {
        foreach (var col in m)
        {
            bool ha = col.Contains(a), hb = col.Contains(b);
            if (ha != hb) { if (ha) { col.Remove(a); col.Add(b); } else { col.Remove(b); col.Add(a); } }
        }
    }

    // death(positive p) = position of the negative column whose pivot is p, or -1.
    int Death(int p)
    {
        for (int q = 0; q < _r.Count; q++) if (_r[q].Count > 0 && _r[q].Max == p) return q;
        return -1;
    }
    int Birth(int q) => _r[q].Max;

    /// <summary>Transpose the cells at positions <paramref name="i"/> and i+1 (must have no face/coface
    /// relation). Returns true iff the barcode changed. Full <c>vine_swap</c> (handles trivial swaps).</summary>
    public bool VineSwap(int i)
    {
        bool iPos = IsZeroColumn(i), iiPos = IsZeroColumn(i + 1);
        bool diffDim = _dim[i] != _dim[i + 1];

        if (iPos && iiPos)
        {
            if (diffDim) { SwapAtIndex(i); return true; }
            if (!IsZeroEntryU(i, i + 1)) ZeroEntryU(i, i + 1);
            return PositiveVineSwap(i);
        }
        if (!iPos && !iiPos)
        {
            if (diffDim || IsZeroEntryU(i, i + 1)) { SwapAtIndex(i); return true; }
            return NegativeVineSwap(i);
        }
        if (iPos && !iiPos)
        {
            if (diffDim || IsZeroEntryU(i, i + 1)) { SwapAtIndex(i); return true; }
            return PositiveNegativeVineSwap(i);
        }
        if (diffDim || IsZeroEntryU(i, i + 1)) { SwapAtIndex(i); return true; }
        return NegativePositiveVineSwap(i);
    }

    bool PositiveVineSwap(int i)
    {
        int iD = Death(i), iiD = Death(i + 1);
        if (iD != -1 && iiD != -1 && !IsZeroEntryR(iiD, i))
        {
            if (iD < iiD) { SwapAtIndex(i); AddTo(iD, iiD); return true; }
            SwapAtIndex(i); AddTo(iiD, iD); return false;
        }
        SwapAtIndex(i);
        if (iD != -1 || iiD == -1 || IsZeroEntryR(iiD, i + 1)) return true;
        return false;
    }

    bool NegativeVineSwap(int i)
    {
        int iB = Birth(i), iiB = Birth(i + 1);
        AddTo(i, i + 1);
        SwapAtIndex(i);
        if (iB < iiB) return true;
        AddTo(i, i + 1);
        return false;
    }

    bool PositiveNegativeVineSwap(int i)
    {
        ZeroEntryU(i, i + 1);
        SwapAtIndex(i);
        return true;
    }

    bool NegativePositiveVineSwap(int i)
    {
        AddTo(i, i + 1);
        SwapAtIndex(i);
        AddTo(i, i + 1);
        return false;
    }
}
