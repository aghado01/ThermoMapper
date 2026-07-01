#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// ReflectionZigzag — a third, genuinely independent general zigzag oracle: the Maria–Oudot algorithm,
/// "zigzag persistence via reflections and transpositions" (the algorithm behind GUDHI's
/// <c>zigzag_persistence.h</c>, Clément Maria / Hannah Schreiber), on a compact Z2 matrix. Distinct in
/// mechanism from the rank/inclusion–exclusion <see cref="ZigzagBarcodeNaive"/> (Z1) and the cone
/// reduce-to-standard <see cref="FastZigzag"/> (Z2) — so it guards the general engine against a *shared*
/// blind spot (the kind that let the Z2 cell-re-entry bug hide behind a 2-way check).
///
/// <para>The engine is a streaming Z2 boundary reduction. A <b>forward arrow</b> reduces the inserted cell's
/// boundary: empty ⇒ a new class opens (<i>birth</i>); non-empty ⇒ it closes the class born at the pivot
/// (<i>death</i>). A <b>backward arrow</b> removes a maximal cell: if its (reduced) column is empty it is an
/// unpaired positive cell whose class <i>dies</i>; if non-empty it is a negative cell and its pivot partner is
/// <i>reborn</i>. Bar ends follow the arrow kind: BirthEnd Closed iff born at an Add; DeathEnd Closed iff died
/// at a Delete.</para>
///
/// <para><b>Stages.</b> Forward arrows (stage 1) and <b>reverse-teardown</b> backward arrows (stage 2b-i, the
/// removed cell already at the last position) are implemented and validated against Z1/Z2 on increasing and
/// up-down filtrations. Arbitrary removal — transpose the maximal cell to the last position via the vineyard
/// (<see cref="RuVineyard"/>) — is stage 2b-ii; a <see cref="ZigzagDirection.Delete"/> of a non-last cell throws
/// until then. Full ∅→∅ churn as a 5th cross-check corner is stage 3.</para>
/// </summary>
public static class ReflectionZigzag
{
    public static Barcode Compute(ZigzagFiltration f, int maxDimension = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(f);
        int m = f.Count;
        var bars = new List<Bar>();
        if (m == 0) return new Barcode(bars, "ReflectionZigzag");

        var r = new List<SortedSet<int>>();          // reduced boundary column by position
        var dim = new List<int>();                    // dimension by position
        var posOf = new Dictionary<int, int>();       // cell id -> current position
        var cellAt = new List<int>();                 // position -> cell id
        var pivotToCol = new Dictionary<int, int>();  // pivot position -> negative column position
        var births = new Dictionary<int, (int Arrow, bool ByAdd)>();  // alive positive position -> birth

        void Emit(int bArrow, bool bByAdd, int dArrow, int d, bool dByDelete)
        {
            if (maxDimension != int.MaxValue && d > maxDimension) return;
            bars.Add(new Bar(bArrow, dArrow, d, null, null, null,
                bByAdd ? IntervalEnd.Closed : IntervalEnd.Open,
                dByDelete ? IntervalEnd.Closed : IntervalEnd.Open));
        }

        for (int arrow = 0; arrow < m; arrow++)
        {
            var s = f[arrow];
            int cell = s.GlobalCellId;

            if (s.Direction == ZigzagDirection.Add)
            {
                var bnd = s.BoundaryAtAdd!;
                int pos = r.Count;
                int d = bnd.Length == 0 ? 0 : dim[posOf[bnd[0]]] + 1;
                var col = new SortedSet<int>(bnd.Select(b => posOf[b]));
                r.Add(col); dim.Add(d); cellAt.Add(cell); posOf[cell] = pos;

                while (col.Count > 0)
                {
                    int low = col.Max;
                    if (pivotToCol.TryGetValue(low, out int j)) Xor(col, r[j]);
                    else break;
                }

                if (col.Count == 0)
                {
                    births[pos] = (arrow, true);                                 // birth (forward)
                }
                else
                {
                    int p = col.Max;
                    pivotToCol[p] = pos;
                    var b = births[p];
                    births.Remove(p);
                    Emit(b.Arrow, b.ByAdd, arrow, dim[p], dByDelete: false);      // death at an Add -> Open
                }
            }
            else  // Delete
            {
                int pos = posOf[cell];
                if (pos != r.Count - 1)
                    throw new NotSupportedException(
                        "ReflectionZigzag (stage 2b-i) removes the maximal cell at the last position " +
                        "(reverse teardown). Arbitrary removal needs the vineyard transpose-to-last (stage 2b-ii).");

                if (r[pos].Count == 0)
                {
                    var b = births[pos];                                          // unpaired positive -> death
                    births.Remove(pos);
                    Emit(b.Arrow, b.ByAdd, arrow, dim[pos], dByDelete: true);     // death at a Delete -> Closed
                }
                else
                {
                    int p = r[pos].Max;                                           // negative -> rebirth its partner
                    pivotToCol.Remove(p);
                    births[p] = (arrow, false);                                   // birth (backward) -> BirthEnd Open
                }

                r.RemoveAt(pos); dim.RemoveAt(pos); cellAt.RemoveAt(pos); posOf.Remove(cell);
            }
        }

        // Survivors -> infinite bars (run to the end of the filtration; death end Open).
        foreach (var kv in births)
            Emit(kv.Value.Arrow, kv.Value.ByAdd, m, dim[kv.Key], dByDelete: false);

        return new Barcode(bars, "ReflectionZigzag");
    }

    static void Xor(SortedSet<int> col, SortedSet<int> other)
    {
        foreach (int x in other) if (!col.Add(x)) col.Remove(x);
    }
}
