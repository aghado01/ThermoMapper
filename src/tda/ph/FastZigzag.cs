#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// FastZigzag (Dey–Hou, arXiv:2204.11080): zigzag persistence computed by reducing the zigzag
/// filtration to a single coned standard (non-zigzag) Δ-complex filtration, running standard
/// persistence, and mapping each finite pair back to a 4-type zigzag bar via the creator/destroyer
/// correspondence. Validated against the <see cref="ZigzagBarcodeNaive"/> oracle.
/// <para>With representatives, each bar's <see cref="Bar.Cycle"/> is a cycle of the complex at the bar's first level,
/// as global cell ids, whose class there is the bar's own. A later level may have lost its cells: this names the class
/// where it begins. <see cref="ComputeCycles"/> gives the compatible cycle at every level.</para>
/// </summary>
public static class FastZigzag
{
    private enum EventKind { Add, Delete }

    private readonly struct CellEvent
    {
        public readonly int Index;     // F̂ step index of the add/delete (padding deletions are >= numSteps)
        public readonly EventKind Kind;
        public readonly int Dim;       // dimension of the ORIGINAL simplex (not the coned cell)
        public readonly int Incarnation;
        public CellEvent(int index, EventKind kind, int dim, int incarnation) { Index = index; Kind = kind; Dim = dim; Incarnation = incarnation; }
    }

    public static Barcode Compute(ZigzagFiltration f, int maxDimension = int.MaxValue, bool representatives = false) =>
        Run(f, maxDimension, representatives, keepChains: false).Barcode;

    /// <summary>The barcode with each bar's compatible cycle at every level it covers, read from the coned run's <c>V</c>, at
    /// about the allocations of <see cref="Compute"/> with representatives; figures in <c>docs/TDA/ZigzagPersistence.md</c>.</summary>
    public static ZigzagCycles ComputeCycles(ZigzagFiltration f, int maxDimension = int.MaxValue) =>
        Run(f, maxDimension, representatives: true, keepChains: true);

    static ZigzagCycles Run(ZigzagFiltration f, int maxDimension, bool representatives, bool keepChains)
    {
        int numSteps = f.Count;
        if (numSteps == 0)
            return new ZigzagCycles(new Barcode(Array.Empty<Bar>(), "FastZigzag"), Array.Empty<int[]>(), Array.Empty<int>(), Array.Empty<int>(), Array.Empty<int>(), new());

        // ---- Replay F into per-INCARNATION accounting. Each insertion of a cell is a distinct
        // incarnation: a cell deleted and re-added gets two lifetimes, and the coned filtration must see them
        // as separate simplices. Keying by global cell id instead (one slot per cell, last-add/last-delete
        // wins) conflated the lifetimes into one column and merged their bars — the re-entry bug. ----
        int maxCellId = -1;
        foreach (var s in f) if (s.GlobalCellId > maxCellId) maxCellId = s.GlobalCellId;
        int U = maxCellId + 1;

        var current = new int[U];        // cell id -> its currently-present incarnation, or -1
        for (int i = 0; i < U; i++) current[i] = -1;

        var dimI = new List<int>();      // incarnation -> dimension of the original simplex
        var bndI = new List<int[]>();    // incarnation -> boundary, as incarnation ids
        var addStepI = new List<int>();  // incarnation -> F step added
        var delStepI = new List<int>();  // incarnation -> F step deleted (padding deletions are >= numSteps)
        var origCellI = new List<int>(); // incarnation -> original cell id
        var delList = new List<int>();   // incarnations in deletion order (real deletions, in F order)

        for (int i = 0; i < numSteps; i++)
        {
            var s = f[i];
            int cell = s.GlobalCellId;
            if (s.Direction == ZigzagDirection.Add)
            {
                var bnd = s.BoundaryAtAdd ?? Array.Empty<int>();
                var bi = new int[bnd.Length];
                for (int k = 0; k < bnd.Length; k++) bi[k] = current[bnd[k]]; // boundary -> current incarnations
                int inc = dimI.Count;
                dimI.Add(bnd.Length > 0 ? dimI[bi[0]] + 1 : 0);
                bndI.Add(bi);
                addStepI.Add(i);
                delStepI.Add(-1);
                origCellI.Add(cell);
                current[cell] = inc;
            }
            else
            {
                int inc = current[cell];
                delStepI[inc] = i;
                delList.Add(inc);
                current[cell] = -1;
            }
        }
        int numInc = dimI.Count;

        // ---- Pad: cone every survivor (still-present incarnation) by decreasing dim so cofaces precede
        //  faces. Padding steps get indices >= numSteps. ----
        var survivors = new List<int>();
        for (int c = 0; c < U; c++) if (current[c] != -1) survivors.Add(current[c]);
        survivors.Sort((a, b) => dimI[b].CompareTo(dimI[a])); // decreasing dim
        int padIdx = numSteps;
        foreach (int inc in survivors) { delStepI[inc] = padIdx++; delList.Add(inc); }

        // ---- Build the coned filtration Ê (Algorithm 3.1), over incarnations. Birth = column index. ----
        const int omega = 0;
        var cells = new List<(int Dim, double Birth, int[] Boundary)> { (0, 0, Array.Empty<int>()) };
        var cellEvent = new List<CellEvent?> { null }; // ω has no event
        var origOf = new List<int> { -1 };             // Ê column -> original cell id (-1 = ω / cone)
        var cid = new int[numInc];     // incarnation -> column of its added cell
        var coneId = new int[numInc];  // incarnation -> column of its coned cell
        for (int i = 0; i < numInc; i++) { cid[i] = -1; coneId[i] = -1; }

        // Added incarnations, in F add-order (= incarnation-id order).
        for (int inc = 0; inc < numInc; inc++)
        {
            var bnd = bndI[inc];
            var col = new int[bnd.Length];
            for (int k = 0; k < bnd.Length; k++) col[k] = cid[bnd[k]];
            int colIdx = cells.Count;
            cells.Add((dimI[inc], colIdx, col));
            cellEvent.Add(new CellEvent(addStepI[inc], EventKind.Add, dimI[inc], inc));
            origOf.Add(origCellI[inc]);
            cid[inc] = colIdx;
        }

        // Coned incarnations, in REVERSE deletion order (so face-cones precede coface-cones).
        for (int di = delList.Count - 1; di >= 0; di--)
        {
            int inc = delList[di];
            var bnd = bndI[inc];
            var col = new List<int>(bnd.Length + 1) { cid[inc] };
            if (bnd.Length == 0) col.Add(omega);                 // vertex cone ω·v = edge {v, ω}
            else foreach (int tau in bnd) col.Add(coneId[tau]);  // ω·σ boundary = {σ} ∪ {ω·τ}
            int colIdx = cells.Count;
            cells.Add((dimI[inc] + 1, colIdx, col.ToArray()));
            cellEvent.Add(new CellEvent(delStepI[inc], EventKind.Delete, dimI[inc], inc));
            origOf.Add(-1);
            coneId[inc] = colIdx;
        }

        // ---- Standard persistence on Ê (cones raise dimension by one) ----
        var eHat = new CellFiltration(cells, "FastZigzag-coned");
        int stdMax = maxDimension == int.MaxValue ? int.MaxValue : maxDimension + 1;
        ReductionChains? chains = null;
        Barcode std = keepChains
            ? PersistentHomology.ComputeWithChains(eHat, stdMax, out chains)
            : PersistentHomology.Compute(eHat, stdMax, representatives);

        // ---- Map finite pairs back via creator/destroyer ----
        var bars = new List<Bar>();
        var chainsByBirthStep = new Dictionary<int, (bool, int[], int[])>();
        foreach (var bar in std.Bars)
        {
            if (bar.IsInfinite) continue;
            CellEvent? eb = cellEvent[(int)bar.Birth];
            CellEvent? ed = cellEvent[(int)bar.Death];
            if (eb is null || ed is null) continue; // the ω (cone-apex) class

            CellEvent e1 = eb.Value, e2 = ed.Value;
            CellEvent lo = e1.Index <= e2.Index ? e1 : e2; // creator = earlier F̂ event
            CellEvent hi = e1.Index <= e2.Index ? e2 : e1; // destroyer = later F̂ event

            if (lo.Index >= numSteps) continue; // born in padding => artifact

            int dim = lo.Kind == EventKind.Add ? lo.Dim : lo.Dim - 1;
            if (dim < 0) continue;
            if (maxDimension != int.MaxValue && dim > maxDimension) continue;

            double birth = lo.Index;
            IntervalEnd bEnd = lo.Kind == EventKind.Add ? IntervalEnd.Closed : IntervalEnd.Open;

            double death;
            IntervalEnd dEnd;
            if (hi.Index >= numSteps) { death = numSteps; dEnd = IntervalEnd.Open; } // survivor
            else { death = hi.Index; dEnd = hi.Kind == EventKind.Delete ? IntervalEnd.Closed : IntervalEnd.Open; }

            int[]? cyc = null;
            if (representatives)
            {
                // Created by an addition: the pair's reduced column, which holds only incarnations alive at that step
                // (docs/TDA/ZigzagPersistence.md), and ω in H₀, not a cell of the complex. Created by a deletion: the
                // deleted cell's boundary, the one non-zero class in the backward map's kernel (Dey–Hou's Definition 5).
                var oc = new SortedSet<int>();
                if (lo.Kind == EventKind.Add)
                {
                    foreach (int c in bar.Cycle!) if (origOf[c] >= 0) oc.Add(origOf[c]);
                }
                else
                {
                    foreach (int face in bndI[lo.Incarnation]) oc.Add(origCellI[face]);
                }
                cyc = new int[oc.Count];
                oc.CopyTo(cyc);
            }
            bars.Add(new Bar(birth, death, dim, null, null, cyc, bEnd, dEnd));

            if (chains is not null)
            {
                var added = new List<int>();
                var coned = new List<int>();
                foreach (int c in chains.ChainOf((int)bar.Death))
                {
                    CellEvent e = cellEvent[c]!.Value;
                    (e.Kind == EventKind.Add ? added : coned).Add(e.Incarnation);
                }

                chainsByBirthStep[lo.Index] = (lo.Kind == EventKind.Add, added.ToArray(), coned.ToArray());
            }
        }

        return new ZigzagCycles(new Barcode(bars, "FastZigzag"), bndI, addStepI, delStepI, origCellI, chainsByBirthStep);
    }
}
