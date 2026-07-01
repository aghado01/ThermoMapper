#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Caller-supplied dual graph of an ℝᵖ-embedded complex (Dey–Hou §5, Def 2). The A1 path: defers the
/// geometric void reconstruction (ref [10]) and accepts the dual structure directly. Dual vertices =
/// K-voids ∪ p-simplices; a dual edge per (p−1)-simplex joins the two cells incident across it.
/// </summary>
public sealed class DualGraphSpec
{
    /// <summary>Embedding dimension p (≥ 2); the engine computes the (p−1)-th barcode.</summary>
    public int P { get; }

    /// <summary>Dual-vertex ids of the K-voids — present in every G_i. Disjoint from the p-simplex dual ids.</summary>
    public IReadOnlyList<int> VoidDualVertexIds { get; }

    /// <summary>K cell-id of each p-simplex → its dual-vertex id (in G_i iff that p-simplex ∉ K_i).</summary>
    public IReadOnlyDictionary<int, int> PSimplexDualVertex { get; }

    /// <summary>K cell-id of each (p−1)-simplex → the two dual-vertex ids its dual edge joins (in G_i iff ∉ K_i).</summary>
    public IReadOnlyDictionary<int, (int A, int B)> PMinus1SimplexDualEdge { get; }

    public DualGraphSpec(
        int p,
        IReadOnlyList<int> voidDualVertexIds,
        IReadOnlyDictionary<int, int> pSimplexDualVertex,
        IReadOnlyDictionary<int, (int A, int B)> pMinus1SimplexDualEdge)
    {
        if (p < 2) throw new ArgumentOutOfRangeException(nameof(p), "Codimension-one duality requires p ≥ 2.");
        P = p;
        VoidDualVertexIds = voidDualVertexIds ?? throw new ArgumentNullException(nameof(voidDualVertexIds));
        PSimplexDualVertex = pSimplexDualVertex ?? throw new ArgumentNullException(nameof(pSimplexDualVertex));
        PMinus1SimplexDualEdge = pMinus1SimplexDualEdge ?? throw new ArgumentNullException(nameof(pMinus1SimplexDualEdge));
    }
}

/// <summary>
/// Z5d (A1 skeleton) — the (p−1)-th zigzag barcode of an ℝᵖ-embedded complex via Alexander duality
/// (Dey–Hou §5, Prop 22): <c>Pers(H_{p−1}(F)) = Pers(H̃₀(E))</c> for the dual filtration <c>E</c> of the
/// dual graph, reducing a (p−1)-dimensional problem to the existing H0 graph engine (<see cref="GraphZigzag"/>).
///
/// <para>A1 reduction core: the dual graph is supplied (<see cref="DualGraphSpec"/>) rather than
/// reconstructed from an embedding (ref [10], deferred). One connected dual graph for the whole complex —
/// the (p−1)-connected decomposition (Algorithm 3 item C) and per-component union (item F) are not yet
/// applied.</para>
///
/// <para>Three conventions carry the duality (z5d brief):
/// <list type="number">
/// <item>E reverses F's arrows (F-Add of a simplex ⇒ E-Delete of its dual), which <b>flips</b> each bar's
/// Closed/Open ends — the silent-parity trap.</item>
/// <item>E starts non-empty (G₀ is the full dual graph, since K₀ = ∅): a synthetic Add-prefix builds G₀ from ∅
/// (and a teardown suffix keeps E ∅→∅ for the engine). Bars born in the prefix — the spanning <c>[0,m]</c>
/// reduced-homology component and transient build-up — are dropped.</item>
/// <item>Collapsing F's no-op arrows (dims ∉ {p−1, p}) is undone by remapping each surviving bar's birth/death
/// step back onto F's axis.</item>
/// </list>
/// Output bars live in F's step axis with dimension p−1, matching the direct oracles at that dimension.</para>
/// </summary>
public static class CodimensionOneZigzag
{
    public static Barcode Compute(ZigzagFiltration f, DualGraphSpec dual,
                                  GraphZigzagAlgorithm algorithm = GraphZigzagAlgorithm.Reference)
    {
        ArgumentNullException.ThrowIfNull(f);
        ArgumentNullException.ThrowIfNull(dual);

        int p = dual.P;
        int m = f.Count;

        // E-cell-id allocation: dual vertices (voids, then p-simplex duals), then dual edges.
        var vEId = new Dictionary<int, int>();    // dual-vertex id -> E cell id
        var eEId = new Dictionary<int, int>();    // (p−1) K cell-id  -> E edge cell id
        int nextE = 0;
        foreach (int vd in dual.VoidDualVertexIds) vEId[vd] = nextE++;
        foreach (var kv in dual.PSimplexDualVertex) vEId[kv.Value] = nextE++;
        foreach (var kv in dual.PMinus1SimplexDualEdge) eEId[kv.Key] = nextE++;

        var e = new ZigzagFiltration();
        var fIndexOf = new List<int>();           // per E-step -> F-step (−1 for synthetic prefix/suffix)

        var presentE = new bool[nextE];           // E-cell currently in G_i — the endpoint-presence guard
        void EAdd(int cell, int[] bnd, int fIdx)
        {
            foreach (int b in bnd)
                if (!presentE[b])
                    throw new InvalidOperationException(
                        $"Dual edge (E-cell {cell}) re-enters at F-step {fIdx} but endpoint E-cell {b} is " +
                        "absent: malformed zigzag — a (p−1)-simplex was removed while an incident p-simplex was " +
                        "still present (filtration validity forbids deleting a face under a present coface).");
            presentE[cell] = true;
            e.Add(cell, bnd);
            fIndexOf.Add(fIdx);
        }
        void EDel(int cell, int fIdx) { presentE[cell] = false; e.Delete(cell); fIndexOf.Add(fIdx); }

        // (2) prefix: build G₀ = full dual graph from ∅ (all vertices, then all edges).
        foreach (int vd in dual.VoidDualVertexIds) EAdd(vEId[vd], Array.Empty<int>(), -1);
        foreach (var kv in dual.PSimplexDualVertex) EAdd(vEId[kv.Value], Array.Empty<int>(), -1);
        foreach (var kv in dual.PMinus1SimplexDualEdge)
        {
            var (a, b) = kv.Value;
            EAdd(eEId[kv.Key], new[] { vEId[a], vEId[b] }, -1);
        }
        int prefixLen = e.Count;

        // (1) real steps: walk F, emit the reversed dual op (F-Add ⇒ E-Delete, F-Delete ⇒ E-Add).
        for (int i = 0; i < m; i++)
        {
            var step = f[i];
            int c = step.GlobalCellId;
            bool fAdd = step.Direction == ZigzagDirection.Add;

            if (dual.PSimplexDualVertex.TryGetValue(c, out int dv))
            {
                if (fAdd) EDel(vEId[dv], i);
                else EAdd(vEId[dv], Array.Empty<int>(), i);
            }
            else if (dual.PMinus1SimplexDualEdge.TryGetValue(c, out var ends))
            {
                if (fAdd) EDel(eEId[c], i);
                else EAdd(eEId[c], new[] { vEId[ends.A], vEId[ends.B] }, i);
            }
            // else: dim ∉ {p−1, p} -> no-op arrow, collapsed (no E-step).
        }
        int realEnd = e.Count;

        // (2) suffix: tear G_m down so E is ∅→∅ for the engine (edges first, then vertices).
        // G_m is the full dual graph for an ∅→∅ F; the resulting bars are spanning artifacts, dropped below.
        foreach (var kv in dual.PMinus1SimplexDualEdge) EDel(eEId[kv.Key], -1);
        foreach (var kv in dual.PSimplexDualVertex) EDel(vEId[kv.Value], -1);
        foreach (int vd in dual.VoidDualVertexIds) EDel(vEId[vd], -1);

        // (B) reduce: H0 of E via the existing engine. The dual graph is a multigraph (parallel
        // (p−1)-simplices across one cell pair); both engines are now multigraph-safe — Reference via
        // edge-multiset adjacency, Fast via per-vertex-pair multiplicity over DynamicConnectivity — so either
        // path is valid. Default Reference (the oracle the A1/A2 fixtures pin); pass Fast for near-linear.
        var h0 = GraphZigzag.Compute(e, 0, algorithm: algorithm);
        int totalE = e.Count;

        var outBars = new List<Bar>();
        foreach (var bar in h0.Bars)
        {
            int birthStep = (int)bar.Birth;       // creating E-step (b−1)
            // (2) keep only features born at a real E-step; drops spanning [0,m] + prefix/suffix build-up.
            if (birthStep < prefixLen || birthStep >= realEnd) continue;

            int deathStep = (int)bar.Death;       // killing E-step, or totalE if it survives E
            int fBirth = fIndexOf[birthStep];
            int fDeath = (deathStep < realEnd) ? fIndexOf[deathStep] : m;

            // (1) the reversal flips each end.
            outBars.Add(new Bar(fBirth, fDeath, p - 1, null, null, null,
                Flip(bar.BirthEnd), Flip(bar.DeathEnd)));
        }

        return new Barcode(outBars, "Zigzag Step");
    }

    static IntervalEnd Flip(IntervalEnd end) =>
        end == IntervalEnd.Closed ? IntervalEnd.Open : IntervalEnd.Closed;
}
