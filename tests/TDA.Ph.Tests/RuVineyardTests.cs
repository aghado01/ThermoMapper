#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Stage 2a — isolated validation of the RU vineyard transposition (<see cref="RuVineyard"/>). The invariant:
/// a <c>VineSwap(i)</c> must leave the matrix reduced to exactly the pairing that a from-scratch standard
/// reduction of the reordered filtration produces. The from-scratch reduction here is an independent
/// implementation (R only, no U), so this cross-checks the vineyard's R+U bookkeeping — especially the
/// U = V⁻¹ update direction — locally, before it is wired to ReflectionZigzag's backward arrow (stage 2b).
/// </summary>
public sealed class RuVineyardTests
{
    sealed class Complex
    {
        public readonly Dictionary<int, int> Dim = new();
        public readonly Dictionary<int, int[]> Bnd = new();   // boundary (codim-1 faces) by cell id
        public readonly List<int> Order = new();              // a valid filtration order (faces before cofaces)
    }

    static List<(int Dim, int[] Boundary)> Positional(Complex k, List<int> order)
    {
        var pos = new Dictionary<int, int>();
        for (int p = 0; p < order.Count; p++) pos[order[p]] = p;
        return order.Select(id => (k.Dim[id], k.Bnd[id].Select(f => pos[f]).ToArray())).ToList();
    }

    // Independent standard R = DV reduction (R only) -> pairs in positions.
    static List<(int Birth, int Death, int Dim)> StandardPairs(List<(int Dim, int[] Boundary)> cells)
    {
        int n = cells.Count;
        var r = new SortedSet<int>?[n];
        var pivotToCol = new Dictionary<int, int>();
        var negPivot = new int[n];
        for (int i = 0; i < n; i++) negPivot[i] = -1;
        for (int j = 0; j < n; j++)
        {
            var col = new SortedSet<int>(cells[j].Boundary);
            while (col.Count > 0)
            {
                int low = col.Max;
                if (pivotToCol.TryGetValue(low, out int j2)) foreach (int x in r[j2]!) { if (!col.Add(x)) col.Remove(x); }
                else break;
            }
            if (col.Count > 0) { r[j] = col; pivotToCol[col.Max] = j; negPivot[j] = col.Max; }
        }
        var res = new List<(int, int, int)>();
        var paired = new HashSet<int>();
        for (int j = 0; j < n; j++) if (negPivot[j] != -1) { res.Add((negPivot[j], j, cells[negPivot[j]].Dim)); paired.Add(negPivot[j]); }
        for (int p = 0; p < n; p++) if (negPivot[p] == -1 && !paired.Contains(p)) res.Add((p, -1, cells[p].Dim));
        return res;
    }

    // Position pairs -> cell-id pairs (so they compare across orderings).
    static HashSet<(int, int, int)> CellPairs(List<(int Birth, int Death, int Dim)> posPairs, List<int> order)
        => posPairs.Select(t => (order[t.Birth], t.Death == -1 ? -1 : order[t.Death], t.Dim)).ToHashSet();

    // Two adjacent cells may be transposed iff neither is a (codim-1) face of the other.
    static bool ValidSwap(Complex k, List<int> order, int i)
        => !k.Bnd[order[i + 1]].Contains(order[i]) && !k.Bnd[order[i]].Contains(order[i + 1]);

    static Complex Build((int Id, int Dim, int[] Bnd)[] cells)
    {
        var k = new Complex();
        foreach (var (id, d, b) in cells) { k.Dim[id] = d; k.Bnd[id] = b; k.Order.Add(id); }
        return k;
    }

    static Complex Triangle() => Build(new (int, int, int[])[]
    {
        (0, 0, new int[0]), (1, 0, new int[0]), (2, 0, new int[0]),
        (3, 1, new[] { 0, 1 }), (4, 1, new[] { 1, 2 }), (5, 1, new[] { 0, 2 }),
        (6, 2, new[] { 3, 4, 5 }),
    });

    static Complex HollowTetra() => Build(new (int, int, int[])[]
    {
        (0, 0, new int[0]), (1, 0, new int[0]), (2, 0, new int[0]), (3, 0, new int[0]),
        (4, 1, new[] { 0, 1 }), (5, 1, new[] { 0, 2 }), (6, 1, new[] { 0, 3 }),
        (7, 1, new[] { 1, 2 }), (8, 1, new[] { 1, 3 }), (9, 1, new[] { 2, 3 }),
        (10, 2, new[] { 4, 5, 7 }), (11, 2, new[] { 4, 6, 8 }), (12, 2, new[] { 5, 6, 9 }), (13, 2, new[] { 7, 8, 9 }),
    });

    static void AssertSingleTranspositions(Complex k)
    {
        for (int i = 0; i < k.Order.Count - 1; i++)
        {
            if (!ValidSwap(k, k.Order, i)) continue;
            var vy = new RuVineyard(Positional(k, k.Order));
            vy.VineSwap(i);
            var op = new List<int>(k.Order); (op[i], op[i + 1]) = (op[i + 1], op[i]);
            var expected = CellPairs(StandardPairs(Positional(k, op)), op);
            Assert.Equal(expected, CellPairs(vy.Pairs(), op));
        }
    }

    static void AssertRandomSequences(Complex k, int seeds, int steps)
    {
        for (int seed = 0; seed < seeds; seed++)
        {
            var rng = new Random(seed);
            var order = new List<int>(k.Order);
            var vy = new RuVineyard(Positional(k, order));
            for (int s = 0; s < steps; s++)
            {
                var valid = Enumerable.Range(0, order.Count - 1).Where(i => ValidSwap(k, order, i)).ToList();
                if (valid.Count == 0) break;
                int i = valid[rng.Next(valid.Count)];
                vy.VineSwap(i);
                (order[i], order[i + 1]) = (order[i + 1], order[i]);
                var expected = CellPairs(StandardPairs(Positional(k, order)), order);
                Assert.Equal(expected, CellPairs(vy.Pairs(), order));
            }
        }
    }

    [Fact] public void Triangle_SingleTranspositions() => AssertSingleTranspositions(Triangle());
    [Fact] public void HollowTetra_SingleTranspositions() => AssertSingleTranspositions(HollowTetra());
    [Fact] public void Triangle_RandomSequences() => AssertRandomSequences(Triangle(), 20, 15);
    [Fact] public void HollowTetra_RandomSequences() => AssertRandomSequences(HollowTetra(), 20, 30);
}
