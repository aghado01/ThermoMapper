#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// ReflectionZigzag (Maria–Oudot) — stage 1 (forward arrows) + stage 2b-i (reverse-teardown backward arrows).
/// A forward-only zigzag is ordinary persistence; an up-down filtration additionally exercises the
/// backward-arrow birth/death rule (the removed cell is at the last position, no transposition). The new oracle
/// must agree with the trusted general oracles <see cref="ZigzagBarcodeNaive"/> (Z1) and <see cref="FastZigzag"/>
/// (Z2) across every dimension on both. (<c>AssertForwardParity</c> is the generic Z3==Z1==Z2 check despite its
/// name.) Arbitrary removal (vineyard transpose-to-last) is stage 2b-ii; full ∅→∅ churn is stage 3.
/// </summary>
public sealed class ReflectionZigzagTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc, int dim) =>
        bc.Bars.Where(b => b.Dimension == dim)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertForwardParity(ZigzagFiltration f, int maxDim)
    {
        for (int d = 0; d <= maxDim; d++)
        {
            var z1 = Sig(ZigzagBarcodeNaive.Compute(f, d), d).ToList();
            var z2 = Sig(FastZigzag.Compute(f, d), d).ToList();
            var refl = Sig(ReflectionZigzag.Compute(f), d).ToList();
            Assert.Equal(z1, z2);     // the two trusted oracles agree (sanity)
            Assert.Equal(z1, refl);   // the new Maria–Oudot oracle matches them
        }
    }

    static ZigzagFiltration TriangleBuildup()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Add(6, new[] { 3, 4, 5 });
        return f;
    }

    static ZigzagFiltration HollowTetraBuildup()
    {
        var f = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 0, 3 });
        f.Add(7, new[] { 1, 2 }); f.Add(8, new[] { 1, 3 }); f.Add(9, new[] { 2, 3 });
        f.Add(10, new[] { 4, 5, 7 }); f.Add(11, new[] { 4, 6, 8 });
        f.Add(12, new[] { 5, 6, 9 }); f.Add(13, new[] { 7, 8, 9 });
        return f;
    }

    static ZigzagFiltration FilledTetraBuildup()
    {
        var f = HollowTetraBuildup();
        f.Add(14, new[] { 10, 11, 12, 13 });   // 3-simplex fills the void (kills H2)
        return f;
    }

    // Up-down filtrations (build up, tear down in REVERSE order) — each removed cell is at the last position,
    // so stage 2b-i handles them without transpositions. These exercise the backward-arrow birth/death rule.
    static ZigzagFiltration TriangleFillUpDown()
    {
        var f = TriangleBuildup();
        f.Delete(6); f.Delete(5); f.Delete(4); f.Delete(3); f.Delete(2); f.Delete(1); f.Delete(0);
        return f;
    }

    static ZigzagFiltration HollowTetraFillUpDown()
    {
        var f = FilledTetraBuildup();
        f.Delete(14);
        f.Delete(13); f.Delete(12); f.Delete(11); f.Delete(10);
        f.Delete(9); f.Delete(8); f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        for (int v = 3; v >= 0; v--) f.Delete(v);
        return f;
    }

    [Fact] public void TriangleBuildup_MatchesOracles() => AssertForwardParity(TriangleBuildup(), 1);
    [Fact] public void HollowTetraBuildup_MatchesOracles() => AssertForwardParity(HollowTetraBuildup(), 2);
    [Fact] public void FilledTetraBuildup_MatchesOracles() => AssertForwardParity(FilledTetraBuildup(), 2);
    [Fact] public void TriangleFillUpDown_MatchesOracles() => AssertForwardParity(TriangleFillUpDown(), 1);
    [Fact] public void HollowTetraFillUpDown_MatchesOracles() => AssertForwardParity(HollowTetraFillUpDown(), 2);

    [Fact]
    public void ArbitraryRemoval_ThrowsUntilStage2b2()
    {
        // Remove v0 while v1 present: v0 is maximal (no coface) but NOT at the last position -> needs the
        // vineyard transpose-to-last (stage 2b-ii), so stage 2b-i throws.
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Delete(0);
        Assert.Throws<NotSupportedException>(() => ReflectionZigzag.Compute(f));
    }
}
