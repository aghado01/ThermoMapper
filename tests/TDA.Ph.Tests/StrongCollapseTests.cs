#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>Strong collapse must reach the core and preserve homology (oracle: PH before == after).</summary>
public sealed class StrongCollapseTests
{
    // Full simplicial complex (every face) from a set of maximal simplices, all at filtration value 0.
    static SimplicialFiltration Complex(int[][] maximal)
    {
        var seen = new HashSet<string>();
        var simplices = new List<Simplex>();
        foreach (var m in maximal)
        {
            int n = m.Length;
            for (int mask = 1; mask < (1 << n); mask++)
            {
                var verts = new List<int>();
                for (int b = 0; b < n; b++) if ((mask & (1 << b)) != 0) verts.Add(m[b]);
                verts.Sort();
                if (seen.Add(string.Join(",", verts)))
                    simplices.Add(new Simplex(0.0, verts.ToArray()));
            }
        }
        return new SimplicialFiltration(simplices);
    }

    static int[] Betti(int[][] maximal, int maxDim)
    {
        var bc = PersistentHomology.Compute(Complex(maximal), maxDim);
        var b = new int[maxDim + 1];
        foreach (var bar in bc.Bars)
            if (bar.IsInfinite && bar.Dimension <= maxDim) b[bar.Dimension]++;
        return b;
    }

    static void AssertBettiPreserved(int maxDim, int[][] maximal) =>
        Assert.Equal(Betti(maximal, maxDim), Betti(StrongCollapse.Core(maximal), maxDim));

    [Fact]
    public void FilledTriangle_CollapsesToPoint()
    {
        var core = StrongCollapse.Core(new[] { new[] { 0, 1, 2 } });
        Assert.Single(core);
        Assert.Single(core[0]); // one vertex
    }

    [Fact]
    public void HollowTriangle_IsMinimal()
    {
        var core = StrongCollapse.Core(new[] { new[] { 0, 1 }, new[] { 1, 2 }, new[] { 0, 2 } });
        Assert.Equal(3, core.Length); // a circle: nothing is dominated
    }

    [Fact] public void Betti_FilledTriangle() => AssertBettiPreserved(2, new[] { new[] { 0, 1, 2 } });
    [Fact] public void Betti_HollowTriangle() => AssertBettiPreserved(2, new[] { new[] { 0, 1 }, new[] { 1, 2 }, new[] { 0, 2 } });
    [Fact] public void Betti_TwoTrianglesDisk() => AssertBettiPreserved(2, new[] { new[] { 0, 1, 2 }, new[] { 1, 2, 3 } });
    [Fact] public void Betti_HollowSquare() => AssertBettiPreserved(2, new[] { new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 0, 3 } });
    [Fact] public void Betti_FilledTetrahedron() => AssertBettiPreserved(3, new[] { new[] { 0, 1, 2, 3 } });

    [Fact]
    public void Retraction_MapsEveryVertexToTheCore()
    {
        var (core, r) = StrongCollapse.CoreWithRetraction(new[] { new[] { 0, 1, 2 } });
        var coreVerts = core.SelectMany(s => s).ToHashSet();
        Assert.Single(coreVerts);                 // filled triangle -> one core vertex
        int c = coreVerts.First();
        Assert.Equal(c, r[0]);
        Assert.Equal(c, r[1]);
        Assert.Equal(c, r[2]);
    }

    [Fact]
    public void Retraction_IsIdentityOnMinimalComplex()
    {
        var (_, r) = StrongCollapse.CoreWithRetraction(new[] { new[] { 0, 1 }, new[] { 1, 2 }, new[] { 0, 2 } });
        Assert.Equal(0, r[0]);
        Assert.Equal(1, r[1]);
        Assert.Equal(2, r[2]);
    }
}
