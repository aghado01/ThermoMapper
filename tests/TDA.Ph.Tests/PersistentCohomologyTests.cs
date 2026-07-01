#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>
/// Parity tests: <see cref="PersistentCohomology"/> vs <see cref="PersistentHomology"/>.
/// </summary>
public sealed class PersistentCohomologyTests
{
    static IEnumerable<(double Birth, double Death, int Dim)> Signatures(Barcode bc) =>
        bc.Bars
            .Select(b => (b.Birth, b.Death, b.Dimension))
            .OrderBy(x => x.Dimension)
            .ThenBy(x => x.Birth)
            .ThenBy(x => x.Death);

    static void AssertBarcodeParity(SimplicialFiltration filtration, int maxDimension = int.MaxValue)
    {
        Barcode ph = PersistentHomology.Compute(filtration, maxDimension);
        Barcode pcoh = PersistentCohomology.Compute(filtration, maxDimension);
        Assert.Equal(Signatures(ph), Signatures(pcoh));
    }

    [Fact]
    public void SingleEdge_ParityWithHomology()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(1.0, 0, 1),
        }, "t");

        AssertBarcodeParity(filtration);
    }

    [Fact]
    public void Circle_ParityWithHomology()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
        }, "t");

        AssertBarcodeParity(filtration);
    }

    [Fact]
    public void FilledTriangle_ParityWithHomology()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
            new Simplex(2.0, 0, 1, 2),
        }, "t");

        AssertBarcodeParity(filtration);
    }

    [Fact]
    public void FilledTriangle_MaxDim0_ParityWithHomology()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
            new Simplex(2.0, 0, 1, 2),
        }, "t");

        AssertBarcodeParity(filtration, maxDimension: 0);
    }

    [Fact]
    public void TetrahedronBoundary_H2_ParityWithHomology()
    {
        AssertBarcodeParity(SimplicialFiltrationFixtures.TetrahedronBoundary());
    }

    [Fact]
    public void GridTorus_H1AndH2_ParityWithHomology()
    {
        AssertBarcodeParity(SimplicialFiltrationFixtures.GridTorus(n: 3, m: 3));
    }

    [Fact]
    public void TwoLoopWedge_H1_ParityWithHomology()
    {
        AssertBarcodeParity(SimplicialFiltrationFixtures.TwoLoopWedge());
    }

    [Fact]
    public void Circle_InfiniteH1_CocycleIsBirthEdge()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
        }, "t");

        Barcode barcode = PersistentCohomology.Compute(filtration, representatives: true);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        Assert.True(bar.IsInfinite);
        Assert.NotNull(bar.Cocycle);
        Assert.Single(bar.Cocycle!);
        Assert.Equal(1, filtration.Simplices[bar.Cocycle![0]].Dimension);
        Assert.Equal(bar.Generator, bar.Cocycle![0]);
    }

    [Fact]
    public void FilledTriangle_FiniteH1_CocycleContainsDeathSimplex()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
            new Simplex(2.0, 0, 1, 2),
        }, "t");

        Barcode barcode = PersistentCohomology.Compute(filtration, representatives: true);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        Assert.False(bar.IsInfinite);
        Assert.NotNull(bar.Cocycle);
        Assert.NotEmpty(bar.Cocycle!);
        Assert.Contains(bar.Cocycle!, i => filtration.Simplices[i].Dimension == 2);
    }

    [Fact]
    public void Circle_WithoutRepresentatives_CocycleIsNull()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
        }, "t");

        Barcode barcode = PersistentCohomology.Compute(filtration, representatives: false);

        Assert.All(barcode.Bars, b => Assert.Null(b.Cocycle));
    }
}
