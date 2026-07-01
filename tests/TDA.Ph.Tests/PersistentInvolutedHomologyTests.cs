#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>
/// Parity and cycle-representative tests for <see cref="PersistentInvolutedHomology"/>.
/// </summary>
public sealed class PersistentInvolutedHomologyTests
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
        Barcode inv = PersistentInvolutedHomology.Compute(filtration, maxDimension);
        Assert.Equal(Signatures(ph), Signatures(inv));
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
    public void TetrahedronBoundary_ParityWithHomology()
    {
        AssertBarcodeParity(SimplicialFiltrationFixtures.TetrahedronBoundary());
    }

    [Fact]
    public void TwoLoopWedge_ParityWithHomology()
    {
        AssertBarcodeParity(SimplicialFiltrationFixtures.TwoLoopWedge());
    }

    [Fact]
    public void Circle_InfiniteH1_CycleContainsEdge()
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

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        Assert.True(bar.IsInfinite);
        Assert.NotNull(bar.Cycle);
        Assert.Contains(bar.Cycle!, i => filtration.Simplices[i].Dimension == 1);
    }

    [Fact]
    public void FilledTriangle_FiniteH1_CycleContainsLoopEdge()
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

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        Assert.False(bar.IsInfinite);
        Assert.NotNull(bar.Cycle);
        Assert.Contains(bar.Cycle!, i => filtration.Simplices[i].Dimension == 1);
        Assert.DoesNotContain(bar.Cycle!, i => filtration.Simplices[i].Dimension == 2);
    }

    [Fact]
    public void WithoutRepresentatives_CycleIsNull()
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

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: false);

        Assert.All(barcode.Bars, b => Assert.Null(b.Cycle));
    }
}
