#nullable enable
using System.Linq;
using Maths.Topology;
using Xunit;

using TDA.Ph;
namespace TDA.Mapper.Tests;

/// <summary>
/// Tests for <see cref="PersistentHomology"/> (standard Z/2Z column reduction).
/// Fixtures: single edge, unfilled triangle (circle), filled triangle.
/// </summary>
public sealed class PersistentHomologyTests
{
    // ── H0: single edge ───────────────────────────────────────────────────────

    // Two vertices at t=0, one edge at t=1.
    // Expected H0 barcode: 2 bars born at 0 — one dies at 1 (vertex killed by edge),
    // one persists. No H1 bars.

    [Fact]
    public void SingleEdge_H0Barcode_OneFiniteOneInfinite()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(1.0, 0, 1),
        }, "t");

        Barcode barcode = PersistentHomology.Compute(filtration);

        var h0 = barcode.Bars.Where(b => b.Dimension == 0).ToList();
        Assert.Equal(2, h0.Count);
        Assert.Single(h0, b => b.IsInfinite);
        Assert.Single(h0, b => !b.IsInfinite && b.Death == 1.0);
        Assert.DoesNotContain(barcode.Bars, b => b.Dimension >= 1);
    }

    // ── H0 + H1: circle (3 vertices, 3 edges, no triangle) ───────────────────

    // Vertices at t=0, edges at t=1, no filler.
    // Expected H0: 3 bars born at 0 — 2 die at 1, 1 persists.
    // Expected H1: 1 infinite bar (loop born at t=1 when the circle closes).

    [Fact]
    public void Circle_H0AndH1Barcode_CorrectCounts()
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

        Barcode barcode = PersistentHomology.Compute(filtration);

        var h0 = barcode.Bars.Where(b => b.Dimension == 0).ToList();
        Assert.Equal(3, h0.Count);
        Assert.Single(h0, b => b.IsInfinite);
        Assert.Equal(2, h0.Count(b => !b.IsInfinite && b.Death == 1.0));

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        Assert.Single(h1);
        Assert.True(h1[0].IsInfinite);
        Assert.Equal(1.0, h1[0].Birth);
    }

    // ── H0 + H1 finite: filled triangle ──────────────────────────────────────

    // Vertices at t=0, edges at t=1, triangle at t=2.
    // Expected H0: same as circle (2 finite, 1 infinite).
    // Expected H1: 1 finite bar — loop born at 1, killed by triangle at 2.

    [Fact]
    public void FilledTriangle_H1FiniteBar_BornAtEdgesDiedAtFill()
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

        Barcode barcode = PersistentHomology.Compute(filtration);

        var h0 = barcode.Bars.Where(b => b.Dimension == 0).ToList();
        Assert.Equal(3, h0.Count);
        Assert.Single(h0, b => b.IsInfinite);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        Assert.False(bar.IsInfinite);
        Assert.Equal(1.0, bar.Birth);
        Assert.Equal(2.0, bar.Death);
    }

    // ── maxDimension filter ───────────────────────────────────────────────────

    [Fact]
    public void FilledTriangle_MaxDim0_OnlyH0Bars()
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

        Barcode barcode = PersistentHomology.Compute(filtration, maxDimension: 0);

        Assert.All(barcode.Bars, b => Assert.Equal(0, b.Dimension));
    }
}
