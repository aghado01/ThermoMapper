#nullable enable
using System;
using Maths.Topology;
using Xunit;

namespace TDA.Ph.Tests;

public sealed class DiagramMetricsTests
{
    [Fact]
    public void IdenticalFiniteDiagrams_ReturnZero()
    {
        var bars = new[] { new Bar(0.0, 2.0, 1), new Bar(1.0, 3.0, 1) };
        var a = new Barcode(bars);
        var b = new Barcode(bars);

        double w2 = DiagramMetrics.Wasserstein(a, b, dimension: 1, p: 2.0);
        Assert.Equal(0.0, w2, precision: 12);
    }

    [Fact]
    public void SinglePointVsEmpty_W2EqualsDiagonalDistance()
    {
        var a = new Barcode(new[] { new Bar(0.0, 2.0, 0) });
        var b = new Barcode(Array.Empty<Bar>());

        // L∞ distance to diagonal = (Death - Birth) / 2 = 1
        double w2 = DiagramMetrics.Wasserstein(a, b, dimension: 0, p: 2.0);
        Assert.Equal(1.0, w2, precision: 12);
    }

    [Fact]
    public void EssentialBars_MatchedInfiniteOnMismatch_ReturnsZero()
    {
        var a = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });
        var b = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });

        double w = DiagramMetrics.Wasserstein(a, b, dimension: 0);
        Assert.Equal(0.0, w);
    }

    [Fact]
    public void EssentialBars_MismatchedInfiniteOnMismatch_ReturnsInfinity()
    {
        var a = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });
        var b = new Barcode(Array.Empty<Bar>());

        double w = DiagramMetrics.Wasserstein(a, b, dimension: 0);
        Assert.True(double.IsPositiveInfinity(w));
    }

    [Fact]
    public void EssentialBars_MismatchedFinitePenalty_ChargesPerBar()
    {
        var a = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });
        var b = new Barcode(Array.Empty<Bar>());
        var policy = DiagramMetrics.EssentialPolicy.FinitePenalty(perBar: 2.0);

        double w2 = DiagramMetrics.Wasserstein(a, b, dimension: 0, p: 2.0, essential: policy);
        Assert.Equal(2.0, w2, precision: 12);
    }

    [Fact]
    public void Wasserstein_RejectsPBelowOne()
    {
        var a = new Barcode(Array.Empty<Bar>());
        var b = new Barcode(Array.Empty<Bar>());

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            DiagramMetrics.Wasserstein(a, b, dimension: 0, p: 0.5));
    }

    [Fact]
    public void Wasserstein_RejectsPositiveInfinityP()
    {
        var a = new Barcode(Array.Empty<Bar>());
        var b = new Barcode(Array.Empty<Bar>());

        var ex = Assert.Throws<ArgumentOutOfRangeException>(() =>
            DiagramMetrics.Wasserstein(a, b, dimension: 0, p: double.PositiveInfinity));
        Assert.Contains("Bottleneck", ex.Message);
    }

    [Fact]
    public void Bottleneck_ThrowsNotImplemented()
    {
        var a = new Barcode(Array.Empty<Bar>());
        var b = new Barcode(Array.Empty<Bar>());

        Assert.Throws<NotImplementedException>(() =>
            DiagramMetrics.Bottleneck(a, b, dimension: 0));
    }
}
