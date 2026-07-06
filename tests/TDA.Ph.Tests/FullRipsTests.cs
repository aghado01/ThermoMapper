#nullable enable
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

public sealed class FullRipsTests
{
    [Fact]
    public void Build_UnboundedTriangle_MaterializesFullTwoSkeleton()
    {
        double[][] points =
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ];

        SimplicialFiltration filtration = FullRips.Build(points);

        Assert.Equal("FullRips", filtration.Label);
        Assert.Equal(3, filtration.Simplices.Count(s => s.Dimension == 0));
        Assert.Equal(3, filtration.Simplices.Count(s => s.Dimension == 1));
        Assert.Equal(1, filtration.Simplices.Count(s => s.Dimension == 2));
    }

    [Fact]
    public void Build_ThresholdOmittedDiagonal_LeavesPathOnly()
    {
        double[][] points =
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
        ];

        SimplicialFiltration filtration = FullRips.Build(points, threshold: 1.0);

        Assert.Equal(3, filtration.Simplices.Count(s => s.Dimension == 0));
        Assert.Equal(2, filtration.Simplices.Count(s => s.Dimension == 1));
        Assert.DoesNotContain(filtration.Simplices, s => s.Dimension == 2);

        Barcode barcode = PersistentHomology.Compute(filtration, maxDimension: 1);
        Assert.Equal(1, barcode.Bars.Count(b => b.Dimension == 0 && b.IsInfinite));
        Assert.DoesNotContain(barcode.Bars, b => b.Dimension == 1 && b.IsInfinite);
    }
}
