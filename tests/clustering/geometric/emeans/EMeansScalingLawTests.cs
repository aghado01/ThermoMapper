using Xunit;

namespace Clustering.Geometric.EMeans.Tests;

/// <summary>
/// The published laws (Bessa et al. 2026, Table 2), pinned against values worked by hand from the coefficients.
/// </summary>
public class EMeansScalingLawTests
{
    // B = 0.035 · k^0.13 · n^0.88 · T^0.81 before rounding: 28,790 / 50,475 / 222,709 / 57,334.
    [Theory]
    [InlineData(1_000_000, 1024, 2, 32_768)]
    [InlineData(1_000_000, 1024, 4, 65_536)]
    [InlineData(1_000_000, 1024, 25, 262_144)]
    [InlineData(1_000_000, 16_384, 3, 65_536)]
    public void BatchLaw_RoundsToTheNearestPowerOfTwo(int n, int k, int epochs, int expected)
        => Assert.Equal(expected, EMeansScalingLaw.BatchSize(n, k, epochs));

    // Above 2M points the law is divided by 8: 30,659 and 65,101 before rounding.
    [Theory]
    [InlineData(4_000_000, 4096, 5, 32_768)]
    [InlineData(10_000_000, 65_536, 3, 65_536)]
    public void BatchLaw_IsDividedByEightAboveTwoMillionPoints(int n, int k, int epochs, int expected)
        => Assert.Equal(expected, EMeansScalingLaw.BatchSize(n, k, epochs));

    [Fact]
    public void BatchLaw_IsCappedAtN()
        => Assert.Equal(100, EMeansScalingLaw.BatchSize(100, 2, 1000));

    [Fact]
    public void StepLaw_MatchesTheHandWorkedValue()
        => Assert.Equal(0.362779, EMeansScalingLaw.PeakStep(1_000_000, 16_384, 3, 65_536), 5);

    [Fact]
    public void StepLaw_IsCappedAtOne()
    {
        Assert.Equal(1.258330, EMeansScalingLaw.UnclippedPeakStep(1_000_000, 1024, 4, 65_536), 5);
        Assert.Equal(1.0, EMeansScalingLaw.PeakStep(1_000_000, 1024, 4, 65_536));
    }

    [Fact]
    public void Extrapolations_AreEmptyInsideTheFittingGrid()
        => Assert.Empty(EMeansScalingLaw.Extrapolations(1_000_000, 1024, 4, 65_536));

    [Fact]
    public void Extrapolations_NameEachBoundASmallProblemCrosses()
    {
        var crossed = EMeansScalingLaw.Extrapolations(1000, 5, 5, 64);
        Assert.Contains(crossed, entry => entry.StartsWith("n="));
        Assert.Contains(crossed, entry => entry.StartsWith("k="));
        Assert.Contains(crossed, entry => entry.StartsWith("B="));
        Assert.DoesNotContain(crossed, entry => entry.StartsWith("T="));
    }
}
