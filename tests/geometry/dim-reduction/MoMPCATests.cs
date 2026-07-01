using System;
using Maths.LinAlg;
using Maths.Geometry.DimReduction;
using Xunit;

namespace Maths.Geometry.Tests;

public sealed class MoMPCATests
{
    /// <summary>
    /// Three clean blocks share a mean ≈ (5, 0) and a high-variance axis ≈ (1, 0); two corrupted
    /// "nodes" sit at mean (5, 20) along the (1, 1) axis. The scale-calibrated product median on
    /// Rᵖ × Gr(r, p) resists the minority of corrupted blocks in BOTH factors — recovering the clean
    /// axis and a mean whose y-coordinate is pulled back toward 0, not toward the outliers' 20.
    /// </summary>
    [Fact]
    public void ComputeMoM_CorruptedBlocks_RecoverCleanMeanAndSubspace()
    {
        const int perBlock = 12;
        double[][] data = new double[5 * perBlock][];
        int idx = 0;
        idx = AddBlock(data, idx, perBlock, 5.0, 0.0, 0.0);    // clean
        idx = AddBlock(data, idx, perBlock, 5.2, 0.3, 0.04);   // clean
        idx = AddBlock(data, idx, perBlock, 4.8, -0.3, -0.04); // clean
        idx = AddBlock(data, idx, perBlock, 5.0, 20.0, 1.0);   // corrupted node
        idx = AddBlock(data, idx, perBlock, 5.0, 20.0, 1.0);   // corrupted node

        PcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 5, numComponents: 1);

        double[] axis = mom.Components[0];
        Assert.True(
            Math.Abs(axis[0]) > 0.9,
            $"Expected the robust axis near (1,0); |x|={Math.Abs(axis[0]):F4} (corrupted (1,1) would give ~0.707).");
        Assert.InRange(mom.Mean[1], -3.0, 3.0);   // robust mean resists the y=20 outlier nodes
    }

    [Fact]
    public void ComputeMoM_SingleBlock_DelegatesToPca()
    {
        var rng = new Random(99);
        double[][] data = new double[24][];
        for (int i = 0; i < data.Length; i++)
            data[i] = new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble() };

        PcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 1, numComponents: 2);

        Assert.Equal(2, mom.Components.Length);
    }

    // A contiguous block of `count` points with mean (meanX, meanY), spread along (1, slope).
    private static int AddBlock(double[][] data, int start, int count, double meanX, double meanY, double slope)
    {
        double mid = (count - 1) / 2.0;
        for (int i = 0; i < count; i++)
        {
            double t = i - mid;
            data[start + i] = new[] { meanX + t, meanY + slope * t };
        }
        return start + count;
    }
}
