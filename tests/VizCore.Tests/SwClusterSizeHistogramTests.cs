using System;
using System.IO;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Maths.Information;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The log-binned cluster-size histogram counts in 64 bits. At high T nearly every node is a singleton, so bin 0 gains
/// about N a cycle, and a 32-bit count wrapped after about 2³¹/N measured cycles (about 21,000 at N = 100,000); pooling
/// replicas wrapped sooner. Each test carries counts past <see cref="int.MaxValue"/> through one leg of the data path.
/// </summary>
public sealed class SwClusterSizeHistogramTests
{
    [Fact]
    public void Serialization_RoundTripsCountsPastInt32()
    {
        Accumulator acc = Frame(1.0, new[] { int.MaxValue + 10L, 3L });

        using var ms = new MemoryStream();
        AccumulatorSerializer.Instance.WriteTo(acc, ms);
        ms.Position = 0;
        Accumulator round = AccumulatorSerializer.Instance.ReadFrom(ms);

        Assert.Equal(acc.ClusterSizeHistogram, round.ClusterSizeHistogram);
    }

    [Fact]
    public void LabelEntropyCurve_PoolsReplicasPastInt32()
    {
        // Each replica's bins fit in 32 bits; their pooled bin 0, 3·2³⁰, does not. Pooled shares are 3/4 and 1/4.
        long[] replica = [3L << 29, 1L << 29];
        LabelEntropyCurve curve = LabelEntropyCurve.From([Frame(1.0, replica), Frame(1.0, replica)]);

        double expected = -(0.75 * Math.Log(0.75) + 0.25 * Math.Log(0.25));
        Assert.InRange(Math.Abs(curve.Entropy[0] - expected), 0.0, 1e-12);
    }

    [Fact]
    public void Shannon_EntropyNats_CountsPastInt32()
    {
        double expected = -(0.75 * Math.Log(0.75) + 0.25 * Math.Log(0.25));

        Assert.InRange(Math.Abs(Shannon.EntropyNats(new[] { 3L << 31, 1L << 31 }) - expected), 0.0, 1e-12);
    }

    private static Accumulator Frame(double temperature, long[] histogram) => new()
    {
        Temperature = temperature, Q = 4, DrawCount = 1,
        Spins = new int[2], ClusterSizeHistogram = histogram,
        RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
        RunningSumSqClusterSizes = 0, RunningSumSqClusterSizesExcl = 0,
        RunningSumEnergy = 0, RunningSumEnergySq = 0,
        RunningSumMag = 0, RunningSumMagSq = 0,
    };
}
