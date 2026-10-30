using System;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The REPL's <c>--accumulation</c> tokens: each sets its own <see cref="AccumulationSpec"/> flag and no other, tokens
/// combine by comma, and an unknown token names every valid one.
/// </summary>
public sealed class SpcAccumulationTokensTests
{
    [Theory]
    [InlineData("affinities")]
    [InlineData("alignments")]
    [InlineData("comembership")]
    [InlineData("cluster-size-landscape")]
    [InlineData("order-landscape")]
    [InlineData("burnin-trace")]
    [InlineData("exact-cluster-size-histogram")]
    public void EachToken_SetsItsFlagAlone(string token)
    {
        AccumulationSpec spec = SpcCommand.ParseAccumulationSpec(token);

        AccumulationSpec expected = token switch
        {
            "affinities"                   => new AccumulationSpec { Affinities = true },
            "alignments"                   => new AccumulationSpec { Alignments = true },
            "comembership"                 => new AccumulationSpec { CoMembership = true },
            "cluster-size-landscape"       => new AccumulationSpec { ClusterSizeLandscape = true },
            "order-landscape"              => new AccumulationSpec { OrderLandscape = true },
            "burnin-trace"                 => new AccumulationSpec { BurnInTrace = true },
            "exact-cluster-size-histogram" => new AccumulationSpec { ExactClusterSizeHistogram = true },
            _ => throw new ArgumentOutOfRangeException(nameof(token)),
        };
        Assert.Equal(expected, spec);
    }

    [Fact]
    public void Tokens_CombineByComma()
    {
        AccumulationSpec spec = SpcCommand.ParseAccumulationSpec("comembership, Burnin-Trace,exact-cluster-size-histogram");

        Assert.Equal(new AccumulationSpec { CoMembership = true, BurnInTrace = true, ExactClusterSizeHistogram = true }, spec);
    }

    [Fact]
    public void None_IsTheDefault()
        => Assert.Equal(AccumulationSpec.None, SpcCommand.ParseAccumulationSpec("none"));

    [Fact]
    public void UnknownToken_NamesEveryValidToken()
    {
        var ex = Assert.Throws<ArgumentException>(() => SpcCommand.ParseAccumulationSpec("burn-in-trace"));

        foreach (string token in new[] { "none", "affinities", "alignments", "comembership", "cluster-size-landscape",
                                         "order-landscape", "burnin-trace", "exact-cluster-size-histogram" })
            Assert.Contains(token, ex.Message);
    }
}
