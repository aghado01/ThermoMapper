using System;
using System.Text.Json;
using System.Text.Json.Serialization;
using Maths.Regression.GPrior;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The σ² prior as a value (F-C2): it round-trips through JSON under the repo's conventions, compares by value, refuses
/// a parameter its kind does not take whether built or read back, and keeps the spelling BARS's fingerprint reads.
/// </summary>
public sealed class Sigma2PriorTests
{
    private static readonly Sigma2Prior[] Priors =
        [Sigma2Prior.ScaledToResponse(), Sigma2Prior.ScaledToResponse(0.5), Sigma2Prior.Jeffreys, Sigma2Prior.InverseGamma(2.0, 0.05)];

    // Plain options, and the boundary's string enums and snake case with named floats refused.
    private static readonly JsonSerializerOptions[] Conventions =
    [
        new(),
        new() { PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower, Converters = { new JsonStringEnumConverter() } },
    ];

    [Fact]
    public void RoundTripsThroughJson_ByValue()
    {
        foreach (JsonSerializerOptions options in Conventions)
            foreach (Sigma2Prior prior in Priors)
            {
                string json = JsonSerializer.Serialize(prior, options);
                Assert.Equal(prior, JsonSerializer.Deserialize<Sigma2Prior>(json, options));
                // The derived flags are read from the kind, not stored beside it.
                Assert.DoesNotContain("proper", json, StringComparison.OrdinalIgnoreCase);
                Assert.DoesNotContain("reads", json, StringComparison.OrdinalIgnoreCase);
            }
        Assert.Equal("""{"kind":"ScaledToResponse","shape":0.001,"scale":null}""",
                     JsonSerializer.Serialize(Sigma2Prior.ScaledToResponse(), Conventions[1]));
    }

    [Fact]
    public void ComparesByValue()
    {
        Assert.Equal(Sigma2Prior.ScaledToResponse(1e-3), Sigma2Prior.ScaledToResponse());
        Assert.Equal(Sigma2Prior.ScaledToResponse(1e-3).GetHashCode(), Sigma2Prior.ScaledToResponse().GetHashCode());
        Assert.NotEqual(Sigma2Prior.ScaledToResponse(), Sigma2Prior.ScaledToResponse(0.5));
        Assert.NotEqual(Sigma2Prior.InverseGamma(2.0, 0.05), Sigma2Prior.InverseGamma(2.0, 0.5));
        Assert.Equal(Sigma2Prior.Jeffreys, new Sigma2Prior(Sigma2PriorKind.Jeffreys));
    }

    [Fact]
    public void RefusesParametersItsKindDoesNotTake_BuiltOrRead()
    {
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.Jeffreys, shape: 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Sigma2Prior(Sigma2PriorKind.ScaledToResponse));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.ScaledToResponse(-1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.ScaledToResponse(double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.ScaledToResponse(double.PositiveInfinity));
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.ScaledToResponse, 1e-3, 1.0));
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.InverseGamma, 2.0));
        Assert.ThrowsAny<ArgumentException>(() => Sigma2Prior.InverseGamma(2.0, -1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Sigma2Prior((Sigma2PriorKind)7));

        // A manifest is checked on load as a built prior is.
        Assert.Throws<ArgumentOutOfRangeException>(() => JsonSerializer.Deserialize<Sigma2Prior>("""{"Kind":0,"Shape":-1}"""));
        Assert.Throws<ArgumentException>(() => JsonSerializer.Deserialize<Sigma2Prior>("""{"Kind":1,"Shape":0.5}"""));
    }

    [Fact]
    public void ToString_IsTheSpellingTheFingerprintReads()
    {
        Assert.Equal("ScaledToResponse(0.001)", Sigma2Prior.ScaledToResponse().ToString());
        Assert.Equal("Jeffreys", Sigma2Prior.Jeffreys.ToString());
        Assert.Equal("InverseGamma(2, 0.05)", Sigma2Prior.InverseGamma(2.0, 0.05).ToString());
    }

    [Fact]
    public void ResolvesAsBefore()
    {
        Assert.Null(Sigma2Prior.Jeffreys.For(2.0));
        Assert.Null(Sigma2Prior.ScaledToResponse().For(0.0));
        Assert.Equal(1e-3 * 2.0, Sigma2Prior.ScaledToResponse().For(2.0)!.Value.Scale, 15);
        Assert.Equal(0.05, Sigma2Prior.InverseGamma(2.0, 0.05).For(double.NaN)!.Value.Scale);
        Assert.False(Sigma2Prior.Jeffreys.IsProper);
        Assert.True(Sigma2Prior.ScaledToResponse().ReadsResponse && !Sigma2Prior.InverseGamma(2.0, 0.05).ReadsResponse);
    }
}
