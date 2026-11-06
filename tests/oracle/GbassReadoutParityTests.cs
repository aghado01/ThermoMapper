using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// GBASS's posterior readouts (<c>gbass_readouts_oracle.R</c>; E1's build brief, slice 0): <c>tbass</c>, <c>qbass</c>, <c>nwbass</c> and
/// <c>hbass</c> on one data set, the fitted surface over a grid with w, β, γ, τ and M. GBASS's chain is not the engine's, so the engine
/// meets these readouts, never their traces, each difference stated beside the comparison; the slice that builds each noise holds the
/// engine against its wrapper here. Slice 0 checks that each mode runs from <c>r/</c> and reads what its noise means: the q-quantile
/// surface has about q of the responses below it, and t's w is the noise variance.
/// </summary>
public sealed class GbassReadoutParityTests(ITestOutputHelper output)
{
    internal const int Observations = 200;

    /// <summary>The readouts' data: 200 cases on [0, 1]², f = sin(2πx₁) + x₂, noise 0.3·t(3), whose variance is 0.27.</summary>
    internal static double[][] Data()
    {
        var rng = new Xoshiro256PlusPlus(20261107);
        var t3 = StudentT.Standard(3.0);
        return [.. Enumerable.Range(0, Observations).Select(_ =>
        {
            double x1 = rng.NextDouble(), x2 = rng.NextDouble();
            return new[] { x1, x2, Math.Sin(2 * Math.PI * x1) + x2 + 0.3 * t3.Sample(rng) };
        })];
    }

    /// <summary>The surface's grid: 5 × 3 points on [0.05, 0.95]².</summary>
    internal static double[][] Grid() => [.. from b in Enumerable.Range(0, 3) from a in Enumerable.Range(0, 5) select new[] { 0.05 + 0.225 * a, 0.05 + 0.45 * b }];

    internal static JsonElement? Readout(string model, int nmcmc, int nburn, params string[] args)
    {
        if (!ROracle.IsAvailable) return null;
        string grid = GbassOracle.WriteCsv(Grid());
        try { return GbassOracle.Run("gbass_readouts_oracle.R", Data(), [model, grid, nmcmc.ToString(), nburn.ToString(), "11", .. args]); }
        finally { System.IO.File.Delete(grid); }
    }

    private static readonly Dictionary<string, Lazy<JsonElement?>> Fits = new()
    {
        ["tbass"] = new(() => Readout("tbass", 2000, 1000, "df=3")),
        ["qbass"] = new(() => Readout("qbass", 2000, 1000, "q=0.9")),
        ["nwbass"] = new(() => Readout("nwbass", 2000, 1000, "s_beta=10", "m_gamma=1", "s_gamma=3", "lag_beta=500")),
        ["hbass"] = new(() => Readout("hbass", 2000, 1000)),
    };

    [Theory]
    [InlineData("tbass")]
    [InlineData("qbass")]
    [InlineData("nwbass")]
    [InlineData("hbass")]
    public void EachMode_RunsFromR_AndEmitsItsReadouts(string model)
    {
        if (Fits[model].Value is not { } r) return;
        Assert.Equal(model, r.GetProperty("model").GetString());
        Assert.Equal(Observations, r.GetProperty("N").GetInt32());
        JsonElement surface = r.GetProperty("surface");
        int points = Grid().Length;
        foreach (string k in new[] { "mean", "q05", "q50", "q95" })
            Assert.Equal(points, GbassOracle.Numbers(surface.GetProperty(k)).Count(double.IsFinite));
        double[] lo = GbassOracle.Numbers(surface.GetProperty("q05")), hi = GbassOracle.Numbers(surface.GetProperty("q95"));
        Assert.All(Enumerable.Range(0, points), i => Assert.True(lo[i] <= hi[i]));
        double[] proposed = GbassOracle.Numbers(r.GetProperty("cnt1")), accepted = GbassOracle.Numbers(r.GetProperty("cnt2"));
        Assert.True(accepted.Sum() > 0 && proposed.Sum() >= accepted.Sum());
        output.WriteLine($"{model}: {r.GetProperty("seconds").GetDouble():F1} s, M {r.GetProperty("M").GetProperty("mean").GetDouble():F2}, " +
                         $"w {r.GetProperty("w").GetProperty("mean").GetDouble():G4}, below the fit {r.GetProperty("below_fit").GetDouble():F3}, " +
                         $"acceptance B/D/M {string.Join("/", accepted.Zip(proposed, (a, p) => (a / p).ToString("F3")))}");
    }

    [Fact]
    public void Qbass_PutsAboutQOfTheResponsesBelowItsSurface()
    {
        if (Fits["qbass"].Value is not { } r) return;
        double below = r.GetProperty("below_fit").GetDouble();
        Assert.InRange(below, 0.9 - 0.06, 0.9 + 0.06);
        Assert.Equal(1.0 / 0.9 - 1.0 / 0.1, r.GetProperty("beta").GetProperty("mean").GetDouble(), 12);   // β held at 1/q − 1/(1 − q)
    }

    [Fact]
    public void Tbass_ReadsWAsTheNoiseVariance()
    {
        if (Fits["tbass"].Value is not { } r) return;
        JsonElement w = r.GetProperty("w");
        output.WriteLine($"tbass w: mean {w.GetProperty("mean").GetDouble():G4}, 90% [{w.GetProperty("q05").GetDouble():G4}, {w.GetProperty("q95").GetDouble():G4}]; the noise's variance 0.27");
        Assert.InRange(0.27, 0.5 * w.GetProperty("q05").GetDouble(), 2.0 * w.GetProperty("q95").GetDouble());
    }
}
