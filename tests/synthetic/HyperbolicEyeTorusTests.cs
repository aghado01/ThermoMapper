using System;
using Synthetic.Euclidean;
using Synthetic.Manifolds;
using Xunit;

namespace Synthetic.Tests;

public sealed class HyperbolicEyeTorusTests
{
    private static EyeTorusToy.EyeTorusToyConfig SmallConfig(double? rhoMax = null, double warp = 1.0) => new()
    {
        CentralPoints = 800, UpperPoints = 500, LowerPoints = 700,
        BackgroundDensityRatio = 0.15, NoiseGridSize = 16,
        MaxGeodesicRadius = rhoMax ?? double.PositiveInfinity,
        WarpStrength = warp, Seed = 5,
    };

    [Fact]
    public void Generate_AllPointsInsideUnitBall()
    {
        var ds = HyperbolicEyeTorus.Generate(SmallConfig());
        foreach (var p in ds.Features)
        {
            double r = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            Assert.True(r < 1.0, $"point at ball-radius {r} is outside the unit ball");
            Assert.Equal(3, p.Length);
            Assert.All(p, c => Assert.True(double.IsFinite(c)));
        }
    }

    [Fact]
    public void Generate_StructureStaysWithinTheGeodesicShell()
    {
        double rhoMax = 2.0;
        double outerBallR = Math.Tanh(rhoMax / 2.0);
        var ds = HyperbolicEyeTorus.Generate(SmallConfig(rhoMax: rhoMax, warp: 0.5));
        int bgLabel = (int)ds.Parameters["backgroundLabel"];

        for (int i = 0; i < ds.Features.Length; i++)
        {
            if (ds.Labels[i] == bgLabel) continue; // background is allowed out to the shell too
            var p = ds.Features[i];
            double r = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            Assert.True(r <= outerBallR + 1e-9, $"structure point at {r} exceeds shell {outerBallR}");
        }
    }

    [Fact]
    public void Generate_IsDeterministicForSameSeed()
    {
        var a = HyperbolicEyeTorus.Generate(SmallConfig(rhoMax: 2.5));
        var b = HyperbolicEyeTorus.Generate(SmallConfig(rhoMax: 2.5));
        Assert.Equal(a.Features.Length, b.Features.Length);
        for (int i = 0; i < a.Features.Length; i++)
            for (int d = 0; d < 3; d++)
                Assert.Equal(a.Features[i][d], b.Features[i][d], 12);
    }

    [Fact]
    public void WarpStrength_PushesInnerStructureOutwardWhenFaithful()
    {
        // The central ring is inner structure; conformal warp (warp=1) places it
        // at a larger apparent ball-radius than the cosmetically-Euclidean (warp=0)
        // realization — that radial bunching is the distortion the dial exposes.
        double cosmetic = MeanCentralBallRadius(warp: 0.0);
        double faithful = MeanCentralBallRadius(warp: 1.0);
        Assert.True(faithful > cosmetic + 1e-6,
            $"faithful central radius ({faithful:G4}) should exceed cosmetic ({cosmetic:G4})");
    }

    [Fact]
    public void FlattenToPlane_NaivelyProjectsToTwoDimensions()
    {
        var ds = HyperbolicEyeTorus.Generate(SmallConfig(rhoMax: 2.5));
        var flat = EyeTorusToy.FlattenToPlane(ds);
        Assert.Equal(ds.Features.Length, flat.Features.Length);
        Assert.All(flat.Features, p => Assert.Equal(2, p.Length));
        Assert.Equal(2, flat.Metadata!.AmbientDimensionality);
    }

    private static double MeanCentralBallRadius(double warp)
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            CentralPoints = 2000, UpperPoints = 1, LowerPoints = 1,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16,
            MaxGeodesicRadius = 2.5, WarpStrength = warp, Seed = 7,
        };
        var ds = HyperbolicEyeTorus.Generate(cfg);
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < ds.Features.Length; i++)
        {
            if (ds.Labels[i] != 0) continue; // central ring
            var p = ds.Features[i];
            sum += Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            count++;
        }
        return sum / count;
    }
}
