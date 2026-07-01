using System;
using Synthetic.Euclidean;
using Xunit;

namespace Synthetic.Tests;

public sealed class EyeSkeletonTests
{
    [Fact]
    public void BuildSkeleton_CarriesKnobsAndThreeStrokes()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            WarpStrength = 0.3,
            MaxGeodesicRadius = 4.0,
            Center = new[] { 1.0, 2.0, 3.0 },
            CentralPoints = 10, UpperPoints = 20, LowerPoints = 30,
        };
        var sk = EyeTorusToy.BuildSkeleton(cfg);

        Assert.Equal(0.3, sk.WarpStrength);
        Assert.Equal(4.0, sk.MaxGeodesicRadius);
        Assert.Equal(new[] { 1.0, 2.0, 3.0 }, sk.Center);
        Assert.Equal(3, sk.Strokes.Count);
        Assert.Equal(new[] { 0, 1, 2 },
            new[] { sk.Strokes[0].Label, sk.Strokes[1].Label, sk.Strokes[2].Label });
        Assert.Equal(10, sk.Strokes[0].PointCount);
        Assert.Equal(2.0 * Math.PI, sk.Strokes[0].ArcEnd, 12); // central is a full sweep
    }

    [Fact]
    public void Config_HyperbolicKnobsDefaultToFlatNeutral()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig();
        Assert.Equal(1.0, cfg.WarpStrength);
        Assert.Equal(double.PositiveInfinity, cfg.MaxGeodesicRadius);
        Assert.Equal(new[] { 0.0, 0.0, 0.0 }, cfg.Center);
    }

    [Fact]
    public void MaxGeodesicRadius_ConfinesStructureWithinShell()
    {
        const double cap = 3.0;
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            MaxGeodesicRadius = cap,
            CentralPoints = 500, UpperPoints = 300, LowerPoints = 500, // lower R=3.8 exceeds cap
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 4,
        };
        var ds = EyeTorusToy.Generate(cfg);

        foreach (var p in ds.Features) // bg ratio 0 → every point is structure
        {
            double r = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            Assert.True(r <= cap + 1e-9, $"point at radius {r} exceeds cap {cap}");
        }
    }

    [Fact]
    public void Center_TranslatesStructure()
    {
        var center = new[] { 5.0, -4.0, 2.0 };
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            Center = center,
            CentralPoints = 2000, UpperPoints = 1, LowerPoints = 1,
            StructureNoiseSigma = 0.0, ZThickness = 0.0, DensityGradientStrength = 0.0,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 6,
        };
        var ds = EyeTorusToy.Generate(cfg);

        // The full central ring averages to its center; z is exact (no jitter / thickness).
        double mx = 0, my = 0, mz = 0;
        for (int i = 0; i < 2000; i++)
        {
            mx += ds.Features[i][0];
            my += ds.Features[i][1];
            mz += ds.Features[i][2];
        }
        mx /= 2000; my /= 2000; mz /= 2000;

        Assert.True(Math.Abs(mx - center[0]) < 0.25, $"mean x {mx} not near {center[0]}");
        Assert.True(Math.Abs(my - center[1]) < 0.25, $"mean y {my} not near {center[1]}");
        Assert.Equal(center[2], mz, 9);
    }
}
