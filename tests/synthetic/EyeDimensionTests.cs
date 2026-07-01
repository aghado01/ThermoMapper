using System;
using Synthetic.Euclidean;
using Synthetic.Manifolds;
using Xunit;

namespace Synthetic.Tests;

public sealed class EyeDimensionTests
{
    [Fact]
    public void FourD_StrokesOccupyIndependentPlanes()
    {
        // central sweeps (0,1)/axial 2, upper (0,2)/axial 3, lower (0,3)/axial 1.
        // With no jitter and flat tubes, each stroke's "unused" axis stays exactly 0
        // — the rings live in independent 2-planes (unlinked in 4-D).
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            Dimension = 4,
            CentralPoints = 600, UpperPoints = 400, LowerPoints = 500,
            StructureNoiseSigma = 0.0, ZThickness = 0.0,
            BackgroundDensityRatio = 0.1, NoiseGridSize = 16, Seed = 4,
        };
        var ds = EyeTorusToy.Generate(cfg);

        Assert.Equal(4, (int)ds.Parameters["dimension"]);
        Assert.All(ds.Features, p => Assert.Equal(4, p.Length));

        for (int i = 0; i < ds.Features.Length; i++)
        {
            var p = ds.Features[i];
            switch (ds.Labels[i])
            {
                case 0: Assert.Equal(0.0, p[3], 9); break; // central: dim 3 unused
                case 1: Assert.Equal(0.0, p[1], 9); break; // upper: dim 1 unused
                case 2: Assert.Equal(0.0, p[2], 9); break; // lower: dim 2 unused
            }
        }
    }

    [Fact]
    public void FourD_Hyperbolic_AllPointsInsideUnitBall()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            Dimension = 4,
            CentralPoints = 500, UpperPoints = 300, LowerPoints = 400,
            BackgroundDensityRatio = 0.1, NoiseGridSize = 16,
            MaxGeodesicRadius = 2.0, Seed = 5,
        };
        var ds = HyperbolicEyeTorus.Generate(cfg);

        Assert.Equal(4, (int)ds.Parameters["dimension"]);
        foreach (var p in ds.Features)
        {
            Assert.Equal(4, p.Length);
            double r2 = 0.0;
            foreach (var c in p) { Assert.True(double.IsFinite(c)); r2 += c * c; }
            Assert.True(Math.Sqrt(r2) < 1.0, $"4-D point at ball-radius {Math.Sqrt(r2)} escaped the unit ball");
        }
    }

    [Fact]
    public void Dimension_DefaultsToThree()
    {
        var ds = EyeTorusToy.Generate(new EyeTorusToy.EyeTorusToyConfig
        {
            CentralPoints = 50, UpperPoints = 50, LowerPoints = 50,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16,
        });
        Assert.Equal(3, (int)ds.Parameters["dimension"]);
        Assert.All(ds.Features, p => Assert.Equal(3, p.Length));
    }
}
