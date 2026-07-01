using System;
using Synthetic.Euclidean;
using Xunit;

namespace Synthetic.Tests;

public sealed class EyeTorusToyTests
{
    [Fact]
    public void Generate_ProducesExactCountsAndLabelInvariants()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            CentralPoints = 400,
            UpperPoints = 250,
            LowerPoints = 350,
            BackgroundDensityRatio = 0.1,
            NoiseGridSize = 16,
            Seed = 3,
        };
        var ds = EyeTorusToy.Generate(cfg);

        int structure = 400 + 250 + 350;
        int bg = (int)(structure * 0.1);
        Assert.Equal(structure + bg, ds.Features.Length);
        Assert.Equal(structure + bg, ds.Labels.Length);

        foreach (var p in ds.Features)
        {
            Assert.Equal(3, p.Length);
            Assert.All(p, c => Assert.True(double.IsFinite(c)));
        }

        // Two-level ground truth: level 0 = signal/background, level 1 = fine.
        Assert.NotNull(ds.LabelsByLevel);
        Assert.Equal(2, ds.LabelsByLevel!.Length);
        var coarse = ds.LabelsByLevel[0];
        var fine = ds.LabelsByLevel[1];
        Assert.Same(ds.Labels, fine); // primary labels are the fine level

        int bgFine = 0, bgCoarse = 0;
        for (int i = 0; i < fine.Length; i++)
        {
            if (fine[i] == 3) bgFine++;
            if (coarse[i] == 1) bgCoarse++;
        }
        Assert.Equal(bg, bgFine);
        Assert.Equal(bg, bgCoarse);

        Assert.Equal(3, ds.ClusterCount);
        Assert.Equal(3, ds.Metadata!.GTNumClusters);
        Assert.Equal(3, ds.Metadata.AmbientDimensionality);
    }

    [Fact]
    public void Generate_IsDeterministicForSameSeed()
    {
        EyeTorusToy.EyeTorusToyConfig Cfg() => new()
        {
            Seed = 99, NoiseGridSize = 16,
            CentralPoints = 200, UpperPoints = 150, LowerPoints = 150,
        };
        var a = EyeTorusToy.Generate(Cfg());
        var b = EyeTorusToy.Generate(Cfg());

        Assert.Equal(a.Features.Length, b.Features.Length);
        for (int i = 0; i < a.Features.Length; i++)
            for (int d = 0; d < 3; d++)
                Assert.Equal(a.Features[i][d], b.Features[i][d], 12);
    }

    [Fact]
    public void OverrideSeed_ChangesPointCloud()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            NoiseGridSize = 16, CentralPoints = 200, UpperPoints = 100, LowerPoints = 100,
        };
        var a = EyeTorusToy.Generate(cfg, overrideSeed: 1);
        var b = EyeTorusToy.Generate(cfg, overrideSeed: 2);

        bool anyDiff = false;
        for (int i = 0; i < a.Features.Length && !anyDiff; i++)
            if (Math.Abs(a.Features[i][0] - b.Features[i][0]) > 1e-9)
                anyDiff = true;
        Assert.True(anyDiff);
    }

    [Fact]
    public void FlattenToPlane_DropsToTwoDimensionsPreservingLabels()
    {
        var ds = EyeTorusToy.Generate(new EyeTorusToy.EyeTorusToyConfig
        {
            NoiseGridSize = 16, CentralPoints = 300, UpperPoints = 200, LowerPoints = 200,
        });
        var flat = EyeTorusToy.FlattenToPlane(ds);

        Assert.Equal(ds.Features.Length, flat.Features.Length);
        Assert.All(flat.Features, p => Assert.Equal(2, p.Length));
        Assert.Same(ds.Labels, flat.Labels);
        Assert.Equal(2, flat.Metadata!.AmbientDimensionality);
        Assert.Equal("orthographic-xy", (string)flat.Parameters["projection"]);

        for (int i = 0; i < ds.Features.Length; i++)
        {
            Assert.Equal(ds.Features[i][0], flat.Features[i][0]);
            Assert.Equal(ds.Features[i][1], flat.Features[i][1]);
        }
    }

    [Fact]
    public void CentralRing_PointsLieNearMajorRadius()
    {
        // With a thin tube, zero jitter, and no z-extent, central-ring points sit
        // at planar radius ≈ CentralMajorR ± CentralMinorR from the z-axis.
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            CentralPoints = 1000, UpperPoints = 1, LowerPoints = 1,
            CentralMajorR = 2.5, CentralMinorR = 0.1,
            StructureNoiseSigma = 0.0, DensityGradientStrength = 0.0, ZThickness = 0.0,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 5,
        };
        var ds = EyeTorusToy.Generate(cfg);

        for (int i = 0; i < 1000; i++) // central points are emitted first
        {
            var p = ds.Features[i];
            double rho = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);
            Assert.InRange(rho, 2.5 - 0.11, 2.5 + 0.11);
        }
    }
}
