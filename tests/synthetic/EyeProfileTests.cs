using System;
using Synthetic.Euclidean;
using Xunit;

namespace Synthetic.Tests;

public sealed class EyeProfileTests
{
    [Fact]
    public void Pupil_AddsCenteredBallAndShiftsBackgroundLabel()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            CentralPoints = 100, UpperPoints = 100, LowerPoints = 100,
            PupilPoints = 500, PupilDilation = 0.5,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 8,
        };
        var ds = EyeTorusToy.Generate(cfg);

        Assert.Equal(800, ds.Features.Length);          // 300 strokes + 500 pupil, no bg
        Assert.Equal(4, ds.ClusterCount);               // 3 strokes + pupil
        Assert.Equal(4, (int)ds.Parameters["backgroundLabel"]);

        double bore = 2.5 - 0.6;
        double maxR = 0.5 * bore;                        // dilation 0.5
        int pupilCount = 0;
        for (int i = 0; i < ds.Features.Length; i++)
        {
            if (ds.Labels[i] != 3) continue;             // pupil label = 3
            pupilCount++;
            var p = ds.Features[i];
            double r = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            Assert.True(r <= maxR + 1e-9, $"pupil point at {r} exceeds bore-dilation {maxR}");
        }
        Assert.Equal(500, pupilCount);
    }

    [Fact]
    public void SolidCrossSection_FillsTheTubeInterior()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            CentralCrossSection = CrossSectionShape.Solid,
            CentralMajorR = 2.5, CentralMinorR = 0.6, ZThickness = 1.0,
            CentralPoints = 3000, UpperPoints = 1, LowerPoints = 1,
            StructureNoiseSigma = 0.0, DensityGradientStrength = 0.0,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 9,
        };
        var ds = EyeTorusToy.Generate(cfg);

        double sum = 0.0, min = double.MaxValue;
        int count = 0;
        for (int i = 0; i < 3000; i++) // central points emitted first
        {
            var p = ds.Features[i];
            double rho = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);
            double crossR = Math.Sqrt((rho - 2.5) * (rho - 2.5) + p[2] * p[2]);
            sum += crossR;
            min = Math.Min(min, crossR);
            count++;
        }
        double mean = sum / count;
        Assert.True(mean < 0.5, $"solid mean cross-radius {mean} should sit well inside the shell (0.6)");
        Assert.True(min < 0.1, $"solid tube should be filled to its center; min cross-radius was {min}");
    }

    [Fact]
    public void RibbonCrossSection_IsThinOutOfPlaneAndWideInPlane()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            CrossSection = CrossSectionShape.Ribbon, RibbonThickness = 0.04,
            CentralMajorR = 2.5, CentralMinorR = 0.6, ZThickness = 1.0,
            CentralPoints = 3000, UpperPoints = 1, LowerPoints = 1,
            StructureNoiseSigma = 0.0,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 10,
        };
        var ds = EyeTorusToy.Generate(cfg);

        double maxAbsZ = 0.0, maxRadial = 0.0;
        for (int i = 0; i < 3000; i++)
        {
            var p = ds.Features[i];
            double rho = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);
            maxAbsZ = Math.Max(maxAbsZ, Math.Abs(p[2]));
            maxRadial = Math.Max(maxRadial, Math.Abs(rho - 2.5));
        }
        Assert.True(maxAbsZ < 0.4, $"ribbon should be thin out-of-plane; max |z| was {maxAbsZ}");
        Assert.True(maxRadial > 0.45, $"ribbon should be wide in-plane; max radial offset was {maxRadial}");
    }

    [Fact]
    public void HalfArcTaper_IsAppliedToHalvesNotCentral()
    {
        var sk = EyeTorusToy.BuildSkeleton(new EyeTorusToy.EyeTorusToyConfig { HalfArcTaper = 0.7 });
        Assert.Equal(0.0, sk.Strokes[0].EndTaper); // central never tapers
        Assert.Equal(0.7, sk.Strokes[1].EndTaper); // upper
        Assert.Equal(0.7, sk.Strokes[2].EndTaper); // lower
    }

    [Fact]
    public void HalfArcTaper_ThinsTheArcEndsRelativeToTheMiddle()
    {
        var cfg = new EyeTorusToy.EyeTorusToyConfig
        {
            HalfArcTaper = 1.0,
            UpperMajorR = 2.0, UpperMinorR = 0.4, ZThickness = 1.0,
            CentralPoints = 1, UpperPoints = 4000, LowerPoints = 1,
            StructureNoiseSigma = 0.0, DensityGradientStrength = 0.0,
            BackgroundDensityRatio = 0.0, NoiseGridSize = 16, Seed = 11,
        };
        var ds = EyeTorusToy.Generate(cfg);

        double arcStart = -Math.PI * 0.6, arcSpan = Math.PI * 1.2;
        double endSum = 0, midSum = 0;
        int endN = 0, midN = 0;
        for (int i = 0; i < ds.Features.Length; i++)
        {
            if (ds.Labels[i] != 1) continue; // upper stroke only
            var p = ds.Features[i];
            double theta = Math.Atan2(p[1], p[0]);
            double t = (theta - arcStart) / arcSpan;
            double rho = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);
            double crossR = Math.Sqrt((rho - 2.0) * (rho - 2.0) + p[2] * p[2]);
            if (t < 1.0 / 3.0 || t > 2.0 / 3.0) { endSum += crossR; endN++; }
            else { midSum += crossR; midN++; }
        }
        Assert.True(endN > 0 && midN > 0);
        Assert.True(endSum / endN < midSum / midN,
            $"tapered arc ends ({endSum / endN:G4}) should be thinner than the middle ({midSum / midN:G4})");
    }
}
