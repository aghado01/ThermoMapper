using System;
using Maths.Rng;
using Maths.Samplers;
using Xunit;

namespace Synthetic.Tests;

public sealed class SpectralNoiseFieldTests
{
    private static readonly double[] BoxMin = { -5.0, -5.0, -5.0 };
    private static readonly double[] BoxMax = { 5.0, 5.0, 5.0 };

    [Fact]
    public void Generate_FieldIsFiniteAndStraddlesZero()
    {
        var f = SpectralNoiseField.Generate(
            new Xoshiro256PlusPlus(42), 16, SpectralNoiseField.Pink, BoxMin, BoxMax);

        Assert.True(double.IsFinite(f.Max) && double.IsFinite(f.Min));
        Assert.True(f.Max > 0.0, "standardized field should reach above its mean");
        Assert.True(f.Min < 0.0, "standardized field should reach below its mean");

        for (int i = 0; i < 50; i++)
        {
            double v = f.Sample(BoxMin[0] + i * 0.2, 0.1 * i, -0.05 * i);
            Assert.True(double.IsFinite(v));
        }
    }

    [Fact]
    public void Generate_IsDeterministicForSameSeed()
    {
        var a = SpectralNoiseField.Generate(
            new Xoshiro256PlusPlus(7), 16, SpectralNoiseField.Brown, BoxMin, BoxMax);
        var b = SpectralNoiseField.Generate(
            new Xoshiro256PlusPlus(7), 16, SpectralNoiseField.Brown, BoxMin, BoxMax);

        for (double x = -4.0; x <= 4.0; x += 1.3)
            Assert.Equal(a.Sample(x, x, x), b.Sample(x, x, x), 12);
    }

    [Fact]
    public void Generate_BrownIsSmootherThanWhite()
    {
        // Higher β concentrates power at low frequencies → smoother field →
        // smaller mean absolute difference between adjacent grid nodes.
        double rough = MeanAdjacentDiff(SpectralNoiseField.White);
        double smooth = MeanAdjacentDiff(SpectralNoiseField.Brown);
        Assert.True(smooth < rough, $"brown ({smooth:G4}) should be smoother than white ({rough:G4})");
    }

    private static double MeanAdjacentDiff(double beta)
    {
        const int g = 32;
        var f = SpectralNoiseField.Generate(new Xoshiro256PlusPlus(11), g, beta, BoxMin, BoxMax);
        double step = (BoxMax[0] - BoxMin[0]) / (g - 1);
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < g - 1; i++)
        {
            double x0 = BoxMin[0] + i * step;
            sum += Math.Abs(f.Sample(x0 + step, 0.0, 0.0) - f.Sample(x0, 0.0, 0.0));
            count++;
        }
        return sum / count;
    }

    [Fact]
    public void Generate_FourD_FiniteAndDeterministic()
    {
        double[] min = { -3.0, -3.0, -3.0, -3.0 };
        double[] max = { 3.0, 3.0, 3.0, 3.0 };
        var a = SpectralNoiseField.Generate(new Xoshiro256PlusPlus(21), 8, SpectralNoiseField.Pink, min, max);
        var b = SpectralNoiseField.Generate(new Xoshiro256PlusPlus(21), 8, SpectralNoiseField.Pink, min, max);

        Assert.Equal(4, a.Dimension);
        Assert.True(double.IsFinite(a.Max) && double.IsFinite(a.Min));
        Assert.True(a.Max > 0.0 && a.Min < 0.0);

        var p = new[] { 0.5, -1.0, 2.0, 0.0 };
        Assert.True(double.IsFinite(a.Sample(p)));
        Assert.Equal(a.Sample(p), b.Sample(p), 12); // deterministic for same seed
    }

    [Fact]
    public void Generate_NonPowerOfTwoGrid_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            SpectralNoiseField.Generate(
                new Xoshiro256PlusPlus(1), 20, SpectralNoiseField.Pink, BoxMin, BoxMax));
    }
}
