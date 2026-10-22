using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The birth density of step 3.9: ε·b(t) + (1 − ε)·(local mixture), and its two bases. Each is a probability density on
/// (0, 1) whose samples follow it; ε = 0 is exactly DMGK's local mixture; and with no knots the birth comes from the base.
/// </summary>
public sealed class KnotBirthDensityTests
{
    // A deliberately uneven design: 40 points packed into [0.6, 0.9].
    private static readonly double[] Clustered = Enumerable.Range(0, 40).Select(i => 0.6 + 0.3 * (i + 0.5) / 40).ToArray();

    private static double Integral(Func<double, double> logDensity)
    {
        const int cells = 200_000;   // midpoint rule on (0, 1)
        double s = 0.0;
        for (int i = 0; i < cells; i++) s += Math.Exp(logDensity((i + 0.5) / cells));
        return s / cells;
    }

    public static TheoryData<double, string> Mixtures => new()
    {
        { 0.0, "uniform" }, { 0.1, "uniform" }, { 0.3, "design" }, { 1.0, "design" },
    };

    private static KnotBirthDensity Density(double exploration, string baseName)
        => new(new LocalBetaKernel(50.0), exploration, baseName == "design" ? new DesignHistogramBase(Clustered) : UniformKnotBase.Instance);

    [Theory]
    [MemberData(nameof(Mixtures))]
    public void Density_IntegratesToOne_AndItsSamplesFollowIt(double exploration, string baseName)
    {
        KnotBirthDensity q = Density(exploration, baseName);
        double[] knots = [0.2, 0.25, 0.7];
        Assert.Equal(1.0, Integral(t => q.LogDensity(t, knots)), 3);
        Assert.Equal(1.0, Integral(t => q.LogDensity(t, [])), 3);

        // Twenty equal cells: sample frequencies against the density's mass in each, within five standard errors.
        var rng = new Xoshiro256PlusPlus(seed: 9);
        const int draws = 200_000, cellsCount = 20;
        var counts = new double[cellsCount];
        for (int d = 0; d < draws; d++) counts[Math.Min((int)(q.Sample(knots, rng) * cellsCount), cellsCount - 1)] += 1.0;
        for (int c = 0; c < cellsCount; c++)
        {
            double mass = 0.0;
            for (int i = 0; i < 2000; i++) mass += Math.Exp(q.LogDensity((c + (i + 0.5) / 2000) / cellsCount, knots));
            mass /= 2000 * cellsCount;
            double observed = counts[c] / draws, se = Math.Sqrt(mass * (1 - mass) / draws);
            Assert.True(Math.Abs(observed - mass) <= 5 * se + 1e-4, $"cell {c}: {observed:F4} against {mass:F4}");
        }
    }

    [Fact]
    public void LocalOnly_IsDmgksMixture_AndEmpty_IsTheBase()
    {
        var kernel = new LocalBetaKernel(50.0);
        var q = new KnotBirthDensity(kernel);
        double[] knots = [0.3, 0.6];
        foreach (double t in new[] { 0.05, 0.31, 0.5, 0.95 })
            Assert.Equal(ProposalMath.LogMixtureDensity(t, knots, kernel), q.LogDensity(t, knots));
        Assert.Equal(0.0, q.LogDensity(0.4, []));

        var design = new DesignHistogramBase(Clustered);
        var qd = new KnotBirthDensity(kernel, 0.1, design);
        Assert.Equal(design.LogDensity(0.75), qd.LogDensity(0.75, []));
        Assert.Equal(Math.Log(design.Floor), design.LogDensity(0.2), 12);   // an empty bin keeps the floor
    }

    [Fact]
    public void Exploration_BoundsTheReverseDensity_OfAnIsolatedKnot()
    {
        // OQ11's trap: from knots at 0.39 and 0.96, the local reverse density at 0.39 is about e^-37. With ε = 0.1 and the
        // uniform base it is at least 0.1.
        double[] remaining = [0.958];
        double local = new KnotBirthDensity(new LocalBetaKernel(50.0)).LogDensity(0.3888, remaining);
        double mixed = new KnotBirthDensity(new LocalBetaKernel(50.0), 0.1).LogDensity(0.3888, remaining);
        Assert.InRange(local, -38.0, -36.0);
        Assert.True(mixed >= Math.Log(0.1));
    }

    [Fact]
    public void Construction_RefusesOutOfRangeSettings()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new KnotBirthDensity(new LocalBetaKernel(50.0), -0.01));
        Assert.Throws<ArgumentOutOfRangeException>(() => new KnotBirthDensity(new LocalBetaKernel(50.0), 1.01));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DesignHistogramBase(Clustered, floor: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DesignHistogramBase([0.5, 1.5]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(),
            new PoissonPrior(5.0), new LocalBetaKernel(50.0), birthExploration: 2.0));
    }
}
