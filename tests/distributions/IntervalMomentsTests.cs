using System;
using Maths.Distributions;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// Interval moments ∫ (x − c)ᵏ dF over (a, b] (<see cref="IIntervalMoments"/>): each measure against a reference that does
/// not share its method. The uniform and the full-line normal against closed forms; the normal's series against a
/// composite Gauss–Legendre rule, on wide, narrow and far-tail intervals; Student t against the Cauchy's closed forms and
/// a fine composite rule at ν = 2 and 3; truncation, mixtures under both truncation conventions, and the empirical measure
/// by enumeration, its atoms counted once across a partition.
/// </summary>
public sealed class IntervalMomentsTests
{
    private const int Powers = 7;

    [Fact]
    public void GaussLegendre_IntegratesPolynomialsToDegreeTwoNMinusOne()
    {
        foreach (int n in new[] { 1, 2, 5, 20, 40 })
        {
            GaussLegendre rule = GaussLegendre.Of(n);
            double total = 0.0;
            foreach (double w in rule.Weights) total += w;
            Assert.Equal(2.0, total, 13);
            for (int k = 0; k <= 2 * n - 1; k++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++) sum += rule.Weights[i] * Math.Pow(rule.Nodes[i], k);
                double exact = k % 2 == 0 ? 2.0 / (k + 1) : 0.0;
                Assert.True(Math.Abs(sum - exact) <= 1e-14, $"n = {n}, x^{k}: {sum} against {exact}");
            }
        }
        Assert.Equal(1.0 / Math.Sqrt(3.0), GaussLegendre.Of(2).Nodes[1], 15);
        Assert.Same(GaussLegendre.Of(20), GaussLegendre.Of(20));
    }

    [Fact]
    public void Uniform_MatchesItsClosedFormClippedToTheSupport()
    {
        Uniform u = Uniform.On(0.25, 0.75);
        var m = new double[Powers];
        Assert.Equal(0.0, u.IntervalMoments(-1.0, 0.4, 0.1, 0.0, m));
        // ∫ over [0.25, 0.4] of (x − 0.1)ᵏ / 0.5.
        for (int k = 0; k < Powers; k++)
            Assert.Equal((Math.Pow(0.3, k + 1) - Math.Pow(0.15, k + 1)) / ((k + 1) * 0.5), m[k], 15);
        u.IntervalMoments(0.8, 2.0, 0.0, 0.0, m);
        Assert.All(m, v => Assert.Equal(0.0, v));
        u.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.5, Math.Log(2.0), m);
        Assert.Equal(0.5, m[0], 15);   // the scale divides by e^logScale
        Assert.Equal(0.3, u.IntervalProbability(0.1, 0.4), 15);
        Assert.Throws<ArgumentException>(() => Uniform.On(1.0, 1.0));
    }

    [Fact]
    public void Normal_OnTheWholeLine_HasTheGaussianMoments()
    {
        Normal f = Normal.MeanSd(0.7, 0.2);
        var m = new double[Powers];
        f.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.7, 0.0, m);
        // Central moments: σᵏ(k − 1)!! at even k, 0 at odd.
        double[] central = [1.0, 0.0, 0.04, 0.0, 3 * 0.0016, 0.0, 15 * 0.000064];
        for (int k = 0; k < Powers; k++) Assert.True(Math.Abs(m[k] - central[k]) <= 1e-15, $"k = {k}: {m[k]}");
        f.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.0, 0.0, m);
        Assert.Equal(0.7, m[1], 14);
        Assert.Equal(0.49 + 0.04, m[2], 14);
    }

    [Theory]
    [InlineData(0.3, 0.1, 0.0, 1.0, 0.5)]       // wide, both sides of the mean
    [InlineData(0.5, 0.3, 0.4, 0.6, 0.5)]       // narrow truncation around the mean
    [InlineData(0.7, 0.2, 0.31, 0.3115, 0.31075)] // knots 0.0015 apart: the recurrence's failure case
    [InlineData(0.0, 1.0, 3.0, 5.0, 4.0)]       // tail
    [InlineData(0.0, 1.0, 6.0, 6.02, 6.01)]     // narrow and in the tail
    [InlineData(5.0, 0.1, 0.0, 1.0, 0.5)]       // 40-50 sds away
    public void Normal_Series_MatchesCompositeQuadrature(double mean, double sd, double a, double b, double centre)
    {
        Normal f = Normal.MeanSd(mean, sd);
        var m = new double[Powers];
        // Scaled by the interval's log mass so the far cases compare on a unit scale.
        double logMass = f.LogIntervalProbability(a, b);
        Assert.Equal(0.0, f.IntervalMoments(a, b, centre, logMass, m));
        double[] reference = Composite(x => Math.Exp(f.LogPdf(x) - logMass), a, b, centre, panels: 400);
        // Far out, the log density −z²/2 itself carries a rounding of z²/2 ulps, which both sides inherit.
        double near = Math.Min(Math.Abs(a - mean), Math.Abs(b - mean)) / sd, tolerance = 1e-13 + 4e-16 * near * near;
        for (int k = 0; k < Powers; k++)
            Assert.True(Math.Abs(m[k] - reference[k]) <= tolerance * Math.Max(Math.Abs(reference[k]), Math.Pow(b - a, k)),
                $"k = {k}: {m[k]} against {reference[k]}");
        Assert.Equal(1.0, m[0], 12);
    }

    [Fact]
    public void Normal_TruncatedFarInTheTail_KeepsItsDigits()
    {
        // N(0, 1) on [40, 41]: the unscaled moments underflow, and the conditional mean is 40 + 1/40 − 2/40³ + 10/40⁵ − ….
        var t = Truncated.Of(Normal.Standard, 40.0, 41.0);
        var m = new double[3];
        t.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 40.0, 0.0, m);
        Assert.Equal(1.0, m[0], 12);
        double mills = 1.0 / 40 - 2.0 / Math.Pow(40, 3) + 10.0 / Math.Pow(40, 5);   // E[z − 40 | z > 40], [41, ∞) negligible
        Assert.Equal(mills, m[1], 8);
        Assert.Equal(0.0, Normal.Standard.IntervalProbability(40.0, 41.0));
    }

    [Fact]
    public void Cauchy_MatchesItsClosedForms()
    {
        // Standard Cauchy on (a, b]: ∫ 1 = Δatan/π, ∫ t = ln((1 + b²)/(1 + a²))/(2π), ∫ t² = (b − a − Δatan)/π.
        StudentT f = StudentT.Standard(1.0);
        foreach ((double a, double b) in new[] { (-3.0, 4.0), (0.2, 0.9), (-50.0, -2.0), (10.0, 1000.0) })
        {
            var m = new double[3];
            double error = f.IntervalMoments(a, b, 0.0, 0.0, m);
            double atan = Math.Atan(b) - Math.Atan(a);
            Assert.Equal(atan / Math.PI, m[0], 14);
            Assert.Equal(Math.Log((1 + b * b) / (1 + a * a)) / (2 * Math.PI), m[1], 13);
            Assert.True(Math.Abs(m[2] - (b - a - atan) / Math.PI) <= 1e-13 * Math.Max(1.0, b - a), $"({a}, {b}): {m[2]}");
            Assert.InRange(error, 0.0, 1e-12 * Math.Max(1.0, b * b));
        }
    }

    [Theory]
    [InlineData(2.0, 0.5, 0.15, 0.0, 1.0, 0.5)]
    [InlineData(3.0, 0.5, 0.15, 0.0, 0.4, 0.2)]
    [InlineData(3.0, 0.3, 0.1, 0.62, 0.63, 0.625)]
    [InlineData(1.0, 0.7, 0.2, 0.0, 1.0, 0.5)]
    [InlineData(30.0, 0.5, 0.01, 0.0, 1.0, 0.5)]
    public void StudentT_MatchesAFineCompositeRule(double nu, double location, double scale, double a, double b, double centre)
    {
        StudentT f = StudentT.LocationScale(nu, location, scale);
        var m = new double[Powers];
        double error = f.IntervalMoments(a, b, centre, 0.0, m);
        double[] reference = Composite(x => Math.Exp(f.LogPdf(x)), a, b, centre, panels: 4000);
        for (int k = 0; k < Powers; k++)
            Assert.True(Math.Abs(m[k] - reference[k]) <= 1e-13 * Math.Max(Math.Abs(reference[k]), Math.Pow(b - a, k) * m[0]),
                $"k = {k}: {m[k]} against {reference[k]}");
        Assert.InRange(error, 0.0, 1e-13);
        Assert.Equal(f.IntervalProbability(a, b), m[0], 13);
    }

    [Fact]
    public void StudentT_RefusesAnInfiniteInterval()
    {
        var m = new double[2];
        Assert.Throws<ArgumentException>(() => StudentT.Standard(3.0).IntervalMoments(0.0, double.PositiveInfinity, 0.0, 0.0, m));
        // Truncated to a finite interval, it has them.
        var t = Truncated.Of(StudentT.Standard(1.0), -1.0, 2.0);
        t.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.0, 0.0, m);
        Assert.Equal(1.0, m[0], 14);
        Assert.Equal(Math.Log(5.0 / 2.0) / (2 * Math.PI) / t.Mass, m[1], 13);
    }

    [Fact]
    public void Mixture_ConventionsDiffer_AndEachMatchesItsDefinition()
    {
        // weights (0.4, 0.6), N(0.3, 0.1²) and N(0.7, 0.2²), on [0, 1].
        Normal f1 = Normal.MeanSd(0.3, 0.1), f2 = Normal.MeanSd(0.7, 0.2);
        Mixture mix = Mixture.Of([f1, f2], [2.0, 3.0]);
        Assert.Equal(0.4, mix.Weights[0], 15);
        Assert.Equal(Math.Log(0.4 * Math.Exp(f1.LogPdf(0.5)) + 0.6 * Math.Exp(f2.LogPdf(0.5))), mix.LogPdf(0.5), 14);
        Assert.Equal(0.4 * f1.IntervalProbability(0.0, 1.0) + 0.6 * f2.IntervalProbability(0.0, 1.0), mix.IntervalProbability(0.0, 1.0), 14);

        double z1 = f1.IntervalProbability(0.0, 1.0), z2 = f2.IntervalProbability(0.0, 1.0);
        var before = Truncated.Of(mix, 0.0, 1.0);
        Mixture conditional = Mixture.Of([Truncated.Of(f1, 0.0, 1.0), Truncated.Of(f2, 0.0, 1.0)], [0.4, 0.6]);
        var m1 = new double[Powers];
        var m2 = new double[Powers];
        var mb = new double[Powers];
        var mc = new double[Powers];
        f1.IntervalMoments(0.0, 1.0, 0.5, 0.0, m1);
        f2.IntervalMoments(0.0, 1.0, 0.5, 0.0, m2);
        before.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.5, 0.0, mb);
        conditional.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.5, 0.0, mc);
        for (int k = 0; k < Powers; k++)
        {
            Assert.Equal((0.4 * m1[k] + 0.6 * m2[k]) / (0.4 * z1 + 0.6 * z2), mb[k], 14);
            Assert.Equal(0.4 * m1[k] / z1 + 0.6 * m2[k] / z2, mc[k], 14);
        }
        Assert.True(Math.Abs(mb[1] - mc[1]) > 1e-3);
        Assert.Throws<ArgumentException>(() => Mixture.Of([f1], [0.0]));
        Assert.Throws<InvalidOperationException>(() => default(Mixture).LogPdf(0.0));
    }

    [Fact]
    public void Empirical_SumsItsAtoms_CountingEachOnceAcrossAPartition()
    {
        double[] sample = [0.0, 0.2, 0.2, 0.5, 1.0];
        Empirical e = Empirical.Of(sample);
        Assert.Equal(0.4, e.IntervalProbability(0.1, 0.2), 15);     // (0.1, 0.2] holds both 0.2s
        var whole = new double[Powers];
        var left = new double[Powers];
        var right = new double[Powers];
        e.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, 0.3, 0.0, whole);
        e.IntervalMoments(double.NegativeInfinity, 0.2, 0.3, 0.0, left);      // splits at an atom
        e.IntervalMoments(0.2, double.PositiveInfinity, 0.3, 0.0, right);
        for (int k = 0; k < Powers; k++)
        {
            double sum = 0.0;
            foreach (double x in sample) sum += Math.Pow(x - 0.3, k);
            Assert.Equal(sum / sample.Length, whole[k], 15);
            Assert.Equal(whole[k], left[k] + right[k], 15);
        }
        Assert.Equal(0.6, left[0], 15);
        Assert.Throws<ArgumentException>(() => Empirical.Of([]));
    }

    // Composite 20-node Gauss–Legendre of (x − c)ᵏ f(x) over [a, b] on equal panels.
    private static double[] Composite(Func<double, double> f, double a, double b, double centre, int panels)
    {
        var sums = new double[Powers];
        GaussLegendre rule = GaussLegendre.Of(20);
        double width = (b - a) / panels;
        for (int p = 0; p < panels; p++)
        {
            double mid = a + (p + 0.5) * width, half = 0.5 * width;
            for (int i = 0; i < rule.Count; i++)
            {
                double x = mid + half * rule.Nodes[i];
                double w = half * rule.Weights[i] * f(x);
                for (int k = 0; k < Powers; k++)
                {
                    sums[k] += w;
                    w *= x - centre;
                }
            }
        }
        return sums;
    }
}
