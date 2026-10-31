using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The conjugate updates of a normal scale as limits of a GIG update (N-X3, <see cref="GigParameters"/>).
/// <list type="bullet">
/// <item>The GIG update itself against Bayes' rule, for a prior with a, b &gt; 0: prior × likelihood and the updated kernel
/// differ by a constant, for the variance form and the precision form.</item>
/// <item>The two limits are exactly the gamma and inverse-gamma updates: the same shapes and rates, bit for bit, and the
/// kernels at b = 0 and a = 0 differ from those densities by a constant.</item>
/// <item>The normalized GIG density, its constant by the integral form of Bessel K, converges to each update's density
/// as b, or a, goes to 0, at the rate of the small parameter.</item>
/// </list>
/// </summary>
public sealed class GigLimitTests
{
    // K_ν(z) = ∫₀^∞ exp(−z·cosh t)·cosh(νt) dt, by the trapezoid rule, which converges geometrically for this analytic
    // integrand; ln K by factoring out the largest exponent.
    private static double LogBesselK(double nu, double z)
    {
        double top = Math.Acosh(1.0 + 745.0 / z);   // z·cosh t − z exceeds 745: the integrand underflows
        const double step = 0.002;
        int count = (int)Math.Ceiling(top / step);
        double shift = double.NegativeInfinity;
        var logs = new double[count + 1];
        for (int k = 0; k <= count; k++)
        {
            double t = k * step;
            logs[k] = -z * Math.Cosh(t) + Math.Abs(nu) * t + Math.Log(0.5 * (1.0 + Math.Exp(-2.0 * Math.Abs(nu) * t)));
            shift = Math.Max(shift, logs[k]);
        }
        double sum = 0.0;
        for (int k = 0; k <= count; k++) sum += (k == 0 ? 0.5 : 1.0) * Math.Exp(logs[k] - shift);
        return shift + Math.Log(sum * step);
    }

    // ln f of the normalized GIG(p, a, b): (p/2)·ln(a/b) − ln(2·K_p(√(ab))) + the kernel.
    private static double LogGig(GigParameters g, double x)
        => 0.5 * g.P * Math.Log(g.A / g.B) - Math.Log(2.0) - LogBesselK(g.P, Math.Sqrt(g.A * g.B)) + g.LogKernel(x);

    [Fact]
    public void TheBesselIntegral_MatchesKnownValues()
    {
        // K_{1/2}(z) = √(π/(2z))·e^(−z); K_{3/2}(z) = √(π/(2z))·e^(−z)·(1 + 1/z).
        foreach (double z in new[] { 1e-3, 0.4, 1.0, 7.5 })
        {
            double half = 0.5 * Math.Log(Math.PI / (2.0 * z)) - z;
            Assert.Equal(half, LogBesselK(0.5, z), 9);
            Assert.Equal(half + Math.Log(1.0 + 1.0 / z), LogBesselK(1.5, z), 9);
            Assert.Equal(LogBesselK(2.3, z), LogBesselK(-2.3, z), 12);
        }
    }

    [Fact]
    public void TheGigUpdate_IsBayesRule_ForAGenuineGigPrior()
    {
        var rng = new Xoshiro256PlusPlus(11);
        var prior = new GigParameters(-0.7, 1.3, 2.1);
        double[] r = [.. Enumerable.Range(0, 6).Select(_ => 2.0 * rng.NextDouble() - 1.0)];
        double s = r.Sum(x => x * x);
        // Variance v: likelihood Πᵢ v^(−1/2)·exp(−rᵢ²/(2v)).
        GigParameters byVariance = ConjugateUpdates.NormalVariance(prior, r.Length, s);
        double[] offsets = [.. new[] { 0.2, 0.9, 1.7, 3.5 }.Select(v => prior.LogKernel(v) - 0.5 * r.Length * Math.Log(v) - 0.5 * s / v - byVariance.LogKernel(v))];
        Assert.All(offsets, o => Assert.Equal(offsets[0], o, 12));
        // Precision φ: likelihood Πᵢ φ^(1/2)·exp(−φ rᵢ²/2).
        GigParameters byPrecision = ConjugateUpdates.NormalPrecision(prior, r.Length, s);
        offsets = [.. new[] { 0.2, 0.9, 1.7, 3.5 }.Select(f => prior.LogKernel(f) + 0.5 * r.Length * Math.Log(f) - 0.5 * f * s - byPrecision.LogKernel(f))];
        Assert.All(offsets, o => Assert.Equal(offsets[0], o, 12));
    }

    [Fact]
    public void TheLimits_AreTheGammaAndInverseGammaUpdates_ExactlyAndBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(5);
        for (int trial = 0; trial < 200; trial++)
        {
            double shape = 0.1 + 5.0 * rng.NextDouble(), rate = 0.01 + 4.0 * rng.NextDouble();
            double count = 1 + (int)(30 * rng.NextDouble()), s = 20.0 * rng.NextDouble();

            Gamma gamma = ConjugateUpdates.NormalPrecision(shape, rate, count, s);
            Gamma viaGig = ConjugateUpdates.NormalPrecision(GigParameters.OfGamma(shape, rate), count, s).ToGamma();
            Assert.Equal(gamma.Shape, viaGig.Shape);
            Assert.Equal(gamma.Rate, viaGig.Rate);

            InverseGamma inverse = ConjugateUpdates.NormalVariance(shape, rate, count, s);
            InverseGamma viaGigInverse = ConjugateUpdates.NormalVariance(GigParameters.OfInverseGamma(shape, rate), count, s).ToInverseGamma();
            Assert.Equal(inverse.Shape, viaGigInverse.Shape);
            Assert.Equal(inverse.Scale, viaGigInverse.Scale);
        }
    }

    [Fact]
    public void TheKernels_AtTheLimits_DifferFromTheGammaAndInverseGammaDensitiesByAConstant()
    {
        GigParameters g = GigParameters.OfGamma(2.5, 1.4);
        double[] x = [0.1, 0.6, 1.0, 2.2, 7.0];
        double[] offsets = [.. x.Select(v => g.LogKernel(v) - g.ToGamma().LogPdf(v))];
        Assert.All(offsets, o => Assert.Equal(offsets[0], o, 12));

        GigParameters h = GigParameters.OfInverseGamma(3.1, 0.8);
        offsets = [.. x.Select(v => h.LogKernel(v) - h.ToInverseGamma().LogPdf(v))];
        Assert.All(offsets, o => Assert.Equal(offsets[0], o, 12));

        Assert.True(g.IsGamma && !g.IsInverseGamma);
        Assert.True(h.IsInverseGamma && !h.IsGamma);
        Assert.Throws<InvalidOperationException>(() => new GigParameters(1.0, 1.0, 1.0).ToGamma());
        Assert.Throws<InvalidOperationException>(() => new GigParameters(-1.0, 1.0, 1.0).ToInverseGamma());
    }

    [Fact]
    public void TheNormalizedGigDensity_ConvergesToTheGammaDensity_AsBGoesToZero()
    {
        Gamma gamma = ConjugateUpdates.NormalPrecision(1.5, 0.7, 12, 9.0);
        double[] x = [0.2, 0.8, 1.5, 3.0];
        double previous = double.PositiveInfinity;
        foreach (double b in new[] { 1e-2, 1e-4, 1e-6, 1e-8 })
        {
            GigParameters g = new(gamma.Shape, 2.0 * gamma.Rate, b);
            double error = x.Max(v => Math.Abs(LogGig(g, v) - gamma.LogPdf(v)));
            Assert.True(error < previous, $"b = {b}: {error} after {previous}");
            Assert.True(error <= 4.0 * b, $"b = {b}: {error}");   // O(b): the b/(2x) term and the constant's O(ab) correction
            previous = error;
        }
        Assert.True(previous < 1e-7);
    }

    [Fact]
    public void TheNormalizedGigDensity_ConvergesToTheInverseGammaDensity_AsAGoesToZero()
    {
        InverseGamma inverse = ConjugateUpdates.NormalVariance(1.5, 0.7, 12, 9.0);
        double[] x = [0.2, 0.8, 1.5, 3.0];
        double previous = double.PositiveInfinity;
        foreach (double a in new[] { 1e-2, 1e-4, 1e-6, 1e-8 })
        {
            GigParameters g = new(-inverse.Shape, a, 2.0 * inverse.Scale);
            double error = x.Max(v => Math.Abs(LogGig(g, v) - inverse.LogPdf(v)));
            Assert.True(error < previous, $"a = {a}: {error} after {previous}");
            Assert.True(error <= 4.0 * a * Math.Max(1.0, 3.0), $"a = {a}: {error}");   // a·x/2 at x up to 3
            previous = error;
        }
        Assert.True(previous < 1e-7);
    }
}
