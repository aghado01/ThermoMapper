using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// Every family against its own log tails, at parameters and intervals drawn over many orders of magnitude, where the other
/// tests choose theirs. A truncated draw reads no tail, so Kolmogorov–Smirnov against the tails is independent of it, and a
/// quantile or a density that disagrees with them is wrong on one side or the other. The seeds are fixed, so a failure repeats.
/// Rounding makes atoms the tests allow for: an x within a few ulps of 1, a value below the double range.
/// </summary>
public sealed class DistributionSweepTests
{
    private const int DrawCases = 600;
    private const int Draws = 2_000;
    private const int PointCases = 3_000;

    // P(√n·D > 2.6) ≈ 2e^(−2·2.6²), about 3 in a million per case.
    private const double KolmogorovCut = 2.6;
    private const double TailTolerance = 1e-6;
    private const double DensityTolerance = 1e-3;
    private const double Slice = 1e-3;

    [Fact]
    public void TruncatedDraws_FollowTheTails()
    {
        var rng = new Xoshiro256PlusPlus(20261102);
        var failures = new List<string>();
        var xs = new double[Draws];
        for (int c = 0; c < DrawCases; c++)
        {
            (string name, IContinuousDistribution d) = Family(rng);
            if (Interval(rng, d) is not { } box) continue;
            double logMass = d.LogIntervalProbability(box.Lower, box.Upper);
            if (!double.IsFinite(logMass)) continue;
            for (int i = 0; i < Draws; i++)
            {
                xs[i] = d.SampleOn(rng, box.Lower, box.Upper);
                Assert.InRange(xs[i], box.Lower, box.Upper);
            }
            Array.Sort(xs);
            double gap = 0.0;
            int first = 0;
            for (int i = 0; i < Draws; i++)
            {
                if (i + 1 < Draws && xs[i + 1] == xs[i]) continue;   // at each distinct value's last rank: rounding makes atoms
                double f = xs[i] <= box.Lower ? 0.0 : Math.Exp(d.LogIntervalProbability(box.Lower, xs[i]) - logMass);
                gap = Math.Max(gap, Math.Abs(f - (i + 1.0) / Draws));
                if (i == first) gap = Math.Max(gap, Math.Abs(f - (double)i / Draws));
                first = i + 1;
            }
            double statistic = Math.Sqrt(Draws) * gap;
            if (statistic > KolmogorovCut) failures.Add($"√n·D = {statistic:F2}: {name} on [{box.Lower:R}, {box.Upper:R}]");
        }
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    [Fact]
    public void Quantiles_InvertTheTails()
    {
        var rng = new Xoshiro256PlusPlus(20261103);
        var failures = new List<string>();
        void Check(string name, IContinuousDistribution d, Func<double, double> invert, Func<double, double> tail, double logP)
        {
            double x = invert(logP);
            if (!double.IsFinite(x) || x == 0.0) return;          // the root is beyond the double range
            if (d is Beta && 1.0 - x < 1e-10) return;             // too few doubles there; the complement holds the digits
            double gap = Math.Abs(tail(x) - logP) / Math.Max(1.0, Math.Abs(logP));
            if (!(gap <= TailTolerance)) failures.Add($"gap {gap:G3}: {name} at ln p = {logP:R}, x = {x:R}");
        }
        for (int c = 0; c < PointCases; c++)
        {
            (string name, IContinuousDistribution d) = Family(rng);
            double logP = -LogUniform(rng, 1e-6, 1e4);
            Check(name + " lower", d, d.InverseLogCdf, d.LogCdf, logP);
            Check(name + " upper", d, d.InverseLogSurvival, d.LogSurvival, logP);
        }
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    // A thin slice of the smaller tail: its mass from the tail, against Simpson's rule on the density over it.
    [Fact]
    public void Densities_AreTheTailsSlopes()
    {
        var rng = new Xoshiro256PlusPlus(20261104);
        var failures = new List<string>();
        for (int c = 0; c < PointCases; c++)
        {
            (string name, IContinuousDistribution d) = Family(rng);
            double logTail = -LogUniform(rng, 1.0, rng.NextDouble() < 0.5 ? 600.0 : 20.0);
            bool upper = rng.NextDouble() < 0.5;
            double a = upper ? d.InverseLogSurvival(logTail) : d.InverseLogCdf(logTail - Slice);
            double b = upper ? d.InverseLogSurvival(logTail - Slice) : d.InverseLogCdf(logTail);
            if (!double.IsFinite(a) || !double.IsFinite(b) || !(b > a)) continue;
            if (b - a < 1e-9 * Math.Max(Math.Abs(a), Math.Abs(b))) continue;   // too few doubles between them
            double pa = d.LogPdf(a), pm = d.LogPdf(0.5 * a + 0.5 * b), pb = d.LogPdf(b);
            double high = Math.Max(pa, Math.Max(pm, pb));
            double simpson = Math.Log(b - a) + high + Math.Log((Math.Exp(pa - high) + 4.0 * Math.Exp(pm - high) + Math.Exp(pb - high)) / 6.0);
            double gap = Math.Abs(simpson - (logTail + Math.Log(1.0 - Math.Exp(-Slice))));
            if (!(gap <= DensityTolerance)) failures.Add($"gap {gap:G3}: {name} on [{a:R}, {b:R}], ln f = [{pa:R}, {pm:R}, {pb:R}]");
        }
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    private static double LogUniform(Xoshiro256PlusPlus rng, double lo, double hi) => Math.Exp(Math.Log(lo) + rng.NextDouble() * (Math.Log(hi) - Math.Log(lo)));

    private static double Signed(Xoshiro256PlusPlus rng, double scale) => (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * LogUniform(rng, 1e-3, scale);

    private static (string Name, IContinuousDistribution Distribution) Family(Xoshiro256PlusPlus rng)
    {
        switch ((int)(rng.NextDouble() * 8.0))
        {
            case 0:
            {
                double mean = Signed(rng, 100.0), sd = LogUniform(rng, 1e-3, 1e3);
                return ($"Normal({mean:R}, {sd:R})", Normal.MeanSd(mean, sd));
            }
            case 1:
            {
                double shape = LogUniform(rng, 0.02, 1e6), rate = LogUniform(rng, 1e-3, 1e3);
                return ($"Gamma({shape:R}, rate {rate:R})", Gamma.ShapeRate(shape, rate));
            }
            case 2:
            {
                double shape = LogUniform(rng, 0.02, 1e6), scale = LogUniform(rng, 1e-3, 1e3);
                return ($"InverseGamma({shape:R}, {scale:R})", InverseGamma.ShapeScale(shape, scale));
            }
            case 3:
            {
                double a = LogUniform(rng, 0.02, 1e5), b = LogUniform(rng, 0.02, 1e5);
                return ($"Beta({a:R}, {b:R})", Beta.Shapes(a, b));
            }
            case 4:
            {
                double nu = LogUniform(rng, 0.1, 1e6), location = Signed(rng, 100.0), scale = LogUniform(rng, 1e-3, 1e3);
                return ($"StudentT({nu:R}, {location:R}, {scale:R})", StudentT.LocationScale(nu, location, scale));
            }
            case 5:
            {
                double shape = LogUniform(rng, 0.1, 20.0), scale = LogUniform(rng, 1e-3, 1e3);
                return ($"Weibull({shape:R}, {scale:R})", Weibull.ShapeScale(shape, scale));
            }
            case 6:
            {
                double scale = LogUniform(rng, 1e-3, 1e3);
                return ($"HalfCauchy({scale:R})", HalfCauchy.WithScale(scale));
            }
            default:
            {
                double m1 = Signed(rng, 10.0), s1 = LogUniform(rng, 1e-2, 10.0), m2 = Signed(rng, 10.0), s2 = LogUniform(rng, 1e-2, 10.0);
                double nu = LogUniform(rng, 0.5, 50.0), w = 0.05 + 0.9 * rng.NextDouble();
                return ($"Mixture({w:R}: N({m1:R}, {s1:R}), t({nu:R}, {m2:R}, {s2:R}))",
                        Mixture.Of([Normal.MeanSd(m1, s1), StudentT.LocationScale(nu, m2, s2)], [w, 1.0 - w]));
            }
        }
    }

    // A half-line, both sides of the bulk, or a narrow interval out in a tail, placed by the distribution's own quantiles.
    private static (double Lower, double Upper)? Interval(Xoshiro256PlusPlus rng, IContinuousDistribution d)
    {
        double Tail() => -LogUniform(rng, 1e-3, rng.NextDouble() < 0.3 ? 600.0 : 20.0);
        double lower, upper;
        switch ((int)(rng.NextDouble() * 5.0))
        {
            case 0:
                lower = d.InverseLogSurvival(Tail());
                upper = double.PositiveInfinity;
                break;
            case 1:
                lower = double.NegativeInfinity;
                upper = d.InverseLogCdf(Tail());
                break;
            case 2:
                lower = d.InverseLogCdf(Tail());
                upper = d.InverseLogSurvival(Tail());
                break;
            case 3:
            {
                double t = Tail();
                lower = d.InverseLogSurvival(t);
                upper = d.InverseLogSurvival(t - LogUniform(rng, 1e-2, 3.0));
                break;
            }
            default:
            {
                double t = Tail();
                upper = d.InverseLogCdf(t);
                lower = d.InverseLogCdf(t - LogUniform(rng, 1e-2, 3.0));
                break;
            }
        }
        if (double.IsNaN(lower) || double.IsNaN(upper) || !(upper > lower)) return null;
        if (upper < 1e-290) return null;   // draws pile up at 0, below the double range
        double largest = Math.Max(double.IsFinite(lower) ? Math.Abs(lower) : 0.0, double.IsFinite(upper) ? Math.Abs(upper) : 0.0);
        if (double.IsFinite(lower) && double.IsFinite(upper) && upper - lower < 1e-9 * largest) return null;   // too few doubles to test against
        return (lower, upper);
    }
}
