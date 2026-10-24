using System;
using Maths.Regression.Changepoint;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The Cauchy segment model's quadrature against composite Simpson on a wide, fine window per segment: the automatic
/// grid resolves the narrow posteriors of long segments under light noise, and keeps the mass a narrow prior puts
/// beyond ν ± 6ρ. A fixed coarse grid reports its coarseness through <see cref="CauchySegmentModel.Resolution"/>.
/// </summary>
public sealed class CauchySegmentModelTests
{
    private static double[] Steps(int seed, int[] lengths, double[] levels, double noise)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int n = 0;
        foreach (int l in lengths) n += l;
        var y = new double[n];
        int i = 0;
        for (int s = 0; s < lengths.Length; s++)
            for (int j = 0; j < lengths[s]; j++, i++)
                y[i] = levels[s] + noise * Math.Sqrt(-2.0 * Math.Log(1.0 - rng.NextDouble())) * Math.Cos(2.0 * Math.PI * rng.NextDouble());
        return y;
    }

    // log ∫ N(μ; ν, ρ²) Π_{i=a..b} Cauchy(y_i | μ, γ) dμ and the level's mean and variance, by composite Simpson on
    // [min(ν, segment) − 10ρ, max(ν, segment) + 10ρ] at a fortieth of the narrowest posterior's width.
    private static (double LogEvidence, double Mean, double Var) Simpson(double[] y, int a, int b, double scale, double nu, double rho2)
    {
        double rho = Math.Sqrt(rho2), lo = nu, hi = nu;
        for (int i = a; i <= b; i++) { lo = Math.Min(lo, y[i - 1]); hi = Math.Max(hi, y[i - 1]); }
        lo -= 10.0 * rho;
        hi += 10.0 * rho;
        double width = 1.0 / Math.Sqrt(1.0 / rho2 + 2.0 * (b - a + 1) / (scale * scale));
        int intervals = 2 * (int)Math.Ceiling((hi - lo) / (width / 40.0) / 2.0);
        double h = (hi - lo) / intervals;
        var l = new double[intervals + 1];
        double top = double.NegativeInfinity;
        for (int k = 0; k <= intervals; k++)
        {
            double m = lo + k * h, v = -0.5 * Math.Log(2.0 * Math.PI * rho2) - 0.5 * (m - nu) * (m - nu) / rho2;
            for (int i = a; i <= b; i++)
            {
                double z = (y[i - 1] - m) / scale;
                v -= Math.Log(Math.PI * scale) + Math.Log(1.0 + z * z);
            }
            l[k] = v;
            top = Math.Max(top, v);
        }
        double s0 = 0.0, s1 = 0.0, s2 = 0.0;
        for (int k = 0; k <= intervals; k++)
        {
            double w = (k == 0 || k == intervals ? 1.0 : k % 2 == 1 ? 4.0 : 2.0) * Math.Exp(l[k] - top), d = lo + k * h - nu;
            s0 += w; s1 += w * d; s2 += w * d * d;
        }
        double offset = s1 / s0;
        return (top + Math.Log(s0 * h / 3.0), nu + offset, s2 / s0 - offset * offset);
    }

    private static (double[,] LogD, double[,] Mean, double[,] Var) Tables(ISegmentModel model)
    {
        int n = model.Count;
        var d = new double[n + 1, n + 1];
        var m = new double[n + 1, n + 1];
        var v = new double[n + 1, n + 1];
        model.Fill(d, m, v);
        return (d, m, v);
    }

    [Fact]
    public void AutomaticGrid_ResolvesNarrowPosteriors_UnderLightNoise()
    {
        // Integer steps with noise 0.01 and a matching Cauchy width: ρ/γ ≈ 100, where a long segment's posterior is about
        // γ/√(2L) wide, far below the former fixed grid's spacing of 12ρ/511.
        const double noise = 0.01;
        double[] y = Steps(7, new[] { 10, 10, 10, 10 }, new[] { 0.0, 1.0, 3.0, 1.0 }, noise);
        (_, double nu, double rho2) = ExactChangepoint.DefaultHyperparameters(y);
        var model = new CauchySegmentModel(y, noise, nu, rho2);
        Assert.True(model.Resolution <= 0.5 + 1e-12);
        (double[,] logD, double[,] mean, double[,] var) = Tables(model);

        foreach ((int a, int b) in new[] { (1, 10), (11, 20), (1, 40), (5, 5), (8, 33) })
        {
            (double refD, double refMean, double refVar) = Simpson(y, a, b, noise, nu, rho2);
            Assert.Equal(refD, logD[a, b], 9);
            Assert.True(Math.Abs(mean[a, b] - refMean) <= 1e-7 * Math.Sqrt(refVar), $"[{a},{b}] mean {mean[a, b]} vs {refMean}");
            Assert.Equal(1.0, Math.Sqrt(var[a, b] / refVar), 6);
        }

        // The former node count on the wider range is far too coarse here, and says so.
        var coarse = new CauchySegmentModel(y, noise, nu, rho2, gridPoints: 512);
        Assert.True(coarse.Resolution > 2.0, $"resolution {coarse.Resolution}");
        (double[,] coarseD, _, _) = Tables(coarse);
        Assert.True(Math.Abs(coarseD[1, 40] - Simpson(y, 1, 40, noise, nu, rho2).LogEvidence) > 0.1);
    }

    [Fact]
    public void Range_KeepsTheMassANarrowPriorPutsBeyondSixRho()
    {
        // A prior N(0, 0.5²) and a level near 5: the posterior sits beyond ν + 6ρ = 3, where the former grid ended.
        double[] y = Steps(13, new[] { 20 }, new[] { 5.0 }, 0.1);
        (double[,] logD, double[,] mean, _) = Tables(new CauchySegmentModel(y, 0.1, 0.0, 0.25));
        (double refD, double refMean, _) = Simpson(y, 1, 20, 0.1, 0.0, 0.25);
        Assert.True(mean[1, 20] > 4.9, $"level {mean[1, 20]}");
        Assert.Equal(refMean, mean[1, 20], 8);
        Assert.Equal(refD, logD[1, 20], 9);
    }

    [Fact]
    public void TooSmallAScale_IsRefused_RatherThanLeftUnresolved()
    {
        double[] y = Steps(3, new[] { 50, 50 }, new[] { 0.0, 10.0 }, 0.1);
        Assert.Throws<ArgumentOutOfRangeException>(() => new CauchySegmentModel(y, 1e-7, 5.0, 25.0));
    }
}
