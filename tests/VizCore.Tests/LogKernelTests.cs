using System;
using Graphs.Coupling;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The kernels' log-space evaluators: ln J agrees with ln of the linear value where that is representable,
/// and stays exact in the tail where the linear value has underflowed to 0. The oracles are the kernels' closed
/// forms, not the evaluators themselves.
/// </summary>
public sealed class LogKernelTests
{
    private static readonly double[] Distances = { 0.0, 0.3, 1.0, 2.5, 5.0, 12.0 };
    private static readonly double[] Deltas = { 0.5, 1.0, 3.0 };

    private static void AgreesWhereRepresentable(Func<double, double, double> linear, Func<double, double, double> log)
    {
        foreach (double delta in Deltas)
            foreach (double d in Distances)
            {
                double j = linear(d, delta);
                if (j <= 1e-300) continue;
                Assert.InRange(log(d, delta), Math.Log(j) - 1e-12 * Math.Max(1.0, Math.Abs(Math.Log(j))),
                                              Math.Log(j) + 1e-12 * Math.Max(1.0, Math.Abs(Math.Log(j))));
            }
    }

    [Fact]
    public void Gaussian_Cauchy_Laplacian_AgreeWithTheLogOfTheLinearValue()
    {
        AgreesWhereRepresentable(GaussianKernel.Evaluate, GaussianKernel.LogEvaluate);
        AgreesWhereRepresentable(CauchyKernel.Evaluate, CauchyKernel.LogEvaluate);
        AgreesWhereRepresentable(LaplacianKernel.Evaluate, LaplacianKernel.LogEvaluate);
    }

    [Fact]
    public void Linear_AgreesInsideItsSupport_AndIsMinusInfinityFromOneScaleOut()
    {
        AgreesWhereRepresentable(LinearKernel.Evaluate, LinearKernel.LogEvaluate);
        Assert.Equal(0.0, LinearKernel.LogEvaluate(0.0, 2.0));
        Assert.Equal(Math.Log(0.5), LinearKernel.LogEvaluate(1.0, 2.0), 12);
        Assert.Equal(double.NegativeInfinity, LinearKernel.LogEvaluate(2.0, 2.0));
        Assert.Equal(double.NegativeInfinity, LinearKernel.LogEvaluate(9.0, 2.0));
    }

    [Fact]
    public void TheTails_StayExactWhereTheLinearValueIsZero()
    {
        // Gaussian at d/δ = 60: e^-1800 underflows, and ln J = −d²/2δ² = −1800 exactly.
        Assert.Equal(0.0, GaussianKernel.Evaluate(60.0, 1.0));
        Assert.Equal(-1800.0, GaussianKernel.LogEvaluate(60.0, 1.0));

        // Laplacian at d/δ = 1000.
        Assert.Equal(0.0, LaplacianKernel.Evaluate(1000.0, 1.0));
        Assert.Equal(-1000.0, LaplacianKernel.LogEvaluate(1000.0, 1.0));

        // Cauchy past the point where d² overflows: ln J = −2 ln(d/δ).
        Assert.Equal(0.0, CauchyKernel.Evaluate(1e160, 1.0));
        Assert.Equal(-2.0 * Math.Log(1e160), CauchyKernel.LogEvaluate(1e160, 1.0), 9);
    }

    [Fact]
    public void Cauchy_KeepsItsDigitsForASmallRatio()
    {
        // ln(1 + 1e-20) is 1e-20; a naive log(1 + x) returns 0.
        double ln = CauchyKernel.LogEvaluate(1e-10, 1.0);
        Assert.InRange(ln / -1e-20, 1.0 - 1e-9, 1.0 + 1e-9);
    }

    [Fact]
    public void LogJ_FallsAsTheDistanceGrows()
    {
        foreach (Func<double, double, double> log in new Func<double, double, double>[]
                 { GaussianKernel.LogEvaluate, CauchyKernel.LogEvaluate, LaplacianKernel.LogEvaluate })
        {
            double previous = double.PositiveInfinity;
            for (double d = 0.0; d <= 80.0; d += 4.0)
            {
                double now = log(d, 1.5);
                Assert.True(now <= previous);
                previous = now;
            }
        }
    }

    [Fact]
    public void Mixture_IsTheLogSumExpOfItsComponents_AndSurvivesEveryComponentUnderflowing()
    {
        var b = new MixtureBandwidth(Gaussian: 1.0, Cauchy: 1.5, Laplacian: 0.8);
        var w = new MixtureWeights(Gaussian: 0.2, Cauchy: 0.5, Laplacian: 0.3);

        foreach (double d in Distances)
            Assert.Equal(Math.Log(MixtureKernel.Evaluate(d, b, w)), MixtureKernel.LogEvaluate(d, b, w), 12);

        // Gaussian-only: the mixture's log is ln w + the Gaussian's, beyond the double range.
        Assert.Equal(-1800.0, MixtureKernel.LogEvaluate(60.0, b, new MixtureWeights(1.0, 0.0, 0.0)));

        // Gaussian and Laplacian, both underflowed at d = 2000: the slower decay dominates.
        var twin = new MixtureBandwidth(1.0, 1.0, 1.0);
        double expected = Math.Log(0.5) - 2000.0 + Math.Log(1.0 + Math.Exp(-(2000.0 * 2000.0 / 2.0) + 2000.0));
        Assert.Equal(0.0, MixtureKernel.Evaluate(2000.0, twin, new MixtureWeights(0.5, 0.0, 0.5)));
        Assert.Equal(expected, MixtureKernel.LogEvaluate(2000.0, twin, new MixtureWeights(0.5, 0.0, 0.5)), 9);

        Assert.Equal(double.NegativeInfinity, MixtureKernel.LogEvaluate(1.0, b, new MixtureWeights(0.0, 0.0, 0.0)));
    }
}
