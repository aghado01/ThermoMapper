using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BAPS's admission after <see cref="Maths.LinAlg.BandCholesky"/> became strict by default (shared-code defects brief,
/// slice 5), on the nine tensor layouts of <c>issues/linalg/spectral-tensor-factorization-report.md</c> and a 1-D penalized
/// spline. A valid fit in small units, which the floored factor refused, is admitted and reads as the same fit. A design that
/// leaves the penalty's null space unidentified is refused at every smoothing of the optimizer's box, where the floored gate
/// admitted some. The public λx = λy = 1e20 points are refused, and λ = 1e12 is still admitted.
/// </summary>
public sealed class BapsAdmissionTests
{
    public static TheoryData<int, int> Layouts => new() { { 4, 4 }, { 6, 8 }, { 6, 16 }, { 8, 8 }, { 8, 16 }, { 16, 8 }, { 16, 16 }, { 16, 32 }, { 32, 16 } };

    // The report's fixture: cubic bases, ν − 4 interior knots at k/(ν − 3), a midpoint grid of max(20, 2ν) per axis, and its
    // deterministic response. With fixY, every observation sits at y = 0.5, so the design leaves the penalty's y-linear part
    // unidentified.
    private static (TensorDesign Td, double[] Y) Tensor(int nuX, int nuY, bool fixY = false)
    {
        int gx = Math.Max(20, 2 * nuX), gy = Math.Max(20, 2 * nuY);
        var xs = new List<double>();
        var ys = new List<double>();
        var y = new List<double>();
        for (int a = 0; a < gx; a++)
            for (int b = 0; b < gy; b++)
            {
                double x = (a + 0.5) / gx, v = fixY ? 0.5 : (b + 0.5) / gy;
                xs.Add(x);
                ys.Add(v);
                y.Add(Math.Sin(1.5 * Math.PI * x) * Math.Sin(5 * Math.PI * v) + 0.1 * Math.Cos(37 * x + 53 * v) + 0.05 * Math.Sin(89 * x - 71 * v));
            }
        KnotConfig Knots(int nu) => new([.. Enumerable.Range(1, nu - 4).Select(k => k / (nu - 3.0))]);
        return (TensorProductDesign.Build(new SplineBasis(3), Knots(nuX), [.. xs], new SplineBasis(3), Knots(nuY), [.. ys]), [.. y]);
    }

    private static double[,] Scaled(double[,] z, double c)
    {
        var s = (double[,])z.Clone();
        for (int i = 0; i < s.GetLength(0); i++) for (int j = 0; j < s.GetLength(1); j++) s[i, j] *= c;
        return s;
    }

    // Each axis's trace balance λ_ref = tr(ZᵀZ)/tr(P_axis), and the optimizer's box λ_ref·10^{−8, −4, 0, 4, 8} per axis.
    private static IEnumerable<(double X, double Y)> Box(TensorDesign td)
    {
        double trace = 0.0;
        for (int i = 0; i < td.Design.GetLength(0); i++) for (int j = 0; j < td.Design.GetLength(1); j++) trace += td.Design[i, j] * td.Design[i, j];
        var penalty = new TensorPenalty(td.NuX, td.NuY, 2, 2);
        double Axis(int axis)
        {
            int nu = td.NuX * td.NuY;
            var band = new double[penalty.Bandwidth + 1, nu];
            penalty.AccumulateComponents(band, nu, axis == 0 ? [1.0, 0.0] : [0.0, 1.0]);
            double t = 0.0;
            for (int j = 0; j < nu; j++) t += band[0, j];
            return t;
        }
        double refX = trace / Axis(0), refY = trace / Axis(1);
        foreach (int a in new[] { -8, -4, 0, 4, 8 })
            foreach (int b in new[] { -8, -4, 0, 4, 8 })
                yield return (refX * Math.Pow(10, a), refY * Math.Pow(10, b));
    }

    // A valid fit in units 2^24 smaller: Z' = cZ and λ' = c²λ give β' = β/c exactly, which the floored factor refused because
    // every pivot of c²(ZᵀZ + λP) sits below its 1e-12.
    [Fact]
    public void AValidFitInSmallUnits_IsAdmitted_AndIsTheSameFit()
    {
        const int n = 300;
        double c = Math.ScaleB(1.0, -24);
        double[] x = [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)];
        double[] y = [.. x.Select(v => Math.Sin(2.0 * Math.PI * v) + 0.1 * Math.Cos(31 * v))];
        double[,] z = new SplineBasis(3).Design(new KnotConfig([.. Enumerable.Range(1, 16).Select(k => k / 17.0)]), x);
        var unit = new PenalizedSpline(z, y, new DifferencePenalty(2));
        var small = new PenalizedSpline(Scaled(z, c), y, new DifferencePenalty(2));
        foreach (double lambda in new[] { 1e-3, 1.0, 1e3 })
        {
            PenalizedSolve a = unit.SolveFactored(lambda), b = small.SolveFactored(lambda * c * c);
            for (int j = 0; j < a.Beta.Length; j++) Assert.Equal(BitConverter.DoubleToInt64Bits(a.Beta[j] / c), BitConverter.DoubleToInt64Bits(b.Beta[j]));
            Assert.Equal(a.RssPen, b.RssPen);
            Assert.Equal(a.Factor.LogDet + a.Beta.Length * Math.Log(c * c), b.Factor.LogDet, 9);
        }

        (TensorDesign td, double[] ty) = Tensor(8, 8);
        var tensor = new AnisotropicTensorSpline(td, ty);
        var smallTensor = new AnisotropicTensorSpline(td with { Design = Scaled(td.Design, c) }, ty);
        foreach ((double lx, double ly) in Box(td))
        {
            double[] a = tensor.Coefficients(lx, ly), b = smallTensor.Coefficients(lx * c * c, ly * c * c);
            for (int j = 0; j < a.Length; j++) Assert.Equal(BitConverter.DoubleToInt64Bits(a[j] / c), BitConverter.DoubleToInt64Bits(b[j]));
        }
    }

    [Theory]
    [MemberData(nameof(Layouts))]
    public void AnUnidentifiedDesign_IsRefusedAcrossTheBox(int nuX, int nuY)
    {
        (TensorDesign identified, _) = Tensor(nuX, nuY);
        (TensorDesign td, double[] y) = Tensor(nuX, nuY, fixY: true);
        var spline = new AnisotropicTensorSpline(td, y);
        foreach ((double lx, double ly) in Box(identified))
            Assert.Throws<InvalidOperationException>(() => spline.RemlLogEvidence(lx, ly));
    }

    [Theory]
    [MemberData(nameof(Layouts))]
    public void AnIdentifiedDesign_IsAdmittedAcrossTheBox_AndAtTheTwelfthPower_ButNotTheTwentieth(int nuX, int nuY)
    {
        (TensorDesign td, double[] y) = Tensor(nuX, nuY);
        var spline = new AnisotropicTensorSpline(td, y);
        foreach ((double lx, double ly) in Box(td)) Assert.True(double.IsFinite(spline.RemlLogEvidence(lx, ly)));
        Assert.True(double.IsFinite(spline.RemlLogEvidence(1e12, 1e12)));
        Assert.Throws<InvalidOperationException>(() => spline.RemlLogEvidence(1e20, 1e20));
    }
}
