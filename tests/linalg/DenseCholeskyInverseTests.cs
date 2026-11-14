using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// <see cref="DenseCholesky"/>'s reads for a Gaussian's covariance: the inverse formed on demand, bit for bit as the retired
/// <c>CholeskyDecomposition</c> formed it, and refused where it overflows; and <c>μ + L·z</c> accumulated in its order
/// (shared-code defects brief, slice 2's follow-up).
/// </summary>
public sealed class DenseCholeskyInverseTests
{
    // L bidiagonal, 2^-20 on its diagonal and 1 below it, so A = LLᵀ is positive definite and factors exactly in powers of
    // two; L⁻¹'s first column is (−1)^i·2^(20(i+1)), past double.MaxValue from row 51.
    private static double[,] OverflowingInverse(int n)
    {
        double d = Math.ScaleB(1.0, -20);
        var a = new double[n, n];
        a[0, 0] = d * d;
        for (int i = 1; i < n; i++)
        {
            a[i, i] = 1.0 + d * d;
            a[i, i - 1] = a[i - 1, i] = d;
        }
        return a;
    }

    // The retired class's L⁻¹ and inverse, verbatim, read from a factor.
    private static double[,] EarlierInverse(DenseCholesky f)
    {
        int n = f.Dimension;
        var l = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) l[i, j] = f.Factor(i, j);
        var lInv = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            lInv[i, i] = 1.0 / l[i, i];
            for (int j = 0; j < i; j++)
            {
                double sum = 0.0;
                for (int k = j; k < i; k++) sum -= l[i, k] * lInv[k, j];
                lInv[i, j] = sum / l[i, i];
                lInv[j, i] = 0.0;
            }
        }
        var target = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
            {
                double sum = 0.0;
                for (int k = j; k < n; k++) sum += lInv[k, i] * lInv[k, j];
                target[i, j] = sum;
                target[j, i] = sum;
            }
        return target;
    }

    [Fact]
    [Trait("Tier", "Pinned")]
    public void TheInverse_IsTheRetiredClassesBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1113);
        for (int t = 0; t < 60; t++)
        {
            int n = 1 + (int)(rng.NextUInt64() % 40);
            double[,] a = DenseSpd(n, rng);
            double scale = Math.ScaleB(1.0, (int)(rng.NextUInt64() % 120) - 80);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i, j] *= scale;
            var f = new DenseCholesky(n);
            f.Decompose(a);
            var inverse = new double[n, n];
            Assert.True(f.TryWriteInverseTo(inverse));
            double[,] earlier = EarlierInverse(f);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    Assert.Equal(BitConverter.DoubleToInt64Bits(earlier[i, j]), BitConverter.DoubleToInt64Bits(inverse[i, j]));
        }
    }

    [Fact]
    public void TheInverse_InvertsTheMatrix()
    {
        const int n = 12;
        double[,] a = DenseSpd(n, new Xoshiro256PlusPlus(seed: 7));
        var f = new DenseCholesky(n);
        f.Decompose(a);
        var inverse = new double[n, n];
        Assert.True(f.TryWriteInverseTo(inverse));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < n; k++) s += a[i, k] * inverse[k, j];
                Assert.Equal(i == j ? 1.0 : 0.0, s, 10);
            }
    }

    [Fact]
    public void AnInverseThatOverflows_IsRefused_ThoughTheFactorHolds()
    {
        // [1e-320]: a finite factor and log-determinant, a reciprocal past double.MaxValue.
        var f = new DenseCholesky(1);
        f.Decompose(new double[,] { { 1e-320 } });
        Assert.Equal(-1, f.FailedAt);
        Assert.True(double.IsFinite(f.LogDet));
        Assert.False(f.TryWriteInverseTo(new double[1, 1]));

        // L bidiagonal, 2^-20 on its diagonal and 1 below it: L⁻¹'s first column passes double.MaxValue from row 51.
        const int n = 60;
        var g = new DenseCholesky(n);
        g.Decompose(OverflowingInverse(n));
        Assert.Equal(-1, g.FailedAt);
        Assert.False(g.TryWriteInverseTo(new double[n, n]));
        Assert.Throws<ArgumentException>(() => g.TryWriteInverseTo(new double[n - 1, n]));
    }

    [Fact]
    public void AddSample_AccumulatesLzIntoTheMeanInColumnOrder()
    {
        const int n = 9;
        var rng = new Xoshiro256PlusPlus(seed: 23);
        double[,] a = DenseSpd(n, rng);
        var f = new DenseCholesky(n);
        f.Decompose(a);
        double[] z = Vector(n, rng), mean = Vector(n, rng);
        double[] x = (double[])mean.Clone();
        f.AddSample(z, x);
        for (int i = 0; i < n; i++)
        {
            double expected = mean[i];
            for (int j = 0; j <= i; j++) expected += f.Factor(i, j) * z[j];
            Assert.Equal(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(x[i]));
        }
        Assert.Throws<ArgumentException>(() => f.AddSample(z, new double[n - 1]));
        Assert.Throws<InvalidOperationException>(() => new DenseCholesky(2).AddSample(new double[2], new double[2]));
    }
}
