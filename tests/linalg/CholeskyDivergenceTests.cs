using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// <see cref="CholeskyDecomposition"/> where its absolute pivot floor cannot carry the matrix. Two cases: a smooth
/// kernel's covariance, whose spectrum falls below rounding, so that floored pivots let the later rows overflow; and a
/// finite matrix whose pivot overflows while L stays finite, or whose L⁻¹ overflows. Each is refused as
/// <see cref="DenseCholesky"/> refuses a failed factor. Every finite factor, floored or not, reads bit for bit as the
/// class read it before (shared-code defects brief, slice 2).
/// </summary>
public sealed class CholeskyDivergenceTests
{
    // A squared-exponential covariance, variance 1, on n evenly spaced points of [0, 1].
    internal static double[,] SquaredExponential(int n, double length)
    {
        var k = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double r = (i - j) / (double)(n - 1) / length;
                k[i, j] = k[j, i] = Math.Exp(-0.5 * r * r);
            }
        return k;
    }

    // Rank one at 1e200, one entry 2ε above it: the second pivot is 0 and floored, so the third row's entry is the
    // rounding over 1e-6, 4.4e190, and its pivot, 1e200 − 1e200 − (4.4e190)², overflows to −∞. L stays finite.
    internal static double[,] OverflowingPivot()
    {
        double big = 1e200, off = big * (1 + 2 * Math.ScaleB(1.0, -52));
        return new[,] { { big, big, big }, { big, big, off }, { big, off, big } };
    }

    // Tridiagonal: 0 then 1 on the diagonal, s = √1e-12 beside it. Every pivot is exactly 0 and floored, so L is exactly
    // bidiagonal, s on its diagonal and 1 below it, and finite; L⁻¹'s first column is (−1)^i/s^(i+1), past
    // double.MaxValue from row 51.
    internal static double[,] OverflowingInverse(int n)
    {
        double s = Math.Sqrt(1e-12);
        var a = new double[n, n];
        for (int i = 1; i < n; i++)
        {
            a[i, i] = 1.0;
            a[i, i - 1] = a[i - 1, i] = s;
        }
        return a;
    }

    private static void AssertRefused(CholeskyDecomposition c, int n)
    {
        Assert.True(double.IsNaN(c.LogDet));
        Assert.Throws<InvalidOperationException>(() => c.Solve(new double[n]));
        Assert.Throws<InvalidOperationException>(() => c.WriteInverseTo(new double[n, n]));
        Assert.Throws<InvalidOperationException>(() => c.WriteLTo(new double[n * n]));
        Assert.Throws<InvalidOperationException>(() => c.Sample(new Random(1), new double[n]));
    }

    [Theory]
    [InlineData(50, 0.1)]
    [InlineData(100, 0.5)]
    [InlineData(200, 0.2)]
    [InlineData(401, 0.05)]
    public void ASmoothKernelsDivergingFactor_IsRefused(int n, double length)
    {
        var c = new CholeskyDecomposition(n);
        c.Decompose(SquaredExponential(n, length));
        Assert.False(c.Status.Finite);
        AssertRefused(c, n);
        Assert.InRange(c.FailedAt, 1, n - 1);
    }

    [Fact]
    public void AFinitePivotThatOverflows_IsRefused_ThoughLStaysFinite()
    {
        var c = new CholeskyDecomposition(3);
        c.Decompose(OverflowingPivot());
        Assert.False(c.Status.Finite);
        Assert.Equal(2, c.Status.FlooredPivots);
        AssertRefused(c, 3);
        Assert.Equal(2, c.FailedAt);
    }

    [Fact]
    public void AnInverseFactorThatOverflows_IsRefused_ThoughLAndThePivotsAreFinite()
    {
        const int n = 60;
        var c = new CholeskyDecomposition(n);
        c.Decompose(OverflowingInverse(n));
        Assert.True(c.Status.Finite);
        Assert.Equal(n, c.Status.FlooredPivots);
        AssertRefused(c, n);
        Assert.Equal(51, c.FailedAt);
    }

    [Fact]
    public void ANonFiniteEntry_IsRefused_AndAnUnfactoredInstanceToo()
    {
        var c = new CholeskyDecomposition(3);
        AssertRefused(c, 3);
        Assert.Equal(-1, c.FailedAt);
        Assert.Contains("No matrix", Assert.Throws<InvalidOperationException>(() => c.Solve(new double[3])).Message);
        c.Decompose(new double[,] { { 2, 0, 0 }, { 0, 2, 0 }, { 0, double.NaN, 2 } });
        AssertRefused(c, 3);
        Assert.Equal(2, c.FailedAt);
        Assert.Contains("row 2", Assert.Throws<InvalidOperationException>(() => c.Solve(new double[3])).Message);
    }

    [Fact]
    public void AFailedInstance_FactorsTheNextMatrixAsAFreshOneDoes()
    {
        const int n = 50;
        double[,] a = DenseSpd(n, new Xoshiro256PlusPlus(seed: 911));
        var reused = new CholeskyDecomposition(n);
        reused.Decompose(SquaredExponential(n, 0.1));
        Assert.True(double.IsNaN(reused.LogDet));
        reused.Decompose(a);
        var fresh = new CholeskyDecomposition(n);
        fresh.Decompose(a);
        Assert.Equal(-1, reused.FailedAt);
        Assert.Equal(fresh.LogDet, reused.LogDet);
        var x = new double[n * n];
        var y = new double[n * n];
        reused.WriteLTo(x);
        fresh.WriteLTo(y);
        Assert.Equal(y, x);
        var p = new double[n, n];
        var q = new double[n, n];
        reused.WriteInverseTo(p);
        fresh.WriteInverseTo(q);
        Assert.Equal(q, p);
    }

    // Floored factors that stay finite are read, not refused: the floor's rule is not this class's to change here.
    [Theory]
    [InlineData(50, 1.0)]
    [InlineData(401, 0.5)]
    public void AFlooredFactorThatStaysFinite_IsRead(int n, double length)
    {
        var c = new CholeskyDecomposition(n);
        c.Decompose(SquaredExponential(n, length));
        Assert.True(c.Status.Finite);
        Assert.True(c.Status.FlooredPivots > 0);
        Assert.Equal(-1, c.FailedAt);
        Assert.True(double.IsFinite(c.LogDet));
        Assert.All(c.Solve(new double[n]), v => Assert.Equal(0.0, v));
    }

    // Every finite factor, floored or not, against the class before 2026-11-13, verbatim: the factor, LogDet, the inverse,
    // a solve and a sample, bit for bit.
    [Fact]
    [Trait("Tier", "Pinned")]
    public void EveryFiniteFactor_ReadsAsTheEarlierClassBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1113);
        var cases = new System.Collections.Generic.List<double[,]>();
        for (int t = 0; t < 60; t++)
        {
            int n = 1 + (int)(rng.NextUInt64() % 40);
            double[,] a = DenseSpd(n, rng);
            double scale = Math.ScaleB(1.0, (int)(rng.NextUInt64() % 120) - 80);   // 2^-80 to 2^39: the floor reached below about 2^-39
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i, j] *= scale;
            cases.Add(a);
        }
        cases.Add(SquaredExponential(50, 1.0));
        cases.Add(SquaredExponential(401, 0.5));
        var identity = new double[40, 40];
        for (int i = 0; i < 40; i++) identity[i, i] = 1e-13;
        cases.Add(identity);
        int floored = 0;
        foreach (double[,] a in cases)
        {
            int n = a.GetLength(0);
            var c = new CholeskyDecomposition(n);
            c.Decompose(a);
            Assert.Equal(-1, c.FailedAt);
            if (c.Status.FlooredPivots > 0) floored++;
            var earlier = new EarlierCholesky(n);
            earlier.Decompose(a);
            Assert.Equal(Bits(earlier.LogDet), Bits(c.LogDet));
            var l = new double[n * n];
            c.WriteLTo(l);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Assert.Equal(Bits(earlier.L[i, j]), Bits(l[i * n + j]));
            var inverse = new double[n, n];
            var earlierInverse = new double[n, n];
            c.WriteInverseTo(inverse);
            earlier.WriteInverseTo(earlierInverse);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Assert.Equal(Bits(earlierInverse[i, j]), Bits(inverse[i, j]));
            double[] b = Vector(n, rng);
            Assert.Equal(Array.ConvertAll(earlier.Solve(b), Bits), Array.ConvertAll(c.Solve(b), Bits));
            Assert.Equal(Array.ConvertAll(earlier.Sample(new Random(n), b), Bits), Array.ConvertAll(c.Sample(new Random(n), b), Bits));
        }
        Assert.True(floored >= 10, $"Only {floored} of the cases reached the floor.");
    }

    private static long Bits(double v) => BitConverter.DoubleToInt64Bits(v);

    // CholeskyDecomposition as it stood before 2026-11-13, verbatim but for its fields' names.
    private sealed class EarlierCholesky(int dim)
    {
        public readonly double[,] L = new double[dim, dim];
        private readonly double[,] _lInv = new double[dim, dim];
        public double LogDet { get; private set; }

        public void Decompose(double[,] sigma)
        {
            double logDetL = 0.0;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j <= i; j++)
                {
                    double sum = sigma[i, j];
                    for (int k = 0; k < j; k++) sum -= L[i, k] * L[j, k];
                    if (i == j)
                    {
                        if (sum <= 1e-12) sum = 1e-12;
                        L[i, i] = Math.Sqrt(sum);
                        logDetL += Math.Log(L[i, i]);
                    }
                    else
                    {
                        L[i, j] = sum / L[j, j];
                        L[j, i] = 0.0;
                    }
                }
            LogDet = 2.0 * logDetL;
            for (int i = 0; i < dim; i++)
            {
                _lInv[i, i] = 1.0 / L[i, i];
                for (int j = 0; j < i; j++)
                {
                    double sum = 0.0;
                    for (int k = j; k < i; k++) sum -= L[i, k] * _lInv[k, j];
                    _lInv[i, j] = sum / L[i, i];
                    _lInv[j, i] = 0.0;
                }
            }
        }

        public void WriteInverseTo(double[,] target)
        {
            for (int i = 0; i < dim; i++)
                for (int j = i; j < dim; j++)
                {
                    double sum = 0.0;
                    for (int k = j; k < dim; k++) sum += _lInv[k, i] * _lInv[k, j];
                    target[i, j] = sum;
                    target[j, i] = sum;
                }
        }

        public double[] Solve(double[] b)
        {
            var x = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                double sum = b[i];
                for (int k = 0; k < i; k++) sum -= L[i, k] * x[k];
                x[i] = sum / L[i, i];
            }
            for (int i = dim - 1; i >= 0; i--)
            {
                double sum = x[i];
                for (int k = i + 1; k < dim; k++) sum -= L[k, i] * x[k];
                x[i] = sum / L[i, i];
            }
            return x;
        }

        public double[] Sample(Random rng, double[] mean)
        {
            double[] z = new double[dim];
            for (int i = 0; i < dim - 1; i += 2)
            {
                double u1 = 1.0 - rng.NextDouble();
                double u2 = 1.0 - rng.NextDouble();
                double mag = Math.Sqrt(-2.0 * Math.Log(u1));
                z[i] = mag * Math.Cos(2.0 * Math.PI * u2);
                z[i + 1] = mag * Math.Sin(2.0 * Math.PI * u2);
            }
            if (dim % 2 == 1)
            {
                double u1 = 1.0 - rng.NextDouble();
                double u2 = 1.0 - rng.NextDouble();
                z[dim - 1] = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            }
            double[] x = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                x[i] = mean[i];
                for (int j = 0; j <= i; j++) x[i] += L[i, j] * z[j];
            }
            return x;
        }
    }
}
