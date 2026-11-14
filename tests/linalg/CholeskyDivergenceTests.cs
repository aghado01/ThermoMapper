using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// <see cref="CholeskyDecomposition"/> without its pivot floor: the first pivot that is not finite and positive stops it at
/// <see cref="DenseCholesky"/>'s row, an L⁻¹ that overflows fails it, and a failed factor is refused. On every input the
/// two agree, failure row, <see cref="FactorStatus"/>, factor and <c>LogDet</c>; where the floor never reached a pivot, the
/// class reads bit for bit as before (shared-code defects brief, slice 2 and its follow-up).
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

    // L bidiagonal, 2^-20 on its diagonal and 1 below it, so A = LLᵀ is positive definite and factors exactly in powers of
    // two; L⁻¹'s first column is (−1)^i·2^(20(i+1)), past double.MaxValue from row 51.
    internal static double[,] OverflowingInverse(int n)
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

    private static void AssertRefused(CholeskyDecomposition c, int n)
    {
        Assert.True(double.IsNaN(c.LogDet));
        Assert.Throws<InvalidOperationException>(() => c.Solve(new double[n]));
        Assert.Throws<InvalidOperationException>(() => c.WriteInverseTo(new double[n, n]));
        Assert.Throws<InvalidOperationException>(() => c.WriteLTo(new double[n * n]));
        Assert.Throws<InvalidOperationException>(() => c.Sample(new Random(1), new double[n]));
    }

    private static IEnumerable<double[,]> Ladder()
    {
        var rng = new Xoshiro256PlusPlus(seed: 2026);
        for (int t = 0; t < 60; t++)
        {
            int n = 1 + (int)(rng.NextUInt64() % 40);
            double[,] a = DenseSpd(n, rng);
            double scale = Math.ScaleB(1.0, (int)(rng.NextUInt64() % 120) - 80);   // 2^-80 to 2^39, past the old floor
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i, j] *= scale;
            yield return a;
        }
        foreach (int n in new[] { 50, 100, 200, 401 })
            foreach (double length in new[] { 0.05, 0.1, 0.5, 1.0 })
                yield return SquaredExponential(n, length);
        var tiny = new double[40, 40];
        for (int i = 0; i < 40; i++) tiny[i, i] = 1e-13;
        yield return tiny;
        yield return new double[,] { { 1, 2 }, { 2, 1 } };
        yield return new double[,] { { 1e-300, 1e300 }, { 1e300, 1 } };
        yield return new double[,] { { 2, 0, 0 }, { 0, 2, 0 }, { 0, double.NaN, 2 } };
    }

    [Fact]
    public void EveryInput_FailsWhereTheStrictFactorFails_AndOtherwiseIsItsFactor()
    {
        int failed = 0, held = 0;
        foreach (double[,] a in Ladder())
        {
            int n = a.GetLength(0);
            var c = new CholeskyDecomposition(n);
            c.Decompose(a);
            var d = new DenseCholesky(n);
            d.Decompose(a);
            Assert.Equal(d.FailedAt, c.FailedAt);
            Assert.Equal(d.Status, c.Status);
            Assert.Equal(BitConverter.DoubleToInt64Bits(d.LogDet), BitConverter.DoubleToInt64Bits(c.LogDet));
            if (c.FailedAt >= 0)
            {
                AssertRefused(c, n);
                failed++;
                continue;
            }
            held++;
            var l = new double[n * n];
            c.WriteLTo(l);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Assert.Equal(d.Factor(i, j), l[i * n + j]);
            Assert.Equal(0, c.Status.FlooredPivots);
        }
        Assert.True(failed >= 15 && held >= 60, $"{failed} failed, {held} held.");
    }

    [Theory]
    [InlineData(50, 0.1)]
    [InlineData(100, 0.5)]
    [InlineData(401, 0.05)]
    public void ASmoothKernelBelowRounding_IsRefused(int n, double length)
    {
        var c = new CholeskyDecomposition(n);
        c.Decompose(SquaredExponential(n, length));
        Assert.False(c.Status.PositiveDefinite);
        Assert.InRange(c.FailedAt, 1, n - 1);
        AssertRefused(c, n);
    }

    [Fact]
    public void ASmallScaleIdentity_IsFactoredExactly()
    {
        // The old floor read every variance of 1e-13·I as 1e-12, ten times too large.
        const int n = 40;
        var a = new double[n, n];
        for (int i = 0; i < n; i++) a[i, i] = 1e-13;
        var c = new CholeskyDecomposition(n);
        c.Decompose(a);
        Assert.Equal(-1, c.FailedAt);
        Assert.Equal(n * Math.Log(1e-13), c.LogDet, 9);
        double[] x = c.Solve([.. new double[n - 1], 1.0]);
        Assert.True(Math.Abs(x[n - 1] - 1e13) <= 1e-14 * 1e13, $"{x[n - 1]:R}");
        var s = c.Sample(new Random(3), new double[n]);
        Assert.All(s, v => Assert.InRange(Math.Abs(v), 0.0, 1e-5));
    }

    [Fact]
    public void APivotThatIsNotFinite_IsRefused()
    {
        var c = new CholeskyDecomposition(2);
        c.Decompose(new double[,] { { 1e-300, 1e300 }, { 1e300, 1 } });
        Assert.Equal(1, c.FailedAt);
        Assert.False(c.Status.Finite);
        AssertRefused(c, 2);
    }

    [Fact]
    public void AnInverseFactorThatOverflows_IsRefused_ThoughThePivotsArePositive()
    {
        const int n = 60;
        var c = new CholeskyDecomposition(n);
        c.Decompose(OverflowingInverse(n));
        Assert.True(c.Status.PositiveDefinite);
        Assert.Equal(51, c.FailedAt);
        AssertRefused(c, n);
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

    // Where the old floor never reached a pivot, the class against itself before 2026-11-13, verbatim: the factor, LogDet,
    // the inverse, a solve and a sample, bit for bit.
    [Fact]
    [Trait("Tier", "Pinned")]
    public void EveryFactorTheFloorNeverReached_ReadsAsTheEarlierClassBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1113);
        for (int t = 0; t < 60; t++)
        {
            int n = 1 + (int)(rng.NextUInt64() % 40);
            double[,] a = DenseSpd(n, rng);
            double scale = Math.ScaleB(1.0, (int)(rng.NextUInt64() % 78) - 38);   // 2^-38 to 2^39: every pivot above 1e-12
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i, j] *= scale;
            var c = new CholeskyDecomposition(n);
            c.Decompose(a);
            Assert.Equal(-1, c.FailedAt);
            Assert.True(c.Status.MinPivot > 1e-12);
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
