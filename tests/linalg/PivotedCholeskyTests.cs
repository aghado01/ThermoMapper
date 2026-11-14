using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// <see cref="PivotedCholesky"/>: smooth kernels' covariances, singular to working precision, factored to their numerical
/// rank and reproduced; an indefinite remainder refused, whether on the diagonal or off it; rounding-level truncation kept;
/// and the loadings, bit for bit, of the routine BGP's drawn function carried before (shared-code defects brief, slice 2's
/// follow-up).
/// </summary>
public sealed class PivotedCholeskyTests
{
    // A squared-exponential covariance, variance 1, on n evenly spaced points of [0, 1].
    private static double[,] Kernel(int n, double length)
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

    private static double Reconstruction(PivotedCholesky f, double[,] a)
    {
        int n = a.GetLength(0);
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < f.Rank; k++) s += f.Loading(i, k) * f.Loading(j, k);
                num += (a[i, j] - s) * (a[i, j] - s);
                den += a[i, j] * a[i, j];
            }
        return Math.Sqrt(num / den);
    }

    [Theory]
    [InlineData(50, 0.1, 31)]
    [InlineData(100, 0.5, 11)]
    [InlineData(200, 0.2, 19)]
    [InlineData(401, 0.05, 57)]
    [InlineData(401, 1.0, 8)]
    public void ASmoothKernel_FactorsToItsNumericalRank(int n, double length, int rank)
    {
        double[,] a = Kernel(n, length);
        var f = new PivotedCholesky(n);
        f.Decompose(a);
        Assert.True(f.PositiveSemidefinite, $"remainder {f.RelativeRemainder:E3}");
        Assert.Equal(rank, f.Rank);
        Assert.True(Reconstruction(f, a) <= 2e-13, $"{Reconstruction(f, a):E3}");
        Assert.True(f.RelativeRemainder <= 1e-12);
    }

    [Theory]
    [InlineData(new[] { 1.0, 2.0, 2.0, 1.0 })]       // the remainder's diagonal is −3
    [InlineData(new[] { 0.0, 1.0, 1.0, 0.0 })]       // diagonals 0, an off-diagonal 1: no pivot, and not semidefinite
    [InlineData(new[] { 4.0, 0.0, 0.0, -1e-3 })]     // a negative diagonal past rounding
    [InlineData(new[] { 1.0, double.NaN, double.NaN, 1.0 })]
    public void AnIndefiniteOrNonFiniteMatrix_IsRefused(double[] entries)
    {
        var a = new double[,] { { entries[0], entries[1] }, { entries[2], entries[3] } };
        var f = new PivotedCholesky(2);
        f.Decompose(a);
        Assert.False(f.PositiveSemidefinite);
        Assert.Throws<InvalidOperationException>(() => f.Loading(0, 0));
        Assert.Throws<InvalidOperationException>(() => f.Map(new double[f.Rank], new double[2]));
    }

    [Fact]
    public void TheRemainderOfOneTwoTwoOne_IsReported()
    {
        var f = new PivotedCholesky(2);
        f.Decompose(new double[,] { { 1, 2 }, { 2, 1 } });
        Assert.Equal(1, f.Rank);
        Assert.Equal(3.0, f.RelativeRemainder, 12);
        Assert.Contains("materially indefinite", Assert.Throws<InvalidOperationException>(() => f.Loading(0, 0)).Message);
    }

    [Fact]
    public void SemidefiniteMatrices_AreKept_AtRoundingLevelTruncation()
    {
        // The zero matrix: rank 0, nothing left.
        var f = new PivotedCholesky(3);
        f.Decompose(new double[3, 3]);
        Assert.True(f.PositiveSemidefinite);
        Assert.Equal(0, f.Rank);

        // A Gram matrix of rank 3 in 12 dimensions, formed in floating point: its zero eigenvalues are rounding.
        var rng = new Xoshiro256PlusPlus(seed: 31);
        const int n = 12;
        var v = new double[n, 3];
        for (int i = 0; i < n; i++) for (int c = 0; c < 3; c++) v[i, c] = rng.NextDouble() - 0.5;
        var g = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int c = 0; c < 3; c++) g[i, j] += v[i, c] * v[j, c];
        var h = new PivotedCholesky(n);
        h.Decompose(g);
        Assert.True(h.PositiveSemidefinite, $"remainder {h.RelativeRemainder:E3}");
        Assert.Equal(3, h.Rank);
        Assert.True(Reconstruction(h, g) <= 1e-14);

        // A remainder [[0, δ], [δ, 0]] at δ = 5·10⁻¹¹ of the largest diagonal: indefinite only at the rounding a computed
        // covariance carries, past the 10⁻¹² stop. Kept at the default, √ε; refused when the caller asks for 10⁻¹².
        const double delta = 5e-11;
        var noisy = new double[,] { { 1, 0, 0 }, { 0, 0, delta }, { 0, delta, 0 } };
        var m = new PivotedCholesky(3);
        m.Decompose(noisy);
        Assert.True(m.PositiveSemidefinite, $"remainder {m.RelativeRemainder:E3}");
        Assert.Equal(1, m.Rank);
        Assert.Equal(delta, m.RelativeRemainder);
        m.Decompose(noisy, indefiniteTolerance: 1e-12);
        Assert.False(m.PositiveSemidefinite);

        // Positive definite: full rank, and the map is L·z.
        double[,] a = DenseSpd(n, rng);
        var full = new PivotedCholesky(n);
        full.Decompose(a);
        Assert.True(full.PositiveSemidefinite);
        Assert.Equal(n, full.Rank);
        double[] z = Vector(n, rng), x = new double[n];
        full.Map(z, x);
        for (int i = 0; i < n; i++)
        {
            double s = 0.0;
            for (int k = 0; k < n; k++) s += full.Loading(i, k) * z[k];
            Assert.Equal(s, x[i]);
        }
    }

    [Fact]
    public void AnUnfactoredInstance_IsRefused_AndArgumentsAreChecked()
    {
        var f = new PivotedCholesky(3);
        Assert.Contains("No matrix", Assert.Throws<InvalidOperationException>(() => f.Loading(0, 0)).Message);
        Assert.Throws<ArgumentException>(() => f.Decompose(new double[4, 4]));
        Assert.Throws<ArgumentException>(() => f.Decompose(new double[2, 3]));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.Decompose(new double[2, 2], 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.Decompose(new double[2, 2], 1e-12, -1.0));
        f.Decompose(new double[,] { { 4, 2 }, { 2, 2 } });
        Assert.Throws<ArgumentException>(() => f.Map(new double[1], new double[2]));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.Loading(0, 2));
    }

    // The routine BGP's drawn function carried before 2026-11-13, verbatim: its rank and loadings, bit for bit.
    [Fact]
    [Trait("Tier", "Pinned")]
    public void TheLoadings_AreTheEarlierRoutinesBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1113);
        for (int t = 0; t < 40; t++)
        {
            int n = 2 + (int)(rng.NextUInt64() % 60);
            double[,] a = t % 2 == 0 ? Kernel(n, 0.05 + rng.NextDouble()) : LowRank(n, 1 + (int)(rng.NextUInt64() % (ulong)n), rng);
            var f = new PivotedCholesky(n);
            f.Decompose(a);
            var earlier = new double[n * n];
            int rank = EarlierPivotedFactor(a, earlier);
            Assert.True(f.PositiveSemidefinite, $"case {t}: remainder {f.RelativeRemainder:E3}");
            Assert.Equal(rank, f.Rank);
            for (int i = 0; i < n; i++)
                for (int k = 0; k < rank; k++)
                    Assert.Equal(BitConverter.DoubleToInt64Bits(earlier[i * n + k]), BitConverter.DoubleToInt64Bits(f.Loading(i, k)));
        }
    }

    private static double[,] LowRank(int n, int r, Xoshiro256PlusPlus rng)
    {
        var v = new double[n, r];
        for (int i = 0; i < n; i++) for (int c = 0; c < r; c++) v[i, c] = rng.NextDouble() - 0.5;
        var g = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int c = 0; c < r; c++) g[i, j] += v[i, c] * v[j, c];
        return g;
    }

    private static int EarlierPivotedFactor(double[,] lower, double[] factor)
    {
        int p = lower.GetLength(0);
        Array.Clear(factor, 0, p * p);
        var remaining = new double[p];
        var order = new int[p];
        double largest = 0.0;
        for (int a = 0; a < p; a++)
        {
            remaining[a] = lower[a, a];
            order[a] = a;
            largest = Math.Max(largest, remaining[a]);
        }
        double tolerance = 1e-12 * largest;
        int k = 0;
        for (; k < p; k++)
        {
            int best = k;
            for (int i = k + 1; i < p; i++) if (remaining[order[i]] > remaining[order[best]]) best = i;
            double pivot = remaining[order[best]];
            if (!(pivot > tolerance)) break;
            (order[k], order[best]) = (order[best], order[k]);
            int at = order[k];
            double root = Math.Sqrt(pivot);
            factor[at * p + k] = root;
            for (int i = k + 1; i < p; i++)
            {
                int q = order[i];
                double s = q >= at ? lower[q, at] : lower[at, q];
                for (int c = 0; c < k; c++) s -= factor[q * p + c] * factor[at * p + c];
                double v = s / root;
                factor[q * p + k] = v;
                remaining[q] -= v * v;
            }
        }
        return k;
    }
}
