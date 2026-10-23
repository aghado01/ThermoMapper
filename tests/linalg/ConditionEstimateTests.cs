using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The equilibrated condition estimate (<see cref="ISpdFactor.EquilibratedCondition"/>): checked against κ₁ of
/// D^(−½)·A·D^(−½) from a dense Gauss–Jordan inverse. It never exceeds that value, stays within the factor of 3 the
/// method aims for, agrees between the band and dense factors, ignores column scaling, and measures the matrix whose
/// pivots look healthy while its condition number is about 10¹⁵, which the relative pivot cannot see.
/// </summary>
public sealed class ConditionEstimateTests
{
    // A = L·Lᵀ for a random lower-banded L: symmetric positive definite with half-bandwidth p (p = n−1 is dense).
    private static double[,] Spd(int n, int p, double diagonal, Xoshiro256PlusPlus rng)
    {
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = Math.Max(0, i - p); j <= i; j++)
                l[i, j] = i == j ? diagonal + rng.NextDouble() : rng.NextDouble() - 0.5;
        return Product(l);
    }

    private static double[,] Product(double[,] l)
    {
        int n = l.GetLength(0);
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < n; k++) s += l[i, k] * l[j, k];
                a[i, j] = s;
            }
        return a;
    }

    // κ₁ of the equilibrated matrix, from its explicit inverse.
    private static double ExactEquilibratedCondition(double[,] a)
    {
        int n = a.GetLength(0);
        var e = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) e[i, j] = a[i, j] / Math.Sqrt(a[i, i] * a[j, j]);
        return Norm1(e) * Norm1(Invert(e));
    }

    private static double Norm1(double[,] m)
    {
        int n = m.GetLength(0);
        double best = 0.0;
        for (int j = 0; j < n; j++)
        {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += Math.Abs(m[i, j]);
            best = Math.Max(best, s);
        }
        return best;
    }

    private static double[,] Invert(double[,] a0)
    {
        int n = a0.GetLength(0);
        var a = (double[,])a0.Clone();
        var inv = new double[n, n];
        for (int i = 0; i < n; i++) inv[i, i] = 1.0;
        for (int c = 0; c < n; c++)
        {
            int piv = c;
            for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
            for (int k = 0; k < n; k++) { (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]); (inv[c, k], inv[piv, k]) = (inv[piv, k], inv[c, k]); }
            double d = a[c, c];
            for (int k = 0; k < n; k++) { a[c, k] /= d; inv[c, k] /= d; }
            for (int r = 0; r < n; r++)
            {
                if (r == c) continue;
                double f = a[r, c];
                for (int k = 0; k < n; k++) { a[r, k] -= f * a[c, k]; inv[r, k] -= f * inv[c, k]; }
            }
        }
        return inv;
    }

    private static BandCholesky Band(double[,] a, int p, BandFactorization mode)
    {
        var f = new BandCholesky(a.GetLength(0), p, mode, strict: true);
        f.Decompose(a);
        return f;
    }

    private static DenseCholesky Dense(double[,] a)
    {
        var f = new DenseCholesky(a.GetLength(0));
        f.Decompose(a);
        return f;
    }

    public static TheoryData<int, int, double> Matrices => new()
    {
        { 6, 1, 1.0 }, { 12, 3, 1.0 }, { 12, 3, 0.05 }, { 30, 3, 0.02 }, { 20, 19, 0.3 }, { 8, 7, 0.01 },
    };

    [Theory]
    [MemberData(nameof(Matrices))]
    public void Estimate_IsALowerBound_WithinAFactorOfThree(int n, int p, double diagonal)
    {
        var rng = new Xoshiro256PlusPlus(seed: 1000 + n * 10 + p);
        for (int trial = 0; trial < 20; trial++)
        {
            double[,] a = Spd(n, p, diagonal, rng);
            double exact = ExactEquilibratedCondition(a);
            double slack = 1e-10 + 1e-14 * exact;   // the reference inverse itself loses about log₁₀κ digits
            foreach (ISpdFactor f in new ISpdFactor[] { Band(a, p, BandFactorization.Cholesky), Band(a, p, BandFactorization.Ldlt), Dense(a) })
            {
                double estimate = f.EquilibratedCondition();
                Assert.True(estimate <= exact * (1.0 + slack), $"estimate {estimate:R} above the exact {exact:R}");
                Assert.True(estimate >= exact / 3.0, $"estimate {estimate:G6} below a third of the exact {exact:G6}");
            }
        }
    }

    [Fact]
    public void BandAndDense_Agree_AndIgnoreColumnScaling()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        double[,] a = Spd(15, 3, 0.1, rng);
        double band = Band(a, 3, BandFactorization.Ldlt).EquilibratedCondition();
        Assert.Equal(band, Dense(a).EquilibratedCondition(), band * 1e-10);

        // S·A·S for a diagonal S spanning twelve orders of magnitude has the same equilibrated matrix.
        int n = a.GetLength(0);
        var scaled = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) scaled[i, j] = a[i, j] * Math.Pow(10.0, (i - 7) * 0.8) * Math.Pow(10.0, (j - 7) * 0.8);
        Assert.Equal(band, Band(scaled, 3, BandFactorization.Ldlt).EquilibratedCondition(), band * 1e-8);
    }

    [Fact]
    public void HealthyPivots_DoNotHideAnIllConditionedMatrix()
    {
        // The review's counterexample: A = L·Lᵀ with L lower bidiagonal, diagonal 1 and subdiagonal 2. Every pivot is 1,
        // so MinRelativePivot = 0.2, yet (−½)^j is nearly a null vector and κ(E) is about 2·10¹⁵.
        const int n = 25;
        var l = new double[n, n];
        for (int i = 0; i < n; i++) { l[i, i] = 1.0; if (i > 0) l[i, i - 1] = 2.0; }
        double[,] a = Product(l);
        foreach (ISpdFactor f in new ISpdFactor[] { Band(a, 1, BandFactorization.Ldlt), Dense(a) })
        {
            Assert.Equal(0.2, f.Status.MinRelativePivot, 12);
            Assert.InRange(f.EquilibratedCondition(), 1e14, 1e17);
        }
    }

    [Fact]
    public void FailedStrictFactor_Refuses()
    {
        var a = new double[,] { { 1.0, 2.0 }, { 2.0, 1.0 } };   // indefinite
        var band = new BandCholesky(2, 1, BandFactorization.Ldlt, strict: true);
        band.Decompose(a);
        Assert.Throws<InvalidOperationException>(() => band.EquilibratedCondition());
        var dense = new DenseCholesky(2);
        dense.Decompose(a);
        Assert.Throws<InvalidOperationException>(() => dense.EquilibratedCondition());
    }

    [Fact]
    public void FlooredFactor_Refuses()
    {
        // [1e-20] is clamped to the 1e-12 floor, so the solves describe [1e-12] while the scaling reads 1e-20; the mix
        // reported 1e-8 for a matrix whose condition number is 1.
        var tiny = new BandCholesky(1, 0, BandFactorization.Ldlt);
        tiny.Decompose(new double[,] { { 1e-20 } });
        Assert.Equal(1, tiny.Status.FlooredPivots);
        Assert.Throws<InvalidOperationException>(() => tiny.EquilibratedCondition());

        var strict = new BandCholesky(1, 0, BandFactorization.Ldlt, strict: true);
        strict.Decompose(new double[,] { { 1e-20 } });
        Assert.Equal(1.0, strict.EquilibratedCondition(), 15);
    }

    [Fact]
    public void OrderOne_IsExactlyOne()
    {
        var a = new double[,] { { 4.0 } };
        Assert.Equal(1.0, Band(a, 0, BandFactorization.Cholesky).EquilibratedCondition(), 15);
        Assert.Equal(1.0, Dense(a).EquilibratedCondition(), 15);
    }
}
