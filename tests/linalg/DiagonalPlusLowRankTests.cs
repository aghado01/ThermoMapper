using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// Σ = D + U·diag(s)·Uᵀ against its dense form, factored strictly: log-determinant, solves and bᵀΣ⁻¹b, including
/// zero weights, rank-deficient and wide factors; both draw maps have exactly the covariances they claim; and along a
/// heavy direction the quadratic form keeps its digits where the Woodbury subtraction loses all of them.
/// </summary>
public sealed class DiagonalPlusLowRankTests
{
    [Theory]
    [InlineData(1, 0)]
    [InlineData(5, 2)]
    [InlineData(40, 0)]
    [InlineData(60, 4)]
    [InlineData(30, 6)]
    [InlineData(3, 5)]
    public void MatchesTheDenseCovariance(int m, int r)
    {
        var rng = new Xoshiro256PlusPlus(seed: 1000 + 10 * m + r);
        (double[] d, double[,] u, double[] s) = Components(m, r, rng);
        if (r >= 2)
        {
            s[1] = 0.0;                                              // a zero weight
            for (int i = 0; i < m; i++) u[i, r - 1] = u[i, 0];       // a repeated column
        }
        var sigma = new DiagonalPlusLowRank(d, u, s);
        double[,] dense = Dense(d, u, s);
        var reference = new DenseCholesky(m);
        reference.Decompose(dense);
        Assert.True(reference.Status.PositiveDefinite);

        Assert.True(Math.Abs(sigma.LogDet - reference.LogDet) <= 1e-11 * Math.Max(1.0, Math.Abs(reference.LogDet)),
                    $"ln|Σ| {sigma.LogDet:R} vs {reference.LogDet:R}");
        double[] b = Vector(m, rng);
        AssertClose(reference.Solve(b), sigma.Solve(b), 1e-10);
        double form = reference.InverseQuadraticForm(b);
        Assert.True(Math.Abs(sigma.InverseQuadraticForm(b) - form) <= 1e-11 * form, $"bᵀΣ⁻¹b {sigma.InverseQuadraticForm(b):R} vs {form:R}");
    }

    [Fact]
    public void Draws_HaveExactlyTheirCovariances()
    {
        // The sample map M and innovation map N are linear: Σ_k M e_k (M e_k)ᵀ = Σ and Σ·Σ_k N e_k (N e_k)ᵀ = I.
        const int m = 12, r = 3;
        var rng = new Xoshiro256PlusPlus(seed: 1011);
        (double[] d, double[,] u, double[] s) = Components(m, r, rng);
        var sigma = new DiagonalPlusLowRank(d, u, s);
        double[,] dense = Dense(d, u, s);

        var covariance = new double[m, m];
        var precision = new double[m, m];
        var e = new double[m];
        var x = new double[m];
        var v = new double[m];
        for (int k = 0; k < m; k++)
        {
            Array.Clear(e);
            e[k] = 1.0;
            sigma.Sample(e, x);
            sigma.SampleInnovation(e, v);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                {
                    covariance[i, j] += x[i] * x[j];
                    precision[i, j] += v[i] * v[j];
                }
        }
        double scale = 0.0;
        for (int i = 0; i < m; i++) scale = Math.Max(scale, dense[i, i]);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                Assert.True(Math.Abs(covariance[i, j] - dense[i, j]) <= 1e-13 * scale, $"Cov[{i},{j}] {covariance[i, j]:R} vs {dense[i, j]:R}");
        AssertIdentity(Multiply(dense, precision), 1e-12);
    }

    [Fact]
    public void HeavyDirection_KeepsItsQuadraticForm()
    {
        // D = I, one column u with weight 10¹⁶, b = u/‖u‖: bᵀΣ⁻¹b = 1/(1 + 10¹⁶‖u‖²), about 10⁻¹⁷. The Woodbury form,
        // ‖b‖² − zᵀ(1 + λ²)⁻¹z, subtracts two numbers within 10⁻¹⁷ of 1 and keeps none of it.
        const int m = 50;
        var rng = new Xoshiro256PlusPlus(seed: 1012);
        double[] column = Vector(m, rng);
        var u = new double[m, 1];
        double norm2 = 0.0;
        for (int i = 0; i < m; i++)
        {
            u[i, 0] = column[i];
            norm2 += column[i] * column[i];
        }
        const double weight = 1e16;
        var d = new double[m];
        Array.Fill(d, 1.0);
        var sigma = new DiagonalPlusLowRank(d, u, [weight]);

        var b = new double[m];
        for (int i = 0; i < m; i++) b[i] = column[i] / Math.Sqrt(norm2);
        double exact = 1.0 / (1.0 + weight * norm2);
        double q = sigma.InverseQuadraticForm(b);
        Assert.True(Math.Abs(q - exact) <= 1e-10 * exact, $"bᵀΣ⁻¹b {q:R} vs {exact:R}");

        double bb = 0.0, z = 0.0;
        for (int i = 0; i < m; i++)
        {
            bb += b[i] * b[i];
            z += Math.Sqrt(weight) * column[i] * b[i];
        }
        double woodbury = bb - z * z / (1.0 + weight * norm2);
        Assert.True(Math.Abs(woodbury - exact) > 0.5 * exact, $"the Woodbury form kept its digits: {woodbury:R}");
    }

    [Fact]
    public void Construction_RefusesANonPositiveDiagonal_AndNegativeWeights()
    {
        var u = new double[2, 1];
        Assert.Throws<ArgumentException>(() => new DiagonalPlusLowRank([1.0, 0.0], u, [1.0]));
        Assert.Throws<ArgumentException>(() => new DiagonalPlusLowRank([1.0, double.NaN], u, [1.0]));
        Assert.Throws<ArgumentException>(() => new DiagonalPlusLowRank([1.0, 1.0], u, [-1.0]));
        Assert.Throws<ArgumentException>(() => new DiagonalPlusLowRank([1.0, 1.0, 1.0], u, [1.0]));
    }

    private static (double[] D, double[,] U, double[] S) Components(int m, int r, Xoshiro256PlusPlus rng)
    {
        var d = new double[m];
        for (int i = 0; i < m; i++) d[i] = 0.1 + rng.NextDouble();
        var u = new double[m, r];
        for (int i = 0; i < m; i++) for (int k = 0; k < r; k++) u[i, k] = rng.NextDouble() - 0.5;
        var s = new double[r];
        for (int k = 0; k < r; k++) s[k] = Math.Pow(10.0, 3.0 * rng.NextDouble());
        return (d, u, s);
    }

    private static double[,] Dense(double[] d, double[,] u, double[] s)
    {
        int m = d.Length, r = s.Length;
        var a = new double[m, m];
        for (int i = 0; i < m; i++)
        {
            a[i, i] += d[i];
            for (int j = 0; j < m; j++)
                for (int k = 0; k < r; k++) a[i, j] += u[i, k] * s[k] * u[j, k];
        }
        return a;
    }
}
