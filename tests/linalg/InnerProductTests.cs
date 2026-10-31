using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The output-side inner product (F-X6, <see cref="IInnerProduct"/>) in the Khatri–Rao contractions.
/// <list type="bullet">
/// <item>Q = I is today's arithmetic, bit for bit: the identity product is <see cref="PairwiseSum.Dot"/>, and a declared
/// identity gives the contractions the entries a null product gives.</item>
/// <item>A tridiagonal Q (the Markov residual's precision) gives G = (AᵀA)∘(HᵀQH), its border and h_j = a_jᵀ·Y·Q·h_j against
/// the materialized design X and the weight I ⊗ Q, XᵀWX and XᵀW vec(Y).</item>
/// <item>The product and its log-determinant carry the Gaussian log density of a residual with covariance I ⊗ R: −½Σᵢ yᵢᵀQyᵢ +
/// (n/2)·ln|Q| against the dense (nm)² covariance's.</item>
/// </list>
/// </summary>
public sealed class InnerProductTests
{
    // A tridiagonal SPD Q on m locations: the product acts on k stacked vectors, linear in m per vector.
    private sealed class Tridiagonal(double[] diagonal, double[] off) : IInnerProduct
    {
        public bool IsIdentity => false;
        public int Dimension => diagonal.Length;

        public double LogDeterminant
        {
            get
            {
                int m = diagonal.Length;
                var dense = new double[m, m];
                for (int i = 0; i < m; i++)
                {
                    dense[i, i] = diagonal[i];
                    if (i + 1 < m) dense[i + 1, i] = dense[i, i + 1] = off[i];
                }
                var cholesky = new DenseCholesky(m);
                cholesky.Decompose(dense);
                return cholesky.LogDet;
            }
        }

        public double Dot(ReadOnlySpan<double> u, ReadOnlySpan<double> v)
        {
            int m = diagonal.Length;
            if (u.Length != v.Length || u.Length % m != 0) throw new ArgumentException("Stacked vectors of the output length.");
            double sum = 0.0;
            for (int s = 0; s < u.Length; s += m)
                for (int i = 0; i < m; i++)
                {
                    double qv = diagonal[i] * v[s + i];
                    if (i > 0) qv += off[i - 1] * v[s + i - 1];
                    if (i + 1 < m) qv += off[i] * v[s + i + 1];
                    sum += u[s + i] * qv;
                }
            return sum;
        }

        public double[,] Dense()
        {
            int m = diagonal.Length;
            var q = new double[m, m];
            for (int i = 0; i < m; i++)
            {
                q[i, i] = diagonal[i];
                if (i + 1 < m) q[i + 1, i] = q[i, i + 1] = off[i];
            }
            return q;
        }
    }

    private static Tridiagonal MarkovQ(int m)
    {
        // The precision of a stationary AR(1) residual, ρ = 0.6, unit innovation variance.
        const double rho = 0.6;
        var d = new double[m];
        var e = new double[m - 1];
        for (int i = 0; i < m; i++) d[i] = i == 0 || i == m - 1 ? 1.0 : 1.0 + rho * rho;
        for (int i = 0; i < m - 1; i++) e[i] = -rho;
        return new Tridiagonal(d, e);
    }

    private static (double[][] A, double[][] H, double[,] Y) Problem(int n, int m, int terms, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var a = new double[terms][];
        var h = new double[terms][];
        for (int j = 0; j < terms; j++)
        {
            a[j] = new double[n];
            h[j] = new double[m];
            for (int i = 0; i < n; i++) a[j][i] = j == 0 ? 1.0 : Math.Max(rng.NextDouble() - 0.4, 0.0);
            for (int l = 0; l < m; l++) h[j][l] = j == 1 ? 1.0 : rng.NextDouble();
        }
        var y = new double[n, m];
        for (int i = 0; i < n; i++) for (int l = 0; l < m; l++) y[i, l] = 3.0 + rng.NextDouble();
        return (a, h, y);
    }

    private static void AssertSame(double[,] x, double[,] y)
    {
        Assert.Equal(x.GetLength(0), y.GetLength(0));
        for (int i = 0; i < x.GetLength(0); i++)
            for (int j = 0; j < x.GetLength(1); j++) Assert.Equal(x[i, j], y[i, j]);
    }

    [Fact]
    public void TheIdentity_IsPairwiseSum_AndADeclaredIdentityIsTodaysArithmetic()
    {
        IInnerProduct identity = IdentityInnerProduct.Instance;
        Assert.True(identity.IsIdentity);
        Assert.Equal(0, identity.Dimension);
        Assert.Equal(0.0, identity.LogDeterminant);
        var rng = new Xoshiro256PlusPlus(3);
        foreach (int length in new[] { 1, 7, 64, 1000 })
        {
            var u = new double[length];
            var v = new double[length];
            for (int i = 0; i < length; i++) { u[i] = rng.NextDouble() - 0.3; v[i] = 100.0 + rng.NextDouble(); }
            Assert.Equal(PairwiseSum.Dot(u, v), identity.Dot(u, v));
        }

        (double[][] a, double[][] h, double[,] y) = Problem(9, 11, 5, 21);
        AssertSame(KhatriRao.Gram(a, h), KhatriRao.Gram(a, h, identity));
        double[] plain = KhatriRao.Border(a[..4], h[..4], a[4], h[4], out double diagonal);
        double[] declared = KhatriRao.Border(a[..4], h[..4], a[4], h[4], out double declaredDiagonal, identity);
        Assert.Equal(plain, declared);
        Assert.Equal(diagonal, declaredDiagonal);
        Assert.Equal(KhatriRao.CrossProducts(a, h, y), KhatriRao.CrossProducts(a, h, y, identity));
    }

    [Fact]
    public void GramBorderAndCrossProducts_UnderQ_MatchTheWeightedMaterializedDesign()
    {
        const int n = 9, m = 11, terms = 5;
        (double[][] a, double[][] h, double[,] y) = Problem(n, m, terms, 22);
        Tridiagonal q = MarkovQ(m);
        double[,] dense = q.Dense();

        // X is n·m × terms, column j = a_j ⊗ h_j, and W = I_n ⊗ Q.
        var x = new double[n * m, terms];
        for (int j = 0; j < terms; j++)
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) x[i * m + l, j] = a[j][i] * h[j][l];
        double Weighted(int j, int k)
        {
            double sum = 0.0;
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++)
                    for (int r = 0; r < m; r++) sum += x[i * m + l, j] * dense[l, r] * x[i * m + r, k];
            return sum;
        }

        double[,] g = KhatriRao.Gram(a, h, q);
        for (int j = 0; j < terms; j++)
            for (int k = 0; k < terms; k++)
                Assert.True(Math.Abs(g[j, k] - Weighted(j, k)) <= 1e-12 * (1.0 + Math.Abs(g[j, k])), $"G[{j},{k}]: {g[j, k]} against {Weighted(j, k)}");

        double[] u = KhatriRao.Border(a[..(terms - 1)], h[..(terms - 1)], a[terms - 1], h[terms - 1], out double diagonal, q);
        for (int k = 0; k < terms - 1; k++) Assert.Equal(g[terms - 1, k], u[k], 12);
        Assert.Equal(g[terms - 1, terms - 1], diagonal, 12);

        double[] cross = KhatriRao.CrossProducts(a, h, y, q);
        for (int j = 0; j < terms; j++)
        {
            double expected = 0.0;
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++)
                    for (int r = 0; r < m; r++) expected += x[i * m + l, j] * dense[l, r] * y[i, r];
            Assert.True(Math.Abs(cross[j] - expected) <= 1e-11 * (1.0 + Math.Abs(expected)), $"h[{j}]: {cross[j]} against {expected}");
        }
        Assert.Throws<ArgumentException>(() => KhatriRao.Gram(a, h, MarkovQ(m + 1)));
    }

    [Fact]
    public void TheProductAndItsLogDeterminant_CarryTheGaussianLogDensityOfAnIOtimesRResidual()
    {
        const int n = 4, m = 6;
        Tridiagonal q = MarkovQ(m);
        (_, _, double[,] y) = Problem(n, m, 2, 23);
        var stacked = new double[n * m];
        for (int i = 0; i < n; i++) for (int l = 0; l < m; l++) stacked[i * m + l] = y[i, l] - 3.5;

        // The dense (nm)² covariance I ⊗ R, with R = Q⁻¹ from Q's columns.
        var qFactor = new DenseCholesky(m);
        qFactor.Decompose(q.Dense());
        var r = new double[m, m];
        for (int c = 0; c < m; c++)
        {
            var unit = new double[m];
            unit[c] = 1.0;
            double[] col = qFactor.Solve(unit);
            for (int i = 0; i < m; i++) r[i, c] = col[i];
        }
        var sigma = new double[n * m, n * m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
                for (int k = 0; k < m; k++) sigma[i * m + l, i * m + k] = r[l, k];
        var full = new DenseCholesky(n * m);
        full.Decompose(sigma);
        double dense = -0.5 * full.InverseQuadraticForm(stacked) - 0.5 * full.LogDet;
        double viaProduct = -0.5 * q.Dot(stacked, stacked) + 0.5 * n * q.LogDeterminant;
        Assert.Equal(dense, viaProduct, 10);
    }
}
