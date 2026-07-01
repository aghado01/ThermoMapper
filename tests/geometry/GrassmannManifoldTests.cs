using System;
using Maths.Geometry;
using Maths.LinAlg;
using Xunit;

namespace Maths.Geometry.Tests;

public sealed class GrassmannManifoldTests
{
    // A nearby subspace: perturb Y's basis and re-orthonormalize (stays well inside the
    // injectivity radius, so Log is well-defined).
    private static double[] Perturb(double[] y, int n, int r, double eps, int seed)
    {
        var rng = new Random(seed);
        var z = (double[])y.Clone();
        for (int i = 0; i < z.Length; i++) z[i] += eps * (rng.NextDouble() * 2.0 - 1.0);
        MatrixOps.Orthonormalize(z, n, r);
        return z;
    }

    [Fact]
    public void Distance_SameSubspace_IsZero_AndSymmetric()
    {
        int n = 8, r = 3;
        var g = new GrassmannManifold(n, r);
        var y = MatrixOps.RandomOrthonormal(n, r, 11);
        var z = Perturb(y, n, r, 0.05, 22);

        // Self-distance floor is ~1e-8, not 0: principal angles are arccos(σ), and arccos has
        // infinite slope at 1, so machine-eps in the Gram eigenvalues amplifies to √(2·eps).
        Assert.InRange(g.Distance(y, y), 0.0, 1e-5);
        Assert.Equal(g.Distance(y, z), g.Distance(z, y), 9);
        Assert.True(g.Distance(y, z) > 1e-4);
    }

    [Fact]
    public void LogNorm_Equals_Distance()
    {
        int n = 10, r = 4;
        var g = new GrassmannManifold(n, r);
        var y = MatrixOps.RandomOrthonormal(n, r, 7);
        var z = Perturb(y, n, r, 0.1, 9);

        var t = new double[n * r];
        g.LogMap(y, z, t);
        Assert.Equal(g.Distance(y, z), g.Norm(y, t), 6);
    }

    [Fact]
    public void ExpLog_RoundTrip_RecoversSubspace()
    {
        int n = 10, r = 4;
        var g = new GrassmannManifold(n, r);
        var y = MatrixOps.RandomOrthonormal(n, r, 3);
        var z = Perturb(y, n, r, 0.1, 5);

        var t = new double[n * r];
        g.LogMap(y, z, t);
        var zExp = new double[n * r];
        g.ExpMap(y, t, zExp);

        Assert.InRange(g.Distance(zExp, z), 0.0, 1e-5);   // same subspace
    }

    [Fact]
    public void Exp_PreservesOrthonormality()
    {
        int n = 9, r = 3;
        var g = new GrassmannManifold(n, r);
        var y = MatrixOps.RandomOrthonormal(n, r, 4);
        var z = Perturb(y, n, r, 0.2, 6);

        var t = new double[n * r];
        g.LogMap(y, z, t);
        var e = new double[n * r];
        g.ExpMap(y, t, e);

        for (int a = 0; a < r; a++)
            for (int b = 0; b < r; b++)
            {
                double dot = 0.0;
                for (int k = 0; k < n; k++) dot += e[a * n + k] * e[b * n + k];
                double expected = a == b ? 1.0 : 0.0;
                Assert.InRange(dot, expected - 1e-9, expected + 1e-9);
            }
    }

    [Fact]
    public void Distance_IsGaugeInvariant()
    {
        int n = 8, r = 3;
        var g = new GrassmannManifold(n, r);
        var y = MatrixOps.RandomOrthonormal(n, r, 2);
        var z = MatrixOps.RandomOrthonormal(n, r, 8);
        var q = MatrixOps.RandomOrthonormal(r, r, 15);    // r×r orthogonal gauge

        var zq = MatMulColMajor(z, n, r, q, r);           // Z·Q — same span, different basis
        Assert.Equal(g.Distance(y, z), g.Distance(y, zq), 9);
    }

    private static double[] MatMulColMajor(double[] a, int m, int k, double[] b, int bc)
    {
        var c = new double[m * bc];
        for (int j = 0; j < bc; j++)
            for (int l = 0; l < k; l++)
                for (int i = 0; i < m; i++)
                    c[i + j * m] += a[i + l * m] * b[l + j * k];
        return c;
    }
}
