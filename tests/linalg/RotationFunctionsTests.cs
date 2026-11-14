using System;
using Maths.LinAlg;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// <see cref="MatrixFunctions.SkewExponential"/> and <see cref="MatrixFunctions.RotationLogarithm"/> against rotations in known
/// orthogonal planes: L = Σ θ_j (u_j v_jᵀ − v_j u_jᵀ) has exp(L) = I + Σ [(cos θ_j − 1)(u_j u_jᵀ + v_j v_jᵀ) + sin θ_j (u_j v_jᵀ − v_j u_jᵀ)]
/// in closed form, and for |θ_j| &lt; π it is exp(L)'s principal logarithm. Chosen cases and a seeded sweep over sizes and angles, tiny
/// and near-half-turn angles included.
/// </summary>
public sealed class RotationFunctionsTests
{
    [Theory]
    [InlineData(1e-12)]
    [InlineData(1e-6)]
    [InlineData(0.5)]
    [InlineData(2.0)]
    [InlineData(3.1)]
    public void PlaneRotation_BothWays(double theta)
    {
        var omega = new double[,] { { 0.0, -theta }, { theta, 0.0 } };
        double[,] e = MatrixFunctions.SkewExponential(omega);
        AssertClose(new[,] { { Math.Cos(theta), -Math.Sin(theta) }, { Math.Sin(theta), Math.Cos(theta) } }, e, 1e-15);

        // The angle comes back to rounding in absolute terms, so a 1e-12 turn keeps about four digits.
        double[,] log = MatrixFunctions.RotationLogarithm(e);
        Assert.Equal(theta, log[1, 0], 1e-15 * Math.Max(1.0, theta / (Math.PI - theta)));
        Assert.Equal(-log[1, 0], log[0, 1]);
    }

    [Fact]
    public void Sweep_AgainstKnownPlanes()
    {
        var rng = new Random(20261113);
        double worstExp = 0.0, worstLog = 0.0, worstOrth = 0.0;
        for (int trial = 0; trial < 300; trial++)
        {
            int m = 2 + rng.Next(7);
            int planes = 1 + rng.Next(m / 2);
            var theta = new double[planes];
            for (int j = 0; j < planes; j++)
                theta[j] = rng.Next(5) switch
                {
                    0 => 1e-9 * (rng.NextDouble() - 0.5),
                    1 => (Math.PI - 1e-3) * (rng.Next(2) == 0 ? 1 : -1),
                    _ => (2.0 * rng.NextDouble() - 1.0) * (Math.PI - 0.05),
                };
            double[] basis = MatrixOps.RandomOrthonormal(m, 2 * planes, rng.Next(1, int.MaxValue));
            (double[,] l, double[,] rotation) = KnownPlanes(basis, m, theta);

            double[,] e = MatrixFunctions.SkewExponential(l);
            worstExp = Math.Max(worstExp, MaxAbsDifference(rotation, e));
            worstOrth = Math.Max(worstOrth, MaxAbsDifference(Identity(m), Multiply(Transpose(e), e)));

            double[,] log = MatrixFunctions.RotationLogarithm(rotation);
            // Near a half-turn the logarithm's condition grows like 1/(π − θ); here π − θ ≥ 1e-3.
            worstLog = Math.Max(worstLog, MaxAbsDifference(l, log));
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) Assert.Equal(-log[j, i], log[i, j]);
        }
        Assert.InRange(worstExp, 0.0, 2e-14);
        Assert.InRange(worstOrth, 0.0, 2e-14);
        Assert.InRange(worstLog, 0.0, 1e-11);
    }

    [Fact]
    public void Logarithm_RefusesHalfTurnsAndReflections()
    {
        var halfTurn = new double[,] { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
        var reflection = new double[,] { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, -1.0 } };
        Assert.False(MatrixFunctions.TryRotationLogarithm(halfTurn, out _));
        Assert.False(MatrixFunctions.TryRotationLogarithm(reflection, out _));
        Assert.Throws<ArgumentException>(() => MatrixFunctions.RotationLogarithm(reflection));
    }

    [Fact]
    public void Exponential_ReadsTheSkewPart()
    {
        var omega = new double[,] { { 0.3, -0.7 }, { 0.7, -0.1 } };
        var skew = new double[,] { { 0.0, -0.7 }, { 0.7, 0.0 } };
        AssertClose(MatrixFunctions.SkewExponential(skew), MatrixFunctions.SkewExponential(omega), 0.0);
    }

    private static (double[,] L, double[,] Rotation) KnownPlanes(double[] basis, int m, double[] theta)
    {
        var l = new double[m, m];
        double[,] rotation = Identity(m);
        for (int j = 0; j < theta.Length; j++)
        {
            int u = 2 * j * m, v = (2 * j + 1) * m;
            double c = Math.Cos(theta[j]) - 1.0, s = Math.Sin(theta[j]);
            for (int a = 0; a < m; a++)
                for (int b = 0; b < m; b++)
                {
                    double uv = basis[u + a] * basis[v + b] - basis[v + a] * basis[u + b];
                    l[a, b] += theta[j] * uv;
                    rotation[a, b] += c * (basis[u + a] * basis[u + b] + basis[v + a] * basis[v + b]) + s * uv;
                }
        }
        return (l, rotation);
    }

    private static double[,] Identity(int n)
    {
        var i = new double[n, n];
        for (int k = 0; k < n; k++) i[k, k] = 1.0;
        return i;
    }

    private static double[,] Transpose(double[,] a)
    {
        int n = a.GetLength(0), m = a.GetLength(1);
        var t = new double[m, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) t[j, i] = a[i, j];
        return t;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), m = b.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                double s = 0.0;
                for (int l = 0; l < k; l++) s += a[i, l] * b[l, j];
                c[i, j] = s;
            }
        return c;
    }

    private static double MaxAbsDifference(double[,] a, double[,] b)
    {
        double worst = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]));
        return worst;
    }

    private static void AssertClose(double[,] expected, double[,] actual, double tolerance) =>
        Assert.InRange(MaxAbsDifference(expected, actual), 0.0, tolerance);
}
