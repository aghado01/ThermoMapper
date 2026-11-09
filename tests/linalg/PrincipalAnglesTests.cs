using System;
using System.Linq;
using Maths.LinAlg;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The principal angles against constructed subspaces with known angles, from 0 and 1e-14 to π/2, turned by a random
/// rotation of the ambient space and random rotations within each span, so no basis is aligned with the axes. Every
/// angle comes back within about √(p·k)·ε radians, where the arccos of the cosines alone would leave the small ones
/// near 2e-8; the same holds across ranks, and the answer depends neither on the order of the bases nor on the
/// choice of basis within each span.
/// </summary>
public sealed class PrincipalAnglesTests
{
    private const double Epsilon = 2.220446049250313e-16;

    [Theory]
    [InlineData(new[] { 0.0, 1e-14, 1e-10, 1e-7 }, 12, 101)]
    [InlineData(new[] { 1e-3, 0.5, Math.PI / 4, 1.2 }, 12, 102)]
    [InlineData(new[] { Math.PI / 2 - 1e-9, Math.PI / 2, 1e-12, 0.3 }, 20, 103)]
    [InlineData(new[] { 1e-13, 1e-13, 1e-13 }, 9, 104)]
    public void KnownAngles_AreRecovered_ToAFewUlpsOfTheirSize(double[] theta, int p, int seed)
    {
        int k = theta.Length;
        double[] q = MatrixOps.RandomOrthonormal(p, p, seed);           // a rotation of the ambient space
        double[] a = Rotate(q, Within(Axes(p, k, theta, first: true), p, k, seed + 1000), p, k);
        double[] b = Rotate(q, Within(Axes(p, k, theta, first: false), p, k, seed + 2000), p, k);

        double[] angles = PrincipalAngles.Between(a, k, b, k, p);

        double[] expected = theta.OrderBy(t => t).ToArray();
        double tolerance = 20 * Math.Sqrt(p * k) * Epsilon;
        for (int i = 0; i < k; i++)
            Assert.True(Math.Abs(angles[i] - expected[i]) <= tolerance,
                $"angle {i}: {angles[i]:R} against {expected[i]:R} (tolerance {tolerance:G3})");
    }

    [Fact]
    public void MixedRanks_GiveTheSmallerRanksAngles()
    {
        const int p = 10;
        // A spans {e₀ cos a + e₇ sin a, e₁ cos b + e₈ sin b}; B spans {e₀, e₁, e₂, e₃}.
        const double angleA = 3e-11, angleB = 0.7;
        var a = new double[p * 2];
        a[0] = Math.Cos(angleA); a[7] = Math.Sin(angleA);
        a[p + 1] = Math.Cos(angleB); a[p + 8] = Math.Sin(angleB);
        var b = new double[p * 4];
        for (int j = 0; j < 4; j++) b[j * p + j] = 1.0;

        double[] ab = PrincipalAngles.Between(a, 2, b, 4, p);
        double[] ba = PrincipalAngles.Between(b, 4, a, 2, p);

        Assert.Equal(2, ab.Length);
        Assert.True(Math.Abs(ab[0] - angleA) <= 1e-15, $"{ab[0]:R}");
        Assert.True(Math.Abs(ab[1] - angleB) <= 1e-15, $"{ab[1]:R}");
        Assert.Equal(ab, ba);
    }

    [Fact]
    public void Angles_DependOnlyOnTheSpans_AndNotOnTheOrderOfTheBases()
    {
        const int p = 15, k = 3;
        double[] theta = [2e-9, 0.4, 1.1];
        double[] a = Axes(p, k, theta, first: true), b = Axes(p, k, theta, first: false);
        double[] reference = PrincipalAngles.Between(a, k, b, k, p);
        double[] swapped = PrincipalAngles.Between(b, k, a, k, p);
        double[] turned = PrincipalAngles.Between(Within(a, p, k, 7), k, Within(b, p, k, 8), k, p);

        double tolerance = 20 * Math.Sqrt(p * k) * Epsilon;
        for (int i = 0; i < k; i++)
        {
            Assert.True(Math.Abs(swapped[i] - reference[i]) <= tolerance, $"swapped {i}");
            Assert.True(Math.Abs(turned[i] - reference[i]) <= tolerance, $"turned {i}");
        }

        var a2 = new double[p, k];
        var b2 = new double[p, k];
        for (int j = 0; j < k; j++)
            for (int i = 0; i < p; i++) { a2[i, j] = a[i + j * p]; b2[i, j] = b[i + j * p]; }
        Assert.Equal(reference, PrincipalAngles.Between(a2, b2));
    }

    // Columns e_j (first) or e_j cos θ_j + e_{k+j} sin θ_j, p × k column-major: principal angles exactly θ.
    private static double[] Axes(int p, int k, double[] theta, bool first)
    {
        var basis = new double[p * k];
        for (int j = 0; j < k; j++)
        {
            basis[j * p + j] = first ? 1.0 : Math.Cos(theta[j]);
            if (!first) basis[j * p + k + j] = Math.Sin(theta[j]);
        }
        return basis;
    }

    // The basis times a random k × k rotation: the same span.
    private static double[] Within(double[] basis, int p, int k, int seed) =>
        Multiply(basis, p, k, MatrixOps.RandomOrthonormal(k, k, seed), k);

    // Q (p × p) times the basis.
    private static double[] Rotate(double[] q, double[] basis, int p, int k) => Multiply(q, p, p, basis, k);

    // A (m × n) times B (n × c), column-major.
    private static double[] Multiply(double[] a, int m, int n, double[] b, int c)
    {
        var result = new double[m * c];
        for (int j = 0; j < c; j++)
            for (int l = 0; l < n; l++)
                for (int i = 0; i < m; i++) result[i + j * m] += a[i + l * m] * b[l + j * n];
        return result;
    }
}
