using System;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The standardized shifted Legendre polynomials: orthonormality on [0, 1] by Gauss–Legendre quadrature to degree 40; the
/// closed forms at low degree; ψ_α vanishing at the Gauss–Legendre nodes, independent zeros, to degree 64; the endpoint values
/// and slopes; the slope against Bonnet's identity in a seeded sweep; the lanes and the batch against the scalar bit for bit;
/// and the expansion about a centre against the values.
/// </summary>
public sealed class ShiftedLegendreTests
{
    [Fact]
    public void Orthonormal_ByQuadrature_ToDegreeForty()
    {
        const int top = 40;
        GaussLegendre rule = GaussLegendre.Of(top + 2);   // exact to degree 2·42 − 1 ≥ 80
        var values = new double[rule.Count][];
        for (int k = 0; k < rule.Count; k++)
        {
            values[k] = new double[top + 1];
            ShiftedLegendre.Values(0.5 * (rule.Nodes[k] + 1.0), values[k]);
        }
        double worst = 0.0;
        for (int a = 0; a <= top; a++)
            for (int b = 0; b <= a; b++)
            {
                double sum = 0.0;
                for (int k = 0; k < rule.Count; k++) sum += 0.5 * rule.Weights[k] * values[k][a] * values[k][b];
                worst = Math.Max(worst, Math.Abs(sum - (a == b ? 1.0 : 0.0)));
            }
        Assert.True(worst < 1e-12, $"worst Gram entry error {worst:E2}");
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.13)]
    [InlineData(0.5)]
    [InlineData(0.77)]
    [InlineData(1.0)]
    public void LowDegrees_AreTheClosedForms(double x)
    {
        double u = 2.0 * x - 1.0;
        Assert.Equal(1.0, ShiftedLegendre.Value(0, x));
        Assert.Equal(Math.Sqrt(3.0) * u, ShiftedLegendre.Value(1, x), 1e-15);
        Assert.Equal(Math.Sqrt(5.0) * (3 * u * u - 1) / 2, ShiftedLegendre.Value(2, x), 1e-14);
        Assert.Equal(Math.Sqrt(7.0) * (5 * u * u * u - 3 * u) / 2, ShiftedLegendre.Value(3, x), 1e-14);
        Assert.Equal(Math.Sqrt(9.0) * (35 * Math.Pow(u, 4) - 30 * u * u + 3) / 8, ShiftedLegendre.Value(4, x), 1e-13);
        Assert.Equal(0.0, ShiftedLegendre.Derivative(0, x));
        Assert.Equal(2.0 * Math.Sqrt(3.0), ShiftedLegendre.Derivative(1, x), 1e-15);
        Assert.Equal(2.0 * Math.Sqrt(5.0) * 3 * u, ShiftedLegendre.Derivative(2, x), 1e-14);
        Assert.Equal(2.0 * Math.Sqrt(7.0) * (15 * u * u - 3) / 2, ShiftedLegendre.Derivative(3, x), 1e-13);
    }

    [Fact]
    public void ZerosAreTheGaussLegendreNodes_ToDegreeSixtyFour()
    {
        for (int n = 1; n <= 64; n++)
        {
            GaussLegendre rule = GaussLegendre.Of(n);
            double scale = Math.Sqrt(2.0 * n + 1.0);
            for (int k = 0; k < n; k++)
            {
                double x = 0.5 * (rule.Nodes[k] + 1.0);
                // A root's value is its slope times the node's rounding, which is a few ulps of x.
                double bound = 64 * Math.Abs(ShiftedLegendre.Derivative(n, x)) * Math.Max(x, 1.0 - x) * 2.2e-16 + 1e-15 * scale;
                Assert.True(Math.Abs(ShiftedLegendre.Value(n, x)) < bound, $"degree {n}, node {k}: {ShiftedLegendre.Value(n, x):E2} against {bound:E2}");
            }
        }
    }

    [Fact]
    public void EndpointsHoldTheirValuesAndSlopes_ToDegreeSixty()
    {
        for (int a = 0; a <= 60; a++)
        {
            double root = Math.Sqrt(2.0 * a + 1.0), slope = 2.0 * root * a * (a + 1) / 2.0;
            Assert.Equal(root, ShiftedLegendre.Value(a, 1.0), 1e-12 * root);
            Assert.Equal(a % 2 == 0 ? root : -root, ShiftedLegendre.Value(a, 0.0), 1e-12 * root);
            Assert.Equal(slope, ShiftedLegendre.Derivative(a, 1.0), 1e-12 * Math.Max(1.0, slope));
            Assert.Equal(a % 2 == 1 ? slope : -slope, ShiftedLegendre.Derivative(a, 0.0), 1e-12 * Math.Max(1.0, slope));
        }
    }

    [Fact]
    public void SeededSweep_SlopesMeetBonnetsIdentity_AndValuesStayBounded()
    {
        // (1 − u²)P′ₙ = n(P_{n−1} − uPₙ), away from the ends, where the identity's two sides are both of order n.
        var rng = new Xoshiro256PlusPlus(20261107);
        double worst = 0.0;
        for (int c = 0; c < 4000; c++)
        {
            int n = 1 + (int)(rng.NextDouble() * 60);
            double x = 0.02 + 0.96 * rng.NextDouble(), u = 2.0 * x - 1.0;
            double pn = ShiftedLegendre.Value(n, x) / Math.Sqrt(2.0 * n + 1.0);
            double pm = ShiftedLegendre.Value(n - 1, x) / Math.Sqrt(2.0 * n - 1.0);
            double dn = ShiftedLegendre.Derivative(n, x) / (2.0 * Math.Sqrt(2.0 * n + 1.0));
            Assert.True(Math.Abs(pn) <= 1.0 + 1e-13, $"|P_{n}({u})| = {Math.Abs(pn)}");
            worst = Math.Max(worst, Math.Abs((1.0 - u * u) * dn - n * (pm - u * pn)) / n);
        }
        Assert.True(worst < 1e-12, $"worst identity residual {worst:E2}");
    }

    [Fact]
    public void LanesAndBatch_MatchTheScalar_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(11);
        var batch = new double[36];
        for (int n = 1; n <= 17; n++)
        {
            var x = new double[n];
            for (int i = 0; i < n; i++) x[i] = i == 0 ? 0.0 : i == 1 ? 1.0 : rng.NextDouble();
            for (int degree = 0; degree <= 35; degree++)
            {
                var start = new double[n];
                for (int i = 0; i < n; i++) start[i] = 0.5 + rng.NextDouble();
                var column = (double[])start.Clone();
                ShiftedLegendre.Multiply(degree, x, column);
                for (int i = 0; i < n; i++) Assert.Equal(start[i] * ShiftedLegendre.Value(degree, x[i]), column[i]);
            }
            ShiftedLegendre.Values(x[n - 1], batch);
            for (int degree = 0; degree < batch.Length; degree++) Assert.Equal(ShiftedLegendre.Value(degree, x[n - 1]), batch[degree]);
        }
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(0.0)]
    [InlineData(0.3)]
    public void ExpansionAboutACentre_ReproducesTheValues(double centre)
    {
        for (int degree = 0; degree <= 12; degree++)
        {
            var c = new double[degree + 1];
            ShiftedLegendre.Coefficients(degree, centre, c);
            foreach (double x in new[] { 0.0, 0.11, 0.5, 0.62, 1.0 })
            {
                // The sum cancels terms of size Σ|cₖ||x − c|ᵏ, so that, times a few ulps, bounds its error.
                double sum = 0.0, size = 0.0, power = 1.0;
                for (int k = 0; k <= degree; k++, power *= x - centre)
                {
                    sum += c[k] * power;
                    size += Math.Abs(c[k] * power);
                }
                Assert.Equal(ShiftedLegendre.Value(degree, x), sum, 1e-14 * size + 1e-14);
            }
        }
    }
}
