using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 0 gates on its references (<see cref="ProjectionPursuitReferences"/>, the BPPR brief): each closed-form C against
/// <see cref="GradientQuadrature"/> of the gradient to 10⁻¹² of C's largest entry; each gradient against fourth-order central differences;
/// and the declared directions against C: C equals its projection onto their span and is positive definite on it, so the span is C's
/// range, rank included. Friedman's entries are also read against their stated values, and Cin(2π) against quadrature.
/// </summary>
public sealed class ProjectionPursuitReferenceTests(ITestOutputHelper output)
{
    private const double Tolerance = 1e-12;

    public static TheoryData<string> Names()
    {
        var names = new TheoryData<string>();
        foreach (ProjectionPursuitReference r in ProjectionPursuitReferences.All()) names.Add(r.Name);
        return names;
    }

    private static ProjectionPursuitReference Named(string name) => ProjectionPursuitReferences.All().Single(r => r.Name == name);

    [Theory]
    [MemberData(nameof(Names))]
    public void ClosedForm_MatchesQuadrature(string name)
    {
        ProjectionPursuitReference reference = Named(name);
        double[,] closed = reference.Function.C(reference.Law) ?? throw new InvalidOperationException($"{name} has no closed form under its law.");
        ReferenceBox box = reference.Law;
        // An input the gradient does not vary with needs one node: the inert ones, and Friedman's linear x₄, x₅.
        var rules = Enumerable.Range(0, box.Dimension)
            .Select(j => QuadratureRule.Uniform(box.Lower[j], box.Upper[j], Varies(reference.Function, j) ? 24 : 1, Varies(reference.Function, j) ? 2 : 1))
            .ToList();
        double error = RelativeError(closed, GradientQuadrature.C(reference.Function.Gradient, rules));
        output.WriteLine($"{name}: max |C − quadrature| / max |C| = {error:E2}");
        Assert.True(error <= Tolerance, $"{name}: {error:E2}");
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Gradient_MatchesCentralDifferences(string name)
    {
        ProjectionPursuitReference reference = Named(name);
        ReferenceFunction f = reference.Function;
        ReferenceBox box = reference.Law;
        var rng = new Xoshiro256PlusPlus(36);
        int p = f.Dimension;
        var g = new double[p];
        double worst = 0.0;
        for (int point = 0; point < 32; point++)
        {
            double[] x = [.. Enumerable.Range(0, p).Select(j => box.Lower[j] + (box.Upper[j] - box.Lower[j]) * rng.NextDouble())];
            f.Gradient(x, g);
            double scale = g.Max(Math.Abs);
            for (int j = 0; j < p; j++)
            {
                double step = Math.Pow(2.0, Math.Round(Math.Log2(1e-3 * (box.Upper[j] - box.Lower[j]))));
                double At(double offset)
                {
                    var y = (double[])x.Clone();
                    y[j] += offset;
                    return f.Value(y);
                }
                double difference = (-At(2 * step) + 8 * At(step) - 8 * At(-step) + At(-2 * step)) / (12 * step);
                worst = Math.Max(worst, Math.Abs(difference - g[j]) / scale);
            }
        }
        output.WriteLine($"{name}: worst |∂f − difference| / max |∇f| = {worst:E2}");
        Assert.True(worst <= 1e-8, $"{name}: {worst:E2}");
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Directions_SpanTheRangeOfC(string name)
    {
        ProjectionPursuitReference reference = Named(name);
        double[,] c = reference.Function.C(reference.Law)!;
        int p = c.GetLength(0), k = reference.Directions.Length;
        foreach (double[] w in reference.Directions) Assert.Equal(1.0, w.Sum(v => v * v), 14);
        double[][] q = Orthonormal(reference.Directions);

        var projector = new double[p, p];
        foreach (double[] v in q)
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) projector[i, j] += v[i] * v[j];
        double outside = RelativeError(c, Product(Product(projector, c), projector));

        // QᵀCQ by Cholesky: its smallest pivot² against its largest diagonal entry.
        var g = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++) g[a, b] += q[a][i] * c[i, j] * q[b][j];
        double smallest = double.PositiveInfinity, largest = Enumerable.Range(0, k).Max(a => g[a, a]);
        var l = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = g[i, j];
                for (int m = 0; m < j; m++) s -= l[i, m] * l[j, m];
                if (i == j)
                {
                    smallest = Math.Min(smallest, s);
                    l[i, i] = Math.Sqrt(Math.Max(s, 0.0));
                }
                else l[i, j] = s / l[j, j];
            }
        double conditioning = smallest / largest;
        output.WriteLine($"{name}: ‖C − PCP‖/‖C‖ = {outside:E2}; on the span, smallest pivot / largest diagonal = {conditioning:E2} over {k} direction(s)");
        Assert.True(outside <= 1e-13, $"{name}: C reaches outside the directions' span, {outside:E2}");
        Assert.True(conditioning > 1e-6, $"{name}: C is singular on the directions' span, {conditioning:E2}");
    }

    [Fact]
    public void Friedman_EntriesAndCin()
    {
        double[,] c = ProjectionPursuitReferences.Friedman().Function.C(ReferenceBox.Cube(6, 0.0, 1.0))!;
        Assert.Equal(100.0, c[3, 3]);
        Assert.Equal(50.0, c[3, 4]);
        Assert.Equal(25.0, c[4, 4]);
        Assert.Equal(200.0 / Math.PI, c[0, 3]);
        Assert.Equal(400.0 / 3.0, c[2, 2]);
        Assert.All(Enumerable.Range(0, 6), j => Assert.Equal(0.0, c[5, j]));
        for (int i = 0; i < 6; i++) Assert.Equal(0.0, c[i, 3] - 2.0 * c[i, 4], 12);   // C(e₄ − 2e₅) = 0

        // Cin(2π) by 32-node Gauss–Legendre on eight panels of (1 − cos t)/t, read as 2sin²(t/2)/t, the sum compensated; its terms'
        // own rounding leaves it a few ulps from the decimal series.
        GaussLegendre rule = GaussLegendre.Of(32);
        double sum = 0.0, carry = 0.0, width = 2.0 * Math.PI / 8;
        for (int panel = 0; panel < 8; panel++)
            for (int i = 0; i < rule.Nodes.Length; i++)
            {
                double t = (panel + 0.5) * width + 0.5 * width * rule.Nodes[i], s = Math.Sin(0.5 * t);
                double add = 0.5 * width * rule.Weights[i] * 2.0 * s * s / t, next = sum + add;
                carry += Math.Abs(sum) >= Math.Abs(add) ? sum - next + add : add - next + sum;
                sum = next;
            }
        double quadrature = sum + carry, cin = FriedmanReference.CinTwoPi();
        output.WriteLine($"Cin(2π) = {cin:R}, by quadrature {quadrature:R}; C₁₂ = {c[0, 1]:R}");
        Assert.True(Math.Abs(quadrature - cin) <= 4e-15, $"Cin(2π): {cin:R} against {quadrature:R}");
        Assert.Null(new FriedmanReference(6).C(ReferenceBox.Cube(6, 0.0, 2.0)));
    }

    // ── Helpers ─────────────────────────────────────────────────────────────────────────────────────────────────────

    private static bool Varies(ReferenceFunction f, int j) => f switch
    {
        RidgeSum sum => sum.Ridges.Any(r => r.Direction[j] != 0.0),
        FriedmanReference => j < 3,
        _ => true,
    };

    private static double[][] Orthonormal(double[][] directions)
    {
        var q = new List<double[]>();
        foreach (double[] d in directions)
        {
            var v = (double[])d.Clone();
            for (int pass = 0; pass < 2; pass++)
                foreach (double[] u in q)
                {
                    double dot = v.Zip(u, (a, b) => a * b).Sum();
                    for (int i = 0; i < v.Length; i++) v[i] -= dot * u[i];
                }
            double norm = Math.Sqrt(v.Sum(a => a * a));
            q.Add([.. v.Select(a => a / norm)]);
        }
        return [.. q];
    }

    private static double[,] Product(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), m = b.GetLength(1), inner = a.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < m; k++)
                for (int j = 0; j < inner; j++) c[i, k] += a[i, j] * b[j, k];
        return c;
    }

    private static double RelativeError(double[,] reference, double[,] other)
    {
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < reference.GetLength(0); i++)
            for (int k = 0; k < reference.GetLength(1); k++)
            {
                worst = Math.Max(worst, Math.Abs(reference[i, k] - other[i, k]));
                scale = Math.Max(scale, Math.Abs(reference[i, k]));
            }
        return worst / scale;
    }
}
