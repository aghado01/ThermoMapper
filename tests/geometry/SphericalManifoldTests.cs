using System;
using Maths.Distance.Geodesic;
using Xunit;

namespace Maths.Geometry.Tests;

/// <summary>
/// The sphere's angles through atan2: a distance, a Log map and an Exp map that resolve angles far below the arccos
/// floor of about 2e-8, and agree with each other there; the antipode stays the cut locus, where Log returns zero.
/// Near the antipode the Log map's direction comes from a residual of size sin θ, so the vector holds ε·θ/sin θ.
/// </summary>
public sealed class SphericalManifoldTests
{
    [Theory]
    [InlineData(1e-15)]
    [InlineData(3e-13)]
    [InlineData(1e-9)]
    [InlineData(0.7)]
    [InlineData(Math.PI / 2)]
    [InlineData(Math.PI - 1e-6)]
    public void KnownAngles_AreResolved_InDistanceAndLog(double angle)
    {
        const int n = 5;
        var sphere = new SphericalManifold(n);
        double[] p = Unit(n, [0.3, -0.5, 0.2, 0.7, -0.1]);
        double[] t = Unit(n, Horizontal(p, [0.1, 0.4, -0.6, 0.2, 0.5]));    // unit tangent at p
        var q = new double[n];
        for (int i = 0; i < n; i++) q[i] = Math.Cos(angle) * p[i] + Math.Sin(angle) * t[i];

        double tolerance = 8 * 2.220446049250313e-16 * Math.Max(1.0, angle);
        Assert.InRange(sphere.Distance(p, q), angle - tolerance, angle + tolerance);
        Assert.InRange(SphericalGeodesic.Distance(p, q), angle - tolerance, angle + tolerance);

        var log = new double[n];
        sphere.LogMap(p, q, log);
        Assert.InRange(sphere.Norm(p, log), angle - tolerance, angle + tolerance);
        double direction = tolerance * Math.Max(1.0, angle / Math.Sin(angle));
        for (int i = 0; i < n; i++) Assert.InRange(log[i] - angle * t[i], -direction, direction);
    }

    [Fact]
    public void TinySteps_MoveThePoint_AndTheLogMapReadsThemBack()
    {
        const int n = 4;
        var sphere = new SphericalManifold(n);
        double[] p = Unit(n, [0.5, 0.5, -0.5, 0.5]);
        double[] step = Horizontal(p, [1e-14, -3e-14, 2e-14, 0.0]);

        var moved = new double[n];
        sphere.ExpMap(p, step, moved);
        double length = Math.Sqrt(Dot(step, step));
        Assert.InRange(sphere.Distance(p, moved), length - 2e-16, length + 2e-16);
        var back = new double[n];
        sphere.LogMap(p, moved, back);
        for (int i = 0; i < n; i++) Assert.InRange(back[i] - step[i], -2e-16, 2e-16);
    }

    [Fact]
    public void Antipode_IsTheCutLocus_WhereLogIsZero()
    {
        const int n = 3;
        var sphere = new SphericalManifold(n);
        double[] p = Unit(n, [1.0, 2.0, -2.0]);
        double[] q = [-p[0], -p[1], -p[2]];
        Assert.InRange(sphere.Distance(p, q), Math.PI - 1e-15, Math.PI + 1e-15);
        var log = new double[n];
        sphere.LogMap(p, q, log);
        Assert.All(log, v => Assert.Equal(0.0, v));
    }

    private static double[] Unit(int n, double[] v)
    {
        double norm = Math.Sqrt(Dot(v, v));
        var u = new double[n];
        for (int i = 0; i < n; i++) u[i] = v[i] / norm;
        return u;
    }

    private static double[] Horizontal(double[] p, double[] v)
    {
        double d = Dot(p, v);
        var h = new double[v.Length];
        for (int i = 0; i < v.Length; i++) h[i] = v[i] - d * p[i];
        return h;
    }

    private static double Dot(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }
}
