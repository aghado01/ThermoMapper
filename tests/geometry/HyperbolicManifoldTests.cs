using System;
using Maths.Distance.Geodesic;
using Xunit;

namespace Maths.Geometry.Tests;

/// <summary>
/// The hyperbolic models' distances through asinh and their maps without small-step guards, against closed forms: from
/// the Poincaré ball's origin the distance to x is 2·atanh‖x‖, and from the hyperboloid's apex the point
/// (sinh d·û, cosh d) sits at distance d. Distances far below arcosh's resolution of about 2e-8 come back to a few ε,
/// the Log map's norm agrees, and Exp moves a point by a step of 1e-14, which the Log map reads back.
/// </summary>
public sealed class HyperbolicManifoldTests
{
    private const double Epsilon = 2.220446049250313e-16;

    [Theory]
    [InlineData(1e-15)]
    [InlineData(3e-13)]
    [InlineData(1e-9)]
    [InlineData(0.3)]
    [InlineData(0.9)]
    public void Poincare_DistanceAndLogFromTheOrigin_MatchTheClosedForm(double radius)
    {
        const int n = 4;
        var ball = new PoincareBallManifold(n);
        double[] direction = Unit([0.4, -0.2, 0.7, 0.1]);
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = radius * direction[i];
        var origin = new double[n];

        double expected = 2.0 * Math.Atanh(radius);
        double tolerance = 8 * Epsilon * Math.Max(1.0, expected) / (1.0 - radius * radius);
        Assert.InRange(ball.Distance(origin, x), expected - tolerance, expected + tolerance);
        Assert.InRange(Poincare.Distance(origin, x), expected - tolerance, expected + tolerance);

        var log = new double[n];
        ball.LogMap(origin, x, log);
        Assert.InRange(ball.Norm(origin, log), expected - tolerance, expected + tolerance);
    }

    [Theory]
    [InlineData(1e-15)]
    [InlineData(3e-13)]
    [InlineData(1e-9)]
    [InlineData(0.7)]
    [InlineData(3.0)]
    public void Lorentz_DistanceAndLogFromTheApex_MatchTheClosedForm(double distance)
    {
        const int n = 5;
        var hyperboloid = new LorentzHyperboloidManifold(n);
        double[] direction = Unit([0.3, -0.6, 0.2, 0.7]);
        var apex = new double[n];
        apex[n - 1] = 1.0;
        var q = new double[n];
        for (int i = 0; i < n - 1; i++) q[i] = Math.Sinh(distance) * direction[i];
        q[n - 1] = Math.Cosh(distance);

        double tolerance = 16 * Epsilon * Math.Max(1.0, distance) * Math.Cosh(distance);
        Assert.InRange(hyperboloid.Distance(apex, q), distance - tolerance, distance + tolerance);
        Assert.InRange(Lorentz.Distance(apex, q), distance - tolerance, distance + tolerance);

        var log = new double[n];
        hyperboloid.LogMap(apex, q, log);
        Assert.InRange(hyperboloid.Norm(apex, log), distance - tolerance, distance + tolerance);
        for (int i = 0; i < n - 1; i++) Assert.InRange(log[i] - distance * direction[i], -tolerance, tolerance);
        Assert.InRange(log[n - 1], -tolerance, tolerance);
    }

    [Fact]
    public void TinySteps_MoveThePoint_AndTheLogMapReadsThemBack()
    {
        var ball = new PoincareBallManifold(3);
        double[] p = [0.2, -0.1, 0.3];
        double[] step = [1e-14, -2e-14, 0.5e-14];
        var moved = new double[3];
        ball.ExpMap(p, step, moved);
        // Off the origin the moved point's coordinates carry their own rounding, about ε·‖p‖ ≈ 1e-16 against a step
        // of 2e-14, so the step comes back to a few percent; the guard this replaces left the point unmoved.
        double length = ball.Norm(p, step);
        Assert.InRange(ball.Distance(p, moved), length * (1 - 5e-2), length * (1 + 5e-2));
        var back = new double[3];
        ball.LogMap(p, moved, back);
        for (int i = 0; i < 3; i++) Assert.InRange(back[i] - step[i], -1e-15, 1e-15);

        var hyperboloid = new LorentzHyperboloidManifold(3);
        double[] apex = [0.0, 0.0, 1.0];
        double[] tangent = [2e-14, -1e-14, 0.0];
        var lifted = new double[3];
        hyperboloid.ExpMap(apex, tangent, lifted);
        double norm = hyperboloid.Norm(apex, tangent);
        Assert.InRange(hyperboloid.Distance(apex, lifted), norm - 1e-16, norm + 1e-16);
        var returned = new double[3];
        hyperboloid.LogMap(apex, lifted, returned);
        for (int i = 0; i < 3; i++) Assert.InRange(returned[i] - tangent[i], -1e-16, 1e-16);
    }

    private static double[] Unit(double[] v)
    {
        double norm = 0.0;
        foreach (double e in v) norm += e * e;
        norm = Math.Sqrt(norm);
        var u = new double[v.Length];
        for (int i = 0; i < v.Length; i++) u[i] = v[i] / norm;
        return u;
    }
}
