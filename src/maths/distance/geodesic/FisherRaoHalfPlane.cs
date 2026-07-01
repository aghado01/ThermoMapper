using System;
using System.Runtime.CompilerServices;

namespace Maths.Distance.Geodesic;

public static class FisherRaoHalfPlane
{
    private const double Sqrt2 = 1.4142135623730950488;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Distance(double[] p, double[] q)
        => Distance((ReadOnlySpan<double>)p, q);

    public static double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
    {
        if (p.Length != 2 || q.Length != 2)
            throw new ArgumentException(
                "FisherRaoHalfPlane requires 2-element [mu, log_sigma] vectors.");

        double dmu = p[0] - q[0];
        double sig1 = Math.Exp(p[1]);
        double sig2 = Math.Exp(q[1]);
        double dsig = sig1 - sig2;

        double z = 1.0 + (dmu * dmu + 2.0 * dsig * dsig) / (4.0 * sig1 * sig2);
        return Sqrt2 * Math.Acosh(z);
    }
}
