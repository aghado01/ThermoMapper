using System;

namespace Maths.Distance.Geodesic;

public static class Wasserstein1
{
    public static double Distance(double[] p, double[] q)
    {
        if (p is null || q is null || p.Length != q.Length)
            throw new ArgumentException("Invalid PMF vectors for Wasserstein1 distance.");

        return Distance((ReadOnlySpan<double>)p, q);
    }

    public static double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
    {
        if (p.Length != q.Length)
            throw new ArgumentException("Wasserstein1 distance requires equal-length vectors.");

        int len = p.Length;
        if (len == 0) return 0.0;

        double cumP = 0.0;
        double cumQ = 0.0;
        double sum = 0.0;

        for (int i = 0; i < len; i++)
        {
            cumP += p[i];
            cumQ += q[i];
            sum += Math.Abs(cumP - cumQ);
        }

        return sum;
    }
}
