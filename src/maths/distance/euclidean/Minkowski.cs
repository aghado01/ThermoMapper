using System;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;

namespace Maths.Distance.Euclidean;

public static class Minkowski
{
    private const double MinExponent = 1.0;
    private const double MaxExponent = 2.0;

    public static double Distance(double[] a, double[] b, double p)
    {
        if (a is null || b is null || a.Length != b.Length)
            throw new ArgumentException("Invalid feature vectors for Minkowski distance.");

        p = ClampExponent(p);
        if (p == 1.0) return ManhattanDistance(a, b);
        if (p == 2.0) return TensorPrimitives.Distance<double>(a, b);

        return Distance((ReadOnlySpan<double>)a, b, p);
    }

    public static double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b, double p)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Minkowski distance requires equal-length vectors.");

        p = ClampExponent(p);
        if (p == 1.0) return ManhattanDistance(a, b);
        if (p == 2.0) return TensorPrimitives.Distance<double>(a, b);

        int len = a.Length;
        if (len == 0) return 0.0;

        double maxDiff = 0.0;
        for (int i = 0; i < len; i++)
        {
            double d = Math.Abs(a[i] - b[i]);
            if (d > maxDiff) maxDiff = d;
        }

        if (maxDiff <= 0.0)
            return 0.0;

        double sum = 0.0;
        for (int i = 0; i < len; i++)
        {
            double scaled = Math.Abs(a[i] - b[i]) / maxDiff;
            sum += Math.Pow(scaled, p);
        }

        return maxDiff * Math.Pow(sum, 1.0 / p);
    }

    public static double ManhattanDistance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        int len = a.Length;
        double sum = 0.0;
        for (int i = 0; i < len; i++)
            sum += Math.Abs(a[i] - b[i]);

        return sum;
    }

    public static double ManhattanDistance(double[] a, double[] b)
        => ManhattanDistance((ReadOnlySpan<double>)a, b);

    public static double ClampExponent(double p)
    {
        if (p < MinExponent) return MinExponent;
        if (p > MaxExponent) return MaxExponent;
        return p;
    }
}
