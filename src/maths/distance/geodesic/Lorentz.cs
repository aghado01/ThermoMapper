using System;

namespace Maths.Distance.Geodesic;

/// <summary>
/// Geodesic distance on the Lorentz hyperboloid model of hyperbolic space.
/// Points are represented as vectors in R^N where the last coordinate (index N-1)
/// is the time-like component.
/// </summary>
public static class Lorentz
{
    /// <summary>
    /// Computes the Lorentzian inner product of two vectors in R^N:
    /// &lt;a, b&gt;_L = sum_{i=0}^{N-2} (a_i * b_i) - a_{N-1} * b_{N-1}.
    /// </summary>
    public static double InnerProduct(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Vectors must have the same length.");

        int len = a.Length;
        double sum = 0.0;
        for (int i = 0; i < len - 1; i++)
        {
            sum += a[i] * b[i];
        }
        sum -= a[len - 1] * b[len - 1];
        return sum;
    }

    public static double Distance(double[] a, double[] b)
        => Distance((ReadOnlySpan<double>)a, b);

    public static double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Invalid feature vectors for Lorentz distance.");

        double ip = InnerProduct(a, b);
        double arg = -ip;

        // Clamp to prevent numerical precision issues leading to arg < 1.0 (which would make arcosh NaN)
        if (arg < 1.0)
            return 0.0;

        // arcosh(x) = ln(x + sqrt(x^2 - 1))
        return Math.Log(arg + Math.Sqrt(arg * arg - 1.0));
    }

    /// <summary>
    /// Projects a point in R^N onto the upper sheet of the hyperboloid
    /// by preserving the spatial coordinates and setting the time-like coordinate (last coordinate)
    /// to sqrt(1 + sum_{i=0}^{N-2} x_i^2).
    /// </summary>
    public static void ProjectToHyperboloid(Span<double> point)
    {
        if (point.Length < 2)
            throw new ArgumentException("Point must have at least 2 components (1 spatial + 1 time-like).");

        int len = point.Length;
        double spatialNormSq = 0.0;
        for (int i = 0; i < len - 1; i++)
        {
            spatialNormSq += point[i] * point[i];
        }

        point[len - 1] = Math.Sqrt(1.0 + spatialNormSq);
    }
}
