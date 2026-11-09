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

    /// <summary>
    /// The geodesic distance between two points on the hyperboloid, 2·asinh(√⟨a − b, a − b⟩_L / 2). On the hyperboloid
    /// ⟨a − b, a − b⟩_L = 2(cosh d − 1), so this equals arcosh(−⟨a, b⟩_L), but it reads the distance from the difference
    /// of the points and keeps a small distance's digits, where arcosh of a number near 1 resolves only about
    /// √(2ε) ≈ 2e-8.
    /// </summary>
    public static double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Invalid feature vectors for Lorentz distance.");

        return 2.0 * Math.Asinh(0.5 * Math.Sqrt(Math.Max(0.0, DifferenceSquare(a, b))));
    }

    /// <summary>⟨a − b, a − b⟩_L, the Lorentzian square of the difference, which is 2(cosh d − 1) on the hyperboloid.</summary>
    public static double DifferenceSquare(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Vectors must have the same length.");

        int len = a.Length;
        double sum = 0.0;
        for (int i = 0; i < len - 1; i++)
        {
            double d = a[i] - b[i];
            sum += d * d;
        }
        double t = a[len - 1] - b[len - 1];
        return sum - t * t;
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
