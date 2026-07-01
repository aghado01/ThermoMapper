using System;

namespace Maths.Distance.Geodesic;

/// <summary>
/// Geodesic spherical distance on S^(d-1): the angular arc length in radians.
/// Inputs are treated as raw vectors and normalized internally; zero vectors
/// return π/2 by convention.
/// </summary>
public static class SphericalGeodesic
{
    private const double Epsilon = 1e-12;

    public static double Distance(double[] a, double[] b)
    {
        if (a is null || b is null || a.Length != b.Length)
            throw new ArgumentException("Invalid feature vectors for spherical geodesic distance.");

        return Distance((ReadOnlySpan<double>)a, b);
    }

    public static double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Spherical geodesic distance requires equal-length vectors.");

        double dot = 0.0;
        double normA2 = 0.0;
        double normB2 = 0.0;

        for (int i = 0; i < a.Length; i++)
        {
            double ai = a[i];
            double bi = b[i];
            dot += ai * bi;
            normA2 += ai * ai;
            normB2 += bi * bi;
        }

        if (normA2 <= Epsilon || normB2 <= Epsilon)
            return Math.PI / 2.0;

        double cos = dot / Math.Sqrt(normA2 * normB2);
        cos = Math.Clamp(cos, -1.0, 1.0);
        return Math.Acos(cos);
    }
}
