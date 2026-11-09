using System;

namespace Maths.Distance.Geodesic;

/// <summary>
/// Geodesic spherical distance on S^(d-1): the angular arc length in radians.
/// Inputs are treated as raw vectors and normalized internally; zero vectors
/// return π/2 by convention. The angle is Kahan's 2·atan2(‖â − b̂‖, ‖â + b̂‖), accurate at every angle to a few
/// ε, where the arccos of the cosine resolves only about √(2ε) ≈ 2e-8 near zero.
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

        double invA = 1.0 / Math.Sqrt(normA2), invB = 1.0 / Math.Sqrt(normB2);
        double difference = 0.0, sum = 0.0;
        for (int i = 0; i < a.Length; i++)
        {
            double ai = a[i] * invA, bi = b[i] * invB;
            difference += (ai - bi) * (ai - bi);
            sum += (ai + bi) * (ai + bi);
        }
        return 2.0 * Math.Atan2(Math.Sqrt(difference), Math.Sqrt(sum));
    }
}
