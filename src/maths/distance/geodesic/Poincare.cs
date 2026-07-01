using System;

namespace Maths.Distance.Geodesic;

public static class Poincare
{
    /// <summary>
    /// Margin from the unit-ball boundary used by <see cref="ClampToBall"/>.
    /// </summary>
    public const double BoundaryMargin = 1e-5;

    public static double Distance(double[] a, double[] b)
        => Distance((ReadOnlySpan<double>)a, b);

    public static double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Invalid feature vectors for Poincare distance.");

        int dim = a.Length;
        double limitSq = (1.0 - BoundaryMargin) * (1.0 - BoundaryMargin);

        double normA2 = 0.0, normB2 = 0.0;
        for (int i = 0; i < dim; i++)
        {
            normA2 += a[i] * a[i];
            normB2 += b[i] * b[i];
        }

        double scaleA = normA2 > limitSq ? Math.Sqrt(limitSq / normA2) : 1.0;
        double scaleB = normB2 > limitSq ? Math.Sqrt(limitSq / normB2) : 1.0;

        double effA2 = 0.0, effB2 = 0.0, diff2 = 0.0;
        for (int i = 0; i < dim; i++)
        {
            double ai = a[i] * scaleA;
            double bi = b[i] * scaleB;
            effA2 += ai * ai;
            effB2 += bi * bi;
            double d = ai - bi;
            diff2 += d * d;
        }

        double denom = (1.0 - effA2) * (1.0 - effB2);
        if (denom <= 0.0) return 0.0;

        double arg = 1.0 + 2.0 * diff2 / denom;
        if (arg < 1.0) return 0.0;

        double sqrtTerm = Math.Sqrt(arg * arg - 1.0);
        return Math.Log(arg + sqrtTerm);
    }

    public static void ClampToBall(Span<double> point)
    {
        double n2 = 0.0;
        for (int i = 0; i < point.Length; i++) n2 += point[i] * point[i];

        double limit = 1.0 - BoundaryMargin;
        double limitSq = limit * limit;
        if (n2 <= limitSq) return;

        double scale = limit / Math.Sqrt(n2);
        for (int i = 0; i < point.Length; i++) point[i] *= scale;
    }
}
