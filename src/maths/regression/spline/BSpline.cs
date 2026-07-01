using System;

namespace Maths.Regression.Spline;

public static class BSpline
{
    /// <summary>
    /// Evaluates all non-zero B-spline basis functions of the given degree at parameter t
    /// using the Cox–de Boor recurrence. Returns a length-(knots.Length - degree - 1) array;
    /// only the entries near the active knot span are non-zero.
    /// </summary>
    public static double[] EvaluateBasis(double t, double[] knots, int degree = 3)
    {
        int n = knots.Length - degree - 1;
        double[] result = new double[n];

        // Locate the active knot span, clamped to [degree, n-1] so the degree+1 non-zero basis
        // functions land inside the result array even at the right endpoint t = 1.
        int span = degree;
        while (span < n - 1 && t >= knots[span + 1])
            span++;

        // Cox–de Boor / Piegl & Tiller A2.2 on a local length-(degree+1) array.
        double[] nLocal = new double[degree + 1];
        double[] left  = new double[degree + 1];
        double[] right = new double[degree + 1];
        nLocal[0] = 1.0;

        for (int j = 1; j <= degree; j++)
        {
            left[j]  = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            double saved = 0.0;

            for (int r = 0; r < j; r++)
            {
                double denom = right[r + 1] + left[j - r];
                // Coincident knots (multiplicity) give denom = 0; the term is 0 by convention.
                double temp = denom > 0.0 ? nLocal[r] / denom : 0.0;
                nLocal[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            nLocal[j] = saved;
        }

        // The degree+1 non-zero values occupy global indices [span-degree .. span].
        for (int r = 0; r <= degree; r++)
            result[span - degree + r] = nLocal[r];
        return result;
    }

    /// <summary>
    /// Builds a clamped knot vector for a B-spline parameterized on [0, 1].
    /// The first and last knot are repeated degree+1 times to clamp the curve
    /// to its end control points.
    /// </summary>
    public static double[] MakeClampedKnots(double[] interiorKnots, int degree = 3)
    {
        int nInterior = interiorKnots.Length;
        int nControl  = nInterior + degree + 1;
        double[] knots = new double[nControl + degree + 1];

        for (int i = 0; i <= degree; i++) knots[i] = 0.0;
        for (int i = 0; i < nInterior; i++) knots[degree + 1 + i] = interiorKnots[i];
        for (int i = 0; i <= degree; i++) knots[degree + 1 + nInterior + i] = 1.0;
        return knots;
    }
}
