using System;

namespace Maths.LinAlg;

/// <summary>The pivot accuracy test used by a dense symmetric Jacobi solve.</summary>
public enum JacobiStoppingRule
{
    /// <summary>Absolute off-diagonal tolerance; retained while callers migrate.</summary>
    Absolute,
    /// <summary>|a_pq| ≤ tol √|a_pp a_qq|, including zero pivots when a diagonal is zero.</summary>
    ElementWiseRelative,
}

internal static class JacobiPivot
{
    internal static void Validate(JacobiStoppingRule rule, double tol)
    {
        if (rule != JacobiStoppingRule.Absolute && rule != JacobiStoppingRule.ElementWiseRelative)
            throw new ArgumentOutOfRangeException(nameof(rule));
        if (rule == JacobiStoppingRule.ElementWiseRelative && (!double.IsFinite(tol) || tol < 0 || tol >= 1))
            throw new ArgumentOutOfRangeException(nameof(tol), "Relative tolerance must be finite and in [0, 1).");
    }

    internal static bool Passes(double app, double aqq, double apq, double tol)
    {
        if (!double.IsFinite(app) || !double.IsFinite(aqq) || !double.IsFinite(apq)) return false;
        // Square roots separately avoid overflow/underflow in app*aqq. Multiply tol into the larger root first so a
        // representable threshold is not lost to an intermediate underflow on an uneven pair.
        double a = Math.Sqrt(Math.Abs(app)), b = Math.Sqrt(Math.Abs(aqq));
        return Math.Abs(apq) <= Math.Min(a, b) * (tol * Math.Max(a, b));
    }

    internal static bool TryTangent(double app, double aqq, double apq, out double t)
    {
        // Normalize the local pair, not the matrix: neither aqq-app nor 2*apq may overflow. No global floor is used.
        double scale = Math.Max(Math.Abs(apq), Math.Max(Math.Abs(app), Math.Abs(aqq)));
        double delta = 0.5 * (aqq / scale - app / scale), pivot = apq / scale;
        t = delta == 0 ? 1.0 : pivot / (delta + Math.CopySign(Math.Sqrt(delta * delta + pivot * pivot), delta));
        // A failed accuracy test with an unrepresentable rotation remains in the working matrix and cannot report convergence.
        return double.IsFinite(t) && t != 0 && double.IsFinite(scale) && scale > 0;
    }
}
