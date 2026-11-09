using System;

namespace Maths.Regression.Conformal;

/// <summary>A witnessed change of a grid set's decision, ordered along y by the caller's grid.</summary>
public readonly record struct ConformalBoundary(double Accepted, double Rejected);

internal static class ConformalBisection
{
    internal static void CheckTolerance(double tolerance)
    {
        if (!(tolerance > 0.0 && double.IsFinite(tolerance))) throw new ArgumentOutOfRangeException(nameof(tolerance));
    }

    // The endpoints have opposite decisions. Only narrow that observed crossing; other components can lie between grid values.
    internal static GridInterval Narrow(double lo, double hi, bool lowDecision, Func<double, bool> decide, double tolerance)
    {
        while (hi - lo > tolerance)
        {
            double mid = lo < 0.0 && hi > 0.0 ? (lo + hi) * 0.5 : lo + (hi - lo) * 0.5;
            if (!(mid > lo && mid < hi))
                throw new InvalidOperationException("The requested response tolerance is smaller than the spacing of adjacent doubles at the crossing.");
            if (decide(mid) == lowDecision) lo = mid;
            else hi = mid;
        }
        return new(lo, hi);
    }
}
