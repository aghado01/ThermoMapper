using System;

namespace Clustering.Graphical.CumulativeField;

/// <summary>Log-space running sums for the cumulative-energy ladder.</summary>
internal static class LogLadder
{
    /// <summary>ln(e^a + e^b), stable for any a and b; −∞ is a zero term.</summary>
    public static double LogAddExp(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double max = Math.Max(a, b);
        return max + Math.Log(1.0 + Math.Exp(-Math.Abs(a - b)));
    }
}
