using System;

namespace Graphs.Coupling;

/// <summary>Log-space arithmetic the kernels' log evaluators share.</summary>
internal static class LogMath
{
    /// <summary>
    /// ln(1 + x) with full relative precision for small |x| (x &gt; −1), where <c>Math.Log(1 + x)</c> loses every
    /// digit once 1 + x rounds to 1. Kahan's correction: when 1 + x differs from 1, the log's rounding error is
    /// divided out by <c>(1 + x) − 1</c>, which carries the same rounding.
    /// </summary>
    public static double Log1P(double x)
    {
        if (double.IsPositiveInfinity(x)) return double.PositiveInfinity;
        double u = 1.0 + x;
        if (u == 1.0) return x;
        return Math.Log(u) * x / (u - 1.0);
    }
}
