using System;

namespace Maths.Distributions;

/// <summary>A log tail probability in the coordinate its quantile is solved in, with its slope there.</summary>
internal interface ILogTail
{
    double At(double y, out double slope);
}

/// <summary>
/// Reads a quantile off the log tails, reducing an inverse to a root of the distribution's own <c>ln F</c> or
/// <c>ln S</c>: Newton in a coordinate where the density is log-concave, so the log tail is concave there and a step
/// from one side of the root stays on that side, with bisection inside a bracket when a step leaves it. The last step
/// comes back as a correction, which is second-order small, so a solve stops well above rounding and spends no evaluation
/// on the final digits.
/// </summary>
internal static class TailInversion
{
    private const int MaxIterations = 100;
    private const double StepTolerance = 1e-8;

    /// <summary>
    /// The y with <c>tail.At(y) = target</c> in [<paramref name="lo"/>, <paramref name="hi"/>], either possibly infinite, and
    /// the correction to add to it. <paramref name="scale"/> sets the step that ends the solve: a size of y, or 0 for
    /// relative to |y|, where y itself must keep relative digits.
    /// </summary>
    public static double Solve<T>(T tail, double target, double start, double lo, double hi, bool increasing, double scale, out double correction)
        where T : struct, ILogTail
    {
        double y = Math.Clamp(start, lo, hi), best = y, bestResidual = double.PositiveInfinity, expansion = 1.0;
        correction = 0.0;
        for (int i = 0; i < MaxIterations; i++)
        {
            double residual = tail.At(y, out double slope) - target;
            double step = -residual / slope;
            if (Math.Abs(residual) < bestResidual)
            {
                best = y;
                bestResidual = Math.Abs(residual);
                correction = double.IsFinite(step) ? step : 0.0;
            }
            double tolerance = StepTolerance * (scale > 0.0 ? scale : Math.Abs(y));
            if (residual == 0.0 || Math.Abs(step) <= tolerance) return best;

            if (residual > 0.0 == increasing) hi = y;
            else lo = y;
            if (hi - lo <= tolerance) return best;

            double next = y + step;
            if (!(next > lo && next < hi))
                next = double.IsFinite(lo) && double.IsFinite(hi) ? 0.5 * (lo + hi)
                     : double.IsInfinity(hi) ? y + (expansion *= 2.0)
                     : y - (expansion *= 2.0);
            if (next == y) return best;
            y = next;
        }
        return best;
    }
}
