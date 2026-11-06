using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>A density relative to its largest value, ln f(x) − ln f(mode) ≤ 0, and −∞ off its support.</summary>
internal interface IRelativeLogDensity
{
    double Value(double x);
}

/// <summary>
/// Ratio-of-uniforms (Kinderman and Monahan 1977), with a shift by the mode (Dagpunar 1989; Lehner 1989): (U, V) uniform on the rectangle
/// [u₋, u₊] × (0, 1] and X = U/V + m, accepted when V ≤ √(f(X)/f(m)). The rectangle is the minimal one about the region
/// {(u, v): 0 &lt; v ≤ √(f(u/v + m)/f(m))} when u₋ and u₊ are the least and greatest values of (x − m)·√(f(x)/f(m)); then the accepted X
/// has density f exactly. Where f is T₋½-concave, 1/√f convex, the region is convex and at least half the rectangle (Leydold 2000), so
/// the acceptance is at least ½ at any parameter; two uniforms per proposal and one log. A family with the rectangle's sides in closed
/// form passes them; <see cref="Bounds"/> finds them by search for one without. Exact and untruncated; a truncated draw is the tangent
/// hat's (<see cref="LogConcaveEnvelope"/>).
/// </summary>
internal static class RatioOfUniforms
{
    private const int MaxAttempts = 100_000;

    /// <summary>One draw of <paramref name="h"/>'s law: X = U/V + <paramref name="shift"/>, U on [<paramref name="uMinus"/>,
    /// <paramref name="uPlus"/>], V on (0, 1], accepted when 2·ln V ≤ h(X).</summary>
    public static double Draw<T>(Xoshiro256PlusPlus rng, T h, double shift, double uMinus, double uPlus) where T : struct, IRelativeLogDensity
    {
        double width = uPlus - uMinus;
        if (!(width > 0.0) || !double.IsFinite(width))
            throw new ArithmeticException($"The ratio-of-uniforms rectangle [{uMinus}, {uPlus}] is empty or unbounded.");
        for (int attempt = 0; attempt < MaxAttempts; attempt++)
        {
            double u = uMinus + width * rng.NextDouble(), v = 1.0 - rng.NextDouble();
            double x = u / v + shift;
            if (2.0 * Math.Log(v) <= h.Value(x)) return x;
        }
        throw new ArithmeticException($"Ratio-of-uniforms accepted nothing in {MaxAttempts} attempts.");
    }

    /// <summary>
    /// The rectangle's sides about <paramref name="mode"/> on the support [<paramref name="lower"/>, <paramref name="upper"/>]: the least and
    /// greatest of (x − m)·e^(h(x)/2), each by golden section on its side after its bracket is grown from the mode, which finds them where
    /// the product is unimodal on each side, as it is for a T₋½-concave f. Widened by 10⁻⁹ of themselves, so rounding in the search leaves
    /// the region inside.
    /// </summary>
    public static (double Minus, double Plus) Bounds<T>(T h, double mode, double lower, double upper) where T : struct, IRelativeLogDensity
    {
        double G(double x) => (x - mode) * Math.Exp(0.5 * h.Value(x));
        double Side(double end, double direction)
        {
            if (!(direction > 0 ? end > mode : end < mode)) return 0.0;
            double scale = Math.Max(Math.Abs(mode), 1.0) * 1e-3, near = mode, far = mode + direction * scale;
            for (int k = 0; k < 4000 && (direction > 0 ? far < end : far > end) && Math.Abs(G(far)) >= Math.Abs(G(near)); k++)
            {
                near = far;
                far = mode + direction * (scale *= 2.0);
            }
            if (direction > 0 ? far > end : far < end) far = end;
            double lo = Math.Min(mode, far), hi = Math.Max(mode, far);
            const double golden = 0.6180339887498949;
            for (int k = 0; k < 200; k++)
            {
                double c = hi - golden * (hi - lo), d = lo + golden * (hi - lo);
                if (Math.Abs(G(c)) > Math.Abs(G(d))) hi = d;
                else lo = c;
            }
            return G(0.5 * (lo + hi));
        }
        double minus = Side(lower, -1.0), plus = Side(upper, 1.0);
        return (minus * (1.0 + 1e-9), plus * (1.0 + 1e-9));
    }
}
