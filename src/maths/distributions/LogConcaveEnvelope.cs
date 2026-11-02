using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// A density relative to its largest value on an interval, in a coordinate where it is log-concave:
/// <c>h(δ) = ln f(anchor + δ) − ln f(anchor)</c>, concave, with h(0) = 0 its maximum on the interval. Differences from the
/// anchor keep their digits on a narrow interval or far in a tail, where the two log densities would cancel.
/// </summary>
internal interface IConcaveLogRatio
{
    double Value(double delta);
    double Slope(double delta);
}

/// <summary>
/// Rejection under a hat of tangents, the sampler behind every truncated draw that has no closed form. Each side of the
/// peak gets two tangents of <see cref="IConcaveLogRatio"/>: one at the peak, flat when the peak is interior, and one
/// where h first falls by 1, or at the interval's end if it never does. That is the tangent hat of adaptive rejection
/// (Gilks and Wild 1992) with fixed design points, as in transformed density rejection (Hörmann 1995). It reads the
/// density alone, with no normalizer or CDF, and accepts over half its proposals on any interval of any log-concave density,
/// 87% or more on the families here, tending to 1 far in a tail. Three uniforms per proposal.
/// </summary>
internal static class LogConcaveEnvelope
{
    private const int MaxDesignSteps = 32;
    private const int MaxAttempts = 100_000;

    // An exponential this many scales long has no cut to speak of: e⁻⁴⁰ is below the rounding of 1.
    private const double FullExposure = 40.0;

    /// <summary>
    /// One δ on [−<paramref name="left"/>, <paramref name="right"/>], either span possibly infinite, with density ∝ e^h(δ).
    /// A reach is a first guess at the point on that side where h = −1; any guess is valid, and a good one saves a step.
    /// </summary>
    public static double Draw<T>(Xoshiro256PlusPlus rng, T h, double left, double right, double leftReach, double rightReach)
        where T : struct, IConcaveLogRatio
    {
        Side below = Side.Of(h, -1.0, left, leftReach), above = Side.Of(h, 1.0, right, rightReach);
        double total = below.Area + above.Area;
        if (double.IsNaN(total) || double.IsInfinity(total))
            throw new ArithmeticException($"The envelope over [-{left}, {right}] has no finite area.");
        if (total == 0.0) return 0.0;   // bounds that collapse in this coordinate

        for (int attempt = 0; attempt < MaxAttempts; attempt++)
        {
            double pick = rng.NextDouble() * total, position = rng.NextDouble();
            double hat, delta = pick < above.Area ? above.Propose(pick, position, out hat) : -below.Propose(pick - above.Area, position, out hat);
            // ln u ≤ u − 1, so a uniform below 1 + gap is accepted without its logarithm.
            double gap = h.Value(delta) - hat, u = 1.0 - rng.NextDouble();
            if (u <= 1.0 + gap || Math.Log(u) <= gap) return delta;
        }
        throw new ArithmeticException($"The envelope over [-{left}, {right}] accepted nothing in {MaxAttempts} attempts.");
    }

    /// <summary>The t ≥ 0 where −rt − κt²/2 = −1: the reach of a slope r and curvature κ at the anchor.</summary>
    public static double Reach(double rate, double curvature) => 2.0 / (rate + Math.Sqrt(rate * rate + 2.0 * curvature));

    /// <summary>A draw from the exponential of <paramref name="rate"/> cut to [0, width], from a uniform; uniform at rate 0.</summary>
    public static double TruncatedExponential(double uniform, double rate, double width)
    {
        if (rate == 0.0) return uniform * width;
        double exposure = rate * width;
        return exposure > FullExposure ? -Math.Log(1.0 - uniform) / rate
                                       : -SpecialFunctions.Log1p(uniform * SpecialFunctions.Expm1(-exposure)) / rate;
    }

    /// <summary><c>ln(b/a)</c> for 0 ≤ a ≤ b, either end possibly 0 or +∞ where it is +∞, keeping its digits for b close to a.</summary>
    public static double LogRatio(double a, double b)
    {
        if (!(b > a)) return 0.0;
        double r = (b - a) / a;
        return r < 1.0 ? SpecialFunctions.Log1p(r) : Math.Log(b) - Math.Log(a);
    }

    /// <summary><c>x·e^δ</c>, added to x where δ is small, so x keeps its digits; two factors where it is not, so an in-range result is reached.</summary>
    public static double Scaled(double x, double delta)
    {
        if (Math.Abs(delta) < 0.5) return x + x * SpecialFunctions.Expm1(delta);
        double half = Math.Exp(0.5 * delta);
        return x * half * half;
    }

    // One side of the hat in t = |δ|: the peak's tangent −r₀t on [0, k], the far tangent v − r₁(t − k) on [k, end]. A tangent
    // lies above a concave h wherever it is cut, so rounding in the cut k costs efficiency and never validity.
    private readonly struct Side
    {
        private readonly double _nearRate, _kink, _farRate, _farHeight, _farWidth, _nearArea, _farArea;

        public double Area => _nearArea + _farArea;

        private Side(double nearRate, double kink, double farRate, double farHeight, double farWidth)
        {
            _nearRate = nearRate;
            _kink = kink;
            _farRate = farRate;
            _farHeight = farHeight;
            _farWidth = farWidth;
            _nearArea = Segment(nearRate, kink);
            _farArea = farWidth > 0.0 ? Math.Exp(farHeight) * Segment(farRate, farWidth) : 0.0;
        }

        public static Side Of<T>(T h, double sign, double end, double reach) where T : struct, IConcaveLogRatio
        {
            if (!(end > 0.0)) return default;
            double nearRate = Math.Max(0.0, -sign * h.Slope(0.0));
            double t = DesignPoint(h, sign, end, reach);
            double height = h.Value(sign * t), farRate = Math.Max(nearRate, -sign * h.Slope(sign * t));
            if (!(farRate > nearRate)) return new Side(nearRate, end, 0.0, 0.0, 0.0);   // h is straight here: one tangent
            double kink = Math.Clamp((height + farRate * t) / (farRate - nearRate), 0.0, t);
            return new Side(nearRate, kink, farRate, height - farRate * (kink - t), end - kink);
        }

        public double Propose(double pick, double position, out double hat)
        {
            if (pick < _nearArea)
            {
                double t = TruncatedExponential(position, _nearRate, _kink);
                hat = -_nearRate * t;
                return t;
            }
            double beyond = TruncatedExponential(position, _farRate, _farWidth);
            hat = _farHeight - _farRate * beyond;
            return _kink + beyond;
        }

        private static double Segment(double rate, double width)
        {
            if (rate == 0.0) return width;
            double exposure = rate * width;
            return exposure > FullExposure ? 1.0 / rate : -SpecialFunctions.Expm1(-exposure) / rate;
        }

        // Newton for h = −1 along t. Close is enough: the hat's area is flat around the point.
        private static double DesignPoint<T>(T h, double sign, double end, double reach) where T : struct, IConcaveLogRatio
        {
            double t = Math.Min(reach, end);
            for (int i = 0; i < MaxDesignSteps; i++)
            {
                double value = h.Value(sign * t), slope = sign * h.Slope(sign * t);
                if (Math.Abs(value + 1.0) <= 0.25 || (t == end && value > -1.0)) break;
                double next = slope < 0.0 ? t - (value + 1.0) / slope : 2.0 * t;
                if (!(next > 0.0) || double.IsInfinity(next)) next = value < -1.0 ? 0.5 * t : 2.0 * t;
                t = Math.Min(next, end);
            }
            return t;
        }
    }
}
