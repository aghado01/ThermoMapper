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
/// (Gilks and Wild 1992) with fixed design points, as in transformed density rejection (Hörmann 1995). A reach guessed far past
/// the point where h falls by 1, or a step that overshoots it, is searched back from rather than left where h is astronomically
/// low, so any reach is valid. It reads the density alone, with no normalizer or CDF, and accepts over half its proposals on any
/// interval of any log-concave density, 86% or more on the families here at any reach, tending to 1 far in a tail. Three uniforms per
/// proposal.
/// </summary>
internal static class LogConcaveEnvelope
{
    private const int MaxDesignSteps = 32;
    private const int MaxGuardSteps = 64;
    private const int MaxAttempts = 100_000;

    // A design point where h has fallen past this is astronomical: the far tangent there is so steep that rounding in its cut, which
    // grows with h's size and the tangent's slope, can swamp the hat or sink it below h. The guarded search steps back from below it.
    private const double Astronomical = 1048576.0;

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

    /// <summary>The hat's area for <see cref="Draw"/>'s arguments: ∫e^h over the interval divided by it is the share of proposals accepted.</summary>
    internal static double HatArea<T>(T h, double left, double right, double leftReach, double rightReach) where T : struct, IConcaveLogRatio
        => Side.Of(h, -1.0, left, leftReach).Area + Side.Of(h, 1.0, right, rightReach).Area;

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
            double t = DesignPoint(h, sign, end, reach, out double height);
            double farRate = Math.Max(nearRate, -sign * h.Slope(sign * t));
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

        // Newton for h = −1 along t, as it has always run. Close is enough: the hat's area is flat around the point. Only where it runs
        // out of steps short of that, as from a reach far past the root, where Newton steps back about one unit of h's growth at a
        // time, does the guarded search run, so every point this search already found, and every draw made from it, is unchanged bit
        // for bit.
        private static double DesignPoint<T>(T h, double sign, double end, double reach, out double height) where T : struct, IConcaveLogRatio
        {
            double t = Math.Min(reach, end);
            for (int i = 0; i < MaxDesignSteps; i++)
            {
                double value = h.Value(sign * t), slope = sign * h.Slope(sign * t);
                if (Math.Abs(value + 1.0) <= 0.25 || (t == end && value > -1.0))
                {
                    height = value;
                    return t;
                }
                double next = slope < 0.0 ? t - (value + 1.0) / slope : 2.0 * t;
                if (!(next > 0.0) || double.IsInfinity(next)) next = value < -1.0 ? 0.5 * t : 2.0 * t;
                t = Math.Min(next, end);
            }
            height = h.Value(sign * t);
            return Math.Abs(height + 1.0) <= 0.25 || (t == end && height > -1.0) ? t : GuardedPoint(h, sign, end, reach, out height);
        }

        // The search again from the reach: Newton's steps where h is above −2²⁰, and where it is below, or has overflowed, the guard's
        // step (Guarded), which does not count against Newton's. A reach that is no positive number starts at 1, h's own unit.
        private static double GuardedPoint<T>(T h, double sign, double end, double reach, out double height) where T : struct, IConcaveLogRatio
        {
            double t = reach > 0.0 ? Math.Min(Math.Min(reach, end), double.MaxValue) : Math.Min(1.0, end);
            double inside = 0.0;   // the farthest t read with h above −1: the root lies past it
            for (int steps = 0, guarded = 0; ; )
            {
                double value = h.Value(sign * t), slope = sign * h.Slope(sign * t);
                if (!(value >= -Astronomical))
                {
                    if (guarded++ == MaxGuardSteps)
                        throw new ArithmeticException($"The envelope's design point found no h above −2²⁰ in {MaxGuardSteps} steps.");
                    t = Guarded(t, value, slope, inside);
                    continue;
                }
                if (steps++ == MaxDesignSteps || Math.Abs(value + 1.0) <= 0.25 || (t == end && value > -1.0))
                {
                    height = value;
                    return t;
                }
                if (value > -1.0) inside = Math.Max(inside, t);
                double next = slope < 0.0 ? t - (value + 1.0) / slope : 2.0 * t;
                if (!(next > 0.0) || double.IsInfinity(next)) next = value < -1.0 ? 0.5 * t : 2.0 * t;
                t = Math.Min(next, end);
            }
        }

        // A step back from a t where h is astronomically low. Newton on ln(−h) against ln t lands at t·|h|^(−1/k), k = t·h′/h, which is
        // exact where h is a power of t and close where it is an exponential; concavity with h(0) = 0 makes k ≥ 1, so the step stops
        // no nearer than t/|h|, where h ≥ −1. Where h has overflowed, or that step would not land between the last point inside and
        // t, the step halves the distance from the inside point, in ln t beyond 1 (h's own unit) and in t within it.
        private static double Guarded(double t, double value, double slope, double inside)
        {
            if (double.IsFinite(value) && slope < 0.0 && double.IsFinite(slope))
            {
                double next = t * Math.Exp(-Math.Log(-value) * (value / slope) / t);
                if (next > inside && next < t) return next;
            }
            return t > 1.0 ? Math.Sqrt(Math.Max(inside, 1.0) * t) : 0.5 * (inside + t);
        }
    }
}
