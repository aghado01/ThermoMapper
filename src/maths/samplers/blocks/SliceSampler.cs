using System;
using Maths.Rng;

namespace Maths.Samplers.Blocks;

/// <summary>
/// One univariate slice-sampling update (Neal 2003, the stepping-out and shrinkage procedures of his Figs. 3 and 5). It
/// draws the slice level <c>ln y = f(x) + ln(1 − U)</c>, places an interval of width w at a uniform offset around x,
/// steps each end out by w until it falls below the level (at most m steps in all, split between the ends at random),
/// then draws uniformly from the interval, shrinking it towards x on each draw below the level. The update leaves the
/// target invariant and always moves within the slice; w changes only its cost, since stepping out widens an interval
/// that is too narrow and shrinkage narrows one that is too wide. It needs the log density's values only.
/// </summary>
public static class SliceSampler
{
    /// <summary>Update <paramref name="x"/>, whose log density is <paramref name="logTarget"/>; both move to the new
    /// point. Returns the log-density evaluations it spent.</summary>
    /// <param name="width">The interval's initial width w.</param>
    /// <param name="maxSteps">The stepping-out limit m.</param>
    public static int Step(ref double x, ref double logTarget, double width, int maxSteps, Func<double, double> logDensity,
                           Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(logDensity);
        ArgumentNullException.ThrowIfNull(rng);
        if (!(width > 0.0) || double.IsInfinity(width)) throw new ArgumentOutOfRangeException(nameof(width), "The width must be positive and finite.");
        if (maxSteps < 1) throw new ArgumentOutOfRangeException(nameof(maxSteps));
        if (!double.IsFinite(logTarget)) throw new ArgumentException("The current point must have a finite log density.", nameof(logTarget));

        double level = logTarget + Math.Log(1.0 - rng.NextDouble());
        double left = x - width * rng.NextDouble(), right = left + width;
        int toLeft = (int)(maxSteps * rng.NextDouble()), toRight = maxSteps - 1 - toLeft;
        int evaluations = 0;
        while (toLeft > 0 && Above(left)) { left -= width; toLeft--; }
        while (toRight > 0 && Above(right)) { right += width; toRight--; }

        while (true)
        {
            double candidate = left + rng.NextDouble() * (right - left);
            double value = logDensity(candidate);
            evaluations++;
            if (value > level)
            {
                x = candidate;
                logTarget = value;
                return evaluations;
            }
            if (candidate < x) left = candidate;
            else right = candidate;
            // In exact arithmetic shrinkage ends at x, which is in the slice; stop there if rounding closes the interval.
            if (!(right - left > 1e-14 * (1.0 + Math.Abs(x)))) return evaluations;
        }

        bool Above(double point)
        {
            evaluations++;
            return logDensity(point) > level;
        }
    }
}
