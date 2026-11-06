using System;
using System.Collections.Generic;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Quantiles of a predictive mixture over weighted posterior draws (D07's intervals): the quantile of the mixture
/// itself, never a quantile of each draw's own interval endpoints, which R's diagnostic plot takes. Weights need not be
/// normalized; a draw of weight 0 drops out. A mixture known only through its draws takes the shared weighted quantile,
/// <see cref="Maths.Samplers.Ensemble.DrawSample.Quantile"/>.
/// </summary>
public static class PredictiveMixture
{
    /// <summary>
    /// The <paramref name="probability"/> quantile of the Gaussian mixture Σₛ wₛ·N(μₛ, σₛ²): the root of
    /// Σₛ wₛΦ((q − μₛ)/σₛ) = p, which lies between the least and largest component quantiles. The tail beyond ½ is
    /// solved through the survival function, so neither tail loses its digits; bisection to adjacent doubles.
    /// </summary>
    public static double Quantile(ReadOnlySpan<double> means, ReadOnlySpan<double> sds, ReadOnlySpan<double> weights, double probability)
    {
        int s = means.Length;
        if (s == 0) throw new ArgumentException("A mixture needs a draw.", nameof(means));
        if (sds.Length != s || weights.Length != s) throw new ArgumentException("Every draw needs a mean, an sd and a weight.", nameof(sds));
        if (!(probability > 0.0 && probability < 1.0)) throw new ArgumentOutOfRangeException(nameof(probability), "The probability lies in (0, 1).");
        double total = Total(weights);
        double z = Normal.StandardQuantile(probability), lower = double.PositiveInfinity, upper = double.NegativeInfinity;
        for (int k = 0; k < s; k++)
        {
            if (!double.IsFinite(means[k])) throw new ArgumentException("The means must be finite.", nameof(means));
            if (!(sds[k] > 0.0) || double.IsPositiveInfinity(sds[k])) throw new ArgumentException("The sds must be positive and finite.", nameof(sds));
            if (weights[k] == 0.0) continue;
            double q = means[k] + sds[k] * z;
            lower = Math.Min(lower, q);
            upper = Math.Max(upper, q);
        }
        bool upperTail = probability > 0.5;
        double target = (upperTail ? 1.0 - probability : probability) * total;
        while (true)
        {
            double middle = lower + 0.5 * (upper - lower);
            if (middle <= lower || middle >= upper) return middle;
            // The mass below the midpoint (above it, in the upper tail) against the target.
            bool below = upperTail ? Mass(means, sds, weights, middle, true) > target : Mass(means, sds, weights, middle, false) < target;
            if (below) lower = middle;
            else upper = middle;
        }
    }

    /// <summary>
    /// The <paramref name="probability"/> quantile of the mixture Σₛ wₛ·Fₛ(· − μₛ) of shifted noises <paramref name="noises"/>: the root of
    /// Σₛ wₛFₛ(q − μₛ) = p, which lies between the least and largest component quantiles, the tail beyond ½ read through the survival
    /// function, by Newton with the mixture's density as its slope, kept in the bracket by bisection. A noise law's predictive (<see cref="IBassNoiseLaw.At"/>).
    /// </summary>
    public static double Quantile(ReadOnlySpan<double> means, IReadOnlyList<IContinuousDistribution> noises, ReadOnlySpan<double> weights, double probability)
    {
        int s = means.Length;
        if (s == 0) throw new ArgumentException("A mixture needs a draw.", nameof(means));
        ArgumentNullException.ThrowIfNull(noises);
        if (noises.Count != s || weights.Length != s) throw new ArgumentException("Every draw needs a mean, a noise and a weight.", nameof(noises));
        if (!(probability > 0.0 && probability < 1.0)) throw new ArgumentOutOfRangeException(nameof(probability), "The probability lies in (0, 1).");
        double total = Total(weights), lower = double.PositiveInfinity, upper = double.NegativeInfinity;
        for (int k = 0; k < s; k++)
        {
            if (!double.IsFinite(means[k])) throw new ArgumentException("The means must be finite.", nameof(means));
            if (weights[k] == 0.0) continue;
            double q = means[k] + noises[k].Quantile(probability);
            lower = Math.Min(lower, q);
            upper = Math.Max(upper, q);
        }
        bool upperTail = probability > 0.5;
        double target = (upperTail ? 1.0 - probability : probability) * total;
        // Newton on the mass against the target, its slope the mixture's density, kept inside the bracket by bisection: a noise whose CDF
        // is an integral (the generalized hyperbolic's) is read a few times per root, not once per halving.
        double x = lower + 0.5 * (upper - lower);
        for (int iteration = 0; iteration < 400; iteration++)
        {
            if (x <= lower || x >= upper) return x;
            double mass = 0.0, density = 0.0;
            for (int k = 0; k < s; k++)
            {
                if (weights[k] == 0.0) continue;
                mass += weights[k] * (upperTail ? noises[k].Survival(x - means[k]) : noises[k].Cdf(x - means[k]));
                density += weights[k] * Math.Exp(noises[k].LogPdf(x - means[k]));
            }
            double gap = upperTail ? target - mass : mass - target;   // rises with x
            if (gap == 0.0) return x;
            if (gap < 0.0) lower = x;
            else upper = x;
            double next = density > 0.0 ? x - gap / density : double.NaN;
            if (!(next > lower && next < upper)) next = lower + 0.5 * (upper - lower);
            if (Math.Abs(next - x) <= 1e-15 * Math.Max(1.0, Math.Abs(x))) return next;
            x = next;
        }
        return x;
    }

    private static double Mass(ReadOnlySpan<double> means, ReadOnlySpan<double> sds, ReadOnlySpan<double> weights, double q, bool above)
    {
        double mass = 0.0;
        for (int k = 0; k < means.Length; k++)
        {
            if (weights[k] == 0.0) continue;
            double z = (q - means[k]) / sds[k];
            mass += weights[k] * Math.Exp(Normal.StandardLogCdf(above ? -z : z));
        }
        return mass;
    }

    private static double Total(ReadOnlySpan<double> weights)
    {
        double total = 0.0;
        foreach (double w in weights)
        {
            if (!(w >= 0.0) || double.IsPositiveInfinity(w)) throw new ArgumentException("The weights must be finite and not negative.", nameof(weights));
            total += w;
        }
        if (!(total > 0.0)) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
        return total;
    }
}
