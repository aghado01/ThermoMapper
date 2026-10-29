using System;
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
