using System;

namespace Maths.Distributions;

/// <summary>
/// The empirical distribution of a sample: an atom of mass 1/n at each value, repeats kept, built with <see cref="Of"/>.
/// Its interval probabilities are counts and its interval moments exact sums over the atoms in (lower, upper],
/// compensated, at every power.
/// </summary>
public sealed class Empirical : IIntervalMoments
{
    private readonly double[] _sorted;

    private Empirical(double[] sorted) => _sorted = sorted;

    /// <summary>The empirical distribution of <paramref name="sample"/>: nonempty and finite.</summary>
    public static Empirical Of(ReadOnlySpan<double> sample)
    {
        if (sample.IsEmpty) throw new ArgumentException("An empirical distribution needs a value.", nameof(sample));
        double[] sorted = sample.ToArray();
        foreach (double v in sorted)
            if (!double.IsFinite(v)) throw new ArgumentException("An empirical distribution's values must be finite.", nameof(sample));
        Array.Sort(sorted);
        return new Empirical(sorted);
    }

    /// <summary>The sample size n.</summary>
    public int Count => _sorted.Length;

    /// <summary>The values, ascending.</summary>
    public ReadOnlySpan<double> Values => _sorted;

    /// <summary>The smallest value.</summary>
    public double Lower => _sorted[0];

    /// <summary>The largest value.</summary>
    public double Upper => _sorted[^1];

    /// <summary><c>P(lower &lt; X ≤ upper)</c>, the share of values in (lower, upper].</summary>
    public double IntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return 0.0;
        (int first, int end) = Range(lower, upper);
        return (double)(end - first) / _sorted.Length;
    }

    /// <inheritdoc/>
    public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
    {
        moments.Clear();
        if (Tails.IsEmpty(lower, upper)) return 0.0;
        (int first, int end) = Range(lower, upper);
        Span<double> compensation = moments.Length <= 32 ? stackalloc double[moments.Length] : new double[moments.Length];
        compensation.Clear();
        for (int i = first; i < end; i++)
        {
            double d = _sorted[i] - centre, power = 1.0;
            for (int k = 0; k < moments.Length; k++)
            {
                // Neumaier's compensated sum.
                double sum = moments[k] + power;
                compensation[k] += Math.Abs(moments[k]) >= Math.Abs(power) ? moments[k] - sum + power : power - sum + moments[k];
                moments[k] = sum;
                power *= d;
            }
        }
        double mass = Math.Exp(-logScale) / _sorted.Length;
        for (int k = 0; k < moments.Length; k++) moments[k] = (moments[k] + compensation[k]) * mass;
        return 0.0;
    }

    // The atoms in (lower, upper]: indices [first, end).
    private (int First, int End) Range(double lower, double upper)
    {
        int first = UpperBound(lower), end = UpperBound(upper);
        return (first, Math.Max(first, end));
    }

    // The first index whose value exceeds x.
    private int UpperBound(double x)
    {
        int lo = 0, hi = _sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (_sorted[mid] <= x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
}
