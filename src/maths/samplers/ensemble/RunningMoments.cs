using System;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// Streaming count, mean and centred sum of squares of one scalar sequence (Welford's update). A constant sequence
/// leaves <see cref="SumSquaredDeviations"/> exactly zero, which is what lets R̂ tell a constant functional from
/// chains frozen at different values; running Σ and Σ² cannot, because their difference cancels to rounding noise.
/// </summary>
public struct RunningMoments
{
    /// <summary>Values folded in so far.</summary>
    public long Count { get; private set; }

    /// <summary>Running mean (NaN once a non-finite value has been folded in).</summary>
    public double Mean { get; private set; }

    /// <summary>Σ (x − mean)² over the values folded in.</summary>
    public double SumSquaredDeviations { get; private set; }

    /// <summary>Sample variance with the n − 1 denominator; NaN below two values.</summary>
    public readonly double Variance => Count < 2 ? double.NaN : SumSquaredDeviations / (Count - 1);

    /// <summary>Fold one value in.</summary>
    public void Add(double x)
    {
        Count++;
        double delta = x - Mean;
        Mean += delta / Count;
        SumSquaredDeviations += delta * (x - Mean);
    }

    /// <summary>Population variance (n denominator) — the plug-in variance of the values folded in; NaN when empty.</summary>
    public readonly double PopulationVariance => Count < 1 ? double.NaN : SumSquaredDeviations / Count;

    /// <summary>
    /// The moments of two disjoint sequences taken together (Chan, Golub &amp; LeVeque's parallel update). Associative
    /// and commutative, so pooling across chains does not depend on order.
    /// </summary>
    public static RunningMoments Merge(RunningMoments a, RunningMoments b)
    {
        if (a.Count == 0) return b;
        if (b.Count == 0) return a;
        long n = a.Count + b.Count;
        double delta = b.Mean - a.Mean;
        return new RunningMoments
        {
            Count = n,
            Mean = a.Mean + delta * b.Count / n,
            SumSquaredDeviations = a.SumSquaredDeviations + b.SumSquaredDeviations + delta * delta * a.Count * b.Count / n,
        };
    }

    /// <summary>The moments of a whole sequence.</summary>
    public static RunningMoments Of(ReadOnlySpan<double> values)
    {
        var m = new RunningMoments();
        foreach (double x in values) m.Add(x);
        return m;
    }
}
