using System;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// Streaming count, weight, mean and centred sum of squares of one scalar sequence (West's weighted form of Welford's
/// update). A constant sequence leaves <see cref="SumSquaredDeviations"/> exactly zero, which is what lets R̂ tell a
/// constant functional from chains frozen at different values; running Σ and Σ² cannot, because their difference
/// cancels to rounding noise. Unit weights, the MCMC case, reproduce the unweighted update exactly; other weights pool
/// weighted particles. A readout statistic: <see cref="Merge"/> is associative and commutative.
/// </summary>
public struct RunningMoments
{
    /// <summary>Values folded in so far.</summary>
    public long Count { get; private set; }

    /// <summary>The total weight folded in: <see cref="Count"/> under unit weights.</summary>
    public double SumWeights { get; private set; }

    /// <summary>Running (weighted) mean (NaN once a non-finite value has been folded in).</summary>
    public double Mean { get; private set; }

    /// <summary>Σ wᵢ (xᵢ − mean)² over the values folded in.</summary>
    public double SumSquaredDeviations { get; private set; }

    /// <summary>Sample variance with the W − 1 denominator (n − 1 under unit weights, frequency weights otherwise);
    /// NaN below two values.</summary>
    public readonly double Variance => Count < 2 ? double.NaN : SumSquaredDeviations / (SumWeights - 1.0);

    /// <summary>Fold one value in, with weight <paramref name="weight"/> (1 by default).</summary>
    public void Add(double x, double weight = 1.0)
    {
        if (!(weight >= 0.0) || double.IsPositiveInfinity(weight))
            throw new ArgumentOutOfRangeException(nameof(weight), "A fold weight must be finite and non-negative.");
        Count++;
        SumWeights += weight;
        if (SumWeights == 0.0) return;   // only zero weights so far: nothing to average
        double delta = x - Mean;
        Mean += weight * delta / SumWeights;
        SumSquaredDeviations += weight * delta * (x - Mean);
    }

    /// <summary>Population variance (W denominator) — the plug-in variance of the values folded in; NaN when empty.</summary>
    public readonly double PopulationVariance => SumWeights > 0.0 ? SumSquaredDeviations / SumWeights : double.NaN;

    /// <summary>
    /// The moments of two disjoint sequences taken together (Chan, Golub &amp; LeVeque's parallel update, weighted).
    /// Associative and commutative, so pooling across chains does not depend on order.
    /// </summary>
    public static RunningMoments Merge(RunningMoments a, RunningMoments b)
    {
        if (a.SumWeights == 0.0) return Combine(b, a.Count);
        if (b.SumWeights == 0.0) return Combine(a, b.Count);
        double w = a.SumWeights + b.SumWeights;
        double delta = b.Mean - a.Mean;
        return new RunningMoments
        {
            Count = a.Count + b.Count,
            SumWeights = w,
            Mean = a.Mean + delta * b.SumWeights / w,
            SumSquaredDeviations = a.SumSquaredDeviations + b.SumSquaredDeviations + delta * delta * a.SumWeights * b.SumWeights / w,
        };
    }

    /// <summary>The moments of a whole sequence, under unit weights.</summary>
    public static RunningMoments Of(ReadOnlySpan<double> values)
    {
        var m = new RunningMoments();
        foreach (double x in values) m.Add(x);
        return m;
    }

    // A side with no weight adds only its fold count.
    private static RunningMoments Combine(RunningMoments m, long extraCount)
    {
        m.Count += extraCount;
        return m;
    }
}
