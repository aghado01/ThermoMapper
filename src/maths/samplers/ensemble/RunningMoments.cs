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

    /// <summary>The moments of a whole sequence.</summary>
    public static RunningMoments Of(ReadOnlySpan<double> values)
    {
        var m = new RunningMoments();
        foreach (double x in values) m.Add(x);
        return m;
    }
}
