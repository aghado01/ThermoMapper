using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.GPrior;

/// <summary>Huang and Wand's conditional IW rates rₖ = ν/aₖ, with rₖ ~ Gamma(½, rate 1/(νAₖ²)).</summary>
public sealed class ResolvedOutputRates
{
    private readonly double[] _rates;

    internal ResolvedOutputRates(double nu, double[] scales)
    {
        Nu = nu;
        _rates = new double[scales.Length];
        for (int k = 0; k < scales.Length; k++) _rates[k] = 1.0 / (nu * scales[k] * scales[k]);
    }

    public double Nu { get; }
    public int Count => _rates.Length;

    internal double[] Centre()
    {
        var log = new double[Count];
        for (int k = 0; k < Count; k++) log[k] = Math.Log(0.5 / _rates[k]);
        return log;
    }

    internal double[] Sample(Xoshiro256PlusPlus rng)
    {
        var log = new double[Count];
        for (int k = 0; k < Count; k++)
            do log[k] = Math.Log(Gamma.ShapeRate(0.5, _rates[k]).Sample(rng));
            while (!double.IsFinite(log[k]));
        return log;
    }

    internal double LogDensity(ReadOnlySpan<double> log)
    {
        if (log.Length != Count) throw new ArgumentException("One conditional covariance rate per output is required.", nameof(log));
        double result = 0;
        for (int k = 0; k < Count; k++)
        {
            double r = Math.Exp(log[k]);
            if (!(r > 0 && double.IsFinite(r))) return double.NegativeInfinity;
            result += 0.5 * (Math.Log(_rates[k]) + log[k]) - SpecialFunctions.LogGamma(0.5) - _rates[k] * r;
        }
        return result;
    }

    internal InverseWishart PriorAt(ReadOnlySpan<double> log)
    {
        if (log.Length != Count) throw new ArgumentException("One conditional covariance rate per output is required.", nameof(log));
        var scale = new double[Count, Count];
        for (int k = 0; k < Count; k++) scale[k, k] = 2 * Math.Exp(log[k]);
        return InverseWishart.Of(Nu + Count - 1, scale);
    }
}
