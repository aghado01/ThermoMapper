using System;
using System.Collections.Generic;

namespace Maths.Distributions;

/// <summary>
/// A finite mixture Σₖ wₖFₖ of continuous distributions, the weights normalized to sum to 1, built with <see cref="Of"/>.
/// Densities, tails and interval probabilities are weighted sums, formed in log space; interval moments are the
/// components' weighted moments, so every component must have them.
/// <para>Truncation composes two ways. <see cref="Truncated.Of{T}"/> of a mixture keeps the weights as given and
/// renormalizes the whole mixture on the interval, so a component's share follows its mass there; a mixture of
/// truncated components keeps the weights as the components' shares on the interval.</para>
/// </summary>
public readonly struct Mixture : IContinuousDistribution, IIntervalMoments
{
    private readonly IContinuousDistribution[] _components;
    private readonly double[] _weights, _logWeights;

    private Mixture(IContinuousDistribution[] components, double[] weights)
    {
        _components = components;
        _weights = weights;
        _logWeights = Array.ConvertAll(weights, Math.Log);
    }

    /// <summary>
    /// The mixture of <paramref name="components"/> with <paramref name="weights"/>: as many, nonnegative and finite, with
    /// a positive sum, which is normalized away.
    /// </summary>
    public static Mixture Of(IReadOnlyList<IContinuousDistribution> components, IReadOnlyList<double> weights)
    {
        ArgumentNullException.ThrowIfNull(components);
        ArgumentNullException.ThrowIfNull(weights);
        if (components.Count == 0) throw new ArgumentException("A mixture needs a component.", nameof(components));
        if (weights.Count != components.Count) throw new ArgumentException("A mixture needs one weight per component.", nameof(weights));
        var c = new IContinuousDistribution[components.Count];
        var w = new double[c.Length];
        double total = 0.0;
        for (int k = 0; k < c.Length; k++)
        {
            c[k] = components[k] ?? throw new ArgumentException($"Component {k} is null.", nameof(components));
            if (!(weights[k] >= 0.0) || double.IsPositiveInfinity(weights[k]))
                throw new ArgumentOutOfRangeException(nameof(weights), $"Weight {k} must be nonnegative and finite.");
            w[k] = weights[k];
            total += w[k];
        }
        if (!(total > 0.0)) throw new ArgumentException("A mixture's weights need a positive sum.", nameof(weights));
        for (int k = 0; k < w.Length; k++) w[k] /= total;
        return new Mixture(c, w);
    }

    /// <summary>The components.</summary>
    public IReadOnlyList<IContinuousDistribution> Components => Checked()._components;

    /// <summary>The normalized weights.</summary>
    public IReadOnlyList<double> Weights => Checked()._weights;

    public double LogPdf(double x) => LogSum(static (d, v, _) => d.LogPdf(v), x, 0.0);
    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));
    public double LogCdf(double x) => LogSum(static (d, v, _) => d.LogCdf(v), x, 0.0);
    public double LogSurvival(double x) => LogSum(static (d, v, _) => d.LogSurvival(v), x, 0.0);
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return double.NegativeInfinity;
        return LogSum(static (d, a, b) => d.LogIntervalProbability(a, b), lower, upper);
    }

    /// <inheritdoc/>
    public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
    {
        Checked();
        moments.Clear();
        Span<double> each = moments.Length <= 32 ? stackalloc double[moments.Length] : new double[moments.Length];
        double error = 0.0;
        for (int k = 0; k < _components.Length; k++)
        {
            if (_weights[k] == 0.0) continue;
            if (_components[k] is not IIntervalMoments component)
                throw new NotSupportedException($"Mixture component {k} ({_components[k].GetType().Name}) has no interval moments.");
            // wₖ rides the scale, so a component far in its tail keeps its digits.
            error += component.IntervalMoments(lower, upper, centre, logScale - _logWeights[k], each);
            for (int m = 0; m < moments.Length; m++) moments[m] += each[m];
        }
        return error;
    }

    // ln Σₖ wₖ·exp(f(Fₖ, a, b)).
    private double LogSum(Func<IContinuousDistribution, double, double, double> f, double a, double b)
    {
        Checked();
        Span<double> terms = _components.Length <= 32 ? stackalloc double[_components.Length] : new double[_components.Length];
        double max = double.NegativeInfinity;
        for (int k = 0; k < terms.Length; k++)
        {
            terms[k] = _weights[k] == 0.0 ? double.NegativeInfinity : _logWeights[k] + f(_components[k], a, b);
            max = Math.Max(max, terms[k]);
        }
        if (double.IsNegativeInfinity(max)) return double.NegativeInfinity;
        double sum = 0.0;
        foreach (double t in terms) sum += Math.Exp(t - max);
        return max + Math.Log(sum);
    }

    private Mixture Checked()
        => _components is null ? throw new InvalidOperationException("A default Mixture has no components; build one with Mixture.Of.") : this;
}
