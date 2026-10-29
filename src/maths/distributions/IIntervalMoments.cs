using System;

namespace Maths.Distributions;

/// <summary>
/// A measure's interval moments ∫ (x − c)ᵏ dF(x) over <c>(lower, upper]</c>, k = 0, 1, …, about a centre c the caller
/// picks: the currency piecewise-polynomial expectations reduce to, one interval at a time. Either bound may be infinite
/// where the measure allows it. The interval is half-open, as <see cref="IContinuousDistribution.IntervalProbability"/>'s,
/// so consecutive intervals count a discrete measure's atoms once; a continuous measure gives the same moments on the
/// closed interval.
/// </summary>
public interface IIntervalMoments
{
    /// <summary>
    /// Fills <paramref name="moments"/>[k] with e^(−<paramref name="logScale"/>)·∫ (x − <paramref name="centre"/>)ᵏ dF(x)
    /// over (<paramref name="lower"/>, <paramref name="upper"/>] for k below its length, and returns an estimate of the
    /// largest absolute error among them, on the same scale: 0 for a closed form or an exact sum, a quadrature's own
    /// estimate otherwise. Pass 0 for the moments themselves; a truncation passes its log mass, so moments far in a tail
    /// keep their digits where the unscaled values would underflow. A centre inside the interval keeps the powers from
    /// cancelling.
    /// </summary>
    double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments);
}
