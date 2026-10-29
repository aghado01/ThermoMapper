using System;
using System.Collections.Generic;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass.Presentation;

/// <summary>
/// One scalar's posterior summary over weighted draws (stage 9): the weighted mean and sd, the least and largest draw that
/// carries weight, and the quantiles at <see cref="Probabilities"/> through the shared weighted quantile,
/// <see cref="DrawSample.Quantile"/> (type 7 under equal weights). It reduces the draws to values only; formatting is an
/// adapter's. Draws of weight 0 drop out, and with none carrying weight every value is NaN.
/// </summary>
/// <param name="Mean">Σ wₛxₛ / Σ wₛ, summed in draw order.</param>
/// <param name="Sd">The weighted standard deviation about the mean, √(Σ wₛ(xₛ − x̄)² / Σ wₛ).</param>
/// <param name="Min">The least draw that carries weight.</param>
/// <param name="Max">The largest draw that carries weight.</param>
/// <param name="Probabilities">The quantiles' probabilities, in the order given.</param>
/// <param name="Quantiles">One quantile per probability.</param>
/// <param name="Draws">The draws that carry weight.</param>
public sealed record PosteriorSummary(double Mean, double Sd, double Min, double Max, IReadOnlyList<double> Probabilities,
                                      double[] Quantiles, int Draws)
{
    /// <summary>The probabilities a summary reports unless others are asked for: the median, inside the central 50% and 95%
    /// intervals.</summary>
    public static IReadOnlyList<double> DefaultProbabilities { get; } = Array.AsReadOnly([0.025, 0.25, 0.5, 0.75, 0.975]);

    /// <summary>
    /// The summary of <paramref name="values"/> under <paramref name="weights"/> (null for equal weights) at
    /// <paramref name="probabilities"/> (<see cref="DefaultProbabilities"/> when null). Refused: a weight count that is not
    /// the value count, a weight negative or not finite, a value not finite, and a probability outside [0, 1].
    /// </summary>
    public static PosteriorSummary Of(IReadOnlyList<double> values, IReadOnlyList<double>? weights = null, IReadOnlyList<double>? probabilities = null)
        => Over(values, weights, Checked(probabilities));

    // The summary at probabilities already checked, shared by every summary of one call.
    internal static PosteriorSummary Over(IReadOnlyList<double> values, IReadOnlyList<double>? weights, IReadOnlyList<double> probabilities)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (weights is not null && weights.Count != values.Count)
            throw new ArgumentException($"The {weights.Count} weights do not match the {values.Count} values.", nameof(weights));
        var sample = new DrawSample();
        double sum = 0.0, mass = 0.0, min = double.PositiveInfinity, max = double.NegativeInfinity;
        int draws = 0;
        for (int d = 0; d < values.Count; d++)
        {
            double x = values[d], w = weights?[d] ?? 1.0;
            if (!double.IsFinite(x)) throw new ArgumentException("A summary's values must be finite.", nameof(values));
            sample.Add(x, w);   // refuses a negative or infinite weight
            mass += w;
            sum += w * x;
            if (w == 0.0) continue;
            draws++;
            min = Math.Min(min, x);
            max = Math.Max(max, x);
        }
        var quantiles = new double[probabilities.Count];
        if (draws == 0)
        {
            Array.Fill(quantiles, double.NaN);
            return new PosteriorSummary(double.NaN, double.NaN, double.NaN, double.NaN, probabilities, quantiles, 0);
        }
        double mean = sum / mass, squares = 0.0;
        for (int d = 0; d < values.Count; d++)
        {
            double w = weights?[d] ?? 1.0, e = values[d] - mean;
            squares += w * e * e;
        }
        for (int q = 0; q < quantiles.Length; q++) quantiles[q] = sample.Quantile(probabilities[q]);
        return new PosteriorSummary(mean, Math.Sqrt(squares / mass), min, max, probabilities, quantiles, draws);
    }

    // The probabilities, checked and copied: the defaults when null.
    internal static IReadOnlyList<double> Checked(IReadOnlyList<double>? probabilities, bool open = false)
    {
        if (probabilities is null) return DefaultProbabilities;
        var copy = new double[probabilities.Count];
        for (int q = 0; q < copy.Length; q++)
        {
            double p = probabilities[q];
            if (open ? !(p > 0.0 && p < 1.0) : !(p >= 0.0 && p <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(probabilities), $"Probability {p} lies outside {(open ? "(0, 1)" : "[0, 1]")}.");
            copy[q] = p;
        }
        return Array.AsReadOnly(copy);
    }
}
