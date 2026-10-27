using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Samplers.Smc;

/// <summary>
/// Gaussian random-walk Metropolis for real-vector particles, scaled at each stage from the population: each
/// coordinate's step is a factor times its weighted spread across the particles, 2.38/√d by default (Roberts, Gelman
/// and Gilks 1997). The scale is fixed for the whole stage once <see cref="Prepare"/> has read the population, so every
/// move leaves the stage invariant. A move proposes a new vector and never writes to the old one, so particles that
/// share a vector after resampling stay independent.
/// </summary>
public sealed class RandomWalkKernel : ISmcKernel<double[]>
{
    private readonly double _factor;
    private double[] _scale = [];

    /// <param name="factor">The step as a multiple of each coordinate's spread; NaN takes 2.38/√d.</param>
    public RandomWalkKernel(double factor = double.NaN)
    {
        if (!double.IsNaN(factor) && !(factor > 0.0 && double.IsFinite(factor)))
            throw new ArgumentOutOfRangeException(nameof(factor), "The factor must be finite and positive.");
        _factor = factor;
    }

    /// <summary>The per-coordinate steps set by the last <see cref="Prepare"/>.</summary>
    public IReadOnlyList<double> Scale => _scale;

    /// <inheritdoc/>
    public void Prepare(SmcStage<double[]> stage, IReadOnlyList<double[]> particles, ReadOnlySpan<double> weights)
    {
        int d = particles.Count > 0 ? particles[0].Length : 0;
        double factor = double.IsNaN(_factor) ? 2.38 / Math.Sqrt(Math.Max(d, 1)) : _factor;
        var scale = new double[d];
        for (int j = 0; j < d; j++)
        {
            double mean = 0.0;
            for (int i = 0; i < particles.Count; i++) mean += weights[i] * particles[i][j];
            double variance = 0.0;
            for (int i = 0; i < particles.Count; i++)
            {
                double e = particles[i][j] - mean;
                variance += weights[i] * e * e;
            }
            scale[j] = factor * Math.Sqrt(Math.Max(variance, 1e-24 * (1.0 + mean * mean)));
        }
        _scale = scale;
    }

    /// <inheritdoc/>
    public SmcMove<double[]> Move(double[] state, ITargetEvaluation evaluation, SmcStage<double[]> stage, Xoshiro256PlusPlus rng)
    {
        if (_scale.Length != state.Length) throw new InvalidOperationException("Prepare has not read a population of this dimension.");
        var proposal = new double[state.Length];
        for (int j = 0; j < state.Length; j++) proposal[j] = state[j] + _scale[j] * Normal.SampleStandard(rng);
        ITargetEvaluation proposed = stage.Target.Evaluate(proposal);
        double logRatio = stage.LogDensity(proposed) - stage.LogDensity(evaluation);
        return Math.Log(rng.NextDouble()) < logRatio
            ? new SmcMove<double[]>(proposal, proposed, true)
            : new SmcMove<double[]>(state, evaluation, false);
    }
}
