using System;
using Maths.Samplers;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The geometric bridge between two emulator draws' conditionals (D09, stage 8),
/// (1 − t)·log p̃(θ | y, e_a) + t·log p̃(θ | y, e_b), written as a target on D05's likelihood path: its base is e_a's
/// conditional, log π(θ) + ℓ_a, and the difference Δ = ℓ_b − ℓ_a plays the likelihood, so the stage at t has density
/// log π + ℓ_a + t·Δ. At t = 0 it is e_a's conditional and at t = 1 e_b's, and an SMC run tempers t from 0 to 1 on it
/// (<see cref="SmcSampler.TemperToPosterior"/>), each step from t to t′ reweighting by (t′ − t)·Δ. Emulator draws are
/// discrete and their term arrays cannot be interpolated, so the bridge runs through the draws' induced densities, where
/// Mathews et al.'s (2024) linear tempering interpolates a continuous cut parameter. Both ends are proper posteriors, so
/// t = 0 needs no proper prior. The unnormalized densities' constants cancel in the normalized weights.
/// </summary>
public sealed class CalibrationBridge : ISamplingTarget<CalibrationState>
{
    /// <summary>The bridge from <paramref name="from"/>'s conditional (t = 0) to <paramref name="to"/>'s (t = 1), two
    /// draws of one model.</summary>
    public CalibrationBridge(CalibrationTarget from, CalibrationTarget to)
    {
        ArgumentNullException.ThrowIfNull(from);
        ArgumentNullException.ThrowIfNull(to);
        if (!ReferenceEquals(from.Model, to.Model)) throw new ArgumentException("A bridge joins two draws of one model.", nameof(to));
        From = from;
        To = to;
    }

    /// <summary>The draw the bridge leaves, e_a: its conditional is the bridge at t = 0.</summary>
    public CalibrationTarget From { get; }

    /// <summary>The draw the bridge reaches, e_b: its conditional is the bridge at t = 1.</summary>
    public CalibrationTarget To { get; }

    /// <summary>The model both draws belong to.</summary>
    public CalibrationModel Model => From.Model;

    /// <inheritdoc/>
    public ITargetEvaluation Evaluate(CalibrationState state) => EvaluateState(state);

    /// <summary>The evaluation of <paramref name="state"/> under both draws, each reading the state's cache.</summary>
    public BridgeEvaluation EvaluateState(CalibrationState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        return new BridgeEvaluation(From.EvaluateState(state), To.EvaluateState(state));
    }
}

/// <summary>
/// One state's evaluation under a <see cref="CalibrationBridge"/>: its two draws' evaluations, the base
/// log π + ℓ_a as the log prior and t·Δ as the log likelihood at t, Δ = ℓ_b − ℓ_a, separable in t.
/// </summary>
public sealed class BridgeEvaluation : ITargetEvaluation
{
    internal BridgeEvaluation(CalibrationEvaluation from, CalibrationEvaluation to)
    {
        From = from;
        To = to;
        bool off = double.IsNegativeInfinity(from.LogPrior);
        LogPrior = off ? double.NegativeInfinity : from.LogPrior + from.LogLikelihoodValue;
        Difference = off ? double.NegativeInfinity : to.LogLikelihoodValue - from.LogLikelihoodValue;
    }

    /// <summary>The evaluation under the draw the bridge leaves.</summary>
    public CalibrationEvaluation From { get; }

    /// <summary>The evaluation under the draw the bridge reaches.</summary>
    public CalibrationEvaluation To { get; }

    /// <summary>The bridge's base, log π + ℓ_a: e_a's untempered log posterior; −∞ off the prior's support.</summary>
    public double LogPrior { get; }

    /// <summary>Δ = ℓ_b − ℓ_a, the bridge's likelihood; −∞ off the support.</summary>
    public double Difference { get; }

    /// <inheritdoc/>
    public double LogLikelihood(double inverseTemperature)
    {
        if (!(inverseTemperature >= 0.0 && inverseTemperature <= 1.0)) throw new ArgumentOutOfRangeException(nameof(inverseTemperature));
        return double.IsNegativeInfinity(Difference) ? double.NegativeInfinity : inverseTemperature * Difference;
    }

    /// <inheritdoc/>
    public double LogLikelihoodSlope(double inverseTemperature) => Difference;
}
