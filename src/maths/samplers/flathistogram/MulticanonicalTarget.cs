using System;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// The multicanonical target π_w(x) ∝ prior(x)·w(ℓ₁(x)) over any <see cref="ISamplingTarget{TState}"/>: its evaluation
/// keeps the inner log prior and gives ln w(ℓ₁) as the log likelihood at b = 1, ℓ₁ being the inner evaluation's
/// untempered log likelihood. So π_w is an ordinary target sampled at <see cref="Temperature.Posterior"/>, and any kernel
/// that runs on a target, a reversible-jump chain's moves or a block scan, runs on it unchanged. It is the stage density
/// of the flat-histogram engine, met by wrapping the target rather than by generalizing <see cref="Temperature"/>.
/// <para>During learning the weights move between steps (<see cref="SamcLearning"/>): the run installs each new
/// <see cref="MulticanonicalWeights"/> here, and each chain reweighs its cached evaluation
/// (<see cref="MulticanonicalEvaluation.Reweighed"/>), so no state is evaluated again. In production they stay
/// frozen.</para>
/// </summary>
public sealed class MulticanonicalTarget<TState> : ISamplingTarget<TState>
{
    /// <param name="inner">The target whose prior and untempered log likelihood π_w is built from.</param>
    /// <param name="weights">The weight function.</param>
    public MulticanonicalTarget(ISamplingTarget<TState> inner, MulticanonicalWeights weights)
    {
        ArgumentNullException.ThrowIfNull(inner);
        ArgumentNullException.ThrowIfNull(weights);
        Inner = inner;
        Weights = weights;
    }

    /// <summary>The wrapped target.</summary>
    public ISamplingTarget<TState> Inner { get; }

    /// <summary>The weights candidates are evaluated under.</summary>
    public MulticanonicalWeights Weights { get; private set; }

    /// <summary>Evaluate a state: the inner evaluation once, then ln w of its ℓ₁.</summary>
    public MulticanonicalEvaluation Evaluate(TState state) => new(Inner.Evaluate(state), Weights);

    ITargetEvaluation ISamplingTarget<TState>.Evaluate(TState state) => Evaluate(state);

    // The learning's: install the weights the next steps run under. Chains then reweigh their current evaluations.
    internal void Install(MulticanonicalWeights weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        Weights = weights;
    }
}

/// <summary>
/// One state's evaluation under a <see cref="MulticanonicalTarget{TState}"/>: the inner evaluation, its untempered log
/// likelihood ℓ₁ and ln w(ℓ₁) under the weights it was taken at. The log likelihood is defined at b = 1 only, where it
/// is ln w(ℓ₁); other b are refused, since π_w is not on a tempering path. Off the inner support ℓ₁ is −∞ with no
/// likelihood work, and so is ln w.
/// </summary>
public sealed class MulticanonicalEvaluation : ITargetEvaluation
{
    /// <param name="inner">The inner target's evaluation of the state.</param>
    /// <param name="weights">The weights to read ln w from.</param>
    public MulticanonicalEvaluation(ITargetEvaluation inner, MulticanonicalWeights weights)
    {
        ArgumentNullException.ThrowIfNull(inner);
        ArgumentNullException.ThrowIfNull(weights);
        Inner = inner;
        LogPrior = inner.LogPrior;
        UntemperedLogLikelihood = double.IsNegativeInfinity(LogPrior) ? double.NegativeInfinity : inner.LogLikelihood(1.0);
        Weights = weights;
        LogWeight = double.IsNegativeInfinity(LogPrior) ? double.NegativeInfinity : weights.LogWeight(UntemperedLogLikelihood);
    }

    private MulticanonicalEvaluation(MulticanonicalEvaluation source, MulticanonicalWeights weights)
    {
        Inner = source.Inner;
        LogPrior = source.LogPrior;
        UntemperedLogLikelihood = source.UntemperedLogLikelihood;
        Weights = weights;
        LogWeight = double.IsNegativeInfinity(LogPrior) ? double.NegativeInfinity : weights.LogWeight(UntemperedLogLikelihood);
    }

    /// <summary>The inner target's evaluation.</summary>
    public ITargetEvaluation Inner { get; }

    /// <summary>The inner log prior, the base measure.</summary>
    public double LogPrior { get; }

    /// <summary>ℓ₁, the inner evaluation's log likelihood at b = 1: the energy the weights are a function of.</summary>
    public double UntemperedLogLikelihood { get; }

    /// <summary>The weights <see cref="LogWeight"/> was read from.</summary>
    public MulticanonicalWeights Weights { get; }

    /// <summary>ln w(ℓ₁).</summary>
    public double LogWeight { get; }

    /// <summary>ln w(ℓ₁) at b = 1; any other b is refused (<see cref="NotSupportedException"/>).</summary>
    public double LogLikelihood(double inverseTemperature)
    {
        if (inverseTemperature != 1.0)
            throw new NotSupportedException(
                $"A multicanonical evaluation is defined at b = 1 only, not at b = {inverseTemperature}: π_w is not on a tempering path.");
        return LogWeight;
    }

    /// <summary>The same state under other weights, from the cached ℓ₁: no inner evaluation is repeated.</summary>
    public MulticanonicalEvaluation Reweighed(MulticanonicalWeights weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        return ReferenceEquals(weights, Weights) ? this : new MulticanonicalEvaluation(this, weights);
    }
}
