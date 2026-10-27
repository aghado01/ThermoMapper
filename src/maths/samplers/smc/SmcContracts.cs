using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Samplers.Smc;

/// <summary>
/// One distribution in an SMC sequence: a target at a temperature, with density
/// <see cref="Temperature.LogDensity"/> of the target's evaluation. A tempered run walks one target from b = 0 to 1;
/// a cut posterior walks through different conditional targets.
/// </summary>
public readonly record struct SmcStage<TState>(ISamplingTarget<TState> Target, Temperature Temperature)
{
    /// <summary>The stage's log density of an evaluated state, up to the target's additive constant.</summary>
    public double LogDensity(ITargetEvaluation evaluation) => Temperature.LogDensity(evaluation);
}

/// <summary>
/// The consumer's particles: how one is drawn from the first stage's distribution, and how one is copied when
/// resampling duplicates it. The model, its target sequence and its readout stay with the consumer.
/// </summary>
public interface ISmcParticles<TState>
{
    /// <summary>Draws a particle exactly from the first stage's distribution, which must be normalized (a proper prior
    /// or reference), with <paramref name="rng"/> its own stream.</summary>
    TState Draw(Xoshiro256PlusPlus rng);

    /// <summary>An independent copy, which can be mutated without touching the original: resampling's offspring.</summary>
    TState Fork(TState state);
}

/// <summary>
/// A Markov kernel that leaves a stage's density invariant: the mutation that diversifies the population after
/// reweighting and resampling. SMC does not repair an invalid kernel; a reversible-jump consumer supplies its moves
/// at the stage's temperature.
/// </summary>
public interface ISmcKernel<TState>
{
    /// <summary>
    /// Called once per stage before any particle moves, with the population and its normalized weights, so the kernel
    /// can tune itself from them; every particle then moves under the same fixed kernel, which keeps the stage
    /// invariant. With <see cref="SmcOptions.MaxDegreeOfParallelism"/> above 1, <see cref="Move"/> runs concurrently
    /// on distinct particles and must touch only its own state, evaluation and stream.
    /// </summary>
    void Prepare(SmcStage<TState> stage, IReadOnlyList<TState> particles, ReadOnlySpan<double> weights) { }

    /// <summary>One move from <paramref name="state"/>, whose evaluation under the stage's target is
    /// <paramref name="evaluation"/>.</summary>
    SmcMove<TState> Move(TState state, ITargetEvaluation evaluation, SmcStage<TState> stage, Xoshiro256PlusPlus rng);
}

/// <summary>A kernel's result: the new state, its evaluation under the stage's target, and whether it moved.</summary>
public readonly record struct SmcMove<TState>(TState State, ITargetEvaluation Evaluation, bool Accepted);
