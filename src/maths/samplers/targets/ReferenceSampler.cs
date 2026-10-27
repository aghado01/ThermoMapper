using Maths.Rng;

namespace Maths.Samplers;

/// <summary>
/// Exact draws from a target's b = 0 member on one path, normalized and restricted to the target's support: the prior
/// on the likelihood path, the normalized reference measure on the powered path. A tempering ladder whose hot end is
/// b = 0 draws its hot rung from it every scan, so each state that reaches the hot end restarts from an independent
/// exact draw, the setting of Syed et al.'s (2022) round-trip analysis. An SMC run can start from it.
/// </summary>
public interface IReferenceSampler<TState>
{
    /// <summary>The path whose b = 0 member this samples.</summary>
    TemperingPath Path { get; }

    /// <summary>One exact draw, with its evaluation under the target, using <paramref name="rng"/> only.</summary>
    (TState State, ITargetEvaluation Evaluation) Draw(Xoshiro256PlusPlus rng);

    /// <summary>Candidates evaluated so far, accepted ones included: a rejection sampler's cost, which a comparison
    /// at equal work counts in place of one step per draw.</summary>
    long Attempts { get; }
}
