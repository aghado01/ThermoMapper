using Maths.Rng;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// One move type in a reversible-jump hybrid sampler (Green 1995): a Transform on the carrier
/// <typeparamref name="TState"/> that draws auxiliary variables, maps to a candidate by a bijection, and reports
/// the within-move terms of the acceptance ratio. Green's four factors are split by responsibility: the target
/// ratio and the <i>move-selection</i> ratio <c>j_{m'}(x')/j_m(x)</c> are the engine's (it sees the whole
/// palette); the auxiliary-density ratio and the dimension-matching Jacobian are the move's (only it knows its
/// own bijection). A move declares its reverse via <see cref="ReverseKey"/> so the engine can find the reverse
/// selection probability.
/// </summary>
public interface IRjMove<TState>
{
    /// <summary>This move type's stable identity, unique within a chain's palette.</summary>
    string Key { get; }

    /// <summary>The key of the move type that reverses this one (equal to <see cref="Key"/> for a self-reverse move).</summary>
    string ReverseKey { get; }

    /// <summary>
    /// Unnormalized selection weight at <paramref name="state"/> — Green's <c>j_m(·)</c> up to the engine's
    /// per-state normalization. Constant for a fixed-rate palette, state-dependent in general. (Availability may
    /// also be signalled by a null <see cref="Propose"/> under the engine's stay-on-null convention.)
    /// </summary>
    double Weight(TState state);

    /// <summary>
    /// Propose a candidate from <paramref name="current"/>, or <see langword="null"/> if this move is unavailable
    /// here (e.g. a death move at the minimal dimension). The returned
    /// <see cref="Proposal{TState}.LogProposalRatio"/> is the within-move auxiliary-density ratio
    /// <c>log[q_{m'}(u') / q_m(u)]</c> only — the engine adds the move-selection and target ratios.
    /// </summary>
    Proposal<TState>? Propose(TState current, Xoshiro256PlusPlus rng);
}
