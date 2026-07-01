using System;
using Maths.Rng;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// Decorates a move with a state-dependent selection weight, leaving its proposal (auxiliary draw, bijection,
/// Jacobian) and reverse pairing untouched. This is the seam for Green's state-dependent move-selection
/// probabilities <c>j_m(x)</c> — e.g. the DMGK birth/death schedule <c>b_k = c·min{1, p(k+1)/p(k)}</c> — layered
/// on a move without re-deriving its mechanics. Only meaningful because the engine owns the selection ratio.
/// </summary>
public sealed class ScheduledMove<TState> : IRjMove<TState>
{
    private readonly IRjMove<TState> _inner;
    private readonly Func<TState, double> _weight;

    public ScheduledMove(IRjMove<TState> inner, Func<TState, double> weight)
    {
        ArgumentNullException.ThrowIfNull(inner);
        ArgumentNullException.ThrowIfNull(weight);
        _inner = inner;
        _weight = weight;
    }

    public string Key => _inner.Key;
    public string ReverseKey => _inner.ReverseKey;
    public double Weight(TState state) => _weight(state);
    public Proposal<TState>? Propose(TState current, Xoshiro256PlusPlus rng) => _inner.Propose(current, rng);
}
