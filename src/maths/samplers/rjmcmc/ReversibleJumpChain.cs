using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// One reversible-jump Markov chain (Green 1995): a pure, sequential sampler that walks a carrier
/// <typeparamref name="TState"/> by weighted random choice among <see cref="IRjMove{TState}"/> moves, accepting
/// each proposal under the full Green acceptance ratio. The engine assembles all four factors — target ratio,
/// move-selection ratio <c>j_{m'}(x')/j_m(x)</c> (it owns this: it sees the whole palette and each state's
/// normalization), and the move's auxiliary-density ratio and Jacobian. Selection weights may be state-dependent;
/// a move declares its reverse via <see cref="IRjMove{TState}.ReverseKey"/>. The chain is the unit of parallelism
/// — it shares nothing — so an ensemble layer fans out independent chains and pools their reductions.
/// Model-agnostic: the same engine serves free-knot BARS, partition samplers, and any trans-dimensional target.
/// The chain samples its <see cref="ISamplingTarget{TState}"/> at one <see cref="Samplers.Temperature"/>, the
/// posterior unless told otherwise, and caches the current state's evaluation, so each step evaluates only the
/// candidate. As an <see cref="IChain{TState}"/> it plugs straight into the kernel-agnostic ensemble shell.
/// </summary>
public sealed class ReversibleJumpChain<TState> : IChain<TState>
{
    private readonly IReadOnlyList<IRjMove<TState>> _moves;
    private readonly Dictionary<string, IRjMove<TState>> _byKey;
    private ISamplingTarget<TState> _target;
    private readonly Xoshiro256PlusPlus _rng;
    private ITargetEvaluation _currentEvaluation;   // Current's; only the candidate is evaluated per step
    private double _currentLogDensity;              // its log density at Temperature
    private readonly long[] _moveAttempts;
    private readonly long[] _moveAccepts;

    /// <summary>The most recently accepted (or retained) carrier value.</summary>
    public TState Current { get; private set; }

    /// <summary>Proposals attempted (a step that selected an <i>available</i> move) so far.</summary>
    public long Attempts { get; private set; }

    /// <summary>Proposals accepted so far; <c>Accepted / Attempts</c> is the running acceptance rate.</summary>
    public long Accepted { get; private set; }

    /// <param name="moves">The move set; combined by weighted random choice. Keys must be unique, and each
    /// <see cref="IRjMove{TState}.ReverseKey"/> must name a move in the set that names it back (or itself).</param>
    /// <param name="target">The target; the chain samples its log density at <paramref name="temperature"/>.</param>
    /// <param name="start">The chain's initial state.</param>
    /// <param name="rng">This chain's private RNG stream (the ensemble layer owns the seed tree).</param>
    /// <param name="temperature">The tempered density the chain samples; null is the posterior, b = 1.</param>
    public ReversibleJumpChain(
        IReadOnlyList<IRjMove<TState>> moves,
        ISamplingTarget<TState> target,
        TState start,
        Xoshiro256PlusPlus rng,
        Temperature? temperature = null)
    {
        ArgumentNullException.ThrowIfNull(moves);
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(rng);
        if (moves.Count == 0)
            throw new ArgumentException("A reversible-jump chain needs at least one move.", nameof(moves));

        _byKey = new Dictionary<string, IRjMove<TState>>(moves.Count);
        foreach (IRjMove<TState> m in moves)
        {
            ArgumentNullException.ThrowIfNull(m);
            if (!_byKey.TryAdd(m.Key, m))
                throw new ArgumentException($"Duplicate move key '{m.Key}'.", nameof(moves));
        }
        foreach (IRjMove<TState> m in moves)
        {
            if (!_byKey.TryGetValue(m.ReverseKey, out IRjMove<TState>? reverse))
                throw new ArgumentException($"Move '{m.Key}' names a reverse '{m.ReverseKey}' not in the palette.", nameof(moves));
            if (reverse.ReverseKey != m.Key)
                throw new ArgumentException(
                    $"Move '{m.Key}' names '{m.ReverseKey}' as its reverse, but '{m.ReverseKey}' names '{reverse.ReverseKey}'. " +
                    "Each pair must name each other, or a move must name itself.", nameof(moves));
        }

        _moves = moves;
        _target = target;
        _rng = rng;
        Temperature = temperature ?? Temperature.Posterior;
        _moveAttempts = new long[moves.Count];
        _moveAccepts = new long[moves.Count];
        Current = start;
        _currentEvaluation = target.Evaluate(start);
        _currentLogDensity = Temperature.LogDensity(_currentEvaluation);
    }

    /// <summary>The tempered density this chain samples.</summary>
    public Temperature Temperature { get; }

    /// <summary>
    /// Advance the chain one Draw: select a move ∝ weight(Current), propose, and accept with probability
    /// <c>min(1, exp(Δlog-target + log-auxiliary-ratio + log-Jacobian + log-selection-ratio))</c>, where the
    /// selection ratio is <c>[w_{m'}(x')/ΣW(x')] / [w_m(x)/ΣW(x)]</c> and the target term is the log density at
    /// <see cref="Temperature"/>. A move unavailable from the current state (null proposal) leaves the chain in place
    /// and does not count as an attempt; a proposal whose reverse has zero weight is rejected (no detailed balance),
    /// and so is a candidate outside the target's support (log density −∞ or NaN), whatever its proposal ratio.
    /// Returns the (possibly unchanged) <see cref="Current"/>.
    /// </summary>
    public TState Step()
    {
        double totalCur = TotalWeight(Current);
        if (!(totalCur > 0.0))
            return Current;   // no move available here → stay

        int idx = SelectMove(Current, totalCur);
        IRjMove<TState> move = _moves[idx];
        if (move.Propose(Current, _rng) is not Proposal<TState> p)
            return Current;   // move unavailable from here → stay

        Attempts++;
        _moveAttempts[idx]++;

        double totalCand = TotalWeight(p.Candidate);
        double wReverse = _byKey[move.ReverseKey].Weight(p.Candidate);
        if (!(wReverse > 0.0) || !(totalCand > 0.0))
            return Current;   // proposed state can't reverse this move → reject

        // log[ j_{m'}(x') / j_m(x) ] = log[ w_{m'}(x')/ΣW(x') ] − log[ w_m(x)/ΣW(x) ]
        double logSelection = Math.Log(wReverse / totalCand) - Math.Log(move.Weight(Current) / totalCur);

        // Off the support, the proposal ratio must not decide the move: a sampler that underflows onto the
        // support's edge can report an infinite ratio there, and −∞ + ∞ is NaN. The uniform is still drawn.
        ITargetEvaluation candidate = _target.Evaluate(p.Candidate);
        double candidateLogDensity = Temperature.LogDensity(candidate);
        double logAccept = candidateLogDensity > double.NegativeInfinity
            ? (candidateLogDensity - _currentLogDensity) + p.LogProposalRatio + p.LogJacobian + logSelection
            : double.NegativeInfinity;

        if (logAccept >= 0.0 || _rng.NextDouble() < Math.Exp(logAccept))
        {
            Current = p.Candidate;
            _currentEvaluation = candidate;
            _currentLogDensity = candidateLogDensity;
            Accepted++;
            _moveAccepts[idx]++;
        }
        return Current;
    }

    /// <summary>
    /// Per-move attempt/accept counts (cumulative, in move order) — the diagnostic behind acceptance monitoring
    /// and adaptive proposal-scale tuning. <see cref="IRjMove{TState}.Key"/> labels each entry.
    /// </summary>
    public IReadOnlyList<MoveStat> MoveStats()
    {
        var stats = new MoveStat[_moves.Count];
        for (int i = 0; i < _moves.Count; i++)
            stats[i] = new MoveStat(_moves[i].Key, _moveAttempts[i], _moveAccepts[i]);
        return stats;
    }

    /// <summary>The cached evaluation of <see cref="Current"/>. Replica exchange reads it, and queries it at other
    /// temperatures, without re-evaluating the state.</summary>
    public ITargetEvaluation CurrentEvaluation => _currentEvaluation;

    /// <summary>The cached log density of <see cref="Current"/> at <see cref="Temperature"/>.</summary>
    public double CurrentLogDensity => _currentLogDensity;

    /// <summary>
    /// Force the chain to <paramref name="state"/> with its already-known <paramref name="evaluation"/>, for
    /// replica-exchange swaps, where another replica hands over its state and evaluation. The density is recomputed
    /// at this chain's own temperature; the state is not re-evaluated. The caller is responsible for
    /// <paramref name="evaluation"/> being this chain's target's at <paramref name="state"/>.
    /// </summary>
    public void SetState(TState state, ITargetEvaluation evaluation)
    {
        ArgumentNullException.ThrowIfNull(evaluation);
        Current = state;
        _currentEvaluation = evaluation;
        _currentLogDensity = Temperature.LogDensity(evaluation);
    }

    /// <summary>
    /// Replace the target the chain samples, keeping its state, temperature and cumulative counters, and re-evaluate
    /// <see cref="Current"/>. For a conditional update that changes the target between steps: new robustness weights
    /// or a hyperparameter draw.
    /// </summary>
    public void ReplaceTarget(ISamplingTarget<TState> target)
    {
        ArgumentNullException.ThrowIfNull(target);
        _target = target;
        _currentEvaluation = target.Evaluate(Current);
        _currentLogDensity = Temperature.LogDensity(_currentEvaluation);
    }

    private double TotalWeight(TState state)
    {
        double total = 0.0;
        for (int i = 0; i < _moves.Count; i++)
        {
            double w = _moves[i].Weight(state);
            if (w > 0.0) total += w;
        }
        return total;
    }

    private int SelectMove(TState state, double total)
    {
        if (_moves.Count == 1)
            return 0;

        double u = _rng.NextDouble() * total;
        double cumulative = 0.0;
        for (int i = 0; i < _moves.Count; i++)
        {
            double w = _moves[i].Weight(state);
            if (w <= 0.0) continue;
            cumulative += w;
            if (u < cumulative)
                return i;
        }
        return _moves.Count - 1;   // floating-point guard
    }
}
