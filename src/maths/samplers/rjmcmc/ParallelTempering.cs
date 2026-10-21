using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// Parallel tempering (replica exchange) over a ladder of reversible-jump chains at inverse temperatures
/// 1 = β₀ &gt; β₁ &gt; … &gt; β_{L-1}, all on one declared <see cref="TemperingPath"/>: the powered posterior
/// <c>(prior · likelihood)^β</c> or the likelihood-tempered path <c>prior · likelihood^β</c>. After every sweep,
/// adjacent replicas attempt a Metropolis state swap, each replica's density evaluated at the other's state from the
/// cached evaluations — so the hot replicas, which cross modal barriers freely, ferry the cold replica (β = 1, the true
/// posterior) between modes. The cold chain is the sampler; the rest are scaffolding. For multimodal targets (e.g. a
/// multi-peak T_c posterior) where a single chain would lock into one mode.
/// As an <see cref="IChain{TState}"/> the ladder presents as a single chain whose draw is the cold state and whose
/// acceptance is the cold replica's, so the ensemble shell can drive a tempered replica like any other.
/// </summary>
public sealed class ParallelTempering<TState> : IChain<TState>
{
    private readonly ReversibleJumpChain<TState>[] _chains;
    private readonly Xoshiro256PlusPlus _swapRng;

    /// <summary>Adjacent-pair swap attempts so far.</summary>
    public long SwapAttempts { get; private set; }

    /// <summary>Adjacent-pair swaps accepted so far; <c>SwapAccepts / SwapAttempts</c> is the exchange rate.</summary>
    public long SwapAccepts { get; private set; }

    /// <param name="movesForReplica">Builds replica i's own move set (i = 0 is the cold replica). It is called once
    /// per replica, so any mutable move or kernel state is private to that replica.</param>
    /// <param name="target">The target; replica i samples it at <c>Temperature(β_i, path)</c>.</param>
    /// <param name="start">Common initial state for all replicas.</param>
    /// <param name="betas">The ladder: finite, strictly descending, cold value exactly 1 at index 0. The hot end is
    /// above 0 on the powered path, where β = 0 would be the bare reference measure; on the likelihood path β = 0 is
    /// the prior, which a target with an improper nuisance prior refuses. The array is not retained.</param>
    /// <param name="path">The tempering path every replica takes.</param>
    /// <param name="masterSeed">Seeds the per-replica streams and the swap stream via the <see cref="SeedTree"/>.</param>
    public ParallelTempering(Func<int, IReadOnlyList<IRjMove<TState>>> movesForReplica, ISamplingTarget<TState> target,
                             TState start, double[] betas, TemperingPath path, int masterSeed)
    {
        ArgumentNullException.ThrowIfNull(movesForReplica);
        ArgumentNullException.ThrowIfNull(target);
        ValidateLadder(betas, path);

        Path = path;
        int[] seeds = SeedTree.Derive(masterSeed, betas.Length + 1);
        _chains = new ReversibleJumpChain<TState>[betas.Length];
        for (int i = 0; i < betas.Length; i++)
            _chains[i] = new ReversibleJumpChain<TState>(
                movesForReplica(i), target, start, new Xoshiro256PlusPlus(seeds[i]), new Temperature(betas[i], path));
        _swapRng = new Xoshiro256PlusPlus(seeds[betas.Length]);
    }

    /// <summary>The path every replica takes.</summary>
    public TemperingPath Path { get; }

    /// <summary>The current cold-chain (β = 1) state — the posterior sample.</summary>
    public TState Cold => _chains[0].Current;

    /// <summary>The cold replica (β = 1) itself — the chain whose draws are the posterior, exposed so a consumer
    /// can read its acceptance / per-move diagnostics. Its state advances both by its own moves and by swaps.</summary>
    public ReversibleJumpChain<TState> ColdChain => _chains[0];

    /// <summary>Every replica, cold first, for per-rung diagnostics. Each one's <see cref="ReversibleJumpChain{TState}.Temperature"/>
    /// is its rung. Stepping a replica directly takes it out of the ladder's schedule.</summary>
    public IReadOnlyList<ReversibleJumpChain<TState>> Replicas => _chains;

    /// <summary>
    /// Advance every replica one step, then sweep adjacent-replica swaps, cold to hot. A swap of x_i and x_j is
    /// accepted with <c>log α = d_i(x_j) + d_j(x_i) − d_i(x_i) − d_j(x_j)</c>, d_i replica i's log density, which on
    /// the powered path reduces to <c>(β_i − β_j)(L_j − L_i)</c>. Later pairs see earlier swaps. Returns the cold
    /// state.
    /// </summary>
    public TState Step()
    {
        for (int i = 0; i < _chains.Length; i++) _chains[i].Step();

        for (int i = 0; i < _chains.Length - 1; i++)
        {
            SwapAttempts++;
            ReversibleJumpChain<TState> a = _chains[i], b = _chains[i + 1];
            ITargetEvaluation ea = a.CurrentEvaluation, eb = b.CurrentEvaluation;
            double logA = a.Temperature.LogDensity(eb) + b.Temperature.LogDensity(ea)
                        - a.CurrentLogDensity - b.CurrentLogDensity;
            // Drawn on every attempt, so the swap stream does not advance by the decision: between states whose targets
            // agree to rounding, rounding decides a tie's sign, and it must not decide whether a uniform is consumed.
            double u = _swapRng.NextDouble();
            if (logA >= 0.0 || u < Math.Exp(logA))
            {
                SwapAccepts++;
                TState xa = a.Current, xb = b.Current;
                a.SetState(xb, eb);
                b.SetState(xa, ea);
            }
        }
        return _chains[0].Current;
    }

    // IChain view: the ladder advances by one swept Step; its acceptance is the cold replica's own move acceptance
    // (swaps move state via SetState and are tracked separately as SwapAccepts/SwapAttempts).
    long IChain<TState>.Accepted => _chains[0].Accepted;
    long IChain<TState>.Attempts => _chains[0].Attempts;

    /// <summary>Geometric inverse-temperature ladder from 1 (cold) down to <paramref name="betaMin"/>.</summary>
    public static double[] GeometricLadder(int levels, double betaMin)
    {
        if (levels < 2) throw new ArgumentOutOfRangeException(nameof(levels), "Need at least two levels.");
        if (!(betaMin > 0.0 && betaMin < 1.0)) throw new ArgumentOutOfRangeException(nameof(betaMin), "β_min ∈ (0,1).");
        var b = new double[levels];
        for (int i = 0; i < levels; i++)
            b[i] = Math.Pow(betaMin, i / (double)(levels - 1));   // 1 … β_min
        return b;
    }

    private static void ValidateLadder(double[] betas, TemperingPath path)
    {
        ArgumentNullException.ThrowIfNull(betas);
        if (betas.Length < 2) throw new ArgumentException("Tempering needs at least two levels.", nameof(betas));
        if (betas[0] != 1.0) throw new ArgumentException("The cold rung must be exactly β = 1.", nameof(betas));
        for (int i = 1; i < betas.Length; i++)
            if (!(betas[i] < betas[i - 1]))
                throw new ArgumentException("The ladder must be strictly descending.", nameof(betas));
        double hot = betas[^1];
        if (!(path == TemperingPath.Powered ? hot > 0.0 : hot >= 0.0))
            throw new ArgumentException(path == TemperingPath.Powered
                ? "On the powered path every β must be positive."
                : "Every β must be non-negative.", nameof(betas));
    }
}
