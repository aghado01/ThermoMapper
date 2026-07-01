using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// Parallel tempering (replica exchange) over a ladder of reversible-jump chains at inverse temperatures
/// 1 = β₀ &gt; β₁ &gt; … &gt; β_{L-1} &gt; 0. Each replica samples the β-tempered target; after every sweep, adjacent
/// replicas attempt a Metropolis state swap — accept ∝ <c>exp((β_i − β_{i+1})(L_{i+1} − L_i))</c>, L the
/// untempered log-posterior — so the hot replicas, which cross modal barriers freely, ferry the cold replica
/// (β = 1, the true posterior) between modes. The cold chain is the sampler; the rest are scaffolding. For
/// multimodal targets (e.g. a multi-peak T_c posterior) where a single chain would lock into one mode.
/// As an <see cref="IChain{TState}"/> the ladder presents as a single chain whose draw is the cold state and whose
/// acceptance is the cold replica's, so the ensemble shell can drive a tempered replica like any other.
/// </summary>
public sealed class ParallelTempering<TState> : IChain<TState>
{
    private readonly ReversibleJumpChain<TState>[] _chains;
    private readonly double[] _betas;
    private readonly Xoshiro256PlusPlus _swapRng;

    /// <summary>Adjacent-pair swap attempts so far.</summary>
    public long SwapAttempts { get; private set; }

    /// <summary>Adjacent-pair swaps accepted so far; <c>SwapAccepts / SwapAttempts</c> is the exchange rate.</summary>
    public long SwapAccepts { get; private set; }

    /// <param name="moves">Shared move set for every replica.</param>
    /// <param name="baseTarget">The untempered target (β = 1); replica i samples <c>baseTarget^{β_i}</c>.</param>
    /// <param name="start">Common initial state for all replicas.</param>
    /// <param name="betas">Inverse-temperature ladder, descending from 1 (cold, index 0) to β_min &gt; 0.</param>
    /// <param name="masterSeed">Seeds the per-replica streams and the swap stream via the <see cref="SeedTree"/>.</param>
    public ParallelTempering(IReadOnlyList<IRjMove<TState>> moves, IRjTarget<TState> baseTarget, TState start,
                             double[] betas, int masterSeed)
    {
        ArgumentNullException.ThrowIfNull(moves);
        ArgumentNullException.ThrowIfNull(baseTarget);
        ArgumentNullException.ThrowIfNull(betas);
        if (betas.Length < 2) throw new ArgumentException("Tempering needs at least two levels.", nameof(betas));
        foreach (double b in betas)
            if (!(b > 0.0)) throw new ArgumentException("Inverse temperatures must be positive.", nameof(betas));

        _betas = betas;
        int[] seeds = SeedTree.Derive(masterSeed, betas.Length + 1);
        _chains = new ReversibleJumpChain<TState>[betas.Length];
        for (int i = 0; i < betas.Length; i++)
            _chains[i] = new ReversibleJumpChain<TState>(
                moves, new TemperedTarget<TState>(baseTarget, betas[i]), start, new Xoshiro256PlusPlus(seeds[i]));
        _swapRng = new Xoshiro256PlusPlus(seeds[betas.Length]);
    }

    /// <summary>The current cold-chain (β = 1) state — the posterior sample.</summary>
    public TState Cold => _chains[0].Current;

    /// <summary>The cold replica (β = 1) itself — the chain whose draws are the posterior, exposed so a consumer
    /// can read its acceptance / per-move diagnostics. Its state advances both by its own moves and by swaps.</summary>
    public ReversibleJumpChain<TState> ColdChain => _chains[0];

    /// <summary>Advance every replica one step, then sweep adjacent-replica swaps. Returns the cold state.</summary>
    public TState Step()
    {
        for (int i = 0; i < _chains.Length; i++) _chains[i].Step();

        for (int i = 0; i < _chains.Length - 1; i++)
        {
            SwapAttempts++;
            double li = _chains[i].CurrentLogTarget / _betas[i];           // untempered log π at replica i
            double li1 = _chains[i + 1].CurrentLogTarget / _betas[i + 1];
            double logA = (_betas[i] - _betas[i + 1]) * (li1 - li);
            if (logA >= 0.0 || _swapRng.NextDouble() < Math.Exp(logA))
            {
                SwapAccepts++;
                TState xi = _chains[i].Current, xi1 = _chains[i + 1].Current;
                _chains[i].SetState(xi1, _betas[i] * li1);
                _chains[i + 1].SetState(xi, _betas[i + 1] * li);
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
}
