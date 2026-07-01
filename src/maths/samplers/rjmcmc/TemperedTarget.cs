using System;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// A target raised to inverse-temperature β: <c>log π_β(x) = β · log π(x)</c> — the flattened distribution a hot
/// replica samples in parallel tempering. β = 1 is the true target; β → 0 flattens toward uniform, lowering the
/// barriers between modes so a hot chain traverses them and ferries the cold chain across via swaps. Only the
/// target is tempered; the move's proposal/Hastings terms are untempered (the engine adds them at β = 1).
/// </summary>
public sealed class TemperedTarget<TState> : IRjTarget<TState>
{
    private readonly IRjTarget<TState> _base;

    /// <summary>Inverse temperature β &gt; 0.</summary>
    public double Beta { get; }

    public TemperedTarget(IRjTarget<TState> baseTarget, double beta)
    {
        ArgumentNullException.ThrowIfNull(baseTarget);
        if (!(beta > 0.0)) throw new ArgumentOutOfRangeException(nameof(beta), "Inverse temperature must be positive.");
        _base = baseTarget;
        Beta = beta;
    }

    public double LogPosterior(TState state) => Beta * _base.LogPosterior(state);
}
