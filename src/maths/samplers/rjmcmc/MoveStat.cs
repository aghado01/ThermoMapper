namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// Cumulative attempt/accept counts for one move type — the per-move acceptance diagnostic exposed by
/// <see cref="ReversibleJumpChain{TState}.MoveStats"/>. <c>Accepted / Attempts</c> is that move's acceptance
/// rate, the signal an adaptive proposal-scale tuner targets.
/// </summary>
public readonly record struct MoveStat(string Key, long Attempts, long Accepted);
