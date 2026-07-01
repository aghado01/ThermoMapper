namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// One move's proposed transition: the candidate <typeparamref name="TState"/> plus the two log-terms the
/// acceptance ratio needs <i>from the move</i> (the engine supplies the move-selection and target ratios).
/// <c>LogProposalRatio</c> is the within-move auxiliary-density ratio <c>log[q_{m'}(u') / q_m(u)]</c> and
/// <c>LogJacobian</c> is the log dimension-matching Jacobian (Green 1995, §3.3). A move that is unavailable from
/// the current state returns <see langword="null"/> (a <c>Proposal?</c>) rather than a value.
/// </summary>
public readonly record struct Proposal<TState>(TState Candidate, double LogProposalRatio, double LogJacobian);
