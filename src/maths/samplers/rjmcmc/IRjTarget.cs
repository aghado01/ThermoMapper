namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// The (unnormalized) distribution the chain samples, reduced to the one quantity reversible-jump needs from
/// the model: the log target density at a carrier state. For a marginalized BARS chain this is
/// <c>log marginal-likelihood + log-prior</c> with the coefficients integrated out; trans-dimensionality is
/// invisible here because the move carries the Jacobian. The value is needed only up to an additive constant
/// that is the same for every state — it cancels in the acceptance ratio (Green 1995, Remark 3: one shared
/// unknown constant across subspaces suffices). The chain caches the current state's value, so each step
/// evaluates this once (the candidate's), not twice.
/// </summary>
public interface IRjTarget<TState>
{
    /// <summary>Log target density at <paramref name="state"/> (up to a state-independent additive constant).</summary>
    double LogPosterior(TState state);
}
