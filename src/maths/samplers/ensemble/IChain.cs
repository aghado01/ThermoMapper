namespace Maths.Samplers.Ensemble;

/// <summary>
/// Anything the ensemble shell can drive sequentially as one independent replica: a reversible-jump chain, an
/// MH-over-log-t walk, a NUTS trajectory, a tempered ladder (whose <see cref="Step"/> returns the cold state).
/// The kernel — <i>how</i> the state is traversed — lives behind this seam; the shell only advances it and reads
/// its acceptance bookkeeping, so the same orchestrator serves every member of the engine family.
/// </summary>
/// <typeparam name="TDraw">The per-step draw the shell reduces to functionals (the carrier, or the cold state of a ladder).</typeparam>
public interface IChain<out TDraw>
{
    /// <summary>Advance the chain one step and return the resulting draw.</summary>
    TDraw Step();

    /// <summary>Proposals accepted so far (for a tempered ladder, the cold replica's). <c>Accepted/Attempts</c> is the acceptance rate.</summary>
    long Accepted { get; }

    /// <summary>Proposals attempted so far (for a tempered ladder, the cold replica's).</summary>
    long Attempts { get; }
}
