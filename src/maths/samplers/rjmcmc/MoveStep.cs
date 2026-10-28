using Maths.Rng;
using Maths.Samplers.Blocks;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// A reversible-jump chain's own step as an <see cref="IBlock{TState}"/> over the chain: one weighted move on the chain's
/// stream, accepted when the chain accepts it. A <see cref="BlockScan{TState}"/> over a lone chain then composes its
/// moves with blocks that update the chain's state at its <see cref="ReversibleJumpChain{TState}.Temperature"/>, the
/// same blocks a <see cref="ParallelTempering{TState}"/> replica runs between its moves and the swaps. The scan's stream
/// goes to those blocks; this one ignores it.
/// </summary>
public sealed class MoveStep<TState> : IBlock<ReversibleJumpChain<TState>>
{
    /// <inheritdoc/>
    public string Name => "moves";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<TState> state, Xoshiro256PlusPlus rng)
    {
        long accepted = state.Accepted;
        state.Step();
        return state.Accepted > accepted;
    }
}
