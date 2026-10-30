using System;
using System.Collections.Generic;
using System.Globalization;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Samplers.FlatHistogram;

/// <summary>Reduce one draw to the functionals a run monitors and reweighs, filling every entry.</summary>
public delegate void DrawReduction<TState>(in TState draw, Span<double> functionals);

/// <summary>
/// The kernel a <see cref="MulticanonicalRun"/> drives: any chain that samples its <see cref="MulticanonicalTarget{TState}"/>
/// at <see cref="Temperature.Posterior"/>, and that exposes its current state's evaluation and takes it back reweighed
/// when the learning moves the weights. The run steps every chain, learns from their visits during burn-in, then hands
/// them to the ensemble shell for production on the frozen weights.
/// </summary>
public interface IMulticanonicalChain<TState> : IChain<TState>
{
    /// <summary>The current state's cached evaluation, under the weights it was last reweighed to.</summary>
    MulticanonicalEvaluation CurrentEvaluation { get; }

    /// <summary>The target's weights have moved: re-read the current state's ln w under them from its cached
    /// evaluation (<see cref="MulticanonicalEvaluation.Reweighed"/>), without evaluating it again.</summary>
    void Reweigh();

    /// <summary>Reduce a production draw to the model's functionals (all of them) and fold the chain's own sinks.
    /// Called once per production step, right after <see cref="IChain{TDraw}.Step"/>; never during learning.</summary>
    void Accumulate(in TState draw, Span<double> functionals);

    /// <summary>This chain's restart state, independent of the live chain; the default refuses.</summary>
    IChainState SaveState() => throw new NotSupportedException($"{GetType().Name} cannot save its restart state.");
}

/// <summary>
/// The client side of a <see cref="MulticanonicalRun"/>: how to start one chain on the multicanonical target, which
/// functionals its draws reduce to, and, to resume, how to rebuild a chain from its saved state. The analogue of
/// <see cref="IEnsembleModel{TDraw}"/>, whose declarations the production phase hands to the shell.
/// </summary>
public interface IMulticanonicalModel<TState>
{
    /// <summary>The width of the functional vector each draw reduces to.</summary>
    int FunctionalDim { get; }

    /// <summary>How many leading functionals the shell monitors with rank-normalized diagnostics.</summary>
    int MonitoredDim { get; }

    /// <summary>Each functional's name; "f0", "f1", … by default.</summary>
    IReadOnlyList<string> FunctionalNames
    {
        get
        {
            var names = new string[FunctionalDim];
            for (int f = 0; f < names.Length; f++) names[f] = "f" + f.ToString(CultureInfo.InvariantCulture);
            return names;
        }
    }

    /// <summary>The groups D10's verdict reads together; one group of every functional by default.</summary>
    IReadOnlyList<FunctionalGroup> Groups => FunctionalDim == 0 ? [] : [FunctionalGroup.All(FunctionalDim)];

    /// <summary>
    /// Build one chain on <paramref name="target"/>, which it must sample at <see cref="Temperature.Posterior"/>, from
    /// the two decorrelated seeds the run derived for it, as the shell derives them (<see cref="IEnsembleModel{TDraw}.StartChain"/>).
    /// The target's weights change during learning; the chain reads them through the target.
    /// </summary>
    IMulticanonicalChain<TState> StartChain(MulticanonicalTarget<TState> target, int chainSeed, int readoutSeed);

    /// <summary>A fingerprint of the configuration and data a run depends on, checked on resume; empty by default.</summary>
    string Fingerprint => string.Empty;

    /// <summary>Rebuild a chain on <paramref name="target"/> from its <see cref="IMulticanonicalChain{TState}.SaveState"/>;
    /// the default refuses.</summary>
    IMulticanonicalChain<TState> RestoreChain(MulticanonicalTarget<TState> target, IChainState state)
        => throw new NotSupportedException($"{GetType().Name} does not resume from checkpoints.");
}

/// <summary>A <see cref="MulticanonicalChain{TState}"/>'s restart state: its reversible-jump chain's snapshot.</summary>
public sealed record MulticanonicalChainState<TState>(ChainSnapshot<TState> Chain) : IChainState;

/// <summary>
/// A reversible-jump chain on a <see cref="MulticanonicalTarget{TState}"/> as an <see cref="IMulticanonicalChain{TState}"/>:
/// its moves step it, its cached evaluation is reweighed by a swap-free <see cref="ReversibleJumpChain{TState}.SetState"/>
/// on the same state, and a reduction gives its functionals. For kernels with blocks or other state, a model writes its
/// own <see cref="IMulticanonicalChain{TState}"/>.
/// </summary>
public sealed class MulticanonicalChain<TState> : IMulticanonicalChain<TState>
{
    private readonly ReversibleJumpChain<TState> _chain;
    private readonly MulticanonicalTarget<TState> _target;
    private readonly DrawReduction<TState> _reduce;

    /// <param name="chain">A chain over <paramref name="target"/> at <see cref="Temperature.Posterior"/>.</param>
    /// <param name="target">The multicanonical target the chain samples.</param>
    /// <param name="reduce">The draw's functionals.</param>
    public MulticanonicalChain(ReversibleJumpChain<TState> chain, MulticanonicalTarget<TState> target, DrawReduction<TState> reduce)
    {
        ArgumentNullException.ThrowIfNull(chain);
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(reduce);
        if (chain.Temperature != Temperature.Posterior)
            throw new ArgumentException("A multicanonical chain samples π_w at b = 1.", nameof(chain));
        if (chain.CurrentEvaluation is not MulticanonicalEvaluation)
            throw new ArgumentException("The chain does not sample a multicanonical target.", nameof(chain));
        _chain = chain;
        _target = target;
        _reduce = reduce;
    }

    /// <summary>The reversible-jump chain.</summary>
    public ReversibleJumpChain<TState> Chain => _chain;

    public long Accepted => _chain.Accepted;
    public long Attempts => _chain.Attempts;

    public TState Step() => _chain.Step();

    public MulticanonicalEvaluation CurrentEvaluation => (MulticanonicalEvaluation)_chain.CurrentEvaluation;

    public void Reweigh() => _chain.SetState(_chain.Current, CurrentEvaluation.Reweighed(_target.Weights));

    public void Accumulate(in TState draw, Span<double> functionals) => _reduce(draw, functionals);

    public IChainState SaveState() => new MulticanonicalChainState<TState>(_chain.Snapshot());

    /// <summary>A chain that continues exactly where <paramref name="state"/>'s left off, over the moves it ran with,
    /// on <paramref name="target"/> (whose weights the run restores).</summary>
    public static MulticanonicalChain<TState> Restore(IReadOnlyList<IRjMove<TState>> moves, MulticanonicalTarget<TState> target,
                                                      IChainState state, DrawReduction<TState> reduce)
    {
        if (state is not MulticanonicalChainState<TState> saved)
            throw new ArgumentException($"The state is not a {nameof(MulticanonicalChainState<TState>)}.", nameof(state));
        return new MulticanonicalChain<TState>(ReversibleJumpChain<TState>.Restore(moves, target, saved.Chain), target, reduce);
    }
}
