using System;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// The client side of the kernel/shell split: a model knows how to <b>start one independent replica</b> (its
/// kernel, its start policy — overdispersed for a rugged target, warm-started for a benign one) and how to
/// <b>reduce a draw to the scalar functionals</b> the shell monitors for convergence (R̂/ESS run on invariant
/// functionals, never on the raw — possibly trans-dimensional or label-degenerate — state). <see cref="ChainEnsemble"/>
/// owns the fan-out, rounds, pooling, and diagnostics; the model owns the semantics.
/// </summary>
/// <typeparam name="TDraw">The per-step carrier (e.g. a knot configuration, a bandwidth).</typeparam>
public interface IEnsembleModel<TDraw>
{
    /// <summary>Width of the functional vector each <see cref="IChainRun{TDraw}.Accumulate"/> populates — the scalars R̂ (and means) are computed on.</summary>
    int FunctionalDim { get; }

    /// <summary>
    /// How many leading functionals are monitored: their per-chain draws are kept, so they get the rank-normalized
    /// diagnostics (<see cref="ChainDiagnostics.Diagnose"/>: split R̂, bulk and tail ESS), and an R̂ stopping target
    /// waits on their effective sample size. The rest get classic R̂ from streaming moments only, so a model with a
    /// wide map of functionals (BARS's per-grid fit) keeps the monitored set small and puts it first (BARS: knot
    /// count at index 0). 0 ≤ MonitoredDim ≤ <see cref="FunctionalDim"/>.
    /// </summary>
    int MonitoredDim { get; }

    /// <summary>
    /// Build one replica + its sink from the two decorrelated seeds the shell derived for this chain.
    /// <paramref name="chainSeed"/> drives the transition kernel (and, for a tempered client, seeds its sub-ladder
    /// — hence the raw seed, not a pre-built RNG); <paramref name="readoutSeed"/> is the decoupled stream a client
    /// uses for readout-side randomness (BARS draws spline coefficients per sample), ignored by clients that need none.
    /// </summary>
    IChainRun<TDraw> StartChain(int chainSeed, int readoutSeed);

    /// <summary>
    /// A fingerprint of the configuration and data a run depends on, checked when a run resumes from a checkpoint so
    /// that it continues the same problem (<see cref="RunFingerprint"/>). Empty when the model does not resume.
    /// </summary>
    string Fingerprint => string.Empty;

    /// <summary>
    /// Rebuild one chain from the restart state its <see cref="IResumableChainRun{TDraw}.SaveState"/> gave, so that it
    /// continues exactly where it left off. The default refuses: the model does not resume.
    /// </summary>
    IChainRun<TDraw> RestoreChain(IChainState state)
        => throw new NotSupportedException($"{GetType().Name} does not resume from checkpoints.");
}

/// <summary>
/// One replica plus its thread-local accumulators: the kernel (<see cref="Chain"/>), the burn policy, and the
/// draw→functional reduction. The shell drives <see cref="IChain{TDraw}.Step"/> via <see cref="Chain"/> and folds
/// each draw through <see cref="Accumulate"/>; the model's own rich sinks (BARS's peak/span fields, BGP's predicted
/// means) accumulate inside the implementation and are read back by the client off the returned handles.
/// </summary>
public interface IChainRun<TDraw>
{
    /// <summary>The kernel the shell advances. Its acceptance counts feed the pooled acceptance rate.</summary>
    IChain<TDraw> Chain { get; }

    /// <summary>
    /// Run the burn-in (no accumulation). A hook, not just <c>Chain.Step()×n</c>, so a model can tune its proposal
    /// during burn and freeze it before sampling (BARS's adaptive-τ schedule); the default is plain stepping.
    /// </summary>
    void Burn(int steps);

    /// <summary>
    /// Reduce one post-burn <paramref name="draw"/> to the shell's functionals (must populate all
    /// <see cref="IEnsembleModel{TDraw}.FunctionalDim"/> entries) and fold the model's own sinks. Called once per
    /// sampling step, immediately after the <see cref="Chain"/> produced <paramref name="draw"/>.
    /// </summary>
    void Accumulate(in TDraw draw, Span<double> functionals);
}
