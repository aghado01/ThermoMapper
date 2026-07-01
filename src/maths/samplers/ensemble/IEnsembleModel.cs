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
    /// How many leading functionals also get ESS (their per-chain sequences are retained). ESS is O(n²) per
    /// functional, so a model with a wide R̂ map (e.g. BARS's per-grid consensus) keeps this small and orders the
    /// ESS-worthy functionals first (BARS: knot count at index 0). 0 ≤ EssDim ≤ <see cref="FunctionalDim"/>.
    /// </summary>
    int EssDim { get; }

    /// <summary>
    /// Build one replica + its sink from the two decorrelated seeds the shell derived for this chain.
    /// <paramref name="chainSeed"/> drives the transition kernel (and, for a tempered client, seeds its sub-ladder
    /// — hence the raw seed, not a pre-built RNG); <paramref name="readoutSeed"/> is the decoupled stream a client
    /// uses for readout-side randomness (BARS draws spline coefficients per sample), ignored by clients that need none.
    /// </summary>
    IChainRun<TDraw> StartChain(int chainSeed, int readoutSeed);
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
