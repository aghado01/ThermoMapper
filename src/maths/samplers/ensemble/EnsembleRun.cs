namespace Maths.Samplers.Ensemble;

/// <summary>
/// The pooled diagnostics of an ensemble run — everything the shell can manufacture without knowing what the
/// functionals <i>mean</i>. <see cref="FunctionalMean"/>/<see cref="FunctionalRHat"/> span all
/// <see cref="IEnsembleModel{TDraw}.FunctionalDim"/> functionals (R̂ raw — the client clamps for display if it
/// wants); <see cref="FunctionalEss"/> covers the leading <see cref="IEnsembleModel{TDraw}.EssDim"/>.
/// <see cref="ChainSeeds"/> is the run's reproducibility provenance (the resolved children of the master seed).
/// The client assembles its domain result from these plus the per-chain <see cref="IChainRun{TDraw}"/> handles
/// the shell hands back alongside this record.
/// </summary>
public sealed record EnsembleRun(
    double[] FunctionalMean,
    double[] FunctionalRHat,
    double[] FunctionalEss,
    int SamplesUsed,
    double AcceptanceRate,
    int[] ChainSeeds);
