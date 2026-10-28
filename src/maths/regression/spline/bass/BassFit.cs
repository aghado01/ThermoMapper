using System.Collections.Generic;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A BASS fit (the contract shapes): the retained posterior and its account, immutable. Each chain's draws are kept in
/// draw order (<see cref="BassDraw"/>, every draw's generative state, D04), with structures as sorted ids into one term
/// table, so fitted values, means and intervals are views over the draws, computed on demand.
/// <list type="bullet">
/// <item><see cref="Terms"/> merges the chains' own tables in chain order, so ids never depend on how the chains were
/// scheduled.</item>
/// <item><see cref="DrawWeights"/> is null: MCMC's draws weigh equally (SMC's carry weights when it joins).</item>
/// <item><see cref="ResidualDependence"/> is the direct functional fit's residual-dependence check along the output
/// index (<see cref="BassResidualDependence"/>); the run leaves it null.</item>
/// <item><see cref="Checkpoint"/> is the run's restart state when <see cref="BassRunOptions.Checkpoint"/> asked for it,
/// beside the posterior, never inside it; <see cref="Bass.Resume"/> continues it.</item>
/// </list>
/// </summary>
public sealed record BassFit(
    BassConfig Config, BassResolved Resolved, BassInputTransform Transform, RunIdentity Identity,
    BassTermTable Terms,
    IReadOnlyList<DrawStore<BassDraw>> Draws,
    double[]? DrawWeights,
    BassDiagnostics Diagnostics,
    ResidualDependence? ResidualDependence,
    EnsembleCheckpoint? Checkpoint);
