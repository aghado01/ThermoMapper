using System;
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
/// <item><see cref="Evidence"/> is the run's ln p(y) with its Monte Carlo error and the priors it depends on, from a run
/// tempered down to an exact reference rung under proper priors (<see cref="BassEvidence"/>); null otherwise.</item>
/// <item><see cref="ResidualDependence"/> is the direct functional fit's residual-dependence check along its one output
/// variable (<see cref="BassResidualDependence"/>, F11) under the posterior mean, computed at the end of the run; null for a
/// scalar response and for several output variables, where <see cref="ResidualDependenceFor"/> takes a declared
/// index.</item>
/// <item><see cref="PointwiseLogDensity"/> is OQ17's matrix over the draws (S9), from the data the fit was run on.</item>
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
    BassEvidence? Evidence,
    ResidualDependence? ResidualDependence,
    EnsembleCheckpoint? Checkpoint)
{
    /// <summary>
    /// The draws × units log density (OQ17; <see cref="BassPointwise.LogDensity"/>) on <paramref name="data"/>, which
    /// must be the data the fit was run on: the fit holds the config and the draws, not the data, and the terms' knot and
    /// level indices refer to that data's inputs, so data with another fingerprint are refused.
    /// </summary>
    public PointwiseLogDensity PointwiseLogDensity(BassData data, UnitPartition? units = null,
                                                   PointwiseContent content = PointwiseContent.RetainedParameters)
    {
        BassModel model = ModelOn(data);
        return BassPointwise.LogDensity(model, Terms, Draws, units, content);
    }

    /// <summary>
    /// The residual-dependence check (<see cref="BassResidualDependence.Of"/>) along a declared <paramref name="index"/>
    /// on <paramref name="data"/>, which must be the data the fit was run on, summarized at
    /// <paramref name="probabilities"/>: for several output variables, or another index than the default's.
    /// </summary>
    public ResidualDependence ResidualDependenceFor(BassData data, ResidualIndex index, IReadOnlyList<double>? probabilities = null)
    {
        ArgumentNullException.ThrowIfNull(index);
        BassModel model = ModelOn(data);
        return BassResidualDependence.Of(model, Terms, Draws, index, DrawWeights, probabilities)!;
    }

    // The model on `data`, refused when they are not the data the fit was run on.
    private BassModel ModelOn(BassData data)
    {
        ArgumentNullException.ThrowIfNull(data);
        string fingerprint = Bass.Fingerprint(Config, data);
        if (fingerprint != Identity.Fingerprint)
            throw new ArgumentException($"The fit was run on other data: its fingerprint is {Identity.Fingerprint}, theirs {fingerprint}.", nameof(data));
        return BassModel.Create(Config, data);
    }
}
