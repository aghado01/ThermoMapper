using System;
using System.Collections.Generic;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
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
/// <item><see cref="Predict"/> (stage 5, D07), <see cref="Sensitivity"/> and <see cref="Effects"/> (stage 6, D08) are views
/// over the draws too, each on the data the fit was run on.</item>
/// <item><see cref="Traces"/>, <see cref="TrainingFit"/>, <see cref="Metadata"/>, <see cref="Summary"/> and
/// <see cref="SensitivitySummary"/> (stage 9) are the presentation data, values only: what R's <c>plot</c>, <c>print</c>
/// and <c>summary</c> show, with rendering and formatting left to adapters above maths.</item>
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

    /// <summary>
    /// Prediction at <paramref name="cases"/> (D07, stage 5; <see cref="BassPrediction"/>): latent and predictive draws as
    /// draw × case × location arrays in the requested draw order, the inputs through the fit's transform with
    /// extrapolation flagged, and new output locations for a direct functional fit. <paramref name="data"/> must be the
    /// data the fit was run on, whose knots and levels the terms index; data with another fingerprint are refused.
    /// </summary>
    public BassPrediction Predict(BassData data, PredictionCases cases, BassPredictOptions? options = null)
    {
        BassModel model = ModelOn(data);
        return BassDirectPrediction.Predict(model, Terms, Draws, DrawWeights, cases, options ?? new BassPredictOptions());
    }

    /// <summary>
    /// The Sobol decomposition of the draws (D08, stage 6; <see cref="BassSensitivity.Of"/>) under <paramref name="spec"/>'s
    /// product measure, the training range's uniform by default. <paramref name="data"/> must be the data the fit was run
    /// on, whose knots and levels the terms index and whose values the empirical measure reads; data with another
    /// fingerprint are refused.
    /// </summary>
    public SobolAnalysis Sensitivity(BassData data, SensitivitySpec? spec = null) => BassSensitivity.Of(this, data, spec);

    /// <summary>
    /// The main and two-way effect functions of the draws (D08, stage 6; <see cref="BassSensitivity.Effects"/>) under
    /// <paramref name="spec"/>'s measure, on the data the fit was run on; data with another fingerprint are refused.
    /// </summary>
    public EffectAnalysis Effects(BassData data, SensitivitySpec? spec = null, EffectSpec? effects = null)
        => BassSensitivity.Effects(this, data, spec, effects);

    /// <summary>Each chain's traces of M, σ², τ and λ in draw order (stage 9, P09; <see cref="BassTraces"/>).</summary>
    public BassTraces Traces() => BassTraces.Of(this);

    /// <summary>
    /// Fitted against observed at the training cases (stage 9, P10; <see cref="BassTrainingFit"/>): the posterior mean, the
    /// predictive mixture's quantiles at <paramref name="probabilities"/> (each in (0, 1);
    /// <see cref="PosteriorSummary.DefaultProbabilities"/> when null) and the residuals, through <see cref="Predict"/> on
    /// <paramref name="data"/>, the data the fit was run on; data with another fingerprint are refused.
    /// </summary>
    public BassTrainingFit TrainingFit(BassData data, IReadOnlyList<double>? probabilities = null)
        => BassTrainingFit.Of(this, data, probabilities);

    /// <summary>The fit's metadata (stage 9, P12; <see cref="BassMetadata"/>): the data's shape, the config and what it
    /// resolved to, the run's identity, lengths and account.</summary>
    public BassMetadata Metadata() => BassMetadata.Of(this);

    /// <summary>The metadata with the posterior summaries of M, σ², τ and λ at <paramref name="probabilities"/> (stage 9,
    /// P14; <see cref="BassSummary"/>), under the draw weights.</summary>
    public BassSummary Summary(IReadOnlyList<double>? probabilities = null) => BassSummary.Of(this, probabilities);

    /// <summary>
    /// The posterior summary of <see cref="Sensitivity"/>'s indices and variances at <paramref name="probabilities"/>
    /// (stage 9, P11; <see cref="SobolAnalysisSummary"/>), on the data the fit was run on.
    /// </summary>
    public SobolAnalysisSummary SensitivitySummary(BassData data, SensitivitySpec? spec = null, IReadOnlyList<double>? probabilities = null)
        => SobolAnalysisSummary.Of(this, Sensitivity(data, spec), probabilities);

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
