using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers.Evidence;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A tempered run's evidence ln p(y) (D03, D05): the data's marginal likelihood under BASS's priors restricted to F,
/// estimated over each chain's ladder from b = 1 down to its exact reference rung at b = 0, where the prior on F is
/// normalized. A run reports one only with the reference rung, under proper priors (neither the flat intercept nor the
/// Jeffreys σ² prior, whose constants leave it undefined), and when <see cref="BassSampling.ReportEvidence"/> asks.
/// <list type="bullet">
/// <item><see cref="SteppingStone"/> gives <see cref="LogEvidence"/>: it reads ℓ_b at the rungs only and has no
/// discretization bias, so it is the integrated form's estimator.</item>
/// <item><see cref="PathSampling"/> reads the same draws where the slope comes free: on the linear form, whose slope is ℓ₁
/// and whose ℓ_b is linear in b, so Friel, Hurn and Wyse's diagnostic holds and the two estimates' difference reads the
/// trapezoid's error. On the integrated form it is null: its slope at b = 0 carries −(k/2)·E[1/τ], infinite under
/// Zellner–Siow, so the trapezoid's end has no finite mean, and elsewhere its diagnostic misses E_b[∂²ℓ_b/∂b²] (D05's
/// first point).</item>
/// </list>
/// <see cref="Priors"/> says what the value depends on beyond the data.
/// </summary>
public sealed record BassEvidence(BassEvidenceEstimate SteppingStone, BassEvidenceEstimate? PathSampling, BassEvidencePriors Priors)
{
    /// <summary>The run's evidence, by stepping-stone.</summary>
    public double LogEvidence => SteppingStone.LogEvidence;

    /// <summary>Its Monte Carlo error.</summary>
    public double MonteCarloError => SteppingStone.MonteCarloError;

    /// <summary>
    /// The evidence over <paramref name="ladders"/>, one per chain, each traced since burn-in ended
    /// (<see cref="BassLadder.EndBurnIn"/>) down to its reference rung: stepping-stone, and path sampling on the linear form.
    /// </summary>
    internal static BassEvidence Of(BassModel model, IReadOnlyList<BassLadder> ladders)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(ladders);
        if (ladders.Count < 1) throw new ArgumentException("At least one ladder is required.", nameof(ladders));
        BassEvidenceEstimate? pathSampling = model.Target.Form == TemperedForm.Linear ? Pooled(ladders, LadderEstimator.PathSampling) : null;
        return new BassEvidence(Pooled(ladders, LadderEstimator.SteppingStone), pathSampling, BassEvidencePriors.Of(model.Config, model.Resolved));
    }

    // Each chain's own ladder gives its estimate (each adapted its own rungs); the run's is their mean, as BARS pools.
    private static BassEvidenceEstimate Pooled(IReadOnlyList<BassLadder> ladders, LadderEstimator method)
    {
        var chains = new LadderEstimate[ladders.Count];
        double sum = 0.0, variance = 0.0, discretization = 0.0;
        for (int c = 0; c < chains.Length; c++)
        {
            chains[c] = LadderEvidence.FromLadders([ladders[c].Tempering], method);
            sum += chains[c].LogRatio;
            variance += chains[c].MonteCarloError * chains[c].MonteCarloError;
            discretization += chains[c].DiscretizationError;
        }
        int n = chains.Length;
        return new BassEvidenceEstimate(method, sum / n, Math.Sqrt(variance) / n, discretization / n, chains);
    }
}

/// <summary>
/// One estimator's evidence over a run's ladders: each chain's own estimate (<see cref="LadderEvidence.FromLadders{TState}"/>,
/// over the rungs it froze), and their mean, with Monte Carlo error √Σ err_c² / C and the mean of their discretization
/// diagnostics (NaN for stepping-stone, which has none). The chains' spread checks the errors.
/// </summary>
public sealed record BassEvidenceEstimate(
    LadderEstimator Method, double LogEvidence, double MonteCarloError, double DiscretizationError, IReadOnlyList<LadderEstimate> Chains);

/// <summary>
/// What an evidence depends on beyond the data, carried with it because the marginal likelihood is sensitive to prior
/// scales in a way the posterior is not (D05's third point). Two evidences on the same data compare as a Bayes factor
/// only where these agree, or differ in what the comparison is about.
/// <list type="bullet">
/// <item><b>The coefficient prior:</b> <see cref="Coefficients"/> as declared (τ's shape a, or the held τ) and as
/// resolved, τ ~ Gamma(a, rate a·N) (<see cref="TauPrior"/>) or the held τ (<see cref="FixedTau"/>).</item>
/// <item><b>The count prior</b> on M, λ integrated out (D02), truncated at <see cref="MaxBasis"/>.</item>
/// <item><b>The term measure</b> (D01): the categorical subset measure and the interaction caps as resolved.</item>
/// <item><b>The support rule</b> as declared and its resolved thresholds, and D06's admission gates, which with the caps
/// define F, the support the prior is restricted to and normalized on.</item>
/// <item><b>The centre</b> and its resolved level, and <b>σ²'s prior</b> as declared and as resolved, IG(c, c·s²) under the
/// scaled default. Both read the response, so the default evidence is empirical-Bayes in the level and the scale.</item>
/// </list>
/// </summary>
public sealed record BassEvidencePriors(
    BassCoefficientPrior Coefficients, Gamma? TauPrior, double? FixedTau,
    BassCountPrior Count, int MaxBasis,
    CategoricalSubsetMeasure CategoricalSubsets, int MaxInteraction, int MaxCategoricalInteraction, int MaxFunctionalInteraction,
    SupportRule Support, double SupportMinimum, double? SupportMinimumFunctional, double AngleFloor, double? MaxVarianceInflation,
    GPriorCentre Centre, double CentreLevel, Sigma2Prior Sigma2Prior, InverseGamma Sigma2)
{
    /// <summary>The record for <paramref name="config"/> as <paramref name="resolved"/> on its data, which must carry a
    /// proper σ² prior.</summary>
    public static BassEvidencePriors Of(BassConfig config, BassResolved resolved)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(resolved);
        BassStructureSpec s = config.Structure;
        return new BassEvidencePriors(config.Coefficients, resolved.TauPrior, resolved.FixedTau, config.Count, resolved.MaxBasis,
            s.CategoricalSubsets, resolved.MaxInteraction, resolved.MaxCategoricalInteraction, resolved.MaxFunctionalInteraction,
            s.Support, resolved.SupportMinimum, resolved.SupportMinimumFunctional, s.AngleFloor, s.MaxVarianceInflation,
            config.Intercept, resolved.CentreLevel, config.Sigma2,
            resolved.Sigma2 ?? throw new ArgumentException("The Jeffreys σ² prior leaves the evidence undefined (D03).", nameof(resolved)));
    }
}
