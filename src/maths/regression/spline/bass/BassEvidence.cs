using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers.Evidence;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A tempered run's evidence ln p(y) (D03, D05): the data's marginal likelihood under BASS's priors restricted to F,
/// estimated over each chain's ladder from b = 1 down to its exact reference rung at b = 0, where the prior on F is
/// normalized, or over an SMC fit's stages from the prior to b = 1 (<see cref="OfSmc"/>), each replicate's normalizer. A
/// run reports one only with the reference rung or by SMC that reached b = 1, under proper priors (neither the flat
/// intercept nor the Jeffreys σ² prior, whose constants leave it undefined), and when
/// <see cref="BassSampling.ReportEvidence"/> asks.
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
/// <para>The evidence carries its flags, and a flagged one is still reported (decided 2026-10-28): the largest pair
/// rejection over the chains' sampling phases (<see cref="MaxPairRejection"/>, <see cref="Saturated"/> above
/// <see cref="RejectionThreshold"/>), where stepping-stone's importance weights grow heavy-tailed and its delta-method error
/// stops reading them; whether a self-sizing ladder reached its budget (<see cref="AtMaxLevels"/>), where the count the
/// barrier asked for may not have been taken; and the between-to-within ratio of the errors
/// (<see cref="BassEvidenceEstimate.BetweenToWithin"/>), well above 1 where the chains disagree beyond their own errors.</para>
/// </summary>
public sealed record BassEvidence(BassEvidenceEstimate SteppingStone, BassEvidenceEstimate? PathSampling, BassEvidencePriors Priors,
                                  double MaxPairRejection, bool AtMaxLevels)
{
    /// <summary>
    /// The pair rejection above which the evidence is flagged <see cref="Saturated"/>. Under a Gaussian reading of the swap,
    /// a pair rejecting r has stepping-stone weights of relative variance e^{σ²/2} − 1 with r = 1 − 2Φ(−σ/2): about 220 at
    /// 0.9 and 2,200 at 0.95, against a rung's few hundred to few thousand effective draws, so past 0.9 the weights' mean
    /// may not yet have seen its tail and the delta-method error may not read it. An equalized ladder of N ≈ 2Λ pairs
    /// rejects about half; S8's saturated ladder rejected 0.96–1.00.
    /// </summary>
    public const double RejectionThreshold = 0.9;

    /// <summary>The run's evidence, by stepping-stone.</summary>
    public double LogEvidence => SteppingStone.LogEvidence;

    /// <summary>Its Monte Carlo error: the larger of the within- and between-chain errors.</summary>
    public double MonteCarloError => SteppingStone.MonteCarloError;

    /// <summary>Whether a pair's rejection rose above <see cref="RejectionThreshold"/>.</summary>
    public bool Saturated => MaxPairRejection > RejectionThreshold;

    /// <summary>
    /// The evidence over <paramref name="ladders"/>, one per chain, each traced since burn-in ended
    /// (<see cref="BassLadder.EndBurnIn"/>) down to its reference rung, with the flags of their <paramref name="readout"/>:
    /// stepping-stone, and path sampling on the linear form.
    /// </summary>
    internal static BassEvidence Of(BassModel model, IReadOnlyList<BassLadder> ladders, BassTemperingReadout readout)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(ladders);
        ArgumentNullException.ThrowIfNull(readout);
        if (ladders.Count < 1) throw new ArgumentException("At least one ladder is required.", nameof(ladders));
        BassEvidenceEstimate? pathSampling = model.Target.Form == TemperedForm.Linear ? Pooled(ladders, LadderEstimator.PathSampling) : null;
        return new BassEvidence(Pooled(ladders, LadderEstimator.SteppingStone), pathSampling, BassEvidencePriors.Of(model.Config, model.Resolved),
                                readout.MaxPairRejection, readout.AtMaxLevels);
    }

    /// <summary>
    /// The evidence of an SMC fit (stage 8) from its <paramref name="replicates"/>, each run from the prior on F at b = 0 to
    /// b = 1. Each replicate's estimate is its normalizer ln Ẑ with the genealogy error (Lee and Whiteley 2018; NaN unless
    /// every resampling was multinomial), which is stepping-stone over its stages: with traced stages,
    /// <see cref="LadderEvidence.FromSmc{TState}"/> gives it with its rungs, and path sampling joins on the linear form,
    /// whose error reads the stages' weights and not their shared ancestry. The replicates pool by
    /// <see cref="LadderEvidence.Pool"/>'s rule, the larger of the within-run and between-run errors. There is no ladder,
    /// so the rejection flags read NaN and false.
    /// </summary>
    internal static BassEvidence OfSmc(BassModel model, IReadOnlyList<SmcResult<BassParticle>> replicates)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(replicates);
        if (replicates.Count < 1) throw new ArgumentException("At least one replicate is required.", nameof(replicates));
        bool traced = replicates[0].StageTraces is not null;
        var normalizers = new LadderEstimate[replicates.Count];
        for (int c = 0; c < normalizers.Length; c++)
            normalizers[c] = traced
                ? LadderEvidence.FromSmc(replicates[c])
                : new LadderEstimate(LadderEstimator.SteppingStone, replicates[c].LogNormalizerRatio, replicates[c].LogNormalizerError, double.NaN, []);
        BassEvidenceEstimate? pathSampling = null;
        if (traced && model.Target.Form == TemperedForm.Linear)
            pathSampling = Pooled([.. replicates.Select(r => LadderEvidence.FromSmc(r, LadderEstimator.PathSampling))]);
        return new BassEvidence(Pooled(normalizers), pathSampling, BassEvidencePriors.Of(model.Config, model.Resolved), double.NaN, false);
    }

    // Each chain's own ladder gives its estimate (each adapted its own rungs); the run's pools them as BARS's does.
    private static BassEvidenceEstimate Pooled(IReadOnlyList<BassLadder> ladders, LadderEstimator method)
    {
        var chains = new LadderEstimate[ladders.Count];
        for (int c = 0; c < chains.Length; c++) chains[c] = LadderEvidence.FromLadders([ladders[c].Tempering], method);
        return Pooled(chains);
    }

    // Independent runs' estimates by one estimator, pooled.
    private static BassEvidenceEstimate Pooled(IReadOnlyList<LadderEstimate> runs)
    {
        PooledLadderEstimate pooled = LadderEvidence.Pool(runs);
        return new BassEvidenceEstimate(pooled.Method, pooled.LogRatio, pooled.MonteCarloError, pooled.WithinChainError, pooled.BetweenChainError,
                                        pooled.DiscretizationError, pooled.Chains);
    }
}

/// <summary>
/// One estimator's evidence over a run's ladders: each chain's own estimate (<see cref="LadderEvidence.FromLadders{TState}"/>,
/// over the rungs it froze), and their mean (<see cref="LadderEvidence.Pool"/>), with the within-chain error √Σ err_c² / C,
/// the between-chain error from the chains' spread, the larger of the two as its Monte Carlo error (the within-chain error
/// alone for one chain), and the mean of their discretization diagnostics (NaN for stepping-stone, which has none).
/// </summary>
public sealed record BassEvidenceEstimate(
    LadderEstimator Method, double LogEvidence, double MonteCarloError, double WithinChainError, double BetweenChainError,
    double DiscretizationError, IReadOnlyList<LadderEstimate> Chains)
{
    /// <summary>The between-chain error over the within-chain one: well above 1 where the chains' own errors under-cover;
    /// NaN for one chain.</summary>
    public double BetweenToWithin => BetweenChainError / WithinChainError;
}

/// <summary>
/// What an evidence depends on beyond the data, carried with it because the marginal likelihood is sensitive to prior
/// scales in a way the posterior is not (D05's third point). Two evidences on the same data compare as a Bayes factor
/// only where these agree, or differ in what the comparison is about.
/// <list type="bullet">
/// <item><b>The coefficient prior:</b> <see cref="Coefficients"/> as declared, its covariance form and τ's prior, and τ's
/// prior as resolved on the form's unit-information precision τ₀ (<see cref="Scale"/>): 1/N for the g-form, and 1 for the
/// ridge, whose columns are the centred ones at unit root mean square over the N observations, the normalization its τ is
/// stated against.</item>
/// <item><b>The count prior</b> on M, λ integrated out (D02), truncated at <see cref="MaxBasis"/>.</item>
/// <item><b>The term measure</b> (D01): the categorical subset measure and the interaction caps as resolved.</item>
/// <item><b>The support rule</b> as declared and its resolved thresholds, and D06's admission gates, which with the caps
/// define F, the support the prior is restricted to and normalized on.</item>
/// <item><b>The centre</b> and its resolved level, and <b>σ²'s prior</b> as declared and as resolved, IG(c, c·s²) under the
/// scaled default. Both read the response, so the default evidence is empirical-Bayes in the level and the scale.</item>
/// </list>
/// </summary>
public sealed record BassEvidencePriors(
    BassCoefficientPrior Coefficients, ResolvedTauPrior Scale,
    BassCountPrior Count, int MaxBasis,
    CategoricalSubsetMeasure CategoricalSubsets, int MaxInteraction, int MaxCategoricalInteraction, int MaxFunctionalInteraction,
    SupportRule Support, double SupportMinimum, double? SupportMinimumFunctional, double AngleFloor, double? MaxVarianceInflation,
    GPriorCentre Centre, double CentreLevel, Sigma2Prior Sigma2Prior, InverseGamma Sigma2)
{
    /// <summary>The record for <paramref name="config"/> as <paramref name="resolved"/> on its data, which must carry
    /// proper σ² and τ priors.</summary>
    public static BassEvidencePriors Of(BassConfig config, BassResolved resolved)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(resolved);
        BassStructureSpec s = config.Structure;
        if (!resolved.Scale.IsProper)
            throw new ArgumentException($"τ's prior, {resolved.Scale.Declared}, is improper, which leaves the evidence undefined.", nameof(resolved));
        return new BassEvidencePriors(config.Coefficients, resolved.Scale, config.Count, resolved.MaxBasis,
            s.CategoricalSubsets, resolved.MaxInteraction, resolved.MaxCategoricalInteraction, resolved.MaxFunctionalInteraction,
            s.Support, resolved.SupportMinimum, resolved.SupportMinimumFunctional, s.AngleFloor, s.MaxVarianceInflation,
            config.Intercept, resolved.CentreLevel, config.Sigma2,
            resolved.Sigma2 ?? throw new ArgumentException("The Jeffreys σ² prior leaves the evidence undefined (D03).", nameof(resolved)));
    }
}
