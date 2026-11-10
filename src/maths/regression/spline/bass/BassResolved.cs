using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;

namespace Maths.Regression.Spline.Bass;

/// <summary>One auto-pick a fit made at entry: the setting, the value it took, and why (CONTRIBUTING's automation log).</summary>
public sealed record BassAutoPick(string Setting, string Value, string Reason);

/// <summary>
/// What a config resolves to on its data, at entry: the counts, the caps and thresholds left to the data, and the
/// priors read from the response. With <see cref="AutoPicks"/>, the log of every auto-pick and its reason, it is the
/// fit's account of what it decided beyond the config. The data-dependent refusals are made here too: a response the
/// intercept convention cannot fit (D03), and a support threshold no term could meet.
/// </summary>
public sealed record BassResolved
{
    /// <summary>n_x, the cases.</summary>
    public required int Cases { get; init; }

    /// <summary>n_r, the output locations; 1 for a scalar response.</summary>
    public required int Locations { get; init; }

    /// <summary>N = n_x · n_r, the observations τ's unit information counts (D03).</summary>
    public required int Observations { get; init; }

    /// <summary>The cap on M.</summary>
    public required int MaxBasis { get; init; }

    /// <summary>The caps on a term's continuous, categorical and functional factors.</summary>
    public required int MaxInteraction { get; init; }

    /// <summary>The cap on a term's categorical factors.</summary>
    public required int MaxCategoricalInteraction { get; init; }

    /// <summary>The cap on a term's functional factors; 0 for a scalar response.</summary>
    public required int MaxFunctionalInteraction { get; init; }

    /// <summary>The support threshold on a term's design part.</summary>
    public required double SupportMinimum { get; init; }

    /// <summary>The support threshold on a term's functional part; null for a scalar response.</summary>
    public required double? SupportMinimumFunctional { get; init; }

    /// <summary>The level the coefficient prior is centred on: the response's mean, or zero under the zero centre, or a frozen
    /// reference's level. Under weights the mean is under W only when the coefficient prior reads the weighted columns
    /// (<see cref="WeightedPriorGram"/>).</summary>
    public required double CentreLevel { get; init; }

    /// <summary>The frozen level level centring reads (<see cref="BassPriorReference.Level"/>); null otherwise, ȳ having no energy.</summary>
    public required double? PriorLevel { get; init; }

    /// <summary>κ in the coefficient prior's covariance σ²κ(τG)⁻¹ (<see cref="BassPriorReference.CoefficientConstant"/>); 1 otherwise.</summary>
    public required double CoefficientConstant { get; init; }

    /// <summary>ȳ, the response's mean over all N observations under W, whatever the centre; the zero centre's level energy
    /// Σwᵢȳ² reads it.</summary>
    public required double ResponseMean { get; init; }

    /// <summary>s², the response's mean square about <see cref="CentreLevel"/>, which the scaled σ² prior reads: under W only
    /// when the coefficient prior reads the weighted columns (<see cref="WeightedPriorGram"/>), and a frozen reference's when it
    /// states one.</summary>
    public required double ResponseMeanSquare { get; init; }

    /// <summary>The response's mean square under W about <see cref="CentreLevel"/>, whose ‖·‖²_W the exact-fit floor reads;
    /// <see cref="ResponseMeanSquare"/> at unit weights and under a frozen reference's.</summary>
    public required double WeightedMeanSquare { get; init; }

    /// <summary>The count, the mass and Σ ln wᵢ the marginal reads beyond the weighted statistics; (N, N, 0) at unit weights.</summary>
    public required GPriorWeights Weighting { get; init; }

    /// <summary>σ²'s prior as resolved; null for the Jeffreys prior and for a held σ². Under a hierarchical prior, the inverse gamma
    /// at its rate's mean.</summary>
    public required InverseGamma? Sigma2 { get; init; }

    /// <summary>σ²'s rate under a hierarchical prior, on s²; null under every other.</summary>
    public required ResolvedSigma2Rate? Sigma2Rate { get; init; }

    /// <summary>The held σ² (<see cref="Sigma2PriorKind.Held"/>), which the structure step's marginal is given; null when σ² is
    /// integrated out.</summary>
    public required double? GivenSigma2 { get; init; }

    /// <summary>τ's prior on the coefficient form's unit-information precision τ₀ (1/N for the g-form, 1 for the ridge):
    /// sampled, or held at its declared value or τ₀.</summary>
    public required ResolvedTauPrior Scale { get; init; }

    /// <summary>The training points whose fitted values the diagnostics monitor: the config's, a direct response's
    /// missing locations at the middle one, or the defaults.</summary>
    public required IReadOnlyList<BassMonitoredPoint> MonitoredPoints { get; init; }

    /// <summary>The groups an R̂ stop waits on, by name (<see cref="BassSampling.StoppingScope"/>): every group, or the
    /// predictive group alone. The verdict reports every group either way.</summary>
    public required IReadOnlyList<string> StoppingGroups { get; init; }

    /// <summary>
    /// A tempered run's hot end, its lowest positive b: declared, or auto-picked (<see cref="AutoHotEnd"/>); null when the
    /// run is untempered.
    /// </summary>
    public required double? HotEnd { get; init; }

    /// <summary>
    /// The local steps' windows w = max(1, ⌈ω·n_v⌉) over each hinge variable's n_v knots, continuous and functional
    /// (<see cref="LocalSteps.Window"/>); null when the steps are off.
    /// </summary>
    public required BassStepWindows? StepWindows { get; init; }

    /// <summary>The fixed-M factor moves' shares of the change's slot; null when the family is off.</summary>
    public required BassFactorShares? FactorShares { get; init; }

    /// <summary>A tempered run's split threshold (<see cref="SplitRule"/>); null when the rule is off or the run untempered.</summary>
    public required double? SplitThreshold { get; init; }

    /// <summary>Every auto-pick, with its reason.</summary>
    public required IReadOnlyList<BassAutoPick> AutoPicks { get; init; }

    /// <summary>Several outputs' resolution (E34): their means, centre levels, S_Y and Σ's prior; null for one response. With several
    /// outputs the scalar response's pieces (<see cref="ResponseMean"/>, <see cref="CentreLevel"/>, the mean squares) are NaN and
    /// <see cref="Sigma2"/> is null, since nothing reads them.</summary>
    public BassResolvedOutputs? Outputs { get; init; }

    /// <summary>
    /// <paramref name="config"/> resolved on <paramref name="data"/>, after <see cref="BassConfig.Validate"/>. Refused:
    /// a zero response under the zero centre, and a constant one under level centring or the flat intercept, which
    /// leave no spread to fit (D03); a support threshold above its part's size, which no term could meet; a ladder
    /// with a rung where its tempered form is undefined (D05), its hot end, declared or auto-picked, or the reference
    /// rung's b = 0, and an auto-picked hot end outside (0, 1), which tiny data give; and SMC on the integrated form under
    /// the improper options, which has no b = 0 to start from. The response the intercept convention fits and the scaled σ²
    /// prior reads is <paramref name="working"/>'s (F-X1); null is the response itself. A held σ² needs no spread, so it lifts
    /// the first two refusals. Under a frozen reference (<see cref="BassPriorReference"/>) the level and s² are the reference's, and
    /// those two refusals read its spread when it states one: a working response that happens to be constant is no reason to refuse.
    /// Refused too: a reference that leaves out what the config reads, and a level the zero centre contradicts.
    /// </summary>
    public static BassResolved Resolve(BassConfig config, BassData data, WorkingResponse? working = null)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        config.Validate();
        var log = new List<BassAutoPick>();
        BassResponse response = data.Response;
        int cases = response.Cases, locations = response.Locations, n = response.Observations;
        bool functional = response is FunctionalResponse;
        log.Add(new("Observations", Text(n), functional ? $"N = n_x · n_r = {cases} · {locations}" : "N = n_x, a scalar response"));
        if (!config.HasRandomEffects && data.Groupings.Count > 0)
            throw new ArgumentException("The data hold groupings, which only random effects read; declare RandomIntercepts or RandomSlopes, or leave the groupings out.", nameof(data));
        if (config.HasRandomEffects)
        {
            IBassRandomEffect[] blocks = config.EffectSpecifications;
            if (response is not ScalarResponse) throw new ArgumentException("Random effects currently serve a scalar Gaussian response.", nameof(data));
            RandomEffectsBlock effects = working?.Whitening ?? RandomEffectsBlock.Of(data.Groupings, blocks, data.Predictors.OfType<ContinuousPredictor>().ToArray());
            foreach (CategoricalPredictor grouping in data.Groupings)
                if (!blocks.Any(b => b.Grouping == grouping.Name || b.ParentGrouping == grouping.Name))
                    throw new ArgumentException($"Grouping '{grouping.Name}' is not read by any random-effect block or parent convention.", nameof(data));
            for (int k = 0; k < effects.Count; k++)
            {
                IBassRandomEffect block = blocks[k];
                RandomEffectGrouping grouping = effects.Layout(k);
                double[] observed = grouping.Counts.ToArray();
                string convention = grouping.ParentName is { } parent
                    ? $"centred within '{parent}', Σ_{{g in p}} n_g u_g = 0 for each parent level"
                    : "observation-centred, Σ n_g u_g = 0";
                string scaleText = block.ScalePrior is { } scalePrior ? $"s drawn from C⁺(0, {Text(scalePrior.Scale)}), from {Text(block.Scale)}" : $"s = {Text(block.Scale)}";
                string field = block is BassRandomSlope slope ? $"RandomSlopes[{slope.Grouping}:{slope.Input}]" : $"RandomIntercepts[{block.Grouping}]";
                if (block is BassRandomSlope slopeBlock) convention = $"{slopeBlock.InputScale} '{slopeBlock.Input}', observation-centred, Σ z_i b_g = Σ z_i² b_g = 0";
                log.Add(new(field, $"{observed.Length} levels, {scaleText}",
                    $"{convention}, {Text(observed.Min())} to {Text(observed.Max())} cases a level; jointly integrated by whitening the response and every column"));
            }
            log.Add(new("RandomEffectsRank", Text(effects.Rank), "rank of the constrained joint case design; scale updates factor only this reduced covariance"));
            foreach (RandomEffectAlias alias in effects.Aliases)
                log.Add(new("RandomEffectsAlias", $"'{alias.FirstGrouping}:{alias.FirstLevel}' and '{alias.SecondGrouping}:{alias.SecondLevel}'",
                    "identical case sets; separate proper variance components retain the latent uncertainty, with no additional constraint"));
        }
        IEnumerable<ContinuousPredictor> columns = data.Predictors.OfType<ContinuousPredictor>();
        if (response is FunctionalResponse withOutputs) columns = columns.Concat(withOutputs.OutputVariables);
        foreach (ContinuousPredictor c in columns)
            if (c.Domain is { } domain)
                log.Add(new($"Domain[{c.Name}]", domain.ToString(), "declared: the input scales by it onto [0, 1] in place of its training range"));

        // The response the intercept convention can fit, its level and its mean square about it.
        WorkingResponse fitted = working ?? new WorkingResponse(response);
        ReadOnlySpan<double> y = fitted.Values;
        GPriorWeights weights = fitted.Weighting;
        BassPriorReference? reference = fitted.Reference;
        bool held = config.Sigma2.IsHeld;
        if (!fitted.IsUnweighted)
            log.Add(new("Weights", $"{fitted.Meaning}, prior on the {fitted.PriorGram} columns",
                        $"Σwᵢ = {Text(weights.Mass)} over N = {n}; σ²'s shape counts {Text(weights.Count)}"));
        double mean, level, meanSquare, weightedMeanSquare;
        double? priorLevel = null;
        InverseGamma? sigma2;
        ResolvedSigma2Rate? rate;
        BassResolvedOutputs? several = null;
        if (response is MultiOutputResponse multi)
        {
            // Several outputs (E34): each output's level, S_Y and Σ's prior; the scalar σ² pieces are not read.
            several = ResolveOutputs(config, multi, fitted, log);
            mean = level = meanSquare = weightedMeanSquare = double.NaN;
            sigma2 = null;
            rate = null;
        }
        else
        {
            if (config.OutputCovariance is not null)
                throw new ArgumentException("OutputCovariance is Σ's prior for several outputs (E34); a single response reads Sigma2.", "config");
            bool constant = true, zero = y[0] == 0.0;
            for (int i = 1; i < y.Length && constant; i++) constant = y[i] == y[0];
            mean = fitted.ResponseMean;
            if (config.Intercept == GPriorCentre.Zero)
            {
                if (reference?.Level is double stated && stated != 0.0) throw new ArgumentException($"The reference's level is {Text(stated)}, and the zero centre's is 0.", nameof(working));
                if ((reference?.MeanSquare is double about ? about == 0.0 : constant && zero) && !held) throw new ArgumentException("The response is zero, which leaves σ no posterior (D03).", nameof(data));
                if (constant && !config.Sigma2.IsProper)
                    throw new ArgumentException("Under the Jeffreys σ² prior a constant response is fitted exactly by the intercept, which leaves the posterior improper (D03).", nameof(data));
                level = 0.0;
                log.Add(new("CentreLevel", "0", "the zero centre"));
            }
            else
            {
                if ((reference?.MeanSquare is double about ? about == 0.0 : constant) && !held) throw new ArgumentException("The response is constant, which the intercept absorbs whole, leaving no spread to fit (D03).", nameof(data));
                if (constant && !config.Sigma2.IsProper)
                    throw new ArgumentException("Under the Jeffreys σ² prior a constant response is fitted exactly by the intercept, which leaves the posterior improper (D03).", nameof(data));
                if (config.Intercept == GPriorCentre.Mean && reference is not null)
                {
                    level = reference.Level ?? throw new ArgumentException("The reference states no level, which level centring reads.", nameof(working));
                    priorLevel = level;
                    log.Add(new("CentreLevel", Text(level), "the reference's, frozen apart from the working response"));
                }
                else
                {
                    level = config.Intercept == GPriorCentre.Mean && fitted.PriorStandsApart ? PairwiseSum.Of(y) / n : mean;
                    log.Add(new("CentreLevel", Text(level), config.Intercept == GPriorCentre.Flat ? "the response's mean, the flat intercept's estimate" : "the response's mean, level centring"));
                }
            }
            var squares = new double[n];
            for (int i = 0; i < n; i++) squares[i] = (y[i] - level) * (y[i] - level);
            meanSquare = PairwiseSum.Of(squares) / n;
            weightedMeanSquare = meanSquare;
            if (fitted.WeightArray is { } w)
            {
                weightedMeanSquare = PairwiseSum.Dot(w, squares) / weights.Mass;
                if (!fitted.PriorStandsApart) meanSquare = weightedMeanSquare;
            }
            if (reference?.MeanSquare is double frozenSquare) meanSquare = weightedMeanSquare = frozenSquare;
            else if (reference is not null && config.Sigma2.ReadsResponse) throw new ArgumentException($"The reference states no mean square, which {config.Sigma2} reads.", nameof(working));
            sigma2 = config.Sigma2.For(meanSquare);
            rate = config.Sigma2.RateFor(meanSquare);
            if (held) log.Add(new("Sigma2", Text(config.Sigma2.Value!.Value), "held, so the structure step reads the marginal given σ²"));
            if (config.Sigma2.ReadsResponse)
                log.Add(new("ResponseMeanSquare", Text(meanSquare), reference is null
                    ? "s², the response's mean square about the centre level, which the scaled σ² prior reads"
                    : "s², the reference's, frozen apart from the working response"));
            if (rate is not null) log.Add(new("Sigma2Rate", rate.ToString(), "the hierarchical prior's rate, its mean a·s²; the state carries it"));
        }
        if (config.Noise is { SkewDrawn: true, SkewRescale: null } noise)
            log.Add(new("Noise.SkewRescale", noise.RescalesSkew ? "on" : "off", noise.SkewRescaleByDefault
                ? "the skewed horseshoe's default (2026-11-08): the sweep without the move along β·v fails Geweke's test at b = 1, the move passes"
                : "opt-in: on by default under the skewed horseshoe alone (2026-11-08)"));
        double coefficientConstant = reference?.CoefficientConstant ?? 1.0;
        if (coefficientConstant != 1.0) log.Add(new("CoefficientConstant", Text(coefficientConstant), "the reference's: the coefficient prior's covariance is σ²κ(τG)⁻¹"));

        // The caps, by the predictors of each kind.
        int continuous = 0, categorical = 0;
        foreach (BassPredictor p in data.Predictors)
        {
            if (p is ContinuousPredictor) continuous++;
            else categorical++;
            if (p is CategoricalPredictor factor)
                for (int c = 0; c < factor.Levels.Count; c++)
                    if (!factor.Observed[c])
                        log.Add(new($"Levels of '{factor.Name}'", $"'{factor.Levels[c]}' left out", "no training case takes it"));
        }
        int outputs = response is FunctionalResponse f ? f.OutputVariables.Count : 0;
        BassStructureSpec s = config.Structure;
        int maxInteraction = Cap("MaxInteraction", s.MaxInteraction, continuous, "continuous predictors", log);
        int maxCategorical = Cap("MaxCategoricalInteraction", s.MaxCategoricalInteraction, categorical, "categorical predictors", log);
        int maxFunctional = Cap("MaxFunctionalInteraction", s.MaxFunctionalInteraction, outputs, "output variables", log);
        if (s.Legendre is { } legendre)
        {
            // A polynomial-chaos fit (E19): Legendre terms over the continuous predictors alone, a scalar response, q_max the
            // continuous cap; mixed term types are not built.
            if (categorical > 0) throw new ArgumentException($"A Legendre basis reads continuous predictors alone, and the data hold {categorical} categorical ones.", nameof(data));
            if (functional) throw new ArgumentException("A Legendre basis serves a scalar response; a direct functional response's output part is a hinge product.", nameof(data));
            if (maxInteraction < 1) throw new ArgumentException("A Legendre basis needs a continuous predictor.", nameof(data));
            if (legendre.MaxDegree < maxInteraction)
                throw new ArgumentException($"Structure.Legendre.MaxDegree is {legendre.MaxDegree}, below q_max = {maxInteraction}, so the top orders would hold no term.", "config");
            maxCategorical = maxFunctional = 0;
            var space = new LegendreTermSpace(continuous, legendre.MaxDegree, maxInteraction);
            log.Add(new("Legendre", $"q_max {maxInteraction}, d_max {legendre.MaxDegree}", $"KHAOS's term space over {continuous} predictors, |𝒜| = {Text(space.Size)}"
                + (data.Predictors.OfType<ContinuousPredictor>().Any(c => c.Domain is null)
                    ? "; a predictor without a declared domain takes its training range as the box the basis is orthonormal on, which moves with the sample" : "")));
        }

        int bound = n - 1;   // the intercept counts toward the dimension bound
        int maxBasis;
        if (s.MaxBasis is { } declared)
        {
            maxBasis = Math.Min(declared, bound);
            if (declared > bound) log.Add(new("MaxBasis", Text(maxBasis), $"the declared {declared} exceeds N − 1, the bound with the intercept"));
        }
        else
        {
            maxBasis = Math.Min(1000, bound);
            log.Add(new("MaxBasis", Text(maxBasis), "min(1000, N − 1): R's 1,000, with the intercept counting toward the dimension bound"));
        }

        SupportRule rule = s.Support;
        // Weights that count rows count toward a scalar response's support, as the rows they stand for would.
        double caseCount = !functional && fitted is { IsUnweighted: false, Meaning: WeightMeaning.Replication } ? weights.Mass : cases;
        double support = Threshold(rule.Minimum, caseCount, "Structure.Support.Minimum", "SupportMinimum", "n_x", log);
        double? supportFunctional = functional
            ? Threshold(rule.MinimumFunctional, locations, "Structure.Support.MinimumFunctional", "SupportMinimumFunctional", "n_r", log)
            : null;

        IReadOnlyList<BassMonitoredPoint> points = response is MultiOutputResponse multiple
            ? ResolveOutputPoints(config.Sampling.MonitoredPoints, cases, multiple, log)
            : ResolvePoints(config.Sampling.MonitoredPoints, cases, functional ? locations : null, log);

        BassCoefficientPrior coefficients = config.Coefficients;
        if (coefficients.Form is RawRidgeForm && functional)
            throw new ArgumentException("KHAOS's raw ridge reads a scalar response's columns as evaluated; a direct functional response's are centred by construction.", "config");
        if (coefficients.Form is ModifiedGForm && s.Legendre is null)
            throw new ArgumentException("KHAOS's modified g-prior reads each term's order and degree, a Legendre term's; declare Structure.Legendre.", "config");
        if (coefficients.Form is ModifiedGForm && fitted.WeightArray is not null && !fitted.PriorStandsApart)
            throw new ArgumentException("KHAOS's modified g-prior reads the unweighted or a frozen prior Gram; the weighted one, which moves with W, is not built for it.", "config");
        // One observation's information under the prior's own columns: N, or Σwᵢ for a prior on rows the weights count.
        ResolvedTauPrior scale = coefficients.Scale.For(coefficients.UnitObservations(fitted.PriorStandsApart ? n : weights.Count));
        string unit = coefficients.Form switch { RidgeForm => "1, the normalized ridge's", RawRidgeForm => "1, the raw ridge's", _ => "1/N, the g-form's" };
        switch (coefficients.Scale.Kind)
        {
            case TauPriorKind.UnitInformationGamma:
                log.Add(new("TauPrior", $"Gamma({Text(scale.Gamma!.Value.Shape)}, rate {Text(scale.Gamma.Value.Rate)})", $"the unit-information Gamma, rate a/τ₀, mean τ₀ = {unit}"));
                break;
            case TauPriorKind.HyperG or TauPriorKind.HyperGOverN:
                log.Add(new("TauPrior", scale.ToString(), $"τ₀ = {unit}" + (scale.IsProper ? "" : "; improper, a ≤ 2")));
                break;
            case TauPriorKind.Gamma:
                log.Add(new("TauPrior", coefficients.Scale.ToString(), "on τ itself, whatever τ₀: g₀² = 1/τ ~ IG(a, b), KHAOS's"));
                break;
            case TauPriorKind.Held when coefficients.Scale.Tau is null:
                log.Add(new("HeldTau", Text(scale.HeldTau!.Value), $"unit information, τ₀ = {unit}"));
                break;
        }

        double? hotEnd = null;
        if (config.Sampling.Tempering is { } tempering)
        {
            if (tempering.DeclaredLevels() is null)
                log.Add(new("Tempering.Levels", $"{tempering.StartLevels()}, self-sizing within {tempering.MaxLevels}",
                            "Levels null: the ladder sizes itself during burn-in by its barrier, N ≈ 2Λ pairs"));
            hotEnd = tempering.DeclaredHotEnd();
            if (hotEnd is null)
            {
                double? halving = HalvingTau(coefficients.Form, scale, weights.Count, coefficientConstant);
                double picked = AutoHotEnd(tempering, config.Intercept, weights.Count, halving, coefficientConstant);
                if (!(picked > 0.0 && picked < 1.0))
                    throw new ArgumentException($"Sampling.Tempering.HotEnd auto-picks {Text(picked)} on N = {n} observations, outside (0, 1); set it.", "config");
                hotEnd = picked;
                log.Add(new("Tempering.HotEnd", Text(picked), AutoHotEndReason(tempering, config.Intercept, halving, picked)));
            }

            // A ladder's rungs where its tempered form exists (D05): the integrated form loses b = 0 under the improper
            // options, and under the flat intercept every b up to (1 − 2c)/N.
            if (tempering.Path == TemperingPath.Likelihood)
            {
                string improper = config.Intercept == GPriorCentre.Flat ? "the flat intercept" : "the Jeffreys σ² prior";
                string boundary = several is { } matrix
                    ? $"Σ's posterior needs ν₀ + Nb − 1 > K − 1, so b > {Text((matrix.Count - matrix.Prior.DegreesOfFreedom) / weights.Count)}"
                    : $"b > (1 − 2c)/N = {Text((1.0 - 2.0 * (sigma2?.Shape ?? 0.0)) / weights.Count)}";
                bool setHotEnd = tempering.DeclaredHotEnd() is not null;
                string reading = tempering.InverseTemperatures is not null ? "InverseTemperatures reaches" : setHotEnd ? "HotEnd is" : "HotEnd auto-picks";
                if (!BassTarget.Defines(tempering.Form, config.Intercept, sigma2, weights.Count, hotEnd.Value, held, several?.Prior))
                    throw new ArgumentException(
                        $"Sampling.Tempering.{reading} {Text(hotEnd.Value)}, where the integrated form is undefined under " +
                        $"{improper}, which needs {boundary}; " +
                        $"{(setHotEnd ? "raise it" : "set it above")}, or take the linear form (D05).", "config");
                if (tempering.ReferenceRung && !BassTarget.Defines(tempering.Form, config.Intercept, sigma2, weights.Count, 0.0, held, several?.Prior))
                    throw new ArgumentException(
                        $"Sampling.Tempering.ReferenceRung needs b = 0, where the integrated form is undefined under {improper}; the linear form reaches it (D05).", "config");
            }
        }

        if (config.Sampling.Smc is { } smc)
        {
            // SMC starts from the prior at b = 0, where the integrated form exists under the proper priors alone (D05).
            if (!BassTarget.Defines(smc.Form, config.Intercept, sigma2, weights.Count, 0.0, held, several?.Prior))
                throw new ArgumentException(
                    $"Sampling.Smc.Form is integrated, which has no b = 0 member under {(config.Intercept == GPriorCentre.Flat ? "the flat intercept" : "the Jeffreys σ² prior")}, " +
                    "where SMC starts from the prior; the linear form reaches it (D05).", "config");
            log.Add(new("Sampling.Chains", "not read", $"SMC runs Sampling.Smc.Replicates = {smc.Replicates} independent populations"));
            log.Add(new("Sampling.BurnIn, Samples, Thin, BatchSize, MinEssPerChain, StoppingScope", "not read",
                        $"SMC's draws are each replicate's final population of {smc.Particles} particles, and it ends at b = 1, with no R̂ stop"));
        }

        if (config.Sampling.FlatHistogram is { Run: { } flat })
        {
            log.Add(new("Sampling.Chains", "not read", $"the flat histogram runs Sampling.FlatHistogram.Run.Chains = {flat.Chains} chains"));
            log.Add(new("Sampling.BurnIn, Samples, BatchSize, MinEssPerChain, StoppingScope", "not read",
                        $"the flat histogram runs {flat.LearningIterations} learning and {flat.ProductionDraws} production iterations a chain, " +
                        $"each of Sampling.Thin = {config.Sampling.Thin} sweeps, with no R̂ stop"));
        }

        BassStepWindows? windows = null;
        if (config.Sampling.LocalSteps is { } steps)
        {
            // Each hinge variable's knots, as the prior counts them: its distinct scaled values.
            BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
            int[] Windows(int count, Func<int, int> knots)
            {
                var w = new int[count];
                for (int v = 0; v < count; v++) w[v] = steps.Window(knots(v));
                return w;
            }
            windows = new BassStepWindows(Windows(inputs.Cases.ContinuousCount, v => inputs.ContinuousKnots(v).Length),
                                          Windows(inputs.Outputs?.Count ?? 0, v => inputs.FunctionalKnots(v).Length));
            log.Add(new("LocalSteps windows", windows.ToString(), $"max(1, ⌈ω·n_v⌉) over each hinge variable's n_v knots, ω = {Text(steps.WindowFraction)}"));
        }
        BassFactorShares? shares = config.Sampling.FactorMoves is { } family ? BassFactorShares.Of(family) : null;

        return new BassResolved
        {
            Cases = cases, Locations = locations, Observations = n, MaxBasis = maxBasis,
            MaxInteraction = maxInteraction, MaxCategoricalInteraction = maxCategorical, MaxFunctionalInteraction = maxFunctional,
            SupportMinimum = support, SupportMinimumFunctional = supportFunctional,
            MonitoredPoints = points, ResponseMean = mean, CentreLevel = level, PriorLevel = priorLevel, CoefficientConstant = coefficientConstant,
            ResponseMeanSquare = meanSquare, WeightedMeanSquare = weightedMeanSquare,
            Weighting = weights, Sigma2 = sigma2, Sigma2Rate = rate, GivenSigma2 = config.Sigma2.Value, Scale = scale,
            StoppingGroups = BassDiagnostics.StoppingGroups(config.Sampling.StoppingScope),
            HotEnd = hotEnd, StepWindows = windows, FactorShares = shares, SplitThreshold = config.Sampling.Tempering?.Split?.Threshold,
            AutoPicks = log.AsReadOnly(), Outputs = several,
        };
    }

    // Several outputs on one basis (E34's slice 1): each output's mean and centre level, S_Y, Σ's prior carried from Sigma2's
    // ScaledToResponse(c) as decision 2's IW(K − 1 + 2c, 2c·S_Y), and the refusals of what later slices bring.
    private static BassResolvedOutputs ResolveOutputs(BassConfig config, MultiOutputResponse response, WorkingResponse fitted, List<BassAutoPick> log)
    {
        int k = response.Outputs, n = response.Cases;
        if (config.Noise is not null) throw new ArgumentException("Several outputs take Gaussian noise in E34's slice 1; the noise by row (decision 6) arrives in slice 4.", "config");
        if (config.HasRandomEffects) throw new ArgumentException("Several outputs' random effects arrive in E34's slice 4.", "config");
        // Σ's prior: the declared one, or Sigma2's ScaledToResponse(c) carried to K outputs.
        OutputCovariancePrior chosen = config.OutputCovariance
            ?? (config.Sigma2.Kind == Sigma2PriorKind.ScaledToResponse
                ? OutputCovariancePrior.ScaledToResponse(config.Sigma2.Shape!.Value)
                : throw new ArgumentException($"Several outputs carry Sigma2's ScaledToResponse(c) to Σ's prior (E34's decision 2), and Sigma2 is {config.Sigma2}; "
                                              + "declare OutputCovariance (decision 2's, standardized or a declared IW); Huang and Wand's prior arrives in slice 2.", "config"));

        ReadOnlySpan<double> means = fitted.OutputMeans;
        var levels = new double[k];
        for (int l = 0; l < k; l++) levels[l] = config.Intercept == GPriorCentre.Zero ? 0.0 : means[l];
        // S_Y: scalar BASS's s² carried to K outputs, each entry by its arithmetic, the products about the levels over n.
        ReadOnlySpan<double> values = response.Values;
        var squares = new double[k, k];
        var products = new double[n];
        for (int a = 0; a < k; a++)
            for (int d = 0; d <= a; d++)
            {
                for (int i = 0; i < n; i++) products[i] = (values[i * k + a] - levels[a]) * (values[i * k + d] - levels[d]);
                squares[a, d] = squares[d, a] = PairwiseSum.Of(products) / n;
            }
        // Refused: an output with no spread about its level, and under decision 2's form, whose Ψ₀ is S_Y's, outputs whose S_Y is singular or
        // nearly so, read on the correlation scale: a Cholesky pivot 1 − R² of an output on the ones before it under 1e-10.
        var unit = new double[k, k];
        for (int a = 0; a < k; a++)
        {
            if (!(squares[a, a] > 0.0))
                throw new ArgumentException(config.Intercept == GPriorCentre.Zero
                    ? $"Output '{response.Names[a]}' is zero, which leaves Σ no posterior (D03)."
                    : $"Output '{response.Names[a]}' is constant, which the intercept absorbs whole, leaving no spread to fit (D03).", nameof(response));
            for (int d = 0; d <= a; d++) unit[a, d] = unit[d, a] = squares[a, d] / (Math.Sqrt(squares[a, a]) * Math.Sqrt(squares[d, d]));
        }
        var root = new double[k, k];
        for (int a = 0; a < k && chosen.Kind == OutputCovariancePriorKind.ScaledToResponse; a++)
            for (int d = 0; d <= a; d++)
            {
                double sum = unit[a, d];
                for (int r = 0; r < d; r++) sum -= root[a, r] * root[d, r];
                if (a != d) { root[a, d] = sum / root[d, d]; continue; }
                if (!(sum >= 1e-10))
                    throw new ArgumentException($"Output '{response.Names[a]}' is collinear with the outputs before it about their levels (1 − R² = {Text(Math.Max(sum, 0.0))}), "
                                                + "so S_Y is singular or nearly so; collinear outputs need fewer outputs (E34's decision 2).", nameof(response));
                root[a, a] = Math.Sqrt(sum);
            }

        InverseWishart prior = chosen.For(squares);
        if (config.OutputCovariance is not null && config.Sigma2 != Sigma2Prior.ScaledToResponse(1e-3))
            log.Add(new("Sigma2", "not read", $"OutputCovariance declares Σ's prior, {chosen}, in Sigma2's place"));
        log.Add(new("Outputs", $"K = {k} ({string.Join(", ", response.Names.Select(name => $"'{name}'"))})",
                    "one basis and one Σ for every output, Σ integrated out (E34); the intercept and its convention per output"));
        log.Add(new("CentreLevels", string.Join(", ", levels.Select(Text)), config.Intercept switch
        {
            GPriorCentre.Zero => "the zero centre",
            GPriorCentre.Flat => "the outputs' means, the flat intercept's estimates",
            _ => "the outputs' means, level centring",
        }));
        log.Add(new("SigmaPrior", $"{chosen}: IW({Text(prior.DegreesOfFreedom)}, Ψ₀)", chosen.Kind switch
        {
            OutputCovariancePriorKind.ScaledToResponse => (config.OutputCovariance is null ? "Sigma2's ScaledToResponse(c) carried to K outputs, " : "")
                + "E34's decision 2: ν₀ = K − 1 + 2c, Ψ₀ = 2c·S_Y, S_Y the outputs' mean squares and cross-products about their levels; at K = 1 IG(c, c·s²)",
            OutputCovariancePriorKind.Standardized => "ν₀ = K − 1 + 2c, Ψ₀ = 2c·diag S_Y: invariant to rescaling the outputs, not to mixing them",
            OutputCovariancePriorKind.HuangWand => "conditional IW with one gamma rate per output; the marginal standard deviations have half-t laws",
            _ => "declared in the outputs' own units",
        }));
        return new BassResolvedOutputs(response.Names, [.. means], levels, squares, prior, chosen.ResolveRates(squares));
    }

    /// <summary>
    /// The hot end a ladder auto-picks on <paramref name="observations"/> observations (decided 2026-10-28): 1/N, the b at
    /// which τ₀ makes a unit-information column's shrinkage, τ/(τ + b) under the g-form and τ/(τ + bN) under the ridge,
    /// one half, so the tempered likelihood weighs one observation, as the unit-information prior does; 2/N on the
    /// likelihood path's integrated form under the flat intercept, whose pole at (1 − 2c)/N it clears with ν = Nb − 1 = 1;
    /// and under a held τ the smaller of these and the b at which the held τ halves it (<paramref name="halvingTau"/>). A coefficient
    /// prior's constant κ reads τ as τ/κ, so the unit-information b is 1/(κN), the pole's staying where it is.
    /// </summary>
    internal static double AutoHotEnd(TemperingSpec tempering, GPriorCentre intercept, double observations, double? halvingTau, double constant = 1.0)
    {
        double unit = Pole(tempering, intercept) ? 2.0 / observations : 1.0 / (constant * observations);
        return halvingTau is double tau ? Math.Min(tau, unit) : unit;
    }

    /// <summary>
    /// Under a held τ, the b at which it halves a unit-information column's shrinkage: τ under the g-form, τ/N under either
    /// ridge, whose columns carry about N observations' information each, each over the coefficient prior's constant; null when τ is sampled.
    /// </summary>
    internal static double? HalvingTau(CoefficientForm form, ResolvedTauPrior scale, double observations, double constant = 1.0)
        => scale.HeldTau is double tau ? (form is RidgeForm or RawRidgeForm ? tau / observations : tau) / constant : null;

    private static bool Pole(TemperingSpec tempering, GPriorCentre intercept)
        => tempering.Path == TemperingPath.Likelihood && tempering.Form == TemperedForm.Integrated && intercept == GPriorCentre.Flat;

    private static string AutoHotEndReason(TemperingSpec tempering, GPriorCentre intercept, double? halvingTau, double picked)
    {
        string unit = Pole(tempering, intercept)
            ? "2/N, past the integrated form's pole under the flat intercept at (1 − 2c)/N, with ν = Nb − 1 = 1"
            : "1/N, where τ₀ halves a unit-information column's shrinkage: the tempered likelihood weighs one observation, as the unit-information prior does";
        if (halvingTau is not double tau) return unit;
        return picked == tau ? $"the b at which the held τ halves that shrinkage, at most {unit}" : $"{unit}, below the b at which the held τ halves that shrinkage";
    }

    // The declared points, checked against the data, or the defaults; a direct response's missing location is the middle one.
    private static IReadOnlyList<BassMonitoredPoint> ResolvePoints(IReadOnlyList<BassMonitoredPoint>? declared, int cases, int? locations,
                                                                     List<BassAutoPick> log)
    {
        if (declared is null)
        {
            IReadOnlyList<BassMonitoredPoint> defaults = BassDiagnostics.DefaultPoints(cases, locations);
            log.Add(new("MonitoredPoints", string.Join("; ", defaults), "training cases at the quartile positions of the data's order" +
                        (locations is null ? "" : ", at the middle output location")));
            return defaults;
        }
        var points = new BassMonitoredPoint[declared.Count];
        for (int k = 0; k < points.Length; k++)
        {
            BassMonitoredPoint p = declared[k];
            if (p.Case >= cases) throw new ArgumentException($"Sampling.MonitoredPoints holds case {p.Case}; the data have {cases}.", "config");
            if (locations is not int m)
            {
                if (p.Location is not null) throw new ArgumentException("Sampling.MonitoredPoints gives a location, but the response is scalar.", "config");
                points[k] = p;
                continue;
            }
            if (p.Location >= m) throw new ArgumentException($"Sampling.MonitoredPoints holds location {p.Location}; the response has {m}.", "config");
            points[k] = p with { Location = p.Location ?? (m - 1) / 2 };
        }
        if (points.Distinct().Count() != points.Length)
            throw new ArgumentException("Sampling.MonitoredPoints names one point twice once the middle location fills a missing one.", "config");
        return points;
    }

    // Several outputs' monitored points (E34): each case, the defaults' or a declared one's, at every output, unless a declared point names
    // its output in its location.
    private static IReadOnlyList<BassMonitoredPoint> ResolveOutputPoints(IReadOnlyList<BassMonitoredPoint>? declared, int cases, MultiOutputResponse response,
                                                                         List<BassAutoPick> log)
    {
        int k = response.Outputs;
        IReadOnlyList<BassMonitoredPoint> chosen = declared ?? BassDiagnostics.DefaultPoints(cases, null);
        var points = new List<BassMonitoredPoint>();
        foreach (BassMonitoredPoint p in chosen)
        {
            if (p.Case >= cases) throw new ArgumentException($"Sampling.MonitoredPoints holds case {p.Case}; the data have {cases}.", "config");
            if (p.Location is int output)
            {
                if (output >= k) throw new ArgumentException($"Sampling.MonitoredPoints holds output {output}; the response has {k}.", "config");
                points.Add(p);
            }
            else for (int l = 0; l < k; l++) points.Add(p with { Location = l });
        }
        if (points.Distinct().Count() != points.Count)
            throw new ArgumentException("Sampling.MonitoredPoints names one point twice once a case without an output fills in every output.", "config");
        if (declared is null)
            log.Add(new("MonitoredPoints", string.Join("; ", points), "training cases at the quartile positions of the data's order, at every output"));
        return points;
    }

    // A cap at most the count of its variables, as R caps maxInt by pdes, pcat and pfunc.
    private static int Cap(string setting, int declared, int available, string what, List<BassAutoPick> log)
    {
        if (declared <= available) return declared;
        log.Add(new(setting, Text(available), $"the declared {declared} exceeds the {available} {what}"));
        return available;
    }

    // R's default min(20, 0.1·n) for a part of size n, or the declared threshold, which must be one some term can meet.
    private static double Threshold(double? declared, double size, string field, string setting, string sizeName, List<BassAutoPick> log)
    {
        if (declared is { } d)
        {
            if (d > size) throw new ArgumentException($"{field} is {Text(d)}, above the {size} a term's part can reach, so no term could meet it.", "config");
            return d;
        }
        double resolved = Math.Min(20.0, 0.1 * size);
        log.Add(new(setting, Text(resolved), $"min(20, 0.1·{sizeName}), R's default"));
        return resolved;
    }

    private static string Text(double value) => value.ToString("R", CultureInfo.InvariantCulture);
    private static string Text(int value) => value.ToString(CultureInfo.InvariantCulture);
}

/// <summary>The local steps' window w per hinge variable (<see cref="LocalSteps.Window"/>): the continuous predictors', then
/// the output variables'.</summary>
public sealed record BassStepWindows(int[] Continuous, int[] Functional)
{
    /// <inheritdoc/>
    public override string ToString()
        => $"continuous [{string.Join(", ", Continuous)}]" + (Functional.Length > 0 ? $"; functional [{string.Join(", ", Functional)}]" : "");
}

/// <summary>
/// The fixed-M factor moves' shares of the change's slot (<see cref="FactorMoves"/>), each its weight over their sum:
/// a factor's redraw or local step, a prune, a grow (half the prune-or-grow weight each) and a variable change.
/// </summary>
public sealed record BassFactorShares(double Change, double Prune, double Grow, double VariableChange)
{
    /// <summary>The shares <paramref name="family"/>'s weights give.</summary>
    public static BassFactorShares Of(FactorMoves family)
    {
        ArgumentNullException.ThrowIfNull(family);
        double total = family.Change + family.PruneOrGrow + family.VariableChange;
        return new(family.Change / total, 0.5 * family.PruneOrGrow / total, 0.5 * family.PruneOrGrow / total, family.VariableChange / total);
    }
}
