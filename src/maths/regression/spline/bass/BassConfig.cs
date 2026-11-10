using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A BASS fit's configuration (the contract shapes, decided 2026-10-28). Every field is a value, so the record
/// serializes to JSON whole and is its own manifest, and the run's fingerprint hashes that serialization with the data;
/// <c>BarsConfig</c>'s component instances fingerprint their types only (OQ10). The defaults are D01–D06's. The
/// auto-picks (the basis cap, the support thresholds, N, s² and the centre's level) are resolved at entry and reported
/// in <see cref="BassResolved"/>. A setting joins with the code that reads it.
/// </summary>
public sealed record BassConfig
{
    /// <summary>The term space and its support (D01, D06).</summary>
    public BassStructureSpec Structure { get; init; } = new();

    /// <summary>The count prior, λ integrated out (D02).</summary>
    public BassCountPrior Count { get; init; } = new();

    /// <summary>The coefficient prior, a covariance form and τ's prior: the g-form with Zellner–Siow's unit-information
    /// Gamma, by default (D03).</summary>
    public BassCoefficientPrior Coefficients { get; init; } = BassCoefficientPrior.G();

    /// <summary>σ²'s prior: IG(c, c·s²) scaled to the response, c = 10⁻³, by default (D03).</summary>
    public Sigma2Prior Sigma2 { get; init; } = Sigma2Prior.ScaledToResponse(1e-3);

    /// <summary>
    /// Σ's prior for several outputs (E34; <see cref="OutputCovariancePrior"/>): null, the default, carries <see cref="Sigma2"/>'s
    /// ScaledToResponse(c) to K outputs as decision 2's IW(K − 1 + 2c, 2c·S_Y), and refuses Sigma2's other kinds. Declared, it is read in
    /// Sigma2's place. A single response refuses one. Left out of the serialized form when null, so every config without it keeps its
    /// form and fingerprint.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public OutputCovariancePrior? OutputCovariance { get; init; }

    /// <summary>
    /// The intercept convention: level centring (<see cref="GPriorCentre.Mean"/>) by default; Liang et al.'s flat
    /// intercept (<see cref="GPriorCentre.Flat"/>) and DiMatteo, Genovese and Kass's zero centre selectable (D03). The
    /// median centre is for the noise extensions' weights and is refused.
    /// </summary>
    public GPriorCentre Intercept { get; init; } = GPriorCentre.Mean;

    /// <summary>
    /// The noise beyond the Gaussian (E1; <see cref="BassNoise"/>): null, the default, is Gaussian BASS. Left out of the serialized form
    /// when null, so a Gaussian config's form and fingerprint stand as they were. A fit under a declared noise runs chains, tempered ladders
    /// or SMC, each on the likelihood path in the integrated form, whose rungs temper the likelihood given the latents; not the flat
    /// histogram, whose π_w reads the marginal over a and σ² that the latents couple.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public BassNoise? Noise { get; init; }

    /// <summary>
    /// Named random-effect intercept blocks (E21; <see cref="BassRandomIntercept"/>): null, the default, holds none. Left out of
    /// the serialized form when null, so a config without them keeps its form and fingerprint. Refused beside a noise, whose latent weights
    /// move each sweep while the block is integrated under unit weights. The ordered collection is copied and compared by value.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public IReadOnlyList<BassRandomIntercept>? RandomIntercepts
    {
        get => _randomIntercepts;
        init => _randomIntercepts = value is null ? null : new EffectCollection<BassRandomIntercept>(value);
    }

    private readonly EffectCollection<BassRandomIntercept>? _randomIntercepts;

    /// <summary>Independent named slopes, each with an explicitly chosen input scale; null declares none.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public IReadOnlyList<BassRandomSlope>? RandomSlopes
    {
        get => _randomSlopes;
        init => _randomSlopes = value is null ? null : new EffectCollection<BassRandomSlope>(value);
    }

    private readonly EffectCollection<BassRandomSlope>? _randomSlopes;
    internal bool HasRandomEffects => (RandomIntercepts?.Count ?? 0) + (RandomSlopes?.Count ?? 0) > 0;
    internal IBassRandomEffect[] EffectSpecifications => [.. RandomIntercepts ?? Array.Empty<BassRandomIntercept>(), .. RandomSlopes ?? Array.Empty<BassRandomSlope>()];

    private sealed class EffectCollection<T>(IReadOnlyList<T> source)
        : IReadOnlyList<T>, IEquatable<EffectCollection<T>>
    {
        private readonly T[] _items = source.ToArray();
        public int Count => _items.Length;
        public T this[int index] => _items[index];
        public IEnumerator<T> GetEnumerator() => ((IEnumerable<T>)_items).GetEnumerator();
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
        public bool Equals(EffectCollection<T>? other) => ReferenceEquals(this, other) || other is not null && _items.SequenceEqual(other._items);
        public override bool Equals(object? other) => other is EffectCollection<T> collection && Equals(collection);
        public override int GetHashCode()
        {
            var hash = new HashCode();
            hash.Add(Count);
            foreach (T item in _items) hash.Add(item);
            return hash.ToHashCode();
        }
    }

    /// <summary>The run: chains or SMC's replicates, seeds, lengths, the R̂ stop and the move schedule.</summary>
    public BassSampling Sampling { get; init; } = new();

    /// <summary>
    /// The canonical JSON options: enums by name, compact. The fingerprint hashes this serialization. Reading, a field
    /// the config does not have is refused, so a misspelled setting cannot fall back to its default unnoticed; a field
    /// left out takes its default; a variant's <c>Kind</c> may stand anywhere among its fields.
    /// </summary>
    public static JsonSerializerOptions JsonOptions { get; } = CreateJsonOptions();

    private static JsonSerializerOptions CreateJsonOptions()
    {
        var options = new JsonSerializerOptions
        {
            Converters = { new JsonStringEnumConverter() },
            UnmappedMemberHandling = JsonUnmappedMemberHandling.Disallow,
            AllowOutOfOrderMetadataProperties = true,
        };
        options.MakeReadOnly(populateMissingResolver: true);
        return options;
    }

    /// <summary>This config as canonical JSON.</summary>
    public string ToJson() => JsonSerializer.Serialize(this, JsonOptions);

    /// <summary>A config read from <paramref name="json"/>, checked as <see cref="Validate"/> checks a built one.</summary>
    public static BassConfig FromJson(string json)
    {
        BassConfig config = JsonSerializer.Deserialize<BassConfig>(json, JsonOptions)
                            ?? throw new ArgumentException("The JSON holds no config.", nameof(json));
        config.Validate();
        return config;
    }

    /// <summary>
    /// Checks every field against its range, naming the field in the reason. The data-dependent checks (a support
    /// threshold no term could meet, a response the intercept convention refuses) are <see cref="BassResolved"/>'s.
    /// </summary>
    public void Validate()
    {
        BassStructureSpec s = Structure ?? throw Bad("Structure", "is missing");
        if (s.Degree < 1) throw Bad("Structure.Degree", "must be a positive integer");
        if (s.MaxInteraction < 1) throw Bad("Structure.MaxInteraction", "must be at least 1");
        if (s.MaxCategoricalInteraction < 1) throw Bad("Structure.MaxCategoricalInteraction", "must be at least 1");
        if (s.MaxFunctionalInteraction < 1) throw Bad("Structure.MaxFunctionalInteraction", "must be at least 1");
        if (s.MaxBasis is < 1) throw Bad("Structure.MaxBasis", "must be at least 1");
        SupportRule support = s.Support ?? throw Bad("Structure.Support", "is missing");
        if (support.Kind is not (SupportKind.Count or SupportKind.ColumnMass)) throw Bad("Structure.Support.Kind", "is not a known rule");
        if (support.Minimum is { } min && !(min > 0.0 && double.IsFinite(min))) throw Bad("Structure.Support.Minimum", "must be positive and finite");
        if (support.MinimumFunctional is { } minF && !(minF > 0.0 && double.IsFinite(minF)))
            throw Bad("Structure.Support.MinimumFunctional", "must be positive and finite");
        if (!(s.AngleFloor > 0.0 && s.AngleFloor < 1.0)) throw Bad("Structure.AngleFloor", "must lie in (0, 1)");
        if (s.MaxVarianceInflation is { } vif && !(vif > 1.0 && double.IsFinite(vif)))
            throw Bad("Structure.MaxVarianceInflation", "must be finite and above 1");
        if (s.CategoricalSubsets is not (CategoricalSubsetMeasure.SizeUniform or CategoricalSubsetMeasure.Uniform))
            throw Bad("Structure.CategoricalSubsets", "is not a known measure");
        if (!(s.OrderWeight > 0.0 && double.IsFinite(s.OrderWeight)))
            throw Bad("Structure.OrderWeight", "must be positive and finite, or a death could have no reverse birth");
        if (!(s.VariableWeight > 0.0 && double.IsFinite(s.VariableWeight)))
            throw Bad("Structure.VariableWeight", "must be positive and finite, or a death could have no reverse birth");
        if (s.Legendre is { } legendre)
        {
            if (legendre.Problem() is { } problem) throw Bad($"Structure.Legendre.{problem.Setting}", problem.Reason);
            if (s.Degree != 1) throw Bad("Structure.Degree", "is the hinges' power, which a Legendre basis does not read; its degrees are capped by Structure.Legendre.MaxDegree");
            if (s.Support.Kind == SupportKind.ColumnMass)
                throw Bad("Structure.Support.Kind", "is the column mass, which reads hinges normalized to a maximum of one; a Legendre column sums to about zero, so a Legendre basis takes the count");
        }

        BassCountPrior count = Count ?? throw Bad("Count", "is missing");
        if (!(count.Shape > 0.0 && double.IsFinite(count.Shape))) throw Bad("Count.Shape", "must be positive and finite");
        if (!(count.Rate > 0.0 && double.IsFinite(count.Rate))) throw Bad("Count.Rate", "must be positive and finite");

        BassCoefficientPrior coefficients = Coefficients ?? throw Bad("Coefficients", "is missing");
        if (coefficients.Form is not (GForm or RidgeForm or RawRidgeForm or ModifiedGForm)) throw Bad("Coefficients.Form", coefficients.Form is null ? "is missing" : "is not a known form");
        TauPrior scale = coefficients.Scale ?? throw Bad("Coefficients.Scale", "is missing");
        if (coefficients.Form is RidgeForm or RawRidgeForm && scale is { IsHeld: true, Tau: null })
            throw Bad("Coefficients.Scale.Tau", "must be declared under the ridge, whose held scale is part of the model (D03)");

        if (Sigma2 is null) throw Bad("Sigma2", "is missing");
        if (Intercept == GPriorCentre.Median) throw Bad("Intercept", "the median centre is for the noise extensions' weights (D03)");
        if (Intercept is not (GPriorCentre.Mean or GPriorCentre.Flat or GPriorCentre.Zero)) throw Bad("Intercept", "is not a known convention");
        if (!scale.IsProper && Intercept == GPriorCentre.Flat)
            throw Bad("Coefficients.Scale", $"is {scale}, improper, which leaves τ no posterior under the flat intercept, whose intercept-only structure has no column inside the coefficient prior; take a > 2");

        if (Noise is { } noise)
        {
            if (noise.Interweave && Sigma2.IsHeld)
                throw Bad("Noise.Interweave", "moves σ with v under the z parameterization, which a held σ² does not move");
            // The noise's blocks temper the likelihood given its latents, so a, σ², β and each vᵢ keep their conditionals at every b: the
            // likelihood path in the integrated form. The linear form and the powered path temper the marginal over a and σ², which couples
            // every latent, and the flat histogram's π_w and WHAM read that marginal's rung laws, prior·e^{b·ℓ₁}.
            if (Sampling?.Tempering is { } ladder)
            {
                if (ladder.Path != TemperingPath.Likelihood)
                    throw Bad("Sampling.Tempering.Path", "must be the likelihood path under a noise: the powered path powers the marginal over a and σ², which couples every latent");
                if (ladder.Form != TemperedForm.Integrated)
                    throw Bad("Sampling.Tempering.Form", "must be the integrated form under a noise: the linear form tempers the marginal over a and σ², b·ℓ₁, which couples every latent");
            }
            if (Sampling?.Smc is { Form: not TemperedForm.Integrated })
                throw Bad("Sampling.Smc.Form", "must be the integrated form under a noise: the linear form tempers the marginal over a and σ², b·ℓ₁, which couples every latent");
            if (Sampling?.FlatHistogram is not null)
                throw Bad("Noise", "takes no Sampling.FlatHistogram: its walkers' π_w and WHAM read the rung laws prior·e^{b·ℓ₁} of the marginal over a and σ², "
                               + "whose weights couple every latent, so no vᵢ keeps a step of its own; take Sampling.Tempering or Sampling.Smc");
            bool drawsPrior = Sampling?.Smc is not null || Sampling?.Tempering is { ReferenceRung: true };
            if (drawsPrior && noise.Skewed && !Sigma2.IsProper)
                throw Bad("Sigma2", "is Jeffreys, improper, while the noise's skew puts σ² on the record and b = 0 draws it from its prior (an SMC start or a reference rung); take an inverse gamma, a held σ² or the hierarchical prior");
            if (noise.Reference == BassNoiseReference.None && (Intercept == GPriorCentre.Mean || Sigma2.ReadsResponse))
                throw Bad("Noise.Reference", "is None, while the intercept's level or σ²'s prior reads the response; state Moments or Robust, or take priors that read nothing");
        }

        CheckRandomIntercepts(RandomIntercepts);
        CheckRandomSlopes(RandomSlopes);
        if (RandomIntercepts is { Count: > 0 } && Noise is not null)
            throw Bad("RandomIntercepts", "take no Noise: blocks are integrated by whitening under unit weights, and a noise's latent weights move each sweep; not built");
        if (RandomSlopes is { Count: > 0 } && Noise is not null)
            throw Bad("RandomSlopes", "take no Noise: moving row weights arrive in slice 4");
        for (int k = 0; k < (RandomIntercepts?.Count ?? 0); k++)
        {
            if (RandomIntercepts![k].ScalePrior is not null && Sampling?.FlatHistogram is not null)
                throw Bad($"RandomIntercepts[{k}].ScalePrior", "takes no Sampling.FlatHistogram: its weights and WHAM read the rung laws of a marginal that the drawn s moves, which is not gated; take Sampling.Tempering or Sampling.Smc, or hold s");
        }
        if (RandomSlopes?.Any(s => s.ScalePrior is not null) == true && Sampling?.FlatHistogram is not null)
            throw Bad("RandomSlopes.ScalePrior", "takes no Sampling.FlatHistogram: its weights read a marginal that the drawn s moves; take Sampling.Tempering or Sampling.Smc, or hold s");

        BassSampling run = Sampling ?? throw Bad("Sampling", "is missing");
        if (run.Chains < 1) throw Bad("Sampling.Chains", "must be at least 1");
        if (run.BurnIn < 0) throw Bad("Sampling.BurnIn", "must not be negative");
        if (run.Samples < 1) throw Bad("Sampling.Samples", "must be at least 1");
        if (run.Thin < 1) throw Bad("Sampling.Thin", "must be at least 1");
        if (!(run.RHatTarget == 0.0 || (run.RHatTarget > 1.0 && double.IsFinite(run.RHatTarget))))
            throw Bad("Sampling.RHatTarget", "must be 0, no early stop, or a finite target above 1");
        if (run.RHatTarget > 0.0 && run.Chains < 2) throw Bad("Sampling.Chains", "must be at least 2 for the R̂ stop (Sampling.RHatTarget)");
        if (run.BatchSize < 0) throw Bad("Sampling.BatchSize", "must not be negative");
        if (!(run.MinEssPerChain >= 0.0 && double.IsFinite(run.MinEssPerChain))) throw Bad("Sampling.MinEssPerChain", "must be finite and not negative");
        if (run.StoppingScope is not (BassStoppingScope.All or BassStoppingScope.Predictive)) throw Bad("Sampling.StoppingScope", "is not a known scope");
        if (run.MonitoredPoints is { } points)
        {
            if (points.Count == 0) throw Bad("Sampling.MonitoredPoints", "must hold a point, or be null for the defaults");
            if (points.Any(p => p.Case < 0 || p.Location < 0)) throw Bad("Sampling.MonitoredPoints", "must hold indices, not negative");
            if (points.Distinct().Count() != points.Count) throw Bad("Sampling.MonitoredPoints", "must not repeat a point");
        }
        switch (run.Moves)
        {
            case PriorCancelling p when !(p.C > 0.0 && p.C <= 0.5):
                throw Bad("Sampling.Moves.C", "must lie in (0, ½], so a birth and a death leave the change move its share");
            case EqualOverAvailable or PriorCancelling:
                break;
            case null:
                throw Bad("Sampling.Moves", "is missing");
            default:
                throw Bad("Sampling.Moves", "is not a known schedule");
        }
        if (run.LocalSteps is { } steps)
        {
            if (!(steps.Probability > 0.0 && steps.Probability <= 1.0)) throw Bad("Sampling.LocalSteps.Probability", "must lie in (0, 1]; null turns the steps off");
            if (!(steps.WindowFraction > 0.0 && steps.WindowFraction <= 1.0)) throw Bad("Sampling.LocalSteps.WindowFraction", "must lie in (0, 1]");
        }
        if (run.FactorMoves is { } family)
        {
            if (!(family.Change >= 0.0 && double.IsFinite(family.Change))) throw Bad("Sampling.FactorMoves.Change", "must be finite and not negative");
            if (!(family.PruneOrGrow >= 0.0 && double.IsFinite(family.PruneOrGrow))) throw Bad("Sampling.FactorMoves.PruneOrGrow", "must be finite and not negative");
            if (!(family.VariableChange >= 0.0 && double.IsFinite(family.VariableChange)))
                throw Bad("Sampling.FactorMoves.VariableChange", "must be finite and not negative");
            if (!(family.PruneOrGrow + family.VariableChange > 0.0))
                throw Bad("Sampling.FactorMoves", "must weigh a prune or grow or a variable change; null turns the family off");
            if (run.LocalSteps is not null && !(family.Change > 0.0))
                throw Bad("Sampling.FactorMoves.Change", "must be positive with local steps, which run inside the change");
        }
        if (run.Births is { } births)
        {
            if (births.History is null && births.Residual is null) throw Bad("Sampling.Births", "must declare a component; null turns the mixture off");
            double shares = 0.0;
            if (births.History is { } history)
            {
                if (!(history.Share > 0.0 && history.Share < 1.0))
                    throw Bad("Sampling.Births.History.Share", "must lie in (0, 1), so the fit's own birth keeps a positive share");
                if (!(history.Floor > 0.0 && history.Floor <= 1.0))
                    throw Bad("Sampling.Births.History.Floor", "must lie in (0, 1], so every variable keeps a positive weight");
                if (history.Refresh < 1) throw Bad("Sampling.Births.History.Refresh", "must be at least 1 sweep");
                shares += history.Share;
            }
            if (births.Residual is { } residual)
            {
                if (!(residual.Share > 0.0 && residual.Share < 1.0))
                    throw Bad("Sampling.Births.Residual.Share", "must lie in (0, 1), so the fit's own birth keeps a positive share");
                shares += residual.Share;
            }
            if (!(shares < 1.0)) throw Bad("Sampling.Births", "the components' shares must sum below 1, so the fit's own birth keeps a positive share");
        }
        if (s.Legendre is not null)
        {
            if (run.LocalSteps is not null) throw Bad("Sampling.LocalSteps", "steps D01's hinge knots, which a Legendre basis has none of; its changes are KHAOS's mutations");
            if (run.FactorMoves is not null) throw Bad("Sampling.FactorMoves", "is D01's fixed-M family over hinge parts; a Legendre basis changes its terms by KHAOS's mutations");
        }
        if (run.Smc is { } smc)
        {
            if (!scale.IsProper) throw Bad("Coefficients.Scale", $"is {scale}, improper, which SMC cannot start from (D05); take a > 2");
            if (run.Tempering is not null) throw Bad("Sampling.Tempering", "must be null with Sampling.Smc, which tempers its own sequence from the prior (D05)");
            if (run.RHatTarget > 0.0) throw Bad("Sampling.RHatTarget", "must be 0 with Sampling.Smc: an SMC run ends when b reaches 1, with no R̂ stop");
            if (smc.Particles < 2) throw Bad("Sampling.Smc.Particles", "must be at least 2, so the population has weights to read and a normalizer error");
            if (smc.Replicates < 1) throw Bad("Sampling.Smc.Replicates", "must be at least 1");
            if (!(smc.TargetConditionalEss > 0.0 && smc.TargetConditionalEss < 1.0)) throw Bad("Sampling.Smc.TargetConditionalEss", "must lie in (0, 1)");
            if (smc.MutationSteps < 0) throw Bad("Sampling.Smc.MutationSteps", "must not be negative; 0 makes the run sequential importance sampling");
            if (smc.AdaptiveMutation is { } adaptive)
            {
                if (smc.MutationSteps < 1) throw Bad("Sampling.Smc.MutationSteps", "must be at least 1 under Sampling.Smc.AdaptiveMutation, which spends it per stage");
                if (!(adaptive.StayProbability > 0.0 && adaptive.StayProbability < 1.0)) throw Bad("Sampling.Smc.AdaptiveMutation.StayProbability", "must lie in (0, 1)");
            }
            if (smc.StepCap is { } cap && !(cap > 1.0 && double.IsFinite(cap)))
                throw Bad("Sampling.Smc.StepCap", "must exceed 1 and be finite, or the schedule never reaches b = 1; null leaves the conditional-ESS rule alone");
            if (!(smc.ResampleBelow >= 0.0 && smc.ResampleBelow <= 1.0)) throw Bad("Sampling.Smc.ResampleBelow", "must lie in [0, 1]");
            if (smc.Resampling is not (ResamplingScheme.Multinomial or ResamplingScheme.Stratified or ResamplingScheme.Systematic or ResamplingScheme.Residual))
                throw Bad("Sampling.Smc.Resampling", "is not a known scheme");
            if (smc.Form is not (TemperedForm.Integrated or TemperedForm.Linear)) throw Bad("Sampling.Smc.Form", "is not a known form");
            if (smc.MaxStages < 1) throw Bad("Sampling.Smc.MaxStages", "must be at least 1");
        }
        if (run.Tempering is { } t) CheckTempering(t, "Sampling.Tempering", scale);
        if (run.FlatHistogram is { } flat)
        {
            if (run.Smc is not null) throw Bad("Sampling.FlatHistogram", "must be null with Sampling.Smc: a fit runs one sampler");
            if (run.Tempering is not null)
                throw Bad("Sampling.Tempering", "must be null with Sampling.FlatHistogram, whose pilot declares its own ladder (Sampling.FlatHistogram.Pilot.Tempering)");
            if (run.RHatTarget > 0.0) throw Bad("Sampling.RHatTarget", "must be 0 with Sampling.FlatHistogram, which runs to its lengths with no R̂ stop");
            FlatHistogramSpec engine = flat.Run ?? throw Bad("Sampling.FlatHistogram.Run", "is missing");
            try
            {
                engine.Validate();
            }
            catch (ArgumentException e)
            {
                string reason = e.Message, suffix = $" (Parameter '{e.ParamName}')";
                if (reason.EndsWith(suffix, StringComparison.Ordinal)) reason = reason[..^suffix.Length];
                throw new ArgumentException($"Sampling.FlatHistogram.Run.{e.ParamName}: {reason}", "config", e);
            }
            if (flat.MaxIterations is < 1) throw Bad("Sampling.FlatHistogram.MaxIterations", "must be at least 1, or null to run to the lengths");
            if (engine.NeedsPilot && flat.Pilot is null)
                throw Bad("Sampling.FlatHistogram.Pilot", "is missing, and the run starts from WHAM or reads its window, b_lo or K from a pilot");
            if (!engine.NeedsPilot && flat.Pilot is not null)
                throw Bad("Sampling.FlatHistogram.Pilot", "must be null: the run declares its window, b_lo and K and starts flat, so nothing reads a pilot");
            bool lowestAtZero = engine.LowSlope is double low ? low == 0.0 : flat.Pilot?.Tempering?.ReferenceRung == true;
            if (lowestAtZero && !scale.IsProper)
                throw Bad("Coefficients.Scale", $"is {scale}, improper, which leaves the prior over (structure, τ) improper, and π_w below a window from b_lo = 0 is that prior; take a > 2, or a window from b_lo > 0");
            if (flat.Pilot is { } pilot)
            {
                if (pilot.Chains < 1) throw Bad("Sampling.FlatHistogram.Pilot.Chains", "must be at least 1");
                if (pilot.BurnIn < 0) throw Bad("Sampling.FlatHistogram.Pilot.BurnIn", "must not be negative");
                if (pilot.TracedScans < 2) throw Bad("Sampling.FlatHistogram.Pilot.TracedScans", "must be at least 2, so every rung has draws to spread and weigh");
                TemperingSpec ladder = pilot.Tempering ?? throw Bad("Sampling.FlatHistogram.Pilot.Tempering", "is missing");
                CheckTempering(ladder, "Sampling.FlatHistogram.Pilot.Tempering", scale);
                if (ladder.Path != TemperingPath.Likelihood)
                    throw Bad("Sampling.FlatHistogram.Pilot.Tempering.Path", "must be the likelihood path: WHAM reads rung laws prior·e^{b·ℓ₁}");
                if (ladder.Form != TemperedForm.Linear)
                    throw Bad("Sampling.FlatHistogram.Pilot.Tempering.Form", "must be the linear form: WHAM reads rung laws prior·e^{b·ℓ₁}, which the integrated form's are not (D05)");
            }
        }
    }

    // A ladder's declaration, its fields named under `field`.
    internal static void CheckRandomIntercepts(IReadOnlyList<BassRandomIntercept>? blocks)
    {
        if (blocks is null) return;
        var names = new HashSet<string>(StringComparer.Ordinal);
        for (int k = 0; k < blocks.Count; k++)
        {
            string field = $"RandomIntercepts[{k}]";
            BassRandomIntercept block = blocks[k] ?? throw Bad(field, "is null");
            if (string.IsNullOrWhiteSpace(block.Grouping)) throw Bad($"{field}.Grouping", "must name a grouping");
            if (!names.Add(block.Grouping)) throw Bad($"{field}.Grouping", $"duplicates grouping '{block.Grouping}'");
            if (!(block.Scale > 0) || !double.IsFinite(block.Scale)) throw Bad($"{field}.Scale", "must be positive and finite");
            if (block.Convention is not (RandomInterceptConvention.ObservationCentred or RandomInterceptConvention.WithinParent))
                throw Bad($"{field}.Convention", "is not a known convention");
            if (block.Convention == RandomInterceptConvention.WithinParent && string.IsNullOrWhiteSpace(block.ParentGrouping))
                throw Bad($"{field}.ParentGrouping", "must name the parent of a WithinParent block");
            if (block.Convention == RandomInterceptConvention.ObservationCentred && block.ParentGrouping is not null)
                throw Bad($"{field}.ParentGrouping", "is only defined for the WithinParent convention");
            if (block.ParentGrouping == block.Grouping) throw Bad($"{field}.ParentGrouping", "cannot name the block's own grouping");
            if (block.ScalePrior is { } prior && (!(prior.Scale > 0) || !double.IsFinite(prior.Scale)))
                throw Bad($"{field}.ScalePrior.Scale", "must be positive and finite");
        }
        var parents = blocks.ToDictionary(block => block.Grouping, block => block.ParentGrouping, StringComparer.Ordinal);
        foreach (BassRandomIntercept block in blocks)
        {
            var path = new HashSet<string>(StringComparer.Ordinal);
            string? next = block.Grouping;
            while (next is not null && parents.TryGetValue(next, out string? parent))
            {
                if (!path.Add(next)) throw Bad("RandomIntercepts", $"parent groupings form a cycle at '{next}'");
                next = parent;
            }
        }
    }

    internal static void CheckRandomSlopes(IReadOnlyList<BassRandomSlope>? blocks)
    {
        if (blocks is null) return;
        var names = new HashSet<(string Grouping, string Input)>();
        for (int k = 0; k < blocks.Count; k++)
        {
            string field = $"RandomSlopes[{k}]";
            BassRandomSlope block = blocks[k] ?? throw Bad(field, "is null");
            if (string.IsNullOrWhiteSpace(block.Grouping)) throw Bad($"{field}.Grouping", "must name a grouping");
            if (string.IsNullOrWhiteSpace(block.Input)) throw Bad($"{field}.Input", "must name a continuous predictor");
            if (!names.Add((block.Grouping, block.Input))) throw Bad(field, "duplicates a grouping/input pair");
            if (block.InputScale is not (RandomSlopeScale.Unit or RandomSlopeScale.Standardized)) throw Bad($"{field}.InputScale", "is not a known scale");
            if (!(block.Scale > 0) || !double.IsFinite(block.Scale)) throw Bad($"{field}.Scale", "must be positive and finite");
            if (block.ScalePrior is { } prior && (!(prior.Scale > 0) || !double.IsFinite(prior.Scale)))
                throw Bad($"{field}.ScalePrior.Scale", "must be positive and finite");
        }
    }

    private static void CheckTempering(TemperingSpec t, string field, TauPrior scale)
    {
        if (t.Path is not (TemperingPath.Likelihood or TemperingPath.Powered)) throw Bad($"{field}.Path", "is not a known path");
        if (t.Form is not (TemperedForm.Integrated or TemperedForm.Linear)) throw Bad($"{field}.Form", "is not a known form");
        if (t.Path == TemperingPath.Powered && t.Form == TemperedForm.Linear)
            throw Bad($"{field}.Form", "is the likelihood path's; the powered path powers the untempered marginal whole (D05)");
        if (t.ReferenceRung && t.Path != TemperingPath.Likelihood)
            throw Bad($"{field}.ReferenceRung", "needs the likelihood path, whose b = 0 member is the prior; the powered path's is not normalizable (D05)");
        if (t.ReferenceRung && !scale.IsProper)
            throw Bad($"{field}.ReferenceRung", $"needs a proper τ prior, whose exact draw it takes, and Coefficients.Scale is {scale}; take a > 2");
        int minimum = TemperingSpec.MinimumLevels(t.ReferenceRung);
        if (t.ListFault() is { } fault) throw Bad($"{field}.InverseTemperatures", fault);
        if (t.Levels < minimum)
            throw Bad($"{field}.Levels", t.ReferenceRung ? "must be at least 3 with a reference rung: b = 1, the hot end and b = 0" : "must be at least 2");
        if (t.DeclaredLevels() is null && !t.Adapts())
            throw Bad($"{field}.Adaptive", "must be true for a self-sizing ladder (Levels null), which sizes itself as it adapts; a fixed ladder declares Levels or InverseTemperatures");
        if (t.MaxLevels < minimum) throw Bad($"{field}.MaxLevels", $"must be at least {minimum}, the fewest rungs a ladder takes");
        if (t.MaxLevels < t.DeclaredLevels()) throw Bad($"{field}.MaxLevels", "must not lie below the declared rungs");
        if (t.HotEnd is double hot && !(hot > 0.0 && hot < 1.0)) throw Bad($"{field}.HotEnd", "must lie in (0, 1)");
        if (t.Schedule != SwapSchedule.DeterministicEvenOdd)
            throw Bad($"{field}.Schedule", "is DEO, the schedule D05 decided; another joins with its decision");
        if (t.Split is { } split)
        {
            if (!(split.Threshold > 0.0 && split.Threshold < 1.0)) throw Bad($"{field}.Split.Threshold", "must lie in (0, 1)");
            if (t.DeclaredLevels() is not null) throw Bad($"{field}.Split", "needs a self-sizing ladder (Levels and InverseTemperatures null), since it adds rungs");
        }
    }

    private static ArgumentException Bad(string field, string reason) => new($"{field} {reason}.", "config");
}

/// <summary>The term space and its support (D01, D06).</summary>
public sealed record BassStructureSpec
{
    /// <summary>The hinges' power q, a positive integer; 1, R's default.</summary>
    public int Degree { get; init; } = 1;

    /// <summary>The cap on a term's continuous factors; resolved to at most the continuous predictors' count.</summary>
    public int MaxInteraction { get; init; } = 3;

    /// <summary>The cap on a term's categorical factors; resolved to at most the categorical predictors' count.</summary>
    public int MaxCategoricalInteraction { get; init; } = 3;

    /// <summary>The cap on a term's functional factors; resolved to at most the output variables' count.</summary>
    public int MaxFunctionalInteraction { get; init; } = 3;

    /// <summary>The cap on M; null for min(1000, N − 1), since the intercept counts toward the dimension bound.</summary>
    public int? MaxBasis { get; init; }

    /// <summary>The support rule (D01, D06).</summary>
    public SupportRule Support { get; init; } = new();

    /// <summary>D06's s_min: every column's angle to the span of the others must reach it.</summary>
    public double AngleFloor { get; init; } = 1e-10;

    /// <summary>A declared cap on every coefficient's variance inflation, a statistical support restriction; off by default.</summary>
    public double? MaxVarianceInflation { get; init; }

    /// <summary>The categorical subset measure: R's size-uniform by default, the paper's uniform selectable (D01).</summary>
    public CategoricalSubsetMeasure CategoricalSubsets { get; init; } = CategoricalSubsetMeasure.SizeUniform;

    /// <summary>w₁, the pseudo-count on interaction orders in births; positive, so every death has a reverse birth.</summary>
    public double OrderWeight { get; init; } = 5;

    /// <summary>w₂, the pseudo-count on variables in births; positive, for the same reason.</summary>
    public double VariableWeight { get; init; } = 5;

    /// <summary>
    /// A polynomial-chaos basis in place of the hinges (E19): Legendre terms over the continuous predictors, of order at most
    /// <see cref="MaxInteraction"/> and degree at most its <see cref="LegendreBasis.MaxDegree"/>, with KHAOS's proposal. Null, the
    /// default, is BASS's hinges; left out of the serialized form when null, so a hinge config's form and fingerprint stand as
    /// they were.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public LegendreBasis? Legendre { get; init; }
}

/// <summary>What a term's support is measured by (D06).</summary>
public enum SupportKind
{
    /// <summary>Its nonzero points, every source's rule (R's <c>npart</c>): the default.</summary>
    Count,

    /// <summary>Its column mass Σᵢ B(xᵢ), hinges normalized to a maximum of one: BARS's Q6 rule, an option.</summary>
    ColumnMass,
}

/// <summary>
/// The groups of monitored quantities an R̂ stop waits on (D10): every group by default, or the predictive group alone.
/// The verdict reports every group either way (<see cref="BassDiagnostics.Verdict"/>).
/// </summary>
public enum BassStoppingScope
{
    /// <summary>Every group: predictive, level and structural, so the chains agree on the surface and the noise, on the
    /// posterior level, and on the representation. The default.</summary>
    All,

    /// <summary>The predictive group alone, σ² and the fitted values, for an emulator user who needs the surface and the
    /// noise and accepts a reported structural disagreement. A stop on it can end a run whose chains sit at different
    /// posterior levels, so it is a declared choice.</summary>
    Predictive,
}

/// <summary>The measure on a categorical factor's level subset (D01).</summary>
public enum CategoricalSubsetMeasure
{
    /// <summary>The subset's size uniform on 1 … n_v − 1, then the subset uniform among those of that size: R's, the default.</summary>
    SizeUniform,

    /// <summary>Uniform over the 2^{n_v} − 2 proper subsets: Francom et al. 2019's statement, selectable.</summary>
    Uniform,
}

/// <summary>
/// The support rule: a term needs at least <see cref="Minimum"/> by its <see cref="Kind"/> in its design part, and
/// <see cref="MinimumFunctional"/> in its functional part, counted apart (R's <c>npart</c> and <c>npart.func</c>).
/// Null resolves to R's default, min(20, 0.1·n) for that part's n: the cases, or the output locations.
/// </summary>
public sealed record SupportRule(SupportKind Kind = SupportKind.Count, double? Minimum = null, double? MinimumFunctional = null);

/// <summary>D02's count prior: λ ~ Gamma(h₁, rate h₂), integrated out; h₁ = h₂ = 10, R's and Nott, Kuk and Duc's.</summary>
public sealed record BassCountPrior(double Shape = 10, double Rate = 10);

/// <summary>
/// The coefficient prior (D03; its two axes decided 2026-10-29): a covariance <see cref="Form"/> and τ's prior,
/// <see cref="Scale"/>, from the shared family, which samples τ or holds it.
/// <list type="bullet">
/// <item><see cref="GForm"/>, σ²G_c⁻¹/τ, with its unit information at τ₀ = 1/N; <see cref="RidgeForm"/>, σ²I/τ on the
/// centred columns at unit root mean square, with τ₀ = 1.</item>
/// <item><see cref="Scale"/>: the unit-information Gamma (a = ½, Zellner–Siow, the default), hyper-g, hyper-g/n, each on
/// the form's τ₀, or a held τ. A held τ left null is τ₀ under the g-form, unit information (D03's second mode); under
/// the ridge it must be declared, its scale part of the model.</item>
/// </list>
/// The intercept, a unit-RMS column already, sits inside either form under level centring and the zero centre, so its
/// prior precision is τN under the g-form and τ under the ridge, unit information at τ₀ in both.
/// </summary>
public sealed record BassCoefficientPrior(CoefficientForm Form, TauPrior Scale)
{
    /// <summary>The g-form with τ's prior <paramref name="scale"/>; null for Zellner–Siow's unit-information Gamma, the default.</summary>
    public static BassCoefficientPrior G(TauPrior? scale = null) => new(new GForm(), scale ?? TauPrior.UnitInformationGamma());

    /// <summary>The normalized ridge with τ's prior <paramref name="scale"/>, a held one with its τ declared.</summary>
    public static BassCoefficientPrior Ridge(TauPrior scale) => new(new RidgeForm(), scale);

    /// <summary>KHAOS's raw ridge with τ's prior <paramref name="scale"/>, a held one with its τ declared: <c>khaos</c>'s τ² is 1/τ.</summary>
    public static BassCoefficientPrior RawRidge(TauPrior scale) => new(new RawRidgeForm(), scale);

    /// <summary>KHAOS's modified g-prior with complexity exponent <paramref name="zeta"/> (<see cref="ModifiedGForm"/>) and τ's prior
    /// <paramref name="scale"/>, where τ = 1/g₀²: <c>khaos</c>'s g₀² ~ IG(a, b) is <see cref="TauPrior.Gamma"/>(a, b).</summary>
    public static BassCoefficientPrior ModifiedG(TauPrior scale, double zeta = 1.0) => new(new ModifiedGForm(zeta), scale);

    /// <summary>τ₀ = 1/n₀'s n₀ for <paramref name="observations"/> observations: N for the g-form, 1 for either ridge.</summary>
    public double UnitObservations(double observations) => Form is RidgeForm or RawRidgeForm ? 1.0 : observations;
}

/// <summary>The coefficient prior's covariance form (D03, decided 2026-10-29).</summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(GForm), "G")]
[JsonDerivedType(typeof(RidgeForm), "Ridge")]
[JsonDerivedType(typeof(RawRidgeForm), "RawRidge")]
[JsonDerivedType(typeof(ModifiedGForm), "ModifiedG")]
public abstract record CoefficientForm
{
    private protected CoefficientForm() { }
}

/// <summary>The g-form, a | σ², τ ~ N(μ₀, σ²G_c⁻¹/τ) on the centred columns' Gram: shrinkage follows the design's
/// own information, and the marginal depends on the column space alone.</summary>
public sealed record GForm : CoefficientForm;

/// <summary>
/// The normalized ridge, β | σ², τ ~ N(0, σ²I/τ) on the centred columns scaled to unit root mean square over the N
/// observations: τ = 1 is one observation's worth for every term, whatever its knot. In the design's raw centred
/// coordinates it is the generalized ridge τD, D = diag(G_c)/N, so the marginal does not depend on the columns' scale.
/// </summary>
public sealed record RidgeForm : CoefficientForm;

/// <summary>
/// KHAOS's raw ridge, (β₀, β) | σ², τ ~ N(0, σ²I/τ) on the columns as evaluated, uncentred and unscaled, the intercept's
/// included under the zero centre (<c>khaos</c>'s <c>adaptive_khaos_ridge</c>, whose τ² is 1/τ); under level centring β₀'s
/// prior centres on the level instead, and a flat intercept leaves β₀ out. A Legendre basis is orthonormal under the uniform
/// measure, so its columns sit near unit root mean square and the two ridges nearly agree there; on hinges they do not. The
/// design carries it as a prior apart from the Gram (<see cref="Maths.Regression.GPrior.PriorPencil"/>), D = I on the raw
/// coefficients, so each structure costs its own decomposition.
/// </summary>
public sealed record RawRidgeForm : CoefficientForm;

/// <summary>The run's settings.</summary>
public sealed record BassSampling
{
    /// <summary>Independent chains.</summary>
    public int Chains { get; init; } = 4;

    /// <summary>The seed every chain's streams derive from.</summary>
    public int MasterSeed { get; init; }

    /// <summary>Sweeps discarded per chain; R's 9,000 of 10,000, provisional until BASS's own measurement. With
    /// <see cref="Samples"/> a budget, the most a run spends (D10).</summary>
    public int BurnIn { get; init; } = 9000;

    /// <summary>Draws kept per chain; R's 1,000, provisional likewise.</summary>
    public int Samples { get; init; } = 1000;

    /// <summary>Sweeps per kept draw.</summary>
    public int Thin { get; init; } = 1;

    /// <summary>0, no early stop; above 1, stop once R̂ and ESS meet their targets in the groups
    /// <see cref="StoppingScope"/> names, as BARS's stop reads them. The verdict reads this target, or 1.01 when it is 0
    /// (<see cref="ChainEnsemble.ReferenceRHatTarget"/>).</summary>
    public double RHatTarget { get; init; }

    /// <summary>Draws per chain between stopping checks; 0 for one round of <see cref="Samples"/>.</summary>
    public int BatchSize { get; init; }

    /// <summary>The ESS per chain the early stop also needs.</summary>
    public double MinEssPerChain { get; init; } = ChainEnsemble.MinEssPerChain;

    /// <summary>The groups an R̂ stop waits on (D10): every group by default (<see cref="BassStoppingScope"/>).</summary>
    public BassStoppingScope StoppingScope { get; init; } = BassStoppingScope.All;

    /// <summary>The move-type schedule: equal over the available moves by default (D01).</summary>
    public MoveSchedule Moves { get; init; } = new EqualOverAvailable();

    /// <summary>Local steps inside a change (F-M2's lever 2); null, the default, keeps D01's redraw alone.</summary>
    public LocalSteps? LocalSteps { get; init; }

    /// <summary>The fixed-M factor moves in the change's slot (F-M2's Phase 1); null, the default, keeps the change
    /// alone.</summary>
    public FactorMoves? FactorMoves { get; init; }

    /// <summary>
    /// A mixture on the births (E4; <see cref="BirthMixture"/>): components beside the fit's own birth, which keeps a
    /// positive share, so every structure it reaches stays reachable. Null, the default, runs the fit's own birth alone.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public BirthMixture? Births { get; init; }

    /// <summary>Parallel tempering within each chain; null, the default, runs untempered (D05).</summary>
    public TemperingSpec? Tempering { get; init; }

    /// <summary>
    /// Fitting by SMC from the prior to the posterior (stage 8, D05) in place of MCMC; null, the default, runs MCMC. It
    /// refuses <see cref="Tempering"/>, since SMC tempers its own sequence, and an <see cref="RHatTarget"/>, since it ends
    /// at b = 1; <see cref="Chains"/>, <see cref="BurnIn"/>, <see cref="Samples"/>, <see cref="Thin"/>,
    /// <see cref="BatchSize"/>, <see cref="MinEssPerChain"/> and <see cref="StoppingScope"/> are not read under it, as the
    /// resolved configuration logs.
    /// </summary>
    public SmcSpec? Smc { get; init; }

    /// <summary>
    /// Fitting by the flat-histogram sampler (N17) in place of MCMC: SAMC learns multicanonical weights over ℓ₁ in burn-in,
    /// independent chains run on the frozen weights, and their draws weighted at b = 1 are the fit; null, the default, runs
    /// MCMC. It refuses <see cref="Smc"/>, <see cref="Tempering"/> (its pilot declares its own ladder) and an
    /// <see cref="RHatTarget"/>; <see cref="Chains"/>, <see cref="BurnIn"/>, <see cref="Samples"/>,
    /// <see cref="BatchSize"/>, <see cref="MinEssPerChain"/> and <see cref="StoppingScope"/> are not read under it, the
    /// spec's own lengths and chains standing in. <see cref="Thin"/> is read: each iteration of the run, a learning update
    /// or a production draw, follows that many sweeps of every chain.
    /// </summary>
    public BassFlatHistogramSpec? FlatHistogram { get; init; }

    /// <summary>Whether a run that can give an evidence reports it (D03, D05).</summary>
    public bool ReportEvidence { get; init; } = true;

    /// <summary>
    /// The training points whose fitted values the diagnostics monitor (R̂ and ESS): a case and, for a direct functional
    /// response, an output location (null takes the middle one). Null takes the defaults
    /// (<see cref="BassDiagnostics.DefaultPoints(int, int?)"/>), logged as an auto-pick (decided 2026-10-28).
    /// </summary>
    public IReadOnlyList<BassMonitoredPoint>? MonitoredPoints { get; init; }
}

/// <summary>
/// Parallel tempering within each chain (D05, D10): a ladder of replicas of D04's sweep, the b = 1 replica the sampler
/// (<see cref="BassLadder"/>).
/// <list type="bullet">
/// <item><see cref="Path"/>: the likelihood path, prior · p_b, by default; the powered path, (prior · p₁)^b with τ powered
/// in log τ, selectable.</item>
/// <item><see cref="Form"/>: the likelihood path's tempered marginal, the integrated form p_b by default or the linear
/// form b·ln p₁. The linear form reaches b = 0 under the improper options too.</item>
/// <item><see cref="Levels"/> and <see cref="HotEnd"/>: the replicas, geometric from 1 down to the hot end, then b = 0
/// with the reference rung. An <see cref="Adaptive"/> ladder (the default) re-spaces its interior rungs during burn-in to
/// equalize the swap rejections, from that start, and is frozen after it (D10). Its ends never move: b = 1 and the hot
/// end, or b = 0 with the reference rung, above which the hot end is an interior rung and moves too.</item>
/// <item>By default the ladder sizes itself (decided 2026-10-28): from <see cref="InitialLevels"/> it grows or shrinks
/// during burn-in to Syed et al.'s N ≈ 2Λ pairs by the barrier it reads (<see cref="TemperingLadder.Size"/>), within
/// <see cref="MaxLevels"/>, and settles a round before burn-in ends, so the last round equalizes it; and its hot end is
/// auto-picked on the data (<see cref="BassResolved.HotEnd"/>).</item>
/// <item><see cref="InverseTemperatures"/>: the rungs themselves, R's <c>temp.ladder</c>, in place of the count and the hot
/// end. A list declares its spacing, so it runs as given, as R's does, unless <see cref="Adaptive"/> is set true, which
/// re-spaces its interior from there.</item>
/// <item><see cref="Schedule"/>: DEO, the non-reversible even-odd swaps (Syed et al. 2022), swapped after every sweep.</item>
/// <item><see cref="Split"/>: the split rule (F-M2's lever 4b) for a self-sizing ladder, opt-in.</item>
/// <item><see cref="ReferenceRung"/>: the hot replica draws exactly from the prior on F every scan (OQ2), opt-in. A run
/// reports its evidence only with it (D03), under the proper priors.</item>
/// </list>
/// </summary>
public sealed record TemperingSpec
{
    /// <summary>The rungs a self-sizing ladder starts from, b = 1 and, with the reference rung, b = 0 included.</summary>
    public const int InitialLevels = 8;

    /// <summary>The tempering path: likelihood by default (D05).</summary>
    public TemperingPath Path { get; init; } = TemperingPath.Likelihood;

    /// <summary>The likelihood path's tempered form: integrated by default (D05).</summary>
    public TemperedForm Form { get; init; } = TemperedForm.Integrated;

    /// <summary>The replicas, b = 1 and, with the reference rung, b = 0 included: a count fixes them, and null (the default)
    /// sizes the ladder during burn-in, from <see cref="InitialLevels"/> within <see cref="MaxLevels"/>, which needs
    /// <see cref="Adaptive"/>.</summary>
    public int? Levels { get; init; }

    /// <summary>The most replicas a self-sizing ladder takes, a compute budget: 64 by default. A ladder that reaches it is
    /// flagged (<see cref="BassTemperingReadout.AtMaxLevels"/>).</summary>
    public int MaxLevels { get; init; } = 64;

    /// <summary>The ladder's lowest positive b, in (0, 1): its hot end, or with the reference rung the rung above b = 0, where
    /// an adaptive ladder starts it. Null (the default) auto-picks it on the data (<see cref="BassResolved.HotEnd"/>).</summary>
    public double? HotEnd { get; init; }

    /// <summary>Whether the ladder re-spaces its interior rungs during burn-in; it is frozen after it either way (D10). Null
    /// (the default) reads by what is declared (<see cref="Adapts"/>): a built ladder adapts, and a list runs as given.</summary>
    public bool? Adaptive { get; init; }

    /// <summary>The swap schedule: DEO (D05).</summary>
    public SwapSchedule Schedule { get; init; } = SwapSchedule.DeterministicEvenOdd;

    /// <summary>Whether the hot replica sits at b = 0 and draws exactly from the prior on F; off by default (D05).</summary>
    public bool ReferenceRung { get; init; }

    /// <summary>A self-sizing ladder's split rule (F-M2's lever 4b), which adds a rung inside each pair that rejects past
    /// its threshold; null, the default, sizes by the barrier alone.</summary>
    public SplitRule? Split { get; init; }

    /// <summary>The ladder's rungs themselves, strictly descending from exactly b = 1 and ending above b = 0, or at exactly
    /// b = 0 with the reference rung; null (the default) builds them from <see cref="Levels"/> and <see cref="HotEnd"/>, which
    /// it replaces. R's <c>temp.ladder</c>.</summary>
    public double[]? InverseTemperatures { get; init; }

    /// <summary>The fewest rungs a ladder takes: b = 1 and the hot end, and b = 0 with the reference rung.</summary>
    public static int MinimumLevels(bool referenceRung) => referenceRung ? 3 : 2;

    /// <summary>The replicas declared: the list's, or <see cref="Levels"/>; null when the ladder sizes itself.</summary>
    public int? DeclaredLevels() => InverseTemperatures?.Length ?? Levels;

    /// <summary>The lowest positive b declared: the list's, or <see cref="HotEnd"/>; null when the data picks it.</summary>
    public double? DeclaredHotEnd() => InverseTemperatures is { } b ? b[ReferenceRung ? ^2 : ^1] : HotEnd;

    /// <summary>Whether the ladder re-spaces during burn-in: <see cref="Adaptive"/>, or when null, false for a list, which
    /// declares its spacing, and true for a built ladder.</summary>
    public bool Adapts() => Adaptive ?? InverseTemperatures is null;

    /// <summary>The rungs a chain's ladder starts on: those declared, or a self-sizing ladder's
    /// <see cref="InitialLevels"/> within [<see cref="MinimumLevels"/>, <see cref="MaxLevels"/>].</summary>
    public int StartLevels() => DeclaredLevels() ?? Math.Clamp(InitialLevels, MinimumLevels(ReferenceRung), MaxLevels);

    // Why the list cannot run, or null: BassConfig.Validate and the ladder both refuse it.
    internal string? ListFault()
    {
        if (InverseTemperatures is not { } b) return null;
        if (Levels is not null) return "replaces Levels, so it takes none";
        if (HotEnd is not null) return "replaces HotEnd, so it takes none";
        int minimum = MinimumLevels(ReferenceRung);
        if (b.Length < minimum) return $"must hold at least {minimum} rungs";
        if (b[0] != 1.0) return "must start at exactly b = 1";
        for (int i = 1; i < b.Length; i++)
            if (!(b[i] < b[i - 1])) return "must descend strictly";
        if (ReferenceRung) return b[^1] == 0.0 ? null : "must end at exactly b = 0 with the reference rung";
        return b[^1] > 0.0 ? null : "must end above b = 0 without the reference rung";
    }
}

/// <summary>
/// Fitting by SMC (stage 8; D05, decided 2026-10-27; its settings decided 2026-10-29): <see cref="Replicates"/> independent
/// populations of <see cref="Particles"/>, each tempered from exact draws of the prior on F at b = 0 to the posterior at
/// b = 1 on the likelihood path. Each next b is placed so the reweighting's conditional ESS is
/// <see cref="TargetConditionalEss"/>·N (Zhou, Johansen and Aston 2016); each particle is reweighted by ℓ_{b′} − ℓ_b,
/// the population resampled when its ESS/N falls below <see cref="ResampleBelow"/>, and each particle mutated by
/// <see cref="MutationSteps"/> of D04's sweeps at b (<see cref="BassSmcKernel"/>). Each replicate's final population, with
/// its weights, is the fit's draws, and under D03's proper defaults its normalizer is the evidence.
/// <list type="bullet">
/// <item>The values are provisional per-engine defaults, read in the stage 8 qualification and not tuned.</item>
/// <item><see cref="Resampling"/>: multinomial (D05), whose offspring the normalizer's genealogy error describes; the
/// lower-variance schemes stay selectable and leave that error NaN.</item>
/// <item><see cref="Form"/>: the likelihood path's tempered form, integrated by default. Under D03's improper options (the
/// Jeffreys σ² prior, the flat intercept) the integrated form has no b = 0 member, so the run needs the linear form, and
/// its normalizer is then no evidence.</item>
/// <item><see cref="MaxStages"/>: each replicate's stage budget. A run that spends it before b = 1 is returned flagged
/// (<see cref="BassSmcReadout.ReachedPosterior"/>), its posterior views refused, and its checkpoint continues it.</item>
/// <item><see cref="TraceStages"/>: record each stage's population, for stepping-stone over the stages and, on the linear
/// form, path sampling; off by default, when the evidence is the normalizer alone.</item>
/// <item><see cref="AdaptiveMutation"/>: each stage's sweeps chosen from the fraction of particles whose structure the previous
/// stage's sweeps changed, within <see cref="MutationSteps"/> per stage taken; null keeps the fixed count.</item>
/// <item><see cref="StepCap"/>: a bound on each step after the first, b′ ≤ r·b, so the schedule cannot jump from a collapsed
/// population to b = 1 (the trap measured on KHAOS's raw ridge); null leaves the conditional-ESS rule alone.</item>
/// </list>
/// </summary>
public sealed record SmcSpec
{
    /// <summary>The population per replicate, N.</summary>
    public int Particles { get; init; } = 1000;

    /// <summary>Independent SMC runs, run in parallel within the run's bound (<see cref="BassRunOptions.MaxDegreeOfParallelism"/>):
    /// their spread is the between-run reading, and the evidence's between-run error.</summary>
    public int Replicates { get; init; } = 4;

    /// <summary>The adaptive schedule's conditional-ESS fraction.</summary>
    public double TargetConditionalEss { get; init; } = 0.99;

    /// <summary>D04's sweeps per particle at each stage; 0 makes the run sequential importance sampling.</summary>
    public int MutationSteps { get; init; } = 5;

    /// <summary>Resample when ESS/N falls below this fraction.</summary>
    public double ResampleBelow { get; init; } = 0.5;

    /// <summary>The resampling scheme: multinomial (D05).</summary>
    public ResamplingScheme Resampling { get; init; } = ResamplingScheme.Multinomial;

    /// <summary>The likelihood path's tempered form: integrated by default (D05).</summary>
    public TemperedForm Form { get; init; } = TemperedForm.Integrated;

    /// <summary>The stages a replicate may take after its initial draw.</summary>
    public int MaxStages { get; init; } = 1000;

    /// <summary>Whether each stage's population is recorded, for the evidence over the stages.</summary>
    public bool TraceStages { get; init; }

    /// <summary>Each stage's sweeps chosen from the previous stage's moved fraction, within <see cref="MutationSteps"/> per stage
    /// taken (<see cref="Maths.Samplers.Smc.AdaptiveMutation"/>); null keeps the fixed count. Left out of the serialized form when
    /// null, so a fixed-count config's form and fingerprint stand as they were.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public AdaptiveMutation? AdaptiveMutation { get; init; }

    /// <summary>
    /// r, a bound on each step of the schedule from b &gt; 0: b′ = min(the conditional-ESS rule's b′, r·b, 1), r &gt; 1. The first step,
    /// from b = 0, is the rule's alone. Null, the default, leaves the rule alone. Where it binds, the fit's draws differ from the
    /// uncapped fit's, and it takes at least ln(1/b₁)/ln r stages, b₁ the first step's b.
    /// <para>Opt-in, built 2026-11-09 for KHAOS's SMC trap: there the rule jumps to b = 1 from a population collapsed at near-empty
    /// structures, and r = 1.25 with E4's residual birth reached the posterior in 16 of 16 replicates (the proposal-kernels brief's
    /// slice 2), the cap alone in 11 of 16. Their evidence still spread over 21 nats. Left out of the serialized form when null.</para>
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? StepCap { get; init; }
}

/// <summary>
/// Fitting by the flat-histogram sampler (N17; decided 2026-10-29, Phase 2 2026-10-30): the engine's
/// <see cref="Maths.Samplers.FlatHistogram.MulticanonicalRun"/> over D04's sweep on the multicanonical target
/// π_w ∝ prior(structure, τ)·w(ℓ₁), ℓ₁ the collapsed log marginal at b = 1, with the prior over (structure, τ) as the base
/// measure. SAMC learns the weights in burn-in; then <see cref="FlatHistogramSpec.Chains"/> independent chains run on the
/// frozen weights, and their draws weighted by e^{ℓ₁}/w(ℓ₁) are the fit (<see cref="BassFit.DrawWeights"/>).
/// <list type="bullet">
/// <item><see cref="Run"/>: the engine's spec, its chains, lengths, start, window, b_lo, K and gain. Its lengths are totals,
/// in iterations of <see cref="BassSampling.Thin"/> sweeps per chain.</item>
/// <item><see cref="Pilot"/>: the tempering pilot a WHAM start reads, and the window, b_lo or K the spec leaves to it,
/// declared exactly when the spec needs one (<see cref="FlatHistogramSpec.NeedsPilot"/>). With its reference rung the
/// pilot's hot end is b = 0, and so is b_lo.</item>
/// <item><see cref="MaxIterations"/>: the iterations a call takes, learning updates and production draws per chain
/// together, as SMC's stage budget; null runs to the lengths. A fit short of them is flagged
/// (<see cref="BassFlatHistogramReadout.Complete"/>), its posterior views refused, and its checkpoint continues it.</item>
/// </list>
/// The evidence is ln Z(1) by the engine's ratio estimator, reported only from b_lo = 0 under proper priors; above b_lo = 0
/// the fit reports ln Z(1) − ln Z(b_lo), flagged as a ratio (<see cref="BassFlatHistogramEvidence"/>).
/// </summary>
public sealed record BassFlatHistogramSpec
{
    /// <summary>The engine's spec.</summary>
    public required FlatHistogramSpec Run { get; init; }

    /// <summary>The tempering pilot; null when <see cref="Run"/> reads none.</summary>
    public BassFlatHistogramPilot? Pilot { get; init; }

    /// <summary>The iterations one call takes at most, ≥ 1; null runs to the lengths.</summary>
    public long? MaxIterations { get; init; }
}

/// <summary>
/// The flat-histogram run's tempering pilot (N17 Phase 2): <see cref="Chains"/> ladders of D04's sweep, each from an exact
/// draw of the prior on F, burned as a tempered fit burns them (<see cref="BassLadder"/>, adapted and, self-sizing, sized),
/// then frozen and traced for <see cref="TracedScans"/> scans. Each ladder's rungs combine by WHAM
/// (<see cref="Maths.Samplers.FlatHistogram.TemperingPilot.FromLadders"/>) into ln ĝ and the window. Its ladder is on the
/// likelihood path's linear form, ℓ_b = b·ℓ₁, whose rung laws prior·e^{b·ℓ₁} WHAM needs; by default it sizes itself and
/// has the exact reference rung at b = 0.
/// </summary>
public sealed record BassFlatHistogramPilot
{
    /// <summary>The ladder: the linear form with the reference rung, self-sizing, by default.</summary>
    public TemperingSpec Tempering { get; init; } = new() { Form = TemperedForm.Linear, ReferenceRung = true };

    /// <summary>The ladders, one per chain.</summary>
    public int Chains { get; init; } = 4;

    /// <summary>Each ladder's burn-in scans, ≥ 0, in which it adapts and sizes itself.</summary>
    public required int BurnIn { get; init; }

    /// <summary>Each ladder's traced scans after burn-in, ≥ 2: every rung's draws.</summary>
    public required int TracedScans { get; init; }
}

/// <summary>How a state's move types are chosen (D01).</summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(EqualOverAvailable), "EqualOverAvailable")]
[JsonDerivedType(typeof(PriorCancelling), "PriorCancelling")]
public abstract record MoveSchedule
{
    private protected MoveSchedule() { }
}

/// <summary>Equal probability over the moves available at M (Nott, Kuk and Duc): the default.</summary>
public sealed record EqualOverAvailable : MoveSchedule;

/// <summary>
/// Green's prior-cancelling schedule, b_M = C·min(1, r_M) and d_{M+1} = C·min(1, 1/r_M), with r_M D02's count ratio, so
/// the count prior cancels from births and deaths; an option for F-M1's comparison.
/// </summary>
public sealed record PriorCancelling(double C = 0.4) : MoveSchedule;

/// <summary>
/// Local steps inside a change (F-M2's lever 2, decided 2026-10-29): with probability <see cref="Probability"/> (ρ) a
/// change steps its factor locally instead of redrawing it from the prior. A hinge, continuous or functional, moves its
/// knot j′ = j + δ over its variable's distinct values, its sign kept, with |δ| uniform on {1, …, w} and the sign of δ
/// ±1, w = max(1, ⌈ω·n_v⌉) for the fraction ω = <see cref="WindowFraction"/>; a categorical factor toggles one level,
/// uniform among its n_v. A step off the knots, or to an empty or full subset, is a null move. The proposal is symmetric,
/// so the ratio is Δℓ and, for a toggle, the subset measure's prior ratio.
/// </summary>
public sealed record LocalSteps(double Probability = 0.5, double WindowFraction = 0.02)
{
    /// <summary>
    /// w = max(1, ⌈ω·n⌉) over a variable of <paramref name="values"/> distinct values, ω·n taken in decimal so that a
    /// multiple of ω's decimal places does not round up a step.
    /// </summary>
    public int Window(int values)
    {
        if (values < 1) throw new ArgumentOutOfRangeException(nameof(values));
        return Math.Max(1, (int)Math.Ceiling((decimal)WindowFraction * values));
    }
}

/// <summary>
/// The fixed-M factor moves (F-M2's Phase 1, decided 2026-10-29), sharing the change's slot of the move schedule by these
/// relative weights, thirds by default: <see cref="Change"/>, a factor's redraw or local step; <see cref="PruneOrGrow"/>,
/// a factor pruned from a term's part or grown into it, ½ each; <see cref="VariableChange"/>, a factor's variable
/// exchanged for one the part does not use. The births' and deaths' shares are the schedule's alone.
/// </summary>
public sealed record FactorMoves(double Change = 1, double PruneOrGrow = 1, double VariableChange = 1);

/// <summary>
/// A mixture on the births (E4; the proposal-kernels brief): q = (1 − ω_H)·q_base + ω_H·q_H, the fit's own birth q_base
/// keeping the rest of the share. Each component's probability is exact, so the mixture's is, and the target is unchanged.
/// Deaths stay uniform, and changes and mutations stay the fit's own.
/// </summary>
public sealed record BirthMixture
{
    /// <summary>The history-based birth (<see cref="HistoryBirth"/>); null leaves it out.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public HistoryBirth? History { get; init; }

    /// <summary>The residual-informed birth (<see cref="ResidualBirth"/>); null leaves it out.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public ResidualBirth? Residual { get; init; }
}

/// <summary>
/// The history-based birth (E4; Nott and Kohn 2005, <i>Biometrika</i> 92): the fit's own birth with each variable weighted
/// by max(ĥ_u, <see cref="Floor"/>), ĥ_u the share of the chain's past states in which some term uses u, in place of the
/// current terms' counts (a hinge birth's z_u at every order, one included; a Legendre birth's (aⱼ^α + ε) in
/// <c>make_weights</c>). D10's policy: ĥ is refreshed every <see cref="Refresh"/> sweeps of burn-in, from the past only, and
/// frozen at its end, so production runs a fixed kernel. Under SMC each stage reads it from its weighted population.
/// </summary>
/// <param name="Share">ω_H, the component's share of the births, in (0, 1).</param>
/// <param name="Floor">δ, each weight's floor, in (0, 1]: Nott and Kohn's truncation.</param>
/// <param name="Refresh">Sweeps between refreshes in burn-in, at least 1.</param>
public sealed record HistoryBirth(double Share = 0.25, double Floor = 0.01, int Refresh = 100);

/// <summary>
/// The residual-informed birth (E4; the proposal-kernels brief's slice 2): a parent drawn uniformly over the intercept and the
/// structure's terms that can take one more continuous factor (a hinge, or a Legendre factor), then the new factor ∝ exp(½·b·Δℓ), Δℓ
/// its single column's profile-likelihood gain against the structure's least-squares residual, at the move's inverse temperature b.
/// Its probability is exact, a death's reverse reads the reduced structure's residual, and the fit's own birth keeps its share, so
/// every structure stays reachable. A scalar response only.
/// </summary>
/// <param name="Share">ω_R, the component's share of the births, in (0, 1).</param>
public sealed record ResidualBirth(double Share = 0.25);

/// <summary>
/// A self-sizing ladder's split rule (F-M2's lever 4b, decided 2026-10-29): each adaptation round, every adjacent pair
/// whose swap rejection exceeds <see cref="Threshold"/> gains a rung halfway along its barrier, within
/// <see cref="TemperingSpec.MaxLevels"/>.
/// </summary>
public sealed record SplitRule(double Threshold = 0.9);
