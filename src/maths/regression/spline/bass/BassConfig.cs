using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;
using Maths.Regression.GPrior;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A BASS fit's configuration (the contract shapes, decided 2026-10-28). Every field is a value, so the record
/// serializes to JSON whole and is its own manifest, and the run's fingerprint hashes that serialization with the data;
/// <c>BarsConfig</c>'s component instances fingerprint their types only (OQ10). The defaults are D01–D06's. The
/// auto-picks (the basis cap, the support thresholds, N, s² and the centre's level) are resolved at entry and reported
/// in <see cref="BassResolved"/>. A setting joins with the code that reads it: tempering with stage 4's S8, SMC with
/// stage 8, the ridge priors when they are built.
/// </summary>
public sealed record BassConfig
{
    /// <summary>The term space and its support (D01, D06).</summary>
    public BassStructureSpec Structure { get; init; } = new();

    /// <summary>The count prior, λ integrated out (D02).</summary>
    public BassCountPrior Count { get; init; } = new();

    /// <summary>The coefficient prior: the g-prior with τ sampled, by default, or with τ fixed (D03).</summary>
    public BassCoefficientPrior Coefficients { get; init; } = new SampledG();

    /// <summary>σ²'s prior: IG(c, c·s²) scaled to the response, c = 10⁻³, by default (D03).</summary>
    public Sigma2Prior Sigma2 { get; init; } = Sigma2Prior.ScaledToResponse(1e-3);

    /// <summary>
    /// The intercept convention: level centring (<see cref="GPriorCentre.Mean"/>) by default; Liang et al.'s flat
    /// intercept (<see cref="GPriorCentre.Flat"/>) and DiMatteo, Genovese and Kass's zero centre selectable (D03). The
    /// median centre is for the noise extensions' weights and is refused.
    /// </summary>
    public GPriorCentre Intercept { get; init; } = GPriorCentre.Mean;

    /// <summary>The run: chains, seeds, lengths, the R̂ stop and the move schedule.</summary>
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

        BassCountPrior count = Count ?? throw Bad("Count", "is missing");
        if (!(count.Shape > 0.0 && double.IsFinite(count.Shape))) throw Bad("Count.Shape", "must be positive and finite");
        if (!(count.Rate > 0.0 && double.IsFinite(count.Rate))) throw Bad("Count.Rate", "must be positive and finite");

        switch (Coefficients)
        {
            case SampledG g when !(g.TauShape > 0.0 && double.IsFinite(g.TauShape)):
                throw Bad("Coefficients.TauShape", "must be positive and finite");
            case FixedG { Tau: { } tau } when !(tau > 0.0 && double.IsFinite(tau)):
                throw Bad("Coefficients.Tau", "must be positive and finite");
            case SampledG or FixedG:
                break;
            case null:
                throw Bad("Coefficients", "is missing");
            default:
                throw Bad("Coefficients", "is not a known prior");
        }

        if (Sigma2 is null) throw Bad("Sigma2", "is missing");
        if (Intercept == GPriorCentre.Median) throw Bad("Intercept", "the median centre is for the noise extensions' weights (D03)");
        if (Intercept is not (GPriorCentre.Mean or GPriorCentre.Flat or GPriorCentre.Zero)) throw Bad("Intercept", "is not a known convention");

        BassSampling run = Sampling ?? throw Bad("Sampling", "is missing");
        if (run.Chains < 1) throw Bad("Sampling.Chains", "must be at least 1");
        if (run.BurnIn < 0) throw Bad("Sampling.BurnIn", "must not be negative");
        if (run.Samples < 1) throw Bad("Sampling.Samples", "must be at least 1");
        if (run.Thin < 1) throw Bad("Sampling.Thin", "must be at least 1");
        if (!(run.RHatTarget == 0.0 || (run.RHatTarget > 1.0 && double.IsFinite(run.RHatTarget))))
            throw Bad("Sampling.RHatTarget", "must be 0, no early stop, or a finite target above 1");
        if (run.BatchSize < 0) throw Bad("Sampling.BatchSize", "must not be negative");
        if (!(run.MinEssPerChain >= 0.0 && double.IsFinite(run.MinEssPerChain))) throw Bad("Sampling.MinEssPerChain", "must be finite and not negative");
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
}

/// <summary>What a term's support is measured by (D06).</summary>
public enum SupportKind
{
    /// <summary>Its nonzero points, every source's rule (R's <c>npart</c>): the default.</summary>
    Count,

    /// <summary>Its column mass Σᵢ B(xᵢ), hinges normalized to a maximum of one: BARS's Q6 rule, an option.</summary>
    ColumnMass,
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

/// <summary>The coefficient prior (D03): the g-prior with τ sampled, or with τ fixed.</summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(SampledG), "SampledG")]
[JsonDerivedType(typeof(FixedG), "FixedG")]
public abstract record BassCoefficientPrior
{
    private protected BassCoefficientPrior() { }
}

/// <summary>
/// The g-prior with τ sampled from the unit-information family, τ ~ Gamma(a, rate a·N), whose mean is 1/N at every
/// shape; a = ½ is Zellner–Siow, the default (D03).
/// </summary>
public sealed record SampledG(double TauShape = 0.5) : BassCoefficientPrior;

/// <summary>The g-prior with τ held at <paramref name="Tau"/>; null for unit information, 1/N (D03's second mode).</summary>
public sealed record FixedG(double? Tau = null) : BassCoefficientPrior;

/// <summary>The run's settings.</summary>
public sealed record BassSampling
{
    /// <summary>Independent chains.</summary>
    public int Chains { get; init; } = 4;

    /// <summary>The seed every chain's streams derive from.</summary>
    public int MasterSeed { get; init; }

    /// <summary>Sweeps discarded per chain; R's 9,000 of 10,000, provisional until BASS's own measurement.</summary>
    public int BurnIn { get; init; } = 9000;

    /// <summary>Draws kept per chain; R's 1,000, provisional likewise.</summary>
    public int Samples { get; init; } = 1000;

    /// <summary>Sweeps per kept draw.</summary>
    public int Thin { get; init; } = 1;

    /// <summary>0, no early stop; above 1, stop once R̂ and ESS meet their targets, as BARS's stop reads them.</summary>
    public double RHatTarget { get; init; }

    /// <summary>Draws per chain between stopping checks; 0 for one round of <see cref="Samples"/>.</summary>
    public int BatchSize { get; init; }

    /// <summary>The ESS per chain the early stop also needs.</summary>
    public double MinEssPerChain { get; init; } = ChainEnsemble.MinEssPerChain;

    /// <summary>The move-type schedule: equal over the available moves by default (D01).</summary>
    public MoveSchedule Moves { get; init; } = new EqualOverAvailable();

    /// <summary>Whether a run that can give an evidence reports it (D03, D05).</summary>
    public bool ReportEvidence { get; init; } = true;

    /// <summary>
    /// The training points whose fitted values the diagnostics monitor (R̂ and ESS): a case and, for a direct functional
    /// response, an output location (null takes the middle one). Null takes the defaults
    /// (<see cref="BassDiagnostics.DefaultPoints(int, int?)"/>), logged as an auto-pick (decided 2026-10-28).
    /// </summary>
    public IReadOnlyList<BassMonitoredPoint>? MonitoredPoints { get; init; }
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
