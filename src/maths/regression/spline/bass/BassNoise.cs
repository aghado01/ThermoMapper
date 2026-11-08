using System;
using System.Collections.Generic;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>Which noise a <see cref="BassNoise"/> declares; there, the parameters each kind takes.</summary>
public enum BassNoiseKind
{
    /// <summary>Student t (robust BASS; GBASS's <c>tbass</c>): v ~ IG(ν/2, ν/2), so ε = σ·√(c·v)·ζ is σ√c·t_ν.</summary>
    StudentT,

    /// <summary>The asymmetric Laplace (quantile regression; GBASS's <c>qbass</c>): v ~ Exp(1), c = 2/(q(1 − q)) and β = 1/q − 1/(1 − q),
    /// so ε is ALD(0, σ, q) and the surface f is the q-quantile of y (Kozumi and Kobayashi 2011).</summary>
    AsymmetricLaplace,

    /// <summary>The Normal-Wald, or normal inverse Gaussian (GBASS's <c>nwbass</c>): v ~ GIG(−½, γ², 1), the inverse Gaussian of mean 1/γ and
    /// shape 1, with the tail γ and the skew β each held or drawn, so the data can choose the noise's skew and tails.</summary>
    NormalWald,

    /// <summary>The generalized hyperbolic: v ~ GIG(p, a, b) as declared (GBASS's <c>gbass</c> with a GIG <c>v_prior</c>), its skew held or
    /// drawn. The Normal-Wald is p = −½, b = 1; the skewed t a = 0; the variance gamma b = 0.</summary>
    GeneralizedHyperbolic,

    /// <summary>The generalized beta prime mixture: v ~ GBP(p, a, b) as declared (GBASS's <c>gbass</c> with a GBP <c>v_prior</c>), its skew
    /// held or drawn; the horseshoe at p = 1, a = b = ½ (GBASS's <c>hbass</c>), the Strawderman–Berger at p = ½, a = 1, b = ½. Heavy enough
    /// at both ends of v that an observation is fitted closely or let go ("interpolate or ignore").</summary>
    GeneralizedBetaPrime,
}

/// <summary>
/// A normal prior N(<see cref="Mean"/>, <see cref="Sd"/>²) on one of a noise's scalars, the skew β or the Normal-Wald's tail γ; Sd = 0 holds
/// the scalar at the mean. γ's is the normal's part on γ &gt; 0, which is the conditional's support too (GBASS's draw takes the γ &gt; 0
/// piece's normal untruncated, E1's slice 0).
/// </summary>
public sealed record BassNoisePrior
{
    /// <summary>N(<paramref name="mean"/>, <paramref name="sd"/>²), the mean finite and the sd non-negative and finite.</summary>
    public BassNoisePrior(double mean, double sd)
    {
        if (!double.IsFinite(mean)) throw new ArgumentOutOfRangeException(nameof(mean), "The mean must be finite.");
        if (!(sd >= 0.0 && double.IsFinite(sd))) throw new ArgumentOutOfRangeException(nameof(sd), "The sd must be non-negative and finite.");
        Mean = mean;
        Sd = sd;
    }

    /// <summary>The scalar held at <paramref name="value"/>.</summary>
    public static BassNoisePrior Held(double value) => new(value, 0.0);

    /// <summary>The mean, and the held value at Sd = 0.</summary>
    public double Mean { get; }

    /// <summary>The sd; 0 holds the scalar.</summary>
    public double Sd { get; }

    /// <summary>Whether the scalar is held at the mean.</summary>
    [JsonIgnore]
    public bool IsHeld => Sd == 0.0;
}

/// <summary>What the priors read of the response under a noise whose working response moves (<see cref="BassPriorReference"/>).</summary>
public enum BassNoiseReference
{
    /// <summary>The mean and the mean square about it, D03's read: at β = 0, under the unweighted prior Gram, what Gaussian BASS reads.</summary>
    Moments,

    /// <summary>The median and the normal-consistent MAD², for heavy tails: an outlier cannot move the level or the scale.</summary>
    Robust,

    /// <summary>Nothing: for priors that read no response, a flat or zero intercept with a σ² prior that reads no s². Refused otherwise,
    /// at validation.</summary>
    None,

    /// <summary>The sample q-quantile and the squared mean check loss: the asymmetric Laplace's level and scale by maximum likelihood,
    /// whose f is the q-quantile and whose ȳ is not its level. The asymmetric Laplace's alone.</summary>
    Quantile,
}

/// <summary>
/// The noise a BASS fit declares beyond the Gaussian (E1): a normal variance-mean mixture, εᵢ = σ(β·vᵢ + √(c·vᵢ)·ζᵢ), ζᵢ ~ N(0, 1), with
/// vᵢ from the kind's mixing law, σ² the core's (GBASS's w) and c <see cref="VarianceScale"/>. The core fits y − σβv under weights
/// 1/(c·v) at the current latents (<see cref="VarianceMeanMixture"/>). A value, serialized as it stands, so a config that holds it
/// fingerprints whole; a config without one is Gaussian BASS and serializes as it did.
/// </summary>
public sealed record BassNoise
{
    /// <summary>
    /// The noise of kind <paramref name="kind"/>: <see cref="BassNoiseKind.StudentT"/> takes ν = <paramref name="degreesOfFreedom"/>,
    /// positive and finite, <see cref="BassNoiseKind.AsymmetricLaplace"/> the level q = <paramref name="quantile"/> in (0, 1), and
    /// <see cref="BassNoiseKind.GeneralizedHyperbolic"/> v's GIG <paramref name="mixingPrior"/>, a proper one.
    /// <paramref name="varianceScale"/> is c, positive and finite: under t, 1 makes σ the t's scale, BARS's <c>StudentTWeights</c>' t, and
    /// (ν − 2)/ν makes σ² the noise variance, GBASS's <c>tbass</c> for ν &gt; 2; the asymmetric Laplace's is 2/(q(1 − q)), which its q
    /// fixes. <paramref name="reference"/> is what the priors read of the response; <see cref="BassNoiseReference.Quantile"/> is the
    /// asymmetric Laplace's alone. <paramref name="skewPrior"/> is β's prior for t, the Normal-Wald and the generalized hyperbolic: null
    /// holds β at 0, and t takes one only for a skew (the skewed t). <paramref name="gammaPrior"/> is the Normal-Wald's γ's, held at 1 when
    /// null as GBASS's default holds it; a held γ is positive. <paramref name="betaPrimePrior"/> is the generalized beta prime's
    /// v ~ GBP(p, a, b), each positive and finite.
    /// </summary>
    public BassNoise(BassNoiseKind kind, double? degreesOfFreedom = null, double varianceScale = 1.0,
                     BassNoiseReference reference = BassNoiseReference.Moments, double? quantile = null,
                     BassNoisePrior? skewPrior = null, BassNoisePrior? gammaPrior = null, GigParameters? mixingPrior = null,
                     GbpParameters? betaPrimePrior = null)
    {
        if (kind != BassNoiseKind.StudentT && degreesOfFreedom is not null) throw new ArgumentException($"{kind} takes no ν.", nameof(degreesOfFreedom));
        if (kind != BassNoiseKind.AsymmetricLaplace && quantile is not null) throw new ArgumentException($"{kind} takes no quantile.", nameof(quantile));
        if (kind != BassNoiseKind.AsymmetricLaplace && reference == BassNoiseReference.Quantile)
            throw new ArgumentException("The quantile read is the asymmetric Laplace's.", nameof(reference));
        if (kind != BassNoiseKind.NormalWald && gammaPrior is not null) throw new ArgumentException("γ is the Normal-Wald's.", nameof(gammaPrior));
        if (kind != BassNoiseKind.GeneralizedHyperbolic && mixingPrior is not null)
            throw new ArgumentException("A declared GIG is the generalized hyperbolic's.", nameof(mixingPrior));
        if (kind != BassNoiseKind.GeneralizedBetaPrime && betaPrimePrior is not null)
            throw new ArgumentException("A declared GBP is the generalized beta prime's.", nameof(betaPrimePrior));
        switch (kind)
        {
            case BassNoiseKind.StudentT:
                if (degreesOfFreedom is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "Student t takes ν, positive and finite.");
                if (skewPrior is { IsHeld: true, Mean: 0.0 }) throw new ArgumentException("A skew held at 0 is the symmetric t: leave it out.", nameof(skewPrior));
                break;
            case BassNoiseKind.AsymmetricLaplace:
                if (quantile is not (> 0.0 and < 1.0)) throw new ArgumentOutOfRangeException(nameof(quantile), "The asymmetric Laplace takes q in (0, 1).");
                if (varianceScale != LaplaceScale(quantile.Value))
                    throw new ArgumentException($"The asymmetric Laplace's c is 2/(q(1 − q)) = {LaplaceScale(quantile.Value)} at q = {quantile.Value}.", nameof(varianceScale));
                if (skewPrior is not null) throw new ArgumentException("The asymmetric Laplace's skew is its level's.", nameof(skewPrior));
                break;
            case BassNoiseKind.NormalWald:
                gammaPrior ??= BassNoisePrior.Held(1.0);
                if (gammaPrior.IsHeld && !(gammaPrior.Mean > 0.0)) throw new ArgumentOutOfRangeException(nameof(gammaPrior), "A held γ is positive.");
                break;
            case BassNoiseKind.GeneralizedHyperbolic:
                if (mixingPrior is not { } gig) throw new ArgumentException("The generalized hyperbolic takes v's GIG.", nameof(mixingPrior));
                GeneralizedInverseGaussian.Of(gig.P, gig.A, gig.B);
                break;
            case BassNoiseKind.GeneralizedBetaPrime:
                if (betaPrimePrior is not { } gbp) throw new ArgumentException("The generalized beta prime takes v's GBP.", nameof(betaPrimePrior));
                gbp.ToDistribution();
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(kind));
        }
        if (!(varianceScale > 0.0 && double.IsFinite(varianceScale))) throw new ArgumentOutOfRangeException(nameof(varianceScale), "c must be positive and finite.");
        if (reference is not (BassNoiseReference.Moments or BassNoiseReference.Robust or BassNoiseReference.None or BassNoiseReference.Quantile))
            throw new ArgumentOutOfRangeException(nameof(reference));
        Kind = kind;
        DegreesOfFreedom = degreesOfFreedom;
        VarianceScale = varianceScale;
        Reference = reference;
        Quantile = quantile;
        SkewPrior = skewPrior;
        GammaPrior = gammaPrior;
        MixingPrior = mixingPrior;
        BetaPrimePrior = betaPrimePrior;
    }

    /// <summary>Student t on <paramref name="degreesOfFreedom"/> ν, with c = <paramref name="varianceScale"/>.</summary>
    public static BassNoise StudentT(double degreesOfFreedom, double varianceScale = 1.0, BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.StudentT, degreesOfFreedom, varianceScale, reference);

    /// <summary>The asymmetric Laplace at level <paramref name="quantile"/>, its priors reading the sample quantile by default.</summary>
    public static BassNoise AsymmetricLaplace(double quantile, BassNoiseReference reference = BassNoiseReference.Quantile)
        => new(BassNoiseKind.AsymmetricLaplace, varianceScale: LaplaceScale(quantile), reference: reference, quantile: quantile);

    /// <summary>The skewed t on <paramref name="degreesOfFreedom"/> ν, β under <paramref name="skew"/>, with c = <paramref name="varianceScale"/>.</summary>
    public static BassNoise SkewedT(double degreesOfFreedom, BassNoisePrior skew, double varianceScale = 1.0,
                                    BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.StudentT, degreesOfFreedom, varianceScale, reference, skewPrior: skew);

    /// <summary>The Normal-Wald, γ under <paramref name="gamma"/> (held at 1 when null) and β under <paramref name="skew"/> (held at 0 when
    /// null), with c = <paramref name="varianceScale"/>, GBASS's 1.</summary>
    public static BassNoise NormalWald(BassNoisePrior? gamma = null, BassNoisePrior? skew = null, double varianceScale = 1.0,
                                       BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.NormalWald, varianceScale: varianceScale, reference: reference, skewPrior: skew, gammaPrior: gamma);

    /// <summary>The generalized hyperbolic, v ~ <paramref name="mixing"/>, β under <paramref name="skew"/> (held at 0 when null).</summary>
    public static BassNoise GeneralizedHyperbolic(GigParameters mixing, BassNoisePrior? skew = null, double varianceScale = 1.0,
                                                  BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.GeneralizedHyperbolic, varianceScale: varianceScale, reference: reference, skewPrior: skew, mixingPrior: mixing);

    /// <summary>The generalized beta prime mixture, v ~ <paramref name="mixing"/>, β under <paramref name="skew"/> (held at 0 when null).</summary>
    public static BassNoise GeneralizedBetaPrime(GbpParameters mixing, BassNoisePrior? skew = null, double varianceScale = 1.0,
                                                 BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.GeneralizedBetaPrime, varianceScale: varianceScale, reference: reference, skewPrior: skew, betaPrimePrior: mixing);

    /// <summary>The horseshoe, v ~ GBP(1, ½, ½), β under <paramref name="skew"/> (held at 0 when null). GBASS's <c>hbass</c> is c = 1 here
    /// with σ² under <see cref="Sigma2Prior.HalfCauchy"/>, σ ~ C⁺(0, s) at the reference's mean square s², where <c>hbass</c> takes
    /// c = var(y) and √w ~ C⁺(0, 1).</summary>
    public static BassNoise Horseshoe(BassNoisePrior? skew = null, double varianceScale = 1.0, BassNoiseReference reference = BassNoiseReference.Moments)
        => GeneralizedBetaPrime(GbpParameters.Horseshoe, skew, varianceScale, reference);

    /// <summary>The Strawderman–Berger, v ~ GBP(½, 1, ½), β under <paramref name="skew"/> (held at 0 when null).</summary>
    public static BassNoise StrawdermanBerger(BassNoisePrior? skew = null, double varianceScale = 1.0,
                                              BassNoiseReference reference = BassNoiseReference.Moments)
        => GeneralizedBetaPrime(GbpParameters.StrawdermanBerger, skew, varianceScale, reference);

    // 2/(q(1 − q)), the asymmetric Laplace's c; q outside (0, 1) is the constructor's to refuse.
    private static double LaplaceScale(double q) => q > 0.0 && q < 1.0 ? 2.0 / (q * (1.0 - q)) : double.NaN;

    /// <summary>Which noise this is.</summary>
    public BassNoiseKind Kind { get; }

    /// <summary>Student t's ν; null for a kind that takes none.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? DegreesOfFreedom { get; }

    /// <summary>c, the variance of the Gaussian part per unit σ²·v.</summary>
    public double VarianceScale { get; }

    /// <summary>What the priors read of the response.</summary>
    public BassNoiseReference Reference { get; }

    /// <summary>
    /// Whether the sweep interweaves the noise's scale (E1's slice 7; <see cref="BassInterweavingBlock"/>): after the noise's blocks, σ is drawn
    /// again under the z parameterization, zᵢ = σ·vᵢ held, so σ and v move together and the offset σβ·v stands. Off by default, so a fit's
    /// draws stand as they were; left out of the serialized form when off. Refused with a held σ², which has no scale to move.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
    public bool Interweave { get; init; }

    /// <summary>
    /// Whether the latent slot adds the skew's rescaling move (<see cref="BassSkewRescaleBlock"/>): a Metropolis move of β and every vᵢ along
    /// the ridge β·vᵢ holds, β → β/κ and v → κ·v, the weights refitted once. Null takes the noise's default
    /// (<see cref="SkewRescaleByDefault"/>), on under the skewed horseshoe alone; true or false declares it. Left out of the serialized form
    /// when null, so a config that declares nothing serializes and fingerprints as it did. True needs a drawn skew, which it moves.
    /// </summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public bool? SkewRescale
    {
        get;
        init => field = value == true && !SkewDrawn
            ? throw new ArgumentException("The skew's rescaling moves β, so it needs a drawn skew.", nameof(SkewRescale))
            : value;
    }

    /// <summary>
    /// The move's default: on under the skewed horseshoe, v ~ GBP(1, ½, ½) with β drawn, off under every other law (the user's decision,
    /// 2026-11-08). There the sweep without it fails Geweke's test at b = 1 and the move passes; under the generalized hyperbolic, the
    /// Normal-Wald and the skewed t the sweep passes without it, so its 1.2 to 1.6 times a run's time buys nothing.
    /// </summary>
    [JsonIgnore]
    public bool SkewRescaleByDefault => Kind == BassNoiseKind.GeneralizedBetaPrime && BetaPrimePrior == GbpParameters.Horseshoe && SkewDrawn;

    /// <summary>Whether the sweep runs the move: <see cref="SkewRescale"/> where declared, the noise's default where not.</summary>
    [JsonIgnore]
    public bool RescalesSkew => SkewRescale ?? SkewRescaleByDefault;

    /// <summary>The asymmetric Laplace's q; null for a kind that takes none.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? Quantile { get; }

    /// <summary>β's prior; null where β is the kind's: 0, or the asymmetric Laplace's level's.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public BassNoisePrior? SkewPrior { get; }

    /// <summary>The Normal-Wald's γ's prior; null for another kind.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public BassNoisePrior? GammaPrior { get; }

    /// <summary>The generalized hyperbolic's v ~ GIG(p, a, b); null for another kind.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public GigParameters? MixingPrior { get; }

    /// <summary>The generalized beta prime's v ~ GBP(p, a, b); null for another kind.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public GbpParameters? BetaPrimePrior { get; }

    /// <summary>β, the drift per unit v the kind fixes or a chain starts at: 1/q − 1/(1 − q) for the asymmetric Laplace, the prior's mean
    /// where β has one, 0 otherwise. Read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public double Skew => Kind == BassNoiseKind.AsymmetricLaplace ? 1.0 / Quantile!.Value - 1.0 / (1.0 - Quantile.Value) : SkewPrior?.Mean ?? 0.0;

    /// <summary>Whether β is drawn.</summary>
    [JsonIgnore]
    public bool SkewDrawn => SkewPrior is { IsHeld: false };

    /// <summary>Whether the noise can be skewed, so σ² sits in the offset and is drawn apart from the structure: a skew drawn, or held away
    /// from 0.</summary>
    [JsonIgnore]
    public bool Skewed => Skew != 0.0 || SkewDrawn;

    /// <summary>The Normal-Wald's γ a chain starts at: the held value, or the mean of the prior's part on γ &gt; 0,
    /// m + s·φ(m/s)/Φ(m/s); null for another kind.</summary>
    [JsonIgnore]
    public double? GammaStart => GammaPrior switch
    {
        null => null,
        { IsHeld: true } held => held.Mean,
        { } prior => prior.Mean + prior.Sd * Math.Exp(-0.5 * Math.Pow(prior.Mean / prior.Sd, 2.0) - 0.5 * Math.Log(2.0 * Math.PI)
                                                     - Normal.StandardLogCdf(prior.Mean / prior.Sd)),
    };

    /// <summary>The kind's mixing law at the start, read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public IMixingLaw Mixing => MixingAt(null);

    /// <summary>The mixing law at γ = <paramref name="gamma"/> for the Normal-Wald (its start when null); the kind's own otherwise.</summary>
    internal IMixingLaw MixingAt(double? gamma) => Kind switch
    {
        BassNoiseKind.StudentT => new StudentTMixing(DegreesOfFreedom!.Value),
        BassNoiseKind.AsymmetricLaplace => ExponentialMixing.Instance,
        BassNoiseKind.NormalWald => GigMixing.NormalWald(gamma ?? GammaStart!.Value),
        BassNoiseKind.GeneralizedHyperbolic => new GigMixing(MixingPrior!.Value),
        BassNoiseKind.GeneralizedBetaPrime => new GbpMixing(BetaPrimePrior!.Value),
        _ => throw new InvalidOperationException($"No mixing law for {Kind}."),
    };
}

/// <summary>
/// An observation model's latents at one state, carried by the <see cref="WorkingResponse"/> they make (E1's decision 1): one value per
/// observation, a variance-mean mixture's scale v (a probit's z or a Pólya–Gamma ω would ride the same slot); the skew β, 0 for a
/// symmetric law; σ² when it is drawn apart from the structure, as under a skew, whose offset carries it (GBASS's w; null when the
/// structure step integrates it); and the Normal-Wald's tail γ, which sets v's law (null for a law without one). The target reads that σ²
/// from the state, as it reads its statistics. A retained draw keeps these alone and rebuilds its working response from them
/// (<see cref="IBassObservationModel.Rebuild"/>). Immutable; the values are copied.
/// </summary>
public sealed class ObservationLatents
{
    private readonly double[] _values;

    /// <summary>The latents <paramref name="values"/>, one per observation and each finite, with skew <paramref name="skew"/>, σ² =
    /// <paramref name="sigma2"/> drawn apart, and the Normal-Wald's <paramref name="gamma"/>, each positive and finite.</summary>
    public ObservationLatents(ReadOnlySpan<double> values, double skew = 0.0, double? sigma2 = null, double? gamma = null)
        : this(values.ToArray(), skew, sigma2, gamma) { }

    // Takes `values` as its own.
    internal ObservationLatents(double[] values, double skew, double? sigma2, double? gamma = null)
    {
        foreach (double v in values)
            if (!double.IsFinite(v)) throw new ArgumentException("A latent is finite.", nameof(values));
        if (!double.IsFinite(skew)) throw new ArgumentOutOfRangeException(nameof(skew), "The skew is finite.");
        if (sigma2 is double s && !(s > 0.0 && double.IsFinite(s))) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive and finite.");
        if (gamma is double g && !(g > 0.0 && double.IsFinite(g))) throw new ArgumentOutOfRangeException(nameof(gamma), "γ must be positive and finite.");
        _values = values;
        Skew = skew;
        Sigma2 = sigma2;
        Gamma = gamma;
    }

    /// <summary>The Normal-Wald's γ, v ~ GIG(−½, γ², 1); null for a law without one.</summary>
    public double? Gamma { get; }

    /// <summary>σ² drawn apart from the structure, the state's own; null when the structure step integrates it.</summary>
    public double? Sigma2 { get; }

    /// <summary>One latent per observation, case by case and location by location.</summary>
    public ReadOnlySpan<double> Values => _values;

    /// <summary>β, the drift per unit v; 0 for a symmetric law.</summary>
    public double Skew { get; }

    /// <summary>These latents at σ² = <paramref name="sigma2"/>, the values shared.</summary>
    internal ObservationLatents WithSigma2(double sigma2) => new(_values, Skew, sigma2, Gamma);

    /// <summary>These latents at skew <paramref name="skew"/> and γ = <paramref name="gamma"/>, the values shared.</summary>
    internal ObservationLatents WithShape(double skew, double? gamma) => new(_values, skew, Sigma2, gamma);

    /// <summary>The values' own array, not to be written.</summary>
    internal double[] ValueArray => _values;
}

/// <summary>
/// The noise's law at a retained draw, each observation's latent integrated out, which the readers of a fit take in place of the
/// Gaussian's (E1): the pointwise density (OQ17's scope rule), the predictive noise, the mean's shift and the variance. A model that
/// states none is Gaussian, N(0, σ²/w), read by the readers' own code as before.
/// </summary>
public interface IBassNoiseLaw
{
    /// <summary>ln p(ε = <paramref name="residual"/> | σ², the draw's latents' globals), v integrated out.</summary>
    double LogDensity(double residual, double sigma2, ObservationLatents? latents);

    /// <summary>The noise's distribution at a draw, v integrated out, for the predictive mixture's quantiles.</summary>
    IContinuousDistribution At(double sigma2, ObservationLatents? latents);

    /// <summary>Whether <see cref="LogDensity"/> is an integral over each observation's latent, hundreds of nodes a read, so that the
    /// pointwise density, N × draws reads, is opt-in.</summary>
    bool DensityByQuadrature => false;

    /// <summary>E[ε], σβE[v]: the mean of y at x sits this far from the surface f.</summary>
    double MeanShift(double sigma2, ObservationLatents? latents);

    /// <summary>Var[ε], σ²(β²Var[v] + c·E[v]); +∞ where it is infinite.</summary>
    double Variance(double sigma2, ObservationLatents? latents);

    /// <summary>A new observation's noise: v from its prior, then ζ, from <paramref name="rng"/> in that order.</summary>
    double Sample(double sigma2, ObservationLatents? latents, Xoshiro256PlusPlus rng);
}

/// <summary>
/// The normal variance-mean mixture observation model (E1): what the Gaussian core fits at the current latents, the blocks that move
/// them, and the noise's law. The working response is y − σβ·v under precision weights 1/(c·v), the latents riding it
/// (<see cref="ObservationLatents"/>), the prior's Gram unweighted (D03's robust row) and the priors reading the
/// <see cref="BassNoise.Reference"/> frozen at the start. A chain starts at v = 1. The sweep (D04's robust row): the structure with the
/// coefficients integrated out, τ, then σ² and a, or a given σ², from their conditional (<see cref="BassCoefficientsBlock"/>), then each vᵢ
/// by its law's step (<see cref="BassLatentScaleBlock"/>), then σ² (<see cref="BassScaleBlock"/>), then a skew and the Normal-Wald's γ where they
/// are drawn (<see cref="BassNoiseShapeBlock"/>). Under a pure scale mixture (t, or a skew held at 0) the structure step integrates σ² too
/// and redraws it with a each sweep, so σ² has no block of its own. Under a skew, held away from 0 or drawn, σ² sits in the offset, so it
/// is drawn apart and held on the record, where the target reads it; its draw is the modified half-normal on σ⁻¹, and a held σ² is the
/// config's own. v's law is read at the record's γ. On a tempered rung and an SMC stage the same blocks run at the rung's power b, the
/// likelihood given the latents tempered (<see cref="BassCoefficientsBlock.Power"/>), and the prior's b = 0 draw takes the latents from
/// theirs (<see cref="PriorLatents"/>). One model serves every chain: it holds no chain's state.
/// </summary>
public sealed class VarianceMeanMixture : IBassObservationModel, IBassNoiseLaw
{
    private readonly IBlock<ReversibleJumpChain<BassState>>[] _coefficients, _latent, _scale;
    private readonly double? _held;   // σ² as the config holds it, under a skew

    /// <summary>The model of <paramref name="noise"/> under σ²'s prior <paramref name="sigma2"/> (D03's default when null), which a
    /// skewed noise draws σ² against, or holds it at.</summary>
    public VarianceMeanMixture(BassNoise noise, Sigma2Prior? sigma2 = null)
    {
        ArgumentNullException.ThrowIfNull(noise);
        Noise = noise;
        Mixing = noise.Mixing;
        Skew = noise.Skew;
        Skewed = noise.Skewed;
        _held = sigma2?.IsHeld == true ? sigma2.Value : null;
        _coefficients = [new BassCoefficientsBlock()];
        _latent = noise.RescalesSkew ? [new BassLatentScaleBlock(this), new BassSkewRescaleBlock(this)] : [new BassLatentScaleBlock(this)];
        var scale = new List<IBlock<ReversibleJumpChain<BassState>>>();
        if (Skewed && _held is null) scale.Add(new BassScaleBlock(this));
        if (noise.SkewDrawn || noise.GammaPrior is { IsHeld: false }) scale.Add(new BassNoiseShapeBlock(this));
        if (noise.Interweave && _held is null) scale.Add(new BassInterweavingBlock(this));
        _scale = [.. scale];
    }

    /// <summary>The noise it models.</summary>
    public BassNoise Noise { get; }

    /// <summary>The mixing law of v at the start; the Normal-Wald's moves with γ (<see cref="MixingAt"/>).</summary>
    public IMixingLaw Mixing { get; }

    /// <summary>v's mixing law at <paramref name="latents"/>: the Normal-Wald's at their γ, the kind's own otherwise.</summary>
    public IMixingLaw MixingAt(ObservationLatents? latents)
        => Noise.Kind == BassNoiseKind.NormalWald && latents?.Gamma is double gamma ? Noise.MixingAt(gamma) : Mixing;

    /// <summary>c.</summary>
    public double VarianceScale => Noise.VarianceScale;

    /// <summary>β, the drift per unit v, the kind's or where a chain starts; 0 for a pure scale mixture.</summary>
    public double Skew { get; }

    /// <summary>Whether σ² is drawn apart and held on the record: a skew drawn, or held away from 0.</summary>
    public bool Skewed { get; }

    /// <summary>
    /// The working response at v = 1, β and γ at their starts: under a skew, σ² starts at the response's variance over the noise's per unit
    /// σ², β²Var[v] + c·E[v] (c where v's law has no finite moments), or at the held value; the reference is read from the response. A
    /// direct functional response takes a latent per observation, case by case and location by location, as a scalar one does.
    /// </summary>
    public WorkingResponse Working(BassResponse response)
    {
        ArgumentNullException.ThrowIfNull(response);
        int n = response.Observations;
        var ones = new double[n];
        Array.Fill(ones, 1.0);
        double? sigma2 = null, gamma = Noise.GammaStart;
        if (Skewed)
        {
            ReadOnlySpan<double> y = response.Values;
            double mean = 0.0, squares = 0.0;
            for (int i = 0; i < n; i++) mean += y[i];
            mean /= n;
            for (int i = 0; i < n; i++) squares += (y[i] - mean) * (y[i] - mean);
            IMixingLaw start = Noise.MixingAt(gamma);
            double perUnit = Skew * Skew * start.Variance + VarianceScale * start.Mean;
            if (!(perUnit > 0.0 && double.IsFinite(perUnit))) perUnit = VarianceScale;
            sigma2 = _held ?? Math.Max(squares / Math.Max(n - 1, 1), double.Epsilon) / perUnit;
        }
        BassPriorReference? reference = Noise.Reference switch
        {
            BassNoiseReference.Moments => BassPriorReference.Moments(response.Values),
            BassNoiseReference.Robust => BassPriorReference.Robust(response.Values),
            BassNoiseReference.Quantile => BassPriorReference.Quantile(response.Values, Noise.Quantile!.Value),
            _ => null,
        };
        var latents = new ObservationLatents(ones, Skew, sigma2, gamma);
        (double[]? offset, double[] weights) = Reading(latents);
        return new WorkingResponse(response, offset, weights, reference: reference, latents: latents);
    }

    /// <summary>
    /// Latents for <paramref name="observations"/> observations from their prior, the b = 0 member of the likelihood path that an exact
    /// reference rung and SMC's start draw: β from its prior where it is drawn (the kind's value otherwise), γ likewise, each vᵢ from its law
    /// at that γ, and under a skew σ² from <paramref name="sigma2Prior"/>, the prior at the state (a held σ² as the config holds it).
    /// Refused under a skew whose σ² prior is improper, which has no draw.
    /// </summary>
    internal ObservationLatents PriorLatents(int observations, InverseGamma? sigma2Prior, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double skew = Noise.SkewPrior is { IsHeld: false } b ? b.Mean + b.Sd * Normal.SampleStandard(rng) : Skew;
        double? gamma = Noise.GammaPrior switch
        {
            { IsHeld: false } g => Normal.MeanSd(g.Mean, g.Sd).SampleOn(rng, 0.0, double.PositiveInfinity),
            { } g => g.Mean,
            null => null,
        };
        IMixingLaw law = Noise.MixingAt(gamma);
        var v = new double[observations];
        for (int i = 0; i < observations; i++) v[i] = law.SamplePrior(rng);
        double? sigma2 = null;
        if (Skewed)
            sigma2 = _held ?? (sigma2Prior is { } prior
                ? prior.Sample(rng)
                : throw new InvalidOperationException("Under a skew σ² sits on the record, and an exact draw of it needs a proper prior: take an inverse gamma, a held σ² or the hierarchical prior."));
        return new ObservationLatents(v, skew, sigma2, gamma);
    }

    /// <summary>
    /// The working response at <paramref name="latents"/>, from <paramref name="from"/>, a reading of the same response: the offset σβ·v
    /// under a skew and the weights 1/(c·v) from the latents alone, so a draw's view rebuilt from them is the one the chain held, bit for
    /// bit. <paramref name="weightsMoved"/> false keeps <paramref name="from"/>'s weights, shared, for a move of σ² alone.
    /// </summary>
    internal WorkingResponse At(WorkingResponse from, ObservationLatents latents, bool weightsMoved = true)
    {
        (double[]? offset, double[] weights) = Reading(latents);
        return from.Moved(offset, weightsMoved ? weights : default, latents);
    }

    // The offset σβ·v (none at β = 0) and the weights 1/(c·v) at `latents`.
    private (double[]? Offset, double[] Weights) Reading(ObservationLatents latents)
    {
        ReadOnlySpan<double> v = latents.Values;
        double c = VarianceScale;
        var weights = new double[v.Length];
        for (int i = 0; i < v.Length; i++) weights[i] = 1.0 / (c * v[i]);
        if (latents.Skew == 0.0) return (null, weights);
        double drift = Math.Sqrt(latents.Sigma2 ?? throw new InvalidOperationException("A skewed noise's latents hold σ².")) * latents.Skew;
        var offset = new double[v.Length];
        for (int i = 0; i < v.Length; i++) offset[i] = drift * v[i];
        return (offset, weights);
    }

    /// <summary>True: its retained draws keep their latents alone.</summary>
    public bool Rebuilds => true;

    /// <summary>The working response a retained draw's latents read on <paramref name="model"/>, the model's working response: its
    /// rebuild, equal to the one the chain held.</summary>
    public WorkingResponse? Rebuild(WorkingResponse model, ObservationLatents latents)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(latents);
        return latents.Values.Length == model.Response.Observations ? At(model, latents) : null;
    }

    /// <summary>a's block in the coefficients' slot, v's in the latent slot; in the scale slot σ²'s under a skew with σ² not held, then the
    /// noise's shape where a skew or γ is drawn.</summary>
    public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => slot switch
    {
        BassSweepSlot.Coefficients => _coefficients,
        BassSweepSlot.Latent => _latent,
        BassSweepSlot.Scale => _scale,
        _ => [],
    };

    /// <summary>The noise's own law.</summary>
    public IBassNoiseLaw Law => this;

    /// <inheritdoc/>
    public double LogDensity(double residual, double sigma2, ObservationLatents? latents)
        => MixingAt(latents).MarginalLogDensity(residual, sigma2 * VarianceScale, Math.Sqrt(sigma2) * (latents?.Skew ?? 0.0));

    /// <summary>The noise at a draw: t's and the asymmetric Laplace's closed forms; under a skewed t, the Normal-Wald, the generalized
    /// hyperbolic and the generalized beta prime, whose CDF has no closed form, that distribution with its tails integrated over v
    /// (<see cref="GigMixing.Noise"/>, <see cref="GbpMixing.Noise"/>).</summary>
    public IContinuousDistribution At(double sigma2, ObservationLatents? latents)
    {
        double skew = latents?.Skew ?? 0.0;
        if (Noise.Kind == BassNoiseKind.StudentT && skew == 0.0)
            return StudentT.LocationScale(Noise.DegreesOfFreedom!.Value, 0.0, Math.Sqrt(sigma2 * VarianceScale));
        if (Noise.Kind == BassNoiseKind.AsymmetricLaplace) return AsymmetricLaplace.Of(0.0, Math.Sqrt(sigma2), Noise.Quantile!.Value);
        double unitVariance = sigma2 * VarianceScale, drift = Math.Sqrt(sigma2) * skew;
        return MixingAt(latents) switch
        {
            GigMixing mixing => mixing.Noise(unitVariance, drift),
            StudentTMixing t => new GigMixing(t.Prior).Noise(unitVariance, drift),
            GbpMixing mixing => mixing.Noise(unitVariance, drift),
            _ => throw new InvalidOperationException($"No distribution for {Noise.Kind}."),
        };
    }

    /// <summary>Whether the density is the mixing law's integral over v: the generalized beta prime's.</summary>
    public bool DensityByQuadrature => Mixing.MarginalByQuadrature;

    /// <inheritdoc/>
    public double MeanShift(double sigma2, ObservationLatents? latents)
    {
        double skew = latents?.Skew ?? 0.0;
        return skew == 0.0 ? 0.0 : Math.Sqrt(sigma2) * skew * MixingAt(latents).Mean;
    }

    /// <inheritdoc/>
    public double Variance(double sigma2, ObservationLatents? latents)
    {
        double skew = latents?.Skew ?? 0.0;
        IMixingLaw mixing = MixingAt(latents);
        return sigma2 * ((skew == 0.0 ? 0.0 : skew * skew * mixing.Variance) + VarianceScale * mixing.Mean);
    }

    /// <inheritdoc/>
    public double Sample(double sigma2, ObservationLatents? latents, Xoshiro256PlusPlus rng)
    {
        double v = MixingAt(latents).SamplePrior(rng);
        return Math.Sqrt(sigma2) * ((latents?.Skew ?? 0.0) * v + Math.Sqrt(VarianceScale * v) * Normal.SampleStandard(rng));
    }
}

/// <summary>
/// The coefficients' block (D04's robust row, the coefficients' slot): σ² and a from their conditional given the structure, τ and the
/// working response at the rung's power b, exact, at the cost of a retained draw's composition (<see cref="BassStructureMoves.Compose"/>),
/// with the fit f at every observation. They ride the chain's own moves to the blocks after it in the same sweep, which read the
/// residuals y − f; being redrawn each sweep they are not state. A σ² held or given is taken as it stands.
/// <para><b>On a rung.</b> The noise's blocks temper the likelihood given the latents, p(θ)·p(v)·N(y | f + σβv, σ²c·v)^b, the integrated
/// form on the likelihood path: a and σ² are integrated against the tempered likelihood (the structure step's ℓ_b), so their conditional
/// is the same block's at b, and β and each vᵢ keep theirs at b (<see cref="Power"/>). At b = 0 every block draws from its prior.</para>
/// </summary>
public sealed class BassCoefficientsBlock : IBlock<ReversibleJumpChain<BassState>>
{
    /// <inheritdoc/>
    public string Name => "coefficients";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        moves.DrawSweepFit(state.Current, rng, state.Target as BassTarget, Power(state));
        return true;
    }

    /// <summary>
    /// The likelihood's power b at <paramref name="state"/>'s rung, which every latent block reads. Refused off the likelihood path or in the
    /// linear form: the powered path and the linear form temper the marginal over a and σ², b·ℓ₁, whose weights 1/(c·vᵢ) couple every
    /// latent, so no vᵢ keeps a conditional of its own there.
    /// </summary>
    internal static double Power(ReversibleJumpChain<BassState> state)
    {
        Temperature temperature = state.Temperature;
        double b = temperature.InverseTemperature;
        if (b == 1.0) return 1.0;
        if (temperature.Path != TemperingPath.Likelihood || (state.Target as BassTarget)?.Form != TemperedForm.Integrated)
            throw new NotSupportedException("A latent noise's blocks temper the likelihood given its latents, on the likelihood path in the integrated form; the "
                                            + "linear form and the powered path temper the marginal over a and σ², which couples every latent.");
        return b;
    }
}

/// <summary>
/// The latent scales' block (the latent slot): each vᵢ by its law's step under its conditional given the residual rᵢ = yᵢ − fᵢ of this
/// sweep's coefficients and σ² (<see cref="BassCoefficientsBlock"/>, which runs first), under the drift σβ and the unit variance σ²c
/// (<see cref="IMixingLaw.StepConditional"/>) at the rung's power b: an exact draw where the law has one, else a step from the current vᵢ;
/// then the chain follows the working response at the new v (<see cref="BassStructureMoves.Follow"/>), its offset and weights moved
/// together. A draw that puts the current structure outside F is refused and the chain keeps its state, the exact step under the gate
/// that reads W (decision 4); the scan's counters read its refusals.
/// </summary>
public sealed class BassLatentScaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _model;

    /// <summary>The block of <paramref name="model"/>'s latents.</summary>
    public BassLatentScaleBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "latent scales";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassSweepFit fit = moves.SweepFit is { } drawn && ReferenceEquals(drawn.State, state.Current)
            ? drawn
            : throw new InvalidOperationException("The latent scales read this sweep's coefficients: the coefficients' block runs first, at the current state.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        ReadOnlySpan<double> y = working.Response.Values, current = latents.Values;
        double unitVariance = fit.Sigma2 * _model.VarianceScale, drift = Math.Sqrt(fit.Sigma2) * latents.Skew;
        IMixingLaw mixing = _model.MixingAt(latents);
        int n = y.Length;
        var v = new double[n];
        double power = BassCoefficientsBlock.Power(state);
        for (int i = 0; i < n; i++) v[i] = mixing.StepConditional(current[i], y[i] - fit.Fitted[i], unitVariance, drift, power, rng);
        WorkingResponse moved = _model.At(working, new ObservationLatents(v, latents.Skew, latents.Sigma2, latents.Gamma));
        return BassStructureMoves.Follow(state, moved) is null;
    }
}

/// <summary>
/// The skew's rescaling move (the latent slot, after v's block, under a drawn skew): a Metropolis move of β and every vᵢ along the ridge
/// β·vᵢ holds, (β, v) → (β/κ, κ·v), slice 5's and slice 7's follow-up. The offset σβ·v stands and the weights 1/(c·vᵢ) move to W/κ, refitted
/// once. On the sweep's joint with the coefficients integrated, π ∝ p(structure, τ)·p(σ²)·p(β)·Πp_v(vᵢ)·ℓ_b, ℓ_b the target's marginal given σ²
/// on the working response at the rung's power (the weights' normalizer in it), the map has Jacobian κ^(N−1) and its reverse is 1/κ, so with
/// t = ln κ the move accepts with
/// <para>ln A = Δℓ_b + Σᵢ [ln p_v(κ·vᵢ) − ln p_v(vᵢ)] + ln p(β/κ) − ln p(β) + (N − 1)·t + ln q(−t | x′) − ln q(t | x),</para>
/// <para>declining a move that leaves F, as π is 0 there. It is a step on the coefficients' marginal, so after an accepted move they are drawn
/// afresh at the new state before a block reads them (van Dyk and Park 2008), and its proposal reads the state alone: t ~ N(0, h(x)²) with
/// h = 2.4/√C(x), C the log target's curvature along t, the likelihood's b·N/2 less the second derivatives of v's prior and β's along t, at
/// least 1. The scale moves with the state, so the reverse's density is read at x′ (Hastings). A form holding the coefficients, accepting on
/// the likelihood given them at O(N), was measured beside it: it mixed v's mean slower per second under the skewed horseshoe and failed
/// its Geweke there at b = 1 (the brief's report).</para>
/// </summary>
public sealed class BassSkewRescaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private const double Step = 2.4;       // the random walk's scale over the target's sd along t
    private const double Difference = 1e-2;   // the finite difference v's prior's curvature is read at
    private readonly VarianceMeanMixture _model;

    internal BassSkewRescaleBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (model.Noise.SkewPrior is not { IsHeld: false }) throw new ArgumentException("The skew's rescaling needs a drawn skew.", nameof(model));
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "skew rescale";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        Proposal x = Read(state);
        double t = x.Scale * Normal.SampleStandard(rng);
        if (Propose(state, x, t) is not { } move) return false;
        if (!(Math.Log(1.0 - rng.NextDouble()) < move.LogAcceptance)) return false;
        state.SetState(move.Candidate, move.Evaluation);
        // A step on the coefficients' marginal: they are drawn afresh at the new state, so the blocks after it condition on a draw of theirs.
        BassStructureMoves.Of(state).DrawSweepFit(state.Current, rng, state.Target as BassTarget, x.Power);
        return true;
    }

    /// <summary>The state's reading the move needs: v, β, σ², v's law, β's prior, the power b and the proposal's scale h(x).</summary>
    internal readonly record struct Proposal(double[] V, double Skew, double Sigma2, IMixingLaw Mixing, BassNoisePrior SkewPrior, double Power,
                                             double Scale, WorkingResponse Working, ObservationLatents Latents);

    /// <summary>A proposed move: the candidate state, its evaluation, and ln A.</summary>
    internal readonly record struct Move(BassState Candidate, ITargetEvaluation Evaluation, double LogAcceptance);

    // The current state's reading, which reads no coefficients.
    internal Proposal Read(ReversibleJumpChain<BassState> state)
    {
        WorkingResponse working = BassStructureMoves.Of(state).DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        double sigma2 = latents.Sigma2 ?? throw new InvalidOperationException("A drawn skew's latents hold σ².");
        IMixingLaw mixing = _model.MixingAt(latents);
        BassNoisePrior prior = _model.Noise.SkewPrior!;
        double b = BassCoefficientsBlock.Power(state);
        double[] values = latents.Values.ToArray();
        return new Proposal(values, latents.Skew, sigma2, mixing, prior, b, Scale(values, latents.Skew, 0.5 * b * values.Length, mixing, prior),
                            working, latents);
    }

    // h = Step/√C, C = the likelihood's curvature − ∂²/∂t² Σᵢ ln p_v(vᵢ·eᵗ) − ∂²/∂t² ln p(β·e⁻ᵗ) at t = 0, at least 1.
    internal static double Scale(double[] v, double beta, double likelihood, IMixingLaw mixing, BassNoisePrior prior)
    {
        double up = Math.Exp(Difference), down = Math.Exp(-Difference), centre = 0.0, above = 0.0, below = 0.0;
        foreach (double vi in v)
        {
            centre += mixing.LogPrior(vi);
            above += mixing.LogPrior(vi * up);
            below += mixing.LogPrior(vi * down);
        }
        double latent = (above - 2.0 * centre + below) / (Difference * Difference);
        double skew = -(2.0 * beta * beta - prior.Mean * beta) / (prior.Sd * prior.Sd);
        double curvature = likelihood - latent - skew;
        return Step / Math.Sqrt(double.IsFinite(curvature) ? Math.Max(curvature, 1.0) : 1.0);
    }

    // The move at t from x, or null when it leaves F: the candidate, its evaluation at the chain's temperature, and ln A.
    internal Move? Propose(ReversibleJumpChain<BassState> state, Proposal x, double t)
    {
        double kappa = Math.Exp(t), inverse = Math.Exp(-t);
        int n = x.V.Length;
        var moved = new double[n];
        double latentRatio = 0.0;
        for (int i = 0; i < n; i++)
        {
            moved[i] = kappa * x.V[i];
            latentRatio += x.Mixing.LogPrior(moved[i]) - x.Mixing.LogPrior(x.V[i]);
        }
        double skew = x.Skew * inverse;
        double sd = x.SkewPrior.Sd, m = x.SkewPrior.Mean;
        double skewRatio = -((skew - m) * (skew - m) - (x.Skew - m) * (x.Skew - m)) / (2.0 * sd * sd);
        double reverse = Scale(moved, skew, 0.5 * x.Power * n, x.Mixing, x.SkewPrior);
        double proposalRatio = Math.Log(x.Scale) - Math.Log(reverse) - t * t / (2.0 * reverse * reverse) + t * t / (2.0 * x.Scale * x.Scale);
        WorkingResponse working = _model.At(x.Working, new ObservationLatents(moved, skew, x.Sigma2, x.Latents.Gamma));
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassState candidate = moves.ProposeFollow(state.Current, working);
        if ((state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.")).Refusal(candidate.Fit) is not null) return null;
        ITargetEvaluation evaluation = state.Target.Evaluate(candidate);
        double target = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity;
        return new Move(candidate, evaluation, target + latentRatio + skewRatio + (n - 1) * t + proposalRatio);
    }
}

/// <summary>
/// σ²'s block under a skew (the scale slot; E1's decision 2): σ² from its conditional given the coefficients of this sweep and the latents.
/// With rᵢ = yᵢ − fᵢ, the coefficients' prior variance scaling with σ², K of them inside the prior with energy Q = τ(a − a₀)ᵀD(a − a₀)
/// (<see cref="GPriorConditional.PriorEnergy"/>), and σ²'s prior IG(α₀, β₀) at the state (the hierarchical prior's at its rate; Jeffreys
/// at α₀ = β₀ = 0), u = σ⁻¹ is the modified half-normal MHN(bN + K + 2α₀, ½(bΣᵢ rᵢ²/(c·vᵢ) + Q + 2β₀), bβΣᵢ rᵢ/c) at the rung's power b
/// (Sun, Kong and Pal 2021; GBASS's w step with w's prior a kept, which GBASS drops, E1's slice 0): the coefficients' prior is untempered.
/// The chain follows the offset at the new σ², the weights shared; a draw that puts the structure outside F is refused as the latent
/// block's is.
/// </summary>
public sealed class BassScaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _model;

    /// <summary>The block of <paramref name="model"/>'s σ².</summary>
    public BassScaleBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "scale";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassSweepFit fit = moves.SweepFit is { } drawn && ReferenceEquals(drawn.Structure, state.Current.TermArray)
            ? drawn
            : throw new InvalidOperationException("σ² reads this sweep's coefficients: the coefficients' block runs first, at the current structure.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        ReadOnlySpan<double> y = working.Response.Values, weights = working.Weights;
        double c = _model.VarianceScale, squares = 0.0, sum = 0.0;
        for (int i = 0; i < y.Length; i++)
        {
            double r = y[i] - fit.Fitted[i];
            squares += r * r * weights[i];
            sum += r;
        }
        (double quadratic, int count) = fit.Conditional.PriorEnergy(fit.Centred);
        InverseGamma? prior = (state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.")).Sigma2At(state.Current);
        double u = VarianceMeanConditionals.SampleInverseScale(y.Length, squares, sum, quadratic, count, prior, latents.Skew, c,
                                                               BassCoefficientsBlock.Power(state), rng);
        WorkingResponse moved = _model.At(working, latents.WithSigma2(1.0 / (u * u)), weightsMoved: false);
        return BassStructureMoves.Follow(state, moved) is null;
    }
}

/// <summary>
/// The noise's shape (the scale slot, after σ²): the skew β and the Normal-Wald's γ, each where it is drawn, from its conditional; GBASS's
/// <c>nwbass</c> steps, γ's on its support. With rᵢ = yᵢ − fᵢ of this sweep's coefficients, σ² and c, rᵢ/σ = βvᵢ + √(c·vᵢ)ζᵢ makes β | ·
/// normal with precision bΣvᵢ/c + 1/s² and mean (bΣrᵢ/(σc) + m/s²) over it, under β ~ N(m, s²), at the rung's power b; γ reads v's prior
/// alone, untempered. v ~ GIG(−½, γ², 1), whose normalizer is
/// e^γ/√(2π), makes Πp(vᵢ | γ) ∝ exp(Nγ − γ²Σvᵢ/2), so γ | v is normal with precision Σvᵢ + 1/s² and mean (N + m/s²) over it, on γ &gt; 0,
/// drawn there (<see cref="Normal.SampleOn"/>): the correction of GBASS's draw, which takes that normal on the whole line (E1's slice 0).
/// The target reads neither, and neither prior enters it: no move of the structure changes them. The chain follows the offset at the new
/// β, the weights shared; a draw that puts the structure outside F is refused as the latent block's is.
/// </summary>
public sealed class BassNoiseShapeBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _model;

    /// <summary>The block of <paramref name="model"/>'s skew and γ.</summary>
    public BassNoiseShapeBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "noise shape";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassSweepFit fit = moves.SweepFit is { } drawn && ReferenceEquals(drawn.Structure, state.Current.TermArray)
            ? drawn
            : throw new InvalidOperationException("The noise's shape reads this sweep's coefficients: the coefficients' block runs first, at the current structure.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        ReadOnlySpan<double> y = working.Response.Values, v = latents.Values;
        double sumV = 0.0;
        foreach (double vi in v) sumV += vi;
        double skew = latents.Skew;
        if (_model.Noise.SkewPrior is { IsHeld: false } skewPrior)
        {
            double sigma = Math.Sqrt(latents.Sigma2 ?? throw new InvalidOperationException("A drawn skew's latents hold σ².")), c = _model.VarianceScale;
            double sumR = 0.0;
            for (int i = 0; i < y.Length; i++) sumR += y[i] - fit.Fitted[i];
            skew = VarianceMeanConditionals.SampleSkew(sumV, sumR, sigma, c, skewPrior.Mean, skewPrior.Sd, BassCoefficientsBlock.Power(state), rng);
        }
        double? gamma = latents.Gamma;
        if (_model.Noise.GammaPrior is { IsHeld: false } gammaPrior)
        {
            double inverse = 1.0 / (gammaPrior.Sd * gammaPrior.Sd), precision = sumV + inverse;
            gamma = Normal.MeanSd((v.Length + gammaPrior.Mean * inverse) / precision, 1.0 / Math.Sqrt(precision))
                          .SampleOn(rng, 0.0, double.PositiveInfinity);
        }
        WorkingResponse moved = _model.At(working, latents.WithShape(skew, gamma), weightsMoved: false);
        return BassStructureMoves.Follow(state, moved) is null;
    }
}

/// <summary>
/// The interweaving step (E1's slice 7; the ancillarity-sufficiency interweaving of Yu and Meng 2011, on the noise's scale): σ drawn again
/// under the z parameterization, zᵢ = σ·vᵢ held, in which εᵢ = β·zᵢ + √(c·σ·zᵢ)·ξᵢ (the user's derivation, E1, 2026-10-26). v is the scale's
/// ancillary augmentation, its prior free of σ, and z the other, so a sweep that draws σ under each moves it along both.
/// <para>With u = σ, rᵢ = yᵢ − fᵢ at this sweep's coefficients and A = Σᵢ (rᵢ − β·zᵢ)²/(c·zᵢ), u's conditional given z multiplies the
/// likelihood at power b, u^(−bN/2)·exp(−b·A/(2u)), v's prior at zᵢ/u with its Jacobian, Πᵢ p_v(zᵢ/u)/u, σ²'s prior IG(α₀, β₀) carried to u,
/// and the coefficients' prior, whose variance scales with σ² (the equivariant prior), u^(−K)·exp(−Q/(2u²)). In t = ln u that is
/// −(bN/2 + N + 2α₀ + K)·t − (bA/2)·e^(−t) − (β₀ + Q/2)·e^(−2t) + Σᵢ ln p_v(zᵢ·e^(−t)): log-concave under GIG mixing, and drawn here by one
/// slice step for every law. A prior on u in the GIG family with the coefficients' variance scaling as u, or not at all, would make it a GIG
/// (the derivation's exact route); BASS's priors are equivariant, so the slice step is the route here.</para>
/// <para>Then vᵢ = zᵢ/u: the offset σβ·v = β·z stands and the weights 1/(c·vᵢ) move. Under a skew σ² goes on the record; under a pure scale
/// mixture σ² is the sweep's draw with the coefficients, so only v moves. A move that puts the structure outside F is refused and the chain
/// keeps its state, as the latent block's is.</para>
/// </summary>
public sealed class BassInterweavingBlock : IBlock<ReversibleJumpChain<BassState>>
{
    // The slice's width in ln σ and its stepping-out limit: σ's conditional given z is a fraction of a unit wide in ln σ.
    private const double Width = 1.0;
    private const int MaxSteps = 32;
    private readonly VarianceMeanMixture _model;

    internal BassInterweavingBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "interweaving";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassSweepFit fit = moves.SweepFit is { } drawn && ReferenceEquals(drawn.Structure, state.Current.TermArray)
            ? drawn
            : throw new InvalidOperationException("The interweaving step reads this sweep's coefficients: the coefficients' block runs first, at the current structure.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        ReadOnlySpan<double> y = working.Response.Values, v = latents.Values;
        int n = y.Length;
        double sigma = Math.Sqrt(latents.Sigma2 ?? fit.Sigma2), c = _model.VarianceScale, beta = latents.Skew, b = BassCoefficientsBlock.Power(state);
        var z = new double[n];
        double a = 0.0;
        for (int i = 0; i < n; i++)
        {
            z[i] = sigma * v[i];
            double e = y[i] - fit.Fitted[i] - beta * z[i];
            a += e * e / (c * z[i]);
        }
        (double quadratic, int count) = fit.Conditional.PriorEnergy(fit.Centred);
        InverseGamma? prior = (state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.")).Sigma2At(state.Current);
        double shape = prior?.Shape ?? 0.0, scale = prior?.Scale ?? 0.0;
        var density = new ScaleGivenZ(z, _model.MixingAt(latents), -(0.5 * b * n + n + 2.0 * shape + count), 0.5 * b * a, scale + 0.5 * quadratic);
        double t = Math.Log(sigma), logTarget = density.At(t);
        SliceSampler.Step(ref t, ref logTarget, Width, MaxSteps, density, rng);
        double u = Math.Exp(t);
        var moved = new double[n];
        for (int i = 0; i < n; i++) moved[i] = z[i] / u;
        var next = new ObservationLatents(moved, beta, latents.Sigma2 is null ? null : u * u, latents.Gamma);
        return BassStructureMoves.Follow(state, _model.At(working, next)) is null;
    }

    // σ's conditional given z in t = ln σ, up to a constant: slope·t − h·e^(−t) − q·e^(−2t) + Σᵢ ln p_v(zᵢ·e^(−t)).
    private readonly struct ScaleGivenZ(double[] z, IMixingLaw mixing, double slope, double h, double q) : IScalarLogDensity
    {
        public double At(double t)
        {
            double inverse = Math.Exp(-t), sum = 0.0;
            foreach (double zi in z) sum += mixing.LogPrior(zi * inverse);
            return slope * t - h * inverse - q * inverse * inverse + sum;
        }
    }
}

/// <summary>A sweep's draw of σ² and the coefficients at <paramref name="State"/>: the coefficients (α, β_c) in the design's centred order,
/// the conditional they came from, and the fit f at every observation. Not state: the blocks after the coefficients' read it within the
/// sweep. The structure is <paramref name="State"/>'s terms' array, which a later block's state shares while the structure stands.</summary>
internal sealed record BassSweepFit(BassState State, double Sigma2, double[] Centred, double[] Fitted, GPriorConditional Conditional)
{
    /// <summary>The structure the coefficients belong to.</summary>
    public BassTerm[] Structure => State.TermArray;
}
