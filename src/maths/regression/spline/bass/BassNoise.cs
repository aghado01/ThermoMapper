using System;
using System.Collections.Generic;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.Regression.GPrior;
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
    /// positive and finite, and <see cref="BassNoiseKind.AsymmetricLaplace"/> the level q = <paramref name="quantile"/> in (0, 1).
    /// <paramref name="varianceScale"/> is c, positive and finite: under t, 1 makes σ the t's scale, BARS's <c>StudentTWeights</c>' t, and
    /// (ν − 2)/ν makes σ² the noise variance, GBASS's <c>tbass</c> for ν &gt; 2; the asymmetric Laplace's is 2/(q(1 − q)), which its q
    /// fixes. <paramref name="reference"/> is what the priors read of the response; <see cref="BassNoiseReference.Quantile"/> is the
    /// asymmetric Laplace's alone.
    /// </summary>
    public BassNoise(BassNoiseKind kind, double? degreesOfFreedom = null, double varianceScale = 1.0,
                     BassNoiseReference reference = BassNoiseReference.Moments, double? quantile = null)
    {
        switch (kind)
        {
            case BassNoiseKind.StudentT:
                if (degreesOfFreedom is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "Student t takes ν, positive and finite.");
                if (quantile is not null) throw new ArgumentException("Student t takes no quantile.", nameof(quantile));
                if (reference == BassNoiseReference.Quantile) throw new ArgumentException("The quantile read is the asymmetric Laplace's.", nameof(reference));
                break;
            case BassNoiseKind.AsymmetricLaplace:
                if (quantile is not (> 0.0 and < 1.0)) throw new ArgumentOutOfRangeException(nameof(quantile), "The asymmetric Laplace takes q in (0, 1).");
                if (degreesOfFreedom is not null) throw new ArgumentException("The asymmetric Laplace takes no ν.", nameof(degreesOfFreedom));
                if (varianceScale != LaplaceScale(quantile.Value))
                    throw new ArgumentException($"The asymmetric Laplace's c is 2/(q(1 − q)) = {LaplaceScale(quantile.Value)} at q = {quantile.Value}.", nameof(varianceScale));
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
    }

    /// <summary>Student t on <paramref name="degreesOfFreedom"/> ν, with c = <paramref name="varianceScale"/>.</summary>
    public static BassNoise StudentT(double degreesOfFreedom, double varianceScale = 1.0, BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.StudentT, degreesOfFreedom, varianceScale, reference);

    /// <summary>The asymmetric Laplace at level <paramref name="quantile"/>, its priors reading the sample quantile by default.</summary>
    public static BassNoise AsymmetricLaplace(double quantile, BassNoiseReference reference = BassNoiseReference.Quantile)
        => new(BassNoiseKind.AsymmetricLaplace, varianceScale: LaplaceScale(quantile), reference: reference, quantile: quantile);

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

    /// <summary>The asymmetric Laplace's q; null for a kind that takes none.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? Quantile { get; }

    /// <summary>β, the drift per unit v the kind fixes: 1/q − 1/(1 − q) for the asymmetric Laplace, 0 for t. Read from the kind, so not
    /// serialized.</summary>
    [JsonIgnore]
    public double Skew => Kind == BassNoiseKind.AsymmetricLaplace ? 1.0 / Quantile!.Value - 1.0 / (1.0 - Quantile.Value) : 0.0;

    /// <summary>The kind's mixing law, read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public IMixingLaw Mixing => Kind switch
    {
        BassNoiseKind.StudentT => new StudentTMixing(DegreesOfFreedom!.Value),
        BassNoiseKind.AsymmetricLaplace => ExponentialMixing.Instance,
        _ => throw new InvalidOperationException($"No mixing law for {Kind}."),
    };
}

/// <summary>
/// An observation model's latents at one state, carried by the <see cref="WorkingResponse"/> they make (E1's decision 1): one value per
/// observation, a variance-mean mixture's scale v (a probit's z or a Pólya–Gamma ω would ride the same slot); the skew β, 0 for a
/// symmetric law; and σ² when it is drawn apart from the structure, as under a skew, whose offset carries it (GBASS's w; null when the
/// structure step integrates it). The target reads that σ² from the state, as it reads its statistics. A retained draw keeps these alone
/// and rebuilds its working response from them (<see cref="IBassObservationModel.Rebuild"/>). Immutable; the values are copied.
/// </summary>
public sealed class ObservationLatents
{
    private readonly double[] _values;

    /// <summary>The latents <paramref name="values"/>, one per observation and each finite, with skew <paramref name="skew"/> and, drawn
    /// apart, σ² = <paramref name="sigma2"/>, positive and finite.</summary>
    public ObservationLatents(ReadOnlySpan<double> values, double skew = 0.0, double? sigma2 = null)
        : this(values.ToArray(), skew, sigma2) { }

    // Takes `values` as its own.
    internal ObservationLatents(double[] values, double skew, double? sigma2)
    {
        foreach (double v in values)
            if (!double.IsFinite(v)) throw new ArgumentException("A latent is finite.", nameof(values));
        if (!double.IsFinite(skew)) throw new ArgumentOutOfRangeException(nameof(skew), "The skew is finite.");
        if (sigma2 is double s && !(s > 0.0 && double.IsFinite(s))) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive and finite.");
        _values = values;
        Skew = skew;
        Sigma2 = sigma2;
    }

    /// <summary>σ² drawn apart from the structure, the state's own; null when the structure step integrates it.</summary>
    public double? Sigma2 { get; }

    /// <summary>One latent per observation, case by case and location by location.</summary>
    public ReadOnlySpan<double> Values => _values;

    /// <summary>β, the drift per unit v; 0 for a symmetric law.</summary>
    public double Skew { get; }

    /// <summary>These latents at σ² = <paramref name="sigma2"/>, the values shared.</summary>
    internal ObservationLatents WithSigma2(double sigma2) => new(_values, Skew, sigma2);

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
/// from its own (<see cref="BassLatentScaleBlock"/>), then σ² (<see cref="BassScaleBlock"/>). Under a pure scale mixture (t) the structure
/// step integrates σ² too and redraws it with a each sweep, so the scale slot is empty. Under a skew (the asymmetric Laplace) σ² sits in
/// the offset, so it is drawn apart and held on the record, where the target reads it; its draw is the modified half-normal on σ⁻¹, and
/// a held σ² is the config's own. One model serves every chain: it holds no chain's state.
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
        _held = sigma2?.IsHeld == true ? sigma2.Value : null;
        _coefficients = [new BassCoefficientsBlock()];
        _latent = [new BassLatentScaleBlock(this)];
        _scale = Skew != 0.0 && _held is null ? [new BassScaleBlock(this)] : [];
    }

    /// <summary>The noise it models.</summary>
    public BassNoise Noise { get; }

    /// <summary>The mixing law of v.</summary>
    public IMixingLaw Mixing { get; }

    /// <summary>c.</summary>
    public double VarianceScale => Noise.VarianceScale;

    /// <summary>β, the drift per unit v, the kind's; 0 for a pure scale mixture.</summary>
    public double Skew { get; }

    /// <summary>
    /// The working response at v = 1: under a skew, σ² starts at the response's variance over the noise's per unit σ², β²Var[v] +
    /// c·E[v], or at the held value; the reference is read from the response.
    /// </summary>
    public WorkingResponse Working(BassResponse response)
    {
        ArgumentNullException.ThrowIfNull(response);
        if (response is FunctionalResponse)
            throw new NotSupportedException("A per-observation latent under a direct functional design is E1's slice 5; a scalar response takes it now.");
        int n = response.Observations;
        var ones = new double[n];
        Array.Fill(ones, 1.0);
        double? sigma2 = null;
        if (Skew != 0.0)
        {
            ReadOnlySpan<double> y = response.Values;
            double mean = 0.0, squares = 0.0;
            for (int i = 0; i < n; i++) mean += y[i];
            mean /= n;
            for (int i = 0; i < n; i++) squares += (y[i] - mean) * (y[i] - mean);
            double perUnit = Skew * Skew * Mixing.Variance + VarianceScale * Mixing.Mean;
            sigma2 = _held ?? Math.Max(squares / Math.Max(n - 1, 1), double.Epsilon) / perUnit;
        }
        BassPriorReference? reference = Noise.Reference switch
        {
            BassNoiseReference.Moments => BassPriorReference.Moments(response.Values),
            BassNoiseReference.Robust => BassPriorReference.Robust(response.Values),
            BassNoiseReference.Quantile => BassPriorReference.Quantile(response.Values, Noise.Quantile!.Value),
            _ => null,
        };
        var latents = new ObservationLatents(ones, Skew, sigma2);
        (double[]? offset, double[] weights) = Reading(latents);
        return new WorkingResponse(response, offset, weights, reference: reference, latents: latents);
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

    /// <summary>a's block in the coefficients' slot, v's in the latent slot, σ²'s in the scale slot under a skew with σ² not held.</summary>
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
        => Mixing.MarginalLogDensity(residual, sigma2 * VarianceScale, Math.Sqrt(sigma2) * (latents?.Skew ?? 0.0));

    /// <inheritdoc/>
    public IContinuousDistribution At(double sigma2, ObservationLatents? latents) => Noise.Kind switch
    {
        BassNoiseKind.StudentT => StudentT.LocationScale(Noise.DegreesOfFreedom!.Value, 0.0, Math.Sqrt(sigma2 * VarianceScale)),
        BassNoiseKind.AsymmetricLaplace => AsymmetricLaplace.Of(0.0, Math.Sqrt(sigma2), Noise.Quantile!.Value),
        _ => throw new InvalidOperationException($"No distribution for {Noise.Kind}."),
    };

    /// <inheritdoc/>
    public double MeanShift(double sigma2, ObservationLatents? latents)
    {
        double skew = latents?.Skew ?? 0.0;
        return skew == 0.0 ? 0.0 : Math.Sqrt(sigma2) * skew * Mixing.Mean;
    }

    /// <inheritdoc/>
    public double Variance(double sigma2, ObservationLatents? latents)
    {
        double skew = latents?.Skew ?? 0.0;
        return sigma2 * ((skew == 0.0 ? 0.0 : skew * skew * Mixing.Variance) + VarianceScale * Mixing.Mean);
    }

    /// <inheritdoc/>
    public double Sample(double sigma2, ObservationLatents? latents, Xoshiro256PlusPlus rng)
    {
        double v = Mixing.SamplePrior(rng);
        return Math.Sqrt(sigma2) * ((latents?.Skew ?? 0.0) * v + Math.Sqrt(VarianceScale * v) * Normal.SampleStandard(rng));
    }
}

/// <summary>
/// The coefficients' block (D04's robust row, the coefficients' slot): σ² and a from their conditional given the structure, τ and the
/// working response, exact, at the cost of a retained draw's composition (<see cref="BassStructureMoves.Compose"/>), with the fit f at
/// every observation. They ride the chain's own moves to the blocks after it in the same sweep, which read the residuals y − f; being
/// redrawn each sweep they are not state. A σ² held or given is taken as it stands. The likelihood's power is 1: a tempered rung's
/// conditional is E1's slice 5.
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
        if (state.Temperature.InverseTemperature != 1.0)
            throw new NotSupportedException("A latent noise's blocks on a tempered rung are E1's slice 5; run a chain at b = 1.");
        BassStructureMoves moves = BassStructureMoves.Of(state);
        moves.DrawSweepFit(state.Current, rng, state.Target as BassTarget);
        return true;
    }
}

/// <summary>
/// The latent scales' block (the latent slot): each vᵢ from its conditional given the residual rᵢ = yᵢ − fᵢ of this sweep's coefficients
/// and σ² (<see cref="BassCoefficientsBlock"/>, which runs first), under the drift σβ and the unit variance σ²c, through the mixing law;
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
        ReadOnlySpan<double> y = working.Response.Values;
        double unitVariance = fit.Sigma2 * _model.VarianceScale, drift = Math.Sqrt(fit.Sigma2) * latents.Skew;
        IMixingLaw mixing = _model.Mixing;
        int n = y.Length;
        var v = new double[n];
        for (int i = 0; i < n; i++) v[i] = mixing.SampleConditional(y[i] - fit.Fitted[i], unitVariance, drift, 1.0, rng);
        WorkingResponse moved = _model.At(working, new ObservationLatents(v, latents.Skew, latents.Sigma2));
        return BassStructureMoves.Follow(state, moved) is null;
    }
}

/// <summary>
/// σ²'s block under a skew (the scale slot; E1's decision 2): σ² from its conditional given the coefficients of this sweep and the latents.
/// With rᵢ = yᵢ − fᵢ, the coefficients' prior variance scaling with σ², K of them inside the prior with energy Q = τ(a − a₀)ᵀD(a − a₀)
/// (<see cref="GPriorConditional.PriorEnergy"/>), and σ²'s prior IG(α₀, β₀) at the state (the hierarchical prior's at its rate; Jeffreys
/// at α₀ = β₀ = 0), u = σ⁻¹ is the modified half-normal MHN(N + K + 2α₀, ½(Σᵢ rᵢ²/(c·vᵢ) + Q + 2β₀), βΣᵢ rᵢ/c) (Sun, Kong and Pal 2021;
/// GBASS's w step with w's prior a kept, which GBASS drops, E1's slice 0). The chain follows the offset at the new σ², the weights
/// shared; a draw that puts the structure outside F is refused as the latent block's is.
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
        double shape = prior?.Shape ?? 0.0, scale = prior?.Scale ?? 0.0;
        double u = ModifiedHalfNormal.Of(y.Length + count + 2.0 * shape, 0.5 * (squares + quadratic + 2.0 * scale), latents.Skew * sum / c).Sample(rng);
        WorkingResponse moved = _model.At(working, latents.WithSigma2(1.0 / (u * u)), weightsMoved: false);
        return BassStructureMoves.Follow(state, moved) is null;
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
