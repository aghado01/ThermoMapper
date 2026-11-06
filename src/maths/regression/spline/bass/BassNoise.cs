using System;
using System.Collections.Generic;
using System.Text.Json.Serialization;
using Maths.Distributions;
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
}

/// <summary>What the priors read of the response under a noise whose working response moves (<see cref="BassPriorReference"/>).</summary>
public enum BassNoiseReference
{
    /// <summary>The mean and the mean square about it, D03's read: at β = 0, under the unweighted prior Gram, what Gaussian BASS reads.</summary>
    Moments,

    /// <summary>The median and the normal-consistent MAD², for heavy tails: an outlier cannot move the level or the scale.</summary>
    Robust,

    /// <summary>Nothing: for priors that read no response, a flat or zero intercept with a σ² prior that reads no s². Refused otherwise,
    /// at the first move (<see cref="BassModel.PriorsReadResponse"/>).</summary>
    None,
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
    /// positive and finite. <paramref name="varianceScale"/> is c, positive and finite: 1 makes σ the t's scale, BARS's
    /// <c>StudentTWeights</c>' t; (ν − 2)/ν makes σ² the noise variance, GBASS's <c>tbass</c> for ν &gt; 2. <paramref name="reference"/>
    /// is what the priors read of the response.
    /// </summary>
    public BassNoise(BassNoiseKind kind, double? degreesOfFreedom = null, double varianceScale = 1.0,
                     BassNoiseReference reference = BassNoiseReference.Moments)
    {
        switch (kind)
        {
            case BassNoiseKind.StudentT:
                if (degreesOfFreedom is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "Student t takes ν, positive and finite.");
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(kind));
        }
        if (!(varianceScale > 0.0 && double.IsFinite(varianceScale))) throw new ArgumentOutOfRangeException(nameof(varianceScale), "c must be positive and finite.");
        if (reference is not (BassNoiseReference.Moments or BassNoiseReference.Robust or BassNoiseReference.None))
            throw new ArgumentOutOfRangeException(nameof(reference));
        Kind = kind;
        DegreesOfFreedom = degreesOfFreedom;
        VarianceScale = varianceScale;
        Reference = reference;
    }

    /// <summary>Student t on <paramref name="degreesOfFreedom"/> ν, with c = <paramref name="varianceScale"/>.</summary>
    public static BassNoise StudentT(double degreesOfFreedom, double varianceScale = 1.0, BassNoiseReference reference = BassNoiseReference.Moments)
        => new(BassNoiseKind.StudentT, degreesOfFreedom, varianceScale, reference);

    /// <summary>Which noise this is.</summary>
    public BassNoiseKind Kind { get; }

    /// <summary>Student t's ν; null for a kind that takes none.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? DegreesOfFreedom { get; }

    /// <summary>c, the variance of the Gaussian part per unit σ²·v.</summary>
    public double VarianceScale { get; }

    /// <summary>What the priors read of the response.</summary>
    public BassNoiseReference Reference { get; }

    /// <summary>The kind's mixing law, read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public IMixingLaw Mixing => Kind switch
    {
        BassNoiseKind.StudentT => new StudentTMixing(DegreesOfFreedom!.Value),
        _ => throw new InvalidOperationException($"No mixing law for {Kind}."),
    };
}

/// <summary>
/// An observation model's latents at one state, carried by the <see cref="WorkingResponse"/> they make (E1's decision 1): one value per
/// observation, a variance-mean mixture's scale v (a probit's z or a Pólya–Gamma ω would ride the same slot), and the skew β, 0 for a
/// symmetric law. A retained draw keeps the working response it was composed on, so it keeps these too. Immutable; the values are copied.
/// </summary>
public sealed class ObservationLatents
{
    private readonly double[] _values;

    /// <summary>The latents <paramref name="values"/>, one per observation and each finite, with skew <paramref name="skew"/>.</summary>
    public ObservationLatents(ReadOnlySpan<double> values, double skew = 0.0)
        : this(values.ToArray(), skew) { }

    // Takes `values` as its own.
    internal ObservationLatents(double[] values, double skew)
    {
        foreach (double v in values)
            if (!double.IsFinite(v)) throw new ArgumentException("A latent is finite.", nameof(values));
        if (!double.IsFinite(skew)) throw new ArgumentOutOfRangeException(nameof(skew), "The skew is finite.");
        _values = values;
        Skew = skew;
    }

    /// <summary>One latent per observation, case by case and location by location.</summary>
    public ReadOnlySpan<double> Values => _values;

    /// <summary>β, the drift per unit v; 0 for a symmetric law.</summary>
    public double Skew { get; }

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
/// <see cref="BassNoise.Reference"/> frozen at the start. A chain starts at v = 1. The sweep (D04's robust row): the structure with σ²
/// and a integrated out, τ, then (σ², a) from their conditional (<see cref="BassCoefficientsBlock"/>), then each vᵢ from its own
/// (<see cref="BassLatentScaleBlock"/>); under a pure scale mixture σ² is redrawn with a each sweep, so the scale slot is empty. Student t
/// is the member built (E1's slice 1); its draws need no new distribution. One model serves every chain: it holds no chain's state.
/// </summary>
public sealed class VarianceMeanMixture : IBassObservationModel, IBassNoiseLaw
{
    private readonly IBlock<ReversibleJumpChain<BassState>>[] _coefficients, _latent;

    /// <summary>The model of <paramref name="noise"/>.</summary>
    public VarianceMeanMixture(BassNoise noise)
    {
        ArgumentNullException.ThrowIfNull(noise);
        Noise = noise;
        Mixing = noise.Mixing;
        _coefficients = [new BassCoefficientsBlock()];
        _latent = [new BassLatentScaleBlock(this)];
    }

    /// <summary>The noise it models.</summary>
    public BassNoise Noise { get; }

    /// <summary>The mixing law of v.</summary>
    public IMixingLaw Mixing { get; }

    /// <summary>c.</summary>
    public double VarianceScale => Noise.VarianceScale;

    /// <summary>The working response at v = 1: the response itself under weights 1/c, the reference read from it.</summary>
    public WorkingResponse Working(BassResponse response)
    {
        ArgumentNullException.ThrowIfNull(response);
        if (response is FunctionalResponse)
            throw new NotSupportedException("A per-observation latent under a direct functional design is E1's slice 5; a scalar response takes it now.");
        int n = response.Observations;
        var weights = new double[n];
        Array.Fill(weights, 1.0 / VarianceScale);
        var ones = new double[n];
        Array.Fill(ones, 1.0);
        BassPriorReference? reference = Noise.Reference switch
        {
            BassNoiseReference.Moments => BassPriorReference.Moments(response.Values),
            BassNoiseReference.Robust => BassPriorReference.Robust(response.Values),
            _ => null,
        };
        return new WorkingResponse(response, weights: weights, reference: reference, latents: new ObservationLatents(ones, 0.0));
    }

    /// <summary>a's block in the coefficients' slot, v's in the latent slot; the scale slot is empty under a pure scale mixture.</summary>
    public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => slot switch
    {
        BassSweepSlot.Coefficients => _coefficients,
        BassSweepSlot.Latent => _latent,
        _ => [],
    };

    /// <summary>The noise's own law.</summary>
    public IBassNoiseLaw Law => this;

    /// <inheritdoc/>
    public double LogDensity(double residual, double sigma2, ObservationLatents? latents)
        => Mixing.MarginalLogDensity(residual, sigma2 * VarianceScale, Math.Sqrt(sigma2) * (latents?.Skew ?? 0.0));

    /// <inheritdoc/>
    public IContinuousDistribution At(double sigma2, ObservationLatents? latents)
    {
        if ((latents?.Skew ?? 0.0) != 0.0) throw new NotSupportedException("A skewed noise's distribution arrives with E1's slices 2 and 3.");
        return Noise.Kind switch
        {
            BassNoiseKind.StudentT => StudentT.LocationScale(Noise.DegreesOfFreedom!.Value, 0.0, Math.Sqrt(sigma2 * VarianceScale)),
            _ => throw new InvalidOperationException($"No distribution for {Noise.Kind}."),
        };
    }

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
/// The latent scales' block (the latent slot): each vᵢ from its conditional given the residual yᵢ − fᵢ of this sweep's coefficients and
/// σ² (<see cref="BassCoefficientsBlock"/>, which runs first), through the mixing law, then the chain follows the working response at the
/// new v (<see cref="BassStructureMoves.Follow"/>). A draw that puts the current structure outside F is refused and the chain keeps its
/// state, the exact step under the gate that reads W (decision 4); the scan's counters read its refusals.
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
        ReadOnlySpan<double> y = working.Values;
        double c = _model.VarianceScale, unitVariance = fit.Sigma2 * c;
        IMixingLaw mixing = _model.Mixing;
        int n = y.Length;
        var v = new double[n];
        var weights = new double[n];
        for (int i = 0; i < n; i++)
        {
            v[i] = mixing.SampleConditional(y[i] - fit.Fitted[i], unitVariance, 0.0, 1.0, rng);
            weights[i] = 1.0 / (c * v[i]);
        }
        WorkingResponse moved = working.Moved(working.Offset, weights, new ObservationLatents(v, 0.0));
        return BassStructureMoves.Follow(state, moved) is null;
    }
}

/// <summary>A sweep's draw of σ² and the coefficients at <paramref name="State"/>: the coefficients (α, β_c) in the design's centred order,
/// and the fit f at every observation. Not state: the blocks after the coefficients' read it within the sweep.</summary>
internal sealed record BassSweepFit(BassState State, double Sigma2, double[] Centred, double[] Fitted);
