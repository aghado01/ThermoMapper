using System;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.GPrior;

/// <summary>Which prior on σ² a <see cref="Sigma2Prior"/> is; there, the parameters each kind takes.</summary>
public enum Sigma2PriorKind
{
    /// <summary>IG(a, a·s²), s² read from the response: the default.</summary>
    ScaledToResponse,

    /// <summary>The improper p(σ²) ∝ 1/σ².</summary>
    Jeffreys,

    /// <summary>IG(a, b) in the response's own units.</summary>
    InverseGamma,

    /// <summary>σ² held at a declared value: a point mass, so the block's marginal is the one given σ².</summary>
    Held,

    /// <summary>σ² | r ~ IG(a, r) with r ~ Gamma(h, rate h/(a·s²)), s² read as the scaled prior reads it; r a coordinate of the
    /// chain's state.</summary>
    Hierarchical,
}

/// <summary>
/// The prior on the noise variance σ² of the collapsed g-prior engines (BARS now, BASS next). It resolves, per response,
/// to the <see cref="Distributions.InverseGamma"/> the block takes, or to null for the improper 1/σ².
/// <list type="bullet">
/// <item><see cref="Jeffreys"/>: p(σ²) ∝ 1/σ², scale-free and improper. The marginal is known up to a constant, so a
/// tempered marginal stops short of b = 0 and there is no evidence; a configuration that fits the response exactly has
/// no posterior.</item>
/// <item><see cref="ScaledToResponse"/>: IG(a, a·s²), with s² the response's mean square about the g-prior's centre
/// level, read from the data alone (never from robust weights) and the same for every design. Proper and scale-free:
/// rescaling y rescales s², and the posterior over designs is unchanged (OQ19's route b). a = 10⁻³ is the vague
/// IG(10⁻³, 10⁻³) on the standardized response that Bach &amp; Klein's (2025) code uses; a = ½ is a unit-information
/// prior, one observation's worth at the null fit's variance, which pulls σ² towards the response's own spread when the
/// fit is close. Normal BARS's default, at a = 10⁻³.</item>
/// <item><see cref="InverseGamma"/>: IG(a, b) in the response's own units, so b is small or large depending on them
/// (OQ19); for parity with an analysis that fixes it.</item>
/// <item><see cref="Held"/>: σ² held at a declared value, as probit's identification holds it at 1 or a known noise
/// states it. BASS's structure step then reads the marginal given σ² (<see cref="GPriorMarginal.LogMarginalGiven"/>);
/// BARS, which integrates σ², refuses it.</item>
/// <item><see cref="Hierarchical"/>: σ² | r ~ IG(a, r) with the rate r ~ Gamma(h, rate h/(a·s²)), whose mean is the scaled prior's
/// a·s² at every h, s² read as the scaled prior reads it. Proper and scale-free, with the family's two limits at its ends: h → ∞
/// is the scaled prior, and h → 0 is p(σ²) ∝ 1/σ² for any a, since IG(a, r) mixed over 1/r gives it. a = h = ½ is σ ~ half-Cauchy(s)
/// (<see cref="HalfCauchy"/>; Wand et al. 2011's two-block form). BASS's state carries r and a block in the sweep's scale slot moves
/// it, the structure's marginal being given r as it is given τ; BARS, which has no such block, refuses it.</item>
/// </list>
/// A proper prior makes the marginal complete: defined at every b in [0, 1], where b = 0 gives exactly 0, finite for an
/// exact fit, and absolute, so it is an evidence. The scaled prior's evidence depends on a and on s², as every proper
/// prior's does on its constants.
/// <para>A value: its kind and that kind's parameters, compared by value and serialized as they stand, so a config that
/// holds it serializes and fingerprints whole (the BASS contract shapes, F-C2). A parameter the kind does not take is
/// null, not NaN, so the JSON conventions that refuse named floats read it too. The constructor checks the parameters
/// against the kind, so a prior read back from JSON is checked as a built one is.</para>
/// </summary>
public sealed record Sigma2Prior
{
    /// <summary>
    /// The prior of kind <paramref name="kind"/>: <see cref="Sigma2PriorKind.Jeffreys"/> takes no parameter,
    /// <see cref="Sigma2PriorKind.ScaledToResponse"/> a shape, and <see cref="Sigma2PriorKind.InverseGamma"/> a shape and
    /// a scale, <see cref="Sigma2PriorKind.Held"/> a value, and <see cref="Sigma2PriorKind.Hierarchical"/> a shape and a hyper-shape,
    /// each positive and finite.
    /// </summary>
    public Sigma2Prior(Sigma2PriorKind kind, double? shape = null, double? scale = null, double? value = null, double? hyperShape = null)
    {
        if (kind != Sigma2PriorKind.Held && value is not null) throw new ArgumentException("Only the held prior takes a value.", nameof(value));
        if (kind != Sigma2PriorKind.Hierarchical && hyperShape is not null) throw new ArgumentException("Only the hierarchical prior takes a hyper-shape.", nameof(hyperShape));
        switch (kind)
        {
            case Sigma2PriorKind.Jeffreys:
                if (shape is not null || scale is not null) throw new ArgumentException("The Jeffreys prior takes no parameter.");
                break;
            case Sigma2PriorKind.ScaledToResponse:
                if (shape is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(shape), "The shape must be positive and finite.");
                if (scale is not null) throw new ArgumentException("The scaled prior reads its scale from the response.", nameof(scale));
                break;
            case Sigma2PriorKind.InverseGamma:
                if (shape is not { } a || scale is not { } b) throw new ArgumentException("The inverse gamma takes a shape and a scale.");
                _ = Distributions.InverseGamma.ShapeScale(a, b);   // validates both
                break;
            case Sigma2PriorKind.Held:
                if (shape is not null || scale is not null) throw new ArgumentException("The held prior takes a value alone.");
                if (value is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(value), "The held σ² must be positive and finite.");
                break;
            case Sigma2PriorKind.Hierarchical:
                if (shape is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(shape), "The shape must be positive and finite.");
                if (hyperShape is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(hyperShape), "The hyper-shape must be positive and finite.");
                if (scale is not null) throw new ArgumentException("The hierarchical prior reads its scale from the response.", nameof(scale));
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(kind));
        }
        Kind = kind;
        Shape = shape;
        Scale = scale;
        Value = value;
        HyperShape = hyperShape;
    }

    /// <summary>The improper p(σ²) ∝ 1/σ².</summary>
    public static Sigma2Prior Jeffreys { get; } = new(Sigma2PriorKind.Jeffreys);

    /// <summary>IG(<paramref name="shape"/>, shape·s²), s² the response's mean square about the centre level.</summary>
    public static Sigma2Prior ScaledToResponse(double shape = 1e-3) => new(Sigma2PriorKind.ScaledToResponse, shape);

    /// <summary>IG(<paramref name="shape"/>, <paramref name="scale"/>) in the response's own units.</summary>
    public static Sigma2Prior InverseGamma(double shape, double scale) => new(Sigma2PriorKind.InverseGamma, shape, scale);

    /// <summary>σ² held at <paramref name="value"/>.</summary>
    public static Sigma2Prior Held(double value) => new(Sigma2PriorKind.Held, value: value);

    /// <summary>σ² | r ~ IG(<paramref name="shape"/>, r), r ~ Gamma(<paramref name="hyperShape"/>, rate h/(a·s²)).</summary>
    public static Sigma2Prior Hierarchical(double shape, double hyperShape) => new(Sigma2PriorKind.Hierarchical, shape, hyperShape: hyperShape);

    /// <summary>σ ~ half-Cauchy(s): the hierarchical prior at a = h = ½.</summary>
    public static Sigma2Prior HalfCauchy() => Hierarchical(0.5, 0.5);

    /// <summary>Which prior this is.</summary>
    public Sigma2PriorKind Kind { get; }

    /// <summary>The shape a; null for <see cref="Jeffreys"/>.</summary>
    public double? Shape { get; }

    /// <summary>The fixed scale b of <see cref="InverseGamma"/>; null otherwise.</summary>
    public double? Scale { get; }

    /// <summary>The held σ²; null otherwise, and left out of the serialized form then, so the other kinds' forms stand as they were.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? Value { get; }

    /// <summary>The rate's shape h; null otherwise, and left out of the serialized form then, as <see cref="Value"/> is.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double? HyperShape { get; }

    /// <summary>Whether σ² is held; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool IsHeld => Kind == Sigma2PriorKind.Held;

    /// <summary>True for every prior but <see cref="Jeffreys"/>; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool IsProper => Kind != Sigma2PriorKind.Jeffreys;

    /// <summary>Whether the state carries σ²'s rate; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool IsHierarchical => Kind == Sigma2PriorKind.Hierarchical;

    /// <summary>
    /// The prior for a response whose mean square about the centre level is <paramref name="meanSquare"/>: null for
    /// <see cref="Jeffreys"/>, the fixed IG, or IG(a, a·s²). The scaled prior has none for a response with no spread
    /// (s² = 0), which has no posterior for σ either; it returns null there and the caller refuses the response. A held σ²
    /// has no inverse gamma either: null, with <see cref="Value"/> for its caller. A hierarchical prior's is the inverse gamma at
    /// its rate's mean, which <see cref="RateFor"/> stands behind.
    /// </summary>
    public Distributions.InverseGamma? For(double meanSquare) => (Kind, Shape, Scale) switch
    {
        (Sigma2PriorKind.InverseGamma, { } a, { } b) => Distributions.InverseGamma.ShapeScale(a, b),
        (Sigma2PriorKind.ScaledToResponse or Sigma2PriorKind.Hierarchical, { } a, _) when meanSquare > 0.0 && double.IsFinite(meanSquare)
            => Distributions.InverseGamma.ShapeScale(a, a * meanSquare),
        _ => null,
    };

    /// <summary>The hierarchical prior's rate on the mean square <paramref name="meanSquare"/>; null for every other kind, and
    /// where <see cref="For"/> has none.</summary>
    public ResolvedSigma2Rate? RateFor(double meanSquare)
        => Kind == Sigma2PriorKind.Hierarchical && meanSquare > 0.0 && double.IsFinite(meanSquare)
            ? new ResolvedSigma2Rate(Shape!.Value, HyperShape!.Value, meanSquare)
            : null;

    /// <summary>Whether <see cref="For"/> reads the response's mean square; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool ReadsResponse => Kind is Sigma2PriorKind.ScaledToResponse or Sigma2PriorKind.Hierarchical;

    /// <summary>The spelling BARS's checkpoint fingerprint reads.</summary>
    public override string ToString() => Kind switch
    {
        Sigma2PriorKind.Jeffreys => "Jeffreys",
        Sigma2PriorKind.ScaledToResponse => FormattableString.Invariant($"ScaledToResponse({Shape:R})"),
        Sigma2PriorKind.Held => FormattableString.Invariant($"Held({Value:R})"),
        Sigma2PriorKind.Hierarchical => FormattableString.Invariant($"Hierarchical({Shape:R}, {HyperShape:R})"),
        _ => FormattableString.Invariant($"InverseGamma({Shape:R}, {Scale:R})"),
    };
}

/// <summary>
/// A hierarchical σ² prior's rate on its mean square s² (<see cref="Sigma2Prior.RateFor"/>): r ~ Gamma(h, rate h/(a·s²)), with
/// the density of v = ln r its Jacobian included, its exact draw and the inverse gamma it makes of σ²'s prior, the three that the
/// rate's block, the reference rung and SMC's start read. The density is the one the block's slice step reads on v, so the step
/// adds no Jacobian of its own, as <see cref="ResolvedTauPrior"/>'s is for τ. All of it stands on v, since under a tiny h most of
/// r's prior lies below the doubles: at h = 10⁻³ about half. Only the inverse gamma needs r itself, and it reads r no lower than the
/// doubles' normal range, with <see cref="LogRemainder"/> for the rest.
/// </summary>
public sealed record ResolvedSigma2Rate
{
    private readonly double _logConstant;    // h ln(rate) − ln Γ(h)
    private readonly double _logGammaShape;  // ln Γ(a)

    /// <summary>The lowest v at which the inverse gamma reads r, ln of the smallest normal double.</summary>
    public static readonly double LogSmallestRate = Math.Log(2.2250738585072014E-308);

    internal ResolvedSigma2Rate(double shape, double hyperShape, double meanSquare)
    {
        Shape = shape;
        Centre = shape * meanSquare;
        Hyper = Gamma.ShapeRate(hyperShape, hyperShape / Centre);
        _logConstant = hyperShape * Math.Log(Hyper.Rate) - SpecialFunctions.LogGamma(hyperShape);
        _logGammaShape = SpecialFunctions.LogGamma(shape);
    }

    /// <summary>a, the inverse gamma's shape.</summary>
    public double Shape { get; }

    /// <summary>r's prior.</summary>
    public Gamma Hyper { get; }

    /// <summary>r's mean, a·s²: the scaled prior's rate, and where a chain starts unless it is told otherwise.</summary>
    public double Centre { get; }

    /// <summary>v = ln of <see cref="Centre"/>.</summary>
    public double CentreLogRate => Math.Log(Centre);

    /// <summary>
    /// The log density of v = ln r: h ln(rate) − ln Γ(h) + h·v − rate·eᵛ, which stands below the doubles. Where r overflows it is
    /// outside the support: −∞, as for τ.
    /// </summary>
    public double LogDensity(double logRate)
    {
        if (double.IsNaN(logRate)) throw new ArgumentOutOfRangeException(nameof(logRate));
        double scaled = Hyper.Rate * Math.Exp(logRate);
        return double.IsFinite(scaled) ? _logConstant + Hyper.Shape * logRate - scaled : double.NegativeInfinity;
    }

    /// <summary>r at v, no lower than the smallest normal double: what a draw retains and what <see cref="PriorAt"/> reads.</summary>
    public double RateAt(double logRate) => Math.Exp(Math.Max(logRate, LogSmallestRate));

    /// <summary>σ²'s prior at v, IG(a, <see cref="RateAt"/>).</summary>
    public InverseGamma PriorAt(double logRate) => InverseGamma.ShapeScale(Shape, RateAt(logRate));

    /// <summary>
    /// What the marginal over σ² gains at v beyond <see cref="PriorAt"/>'s reading, once b &gt; 0: a·(v − ln r_min) below the smallest normal
    /// double r_min, 0 above it. The marginal's r enters as rᵃ(r + S/2)^(−a−ν/2), and S, far above r, cannot tell r from r_min.
    /// </summary>
    public double LogRemainder(double logRate) => Shape * Math.Min(logRate - LogSmallestRate, 0.0);

    /// <summary>ln IG(<paramref name="sigma2"/>; a, eᵛ), which the state's prior carries when the structure's marginal is given σ²,
    /// so the rate's conditional still reads σ².</summary>
    public double LogConditional(double sigma2, double logRate)
        => Shape * logRate - _logGammaShape - (Shape + 1.0) * Math.Log(sigma2) - Math.Exp(logRate) / sigma2;

    /// <summary>An exact draw of v, below the doubles too; a draw whose r overflows is drawn again.</summary>
    public double SampleLogRate(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double v;
        do v = Hyper.SampleLog(rng);
        while (double.IsNegativeInfinity(LogDensity(v)));
        return v;
    }

    /// <summary>The spelling a manifest or a log reads.</summary>
    public override string ToString() => FormattableString.Invariant($"IG({Shape:R}, r), r ~ Gamma({Hyper.Shape:R}, rate {Hyper.Rate:R})");
}
