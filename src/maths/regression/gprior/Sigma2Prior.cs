using System;
using System.Text.Json.Serialization;

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
    /// a scale, each positive and finite.
    /// </summary>
    public Sigma2Prior(Sigma2PriorKind kind, double? shape = null, double? scale = null)
    {
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
            default:
                throw new ArgumentOutOfRangeException(nameof(kind));
        }
        Kind = kind;
        Shape = shape;
        Scale = scale;
    }

    /// <summary>The improper p(σ²) ∝ 1/σ².</summary>
    public static Sigma2Prior Jeffreys { get; } = new(Sigma2PriorKind.Jeffreys);

    /// <summary>IG(<paramref name="shape"/>, shape·s²), s² the response's mean square about the centre level.</summary>
    public static Sigma2Prior ScaledToResponse(double shape = 1e-3) => new(Sigma2PriorKind.ScaledToResponse, shape);

    /// <summary>IG(<paramref name="shape"/>, <paramref name="scale"/>) in the response's own units.</summary>
    public static Sigma2Prior InverseGamma(double shape, double scale) => new(Sigma2PriorKind.InverseGamma, shape, scale);

    /// <summary>Which prior this is.</summary>
    public Sigma2PriorKind Kind { get; }

    /// <summary>The shape a; null for <see cref="Jeffreys"/>.</summary>
    public double? Shape { get; }

    /// <summary>The fixed scale b of <see cref="InverseGamma"/>; null otherwise.</summary>
    public double? Scale { get; }

    /// <summary>True for every prior but <see cref="Jeffreys"/>; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool IsProper => Kind != Sigma2PriorKind.Jeffreys;

    /// <summary>
    /// The prior for a response whose mean square about the centre level is <paramref name="meanSquare"/>: null for
    /// <see cref="Jeffreys"/>, the fixed IG, or IG(a, a·s²). The scaled prior has none for a response with no spread
    /// (s² = 0), which has no posterior for σ either; it returns null there and the caller refuses the response.
    /// </summary>
    public Distributions.InverseGamma? For(double meanSquare) => (Kind, Shape, Scale) switch
    {
        (Sigma2PriorKind.InverseGamma, { } a, { } b) => Distributions.InverseGamma.ShapeScale(a, b),
        (Sigma2PriorKind.ScaledToResponse, { } a, _) when meanSquare > 0.0 && double.IsFinite(meanSquare)
            => Distributions.InverseGamma.ShapeScale(a, a * meanSquare),
        _ => null,
    };

    /// <summary>Whether <see cref="For"/> reads the response's mean square; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool ReadsResponse => Kind == Sigma2PriorKind.ScaledToResponse;

    /// <summary>The spelling BARS's checkpoint fingerprint reads.</summary>
    public override string ToString() => Kind switch
    {
        Sigma2PriorKind.Jeffreys => "Jeffreys",
        Sigma2PriorKind.ScaledToResponse => FormattableString.Invariant($"ScaledToResponse({Shape:R})"),
        _ => FormattableString.Invariant($"InverseGamma({Shape:R}, {Scale:R})"),
    };
}
