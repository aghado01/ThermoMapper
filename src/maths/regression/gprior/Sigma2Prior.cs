using System;

namespace Maths.Regression.GPrior;

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
/// fit is close.</item>
/// <item><see cref="InverseGamma"/>: IG(a, b) in the response's own units, so b is small or large depending on them
/// (OQ19); for parity with an analysis that fixes it.</item>
/// </list>
/// A proper prior makes the marginal complete: defined at every b in [0, 1], where b = 0 gives exactly 0, finite for an
/// exact fit, and absolute, so it is an evidence. The scaled prior's evidence depends on a and on s², as every proper
/// prior's does on its constants.
/// </summary>
public sealed class Sigma2Prior
{
    private enum Kind { Jeffreys, Scaled, Fixed }

    private readonly Kind _kind;

    private Sigma2Prior(Kind kind, double shape, double scale)
    {
        _kind = kind;
        Shape = shape;
        Scale = scale;
    }

    /// <summary>The improper p(σ²) ∝ 1/σ².</summary>
    public static Sigma2Prior Jeffreys { get; } = new(Kind.Jeffreys, double.NaN, double.NaN);

    /// <summary>IG(<paramref name="shape"/>, shape·s²), s² the response's mean square about the centre level.</summary>
    public static Sigma2Prior ScaledToResponse(double shape = 1e-3)
    {
        if (!(shape > 0.0) || double.IsPositiveInfinity(shape))
            throw new ArgumentOutOfRangeException(nameof(shape), "The shape must be positive and finite.");
        return new Sigma2Prior(Kind.Scaled, shape, double.NaN);
    }

    /// <summary>IG(<paramref name="shape"/>, <paramref name="scale"/>) in the response's own units.</summary>
    public static Sigma2Prior InverseGamma(double shape, double scale)
    {
        _ = Distributions.InverseGamma.ShapeScale(shape, scale);   // validates both
        return new Sigma2Prior(Kind.Fixed, shape, scale);
    }

    /// <summary>True for every prior but <see cref="Jeffreys"/>.</summary>
    public bool IsProper => _kind != Kind.Jeffreys;

    /// <summary>The shape a; NaN for <see cref="Jeffreys"/>.</summary>
    public double Shape { get; }

    /// <summary>The fixed scale b of <see cref="InverseGamma"/>; NaN otherwise.</summary>
    public double Scale { get; }

    /// <summary>
    /// The prior for a response whose mean square about the centre level is <paramref name="meanSquare"/>: null for
    /// <see cref="Jeffreys"/>, the fixed IG, or IG(a, a·s²). The scaled prior has none for a response with no spread
    /// (s² = 0), which has no posterior for σ either; it returns null there and the caller refuses the response.
    /// </summary>
    public Distributions.InverseGamma? For(double meanSquare) => _kind switch
    {
        Kind.Fixed => Distributions.InverseGamma.ShapeScale(Shape, Scale),
        Kind.Scaled when meanSquare > 0.0 && double.IsFinite(meanSquare) => Distributions.InverseGamma.ShapeScale(Shape, Shape * meanSquare),
        _ => null,
    };

    /// <summary>Whether <see cref="For"/> reads the response's mean square.</summary>
    public bool ReadsResponse => _kind == Kind.Scaled;

    /// <inheritdoc/>
    public override string ToString() => _kind switch
    {
        Kind.Jeffreys => "Jeffreys",
        Kind.Scaled => FormattableString.Invariant($"ScaledToResponse({Shape:R})"),
        _ => FormattableString.Invariant($"InverseGamma({Shape:R}, {Scale:R})"),
    };
}
