using System;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// KHAOS's modified g-prior (E5; Rumsey, Francom et al., arXiv 2510.25036): (β₀, β) | σ², g₀² ~ N(0, σ²S₀) on the columns as
/// evaluated, the intercept's column of ones among them, with S₀ = g₀²·D_g(ΨᵀΨ)⁻¹D_g and D_g = diag(1, g₁, …, g_M). Each term's
/// complexity factor g_m = (1 + q_m(d_m + q_m − 2))^(−ζ/2) reads its order q_m and degree d_m, so a richer term's coefficient is
/// shrunk harder; the intercept's is 1. g₀² is the coefficient prior's τ inverted, τ = 1/g₀², and <c>khaos</c>'s IG(a, b) on it is
/// <see cref="TauPrior.Gamma"/> on τ. ζ = 0 is the g-prior with the intercept inside it under the zero centre, the g-form's
/// (<see cref="GForm"/>). <c>khaos</c>'s default ζ is 1. A Legendre basis's: the factor reads a Legendre term's order and degree.
/// </summary>
public sealed record ModifiedGForm : CoefficientForm
{
    /// <summary>The prior with complexity exponent <paramref name="zeta"/>, finite and not negative.</summary>
    public ModifiedGForm(double zeta = 1.0)
    {
        if (!(zeta >= 0.0 && double.IsFinite(zeta))) throw new ArgumentOutOfRangeException(nameof(zeta), "ζ must be finite and not negative.");
        Zeta = zeta;
    }

    /// <summary>ζ, the complexity factors' exponent.</summary>
    public double Zeta { get; }

    /// <summary>g_m = (1 + q(d + q − 2))^(−ζ/2) for <paramref name="term"/>'s Legendre part, of order q and degree d.</summary>
    public double Factor(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        int q = term.Legendre.Count;
        if (q == 0) throw new ArgumentException("The modified g-prior's factor reads a Legendre term's order and degree.", nameof(term));
        int d = term.LegendreDegree;
        return Math.Pow(1.0 + q * (d + q - 2.0), -0.5 * Zeta);
    }
}

/// <summary>
/// How the coefficient prior's precision reads a structure's columns (E21): the builder a design holds, one per coefficient form,
/// in place of a flag per form. In the design's coordinates, (α, β) on the columns centred under W, a prior that stands apart
/// from the likelihood's Gram has precision D_α(1, δ)(1, δ)ᵀ + diag(0, D_β), with α_u = α + δᵀβ the intercept on the prior's own
/// centring and δ the prior side's centres less the likelihood's (<see cref="PriorPencil"/>):
/// <list type="bullet">
/// <item><see cref="Gram"/>, the g-form: D_α = N and D_β = G_u, the prior side's centred Gram.</item>
/// <item><see cref="NormalizedRidge"/>: D_α = 1 and D_β = diag(G_u)/N.</item>
/// <item><see cref="RawRidge"/>, KHAOS's: D_α = 1, D_β = I and δ = −z̄_w, its coefficients on the columns as evaluated.</item>
/// <item><see cref="ComplexityGram"/>, KHAOS's modified g-prior: the g-form's on the prior-side columns scaled by 1/g_m. With
/// D_g's inverse applied to the raw coefficients, the intercept's column untouched, the prior reads z_m/g_m, so δ = z̄_u/g − z̄_w and
/// D_β = D_g⁻¹G_uD_g⁻¹, D_α = N.</item>
/// </list>
/// The g-form and the normalized ridge stand apart only under weights whose prior Gram is not the likelihood's
/// (<see cref="WorkingResponse.PriorStandsApart"/>), and keep their own forms otherwise, a factor and a spectrum, whose draws
/// stand as they were. KHAOS's two always stand apart, so each structure costs its own decomposition, after which a τ costs
/// O(p): g₀²'s slice step is τ's on the pencil.
/// </summary>
public abstract class BassPriorPrecision
{
    private protected BassPriorPrecision() { }

    /// <summary>The g-form's.</summary>
    public static BassPriorPrecision Gram { get; } = new GramPrecision();

    /// <summary>The normalized ridge's.</summary>
    public static BassPriorPrecision NormalizedRidge { get; } = new NormalizedRidgePrecision();

    /// <summary>KHAOS's raw ridge's.</summary>
    public static BassPriorPrecision RawRidge { get; } = new RawRidgePrecision();

    /// <summary>KHAOS's modified g-prior's, the factors of <paramref name="form"/>.</summary>
    public static BassPriorPrecision ComplexityGram(ModifiedGForm form) => new ComplexityGramPrecision(form ?? throw new ArgumentNullException(nameof(form)));

    /// <summary>The builder <paramref name="form"/> declares.</summary>
    public static BassPriorPrecision Of(CoefficientForm form) => form switch
    {
        GForm => Gram,
        RidgeForm => NormalizedRidge,
        RawRidgeForm => RawRidge,
        ModifiedGForm modified => ComplexityGram(modified),
        null => throw new ArgumentNullException(nameof(form)),
        _ => throw new ArgumentException($"{form.GetType().Name} has no prior-precision builder.", nameof(form)),
    };

    /// <summary>Whether a design keeps the normalized ridge's spectrum where the prior does not stand apart.</summary>
    internal virtual bool IsRidge => false;

    /// <summary>Whether the prior stands apart from the Gram at any weights.</summary>
    internal virtual bool AlwaysApart => false;

    /// <summary>Whether the angle gate reads the prior side's Gram, which is then a precision to keep well conditioned.</summary>
    internal virtual bool Gated => true;

    /// <summary>Whether the prior side holds no columns, only each column's centre (the raw ridge's D = I).</summary>
    internal virtual bool Unit => false;

    /// <summary>D_α, from the prior side's mass Σw₀ (N without frozen weights).</summary>
    internal virtual double Mass(double priorMass) => priorMass;

    /// <summary>D_β's (i, j) entry from the prior side's Gram <paramref name="u"/>.</summary>
    internal virtual double Entry(double[,] u, int i, int j, double priorMass) => u[i, j];

    /// <summary>The prior side's scale on <paramref name="term"/>'s column: 1, or 1/g_m.</summary>
    internal virtual double Scale(BassTerm term) => 1.0;

    private sealed class GramPrecision : BassPriorPrecision { }

    private sealed class NormalizedRidgePrecision : BassPriorPrecision
    {
        internal override bool IsRidge => true;
        internal override bool Gated => false;
        internal override double Mass(double priorMass) => 1.0;
        internal override double Entry(double[,] u, int i, int j, double priorMass) => i == j ? u[i, i] / priorMass : 0.0;
    }

    private sealed class RawRidgePrecision : BassPriorPrecision
    {
        internal override bool AlwaysApart => true;
        internal override bool Gated => false;
        internal override bool Unit => true;
        internal override double Mass(double priorMass) => 1.0;
        internal override double Entry(double[,] u, int i, int j, double priorMass) => i == j ? 1.0 : 0.0;
    }

    private sealed class ComplexityGramPrecision(ModifiedGForm form) : BassPriorPrecision
    {
        internal override bool AlwaysApart => true;
        internal override double Scale(BassTerm term) => 1.0 / form.Factor(term);
    }
}
