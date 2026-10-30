using System;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.GPrior;

/// <summary>Which prior on τ a <see cref="TauPrior"/> is; there, the parameters each kind takes.</summary>
public enum TauPriorKind
{
    /// <summary>τ ~ Gamma(a, rate a/τ₀), mean τ₀: the unit-information Gamma. a = ½ is Zellner–Siow, the default.</summary>
    UnitInformationGamma,

    /// <summary>Liang et al.'s hyper-g on g = 1/τ.</summary>
    HyperG,

    /// <summary>Liang et al.'s hyper-g/n: the hyper-g law on τ₀/τ.</summary>
    HyperGOverN,

    /// <summary>τ held at a declared value, or at τ₀.</summary>
    Held,
}

/// <summary>
/// The prior on τ, the coefficient prior's precision multiplier, shared by the collapsed g-prior engines (BASS now, BARS's
/// fixed unit-information g later; decided 2026-10-29). Every member is stated on the coefficient form's unit-information
/// precision τ₀ = 1/n₀, where its source states it so: τ₀ = 1/N for the g-form σ²G⁻¹/τ, and 1 for the ridge on centred
/// columns at unit root mean square. <see cref="For"/> resolves it on n₀.
/// <list type="bullet">
/// <item><see cref="UnitInformationGamma"/>: τ ~ Gamma(a, rate a/τ₀), whose mean is τ₀ at every shape; a = ½ is
/// Zellner–Siow (g = 1/τ ~ IG(½, n₀/2) under the g-form), BASS's default (D03). A larger a tightens τ towards τ₀, with
/// the held τ₀ as its limit.</item>
/// <item><see cref="HyperG"/> (Liang, Paulo, Molina, Clyde and Berger 2008): on g = 1/τ, p(g) = ((a − 2)/2)(1 + g)^(−a/2),
/// so g/(1 + g) ~ Beta(1, a/2 − 1). Stated on g itself, not on τ₀.</item>
/// <item><see cref="HyperGOverN"/>: the same law on g/n₀ = τ₀/τ, Liang et al.'s hyper-g/n under the g-form. Under the ridge,
/// where τ₀ = 1, it coincides with <see cref="HyperG"/>.</item>
/// <item><see cref="Held"/>: τ held at <see cref="Tau"/>, or at τ₀ when null: under the g-form unit information, 1/N.</item>
/// </list>
/// The hyper-g members are proper for a &gt; 2 only. Between 1 and 2 they are improper, their density stated up to a
/// constant: the posterior stays proper wherever every structure keeps a column inside the coefficient prior (the
/// intercept inside it), since τ's likelihood then vanishes as τ^(k/2) at τ = 0; and an engine refuses them wherever a
/// proper prior is needed (a b = 0 rung, SMC from the prior, the evidence). Empirical Bayes, which maximizes τ rather than
/// integrating it, is not a member.
/// <para>A value, in <see cref="Sigma2Prior"/>'s style: its kind and that kind's parameters, null where the kind takes
/// none, checked by the constructor, so one read back from JSON is checked as a built one is.</para>
/// </summary>
public sealed record TauPrior
{
    /// <summary>
    /// The prior of kind <paramref name="kind"/>: the Gamma takes a positive, finite shape; a hyper-g member its a, finite
    /// and above 1; the held value a positive, finite τ, or none for τ₀.
    /// </summary>
    public TauPrior(TauPriorKind kind, double? shape = null, double? tau = null)
    {
        switch (kind)
        {
            case TauPriorKind.UnitInformationGamma:
                if (shape is not (> 0.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(shape), "The Gamma's shape must be positive and finite.");
                if (tau is not null) throw new ArgumentException("The Gamma holds no τ.", nameof(tau));
                break;
            case TauPriorKind.HyperG or TauPriorKind.HyperGOverN:
                if (shape is not (> 1.0 and < double.PositiveInfinity))
                    throw new ArgumentOutOfRangeException(nameof(shape),
                        "Hyper-g's a must be finite and above 1; at or below 1 even a structure with a column inside the prior has no posterior for τ.");
                if (tau is not null) throw new ArgumentException("Hyper-g holds no τ.", nameof(tau));
                break;
            case TauPriorKind.Held:
                if (shape is not null) throw new ArgumentException("A held τ takes no shape.", nameof(shape));
                if (tau is { } t && !(t > 0.0 && double.IsFinite(t)))
                    throw new ArgumentOutOfRangeException(nameof(tau), "A held τ must be positive and finite, or null for τ₀.");
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(kind));
        }
        Kind = kind;
        Shape = shape;
        Tau = tau;
    }

    /// <summary>τ ~ Gamma(<paramref name="shape"/>, rate a/τ₀); ½, Zellner–Siow, by default.</summary>
    public static TauPrior UnitInformationGamma(double shape = 0.5) => new(TauPriorKind.UnitInformationGamma, shape);

    /// <summary>Hyper-g with <paramref name="a"/> on g = 1/τ; Liang et al. recommend a in (2, 4].</summary>
    public static TauPrior HyperG(double a) => new(TauPriorKind.HyperG, a);

    /// <summary>Hyper-g/n with <paramref name="a"/>: hyper-g's law on τ₀/τ.</summary>
    public static TauPrior HyperGOverN(double a) => new(TauPriorKind.HyperGOverN, a);

    /// <summary>τ held at <paramref name="tau"/>; null holds it at τ₀.</summary>
    public static TauPrior Held(double? tau = null) => new(TauPriorKind.Held, tau: tau);

    /// <summary>Which prior this is.</summary>
    public TauPriorKind Kind { get; }

    /// <summary>The Gamma's shape a, or a hyper-g member's a; null for a held τ.</summary>
    public double? Shape { get; }

    /// <summary>The held τ; null for τ₀, or for a sampled τ.</summary>
    public double? Tau { get; }

    /// <summary>Whether τ is held, so the state has no τ coordinate; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool IsHeld => Kind == TauPriorKind.Held;

    /// <summary>Whether the prior integrates to one: the Gamma, the held value, and hyper-g for a &gt; 2; read from the
    /// kind, so not serialized.</summary>
    [JsonIgnore]
    public bool IsProper => Kind is not (TauPriorKind.HyperG or TauPriorKind.HyperGOverN) || Shape > 2.0;

    /// <summary>
    /// The prior on a coefficient form whose unit-information precision is τ₀ = 1/<paramref name="unitObservations"/>:
    /// n₀ = N for the g-form, 1 for the normalized ridge.
    /// </summary>
    public ResolvedTauPrior For(double unitObservations) => new(this, unitObservations);

    /// <summary>The spelling a manifest or a log reads.</summary>
    public override string ToString() => Kind switch
    {
        TauPriorKind.Held => Tau is { } t ? FormattableString.Invariant($"Held({t:R})") : "Held(τ₀)",
        _ => FormattableString.Invariant($"{Kind}({Shape:R})"),
    };
}

/// <summary>
/// A <see cref="TauPrior"/> on its form's unit-information precision τ₀ = 1/n₀ (<see cref="TauPrior.For"/>): the
/// density of u = ln τ with its Jacobian, the exact draw of u, and the held value. The density is the one τ's slice step
/// reads on u (D04), so the step adds no Jacobian of its own.
/// </summary>
public sealed record ResolvedTauPrior
{
    private readonly double _logConstant;   // the Gamma's a ln(rate) − ln Γ(a); hyper-g's ln((a − 2)/2), or 0 when improper
    private readonly double _logScale;      // hyper-g's ln of the g its law is on, over 1/τ: 0, or ln τ₀ for hyper-g/n

    internal ResolvedTauPrior(TauPrior declared, double unitObservations)
    {
        ArgumentNullException.ThrowIfNull(declared);
        if (!(unitObservations > 0.0 && double.IsFinite(unitObservations)))
            throw new ArgumentOutOfRangeException(nameof(unitObservations), "n₀ must be positive and finite.");
        Declared = declared;
        UnitObservations = unitObservations;
        switch (declared.Kind)
        {
            case TauPriorKind.UnitInformationGamma:
                double a = declared.Shape!.Value;
                Gamma = Distributions.Gamma.ShapeRate(a, a * unitObservations);
                _logConstant = a * Math.Log(Gamma.Value.Rate) - SpecialFunctions.LogGamma(a);
                break;
            case TauPriorKind.HyperG or TauPriorKind.HyperGOverN:
                double h = declared.Shape!.Value;
                _logConstant = h > 2.0 ? Math.Log(0.5 * (h - 2.0)) : 0.0;
                _logScale = declared.Kind == TauPriorKind.HyperGOverN ? -Math.Log(unitObservations) : 0.0;
                break;
            default:
                HeldTau = declared.Tau ?? 1.0 / unitObservations;
                break;
        }
    }

    /// <summary>The prior as declared.</summary>
    public TauPrior Declared { get; }

    /// <summary>n₀, so τ₀ = 1/n₀.</summary>
    public double UnitObservations { get; }

    /// <summary>τ₀ = 1/n₀, the form's unit-information precision.</summary>
    public double UnitInformation => 1.0 / UnitObservations;

    /// <summary>Gamma(a, rate a·n₀) for the unit-information Gamma; null otherwise.</summary>
    public Gamma? Gamma { get; }

    /// <summary>The held τ, declared or τ₀; null when τ is sampled.</summary>
    public double? HeldTau { get; }

    /// <inheritdoc cref="TauPrior.IsHeld"/>
    public bool IsHeld => HeldTau is not null;

    /// <inheritdoc cref="TauPrior.IsProper"/>
    public bool IsProper => Declared.IsProper;

    /// <summary>
    /// The log density of u = ln τ, the Jacobian included; a sampled prior only.
    /// <list type="bullet">
    /// <item>The Gamma: a ln(rate) − ln Γ(a) + a·u − rate·eᵘ.</item>
    /// <item>Hyper-g, with v = ln g′ for g′ = g (hyper-g) or g/n₀ (hyper-g/n), g = e^(−u): ln((a − 2)/2) − (a/2)·ln(1 + eᵛ)
    /// + v, the constant dropped when a ≤ 2.</item>
    /// </list>
    /// Where τ leaves the doubles, rounding to 0 or overflowing, it is outside the support: −∞, a mass below e⁻³⁷⁰ of the
    /// prior's.
    /// </summary>
    public double LogDensity(double logTau)
    {
        if (IsHeld) throw new InvalidOperationException("A held τ has no density; the state has no τ coordinate.");
        if (double.IsNaN(logTau)) throw new ArgumentOutOfRangeException(nameof(logTau));
        double tau = Math.Exp(logTau);
        if (Gamma is Gamma g)
        {
            double rate = g.Rate * tau;
            return tau > 0.0 && double.IsFinite(rate) ? _logConstant + g.Shape * logTau - rate : double.NegativeInfinity;
        }
        if (!(tau > 0.0 && double.IsFinite(tau))) return double.NegativeInfinity;
        double v = _logScale - logTau;
        return _logConstant - 0.5 * Declared.Shape!.Value * Softplus(v) + v;
    }

    /// <summary>
    /// An exact draw of u = ln τ; a proper sampled prior only. The Gamma's draw is <see cref="Distributions.Gamma"/>'s; a
    /// hyper-g member's inverts g′/(1 + g′) ~ Beta(1, a/2 − 1): 1 − B = U^(1/β) with β = a/2 − 1, so u = ln(1 − B) − ln B
    /// − v₀, v₀ the offset of g′ from 1/τ. A draw whose τ leaves the doubles is drawn again.
    /// </summary>
    public double SampleLogTau(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (IsHeld) throw new InvalidOperationException("A held τ has no draw.");
        if (!IsProper) throw new InvalidOperationException($"{Declared} is improper, so it has no exact draw.");
        double u;
        if (Gamma is Gamma g)
        {
            do u = Math.Log(g.Sample(rng));
            while (double.IsNegativeInfinity(LogDensity(u)));
            return u;
        }
        double beta = 0.5 * Declared.Shape!.Value - 1.0;
        do
        {
            double logSurvival = Math.Log(1.0 - rng.NextDouble()) / beta;   // ln(1 − B), from U in (0, 1]
            u = logSurvival - Log1mExp(logSurvival) + _logScale;
        }
        while (double.IsNaN(u) || double.IsInfinity(u) || double.IsNegativeInfinity(LogDensity(u)));
        return u;
    }

    /// <summary>The spelling a manifest or a log reads: the declared prior with its n₀.</summary>
    public override string ToString() => FormattableString.Invariant($"{Declared} on τ₀ = 1/{UnitObservations:R}");

    private static double Softplus(double x) => LogFunctions.Softplus(x);

    private static double Log1mExp(double u) => LogFunctions.Log1mExp(u);
}

/// <summary>The log-domain functions the g-prior blocks read, accurate in their tails.</summary>
internal static class LogFunctions
{
    /// <summary>ln(1 + x), accurate for small x (Goldberg's correction).</summary>
    public static double Log1p(double x)
    {
        double u = 1.0 + x;
        return u == 1.0 ? x : double.IsPositiveInfinity(u) ? u : Math.Log(u) * x / (u - 1.0);
    }

    /// <summary>ln(1 + eˣ), without overflow and accurate in both tails.</summary>
    public static double Softplus(double x) => x > 0.0 ? x + Log1p(Math.Exp(-x)) : Log1p(Math.Exp(x));

    /// <summary>ln(1 − eᵘ) for u ≤ 0 (Mächler 2012).</summary>
    public static double Log1mExp(double u)
    {
        if (u > -0.6931471805599453)
        {
            double e = Math.Exp(u);   // eᵘ − 1 by Kahan's correction
            double expm1 = e == 1.0 ? u : (e - 1.0) * u / Math.Log(e);
            return Math.Log(-expm1);
        }
        return Log1p(-Math.Exp(u));
    }
}
