using System;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A polynomial-chaos basis (E19): every term a product of standardized shifted Legendre polynomials over the continuous
/// predictors' scaled values, from KHAOS's term space with order at most <see cref="BassStructureSpec.MaxInteraction"/> (q_max)
/// and degree at most <see cref="MaxDegree"/> (d_max), and KHAOS's proposal over it (E4, F-X4):
/// <list type="bullet">
/// <item><b>The birth,</b> by Bernoulli inclusion: q₀ from {1, …, q_max} ∝ q₀^(−s_q) (<see cref="OrderPenalty"/>); each
/// variable's inclusion probability ηⱼ from the current terms' appearance counts by KHAOS's <c>make_weights</c>
/// (<see cref="InclusionExponent"/>, <see cref="InclusionPseudoCount"/>, <see cref="InclusionPasses"/>), summing to about q₀;
/// the inclusions redrawn while none or more than q_max come out; the degree from {q, …, d_max} ∝ (d − q + 1)^(−s_d)
/// (<see cref="DegreePenalty"/>); and a uniform composition of it over the variables. Its forward probability is exact.</item>
/// <item><b>The change,</b> one of KHAOS's two mutations: the degree redrawn with a new composition, with probability
/// <see cref="DegreeMutation"/> (always when p ≤ 3, KHAOS's rule); otherwise one variable swapped for one the term does not
/// use, ∝ <see cref="SwapPseudoCount"/> plus its appearance count.</item>
/// </list>
/// The defaults are <c>khaos</c> 2.1.1's, the mutation's adaptive probability fixed at its middle (the brief's fifth
/// decision). D01's hinge settings (the hinge degree, w₁, w₂, the local steps, the fixed-M family) do not apply.
/// </summary>
public sealed record LegendreBasis
{
    /// <summary>d_max, the largest degree a term's factors sum to; at least q_max. <c>khaos</c>'s default 15.</summary>
    public int MaxDegree { get; init; } = 15;

    /// <summary>s_q: the birth draws q₀ ∝ q₀^(−s_q); <c>khaos</c>'s J(q₀) = 1/q₀.</summary>
    public double OrderPenalty { get; init; } = 1.0;

    /// <summary>s_d: the birth and the degree mutation draw the degree d ∝ (d − q + 1)^(−s_d); 0, uniform, by default.</summary>
    public double DegreePenalty { get; init; }

    /// <summary><c>make_weights</c>' α, the power on each variable's appearance count.</summary>
    public double InclusionExponent { get; init; } = 1.0;

    /// <summary><c>make_weights</c>' ε, the pseudo-count added to each variable's powered count; positive, so every variable can be
    /// included and every death has a reverse birth.</summary>
    public double InclusionPseudoCount { get; init; } = 2.0;

    /// <summary><c>make_weights</c>' passes, the fixed-point steps that bring the inclusion probabilities' sum toward q₀.</summary>
    public int InclusionPasses { get; init; } = 3;

    /// <summary>The variable swap's pseudo-count: a new variable is drawn ∝ this plus its appearance count; positive.</summary>
    public double SwapPseudoCount { get; init; } = 1.0;

    /// <summary>The probability a change mutates the degree rather than swapping a variable, when p &gt; 3.</summary>
    public double DegreeMutation { get; init; } = 0.5;

    /// <summary>The first setting out of range, as its name and the reason; null when every one holds.</summary>
    public (string Setting, string Reason)? Problem()
    {
        if (MaxDegree < 1) return (nameof(MaxDegree), "must be a positive integer");
        if (!(OrderPenalty >= 0.0 && double.IsFinite(OrderPenalty))) return (nameof(OrderPenalty), "must be finite and not negative");
        if (!(DegreePenalty >= 0.0 && double.IsFinite(DegreePenalty))) return (nameof(DegreePenalty), "must be finite and not negative");
        if (!(InclusionExponent >= 0.0 && double.IsFinite(InclusionExponent))) return (nameof(InclusionExponent), "must be finite and not negative");
        if (!(InclusionPseudoCount > 0.0 && double.IsFinite(InclusionPseudoCount)))
            return (nameof(InclusionPseudoCount), "must be positive and finite, or a death could have no reverse birth");
        if (InclusionPasses < 1) return (nameof(InclusionPasses), "must be at least 1");
        if (!(SwapPseudoCount > 0.0 && double.IsFinite(SwapPseudoCount)))
            return (nameof(SwapPseudoCount), "must be positive and finite, or a swap could have no reverse");
        if (!(DegreeMutation >= 0.0 && DegreeMutation <= 1.0)) return (nameof(DegreeMutation), "must lie in [0, 1]");
        return null;
    }
}
