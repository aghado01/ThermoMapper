using Maths.LinAlg;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// The numerical gate every BAPS solve applies to its factored system, <c>ZᵀZ + λP</c> and its kin: the strict factor is
/// positive definite and its estimated equilibrated condition number (<see cref="ISpdFactor.EquilibratedCondition"/>) is at
/// most <see cref="MaxCondition"/>. Both are scale-free, so a valid fit in small units is admitted. A refused system throws
/// where it is solved.
/// </summary>
/// <remarks>
/// At κ(E) = 10¹⁴ a backward-stable solve still keeps about two digits (κ·ε ≈ 0.02). Measured on the nine tensor layouts of
/// <c>issues/linalg/spectral-tensor-factorization-report.md</c>:
/// <list type="bullet">
/// <item>every identified system in the optimizer's boxes reads κ(E) ≤ 7.8·10⁹, and the public λ = 10¹² points
/// 1.3·10¹²–2.5·10¹³;</item>
/// <item>every system whose design leaves the penalty's null space unidentified, and every λ = 10²⁰ point, reads at least
/// 1.2·10¹⁶.</item>
/// </list>
/// The gate it replaces, a floored factor with no floored pivot, refused a valid fit in small units and admitted 55 of 234
/// unidentified systems. BARS's <c>DesignAdmission</c> bounds κ(E) at 10¹⁰ for a sampler that scores a refusal −∞; BAPS
/// throws on one, and its boxes reach 7.8·10⁹.
/// </remarks>
internal static class BapsAdmission
{
    /// <summary>The largest admissible estimated condition number of the equilibrated system.</summary>
    public const double MaxCondition = 1e14;

    /// <summary>Whether a factored system is admissible: positive definite, with κ(E) at most <see cref="MaxCondition"/>.</summary>
    public static bool Admits(ISpdFactor factor) =>
        factor.Status.PositiveDefinite && factor.EquilibratedCondition() <= MaxCondition;
}
