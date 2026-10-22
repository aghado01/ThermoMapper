using Maths.LinAlg;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The numerical gate every BARS observation model applies to its factored Gram or Hessian (Q6 of the BASS preparation
/// brief, step 3.7): the system can be trusted only when its strict factor is positive definite, its smallest
/// equilibrated pivot is at least <see cref="MinRelativePivot"/>, and its estimated equilibrated condition number is at
/// most <see cref="MaxCondition"/>. A refused configuration scores −∞. Whether each basis function has data under it is
/// the support gate's question, which <see cref="SplineTarget"/> asks first.
/// </summary>
/// <remarks>
/// <para>Both thresholds are on the equilibrated matrix E = D^(−½)·A·D^(−½), D = diag(A), whose condition number
/// governs Cholesky's accuracy (van der Sluis 1969; Demmel 1989); above 10¹⁰ a solve may lose ten of sixteen digits.
/// 1/<see cref="FactorStatus.MinRelativePivot"/> only bounds κ(E) from below, so a small pivot proves a bad matrix but a
/// healthy one proves nothing: a bidiagonal-factor matrix has every relative pivot at 0.2 and κ(E) near 10¹⁵. The
/// condition estimate (<see cref="ISpdFactor.EquilibratedCondition"/>) measures κ(E); the pivot check stays as the
/// free first test, and as the per-step test inside an iterative fit (<see cref="PivotsAdmit"/>). Neither changes when
/// a column is rescaled.</para>
/// <para>BARS's factors are strict, so no pivot is ever floored; a floored factor's solves would describe a different
/// matrix, and <see cref="FactorStatus.FlooredPivots"/> must be 0.</para>
/// <para>For the Normal and Whittle models the factored matrix depends only on the inputs, so the gate restricts the
/// model's support. Poisson and GibbsLoss factor a matrix weighted by the fit, so for them it also depends on the
/// response, and it is a numerical truncation of the posterior rather than a prior restriction.</para>
/// </remarks>
public static class DesignAdmission
{
    /// <summary>The smallest admissible <see cref="FactorStatus.MinRelativePivot"/>.</summary>
    public const double MinRelativePivot = 1e-10;

    /// <summary>The largest admissible estimated condition number of the equilibrated matrix.</summary>
    public const double MaxCondition = 1e10;

    /// <summary>Whether the pivots allow admission: positive definite, none floored, and none below
    /// <see cref="MinRelativePivot"/>. Necessary, not sufficient: an iterative fit checks it at each step before solving,
    /// and <see cref="Admits"/> at the end.</summary>
    public static bool PivotsAdmit(in FactorStatus status) =>
        status.PositiveDefinite && status.FlooredPivots == 0 && status.MinRelativePivot >= MinRelativePivot;

    /// <summary>Whether a factor is admissible: its pivots allow it and its estimated equilibrated condition number is at
    /// most <see cref="MaxCondition"/>.</summary>
    public static bool Admits(ISpdFactor factor) =>
        PivotsAdmit(factor.Status) && factor.EquilibratedCondition() <= MaxCondition;
}
