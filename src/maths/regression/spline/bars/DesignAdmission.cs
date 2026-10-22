using Maths.LinAlg;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The numerical gate every BARS observation model applies to its factored Gram or Hessian (Q6 of the BASS preparation
/// brief, step 3.7): the system can be trusted only when its strict factor is positive definite and its smallest
/// equilibrated pivot is at least <see cref="MinRelativePivot"/>. A refused configuration scores −∞. Whether each basis
/// function has data under it is the support gate's question, which <see cref="SplineTarget"/> asks first.
/// </summary>
/// <remarks>
/// <para><see cref="MinRelativePivot"/> = 1e-10 keeps Q4's six-digit rationale on the equilibrated matrix
/// D^(−½)·A·D^(−½), D = diag(A), whose condition number governs Cholesky's accuracy (van der Sluis 1969; Demmel
/// 1989): 1/<see cref="FactorStatus.MinRelativePivot"/> bounds it from below, so a smaller pivot proves more than about
/// ten of sixteen digits lost. Unlike the pivot ratio it replaced, it does not change when a column is rescaled.</para>
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

    /// <summary>Whether a factor with this status is admissible.</summary>
    public static bool Admits(in FactorStatus status) =>
        status.PositiveDefinite && status.FlooredPivots == 0 && status.MinRelativePivot >= MinRelativePivot;
}
