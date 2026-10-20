using Maths.LinAlg;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The admission rule every BARS observation model applies to its factored Gram or Hessian (decision Q4 of the BASS
/// preparation brief): a knot configuration is in the model's support only when the factor is exact (no pivot
/// floored), positive definite, and not proven ill-conditioned. A rejected configuration scores −∞.
/// </summary>
/// <remarks>
/// <para><see cref="MinPivotRatio"/> = 1e-10. Since <see cref="FactorStatus.PivotRatio"/> ≥ 1/κ for a positive-definite
/// matrix, a smaller ratio proves κ &gt; 1e10, beyond which a solve keeps fewer than about 6 of 16 digits.</para>
/// <para>For the Normal and Whittle models the factored matrix depends only on the inputs, so the rule is a declared
/// restriction of the model's support: every basis function needs data under it. Poisson and GibbsLoss factor a
/// matrix weighted by the fit, so for them admission also depends on the response, and the rule is a numerical
/// truncation of the posterior rather than a prior restriction.</para>
/// </remarks>
public static class DesignAdmission
{
    /// <summary>The smallest admissible <see cref="FactorStatus.PivotRatio"/>.</summary>
    public const double MinPivotRatio = 1e-10;

    /// <summary>Whether a factor with this status is admissible.</summary>
    public static bool Admits(in FactorStatus status) =>
        status.PositiveDefinite && status.FlooredPivots == 0 && status.PivotRatio >= MinPivotRatio;
}
