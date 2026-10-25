namespace Maths.Regression.Spline.Baps;

/// <summary>
/// The prior on the noise variance σ² in <see cref="Baps"/>'s Gibbs schemes and <see cref="HierarchicalPSpline"/>. Either
/// keeps σ²'s conditional conjugate: m residual terms with sum of squares S give IG(a + m/2, b + S/2), or IG(m/2, S/2).
/// </summary>
public enum NoiseVariancePrior
{
    /// <summary>IG(priorShape, priorScale), applied in the response's own units, so b = 10⁻³ is small or large depending
    /// on them: below σ² ≈ b it dominates the data. The default of BAPS's isotropic Gibbs schemes and of
    /// <see cref="HierarchicalPSpline"/>.</summary>
    InverseGamma,

    /// <summary>The Jeffreys prior p(σ²) ∝ 1/σ²: scale-free, and proper in the posterior whenever the residual sum of
    /// squares is positive. The collapsed scheme's evidence assumes it, Bach &amp; Klein (2025) use it, and it is the
    /// anisotropic scheme's default.</summary>
    Jeffreys,
}
