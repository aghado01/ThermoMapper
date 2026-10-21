using System;
using Maths.Rng;

using Maths.Regression.GPrior;
using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Exact marginal likelihood for the Normal model under the unit-information prior (DiMatteo, Genovese &amp;
/// Kass 2001, eq. 6; analytic form in He, Yang &amp; Kang 2024, Lemma 1): with β, σ integrated out,
/// <c>p(y|Z) ∝ (m+1)^(−ν/2) · a^(−m/2)</c> where
/// <c>a = yᵀy − (m/(m+1)) · (Zᵀy)ᵀ (ZᵀZ)⁻¹ (Zᵀy)</c>. Returned up to the σ-integration constant that depends only
/// on m and so cancels in ratios.
/// <para>With per-point weights W (robust scale-mixture BARS, <c>y ~ N(Zβ, σ²W⁻¹)</c>) the coefficient prior stays the
/// <i>unweighted</i> g-prior <c>β ~ N(0, mσ²(ZᵀZ)⁻¹)</c>, so the prior on the curve does not depend on the weights
/// (decision Q2). The marginal is then
/// <c>½log|ZᵀZ/m| − ½log|Zᵀ(W + I/m)Z| + ½Σ log wᵢ − (m/2)·log S</c>, with
/// <c>S = yᵀWy − hᵀ(Zᵀ(W + I/m)Z)⁻¹h</c> and <c>h = ZᵀWy</c>; it reduces to the form above at unit weights.</para>
/// </summary>
/// <remarks>
/// Built on the collapsed g-prior block (<see cref="CollapsedGPrior"/>) at τ = 1/m and inverse temperature 1, with
/// <see cref="BandedGPrior"/> factoring the banded Gram matrices (a root-free LDLᵀ; no inverse is formed). The block
/// computes the scale from the explicit residual, <c>‖y − Zβ̂‖² + Q/(m+1)</c>, the cancellation-free form of
/// <c>a</c>. This model adds the m ≤ ν rule and the admission policy, and keeps its marginal up to the
/// σ-integration constant.
/// </remarks>
public sealed class WeightedNormalModel : IObservationModel
{
    /// <inheritdoc/>
    public int EffectiveDimension(double[,] design, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        return design.GetLength(1);
    }

    /// <inheritdoc/>
    public CoefficientUncertainty Uncertainty => CoefficientUncertainty.Exact;

    /// <inheritdoc/>
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);

        // With at most as many observations as parameters the unit-information marginal is improper —
        // an unfittable (over-knotted) configuration. Report it as infinitely unlikely so the chain rejects.
        if (design.GetLength(0) <= design.GetLength(1))
            return double.NegativeInfinity;

        (GPriorConditional c, bool admitted) = Condition(design, y, weights);
        if (!admitted)
            return double.NegativeInfinity;   // rank-deficient or ill-conditioned design ⇒ outside the support
        if (c.S <= 0.0)
            return double.NegativeInfinity;   // an identically zero response: no posterior for σ

        return c.LogDeterminantTerm - 0.5 * c.Rows * Math.Log(c.S);
    }

    /// <summary>
    /// A draw from the exact conditional posterior given the knots: <c>σ² ~ IG(m/2, S/2)</c>, then β from its Normal
    /// conditional through the band factor's precision draw. At unit weights <c>S = a</c> and
    /// <c>β ~ N((m/(m+1))(ZᵀZ)⁻¹Zᵀy, σ²·(m/(m+1))·(ZᵀZ)⁻¹)</c>; with weights, <c>β ~ N(Pa⁻¹h, σ²Pa⁻¹)</c> with
    /// <c>Pa = Zᵀ(W + I/m)Z</c>. Falls back to the mean when the design is inadmissible or the response is zero.
    /// </summary>
    public double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => SampleConditional(design, y, weights, rng).Beta;

    /// <summary>The joint draw behind <see cref="SampleCoefficients"/>: the coefficients and the noise variance they
    /// were drawn with. σ² is NaN when no draw is possible. The robust weight update conditions on both.</summary>
    internal (double[] Beta, double Sigma2) SampleConditional(double[,] design, double[] y, double[]? weights,
                                                             Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (GPriorConditional c, bool admitted) = Condition(design, y, weights);
        double[] beta = c.Mean.ToArray();
        if (!admitted || !(c.S > 0.0))
            return (beta, double.NaN);
        double sigma2 = c.Sample(null, rng, beta);   // the improper 1/σ² prior
        return (beta, sigma2);
    }

    /// <summary>
    /// <c>Var[f(t) | knots] = E[σ²]·z(t)ᵀ Cov z(t)</c> at each readout row z(t), with σ integrated out (the Student-t
    /// predictive's variance): <c>(a/(m−2))·(m/(m+1))·z(t)ᵀ(ZᵀZ)⁻¹z(t)</c> at unit weights, and
    /// <c>(S/(m−2))·z(t)ᵀPa⁻¹z(t)</c> with weights. NaN entries when the design is inadmissible or m ≤ 2.
    /// </summary>
    public double[]? ConditionalVariance(double[,] design, double[] y, double[]? weights, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        (GPriorConditional c, bool admitted) = Condition(design, y, weights);
        if (!admitted || !(c.S > 0.0) || c.Rows <= 2)
            return ReadoutQuadratics.Unavailable(readoutDesign);
        return c.ConditionalVariance(null, readoutDesign);
    }

    /// <summary>
    /// Posterior mean spline coefficients: <c>(m/(m+1)) (ZᵀZ)⁻¹ Zᵀy</c> at unit weights (DMGK eq. 7's
    /// <c>E{f | k, ξ, y}</c>), and <c>Pa⁻¹h</c> with weights. Evaluated against a design at any points this is the
    /// conditional posterior mean curve for one knot draw; averaging over the chain is the Bayes fit.
    /// <paramref name="weights"/> null means unit weights.
    /// </summary>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
        => Condition(design, y, weights).C.Mean.ToArray();

    // One knot configuration's conditional posterior at BARS's τ = 1/m, untempered, and the admission verdict: G's
    // factor alone without weights; G's and P_a's with them.
    private static (GPriorConditional C, bool Admitted) Condition(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        int m = design.GetLength(0);
        if (y.Length != m)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (weights is not null && weights.Length != m)
            throw new ArgumentException("Weights length must match the design row count.", nameof(weights));

        var banded = new BandedDesign(design);
        double tau = 1.0 / m;
        if (weights is null)
        {
            CollapsedGPrior block = BandedGPrior.Fit(banded, y);
            return (block.At(tau), DesignAdmission.Admits(block.GramStatus));
        }
        CollapsedGPrior weighted = BandedGPrior.Fit(banded, y, weights, tau);
        return (weighted.At(tau), DesignAdmission.Admits(weighted.GramStatus) && DesignAdmission.Admits(weighted.PosteriorStatus));
    }
}
