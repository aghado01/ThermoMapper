using System;
using Maths.LinAlg;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Exact marginal likelihood for the Normal model under the unit-information prior (DiMatteo, Genovese &amp;
/// Kass 2001, eq. 6; analytic form in He, Yang &amp; Kang 2024, Lemma 1): with β, σ integrated out,
/// <c>p(y|Z) ∝ (m+1)^(−ν/2) · a^(−m/2)</c> where
/// <c>a = yᵀW y − (m/(m+1)) · (ZᵀW y)ᵀ (ZᵀW Z)⁻¹ (ZᵀW y)</c>. Conditional on the per-point weights W this is
/// the inner marginal of robust scale-mixture BARS; <c>W ≡ I</c> (null weights) is the non-robust case. Reduces
/// the design + response to the log-marginal currency the chain compares; returned up to the σ-integration
/// constant that depends only on m and so cancels in ratios.
/// </summary>
/// <remarks>
/// <c>A = ZᵀW Z</c> is banded at the spline degree (heptadiagonal for cubics, diagonal for step bases), so it is
/// factored by a root-free banded LDLᵀ (<see cref="BandCholesky"/>) and the regression sum of squares
/// <c>bᵀA⁻¹b</c> and posterior coefficients <c>A⁻¹b</c> are read from the band solve alone — A⁻¹ is never formed.
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
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);

        int m = design.GetLength(0);
        int nu = design.GetLength(1);
        if (y.Length != m)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (weights is not null && weights.Length != m)
            throw new ArgumentException("Weights length must match the design row count.", nameof(weights));

        // With at most as many observations as parameters the unit-information marginal is improper —
        // an unfittable (over-knotted) configuration. Report it as infinitely unlikely so the chain rejects.
        if (m <= nu)
            return double.NegativeInfinity;

        var design2 = new BandedDesign(design);
        var band = new double[design2.Bandwidth + 1, nu];
        var b = new double[nu];
        design2.Accumulate(weights, band, y, b);

        double yWy = 0.0;
        for (int i = 0; i < m; i++)
        {
            double w = weights is null ? 1.0 : weights[i];
            yWy += w * y[i] * y[i];
        }

        var chol = new BandCholesky(nu, design2.Bandwidth, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        if (chol.HitFloor)
            return double.NegativeInfinity;   // singular / rank-deficient design ⇒ unfittable, reject

        double[] sol = chol.Solve(b);   // A⁻¹b — no inverse materialized

        // quad = bᵀ A⁻¹ b — the weighted regression sum of squares.
        double quad = 0.0;
        for (int p = 0; p < nu; p++)
            quad += b[p] * sol[p];

        double shrink = (double)m / (m + 1);
        double aScalar = yWy - shrink * quad;
        if (aScalar <= 0.0)
            return double.NegativeInfinity;   // over-interpolating fit ⇒ reject (never reward a singular fit)

        return -0.5 * nu * Math.Log(m + 1) - 0.5 * m * Math.Log(aScalar);
    }

    /// <summary>
    /// Posterior mean spline coefficients under the unit-information prior:
    /// <c>ĉ = (m/(m+1)) (ZᵀW Z)⁻¹ ZᵀW y</c> (DMGK eq. 7's <c>E{f | k, ξ, y}</c>). Evaluated against a design at
    /// any points this is the conditional posterior mean curve for one knot draw; averaging over the chain is
    /// the Bayes fit. <paramref name="weights"/> null means unit weights.
    /// </summary>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        int m = design.GetLength(0);
        int nu = design.GetLength(1);
        if (y.Length != m)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (weights is not null && weights.Length != m)
            throw new ArgumentException("Weights length must match the design row count.", nameof(weights));

        var design2 = new BandedDesign(design);
        var band = new double[design2.Bandwidth + 1, nu];
        var b = new double[nu];
        design2.Accumulate(weights, band, y, b);

        var chol = new BandCholesky(nu, design2.Bandwidth, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        double[] sol = chol.Solve(b);   // A⁻¹b

        double shrink = (double)m / (m + 1);
        var c = new double[nu];
        for (int p = 0; p < nu; p++)
            c[p] = shrink * sol[p];
        return c;
    }
}
