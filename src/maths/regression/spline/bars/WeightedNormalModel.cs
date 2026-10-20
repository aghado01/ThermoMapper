using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

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
    public CoefficientUncertainty Uncertainty => CoefficientUncertainty.Exact;

    /// <inheritdoc/>
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);

        // With at most as many observations as parameters the unit-information marginal is improper —
        // an unfittable (over-knotted) configuration. Report it as infinitely unlikely so the chain rejects.
        if (design.GetLength(0) <= design.GetLength(1))
            return double.NegativeInfinity;

        Conditioned c = Condition(design, y, weights);
        if (!c.Admitted)
            return double.NegativeInfinity;   // rank-deficient or ill-conditioned design ⇒ outside the support
        if (c.S <= 0.0)
            return double.NegativeInfinity;   // an identically zero response: no posterior for σ

        return c.LogDeterminantTerm - 0.5 * c.M * Math.Log(c.S);
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
        Conditioned c = Condition(design, y, weights);
        double[] beta = (double[])c.MeanCoefficients.Clone();
        if (!c.Admitted || !(c.S > 0.0))
            return (beta, double.NaN);

        double sigma2 = InverseGamma.ShapeScale(0.5 * c.M, 0.5 * c.S).Sample(rng);
        var z = new double[c.Nu];
        for (int j = 0; j < c.Nu; j++) z[j] = Normal.SampleStandard(rng);
        var innov = new double[c.Nu];
        c.Factor.SampleInnovation(z, innov);                               // ~ N(0, factored matrix⁻¹)
        double scale = Math.Sqrt(sigma2 * c.FactorScale);
        for (int j = 0; j < c.Nu; j++) beta[j] += scale * innov[j];
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
        Conditioned c = Condition(design, y, weights);
        if (!c.Admitted || !(c.S > 0.0) || c.M <= 2)
            return ReadoutQuadratics.Unavailable(readoutDesign);
        return ReadoutQuadratics.Of(c.Factor, readoutDesign, c.S / (c.M - 2) * c.FactorScale);
    }

    /// <summary>
    /// Posterior mean spline coefficients: <c>(m/(m+1)) (ZᵀZ)⁻¹ Zᵀy</c> at unit weights (DMGK eq. 7's
    /// <c>E{f | k, ξ, y}</c>), and <c>Pa⁻¹h</c> with weights. Evaluated against a design at any points this is the
    /// conditional posterior mean curve for one knot draw; averaging over the chain is the Bayes fit.
    /// <paramref name="weights"/> null means unit weights.
    /// </summary>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
        => (double[])Condition(design, y, weights).MeanCoefficients.Clone();

    // One knot configuration's conditional posterior pieces: a factor of the coefficients' precision matrix (the
    // covariance is σ²·FactorScale·factored⁻¹), the posterior mean, S, the marginal's determinant term and the
    // admission verdict. Null weights factor G = ZᵀZ once; weights factor G for the prior and Pa = Zᵀ(W + I/m)Z.
    private readonly record struct Conditioned(
        BandCholesky Factor, double[] MeanCoefficients, double S, double FactorScale, double LogDeterminantTerm,
        int M, int Nu, bool Admitted);

    private static Conditioned Condition(double[,] design, double[] y, double[]? weights)
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
        double shrink = (double)m / (m + 1);

        if (weights is null)
        {
            var band = new double[design2.Bandwidth + 1, nu];
            var b = new double[nu];
            design2.Accumulate(null, band, y, b);
            double yy = 0.0;
            for (int i = 0; i < m; i++) yy += y[i] * y[i];

            var chol = new BandCholesky(nu, design2.Bandwidth, BandFactorization.Ldlt);
            chol.DecomposeBanded(band);
            double[] sol = chol.Solve(b);   // G⁻¹b — no inverse materialized
            double quad = 0.0;              // bᵀG⁻¹b, the regression sum of squares
            for (int p = 0; p < nu; p++) quad += b[p] * sol[p];

            // a = ‖y − Zβ̂‖² + quad/(m+1) ≥ quad/(m+1): positive for any nonzero response, an interpolating fit
            // included (over-knotting is the m ≤ ν rule). The subtraction loses at most log₁₀(m+1) digits.
            double a = yy - shrink * quad;
            var mean = new double[nu];
            for (int p = 0; p < nu; p++) mean[p] = shrink * sol[p];
            return new Conditioned(chol, mean, a, shrink, -0.5 * nu * Math.Log(m + 1), m, nu,
                                   DesignAdmission.Admits(chol.Status));
        }

        // Weighted: the prior uses the unweighted Gram G = ZᵀZ; the posterior precision is Pa = Zᵀ(W + I/m)Z.
        var gBand = new double[design2.Bandwidth + 1, nu];
        design2.Accumulate(null, gBand, null, null);
        var gChol = new BandCholesky(nu, design2.Bandwidth, BandFactorization.Ldlt);
        gChol.DecomposeBanded(gBand);

        var wPlus = new double[m];
        double yWy = 0.0, sumLogW = 0.0;
        for (int i = 0; i < m; i++)
        {
            wPlus[i] = weights[i] + 1.0 / m;
            yWy += weights[i] * y[i] * y[i];
            sumLogW += Math.Log(weights[i]);
        }
        var pBand = new double[design2.Bandwidth + 1, nu];
        design2.Accumulate(wPlus, pBand, null, null);
        var h = new double[nu];
        design2.AccumulateRhs(weights, y, h);
        var pChol = new BandCholesky(nu, design2.Bandwidth, BandFactorization.Ldlt);
        pChol.DecomposeBanded(pBand);

        double[] beta = pChol.Solve(h);   // Pa⁻¹h, the posterior mean
        double quadW = 0.0;
        for (int p = 0; p < nu; p++) quadW += h[p] * beta[p];

        // S = ‖y − Zβ̃‖²_W + β̃ᵀ(G/m)β̃ ≥ 0, zero only for an identically zero response.
        double s = yWy - quadW;
        double logDet = 0.5 * (gChol.LogDet - nu * Math.Log(m)) - 0.5 * pChol.LogDet + 0.5 * sumLogW;
        bool admitted = DesignAdmission.Admits(gChol.Status) && DesignAdmission.Admits(pChol.Status);
        return new Conditioned(pChol, beta, s, 1.0, logDet, m, nu, admitted);
    }
}
