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
        if (c.A <= 0.0)
            return double.NegativeInfinity;   // an identically zero response: no posterior for σ

        return -0.5 * c.Nu * Math.Log(c.M + 1) - 0.5 * c.M * Math.Log(c.A);
    }

    /// <summary>
    /// A draw from the exact conditional posterior given the knots: <c>σ² ~ IG(m/2, a/2)</c>, then
    /// <c>β ~ N((m/(m+1))A⁻¹b, σ²·(m/(m+1))·A⁻¹)</c> through the band factor's precision draw. Falls back to the
    /// mean when the design is inadmissible or the response is identically zero.
    /// </summary>
    public double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Conditioned c = Condition(design, y, weights);
        double[] beta = c.Mean();
        if (!c.Admitted || !(c.A > 0.0))
            return beta;

        double sigma2 = 1.0 / Gamma.Sample(rng, 0.5 * c.M, 2.0 / c.A);   // IG(m/2, a/2) as 1/Gamma(shape, scale)
        var z = new double[c.Nu];
        for (int j = 0; j < c.Nu; j++) z[j] = Normal.Sample(rng);
        var innov = new double[c.Nu];
        c.Factor.SampleInnovation(z, innov);                               // ~ N(0, A⁻¹)
        double scale = Math.Sqrt(sigma2 * c.Shrink);
        for (int j = 0; j < c.Nu; j++) beta[j] += scale * innov[j];
        return beta;
    }

    /// <summary>
    /// <c>Var[f(t) | knots] = (a/(m−2))·(m/(m+1))·z(t)ᵀA⁻¹z(t)</c> at each readout row z(t), with σ integrated out
    /// (the Student-t predictive's variance). NaN entries when the design is inadmissible or m ≤ 2.
    /// </summary>
    public double[]? ConditionalVariance(double[,] design, double[] y, double[]? weights, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        Conditioned c = Condition(design, y, weights);
        if (!c.Admitted || !(c.A > 0.0) || c.M <= 2)
            return ReadoutQuadratics.Unavailable(readoutDesign);
        return ReadoutQuadratics.Of(c.Factor, readoutDesign, c.A / (c.M - 2) * c.Shrink);
    }

    /// <summary>
    /// Posterior mean spline coefficients under the unit-information prior:
    /// <c>ĉ = (m/(m+1)) (ZᵀW Z)⁻¹ ZᵀW y</c> (DMGK eq. 7's <c>E{f | k, ξ, y}</c>). Evaluated against a design at
    /// any points this is the conditional posterior mean curve for one knot draw; averaging over the chain is
    /// the Bayes fit. <paramref name="weights"/> null means unit weights.
    /// </summary>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
        => Condition(design, y, weights).Mean();

    // One knot configuration's conditional posterior pieces, from a single banded factorization of A = ZᵀWZ:
    // the factor, A⁻¹b, a = yᵀWy − (m/(m+1))·bᵀA⁻¹b, and the admission verdict.
    private readonly record struct Conditioned(
        BandCholesky Factor, double[] Solution, double A, int M, int Nu, bool Admitted)
    {
        public double Shrink => (double)M / (M + 1);

        public double[] Mean()
        {
            var c = new double[Nu];
            for (int p = 0; p < Nu; p++) c[p] = Shrink * Solution[p];
            return c;
        }
    }

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
        double[] sol = chol.Solve(b);   // A⁻¹b — no inverse materialized

        // quad = bᵀ A⁻¹ b — the weighted regression sum of squares.
        double quad = 0.0;
        for (int p = 0; p < nu; p++)
            quad += b[p] * sol[p];

        // a = ‖y − Zβ̂‖²_W + quad/(m+1) ≥ quad/(m+1): positive for any nonzero response, an interpolating fit
        // included (over-knotting is the m ≤ ν rule). The subtraction loses at most log₁₀(m+1) digits.
        double aScalar = yWy - (double)m / (m + 1) * quad;
        return new Conditioned(chol, sol, aScalar, m, nu, DesignAdmission.Admits(chol.Status));
    }
}
