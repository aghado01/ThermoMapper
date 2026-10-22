using System;
using Maths.LinAlg;
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
    public bool CanTemperLikelihood => true;

    /// <inheritdoc/>
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
        => TemperedMarginal(design, y, weights).LogMarginalLikelihood(1.0);

    /// <summary>
    /// The marginal at inverse temperature b, from the block at τ = 1/m: <c>LDT_b − (mb/2)·ln S_b</c>, dropping
    /// the σ-integration constant, which depends on m and b only. Unweighted, any b in (0, 1]; b = 0 is refused,
    /// since the improper 1/σ² prior then leaves σ unidentified. Weighted, b = 1 only (the block is built at one τ
    /// and b). A configuration with m ≤ ν, an inadmissible design or a zero response is −∞ at every b.
    /// </summary>
    public ITemperedMarginal TemperedMarginal(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);

        // With at most as many observations as parameters the unit-information marginal is improper —
        // an unfittable (over-knotted) configuration. Report it as infinitely unlikely so the chain rejects.
        if (design.GetLength(0) <= design.GetLength(1))
            return new ImpossibleMarginal(DesignRejection.TooFewObservations);

        (CollapsedGPrior? block, double tau) = Fit(design, y, weights);
        if (block is null)   // rank-deficient or ill-conditioned design ⇒ outside the support
            return new ImpossibleMarginal(DesignRejection.Inadmissible);
        return new NormalMarginal(block, tau);
    }

    /// <summary>
    /// A draw from the exact conditional posterior given the knots: <c>σ² ~ IG(m/2, S/2)</c>, then β from its Normal
    /// conditional through the band factor's precision draw. At unit weights <c>S = a</c> and
    /// <c>β ~ N((m/(m+1))(ZᵀZ)⁻¹Zᵀy, σ²·(m/(m+1))·(ZᵀZ)⁻¹)</c>; with weights, <c>β ~ N(Pa⁻¹h, σ²Pa⁻¹)</c> with
    /// <c>Pa = Zᵀ(W + I/m)Z</c>. Falls back to the mean when the response is zero; NaN when the design is
    /// inadmissible, which has no posterior.
    /// </summary>
    public double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => SampleConditional(design, y, weights, rng).Beta;

    /// <summary>The joint draw behind <see cref="SampleCoefficients"/>: the coefficients and the noise variance they
    /// were drawn with. σ² is NaN when no draw is possible. The robust weight update conditions on both.</summary>
    internal (double[] Beta, double Sigma2) SampleConditional(double[,] design, double[] y, double[]? weights,
                                                             Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (Condition(design, y, weights) is not GPriorConditional c)
            return (Unavailable(design), double.NaN);
        double[] beta = c.Mean.ToArray();
        if (!(c.S > 0.0))
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
        if (Condition(design, y, weights) is not GPriorConditional c || !(c.S > 0.0) || c.Rows <= 2)
            return ReadoutQuadratics.Unavailable(readoutDesign);
        return c.ConditionalVariance(null, readoutDesign);
    }

    /// <summary>
    /// Posterior mean spline coefficients: <c>(m/(m+1)) (ZᵀZ)⁻¹ Zᵀy</c> at unit weights (DMGK eq. 7's
    /// <c>E{f | k, ξ, y}</c>), and <c>Pa⁻¹h</c> with weights. Evaluated against a design at any points this is the
    /// conditional posterior mean curve for one knot draw; averaging over the chain is the Bayes fit.
    /// <paramref name="weights"/> null means unit weights. NaN when the design is inadmissible.
    /// </summary>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
        => Condition(design, y, weights) is GPriorConditional c ? c.Mean.ToArray() : Unavailable(design);

    // One knot configuration's conditional posterior at BARS's τ = 1/m, untempered; null when it is inadmissible.
    private static GPriorConditional? Condition(double[,] design, double[] y, double[]? weights)
    {
        (CollapsedGPrior? block, double tau) = Fit(design, y, weights);
        return block?.At(tau);
    }

    // The block for one knot configuration at BARS's τ = 1/m, or null when admission refuses a factor: G's alone
    // without weights; G's and P_a's with them.
    private static (CollapsedGPrior? Block, double Tau) Fit(double[,] design, double[] y, double[]? weights)
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
        CollapsedGPrior? block = weights is null
            ? BandedGPrior.Fit(banded, y, Admits)
            : BandedGPrior.Fit(banded, y, weights, tau, Admits);
        return (block, tau);
    }

    private static bool Admits(FactorStatus status) => DesignAdmission.Admits(status);

    // No posterior: one NaN per coefficient.
    private static double[] Unavailable(double[,] design)
    {
        var nan = new double[design.GetLength(1)];
        Array.Fill(nan, double.NaN);
        return nan;
    }

    // An admitted configuration's marginal. The b = 1 value is computed up front, since every chain and swap reads it.
    private sealed class NormalMarginal : ITemperedMarginal
    {
        private readonly CollapsedGPrior _block;
        private readonly double _tau;
        private readonly double _atOne;

        public NormalMarginal(CollapsedGPrior block, double tau)
        {
            _block = block;
            _tau = tau;
            _atOne = At(block.At(tau));
        }

        public DesignRejection Rejection => double.IsNegativeInfinity(_atOne) ? DesignRejection.ZeroResponse : DesignRejection.None;

        public double LogMarginalLikelihood(double inverseTemperature)
        {
            if (inverseTemperature == 1.0) return _atOne;
            if (_block.IsWeighted)
                throw new NotSupportedException("The weighted marginal is untempered: P_a is built at b = 1.");
            if (!(inverseTemperature > 0.0 && inverseTemperature <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(inverseTemperature),
                    "Under the improper 1/σ² prior the tempered marginal needs b in (0, 1].");
            return At(_block.At(_tau, inverseTemperature));
        }

        // An identically zero response has S = 0 at every b: no posterior for σ.
        private static double At(GPriorConditional c)
            => c.S > 0.0 ? c.LogDeterminantTerm - 0.5 * c.Rows * c.InverseTemperature * Math.Log(c.S) : double.NegativeInfinity;
    }
}
