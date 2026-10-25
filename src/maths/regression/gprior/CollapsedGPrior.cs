using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Blocks;

namespace Maths.Regression.GPrior;

/// <summary>Where a g-prior centres its coefficients.</summary>
public enum GPriorCentre
{
    /// <summary>
    /// On the response's median level: <c>μ₀ = med(y)·G⁻¹Zᵀ1</c> projects the median onto the design's span (<c>med(y)·1</c>
    /// for a basis whose functions sum to one). The prior then shrinks the fit's shape toward flat and leaves its level
    /// where the bulk of the data put it, so a constant shift of the response moves only the fit's level. No single
    /// response can drag the median far, so an outlier that robust weights discount cannot move the prior either. The level
    /// is the data's alone, never the sampled weights'.
    /// </summary>
    Median,

    /// <summary>
    /// On the response's mean level, the least-squares null fit: <c>μ₀ = ȳ·G⁻¹Zᵀ1</c>. Location-invariant like the median,
    /// and a proper counterpart of Liang et al.'s (2008) flat intercept. But an outlier of size R moves it by R/m even when
    /// robust weights discount the outlier, which shifts the fit by about R/m².
    /// </summary>
    Mean,

    /// <summary>
    /// At zero, DiMatteo, Genovese and Kass's original. The fit is shrunk toward zero by τ/(τ + b), level included, and
    /// the level's energy enters the scale S, so a response far from zero against its noise is under-fitted.
    /// </summary>
    Zero,
}

/// <summary>
/// The collapsed g-prior block shared by the free-knot engines (BARS now, BASS next): for a fixed design Z with Gram
/// <c>G = ZᵀZ</c>, the closed forms of <c>y ~ N(Zβ, σ²W⁻¹)</c> under <c>β | σ² ~ N(μ₀, σ²(τG)⁻¹)</c>, with the
/// likelihood raised to an inverse temperature b and σ² integrated out or drawn (BASS derivations review §§1–3, 5).
/// The prior mean μ₀ is zero unless the design layer centres it (<see cref="GPriorCentre"/>); the statistics the block
/// holds are then those of the centred response <c>ỹ = y − Zμ₀</c>, and μ₀ is added back to every mean and draw.
/// <list type="bullet">
/// <item><b>Unweighted (W = I).</b> The block holds the statistics of the least-squares fit that depend on neither τ
/// nor b: the factor of G, β_OLS, <c>‖r_OLS‖²</c> and <c>Q = ‖Zβ_OLS‖²</c>. Any (τ, b) then costs O(1) plus the mean
/// vector (<see cref="At"/>), which is what BASS's τ updates and tempered swaps need. The scale
/// <c>S = b·(‖r‖² + τ/(τ+b)·Q)</c> is a sum of non-negative terms from the explicit residual, so a small τ on a close
/// fit keeps its digits where <c>b·yᵀy − b²Q/(τ+b)</c> would not.</item>
/// <item><b>Weighted.</b> The robust BARS path (decision Q2): the prior stays on the unweighted G, and the posterior
/// precision <c>P_a = Zᵀ(W + τI)Z</c> depends on τ, so the block is built at one τ and b = 1.</item>
/// </list>
/// The block reports its factors' statuses and leaves admission to the caller, and it knows nothing about the basis:
/// a design layer such as <c>BandedGPrior</c> builds it.
/// </summary>
public sealed class CollapsedGPrior
{
    private readonly ISpdFactor _gram;
    private readonly ISpdFactor? _posterior;   // weighted only
    private readonly double[] _coefficients;   // of ỹ: β_OLS (unweighted) or the posterior mean P_a⁻¹ZᵀWỹ (weighted)
    private readonly double[]? _priorMean;     // μ₀; null is zero
    private readonly double _rss, _fitted;     // unweighted: ‖r_OLS‖², Q
    private readonly double _scale, _sumLogWeights, _tau;   // weighted: S, Σ ln wᵢ, the τ it was built at

    private CollapsedGPrior(ISpdFactor gram, ISpdFactor? posterior, double[] coefficients, double rss, double fitted,
                            double scale, double sumLogWeights, double tau, int rows, double[]? priorMean)
    {
        if (priorMean is not null && priorMean.Length != gram.Dimension)
            throw new ArgumentException("One prior mean per Gram column is required.", nameof(priorMean));
        _priorMean = priorMean;
        _gram = gram;
        _posterior = posterior;
        _coefficients = coefficients;
        _rss = rss;
        _fitted = fitted;
        _scale = scale;
        _sumLogWeights = sumLogWeights;
        _tau = tau;
        Rows = rows;
        Columns = gram.Dimension;
    }

    /// <summary>
    /// The unweighted block from a least-squares fit of the centred response ỹ = y − Zμ₀: the factor of G,
    /// <c>β_OLS = G⁻¹Zᵀỹ</c>, the residual sum of squares <c>‖ỹ − Zβ_OLS‖²</c> and the fitted sum of squares
    /// <c>Q = ‖Zβ_OLS‖²</c>, over <paramref name="rows"/> observations, with the prior mean
    /// <paramref name="priorMean"/> (null is zero, when ỹ = y). A caller that maintains these incrementally (BASS's
    /// bordered births) builds the block directly.
    /// </summary>
    public static CollapsedGPrior FromLeastSquares(ISpdFactor gram, double[] olsCoefficients, double residualSumOfSquares,
                                                   double fittedSumOfSquares, int rows, double[]? priorMean = null)
    {
        ArgumentNullException.ThrowIfNull(gram);
        ArgumentNullException.ThrowIfNull(olsCoefficients);
        if (olsCoefficients.Length != gram.Dimension) throw new ArgumentException("One coefficient per Gram column is required.", nameof(olsCoefficients));
        if (!(residualSumOfSquares >= 0.0)) throw new ArgumentOutOfRangeException(nameof(residualSumOfSquares));
        if (!(fittedSumOfSquares >= 0.0)) throw new ArgumentOutOfRangeException(nameof(fittedSumOfSquares));
        if (rows < 1) throw new ArgumentOutOfRangeException(nameof(rows));
        return new CollapsedGPrior(gram, null, olsCoefficients, residualSumOfSquares, fittedSumOfSquares,
                                   double.NaN, 0.0, double.NaN, rows, priorMean);
    }

    /// <summary>
    /// The weighted block at prior precision τ, from the centred response ỹ = y − Zμ₀: the factor of G (for
    /// <c>|τG|</c>), the factor of <c>P_a = Zᵀ(W + τI)Z</c>, the posterior mean <c>β̃ = P_a⁻¹ZᵀWỹ</c>, the scale
    /// <c>S = ‖ỹ − Zβ̃‖²_W + τ‖Zβ̃‖²</c> and <c>Σ ln wᵢ</c>, with the prior mean <paramref name="priorMean"/> (null is
    /// zero).
    /// </summary>
    public static CollapsedGPrior FromWeighted(ISpdFactor gram, ISpdFactor posterior, double[] posteriorMean, double scale,
                                               double sumLogWeights, int rows, double tau, double[]? priorMean = null)
    {
        ArgumentNullException.ThrowIfNull(gram);
        ArgumentNullException.ThrowIfNull(posterior);
        ArgumentNullException.ThrowIfNull(posteriorMean);
        if (posterior.Dimension != gram.Dimension || posteriorMean.Length != gram.Dimension)
            throw new ArgumentException("The factors and the mean must have one dimension.");
        if (!(scale >= 0.0)) throw new ArgumentOutOfRangeException(nameof(scale));
        if (rows < 1) throw new ArgumentOutOfRangeException(nameof(rows));
        CheckTau(tau);
        return new CollapsedGPrior(gram, posterior, posteriorMean, double.NaN, double.NaN, scale, sumLogWeights, tau, rows,
                                   priorMean);
    }

    /// <summary>Observations m.</summary>
    public int Rows { get; }

    /// <summary>Coefficients p (basis columns).</summary>
    public int Columns { get; }

    /// <summary>Whether per-point weights enter (the robust path).</summary>
    public bool IsWeighted => _posterior is not null;

    /// <summary>The prior mean μ₀; empty when it is zero.</summary>
    public ReadOnlySpan<double> PriorMean => _priorMean;

    /// <summary>The status of G's factor.</summary>
    public FactorStatus GramStatus => _gram.Status;

    /// <summary>The status of the posterior precision's factor: P_a's when weighted, otherwise G's.</summary>
    public FactorStatus PosteriorStatus => (_posterior ?? _gram).Status;

    /// <summary>
    /// The conditional posterior at prior precision <paramref name="tau"/> and likelihood inverse temperature
    /// <paramref name="inverseTemperature"/>. Unweighted: any τ &gt; 0 and b &gt; 0, at O(1) beyond the mean vector.
    /// Weighted: only the τ the block was built at, and b = 1, since P_a depends on both.
    /// </summary>
    public GPriorConditional At(double tau, double inverseTemperature = 1.0)
    {
        CheckTau(tau);
        double b = inverseTemperature;
        if (!(b > 0.0) || double.IsPositiveInfinity(b))
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature), "The inverse temperature must be positive and finite.");

        int p = Columns;
        if (_posterior is ISpdFactor posterior)
        {
            if (tau != _tau) throw new ArgumentException($"The weighted block was built at τ = {_tau}; P_a depends on τ.", nameof(tau));
            if (b != 1.0) throw new ArgumentException("The weighted block is untempered.", nameof(inverseTemperature));
            double logDet = 0.5 * (_gram.LogDet + p * Math.Log(tau)) - 0.5 * posterior.LogDet + 0.5 * _sumLogWeights;
            var weightedMean = new double[p];
            for (int j = 0; j < p; j++) weightedMean[j] = (_priorMean?[j] ?? 0.0) + _coefficients[j];
            return new GPriorConditional(posterior, weightedMean, _scale, logDet, 1.0, Rows, p, tau, b);
        }

        double shrink = b / (tau + b);
        var mean = new double[p];
        for (int j = 0; j < p; j++) mean[j] = (_priorMean?[j] ?? 0.0) + shrink * _coefficients[j];
        double s = b * (_rss + tau / (tau + b) * _fitted);
        return new GPriorConditional(_gram, mean, s, 0.5 * p * Math.Log(tau / (tau + b)), 1.0 / (tau + b), Rows, p, tau, b);
    }

    private static void CheckTau(double tau)
    {
        if (!(tau > 0.0) || double.IsPositiveInfinity(tau))
            throw new ArgumentOutOfRangeException(nameof(tau), "τ must be positive and finite.");
    }
}

/// <summary>
/// The block's answer at one (τ, b): <c>β | σ², y ~ N(μ, σ²·c·A⁻¹)</c> with A the factored precision (G, or P_a when
/// weighted) and c = <see cref="CovarianceScale"/>, and <c>σ² | y ~ IG(g₁ + mb/2, g₂ + S/2)</c> under a σ² prior
/// IG(g₁, g₂), or the improper <c>π(σ²) ∝ 1/σ²</c> (g₁ = g₂ = 0) that BARS uses.
/// </summary>
public sealed class GPriorConditional
{
    private static readonly double Log2Pi = Math.Log(2.0 * Math.PI);
    private readonly ISpdFactor _precision;
    private readonly double[] _mean;

    internal GPriorConditional(ISpdFactor precision, double[] mean, double s, double logDeterminantTerm, double covarianceScale,
                               int rows, int columns, double tau, double inverseTemperature)
    {
        _precision = precision;
        _mean = mean;
        S = s;
        LogDeterminantTerm = logDeterminantTerm;
        CovarianceScale = covarianceScale;
        Rows = rows;
        Columns = columns;
        Tau = tau;
        InverseTemperature = inverseTemperature;
    }

    /// <summary>Observations m.</summary>
    public int Rows { get; }

    /// <summary>Coefficients p.</summary>
    public int Columns { get; }

    /// <summary>Prior precision τ.</summary>
    public double Tau { get; }

    /// <summary>Likelihood inverse temperature b (1 is untempered).</summary>
    public double InverseTemperature { get; }

    /// <summary>The scale <c>b·S_b</c> that σ²'s posterior adds half of: <c>b(‖r‖² + τ/(τ+b)·Q)</c> unweighted,
    /// <c>‖y − Zβ̃‖²_W + τ‖Zβ̃‖²</c> weighted.</summary>
    public double S { get; }

    /// <summary>The coefficient integral's determinant term: <c>(p/2)·ln(τ/(τ+b))</c> unweighted;
    /// <c>½ln|τG| − ½ln|P_a| + ½Σ ln wᵢ</c> weighted.</summary>
    public double LogDeterminantTerm { get; }

    /// <summary>c in <c>Cov(β | σ²) = σ²·c·A⁻¹</c>: <c>1/(τ+b)</c> unweighted, 1 weighted.</summary>
    public double CovarianceScale { get; }

    /// <summary>The conditional posterior mean μ.</summary>
    public ReadOnlySpan<double> Mean => _mean;

    /// <summary><c>σ² | y ~ IG(g₁ + mb/2, g₂ + S/2)</c>; null <paramref name="sigmaPrior"/> is the improper 1/σ².</summary>
    public InverseGamma Sigma2Posterior(InverseGamma? sigmaPrior)
        => ConjugateUpdates.NormalVariance(sigmaPrior?.Shape ?? 0.0, sigmaPrior?.Scale ?? 0.0, Rows * InverseTemperature, S);

    /// <summary>
    /// The complete <c>ln p_b(y | Z, τ)</c> with β and σ² integrated out (derivations review §3):
    /// <c>−(mb/2)ln 2π + LogDeterminantTerm + g₁ln g₂ + lnΓ(g₁ + mb/2) − lnΓ(g₁) − (g₁ + mb/2)ln(g₂ + S/2)</c>. Under the
    /// improper prior the terms in g₁ and g₂ drop, leaving the result up to the improper prior's constant; that needs
    /// S &gt; 0. Comparable across τ and b, which a tempered swap needs.
    /// </summary>
    public double LogMarginal(InverseGamma? sigmaPrior)
    {
        double nb = 0.5 * Rows * InverseTemperature;
        double half = 0.5 * S;
        double common = -nb * Log2Pi + LogDeterminantTerm;
        if (sigmaPrior is InverseGamma prior)
            return common + prior.Shape * Math.Log(prior.Scale) + SpecialFunctions.LogGamma(prior.Shape + nb)
                 - SpecialFunctions.LogGamma(prior.Shape) - (prior.Shape + nb) * Math.Log(prior.Scale + half);
        if (!(S > 0.0)) throw new InvalidOperationException("The improper σ² prior needs S > 0.");
        return common + SpecialFunctions.LogGamma(nb) - nb * Math.Log(half);
    }

    /// <summary>
    /// A composition draw: σ² from its posterior, then <c>β = μ + √(σ²c)·v</c> with v ~ N(0, A⁻¹) from the factor
    /// (<see cref="GaussianConditional"/>). Writes β into <paramref name="beta"/> and returns σ². Scratch for p ≤ 128 lives
    /// on the stack.
    /// </summary>
    public double Sample(InverseGamma? sigmaPrior, Xoshiro256PlusPlus rng, Span<double> beta)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (beta.Length != Columns) throw new ArgumentException("One slot per coefficient is required.", nameof(beta));
        double sigma2 = Sigma2Posterior(sigmaPrior).Sample(rng);
        int p = Columns;
        Span<double> z = p <= 128 ? stackalloc double[p] : new double[p];
        Span<double> innovation = p <= 128 ? stackalloc double[p] : new double[p];
        GaussianConditional.Perturb(_precision, _mean, rng, beta, Math.Sqrt(sigma2 * CovarianceScale), z, innovation);
        return sigma2;
    }

    /// <summary>
    /// <c>Var[z(t)ᵀβ | y] = E[σ² | y]·c·z(t)ᵀA⁻¹z(t)</c> for each row z(t) of <paramref name="readoutDesign"/>, with σ²
    /// integrated out. Needs a posterior σ² with finite mean: shape above 1.
    /// </summary>
    public double[] ConditionalVariance(InverseGamma? sigmaPrior, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        InverseGamma posterior = Sigma2Posterior(sigmaPrior);
        if (!(posterior.Shape > 1.0)) throw new InvalidOperationException("σ²'s posterior mean is infinite (shape ≤ 1).");
        return _precision.InverseQuadraticForms(readoutDesign, posterior.Mean * CovarianceScale);
    }
}
