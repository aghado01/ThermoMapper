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

    /// <summary>
    /// Liang et al.'s (2008) flat intercept: the intercept takes a flat prior and is integrated out, the other columns are
    /// centred and take the g-prior, and k counts them alone. Shift-invariant with no double use of the data, at the cost
    /// of a proper prior (BASS's D03 trilemma); BASS selects it as an option. It needs the constant split out of the basis,
    /// as the block's centred form takes it (<see cref="CollapsedGPrior.FromCentredLeastSquares"/>); B-splines do not
    /// split it out, so BARS refuses it.
    /// </summary>
    Flat,
}

/// <summary>
/// The collapsed g-prior block shared by the free-knot engines (BARS now, BASS next): for a fixed design Z with Gram
/// <c>G = ZᵀZ</c>, the closed forms of <c>y ~ N(Zβ, σ²W⁻¹)</c> under <c>β | σ² ~ N(μ₀, σ²(τG)⁻¹)</c>, with the
/// likelihood raised to an inverse temperature b and σ² integrated out or drawn (BASS derivations review §§1–3, 5).
/// The prior mean μ₀ is zero unless the design layer centres it (<see cref="GPriorCentre"/>); the statistics the block
/// holds are then those of the centred response <c>ỹ = y − Zμ₀</c>, and μ₀ is added back to every mean and draw.
/// <list type="bullet">
/// <item><b>Unweighted (W = I).</b> The block holds the statistics of the least-squares fit that depend on neither τ
/// nor b: the factor of G, β_OLS, <c>‖r_OLS‖²</c> and <c>Q = ‖Zβ_OLS‖²</c>. Any (τ, b) then costs O(1) for the marginal
/// (<see cref="MarginalAt"/>), which is what BASS's τ updates and tempered swaps need, and the mean vector beyond it for
/// the conditional (<see cref="At"/>). The scale
/// <c>S = b·(‖r‖² + τ/(τ+b)·Q)</c> is a sum of non-negative terms from the explicit residual, so a small τ on a close
/// fit keeps its digits where <c>b·yᵀy − b²Q/(τ+b)</c> would not.</item>
/// <item><b>Weighted.</b> The robust BARS path (decision Q2): the prior stays on the unweighted G, and the posterior
/// precision <c>P_a = Zᵀ(W + τI)Z</c> depends on τ, so the block is built at one τ and b = 1.</item>
/// <item><b>Centred, with the intercept apart.</b> BASS's design (D03): the factor is of the centred columns' Gram, the
/// constant stands orthogonal to them, and its convention (<see cref="Intercept"/>) places it inside the g-prior or
/// under Liang's flat prior (<see cref="FromCentredLeastSquares"/>). Unweighted, any (τ, b) as above.</item>
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
    private readonly double _responseMean, _level;          // centred: ȳ, and the intercept's centre (ȳ under Flat)

    private CollapsedGPrior(ISpdFactor gram, ISpdFactor? posterior, double[] coefficients, double rss, double fitted,
                            double scale, double sumLogWeights, double tau, int rows, double[]? priorMean,
                            GPriorCentre? intercept = null, double responseMean = 0.0, double level = 0.0)
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
        _responseMean = responseMean;
        _level = level;
        Intercept = intercept;
        Rows = rows;
        Columns = gram.Dimension + (intercept is null ? 0 : 1);
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
    /// The unweighted block for a design whose intercept stands apart from its centred columns, as BASS's design holds it
    /// (D03): Z_c centred over all <paramref name="rows"/> observations, so the constant is orthogonal to every column and
    /// its Gram entry N stands apart. It takes the factor of <c>G_c = Z_cᵀZ_c</c> and the least-squares fit of the centred
    /// response ỹ = y − ȳ on Z_c (<c>β̂_c</c>, <c>‖ỹ − Z_cβ̂_c‖²</c>, <c>Q_c = ‖Z_cβ̂_c‖²</c>), with ȳ; an empty factor is
    /// the intercept-only structure. The coefficients are (α, β_c), the intercept first, and
    /// <paramref name="intercept"/> places α:
    /// <list type="bullet">
    /// <item><see cref="GPriorCentre.Mean"/>, level centring: inside the g-prior and centred on ȳ, so k = M + 1. The
    /// explicit design [1, Z_c] with μ₀ = ȳ·e₀ through <see cref="FromLeastSquares"/> is the same block.</item>
    /// <item><see cref="GPriorCentre.Zero"/>: inside the g-prior and centred at 0, k = M + 1, with α's fitted energy Nȳ²
    /// added to Q.</item>
    /// <item><see cref="GPriorCentre.Flat"/>, Liang et al.'s (2008) convention: a flat prior on α, integrated out, and
    /// k = M. At inverse temperature b the integral gives <c>(2πσ²/(Nb))^½</c>, which spends one observation: σ²'s shape
    /// reads <c>(Nb − 1)/2</c>, and b = 0 is out of reach on this form (the linear form reaches it, D05). The marginal is
    /// defined up to the flat prior's constant, which cancels from every within-rung ratio and every swap, so it is no
    /// evidence. Given σ², <c>α ~ N(ȳ, σ²/(Nb))</c> apart from β_c.</item>
    /// </list>
    /// The median centre is the robust path's, on a basis that carries the constant (BARS), and is refused here.
    /// </summary>
    public static CollapsedGPrior FromCentredLeastSquares(ISpdFactor centredGram, double[] centredCoefficients,
                                                          double residualSumOfSquares, double centredFittedSumOfSquares,
                                                          int rows, double responseMean, GPriorCentre intercept)
    {
        ArgumentNullException.ThrowIfNull(centredGram);
        ArgumentNullException.ThrowIfNull(centredCoefficients);
        if (centredCoefficients.Length != centredGram.Dimension)
            throw new ArgumentException("One coefficient per centred column is required.", nameof(centredCoefficients));
        if (!(residualSumOfSquares >= 0.0)) throw new ArgumentOutOfRangeException(nameof(residualSumOfSquares));
        if (!(centredFittedSumOfSquares >= 0.0)) throw new ArgumentOutOfRangeException(nameof(centredFittedSumOfSquares));
        (double fitted, double level) = CentredFit(centredFittedSumOfSquares, rows, responseMean, intercept);
        return new CollapsedGPrior(centredGram, null, centredCoefficients, residualSumOfSquares, fitted,
                                   double.NaN, 0.0, double.NaN, rows, null, intercept, responseMean, level);
    }

    /// <summary>
    /// The marginal of <see cref="FromCentredLeastSquares"/>'s block at (<paramref name="tau"/>,
    /// <paramref name="inverseTemperature"/>), from the fit's statistics alone: <paramref name="centredColumns"/> centred
    /// columns beside the intercept, ‖r‖² and Q_c. A design that weighs a candidate structure has these before it has
    /// a factor or coefficients, and the marginal reads nothing else.
    /// </summary>
    public static GPriorMarginal CentredMarginalAt(int centredColumns, double residualSumOfSquares, double centredFittedSumOfSquares,
                                                   int rows, double responseMean, GPriorCentre intercept, double tau,
                                                   double inverseTemperature = 1.0)
    {
        if (centredColumns < 0) throw new ArgumentOutOfRangeException(nameof(centredColumns));
        if (!(residualSumOfSquares >= 0.0)) throw new ArgumentOutOfRangeException(nameof(residualSumOfSquares));
        (double fitted, _) = CentredFit(centredFittedSumOfSquares, rows, responseMean, intercept);
        return Unweighted(centredColumns + 1, residualSumOfSquares, fitted, rows, intercept, tau, inverseTemperature);
    }

    // The fitted energy the centred block's convention reads, and the intercept's centre.
    private static (double Fitted, double Level) CentredFit(double centredFitted, int rows, double responseMean, GPriorCentre intercept)
    {
        if (!(centredFitted >= 0.0)) throw new ArgumentOutOfRangeException(nameof(centredFitted));
        if (rows < 1) throw new ArgumentOutOfRangeException(nameof(rows));
        if (!double.IsFinite(responseMean)) throw new ArgumentOutOfRangeException(nameof(responseMean));
        return intercept switch
        {
            GPriorCentre.Mean or GPriorCentre.Flat => (centredFitted, responseMean),
            GPriorCentre.Zero => (centredFitted + rows * responseMean * responseMean, 0.0),
            GPriorCentre.Median => throw new ArgumentException(
                "The median centre is the robust path's, on a basis that carries the constant; a design with its intercept "
                + "apart takes the mean, the zero centre or the flat intercept.", nameof(intercept)),
            _ => throw new ArgumentOutOfRangeException(nameof(intercept)),
        };
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

    /// <summary>Coefficients p: the basis columns, and the intercept first when it stands apart.</summary>
    public int Columns { get; }

    /// <summary>Whether per-point weights enter (the robust path).</summary>
    public bool IsWeighted => _posterior is not null;

    /// <summary>The convention of an intercept that stands apart from the centred columns
    /// (<see cref="FromCentredLeastSquares"/>); null when the basis carries the constant itself. Under
    /// <see cref="GPriorCentre.Flat"/> the marginal is relative, so it gives no evidence.</summary>
    public GPriorCentre? Intercept { get; }

    /// <summary>The prior mean μ₀ of the factor's columns; empty when it is zero, as it is for centred columns.</summary>
    public ReadOnlySpan<double> PriorMean => _priorMean;

    /// <summary>The status of G's factor.</summary>
    public FactorStatus GramStatus => _gram.Status;

    /// <summary>The status of the posterior precision's factor: P_a's when weighted, otherwise G's.</summary>
    public FactorStatus PosteriorStatus => (_posterior ?? _gram).Status;

    /// <summary>
    /// The conditional posterior at prior precision <paramref name="tau"/> and likelihood inverse temperature
    /// <paramref name="inverseTemperature"/>. Unweighted: any τ &gt; 0 and b ≥ 0, at O(1) beyond the mean vector. At b = 0
    /// the likelihood is gone (S = 0, the mean is μ₀): a proper σ² prior's marginal there is exactly 0, the start of an
    /// SMC run or an exact reference rung, and the improper prior's is undefined. The flat intercept needs b &gt; 0.
    /// Weighted: only the τ the block was built at, and b = 1, since P_a depends on both. Its marginal is
    /// <see cref="MarginalAt"/>'s.
    /// </summary>
    public GPriorConditional At(double tau, double inverseTemperature = 1.0)
    {
        GPriorMarginal marginal = MarginalAt(tau, inverseTemperature);
        int p = Columns;
        var mean = new double[p];
        if (_posterior is ISpdFactor posterior)
        {
            for (int j = 0; j < p; j++) mean[j] = (_priorMean?[j] ?? 0.0) + _coefficients[j];
            return new GPriorConditional(posterior, mean, marginal, 1.0, p);
        }

        double b = inverseTemperature, shrink = b / (tau + b);
        if (Intercept is GPriorCentre intercept)
        {
            // α is shrunk towards its centre inside the g-prior, and not at all under the flat prior.
            bool flat = intercept == GPriorCentre.Flat;
            mean[0] = flat ? _responseMean : _level + shrink * (_responseMean - _level);
            for (int j = 1; j < p; j++) mean[j] = shrink * _coefficients[j - 1];
            return new GPriorConditional(_gram, mean, marginal, 1.0 / (tau + b), p, flat ? 1.0 / b : 1.0 / (tau + b));
        }
        for (int j = 0; j < p; j++) mean[j] = (_priorMean?[j] ?? 0.0) + shrink * _coefficients[j];
        return new GPriorConditional(_gram, mean, marginal, 1.0 / (tau + b), p);
    }

    /// <summary>
    /// The marginal at (<paramref name="tau"/>, <paramref name="inverseTemperature"/>) without the conditional mean: S and
    /// the determinant term at O(1), and nothing allocated, for the τ updates and tempered swaps that evaluate it many
    /// times per sweep. The same domain as <see cref="At"/>.
    /// </summary>
    public GPriorMarginal MarginalAt(double tau, double inverseTemperature = 1.0)
    {
        if (_posterior is ISpdFactor posterior)
        {
            CheckTau(tau);
            double b = inverseTemperature;
            if (tau != _tau) throw new ArgumentException($"The weighted block was built at τ = {_tau}; P_a depends on τ.", nameof(tau));
            if (b != 1.0) throw new ArgumentException("The weighted block is untempered.", nameof(inverseTemperature));
            double logDet = 0.5 * (_gram.LogDet + Columns * Math.Log(tau)) - 0.5 * posterior.LogDet + 0.5 * _sumLogWeights;
            return new GPriorMarginal(_scale, logDet, Rows, tau, b);
        }
        return Unweighted(Columns, _rss, _fitted, Rows, Intercept, tau, inverseTemperature);
    }

    // The unweighted marginal at (τ, b), with S and the determinant term's slopes in b beside them:
    // ∂S/∂b = ‖r‖² + (τ/(τ+b))²·Q and ∂/∂b (k/2)ln(τ/(τ+b)) = −(k/2)/(τ+b).
    private static GPriorMarginal Unweighted(int columns, double rss, double fitted, int rows, GPriorCentre? intercept,
                                             double tau, double b)
    {
        CheckTau(tau);
        if (!(b >= 0.0) || double.IsPositiveInfinity(b))
            throw new ArgumentOutOfRangeException("inverseTemperature", "The inverse temperature must be non-negative and finite.");
        double shrink = tau / (tau + b);
        double s = b * (rss + shrink * fitted), sSlope = rss + shrink * shrink * fitted;
        if (intercept != GPriorCentre.Flat)
            return new GPriorMarginal(s, 0.5 * columns * Math.Log(shrink), rows, tau, b, 0, sSlope, -0.5 * columns / (tau + b));

        // α's integral, (2πσ²/(Nb))^½: ½ln 2π moves into the count's 2π term, −½ln(Nb) stays here.
        if (!(b > 0.0))
            throw new ArgumentOutOfRangeException("inverseTemperature",
                "b = 0 is out of reach on the flat intercept's integrated form, whose (2πσ²/(Nb))^½ diverges there.");
        double flatLogDet = 0.5 * (columns - 1) * Math.Log(shrink) - 0.5 * Math.Log(rows * b);
        return new GPriorMarginal(s, flatLogDet, rows, tau, b, 1, sSlope, -0.5 * (columns - 1) / (tau + b) - 0.5 / b);
    }

    private static void CheckTau(double tau)
    {
        if (!(tau > 0.0) || double.IsPositiveInfinity(tau))
            throw new ArgumentOutOfRangeException(nameof(tau), "τ must be positive and finite.");
    }
}

/// <summary>
/// The block's marginal at one (τ, b), without the conditional mean: the scale S, the determinant term and the row count,
/// which are all that the log marginal and σ²'s posterior read. A value, so the τ updates and tempered swaps that
/// evaluate it many times per sweep allocate nothing; the <see cref="GPriorConditional"/> at the same (τ, b) carries it.
/// </summary>
public readonly struct GPriorMarginal
{
    private static readonly double Log2Pi = Math.Log(2.0 * Math.PI);

    private readonly int _spent;   // observations a flat intercept's integral spends: 0 or 1
    private readonly double _sSlope, _logDeterminantSlope;   // ∂/∂b of S and of the determinant term; NaN when weighted

    internal GPriorMarginal(double s, double logDeterminantTerm, int rows, double tau, double inverseTemperature,
                            int spent = 0, double sSlope = double.NaN, double logDeterminantSlope = double.NaN)
    {
        S = s;
        LogDeterminantTerm = logDeterminantTerm;
        Rows = rows;
        Tau = tau;
        InverseTemperature = inverseTemperature;
        _spent = spent;
        _sSlope = sSlope;
        _logDeterminantSlope = logDeterminantSlope;
    }

    /// <summary>Observations m.</summary>
    public int Rows { get; }

    /// <summary>ν, the count σ²'s shape adds half of: <c>mb</c>, or <c>mb − 1</c> under the flat intercept, whose integral
    /// spends one observation.</summary>
    public double DegreesOfFreedom => Rows * InverseTemperature - _spent;

    /// <summary>Prior precision τ.</summary>
    public double Tau { get; }

    /// <summary>Likelihood inverse temperature b (1 is untempered).</summary>
    public double InverseTemperature { get; }

    /// <summary>The scale <c>b·S_b</c> that σ²'s posterior adds half of: <c>b(‖r‖² + τ/(τ+b)·Q)</c> unweighted,
    /// <c>‖y − Zβ̃‖²_W + τ‖Zβ̃‖²</c> weighted.</summary>
    public double S { get; }

    /// <summary>The coefficient integral's determinant term: <c>(k/2)·ln(τ/(τ+b))</c> unweighted, k the g-prior's
    /// columns, less <c>½ln(mb)</c> under the flat intercept; <c>½ln|τG| − ½ln|P_a| + ½Σ ln wᵢ</c> weighted.</summary>
    public double LogDeterminantTerm { get; }

    /// <summary><c>σ² | y ~ IG(g₁ + ν/2, g₂ + S/2)</c>, ν = <see cref="DegreesOfFreedom"/>; null
    /// <paramref name="sigmaPrior"/> is the improper 1/σ².</summary>
    public InverseGamma Sigma2Posterior(InverseGamma? sigmaPrior)
        => ConjugateUpdates.NormalVariance(sigmaPrior?.Shape ?? 0.0, sigmaPrior?.Scale ?? 0.0, DegreesOfFreedom, S);

    /// <summary>
    /// The complete <c>ln p_b(y | Z, τ)</c> with β and σ² integrated out (derivations review §3), ν =
    /// <see cref="DegreesOfFreedom"/>:
    /// <c>−(ν/2)ln 2π + LogDeterminantTerm + g₁ln g₂ + lnΓ(g₁ + ν/2) − lnΓ(g₁) − (g₁ + ν/2)ln(g₂ + S/2)</c>. Under the
    /// improper prior the terms in g₁ and g₂ drop, leaving the result up to the improper prior's constant; that needs
    /// S &gt; 0 and ν &gt; 0. The flat intercept's is likewise up to its prior's constant. Comparable across τ and b,
    /// which a tempered swap needs.
    /// </summary>
    public double LogMarginal(InverseGamma? sigmaPrior)
    {
        double nb = 0.5 * Rows * InverseTemperature - 0.5 * _spent;
        double half = 0.5 * S;
        double common = -nb * Log2Pi + LogDeterminantTerm;
        if (sigmaPrior is InverseGamma prior)
        {
            if (!(prior.Shape + nb > 0.0)) throw new InvalidOperationException("σ²'s posterior shape g₁ + ν/2 is not positive.");
            return InverseTemperature == 0.0
                ? 0.0   // no likelihood, and a proper prior integrates to one; the terms below cancel only to rounding
                : common + prior.Shape * Math.Log(prior.Scale) + SpecialFunctions.LogGamma(prior.Shape + nb)
                 - SpecialFunctions.LogGamma(prior.Shape) - (prior.Shape + nb) * Math.Log(prior.Scale + half);
        }
        if (!(S > 0.0)) throw new InvalidOperationException("The improper σ² prior needs S > 0.");
        if (!(nb > 0.0)) throw new InvalidOperationException("The improper σ² prior needs ν > 0.");
        return common + SpecialFunctions.LogGamma(nb) - nb * Math.Log(half);
    }

    /// <summary>
    /// <c>∂/∂b</c> of <see cref="LogMarginal"/> at this b, the integrated form's slope on the likelihood path (path
    /// sampling's integrand), with ν' = m:
    /// <c>−(m/2)ln 2π + ∂LogDeterminantTerm + (m/2)ψ(g₁ + ν/2) − (m/2)ln(g₂ + S/2) − (g₁ + ν/2)(S'/2)/(g₂ + S/2)</c>.
    /// The same domain as <see cref="LogMarginal"/>, b = 0 included under a proper prior, where it carries
    /// −(k/2)/τ; unweighted blocks only.
    /// </summary>
    public double LogMarginalSlope(InverseGamma? sigmaPrior)
    {
        if (double.IsNaN(_sSlope)) throw new NotSupportedException("The weighted block is untempered, so it has no slope in b.");
        double nb = 0.5 * Rows * InverseTemperature - 0.5 * _spent, m = 0.5 * Rows;
        double half = 0.5 * S;
        double common = -m * Log2Pi + _logDeterminantSlope;
        double shape = sigmaPrior?.Shape ?? 0.0, scale = sigmaPrior?.Scale ?? 0.0;
        if (sigmaPrior is null)
        {
            if (!(S > 0.0)) throw new InvalidOperationException("The improper σ² prior needs S > 0.");
            if (!(nb > 0.0)) throw new InvalidOperationException("The improper σ² prior needs ν > 0.");
        }
        else if (!(shape + nb > 0.0)) throw new InvalidOperationException("σ²'s posterior shape g₁ + ν/2 is not positive.");
        return common + m * SpecialFunctions.Digamma(shape + nb) - m * Math.Log(scale + half)
               - (shape + nb) * 0.5 * _sSlope / (scale + half);
    }
}

/// <summary>
/// The block's answer at one (τ, b): <c>β | σ², y ~ N(μ, σ²·c·A⁻¹)</c> with A the factored precision (G, or P_a when
/// weighted) and c = <see cref="CovarianceScale"/>, and <c>σ² | y ~ IG(g₁ + mb/2, g₂ + S/2)</c> under a σ² prior
/// IG(g₁, g₂), or the improper <c>π(σ²) ∝ 1/σ²</c> (g₁ = g₂ = 0) that BARS uses. Its marginal is the
/// <see cref="GPriorMarginal"/> at the same (τ, b). With the intercept apart from centred columns, coefficient 0 is α,
/// independent of the others given σ², with <c>Var(α | σ²) = σ²·c₀/m</c> (<see cref="InterceptCovarianceScale"/>), and
/// A is the centred columns' Gram.
/// </summary>
public sealed class GPriorConditional
{
    private readonly ISpdFactor _precision;
    private readonly double[] _mean;
    private readonly GPriorMarginal _marginal;

    internal GPriorConditional(ISpdFactor precision, double[] mean, GPriorMarginal marginal, double covarianceScale, int columns,
                               double? interceptCovarianceScale = null)
    {
        _precision = precision;
        _mean = mean;
        _marginal = marginal;
        CovarianceScale = covarianceScale;
        Columns = columns;
        InterceptCovarianceScale = interceptCovarianceScale;
    }

    /// <inheritdoc cref="GPriorMarginal.Rows"/>
    public int Rows => _marginal.Rows;

    /// <summary>Coefficients p.</summary>
    public int Columns { get; }

    /// <inheritdoc cref="GPriorMarginal.Tau"/>
    public double Tau => _marginal.Tau;

    /// <inheritdoc cref="GPriorMarginal.InverseTemperature"/>
    public double InverseTemperature => _marginal.InverseTemperature;

    /// <inheritdoc cref="GPriorMarginal.S"/>
    public double S => _marginal.S;

    /// <inheritdoc cref="GPriorMarginal.LogDeterminantTerm"/>
    public double LogDeterminantTerm => _marginal.LogDeterminantTerm;

    /// <summary>c in <c>Cov(β | σ²) = σ²·c·A⁻¹</c>: <c>1/(τ+b)</c> unweighted, 1 weighted.</summary>
    public double CovarianceScale { get; }

    /// <summary>c₀ in <c>Var(α | σ²) = σ²·c₀/m</c> for an intercept apart from centred columns: <c>1/(τ+b)</c> inside the
    /// g-prior, <c>1/b</c> under the flat prior; null when the basis carries the constant.</summary>
    public double? InterceptCovarianceScale { get; }

    /// <summary>The conditional posterior mean μ.</summary>
    public ReadOnlySpan<double> Mean => _mean;

    /// <inheritdoc cref="GPriorMarginal.Sigma2Posterior"/>
    public InverseGamma Sigma2Posterior(InverseGamma? sigmaPrior) => _marginal.Sigma2Posterior(sigmaPrior);

    /// <inheritdoc cref="GPriorMarginal.LogMarginal"/>
    public double LogMarginal(InverseGamma? sigmaPrior) => _marginal.LogMarginal(sigmaPrior);

    /// <summary>
    /// A composition draw: σ² from its posterior, then <c>β = μ + √(σ²c)·v</c> with v ~ N(0, A⁻¹) from the factor
    /// (<see cref="GaussianConditional"/>), then an intercept apart from them from its normal. Writes β into
    /// <paramref name="beta"/> and returns σ². Scratch for p ≤ 128 lives on the stack.
    /// </summary>
    public double Sample(InverseGamma? sigmaPrior, Xoshiro256PlusPlus rng, Span<double> beta)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (beta.Length != Columns) throw new ArgumentException("One slot per coefficient is required.", nameof(beta));
        double sigma2 = Sigma2Posterior(sigmaPrior).Sample(rng);
        int offset = InterceptCovarianceScale is null ? 0 : 1, p = Columns - offset;
        Span<double> z = p <= 128 ? stackalloc double[p] : new double[p];
        Span<double> innovation = p <= 128 ? stackalloc double[p] : new double[p];
        if (p > 0)
            GaussianConditional.Perturb(_precision, _mean.AsSpan(offset), rng, beta[offset..], Math.Sqrt(sigma2 * CovarianceScale),
                                        z, innovation);
        if (InterceptCovarianceScale is double c0)
            beta[0] = _mean[0] + Math.Sqrt(sigma2 * c0 / Rows) * Normal.SampleStandard(rng);
        return sigma2;
    }

    /// <summary>
    /// ln p(y_U | y_rest) for one unit U of the observations this conditional was fitted to (OQ17's collapsed content). β
    /// and σ² are integrated out, and the prior (its centre, G and τ) is held at the full data's. It is the leave-unit-out
    /// predictive of y | σ² ~ N(Zμ₀, σ²Q⁻¹), with Q = W − WZA⁻¹ZᵀW and A the posterior precision (A⁻¹ = c·factor⁻¹).
    /// <list type="bullet">
    ///   <item>Qr = We, with e = y − Zμ the posterior-mean residual, and rᵀQr = S.</item>
    ///   <item>The unit's block is Q_UU = W_U − W_U Z_U A⁻¹ Z_Uᵀ W_U: m solves with the factor.</item>
    ///   <item>The rest follows <see cref="LeaveOutPredictive"/>, a Student t with 2a + n − m degrees of freedom under
    ///   σ²'s prior IG(a, b); null <paramref name="sigmaPrior"/> is the improper 1/σ², a = b = 0.</item>
    /// </list>
    /// With the intercept apart, <paramref name="design"/>'s column 0 is its column, α's variance enters A⁻¹ apart from
    /// the factor's, and a flat intercept spends one of the n. Untempered only. <paramref name="weights"/> are the fit's
    /// (null for W = I).
    /// </summary>
    public double LeaveOutLogDensity(InverseGamma? sigmaPrior, double[,] design, ReadOnlySpan<double> y, double[]? weights,
                                     int[] unit)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(unit);
        if (InverseTemperature != 1.0) throw new InvalidOperationException("The leave-out predictive is untempered.");
        if (design.GetLength(0) != Rows || design.GetLength(1) != Columns || y.Length != Rows)
            throw new ArgumentException("The design and response must be the ones this conditional was fitted to.", nameof(design));
        int m = unit.Length, p = Columns, offset = InterceptCovarianceScale is null ? 0 : 1;
        var solved = new double[m][];
        var row = new double[p];
        var qr = new double[m];
        var w = new double[m];
        for (int a = 0; a < m; a++)
        {
            int i = unit[a];
            double fit = 0.0;
            for (int j = 0; j < p; j++) { row[j] = design[i, j]; fit += row[j] * _mean[j]; }
            solved[a] = _precision.Solve(row.AsSpan(offset));
            w[a] = weights?[i] ?? 1.0;
            qr[a] = w[a] * (y[i] - fit);
        }
        var block = new double[m, m];
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
            {
                double quadratic = 0.0;
                for (int j = offset; j < p; j++) quadratic += design[unit[a], j] * solved[b][j - offset];
                block[a, b] = InterceptCovarianceScale is double c0
                    ? (a == b ? w[a] : 0.0) - w[a] * w[b] * (CovarianceScale * quadratic + c0 / Rows * design[unit[a], 0] * design[unit[b], 0])
                    : (a == b ? w[a] : 0.0) - w[a] * w[b] * CovarianceScale * quadratic;
            }
        return LeaveOutPredictive.StudentT(block, qr, S, _marginal.DegreesOfFreedom, sigmaPrior?.Shape ?? 0.0, sigmaPrior?.Scale ?? 0.0);
    }

    /// <summary>
    /// <c>Var[z(t)ᵀβ | y] = E[σ² | y]·c·z(t)ᵀA⁻¹z(t)</c> for each row z(t) of <paramref name="readoutDesign"/>, with σ²
    /// integrated out; an intercept apart adds <c>E[σ² | y]·c₀/m·z₀(t)²</c>. Needs a posterior σ² with finite mean: shape
    /// above 1.
    /// </summary>
    public double[] ConditionalVariance(InverseGamma? sigmaPrior, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        InverseGamma posterior = Sigma2Posterior(sigmaPrior);
        if (!(posterior.Shape > 1.0)) throw new InvalidOperationException("σ²'s posterior mean is infinite (shape ≤ 1).");
        if (InterceptCovarianceScale is not double c0)
            return _precision.InverseQuadraticForms(readoutDesign, posterior.Mean * CovarianceScale);

        int g = readoutDesign.GetLength(0), p = readoutDesign.GetLength(1);
        if (p != Columns) throw new ArgumentException("One readout column per coefficient is required.", nameof(readoutDesign));
        var centred = new double[g, p - 1];
        for (int t = 0; t < g; t++)
            for (int j = 1; j < p; j++) centred[t, j - 1] = readoutDesign[t, j];
        double[] variance = p > 1 ? _precision.InverseQuadraticForms(centred, posterior.Mean * CovarianceScale) : new double[g];
        for (int t = 0; t < g; t++) variance[t] += posterior.Mean * c0 / Rows * readoutDesign[t, 0] * readoutDesign[t, 0];
        return variance;
    }
}
