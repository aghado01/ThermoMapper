using System;
using Maths.LinAlg;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// The penalized P-spline as a Gaussian mixed model (Eilers &amp; Marx 1996; Wand's REML representation) — BAPS's
/// fixed-basis, measure-side core. On a rich fixed B-spline basis the roughness penalty <c>λ·DᵀD</c> is the
/// precision of an (improper) Gaussian prior on β, so <c>β̂(λ) = (ZᵀZ + λP)⁻¹ Zᵀy</c> is the posterior mean — the
/// Reinsch band solve — and the smoothing λ is chosen by the marginal/REML evidence rather than tuned. The flat
/// prior on the penalty's null space (polynomials of degree &lt; r) is exactly the REML integration of the fixed
/// effects, so no eigen-reparametrization or dense n×n covariance is needed: <c>A = ZᵀZ + λP</c> stays banded at
/// half-bandwidth max(degree, r), is factored once by <see cref="BandCholesky"/>, and the evidence reads off its
/// log-determinant. The (λ, β) sampler — Gibbs over variance components, or λ from this marginal evidence — sits
/// on top of this object.
/// </summary>
public sealed class PenalizedSpline
{
    private readonly BandedDesign _design;
    private readonly IBandPenalty _penalty;
    private readonly double[] _y;
    private readonly int _n;
    private readonly int _nu;
    private readonly int _bw;
    private readonly double _yy;   // yᵀy
    private readonly double[,] _gram;   // ZᵀZ band, accumulated once: every solve copies it
    private readonly double[] _zty;     // Zᵀy
    private readonly Lazy<double[]> _spectrum;   // generalized eigenvalues μ_j of (P, ZᵀZ + P), for edf

    public PenalizedSpline(double[,] design, double[] y, IBandPenalty penalty)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(penalty);
        _design = new BandedDesign(design);
        _penalty = penalty;
        _y = y;
        _n = _design.Rows;
        _nu = _design.Dimension;
        if (y.Length != _n)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        _bw = Math.Max(_design.Bandwidth, penalty.Bandwidth);
        _gram = new double[_bw + 1, _nu];
        _zty = new double[_nu];
        _design.Accumulate(null, _gram, y, _zty);
        double yy = 0.0;
        for (int i = 0; i < _n; i++) yy += y[i] * y[i];
        _yy = yy;
        _spectrum = new Lazy<double[]>(PenaltySpectrum);
    }

    /// <summary>Number of basis functions ν.</summary>
    public int Dimension => _nu;

    /// <summary>Penalty null-space dimension r (= penalty order) — the unpenalized polynomial degrees of freedom.</summary>
    public int Nullity => _penalty.Nullity;

    // (β̂, penalized RSS, log|A|) at smoothing λ; A = ZᵀZ + λP factored once and shared by every read-out.
    private (double[] beta, double rssPen, double logDetA) Solve(double lambda)
    {
        PenalizedSolve solve = SolveFactored(lambda);
        return (solve.Beta, solve.RssPen, solve.Factor.LogDet);
    }

    /// <summary>β̂(λ) = A⁻¹Zᵀy, RSS_pen(λ) = yᵀy − β̂ᵀZᵀy and the factor of A = ZᵀZ + λP — what an exact conditional draw
    /// of (σ², β) given λ needs.</summary>
    internal PenalizedSolve SolveFactored(double lambda)
    {
        var band = (double[,])_gram.Clone();          // ZᵀZ
        _penalty.AccumulateInto(band, _nu, lambda);   // + λP
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        if (chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
            throw new InvalidOperationException($"ZᵀZ + λP is not numerically positive definite at λ = {lambda} ({chol.Status}).");
        double[] beta = chol.Solve(_zty);

        double bBeta = 0.0;                            // RSS_pen = yᵀy − β̂ᵀZᵀy
        for (int p = 0; p < _nu; p++) bBeta += _zty[p] * beta[p];
        return new PenalizedSolve(beta, _yy - bBeta, chol);
    }

    /// <summary>Posterior-mean (penalized) coefficients <c>β̂(λ) = (ZᵀZ + λP)⁻¹ Zᵀy</c> — the Reinsch solve.</summary>
    public double[] Coefficients(double lambda) => Solve(lambda).beta;

    /// <summary>Profiled REML noise variance <c>σ̂²(λ) = RSS_pen / (n − r)</c>.</summary>
    public double ProfiledVariance(double lambda) => Solve(lambda).rssPen / (_n - _penalty.Nullity);

    /// <summary>
    /// Profiled REML log-evidence <c>ℓ_p(λ)</c> (σ² profiled out, β integrated under the improper penalty prior):
    /// <c>−½[(n−r)·log RSS_pen − (ν−r)·log λ + log|A|]</c>, up to a λ-independent constant. The marginal objective
    /// that selects the smoothing — its interior maximum is the REML-optimal λ, and it is the target for the
    /// marginal-evidence λ sampler.
    /// </summary>
    public double RemlLogEvidence(double lambda) => RemlLogEvidence(lambda, out _);

    // The evidence and the solve behind it, which the collapsed sampler keeps for the accepted λ's conditional draw.
    internal double RemlLogEvidence(double lambda, out PenalizedSolve solve)
    {
        if (!(lambda > 0.0)) throw new ArgumentOutOfRangeException(nameof(lambda), "Smoothing λ must be positive.");
        solve = SolveFactored(lambda);
        int r = _penalty.Nullity;
        double minusTwoEll = (_n - r) * Math.Log(solve.RssPen) - (_nu - r) * Math.Log(lambda) + solve.Factor.LogDet;
        return -0.5 * minusTwoEll;
    }

    /// <summary>
    /// Effective degrees of freedom <c>edf(λ) = tr[(ZᵀZ + λP)⁻¹ZᵀZ]</c> (Hodges &amp; Sargent 2001): from the number of
    /// coefficients the data identify as λ → 0 down to the nullity r as λ → ∞. Read off the penalty spectrum in O(ν).
    /// </summary>
    public double EffectiveDegreesOfFreedom(double lambda)
    {
        if (!(lambda > 0.0)) throw new ArgumentOutOfRangeException(nameof(lambda), "Smoothing λ must be positive.");
        double edf = 0.0;
        foreach (double mu in _spectrum.Value) edf += (1.0 - mu) / (1.0 - mu + lambda * mu);
        return edf;
    }

    /// <summary>
    /// <c>|d edf / d ln λ| = Σ_j λμ_j(1 − μ_j)/(1 − μ_j + λμ_j)²</c>: the density on u = ln λ that a uniform prior on
    /// <see cref="EffectiveDegreesOfFreedom"/> induces, up to the constant edf range (a uniform shrinkage prior in the
    /// sense of Natarajan &amp; Kass 2000). It integrates to the edf range and decays like 1/λ, so it makes the
    /// marginal-evidence target proper; each term is positive, so it stays accurate at extreme λ.
    /// </summary>
    public double EdfSlope(double lambda)
    {
        if (!(lambda > 0.0)) throw new ArgumentOutOfRangeException(nameof(lambda), "Smoothing λ must be positive.");
        double slope = 0.0;
        foreach (double mu in _spectrum.Value)
        {
            double denom = 1.0 - mu + lambda * mu;
            slope += lambda * mu * (1.0 - mu) / (denom * denom);
        }
        return slope;
    }

    // Generalized eigenvalues μ_j ∈ [0, 1] of P relative to M = ZᵀZ + P. With W = L⁻ᵀD^{−½} from M's LDLᵀ factor
    // (WWᵀ = M⁻¹), they are the eigenvalues of WᵀPW, and in that basis A(λ) = ZᵀZ + λP is diag(1 − μ + λμ). The r
    // smallest belong to the penalty null space and are set to exactly zero.
    private double[] PenaltySpectrum()
    {
        var band = (double[,])_gram.Clone();          // ZᵀZ
        _penalty.AccumulateInto(band, _nu, 1.0);      // + P
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        if (chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
            throw new InvalidOperationException($"ZᵀZ + P is singular: the design does not identify the penalty null space ({chol.Status}).");

        var w = new double[_nu, _nu];                  // columns of W = L⁻ᵀD^{−½}
        var unit = new double[_nu];
        var col = new double[_nu];
        for (int k = 0; k < _nu; k++)
        {
            Array.Clear(unit);
            unit[k] = 1.0;
            chol.SampleInnovation(unit, col);
            for (int i = 0; i < _nu; i++) w[i, k] = col[i];
        }

        var p = new double[_bw + 1, _nu];              // P in lower-band storage p[d, j] = P(j+d, j)
        _penalty.AccumulateInto(p, _nu, 1.0);
        var pw = new double[_nu, _nu];
        for (int j = 0; j < _nu; j++)
            for (int d = 0; d <= _bw && j + d < _nu; d++)
            {
                double v = p[d, j];
                if (v == 0.0) continue;
                for (int k = 0; k < _nu; k++)
                {
                    pw[j + d, k] += v * w[j, k];
                    if (d > 0) pw[j, k] += v * w[j + d, k];
                }
            }

        var c = new double[_nu, _nu];                  // WᵀPW, symmetrized
        for (int a = 0; a < _nu; a++)
            for (int b = 0; b <= a; b++)
            {
                double s = 0.0, t = 0.0;
                for (int i = 0; i < _nu; i++) { s += w[i, a] * pw[i, b]; t += w[i, b] * pw[i, a]; }
                c[a, b] = c[b, a] = 0.5 * (s + t);
            }

        double[] mu = DenseEigen.DecomposeSymmetric(c).Eigenvalues;
        Array.Sort(mu);
        int r = _penalty.Nullity;
        bool identified = false;
        for (int i = 0; i < mu.Length; i++)
        {
            mu[i] = i < r ? 0.0 : Math.Clamp(mu[i], 0.0, 1.0);
            identified |= mu[i] > 0.0 && mu[i] < 1.0;
        }
        if (!identified)
            throw new InvalidOperationException("No penalized direction is identified by the data: edf is constant in λ.");
        return mu;
    }
}

/// <summary>A penalized solve at one λ: β̂(λ), the penalized residual sum of squares and the factor of A = ZᵀZ + λP.</summary>
internal sealed record PenalizedSolve(double[] Beta, double RssPen, BandCholesky Factor);
