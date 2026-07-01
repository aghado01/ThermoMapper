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
        double yy = 0.0;
        for (int i = 0; i < _n; i++) yy += y[i] * y[i];
        _yy = yy;
    }

    /// <summary>Number of basis functions ν.</summary>
    public int Dimension => _nu;

    /// <summary>Penalty null-space dimension r (= penalty order) — the unpenalized polynomial degrees of freedom.</summary>
    public int Nullity => _penalty.Nullity;

    // (β̂, penalized RSS, log|A|) at smoothing λ; A = ZᵀZ + λP factored once and shared by every read-out.
    private (double[] beta, double rssPen, double logDetA) Solve(double lambda)
    {
        var band = new double[_bw + 1, _nu];
        var b = new double[_nu];
        _design.Accumulate(null, band, _y, b);       // ZᵀZ (band) and Zᵀy
        _penalty.AccumulateInto(band, _nu, lambda);   // + λ·DᵀD
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        double[] beta = chol.Solve(b);

        double bBeta = 0.0;                            // RSS_pen = yᵀy − β̂ᵀZᵀy
        for (int p = 0; p < _nu; p++) bBeta += b[p] * beta[p];
        return (beta, _yy - bBeta, chol.LogDet);
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
    public double RemlLogEvidence(double lambda)
    {
        if (!(lambda > 0.0)) throw new ArgumentOutOfRangeException(nameof(lambda), "Smoothing λ must be positive.");
        var (_, rssPen, logDetA) = Solve(lambda);
        int r = _penalty.Nullity;
        double minusTwoEll = (_n - r) * Math.Log(rssPen) - (_nu - r) * Math.Log(lambda) + logDetA;
        return -0.5 * minusTwoEll;
    }
}
