using System;
using Maths.LinAlg;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// Anisotropic tensor-product P-spline (He, Yang &amp; Kang; Currie–Durban–Eilers): the bivariate smoother with
/// <i>separate</i> smoothing per axis, <c>A = ZᵀZ + λ_x(D_xᵀD_x⊗I) + λ_y(I⊗D_yᵀD_y)</c>, chosen by 2-D REML
/// (<see cref="RemlOptimum"/>). The summed penalty's pseudo-determinant does not factor into the two smoothings, so the
/// smoothing variances have no conjugate Gibbs step; but P_x and P_y share the Kronecker eigenbasis, so the combined
/// penalty's eigenvalues are the pairwise sums <c>λ_x a_i + λ_y b_j</c> of the two small 1-D penalty spectra
/// (<see cref="TensorPenalty.Spectrum"/>). That makes the REML evidence
/// <c>−½[(n−r)·log RSS_pen − Σ_{ij} log(λ_x a_i + λ_y b_j) + log|A|]</c> computable from the spectrum plus the banded
/// factor of A (<see cref="BandCholesky"/>), and it is what lets <see cref="Baps"/>'s anisotropic scheme sample the two
/// smoothing variances by Metropolis (Bach &amp; Klein 2025).
/// </summary>
public sealed class AnisotropicTensorSpline
{
    private readonly BandedDesign _design;
    private readonly TensorPenalty _penalty;
    private readonly double[,] _gram;        // ZᵀZ band — constant
    private readonly double[] _zty;          // Zᵀy — constant
    private readonly int _n;
    private readonly int _nu;
    private readonly int _bw;
    private readonly int _r;                 // combined null-space dim = (#null a_i)·(#null b_j)
    private readonly double _yy;

    public AnisotropicTensorSpline(TensorDesign td, double[] y, int orderX = 2, int orderY = 2)
    {
        ArgumentNullException.ThrowIfNull(td);
        ArgumentNullException.ThrowIfNull(y);
        _design = new BandedDesign(td.Design);
        _n = _design.Rows;
        _nu = _design.Dimension;
        if (y.Length != _n) throw new ArgumentException("Response length must match the design row count.", nameof(y));

        _penalty = new TensorPenalty(td.NuX, td.NuY, orderX, orderY);
        _bw = Math.Max(_design.Bandwidth, _penalty.Bandwidth);
        _r = _penalty.Spectrum.Nullity;

        _gram = new double[_bw + 1, _nu];
        _design.Accumulate(null, _gram, null, null);
        _zty = new double[_nu];
        _design.AccumulateRhs(null, y, _zty);
        double yy = 0.0;
        foreach (double v in y) yy += v * v;
        _yy = yy;
    }

    /// <summary>The anisotropic penalty: component 0 is x, component 1 is y.</summary>
    public TensorPenalty Penalty => _penalty;

    /// <summary>Penalized coefficients <c>β̂(λ_x, λ_y) = A⁻¹ Zᵀy</c> — the banded anisotropic Reinsch solve.</summary>
    public double[] Coefficients(double lambdaX, double lambdaY)
    {
        var chol = Factor(lambdaX, lambdaY, out _);
        return chol.Solve(_zty);
    }

    /// <summary>
    /// Profiled REML log-evidence over both smoothing parameters. Its 2-D maximum is the REML-optimal
    /// <c>(λ_x, λ_y)</c> (<see cref="RemlOptimum"/>).
    /// </summary>
    public double RemlLogEvidence(double lambdaX, double lambdaY)
    {
        if (!(lambdaX > 0.0 && lambdaY > 0.0))
            throw new ArgumentOutOfRangeException(nameof(lambdaX), "Smoothing parameters must be positive.");

        BandCholesky chol = Factor(lambdaX, lambdaY, out double logDetA);
        double[] beta = chol.Solve(_zty);
        double bBeta = 0.0;
        for (int p = 0; p < _nu; p++) bBeta += _zty[p] * beta[p];
        double rssPen = _yy - bBeta;

        double logPriorDet = _penalty.Spectrum.LogPseudoDeterminant([lambdaX, lambdaY]);   // Σ log(λ_x a_i + λ_y b_j)

        double minusTwoEll = (_n - _r) * Math.Log(rssPen) - logPriorDet + logDetA;
        return -0.5 * minusTwoEll;
    }

    /// <summary>
    /// The REML-optimal <c>(λ_x, λ_y)</c>: the maximum of <see cref="RemlLogEvidence"/> over <c>(ln λ_x, ln λ_y)</c>, by
    /// Nelder–Mead from each axis's balance point <c>λ_ref = tr(ZᵀZ)/tr(P_axis)</c>, within <c>λ_ref·[10⁻⁸, 10⁸]</c> per
    /// axis. The evidence levels off as an axis's λ grows when the data have no curvature along that axis (they lie in its
    /// null space), so an axis that ends at its box's edge is flagged rather than reported as an interior optimum. A trial
    /// point whose system is not numerically positive definite counts as infeasible.
    /// </summary>
    /// <param name="tolerance">Convergence on ln λ, and on the evidence relative to its size.</param>
    /// <param name="maxEvaluations">The evidence evaluations allowed.</param>
    public AnisotropicReml RemlOptimum(double tolerance = 1e-8, int maxEvaluations = 400)
    {
        double traceGram = 0.0;
        for (int j = 0; j < _nu; j++) traceGram += _gram[0, j];
        double refX = Math.Log(traceGram / ComponentTrace(0)), refY = Math.Log(traceGram / ComponentTrace(1));
        const double span = 18.420680743952367;   // ln 10⁸
        double[] lower = [refX - span, refY - span];
        double[] upper = [refX + span, refY + span];
        NelderMead.Result r = NelderMead.Minimize(
            u =>
            {
                try { return -RemlLogEvidence(Math.Exp(u[0]), Math.Exp(u[1])); }
                catch (InvalidOperationException) { return double.PositiveInfinity; }   // not numerically positive definite
            },
            [refX, refY], [1.0, 1.0], lower, upper, tolerance, maxEvaluations);
        const double edge = 1e-6;
        return new AnisotropicReml(Math.Exp(r.Point[0]), Math.Exp(r.Point[1]), -r.Value, r.Evaluations, r.Converged,
                                   r.Point[0] - lower[0] < edge || upper[0] - r.Point[0] < edge,
                                   r.Point[1] - lower[1] < edge || upper[1] - r.Point[1] < edge);
    }

    // tr(P_component): the diagonal of that component alone.
    private double ComponentTrace(int component)
    {
        var lambdas = new double[2];
        lambdas[component] = 1.0;
        var band = new double[_penalty.Bandwidth + 1, _nu];
        _penalty.AccumulateComponents(band, _nu, lambdas);
        double trace = 0.0;
        for (int j = 0; j < _nu; j++) trace += band[0, j];
        return trace;
    }

    private BandCholesky Factor(double lambdaX, double lambdaY, out double logDetA)
    {
        var band = (double[,])_gram.Clone();                               // ZᵀZ
        _penalty.AccumulateComponents(band, _nu, [lambdaX, lambdaY]);      // + λx·Px + λy·Py
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        if (chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
            throw new InvalidOperationException(
                $"ZᵀZ + λx·Px + λy·Py is not numerically positive definite at (λx, λy) = ({lambdaX}, {lambdaY}) ({chol.Status}).");
        logDetA = chol.LogDet;
        return chol;
    }
}

/// <summary>
/// The REML-optimal smoothing of an anisotropic tensor spline: <c>(λ_x, λ_y)</c>, the evidence there, the evaluations
/// spent, whether the search converged, and whether each axis ended at its search box's edge (for the upper edge: no
/// curvature along that axis, so its smoothing grows without bound).
/// </summary>
public sealed record AnisotropicReml(double LambdaX, double LambdaY, double LogEvidence, int Evaluations, bool Converged,
                                     bool LambdaXAtBound, bool LambdaYAtBound);
