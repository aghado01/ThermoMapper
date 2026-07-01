using System;
using Maths.LinAlg;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// Anisotropic tensor-product P-spline (He, Yang &amp; Kang; Currie–Durban–Eilers): the bivariate smoother with
/// <i>separate</i> smoothing per axis, <c>A = ZᵀZ + λ_x(D_xᵀD_x⊗I) + λ_y(I⊗D_yᵀD_y)</c>, chosen by 2-D REML.
/// Clean Gibbs is unavailable — the summed-penalty pseudo-determinant doesn't factor — but P_x and P_y commute,
/// so the combined penalty's eigenvalues are the pairwise sums <c>λ_x a_i + λ_y b_j</c> of the two small 1-D
/// penalty spectra. That makes the REML evidence
/// <c>−½[(n−r)·log RSS_pen − Σ_{ij} log(λ_x a_i + λ_y b_j) + log|A|]</c> computable from a one-time
/// eigendecomposition of each 1-D penalty (<see cref="DenseEigen"/>) plus the banded factor of A
/// (<see cref="BandCholesky"/>). The same eigen-structure is the Kronecker/Demmler–Reinsch route that would make
/// the whole λ-grid O(ν) per point.
/// </summary>
public sealed class AnisotropicTensorSpline
{
    private readonly BandedDesign _design;
    private readonly DifferencePenalty _px;
    private readonly DifferencePenalty _py;
    private readonly double[,] _gram;        // ZᵀZ band — constant
    private readonly double[] _zty;          // Zᵀy — constant
    private readonly double[] _ax;           // 1-D x-penalty eigenvalues (null modes clamped to 0)
    private readonly double[] _by;           // 1-D y-penalty eigenvalues
    private readonly int _n;
    private readonly int _nu;
    private readonly int _nuX;
    private readonly int _nuY;
    private readonly int _bw;
    private readonly int _r;                 // combined null-space dim = (#null a_i)·(#null b_j)
    private readonly double _yy;

    public AnisotropicTensorSpline(TensorDesign td, double[] y, int orderX = 2, int orderY = 2)
    {
        ArgumentNullException.ThrowIfNull(td);
        ArgumentNullException.ThrowIfNull(y);
        _nuX = td.NuX;
        _nuY = td.NuY;
        _design = new BandedDesign(td.Design);
        _n = _design.Rows;
        _nu = _design.Dimension;
        if (y.Length != _n) throw new ArgumentException("Response length must match the design row count.", nameof(y));

        _px = new DifferencePenalty(orderX);
        _py = new DifferencePenalty(orderY);
        _bw = Math.Max(_design.Bandwidth, Math.Max(orderX * _nuY, orderY));

        _ax = ClampNull(DenseEigen.DecomposeSymmetric(_px.ToDense(_nuX)).Eigenvalues);
        _by = ClampNull(DenseEigen.DecomposeSymmetric(_py.ToDense(_nuY)).Eigenvalues);
        int nullX = 0, nullY = 0;
        foreach (double a in _ax) if (a == 0.0) nullX++;
        foreach (double b in _by) if (b == 0.0) nullY++;
        _r = nullX * nullY;

        _gram = new double[_bw + 1, _nu];
        _design.Accumulate(null, _gram, null, null);
        _zty = new double[_nu];
        _design.AccumulateRhs(null, y, _zty);
        double yy = 0.0;
        foreach (double v in y) yy += v * v;
        _yy = yy;
    }

    /// <summary>Penalized coefficients <c>β̂(λ_x, λ_y) = A⁻¹ Zᵀy</c> — the banded anisotropic Reinsch solve.</summary>
    public double[] Coefficients(double lambdaX, double lambdaY)
    {
        var chol = Factor(lambdaX, lambdaY, out _);
        return chol.Solve(_zty);
    }

    /// <summary>
    /// Profiled REML log-evidence over both smoothing parameters. Its 2-D maximum is the REML-optimal
    /// <c>(λ_x, λ_y)</c> — the basis for anisotropic smoothing selection.
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

        double logPriorDet = 0.0;                 // Σ log(λ_x a_i + λ_y b_j) over the non-null modes
        for (int i = 0; i < _nuX; i++)
            for (int j = 0; j < _nuY; j++)
            {
                double ev = lambdaX * _ax[i] + lambdaY * _by[j];
                if (ev > 0.0) logPriorDet += Math.Log(ev);
            }

        double minusTwoEll = (_n - _r) * Math.Log(rssPen) - logPriorDet + logDetA;
        return -0.5 * minusTwoEll;
    }

    private BandCholesky Factor(double lambdaX, double lambdaY, out double logDetA)
    {
        var band = new double[_bw + 1, _nu];
        Array.Copy(_gram, band, _gram.Length);                       // ZᵀZ
        for (int j = 0; j < _nuX; j++) _py.AccumulateStrided(band, _nuY, 1, j * _nuY, lambdaY);   // I ⊗ D_yᵀD_y
        for (int k = 0; k < _nuY; k++) _px.AccumulateStrided(band, _nuX, _nuY, k, lambdaX);       // D_xᵀD_x ⊗ I
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);
        logDetA = chol.LogDet;
        return chol;
    }

    private static double[] ClampNull(double[] eig)
    {
        double max = 0.0;
        foreach (double e in eig) if (e > max) max = e;
        double tol = 1e-9 * max;
        var clamped = new double[eig.Length];
        for (int i = 0; i < eig.Length; i++) clamped[i] = eig[i] > tol ? eig[i] : 0.0;
        return clamped;
    }
}
