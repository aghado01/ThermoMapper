using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// Tensor-product difference penalty for bivariate P-splines (He, Yang &amp; Kang 2024; Eilers–Marx–Currie): the
/// roughness penalty <c>P = D_xᵀD_x ⊗ I_y + I_x ⊗ D_yᵀD_y</c> over coefficients flattened as <c>j·ν_y + k</c>
/// (x outer, y inner). In that ordering it is just two strided 1-D difference penalties — the y-penalty runs
/// contiguously within each x-block, the x-penalty couples blocks at stride <c>ν_y</c> — so it assembles straight
/// into the flattened Gram band via <see cref="DifferencePenalty.AccumulateStrided"/> and the whole system stays
/// a banded solve (flattened half-bandwidth <c>order_x·ν_y</c>; put the smaller dimension inner to keep it small).
/// </summary>
/// <remarks>
/// <see cref="IBandPenalty.AccumulateInto"/> is the isotropic penalty: a single λ scales both directions. As an
/// <see cref="IAnisotropicPenalty"/> it is also the anisotropic one, <c>λ_x(D_xᵀD_x ⊗ I) + λ_y(I ⊗ D_yᵀD_y)</c>, with
/// component 0 the x-penalty and component 1 the y-penalty.
/// </remarks>
public sealed class TensorPenalty : IAnisotropicPenalty
{
    private readonly DifferencePenalty _px;
    private readonly DifferencePenalty _py;
    private readonly int _nuX;
    private readonly int _nuY;
    private readonly Lazy<Eigenbasis> _basis;

    /// <param name="nuX">Number of x-direction basis functions (outer index).</param>
    /// <param name="nuY">Number of y-direction basis functions (inner index).</param>
    /// <param name="orderX">x-direction difference order.</param>
    /// <param name="orderY">y-direction difference order.</param>
    public TensorPenalty(int nuX, int nuY, int orderX = 2, int orderY = 2)
    {
        if (nuX <= orderX) throw new ArgumentOutOfRangeException(nameof(nuX), "Need more x basis functions than the x penalty order.");
        if (nuY <= orderY) throw new ArgumentOutOfRangeException(nameof(nuY), "Need more y basis functions than the y penalty order.");
        _nuX = nuX;
        _nuY = nuY;
        _px = new DifferencePenalty(orderX);
        _py = new DifferencePenalty(orderY);
        _basis = new Lazy<Eigenbasis>(BuildBasis);
    }

    /// <summary>Two components: the x-penalty (0) and the y-penalty (1).</summary>
    public int Components => 2;

    /// <summary>
    /// The joint spectrum. <c>D_xᵀD_x ⊗ I</c> and <c>I ⊗ D_yᵀD_y</c> are both diagonal in the Kronecker eigenbasis
    /// <c>Q_x ⊗ Q_y</c>, so mode <c>l = i·ν_y + k</c> carries the x-penalty's eigenvalue <c>a_i</c> and the y-penalty's
    /// <c>b_k</c>, each 1-D penalty's null modes set to 0 (Bach &amp; Klein 2025, Theorem 3.1).
    /// </summary>
    public PenaltySpectrum Spectrum => _basis.Value.Spectrum;

    /// <summary>
    /// One draw from <c>N(0, (w_x D_xᵀD_x⊗I + w_y I⊗D_yᵀD_y)⁺)</c>: <c>(Q_x ⊗ Q_y)z̃</c> with <c>z̃ ~ N(0, 1/(w_x a_i + w_y b_k))</c>
    /// on each non-null mode (i, k) and 0 on the null modes, formed as <c>Q_x Z̃ Q_yᵀ</c> in O(ν_xν_y(ν_x + ν_y)).
    /// </summary>
    public void SamplePrior(ReadOnlySpan<double> weights, Xoshiro256PlusPlus rng, Span<double> beta)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (weights.Length != 2) throw new ArgumentException("Need [w_x, w_y].", nameof(weights));
        if (beta.Length != _nuX * _nuY) throw new ArgumentException($"beta must have nuX·nuY = {_nuX * _nuY} entries.", nameof(beta));
        Eigenbasis e = _basis.Value;
        var t = new double[_nuX, _nuY];   // Z̃ Q_yᵀ
        var z = new double[_nuY];
        for (int i = 0; i < _nuX; i++)
        {
            for (int k = 0; k < _nuY; k++)
            {
                double s = weights[0] * e.A[i] + weights[1] * e.B[k];
                z[k] = s > 0.0 ? Normal.SampleStandard(rng) / Math.Sqrt(s) : 0.0;
            }
            for (int k2 = 0; k2 < _nuY; k2++)
            {
                double sum = 0.0;
                for (int k = 0; k < _nuY; k++) sum += z[k] * e.Qy[k][k2];
                t[i, k2] = sum;
            }
        }
        for (int i2 = 0; i2 < _nuX; i2++)
            for (int k2 = 0; k2 < _nuY; k2++)
            {
                double sum = 0.0;
                for (int i = 0; i < _nuX; i++) sum += e.Qx[i][i2] * t[i, k2];
                beta[i2 * _nuY + k2] = sum;
            }
    }

    /// <summary>Flattened half-bandwidth: <c>order_x·ν_y</c> (the x-penalty's block stride; dominates order_y).</summary>
    public int Bandwidth => Math.Max(_px.Order * _nuY, _py.Order);

    /// <summary>Null space = polynomials of degree &lt; order in each factor: <c>order_x · order_y</c>.</summary>
    public int Nullity => _px.Order * _py.Order;

    /// <summary>Adds <c>λ·(D_xᵀD_x ⊗ I_y + I_x ⊗ D_yᵀD_y)</c> into the flattened lower-band storage.</summary>
    public void AccumulateInto(double[,] band, int dim, double lambda)
    {
        if (dim != _nuX * _nuY)
            throw new ArgumentException($"dim must equal nuX·nuY = {_nuX * _nuY}.", nameof(dim));

        for (int j = 0; j < _nuX; j++)                              // I_x ⊗ D_yᵀD_y — within each x-block
            _py.AccumulateStrided(band, _nuY, stride: 1, offset: j * _nuY, lambda);

        for (int k = 0; k < _nuY; k++)                             // D_xᵀD_x ⊗ I_y — across blocks at stride ν_y
            _px.AccumulateStrided(band, _nuX, stride: _nuY, offset: k, lambda);
    }

    /// <summary>Adds <c>λ_x(D_xᵀD_x ⊗ I_y) + λ_y(I_x ⊗ D_yᵀD_y)</c>, with <c>lambdas = [λ_x, λ_y]</c>.</summary>
    public void AccumulateComponents(double[,] band, int dim, ReadOnlySpan<double> lambdas)
    {
        if (dim != _nuX * _nuY)
            throw new ArgumentException($"dim must equal nuX·nuY = {_nuX * _nuY}.", nameof(dim));
        if (lambdas.Length != 2) throw new ArgumentException("Need [λ_x, λ_y].", nameof(lambdas));
        double lambdaX = lambdas[0], lambdaY = lambdas[1];

        for (int j = 0; j < _nuX; j++)                              // λ_y·(I_x ⊗ D_yᵀD_y)
            _py.AccumulateStrided(band, _nuY, stride: 1, offset: j * _nuY, lambdaY);

        for (int k = 0; k < _nuY; k++)                             // λ_x·(D_xᵀD_x ⊗ I_y)
            _px.AccumulateStrided(band, _nuX, stride: _nuY, offset: k, lambdaX);
    }

    /// <summary>The x-roughness <c>βᵀ(D_xᵀD_x ⊗ I)β</c> (component 0) or the y-roughness <c>βᵀ(I ⊗ D_yᵀD_y)β</c> (1).</summary>
    public double ComponentRoughness(double[] beta, int component)
    {
        ArgumentNullException.ThrowIfNull(beta);
        if (beta.Length != _nuX * _nuY)
            throw new ArgumentException($"beta must have nuX·nuY = {_nuX * _nuY} entries.", nameof(beta));
        double sum = 0.0;
        switch (component)
        {
            case 0:
                for (int k = 0; k < _nuY; k++) sum += _px.RoughnessStrided(beta, _nuX, stride: _nuY, offset: k);
                return sum;
            case 1:
                for (int j = 0; j < _nuX; j++) sum += _py.RoughnessStrided(beta, _nuY, stride: 1, offset: j * _nuY);
                return sum;
            default:
                throw new ArgumentOutOfRangeException(nameof(component), "The tensor penalty has components 0 (x) and 1 (y).");
        }
    }

    /// <summary><c>βᵀ(D_xᵀD_x ⊗ I_y + I_x ⊗ D_yᵀD_y)β</c>: each x-block's y-roughness plus each y-column's x-roughness.</summary>
    public double Roughness(double[] beta)
    {
        ArgumentNullException.ThrowIfNull(beta);
        if (beta.Length != _nuX * _nuY)
            throw new ArgumentException($"beta must have nuX·nuY = {_nuX * _nuY} entries.", nameof(beta));
        double sum = 0.0;
        for (int j = 0; j < _nuX; j++) sum += _py.RoughnessStrided(beta, _nuY, stride: 1, offset: j * _nuY);
        for (int k = 0; k < _nuY; k++) sum += _px.RoughnessStrided(beta, _nuX, stride: _nuY, offset: k);
        return sum;
    }

    // Each 1-D penalty's eigendecomposition, its eigenvectors as rows: the joint spectrum and the Kronecker eigenbasis.
    private sealed record Eigenbasis(PenaltySpectrum Spectrum, double[] A, double[] B, double[][] Qx, double[][] Qy);

    private Eigenbasis BuildBasis()
    {
        EigenResult ex = DenseEigen.DecomposeSymmetric(_px.ToDense(_nuX));
        EigenResult ey = DenseEigen.DecomposeSymmetric(_py.ToDense(_nuY));
        double[] a = ClampNull(ex.Eigenvalues);
        double[] b = ClampNull(ey.Eigenvalues);
        var gx = new double[_nuX * _nuY];
        var gy = new double[_nuX * _nuY];
        for (int i = 0; i < _nuX; i++)
            for (int k = 0; k < _nuY; k++) { gx[i * _nuY + k] = a[i]; gy[i * _nuY + k] = b[k]; }
        return new Eigenbasis(new PenaltySpectrum([gx, gy]), a, b, ex.Eigenvectors, ey.Eigenvectors);
    }

    // A difference penalty's null modes come out of the eigensolver at rounding level; they are set to exactly 0.
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
