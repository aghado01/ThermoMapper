using System;

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
/// Isotropic: a single λ scales both directions. Anisotropic smoothing (separate λ_x, λ_y) is the natural
/// refinement — the assembly already separates the two terms.
/// </remarks>
public sealed class TensorPenalty : IBandPenalty
{
    private readonly DifferencePenalty _px;
    private readonly DifferencePenalty _py;
    private readonly int _nuX;
    private readonly int _nuY;

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
}
