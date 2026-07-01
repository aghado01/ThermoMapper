namespace Maths.Regression.Spline.Baps;

/// <summary>
/// A roughness penalty that assembles its banded precision <c>λ·P</c> into a (shared) Gram band — the swappable
/// smoother behind the penalized P-spline. Implemented by the 1-D <see cref="DifferencePenalty"/> and the
/// tensor-product <see cref="TensorPenalty"/>; consumed by <see cref="PenalizedSpline"/> (and BAPS through it).
/// Keeping <see cref="Bandwidth"/> and <see cref="Nullity"/> distinct is what lets the tensor case work: there a
/// difference order of 2 still gives a flattened bandwidth of <c>order·ν_inner</c> and a null space of
/// <c>order_x·order_y</c>.
/// </summary>
public interface IBandPenalty
{
    /// <summary>Half-bandwidth of P in the (flattened) coefficient ordering.</summary>
    int Bandwidth { get; }

    /// <summary>Dimension of P's null space — the unpenalized polynomial degrees of freedom (the REML fixed effects).</summary>
    int Nullity { get; }

    /// <summary>Adds <c>λ·P</c> into LAPACK lower-band storage <c>band[d, j] = A(j+d, j)</c> over <paramref name="dim"/> coefficients.</summary>
    void AccumulateInto(double[,] band, int dim, double lambda);
}
