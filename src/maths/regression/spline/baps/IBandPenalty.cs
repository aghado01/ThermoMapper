namespace Maths.Regression.Spline.Baps;

/// <summary>
/// A roughness penalty that assembles its banded precision <c>λ·P</c> into a (shared) Gram band — the swappable
/// smoother behind the penalized P-spline. Implemented by the 1-D <see cref="DifferencePenalty"/> and the
/// tensor-product <see cref="TensorPenalty"/>; consumed by <see cref="PenalizedSpline"/>, <see cref="Baps"/> and
/// <see cref="HierarchicalPSpline"/>.
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

    /// <summary>The roughness quadratic <c>βᵀPβ</c>: the scale of the smoothing variance's inverse-gamma conditional.</summary>
    double Roughness(double[] beta);
}

/// <summary>
/// A penalty that is a sum of local terms, <c>P = Dᵀ diag(w) D</c>, with one weight per row of D: what locally adaptive
/// smoothing reweights (<see cref="BapsScheme.LocallyAdaptiveGibbs"/>). D must have full row rank, so the adaptive
/// prior's normalizer factorizes as <c>(Π wᵢ)·|DDᵀ|</c> and each weight's multiplier has a Gamma full conditional; the
/// 1-D difference penalty qualifies, and the tensor penalty, whose stacked rows outnumber its rank, does not.
/// </summary>
public interface ILocallyAdaptivePenalty : IBandPenalty
{
    /// <summary>The number of local terms over <paramref name="dim"/> coefficients: one weight each.</summary>
    int DifferenceCount(int dim);

    /// <summary>Adds <c>Dᵀ diag(weights) D</c> into band storage.</summary>
    void AccumulateInto(double[,] band, int dim, double[] weights);

    /// <summary>Writes each local term's square <c>(Dβ)ᵢ²</c> into <paramref name="into"/>.</summary>
    void SquaredDifferencesInto(double[] beta, double[] into);
}
