using System;
using Maths.Rng;

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

/// <summary>
/// A penalty whose components are smoothed separately, <c>P(λ) = Σⱼ λⱼPⱼ</c>: the anisotropic tensor-product penalty,
/// one component per axis. <see cref="IBandPenalty.AccumulateInto(double[,], int, double)"/> is its isotropic case, every
/// λⱼ = λ. The components share an eigenbasis (<see cref="Spectrum"/>), so P(λ)'s pseudo-determinant costs O(ν) at any λ;
/// that is what makes the per-component smoothing selectable by REML and its variances samplable (Bach &amp; Klein 2025).
/// </summary>
public interface IAnisotropicPenalty : IBandPenalty
{
    /// <summary>The number of components, each with its own smoothing.</summary>
    int Components { get; }

    /// <summary>Adds <c>Σⱼ λⱼPⱼ</c> into band storage, one smoothing per component.</summary>
    void AccumulateComponents(double[,] band, int dim, ReadOnlySpan<double> lambdas);

    /// <summary>Component j's roughness <c>βᵀPⱼβ</c>.</summary>
    double ComponentRoughness(double[] beta, int component);

    /// <summary>The components' joint spectrum.</summary>
    PenaltySpectrum Spectrum { get; }

    /// <summary>One draw from the improper prior <c>N(0, (Σⱼ wⱼPⱼ)⁺)</c>, which has no mass on the null space, at precisions
    /// <paramref name="weights"/> = 1/τⱼ²: what prior scaling (<see cref="WeibullVariancePrior.ScaledTo"/>) draws.</summary>
    void SamplePrior(ReadOnlySpan<double> weights, Xoshiro256PlusPlus rng, Span<double> beta);
}
