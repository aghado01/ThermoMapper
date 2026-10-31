using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// An inner product on an output index, ⟨u, v⟩ = uᵀQv for a symmetric positive-definite Q (the declared output-side
    /// inner product of direct functional BASS, F-X6). A residual with covariance σ²(I ⊗ R), independent across cases and
    /// correlated along the index, has Q = R⁻¹; independent noise has Q = I (<see cref="IdentityInnerProduct"/>). It acts
    /// on stacked vectors, k cases of length m end to end, as I_k ⊗ Q, so one member gives the output-side Gram entries
    /// h_iᵀQh_j (k = 1) and the energy Σᵢ yᵢᵀQyᵢ (k = n). A banded Q keeps each product linear in m.
    /// </summary>
    public interface IInnerProduct
    {
        /// <summary>Whether Q = I, which the contractions take the unweighted arithmetic path for.</summary>
        bool IsIdentity { get; }

        /// <summary>The output length m that Q acts on; 0 for the identity, which acts on any length.</summary>
        int Dimension { get; }

        /// <summary>ln|Q|, whose negative n times is the collapsed marginal's log|R| term per case; 0 for the identity.</summary>
        double LogDeterminant { get; }

        /// <summary>Σᵢ uᵢᵀQvᵢ over the k stacked vectors of <paramref name="u"/> and <paramref name="v"/>, whose lengths are equal
        /// and a multiple of <see cref="Dimension"/>.</summary>
        double Dot(ReadOnlySpan<double> u, ReadOnlySpan<double> v);
    }

    /// <summary>Q = I: independent noise along the output index. Its product is <see cref="PairwiseSum.Dot"/>, so a declared
    /// identity does the arithmetic the unweighted contractions do.</summary>
    public sealed class IdentityInnerProduct : IInnerProduct
    {
        /// <summary>The identity.</summary>
        public static IdentityInnerProduct Instance { get; } = new();

        private IdentityInnerProduct() { }

        /// <inheritdoc/>
        public bool IsIdentity => true;

        /// <inheritdoc/>
        public int Dimension => 0;

        /// <inheritdoc/>
        public double LogDeterminant => 0.0;

        /// <inheritdoc/>
        public double Dot(ReadOnlySpan<double> u, ReadOnlySpan<double> v) => PairwiseSum.Dot(u, v);
    }
}
