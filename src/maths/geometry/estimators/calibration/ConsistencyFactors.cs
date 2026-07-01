// ============================================================================
// Estimators/Shared/ConsistencyFactors.cs
// ============================================================================
// Calibration constants c_D for scatter estimators under reference distributions.
// Independent of any specific estimator — pass the result as consistencyFactor
// to GeometricMean.ComputeWithScatter or GeometricMedian.ComputeWithScatter.
//
// Usage:
//   double cD = ConsistencyFactors.Gaussian(dim);
//   GeometricMedian.ComputeWithScatter(manifold, data, weights,
//       loc, scatter, opts, consistencyFactor: cD);
//
// Derivation (both cases):
//   The scatter accumulator computes Σ̂ = (c_D / Σᵢwᵢ) · Σᵢ wᵢ vᵢvᵢᵀ.
//   For consistency E[Σ̂] → I we need c_D = D · E[1/r] / E[r]
//   where r = ||v|| is drawn from the distance distribution of the reference.
// ============================================================================
using System;

namespace Maths.Geometry.Estimators.Calibration
{
    public static class ConsistencyFactors
    {
        /// <summary>
        /// Consistency factor c_D for the Weiszfeld scatter under a Gaussian
        /// reference distribution in D dimensions.
        /// <para>
        /// Derivation: r ~ χ_D, so E[r] = √2 · Γ((D+1)/2) / Γ(D/2) and
        /// E[1/r] = Γ((D-1)/2) / (√2 · Γ(D/2)). The Γ terms cancel, giving
        /// c_D = D · E[1/r] / E[r] = D/(D−1).
        /// </para>
        /// D=1 is undefined for this estimator; use MAD for the 1D case.
        /// </summary>
        public static double Gaussian(int dim)
        {
            if (dim <= 1) throw new ArgumentOutOfRangeException(nameof(dim),
                "Weiszfeld scatter is undefined in 1D; use MAD for the scalar case.");
            return (double)dim / (dim - 1);
        }

        /// <summary>
        /// Consistency factor c_D for the Weiszfeld scatter under a spherical
        /// Laplace reference distribution in D dimensions.
        /// <para>
        /// Derivation: r ~ Gamma(D,1), so E[r] = D and E[1/r] = 1/(D−1).
        /// c_D = D · E[1/r] / E[r] = 1/(D−1).
        /// </para>
        /// D=1 is undefined; use MAD for the scalar case.
        /// </summary>
        public static double Laplace(int dim)
        {
            if (dim <= 1) throw new ArgumentOutOfRangeException(nameof(dim),
                "Weiszfeld scatter is undefined in 1D; use MAD for the scalar case.");
            return 1.0 / (dim - 1);
        }
    }
}
