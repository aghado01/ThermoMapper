using System;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Standard <see cref="IResponsibilityConstraint"/> with linearly annealed λ.
    /// See docs/gmm.md for the semi-supervised path.
    /// </summary>
    public sealed class AnnealedSoftConstraint : IResponsibilityConstraint
    {
        private readonly double[,] _sik;   // N×K, row-normalised on demand
        private readonly double _lambdaStart;
        private readonly double _lambdaEnd;
        private readonly int? _lambdaHorizon;

        /// <param name="confidenceMatrix">
        /// N×K matrix where sᵢₖ ∈ [0, 1] is the supervisor's confidence that
        /// point i belongs to component k. Copied internally.
        /// </param>
        /// <param name="normalizeRows">
        /// When true, each row is normalised to sum to 1 at construction time.
        /// Zero-sum rows are left as-is; points with a zero-sum row are not blended.
        /// </param>
        /// <param name="lambdaStart">λ at iteration 0. Default 0.8.</param>
        /// <param name="lambdaEnd">λ at iteration <c>lambdaHorizon − 1</c>. Default 0.0.</param>
        /// <param name="lambdaHorizon">
        /// Iteration count over which λ decays from <c>lambdaStart</c> to
        /// <c>lambdaEnd</c>; after this many iterations λ is clamped at
        /// <c>lambdaEnd</c>. Null falls back to <c>maxIterations</c> supplied at
        /// <see cref="Apply"/> time, which couples decay to the iteration cap.
        /// </param>
        public AnnealedSoftConstraint(
            double[,] confidenceMatrix,
            bool normalizeRows = false,
            double lambdaStart = 0.8,
            double lambdaEnd = 0.0,
            int? lambdaHorizon = null)
        {
            if (lambdaStart < 0.0 || lambdaStart > 1.0)
                throw new ArgumentOutOfRangeException(nameof(lambdaStart), "Must be in [0, 1].");
            if (lambdaEnd < 0.0 || lambdaEnd > 1.0)
                throw new ArgumentOutOfRangeException(nameof(lambdaEnd), "Must be in [0, 1].");
            if (lambdaHorizon is int h && h < 1)
                throw new ArgumentOutOfRangeException(nameof(lambdaHorizon), "Must be ≥ 1 when supplied.");

            int n = confidenceMatrix.GetLength(0);
            int k = confidenceMatrix.GetLength(1);
            _sik = new double[n, k];

            for (int i = 0; i < n; i++)
            {
                double rowSum = 0.0;
                for (int ki = 0; ki < k; ki++) rowSum += confidenceMatrix[i, ki];

                double scale = (normalizeRows && rowSum > 1e-12) ? 1.0 / rowSum : 1.0;
                for (int ki = 0; ki < k; ki++)
                    _sik[i, ki] = confidenceMatrix[i, ki] * scale;
            }

            _lambdaStart = lambdaStart;
            _lambdaEnd = lambdaEnd;
            _lambdaHorizon = lambdaHorizon;
        }

        /// <inheritdoc/>
        public void Apply(double[,] responsibilities, int n, int k, int iteration, int maxIterations)
        {
            if (n != _sik.GetLength(0) || k != _sik.GetLength(1))
                throw new InvalidOperationException(
                    $"Constraint matrix is {_sik.GetLength(0)}×{_sik.GetLength(1)} but " +
                    $"E-step produced an {n}×{k} responsibility matrix.");

            double lambda = ComputeLambda(iteration, maxIterations);
            if (lambda < 1e-12) return;
            double oneMinusLambda = 1.0 - lambda;

            for (int i = 0; i < n; i++)
            {
                // Supervisor abstains on zero-sum rows; leave E-step result intact.
                double sikRowSum = 0.0;
                for (int ki = 0; ki < k; ki++) sikRowSum += _sik[i, ki];
                if (sikRowSum < 1e-12) continue;

                double blendedSum = 0.0;
                for (int ki = 0; ki < k; ki++)
                {
                    double blended = oneMinusLambda * responsibilities[i, ki]
                                   + lambda * _sik[i, ki];
                    responsibilities[i, ki] = blended;
                    blendedSum += blended;
                }

                if (blendedSum > 1e-12)
                {
                    double invSum = 1.0 / blendedSum;
                    for (int ki = 0; ki < k; ki++)
                        responsibilities[i, ki] *= invSum;
                }
            }
        }

        /// <summary>λ at the given iteration, for diagnostics and schedule tuning.</summary>
        public double LambdaAt(int iteration, int maxIterations) =>
            ComputeLambda(iteration, maxIterations);

        private double ComputeLambda(int iteration, int maxIterations)
        {
            int horizon = _lambdaHorizon ?? maxIterations;
            if (horizon <= 1) return _lambdaEnd;
            double t = Math.Min(1.0, (double)iteration / (horizon - 1));
            return _lambdaStart + (_lambdaEnd - _lambdaStart) * t;
        }
    }
}
