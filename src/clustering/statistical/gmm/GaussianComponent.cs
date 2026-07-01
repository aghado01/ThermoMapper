using System;
using Maths.LinAlg;
using Maths.Distance;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// One Gaussian component (weight π_k, mean μ_k, covariance Σ_k) plus the cached
    /// quantities (Σ⁻¹, ln|Σ|, ln π_k, Cholesky factor) used in the E-step, inference,
    /// and sampling. Scratchpads support an allocation-free M-step.
    /// </summary>
    public sealed class GaussianComponent
    {
        private readonly int _dim;

        // ── Parameters ──────────────────────────────────────────────────────────

        /// <summary>Mixture weight π_k.</summary>
        public double Weight { get; set; }

        /// <summary>Mean vector μ_k. Write elements in place; do not replace the array reference.</summary>
        public double[] Mean { get; }

        /// <summary>Covariance Σ_k. Write elements in place; do not replace the array reference.</summary>
        public double[,] Covariance { get; }

        // ── Cached quantities (updated by UpdateCache) ───────────────────────────

        /// <summary>Σ⁻¹ — used in Mahalanobis distance computation during the E-Step.</summary>
        public double[,] CovarianceInverse { get; }

        /// <summary>
        /// −½(D·ln(2π) + ln|Σ|). Cached to avoid recomputing the log-determinant on every
        /// PDF evaluation. Set by <see cref="UpdateCache"/>.
        /// </summary>
        public double LogNormalizationFactor { get; private set; }

        /// <summary>
        /// ln(max(<see cref="Weight"/>, 1e-300)). Cached so the E-step and any
        /// downstream weighted-log-density evaluation (mode ascent, prediction) avoid
        /// recomputing the log of the mixture weight per data point. Set by
        /// <see cref="UpdateCache"/>; refresh by calling <see cref="UpdateCache"/>
        /// after any direct assignment to <see cref="Weight"/>.
        /// </summary>
        public double LogWeight { get; private set; }

        // ── M-Step scratchpads (pre-allocated; internal to GaussianMixtureModel) ─

        internal double[] ScratchMean { get; }
        internal double[,] ScratchCov { get; }

        // ── Cholesky work arrays ─────────────────────────────────────────────────

        private readonly CholeskyDecomposition _chol;

        public GaussianComponent(int dim)
        {
            _dim = dim;
            Mean = new double[dim];
            Covariance = new double[dim, dim];
            CovarianceInverse = new double[dim, dim];
            ScratchMean = new double[dim];
            ScratchCov = new double[dim, dim];
            _chol = new CholeskyDecomposition(dim);
        }

        // ── M-Step helpers ───────────────────────────────────────────────────────

        /// <summary>Zeroes scratchpads. Call once per component per EM iteration.</summary>
        internal void ResetScratch()
        {
            Array.Clear(ScratchMean, 0, ScratchMean.Length);
            Array.Clear(ScratchCov, 0, ScratchCov.Length);
        }

        /// <summary>
        /// Commits the M-step update: copies normalised <see cref="ScratchMean"/> into
        /// <see cref="Mean"/>, normalises <see cref="ScratchCov"/> into
        /// <see cref="Covariance"/> with diagonal ridge, updates <see cref="Weight"/>,
        /// and refreshes the cache. ScratchMean must already be divided by
        /// <paramref name="effectiveCount"/>.
        /// </summary>
        internal void CommitScratch(double effectiveCount, int n, double regularization = 1e-6)
        {
            Weight = effectiveCount / n;

            if (effectiveCount > 1e-12)
            {
                Array.Copy(ScratchMean, Mean, _dim);

                double invCount = 1.0 / effectiveCount;
                for (int d1 = 0; d1 < _dim; d1++)
                    for (int d2 = 0; d2 < _dim; d2++)
                    {
                        Covariance[d1, d2] = ScratchCov[d1, d2] * invCount
                                             + (d1 == d2 ? regularization : 0.0);
                    }
            }
            // Dead cluster (effectiveCount ≈ 0): leave Mean/Covariance untouched and
            // let Weight ≈ 0 suppress this component in subsequent E-steps.

            UpdateCache();
        }

        // ── Cache update ─────────────────────────────────────────────────────────

        /// <summary>
        /// Refreshes Σ⁻¹, ln|Σ|, ln π_k, and the Cholesky factor from the current
        /// <see cref="Weight"/> and <see cref="Covariance"/>. Call after any direct
        /// assignment to those fields.
        /// </summary>
        public void UpdateCache()
        {
            _chol.Decompose(Covariance);
            _chol.WriteInverseTo(CovarianceInverse);
            LogNormalizationFactor = -0.5 * (_dim * Math.Log(2.0 * Math.PI) + _chol.LogDet);
            LogWeight = Math.Log(Math.Max(Weight, 1e-300));
        }

        // ── PDF evaluation ───────────────────────────────────────────────────────

        /// <summary>ln N(x | μ, Σ). Numerically safe in all dimensions.</summary>
        public double EvaluateLogPdf(double[] x) =>
            LogNormalizationFactor - 0.5 * ComputeMahalanobisSquared(x);

        /// <summary>N(x | μ, Σ). May underflow in high dimensions; use log form when possible.</summary>
        public double EvaluatePdf(double[] x) =>
            Math.Exp(EvaluateLogPdf(x));

        /// <summary>Squared Mahalanobis distance from <paramref name="x"/> to <see cref="Mean"/>.</summary>
        public double MahalanobisSquared(double[] x)
            => Mahalanobis.DistanceSquared(x, Mean, CovarianceInverse);

        /// <summary>Draws one sample using the cached Cholesky factor.</summary>
        public double[] Sample(Random rng) => _chol.Sample(rng, Mean);

        private double ComputeMahalanobisSquared(double[] x)
            => Mahalanobis.DistanceSquared(x, Mean, CovarianceInverse);
    }
}
