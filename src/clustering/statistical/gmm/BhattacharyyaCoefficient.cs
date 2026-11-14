using System;
using Maths.LinAlg;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Bhattacharyya distance and coefficient between two Gaussians.
    /// D_B = ⅛ (μ₁−μ₂)ᵀ Σ̄⁻¹ (μ₁−μ₂) + ½ ln(|Σ̄| / √(|Σ₁||Σ₂|)) with Σ̄ = (Σ₁+Σ₂)/2;
    /// BC = exp(−D_B) ∈ [0,1] (1 = identical, 0 = no overlap). Useful as a pairwise
    /// overlap test between fitted components.
    /// </summary>
    public static class BhattacharyyaCoefficient
    {
        private static readonly double HalfLn2Pi = 0.5 * Math.Log(2.0 * Math.PI);

        // ── Public API ────────────────────────────────────────────────────────────

        /// <summary>Coefficient BC = exp(−D_B) ∈ [0, 1].</summary>
        public static double Between(GaussianComponent a, GaussianComponent b)
            => Math.Exp(-Distance(a, b));

        /// <summary>Distance D_B ∈ [0, ∞). Zero for identical distributions.</summary>
        public static double Distance(GaussianComponent a, GaussianComponent b)
        {
            int d = a.Mean.Length;

            double[,] sigmaBar = new double[d, d];
            double[,] sigmaBarInv = new double[d, d];

            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    sigmaBar[i, j] = 0.5 * (a.Covariance[i, j] + b.Covariance[i, j]);

            if (!double.IsFinite(a.LogNormalizationFactor)) throw new ArgumentException(NoFactor("Component a's covariance"), nameof(a));
            if (!double.IsFinite(b.LogNormalizationFactor)) throw new ArgumentException(NoFactor("Component b's covariance"), nameof(b));

            var chol = new CholeskyDecomposition(d);
            Factor(chol, sigmaBar, "The averaged covariance (Σ₁ + Σ₂)/2", null);
            chol.WriteInverseTo(sigmaBarInv);
            double logDetSigmaBar = chol.LogDet;

            // ln|Σ| recovered from cached LogNormalizationFactor = −½(D·ln(2π) + ln|Σ|).
            double logDetA = -(2.0 * a.LogNormalizationFactor + d * 2.0 * HalfLn2Pi);
            double logDetB = -(2.0 * b.LogNormalizationFactor + d * 2.0 * HalfLn2Pi);

            double mahalTerm = 0.0;
            for (int i = 0; i < d; i++)
            {
                double row = 0.0;
                for (int j = 0; j < d; j++)
                    row += sigmaBarInv[i, j] * (a.Mean[j] - b.Mean[j]);
                mahalTerm += row * (a.Mean[i] - b.Mean[i]);
            }
            mahalTerm *= 0.125;

            double logDetTerm = 0.5 * (logDetSigmaBar - 0.5 * (logDetA + logDetB));
            return mahalTerm + logDetTerm;
        }

        /// <summary>Coefficient overload for callers holding raw (mean, covariance) arrays.</summary>
        public static double Between(
            double[] mean1, double[,] cov1,
            double[] mean2, double[,] cov2)
            => Math.Exp(-Distance(mean1, cov1, mean2, cov2));

        /// <summary>Distance overload for callers holding raw (mean, covariance) arrays.</summary>
        public static double Distance(
            double[] mean1, double[,] cov1,
            double[] mean2, double[,] cov2)
        {
            int d = mean1.Length;

            double[,] sigmaBar = new double[d, d];
            double[,] sigmaBarInv = new double[d, d];

            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    sigmaBar[i, j] = 0.5 * (cov1[i, j] + cov2[i, j]);

            var chol = new CholeskyDecomposition(d);

            Factor(chol, sigmaBar, "The averaged covariance (Σ₁ + Σ₂)/2", null);
            chol.WriteInverseTo(sigmaBarInv);
            double logDetSigmaBar = chol.LogDet;

            Factor(chol, cov1, "The first covariance", nameof(cov1));
            double logDetCov1 = chol.LogDet;

            Factor(chol, cov2, "The second covariance", nameof(cov2));
            double logDetCov2 = chol.LogDet;

            double mahalTerm = 0.0;
            for (int i = 0; i < d; i++)
            {
                double row = 0.0;
                for (int j = 0; j < d; j++)
                    row += sigmaBarInv[i, j] * (mean1[j] - mean2[j]);
                mahalTerm += row * (mean1[i] - mean2[i]);
            }
            mahalTerm *= 0.125;

            double logDetTerm = 0.5 * (logDetSigmaBar - 0.5 * (logDetCov1 + logDetCov2));

            return mahalTerm + logDetTerm;
        }

        // Factors Σ, refusing one with no factor (CholeskyDecomposition.FailedAt): the distance has no value there.
        private static void Factor(CholeskyDecomposition chol, double[,] sigma, string what, string? parameter)
        {
            chol.Decompose(sigma);
            if (chol.FailedAt >= 0)
                throw new ArgumentException(
                    $"{NoFactor(what)} It is not positive definite to working precision at row {chol.FailedAt}.", parameter);
        }

        private static string NoFactor(string what) => $"{what} has no Cholesky factor, so the distance has no value.";
    }
}
