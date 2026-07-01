using System;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Pre-allocated scratch for repeated <c>ModeAscent.Ascend</c> calls.
    /// </summary>
    public sealed class ModeAscentScratch
    {
        internal readonly double[] X;
        internal readonly double[] XNew;
        internal readonly double[] LogR;
        internal readonly double[] R;
        internal readonly double[] Grad;

        public ModeAscentScratch(int dimension, int components)
        {
            X = new double[dimension];
            XNew = new double[dimension];
            LogR = new double[components];
            R = new double[components];
            Grad = new double[dimension];
        }
    }

    /// <summary>
    /// Gradient ascent on a Gaussian mixture log-density. Used by
    /// <see cref="ModalMergeStrategy"/> to map components to density basins.
    /// Analytic gradient ∇ log p(x) = Σ_k r_k(x) · Σ_k⁻¹ (μ_k − x); backtracking
    /// line search keeps the density monotone.
    /// </summary>
    public static class ModeAscent
    {
        // ── Public API ────────────────────────────────────────────────────────────

        /// <summary>
        /// Ascends from <paramref name="start"/> to the nearest local mode of the
        /// mixture log-density. Components must have up-to-date caches.
        /// </summary>
        public static double[] Ascend(
            double[] start,
            GaussianComponent[] components,
            int maxSteps = 50,
            double tol = 1e-7)
        {
            var scratch = new ModeAscentScratch(start.Length, components.Length);
            double[] result = new double[start.Length];
            Ascend(start, components, result, scratch, maxSteps, tol);
            return result;
        }

        /// <summary>
        /// Allocation-free overload: caller supplies a result buffer and reusable
        /// scratch. The result buffer receives the converged mode.
        /// </summary>
        public static void Ascend(
            double[] start,
            GaussianComponent[] components,
            double[] result,
            ModeAscentScratch scratch,
            int maxSteps = 50,
            double tol = 1e-7)
        {
            int d = start.Length;
            int k = components.Length;

            double[] x = scratch.X;
            double[] xNew = scratch.XNew;
            double[] logR = scratch.LogR;
            double[] r = scratch.R;
            double[] grad = scratch.Grad;

            Array.Copy(start, x, d);

            for (int step = 0; step < maxSteps; step++)
            {
                // Log-space responsibilities r_k(x) for stability at large Mahalanobis.
                double logMax = double.NegativeInfinity;
                for (int ki = 0; ki < k; ki++)
                {
                    logR[ki] = components[ki].LogWeight
                               + components[ki].EvaluateLogPdf(x);
                    if (logR[ki] > logMax) logMax = logR[ki];
                }

                double rSum = 0.0;
                for (int ki = 0; ki < k; ki++)
                {
                    r[ki] = Math.Exp(logR[ki] - logMax);
                    rSum += r[ki];
                }
                for (int ki = 0; ki < k; ki++) r[ki] /= rSum;

                // ∇ log p(x) = Σ_k r_k · Σ_k⁻¹ (μ_k − x).
                Array.Clear(grad, 0, d);
                for (int ki = 0; ki < k; ki++)
                {
                    if (r[ki] < 1e-15) continue;
                    double[] mu = components[ki].Mean;
                    double[,] sigmaInv = components[ki].CovarianceInverse;
                    for (int i = 0; i < d; i++)
                    {
                        double dot = 0.0;
                        for (int j = 0; j < d; j++)
                            dot += sigmaInv[i, j] * (mu[j] - x[j]);
                        grad[i] += r[ki] * dot;
                    }
                }

                // Backtracking line search keeps density monotone (up to 30 halvings).
                double logP0 = LogMixtureDensity(x, components, logR);
                double alpha = 1.0;
                bool improved = false;

                for (int ls = 0; ls < 30; ls++)
                {
                    for (int i = 0; i < d; i++) xNew[i] = x[i] + alpha * grad[i];
                    if (LogMixtureDensity(xNew, components, logR) > logP0)
                    {
                        improved = true;
                        break;
                    }
                    alpha *= 0.5;
                }

                if (!improved) break;  // already at a local mode

                double disp = 0.0;
                for (int i = 0; i < d; i++)
                {
                    double delta = xNew[i] - x[i];
                    disp += delta * delta;
                    x[i] = xNew[i];
                }
                if (Math.Sqrt(disp) < tol) break;
            }

            Array.Copy(x, result, d);
        }

        /// <summary>
        /// Index of the component whose density basin contains <paramref name="point"/>:
        /// ascends to the local mode, then returns the component with the highest
        /// weighted log-pdf at that mode.
        /// </summary>
        public static int GetBasin(
            double[] point,
            GaussianComponent[] components,
            int maxSteps = 50,
            double tol = 1e-7)
        {
            double[] mode = Ascend(point, components, maxSteps, tol);

            int best = 0;
            double bestScore = double.NegativeInfinity;
            for (int ki = 0; ki < components.Length; ki++)
            {
                double score = components[ki].LogWeight
                               + components[ki].EvaluateLogPdf(mode);
                if (score > bestScore) { bestScore = score; best = ki; }
            }
            return best;
        }

        // ── Internals ─────────────────────────────────────────────────────────────

        // Log mixture density at x via log-sum-exp. scratch length ≥ components.Length.
        private static double LogMixtureDensity(
            double[] x,
            GaussianComponent[] components,
            double[] scratch)
        {
            double logMax = double.NegativeInfinity;
            for (int ki = 0; ki < components.Length; ki++)
            {
                scratch[ki] = components[ki].LogWeight
                              + components[ki].EvaluateLogPdf(x);
                if (scratch[ki] > logMax) logMax = scratch[ki];
            }
            double sum = 0.0;
            for (int ki = 0; ki < components.Length; ki++)
                sum += Math.Exp(scratch[ki] - logMax);
            return logMax + Math.Log(sum);
        }
    }
}
