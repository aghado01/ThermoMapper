using System;
using System.Buffers;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Geometry.Estimators.Ambient;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Flat K-component finite Gaussian mixture, fit by EM. See docs/gmm.md for
    /// the user-facing reference (mental model, fit modes, supervision, K-selection,
    /// merge strategies).
    /// </summary>
    public sealed class GaussianMixtureModel
    {
        public int K { get; }
        public int Dimension { get; }
        public GaussianComponent[] Components { get; }

        /// <summary>
        /// Diagonal ridge added to each component's covariance after the M-step and
        /// to the variance floor in <see cref="RandomInitialize"/>. Default 1e-6;
        /// adjust to match feature-scale magnitude when data is not roughly O(1).
        /// </summary>
        public double CovarianceRegularization { get; }

        /// <summary>True if the last <c>Fit</c> call converged within the tolerance.</summary>
        public bool IsConverged { get; private set; }

        /// <summary>Log-likelihood at the end of the last <c>Fit</c> call.</summary>
        public double FinalLogLikelihood { get; private set; }

        /// <summary>Number of EM iterations performed in the last <c>Fit</c> call.</summary>
        public int NumIterations { get; private set; }

        // Persisted from the most recent Fit call. Exposed via GetFinalResponsibilities().
        private double[,]? _responsibilities;

        // True once InitializeWithParameters has been called or RandomInitialize has run.
        private bool _isInitialized;

        public GaussianMixtureModel(int k, int dimension, double covarianceRegularization = 1e-6)
        {
            if (covarianceRegularization < 0.0)
                throw new ArgumentOutOfRangeException(nameof(covarianceRegularization), "Must be ≥ 0.");

            K = k;
            Dimension = dimension;
            CovarianceRegularization = covarianceRegularization;
            Components = new GaussianComponent[K];
            for (int i = 0; i < K; i++)
                Components[i] = new GaussianComponent(dimension);
        }

        // ── Initialisation ───────────────────────────────────────────────────────

        /// <summary>
        /// Warm-start with explicit means, covariances, and weights. Call before
        /// <c>Fit</c> to bypass random-sample initialisation.
        /// </summary>
        /// <param name="means">K vectors of length <see cref="Dimension"/>.</param>
        /// <param name="covariances">K matrices of shape <see cref="Dimension"/>².</param>
        /// <param name="weights">K weights; normalised internally.</param>
        public void InitializeWithParameters(double[][] means, double[][,] covariances, double[] weights)
        {
            double weightSum = 0.0;
            for (int i = 0; i < K; i++) weightSum += weights[i];
            double invWeightSum = weightSum > 0.0 ? 1.0 / weightSum : 1.0;

            for (int i = 0; i < K; i++)
            {
                Array.Copy(means[i], Components[i].Mean, Dimension);
                Array.Copy(covariances[i], Components[i].Covariance, Dimension * Dimension);
                Components[i].Weight = weights[i] * invWeightSum;
                Components[i].UpdateCache();
            }
            _isInitialized = true;
        }

        // ── Fitting ──────────────────────────────────────────────────────────────

        /// <summary>Fits with fully unsupervised EM.</summary>
        public void Fit(double[][] data, int maxIterations = 100, double tolerance = 1e-4)
            => FitCore(data, hardLabels: null, constraint: null, maxIterations, tolerance);

        /// <summary>Fits with a soft-label responsibility constraint.</summary>
        public void Fit(double[][] data, IResponsibilityConstraint constraint,
                        int maxIterations = 100, double tolerance = 1e-4)
            => FitCore(data, hardLabels: null, constraint, maxIterations, tolerance);

        /// <summary>
        /// Fits with hard labels on a subset of points. <paramref name="hardLabels"/>[i]
        /// is either the pinned component index (≥ 0) or −1 (unconstrained).
        /// </summary>
        public void Fit(double[][] data, int[] hardLabels, int maxIterations = 100, double tolerance = 1e-4)
            => FitCore(data, hardLabels, constraint: null, maxIterations, tolerance);

        /// <summary>
        /// Fits with hard labels on a pinned subset and a soft constraint on the rest.
        /// Pinned rows are not perturbed by the constraint.
        /// </summary>
        public void Fit(double[][] data, int[] hardLabels, IResponsibilityConstraint constraint,
                        int maxIterations = 100, double tolerance = 1e-4)
            => FitCore(data, hardLabels, constraint, maxIterations, tolerance);

        private void FitCore(double[][] data, int[]? hardLabels, IResponsibilityConstraint? constraint,
                             int maxIterations, double tolerance)
        {
            int n = data.Length;

            if (!_isInitialized)
                RandomInitialize(data);

            // Allocate / resize responsibility matrix only when data size changes.
            if (_responsibilities == null || _responsibilities.GetLength(0) != n)
                _responsibilities = new double[n, K];

            // Per-point log-likelihood accumulator; allocated once, reused across iterations.
            double[] perPointLL = new double[n];

            double prevLogLikelihood = double.MinValue;
            IsConverged = false;
            NumIterations = 0;

            // Constraint is applied *after* convergence is decided so that on return
            // _responsibilities reflects what the model produced — not the blended
            // signal that drives the M-step.
            for (int iter = 0; iter < maxIterations; iter++)
            {
                double currentLogLikelihood = ExpectationStep(data, hardLabels, _responsibilities, perPointLL);
                NumIterations = iter + 1;

                if (Math.Abs(currentLogLikelihood - prevLogLikelihood) < tolerance)
                {
                    IsConverged = true;
                    FinalLogLikelihood = currentLogLikelihood;
                    return;
                }

                prevLogLikelihood = currentLogLikelihood;
                constraint?.Apply(_responsibilities, n, K, iter, maxIterations);
                if (hardLabels != null && constraint != null)
                    RepinHardLabels(_responsibilities, hardLabels, n);
                MaximizationStep(data, _responsibilities, n);
            }

            // Loop exhausted: re-run E-step so _responsibilities matches final params.
            FinalLogLikelihood = ExpectationStep(data, hardLabels, _responsibilities, perPointLL);
        }

        // ── Random initialisation (randSample) ───────────────────────────────────

        /// <summary>
        /// Random-sample initialisation: K means drawn from data without replacement,
        /// per-dimension sample variance as diagonal covariance, uniform weights. Pass
        /// an explicit <paramref name="rng"/> for reproducibility.
        /// </summary>
        public void RandomInitialize(double[][] data, Random? rng = null)
        {
            rng ??= new Random();
            int n = data.Length;

            // Per-dimension mean and variance over the full dataset.
            double[] mean = new double[Dimension];
            for (int i = 0; i < n; i++)
                for (int d = 0; d < Dimension; d++)
                    mean[d] += data[i][d];
            for (int d = 0; d < Dimension; d++) mean[d] /= n;

            double[] variance = new double[Dimension];
            for (int i = 0; i < n; i++)
                for (int d = 0; d < Dimension; d++)
                {
                    double diff = data[i][d] - mean[d];
                    variance[d] += diff * diff;
                }
            // Bessel-corrected; floor at the configured regularisation value to prevent zero.
            double invNm1 = n > 1 ? 1.0 / (n - 1) : 1.0;
            for (int d = 0; d < Dimension; d++)
                variance[d] = variance[d] * invNm1 + CovarianceRegularization;

            // Partial Fisher-Yates to sample K distinct row indices.
            int[] indices = new int[n];
            for (int i = 0; i < n; i++) indices[i] = i;
            for (int i = 0; i < K; i++)
            {
                int j = i + rng.Next(n - i);
                (indices[i], indices[j]) = (indices[j], indices[i]);
            }

            double uniformWeight = 1.0 / K;
            for (int k = 0; k < K; k++)
            {
                Array.Copy(data[indices[k]], Components[k].Mean, Dimension);
                Components[k].Weight = uniformWeight;
                Array.Clear(Components[k].Covariance, 0, Dimension * Dimension);
                for (int d = 0; d < Dimension; d++)
                    Components[k].Covariance[d, d] = variance[d];
                Components[k].UpdateCache();
            }
            _isInitialized = true;
        }

        // ── Robust initialisation (geometric median + Weiszfeld scatter) ──────────

        /// <summary>
        /// Robust initialization: seeds each component's μ and Σ via
        /// <see cref="EuclideanMedian.ComputeWithScatter"/> over points in the
        /// nearest-seed cluster. Strictly better than <see cref="RandomInitialize"/>
        /// for contaminated or heavy-tailed data; EM then proceeds standard with
        /// monotonic LL improvement and valid BIC intact.
        /// </summary>
        /// <param name="data">Full dataset.</param>
        /// <param name="seeds">
        /// K distinct row indices into <paramref name="data"/> used as warm-start
        /// centroids (one per component). Pass <c>null</c> to draw K seeds via
        /// partial Fisher-Yates (same strategy as <see cref="RandomInitialize"/>).
        /// </param>
        /// <param name="consistencyFactor">
        /// Bias-correction multiplier forwarded to
        /// <see cref="EuclideanMedian.ComputeWithScatter"/>. Default 1.0; see
        /// <c>ConsistencyFactors.cs</c> for principled Gaussian-consistent choices.
        /// </param>
        /// <param name="rng">Random source used when <paramref name="seeds"/> is null.</param>
        public void RobustInitialize(
            double[][] data,
            int[]? seeds = null,
            double consistencyFactor = 1.0,
            Random? rng = null)
        {
            int n = data.Length;
            rng ??= new Random();

            // ── Step 1: K seed indices ────────────────────────────────────────────
            int[] seedIndices;
            if (seeds != null)
            {
                if (seeds.Length != K)
                    throw new ArgumentException(
                        $"seeds.Length ({seeds.Length}) must equal K ({K}).", nameof(seeds));
                seedIndices = seeds;
            }
            else
            {
                int[] pool = new int[n];
                for (int i = 0; i < n; i++) pool[i] = i;
                for (int i = 0; i < K; i++)
                {
                    int j = i + rng.Next(n - i);
                    (pool[i], pool[j]) = (pool[j], pool[i]);
                }
                seedIndices = new int[K];
                Array.Copy(pool, seedIndices, K);
            }

            // ── Step 2: nearest-seed hard assignment ──────────────────────────────
            List<int>[] clusters = new List<int>[K];
            for (int k = 0; k < K; k++) clusters[k] = new List<int>();

            for (int i = 0; i < n; i++)
            {
                int nearest = 0;
                double nearestDist = double.MaxValue;
                for (int k = 0; k < K; k++)
                {
                    double dist = 0.0;
                    double[] seed = data[seedIndices[k]];
                    for (int d = 0; d < Dimension; d++)
                    {
                        double diff = data[i][d] - seed[d];
                        dist += diff * diff;
                    }
                    if (dist < nearestDist) { nearestDist = dist; nearest = k; }
                }
                clusters[nearest].Add(i);
            }

            // ── Step 3: robust mean + scatter per component ───────────────────────
            double[] flatScatter = ArrayPool<double>.Shared.Rent(Dimension * Dimension);
            double uniformWeight = 1.0 / K;

            try
            {
                for (int k = 0; k < K; k++)
                {
                    GaussianComponent comp = Components[k];
                    List<int> cluster = clusters[k];

                    if (cluster.Count == 0)
                    {
                        // Empty cluster: fall back to seed point, identity covariance.
                        Array.Copy(data[seedIndices[k]], comp.Mean, Dimension);
                        Array.Clear(comp.Covariance, 0, Dimension * Dimension);
                        for (int d = 0; d < Dimension; d++)
                            comp.Covariance[d, d] = CovarianceRegularization;
                        comp.Weight = uniformWeight;
                        comp.UpdateCache();
                        continue;
                    }

                    int[] idxArr = cluster.ToArray();
                    double[] wArr = new double[idxArr.Length];
                    Array.Fill(wArr, 1.0); // uniform — geometric median of cluster members

                    // Warm-start destination with cluster Euclidean mean.
                    Array.Clear(comp.Mean, 0, Dimension);
                    for (int ii = 0; ii < idxArr.Length; ii++)
                        for (int d = 0; d < Dimension; d++)
                            comp.Mean[d] += data[idxArr[ii]][d];
                    double invCount = 1.0 / idxArr.Length;
                    for (int d = 0; d < Dimension; d++) comp.Mean[d] *= invCount;

                    Span<double> scatter = flatScatter.AsSpan(0, Dimension * Dimension);
                    scatter.Clear();

                    EuclideanMedian.ComputeWithScatter(
                        data,
                        idxArr.AsSpan(),
                        wArr.AsSpan(),
                        comp.Mean.AsSpan(),
                        scatter,
                        consistencyFactor: consistencyFactor);

                    // Copy flat row-major scatter into comp.Covariance with diagonal ridge.
                    for (int d1 = 0; d1 < Dimension; d1++)
                        for (int d2 = 0; d2 < Dimension; d2++)
                            comp.Covariance[d1, d2] = scatter[d1 * Dimension + d2]
                                                     + (d1 == d2 ? CovarianceRegularization : 0.0);

                    comp.Weight = uniformWeight;
                    comp.UpdateCache();
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(flatScatter);
            }

            _isInitialized = true;
        }

        private void RepinHardLabels(double[,] responsibilities, int[] hardLabels, int n)
        {
            for (int i = 0; i < n; i++)
            {
                int locked = hardLabels[i];
                if (locked < 0) continue;
                for (int k = 0; k < K; k++)
                    responsibilities[i, k] = (k == locked) ? 1.0 : 0.0;
            }
        }

        // ── E-Step ───────────────────────────────────────────────────────────────

        private double ExpectationStep(
            double[][] data,
            int[]? hardLabels,
            double[,] responsibilities,
            double[] perPointLL)
        {
            int n = data.Length;
            int localK = K; // captured by the parallel lambda

            Parallel.For(0, n, i =>
            {
                // Hard-label: pin to one-hot but still contribute the pinned component's
                // log-density so FinalLogLikelihood / BIC reflect the true joint LL.
                if (hardLabels != null && hardLabels[i] >= 0)
                {
                    int lockedComp = hardLabels[i];
                    for (int k = 0; k < localK; k++)
                        responsibilities[i, k] = (k == lockedComp) ? 1.0 : 0.0;
                    perPointLL[i] = Components[lockedComp].LogWeight
                                  + Components[lockedComp].EvaluateLogPdf(data[i]);
                    return;
                }

                // Log-sum-exp keeps the E-step stable in high dimensions. stackalloc
                // is safe per parallel invocation; for K ≫ 100 switch to ArrayPool.
                Span<double> logW = stackalloc double[localK];
                double maxLogW = double.MinValue;

                for (int k = 0; k < localK; k++)
                {
                    logW[k] = Components[k].LogWeight
                               + Components[k].EvaluateLogPdf(data[i]);
                    if (logW[k] > maxLogW) maxLogW = logW[k];
                }

                double sumExp = 0.0;
                for (int k = 0; k < localK; k++)
                {
                    double e = Math.Exp(logW[k] - maxLogW);
                    logW[k] = e;
                    sumExp += e;
                }

                double invSumExp = 1.0 / sumExp;
                for (int k = 0; k < localK; k++)
                    responsibilities[i, k] = logW[k] * invSumExp;

                perPointLL[i] = maxLogW + Math.Log(sumExp);
            });

            // Accumulate outside the parallel block — no lock contention.
            double total = 0.0;
            for (int i = 0; i < n; i++) total += perPointLL[i];
            return total;
        }

        // ── M-Step ───────────────────────────────────────────────────────────────

        private void MaximizationStep(double[][] data, double[,] responsibilities, int n)
        {
            for (int k = 0; k < K; k++)
            {
                GaussianComponent comp = Components[k];
                comp.ResetScratch();

                double effectiveCount = 0.0;

                // Pass 1: accumulate weighted mean (raw sum into ScratchMean).
                for (int i = 0; i < n; i++)
                {
                    double resp = responsibilities[i, k];
                    effectiveCount += resp;
                    for (int d = 0; d < Dimension; d++)
                        comp.ScratchMean[d] += resp * data[i][d];
                }

                // CommitScratch expects ScratchMean already divided by effectiveCount.
                if (effectiveCount > 1e-12)
                {
                    double invCount = 1.0 / effectiveCount;
                    for (int d = 0; d < Dimension; d++)
                        comp.ScratchMean[d] *= invCount;
                }

                // Pass 2: weighted centred outer products. Compute upper triangle and
                // mirror immediately to keep ScratchCov symmetric for CommitScratch.
                for (int i = 0; i < n; i++)
                {
                    double resp = responsibilities[i, k];
                    if (resp < 1e-15) continue; // skip negligible contributions

                    for (int d1 = 0; d1 < Dimension; d1++)
                    {
                        double diff1 = data[i][d1] - comp.ScratchMean[d1];
                        for (int d2 = d1; d2 < Dimension; d2++)
                        {
                            double contrib = resp * diff1 * (data[i][d2] - comp.ScratchMean[d2]);
                            comp.ScratchCov[d1, d2] += contrib;
                            if (d1 != d2) comp.ScratchCov[d2, d1] += contrib;
                        }
                    }
                }

                comp.CommitScratch(effectiveCount, n, CovarianceRegularization);
            }
        }

        // ── Inference ────────────────────────────────────────────────────────────

        /// <summary>Hard component assignment per point (argmax weighted log-pdf).</summary>
        public int[] Predict(double[][] data)
        {
            int n = data.Length;
            int[] predictions = new int[n];

            for (int i = 0; i < n; i++)
            {
                int bestK = 0;
                double maxLogPdf = double.MinValue;

                for (int k = 0; k < K; k++)
                {
                    double logPdf = Components[k].LogWeight
                                    + Components[k].EvaluateLogPdf(data[i]);
                    if (logPdf > maxLogPdf)
                    {
                        maxLogPdf = logPdf;
                        bestK = k;
                    }
                }

                predictions[i] = bestK;
            }

            return predictions;
        }

        /// <summary>Soft posterior P(component k | x_i) per point, shape [N, K].</summary>
        public double[,] PredictProba(double[][] data)
        {
            int n = data.Length;
            double[,] result = new double[n, K];
            double[] perPointLL = new double[n];
            ExpectationStep(data, hardLabels: null, result, perPointLL);
            return result;
        }

        /// <summary>
        /// Responsibility matrix [N, K] from the most recent <c>Fit</c> call, or null
        /// if not yet fit. Returns the live internal buffer — do not modify; copy if it
        /// must outlive the next <c>Fit</c> on this instance.
        /// </summary>
        public double[,]? GetFinalResponsibilities() => _responsibilities;

        // ── Density ──────────────────────────────────────────────────────────────

        /// <summary>
        /// Mixture PDF at each point: Σ_k π_k · N(x_i | μ_k, Σ_k). May underflow to 0
        /// in high dimensions; prefer log-space arithmetic for those regimes.
        /// </summary>
        public double[] Pdf(double[][] data)
        {
            int n = data.Length;
            double[] result = new double[n];
            for (int i = 0; i < n; i++)
                for (int k = 0; k < K; k++)
                    result[i] += Components[k].Weight * Components[k].EvaluatePdf(data[i]);
            return result;
        }

        // ── Distance ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Squared Mahalanobis distance per point to each component mean:
        /// (x_i − μ_k)ᵀ Σ_k⁻¹ (x_i − μ_k).
        /// </summary>
        public double[,] Mahal(double[][] data)
        {
            int n = data.Length;
            double[,] result = new double[n, K];
            for (int i = 0; i < n; i++)
                for (int k = 0; k < K; k++)
                    result[i, k] = Components[k].MahalanobisSquared(data[i]);
            return result;
        }

        // ── Sampling ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Draws <paramref name="n"/> samples. Pass <paramref name="componentIndices"/>
        /// (length n) to receive the component each sample came from.
        /// </summary>
        public double[][] Sample(int n, Random? rng = null, int[]? componentIndices = null)
        {
            rng ??= new Random();
            double[][] result = new double[n][];

            // Weight CDF for O(K) multinomial draw; normalise tail so the fallthrough
            // at K-1 is correct when weights don't sum to exactly 1.
            double[] cdf = new double[K];
            cdf[0] = Components[0].Weight;
            for (int k = 1; k < K; k++)
                cdf[k] = cdf[k - 1] + Components[k].Weight;
            double invTotal = cdf[K - 1] > 0.0 ? 1.0 / cdf[K - 1] : 1.0;
            for (int k = 0; k < K; k++) cdf[k] *= invTotal;

            for (int i = 0; i < n; i++)
            {
                double u = rng.NextDouble();
                int compIdx = K - 1;
                for (int k = 0; k < K - 1; k++)
                    if (u < cdf[k]) { compIdx = k; break; }

                result[i] = Components[compIdx].Sample(rng);
                if (componentIndices != null) componentIndices[i] = compIdx;
            }
            return result;
        }
    }
}
