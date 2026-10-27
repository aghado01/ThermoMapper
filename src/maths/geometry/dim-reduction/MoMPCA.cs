using System;
using Maths.LinAlg;
using Maths.Geometry;
using Maths.Geometry.Estimators.Intrinsic;

namespace Maths.Geometry.DimReduction
{
    public static class MoMPCA
    {
        /// <summary>
        /// Robust distributed PCA by scale-calibrated median-of-means aggregation, after Kisung You,
        /// "Scale-Calibrated Median-of-Means for Robust Distributed PCA" (arXiv 2605.20681). Splits the
        /// data into K blocks, runs local PCA per block, and aggregates the block (mean, subspace)
        /// estimates with a single geometric (spatial) median on the product manifold Rᵖ × Gr(r, p),
        /// via the shared IRLS solver (<see cref="GeometricMedian"/>).
        ///
        /// <para>The product metric is scale-calibrated (paper §5.3): d²_α = α·d_Eucl² + (2−α)·d_Grass²
        /// with α̂ = 2·τ̂_U/(τ̂_μ + τ̂_U) from the per-tangent-dimension robust block scales — so a noisy
        /// subspace (weak eigengap) is automatically down-weighted relative to the mean.</para>
        ///
        /// <para>The paper estimates (μ, U) and no more. Within Û this implementation adds principal axes and
        /// robust variances: the geometric median, in the Frobenius norm, of the blocks' covariances projected into Û,
        /// and its eigendecomposition. The median resists a minority of blocks with inflated variance as the subspace
        /// median resists tilted ones; where its eigenvalues tie, the axes between them are arbitrary, as in any
        /// PCA.</para>
        /// </summary>
        /// <param name="data">The full dataset</param>
        /// <param name="kBlocks">Number of blocks to split the data into</param>
        /// <param name="numComponents">Dimension r of the principal subspace</param>
        public static MoMPcaResult ComputeMoM(double[][] data, int kBlocks, int numComponents)
        {
            if (kBlocks <= 1)
            {
                PcaResult pca = Pca.Compute(data, numComponents, center: true, whiten: false);
                return new MoMPcaResult(pca.Mean, pca.Components, pca.Eigenvalues, double.NaN);
            }
            int n = data.Length;
            int d = data[0].Length;
            int r = numComponents;
            int blockSize = n / kBlocks;

            // 1. Local PCA per block → block mean (Euclidean point) + block subspace (Grassmann point,
            //    the r orthonormal components packed column-major as a d×r frame).
            double[][] blockMeans = new double[kBlocks][];
            double[][] blockFrames = new double[kBlocks][];
            for (int k = 0; k < kBlocks; k++)
            {
                int start = k * blockSize;
                int end = (k == kBlocks - 1) ? n : start + blockSize;
                int localN = end - start;

                double[][] localData = new double[localN][];
                Array.Copy(data, start, localData, 0, localN);

                var localPca = Pca.Compute(localData, numComponents, center: true, whiten: false);
                blockMeans[k] = localPca.Mean;
                blockFrames[k] = PackFrame(localPca.Components, d, r);
            }

            double[] weights = new double[kBlocks];
            Array.Fill(weights, 1.0);

            var euclid = new EuclideanVectorManifold(d);
            var grass = new GrassmannManifold(d, r);

            // 2. Preliminary per-factor medians — the calibration center and the joint warm-start.
            double[] meanPrelim = ArithmeticMean(blockMeans, d);
            GeometricMedian.Compute(euclid, blockMeans, weights, meanPrelim);

            // Warm-start at the medoid frame rather than blockFrames[0]: a contaminated leading
            // block can sit on the Grassmann cut locus of the clean majority, where Weiszfeld
            // stalls on the corrupted initialization (see GeometricMedian.MedoidIndex) — and a
            // stalled framePrelim poisons the scale calibration and the joint warm-start below.
            double[] framePrelim = (double[])blockFrames[GeometricMedian.MedoidIndex(grass, blockFrames)].Clone();
            GeometricMedian.Compute(grass, blockFrames, weights, framePrelim);

            // 3. Scale calibration α̂ = 2·τ̂_U/(τ̂_μ + τ̂_U) from the per-tangent-dimension block scales.
            double alpha = CalibrateScale(blockMeans, meanPrelim, blockFrames, framePrelim, grass, d, r);

            // 4. Calibrated joint median on the scaled product Rᵖ × Gr(r, p). Metric scaling leaves
            //    Log/Exp unchanged (You 2601.10992); ScaledManifold only reweights the factor distances.
            var product = new ProductManifold<ScaledManifold<EuclideanVectorManifold>, ScaledManifold<GrassmannManifold>>(
                new ScaledManifold<EuclideanVectorManifold>(euclid, alpha),
                new ScaledManifold<GrassmannManifold>(grass, 2.0 - alpha));

            int frameDim = d * r;
            int prodDim = d + frameDim;
            double[][] productData = new double[kBlocks][];
            for (int k = 0; k < kBlocks; k++)
            {
                var point = new double[prodDim];
                Array.Copy(blockMeans[k], 0, point, 0, d);
                Array.Copy(blockFrames[k], 0, point, d, frameDim);
                productData[k] = point;
            }
            double[] productMedian = new double[prodDim];
            Array.Copy(meanPrelim, 0, productMedian, 0, d);
            Array.Copy(framePrelim, 0, productMedian, d, frameDim);
            GeometricMedian.Compute(product, productData, weights, productMedian);

            // 5. Unpack the robust mean and subspace.
            double[] robustMean = productMedian[..d];
            double[][] basis = new double[r][];
            for (int c = 0; c < r; c++)
            {
                basis[c] = new double[d];
                Array.Copy(productMedian, d + c * d, basis[c], 0, d);
            }

            // 6. Principal axes within Û (beyond the paper): the median of the blocks' projected covariances,
            //    then its eigenvectors, ordered by eigenvalue and signed as Pca signs its components.
            var projected = new double[kBlocks][];
            for (int k = 0; k < kBlocks; k++)
            {
                int start = k * blockSize;
                int end = (k == kBlocks - 1) ? n : start + blockSize;
                projected[k] = ProjectedCovariance(data, start, end, blockMeans[k], basis, r);
            }
            double[] covarianceMedian = ArithmeticMean(projected, r * r);
            GeometricMedian.Compute(new EuclideanVectorManifold(r * r), projected, weights, covarianceMedian);
            var median = new double[r, r];
            for (int a = 0; a < r; a++)
                for (int b = 0; b < r; b++)
                    median[a, b] = 0.5 * (covarianceMedian[a * r + b] + covarianceMedian[b * r + a]);
            EigenResult eig = DenseEigen.DecomposeSymmetric(median);

            var order = new int[r];
            for (int i = 0; i < r; i++) order[i] = i;
            Array.Sort(order, (i, j) => eig.Eigenvalues[j].CompareTo(eig.Eigenvalues[i]));
            var axes = new double[r][];
            var variances = new double[r];
            for (int i = 0; i < r; i++)
            {
                double[] q = eig.Eigenvectors[order[i]];
                axes[i] = new double[d];
                for (int c = 0; c < r; c++)
                    for (int j = 0; j < d; j++) axes[i][j] += q[c] * basis[c][j];
                Pca.Orient(axes[i]);
                variances[i] = Math.Max(eig.Eigenvalues[order[i]], 0.0);
            }
            return new MoMPcaResult(robustMean, axes, variances, alpha);
        }

        // Block k's covariance about its own mean, projected into the subspace: (1/n_k)·Σ_i z_i z_iᵀ with
        // z_i = Uᵀ(x_i − μ_k), the n_k denominator Pca uses; flattened row-major, r × r.
        private static double[] ProjectedCovariance(double[][] data, int start, int end, double[] mean, double[][] basis, int r)
        {
            int d = mean.Length;
            var s = new double[r * r];
            var z = new double[r];
            for (int i = start; i < end; i++)
            {
                for (int c = 0; c < r; c++)
                {
                    double dot = 0.0;
                    for (int j = 0; j < d; j++) dot += (data[i][j] - mean[j]) * basis[c][j];
                    z[c] = dot;
                }
                for (int a = 0; a < r; a++)
                    for (int b = 0; b < r; b++) s[a * r + b] += z[a] * z[b];
            }
            int count = end - start;
            for (int t = 0; t < s.Length; t++) s[t] /= count;
            return s;
        }

        // PCA components (r vectors of length d) → column-major d×r Grassmann frame [v0 | v1 | …].
        private static double[] PackFrame(double[][] components, int d, int r)
        {
            var frame = new double[d * r];
            for (int c = 0; c < r; c++)
                Array.Copy(components[c], 0, frame, c * d, d);
            return frame;
        }

        private static double[] ArithmeticMean(double[][] points, int d)
        {
            int k = points.Length;
            var mean = new double[d];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < d; j++) mean[j] += points[i][j] / k;
            return mean;
        }

        // α̂ = 2·τ̂_U/(τ̂_μ + τ̂_U), with τ̂_μ = median_k‖μ_k − μ̃‖² / p and
        // τ̂_U = median_k d_Gr(U_k, Ũ)² / [r(p−r)] (the √b block-size standardizer cancels in the ratio).
        private static double CalibrateScale(
            double[][] blockMeans, double[] meanCenter,
            double[][] blockFrames, double[] frameCenter,
            GrassmannManifold grass, int p, int r)
        {
            int gIntrinsic = r * (p - r);
            if (gIntrinsic <= 0) return 1.0;   // full subspace: no Grassmann variation, stay balanced

            int k = blockMeans.Length;
            var dMu = new double[k];
            var dU = new double[k];
            for (int i = 0; i < k; i++)
            {
                double s = 0;
                for (int j = 0; j < p; j++) { double e = blockMeans[i][j] - meanCenter[j]; s += e * e; }
                dMu[i] = Math.Sqrt(s);
                dU[i] = grass.Distance(blockFrames[i], frameCenter);
            }

            double mMu = Median(dMu);
            double mU = Median(dU);
            double tauMu = mMu * mMu / p;
            double tauU = mU * mU / gIntrinsic;
            double denom = tauMu + tauU;
            if (denom <= 0) return 1.0;
            return Math.Clamp(2.0 * tauU / denom, 1e-6, 2.0 - 1e-6);
        }

        private static double Median(double[] values)
        {
            var sorted = (double[])values.Clone();
            Array.Sort(sorted);
            int m = sorted.Length;
            return m % 2 == 1 ? sorted[m / 2] : 0.5 * (sorted[m / 2 - 1] + sorted[m / 2]);
        }
    }

    /// <summary>
    /// The scale-calibrated median-of-means estimate of θ = (μ, U) ∈ ℝᵖ × Gr(r, p), You (arXiv 2605.20681, §2.1),
    /// with principal axes and robust variances within Û. The paper's target is the subspace itself, "not a
    /// particular signed, ordered, or oriented orthonormal basis"; the axes and variances are this implementation's
    /// extension, from the median of the blocks' covariances projected into Û.
    /// </summary>
    /// <param name="Mean">μ̂, the mean factor of the product median.</param>
    /// <param name="Basis">r orthonormal principal axes spanning Û, the subspace factor, in descending order of
    /// <paramref name="Variances"/>, each signed so its largest-magnitude entry is positive (as <see cref="Pca"/>
    /// signs its components).</param>
    /// <param name="Variances">The eigenvalues of the median projected block covariance, descending: each axis's
    /// variance, with the n denominator, robust to a minority of blocks.</param>
    /// <param name="Alpha">The calibrated scale α̂ ∈ (0, 2) of the product metric (§5.3); NaN for a single block,
    /// where there is nothing to aggregate and the estimate is the full-data PCA.</param>
    public sealed record MoMPcaResult(double[] Mean, double[][] Basis, double[] Variances, double Alpha);
}
