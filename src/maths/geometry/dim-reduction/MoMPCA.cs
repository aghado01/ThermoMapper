using System;
using System.Linq;
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
        /// </summary>
        /// <param name="data">The full dataset</param>
        /// <param name="kBlocks">Number of blocks to split the data into</param>
        /// <param name="numComponents">Number of principal components to extract</param>
        /// <returns>A PcaResult representing the robust aggregation</returns>
        public static PcaResult ComputeMoM(double[][] data, int kBlocks, int numComponents)
        {
            if (kBlocks <= 1) return Pca.Compute(data, numComponents, center: true, whiten: false);
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

            double[] framePrelim = (double[])blockFrames[0].Clone();
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
            double[][] robustSubspace = new double[r][];
            for (int c = 0; c < r; c++)
            {
                robustSubspace[c] = new double[d];
                Array.Copy(productMedian, d + c * d, robustSubspace[c], 0, d);
            }

            // 6. Eigenvalues = variance of the centered data captured along each robust component.
            double[] eigenvalues = new double[r];
            for (int c = 0; c < r; c++)
            {
                double varSum = 0;
                for (int i = 0; i < n; i++)
                {
                    double dot = 0;
                    for (int j = 0; j < d; j++) dot += (data[i][j] - robustMean[j]) * robustSubspace[c][j];
                    varSum += dot * dot;
                }
                eigenvalues[c] = varSum / n;
            }

            double totalVar = eigenvalues.Sum(); // variance captured, not total data variance
            double[] explained = new double[r];
            for (int c = 0; c < r; c++) explained[c] = totalVar > 0 ? eigenvalues[c] / totalVar : 0;

            return new PcaResult(robustSubspace, explained, eigenvalues, robustMean, false);
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
}
