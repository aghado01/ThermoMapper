using System;
using Maths.LinAlg;

namespace Maths.Distance
{
    public static class Cosine
    {
        public static double Similarity(ReadOnlySpan<double> u, ReadOnlySpan<double> v)
        {
            if (u.Length != v.Length) throw new ArgumentException("Lengths must match");
            double dot = 0;
            double normU = 0;
            double normV = 0;
            for (int i = 0; i < u.Length; i++)
            {
                dot += u[i] * v[i];
                normU += u[i] * u[i];
                normV += v[i] * v[i];
            }
            if (normU == 0 || normV == 0) return 0;
            return dot / (Math.Sqrt(normU) * Math.Sqrt(normV));
        }

        public static double Distance(ReadOnlySpan<double> u, ReadOnlySpan<double> v)
        {
            return 1.0 - Similarity(u, v);
        }

        /// <summary>
        /// A hybrid measure combining angular similarity (Cosine) with radial distance (Magnitude).
        /// As suggested in 2504.16318, combining norm ratio with cosine addresses frequency bias.
        /// </summary>
        public static double NormAwareSimilarity(ReadOnlySpan<double> u, ReadOnlySpan<double> v, double alpha = 0.5)
        {
            double sim = Similarity(u, v);
            double normU = 0;
            double normV = 0;
            for (int i = 0; i < u.Length; i++)
            {
                normU += u[i] * u[i];
                normV += v[i] * v[i];
            }
            normU = Math.Sqrt(normU);
            normV = Math.Sqrt(normV);
            double maxNorm = Math.Max(normU, normV);
            double minNorm = Math.Min(normU, normV);
            
            // Radial similarity is 1.0 if norms match exactly, decreases as they diverge
            double radialSim = maxNorm > 0 ? (minNorm / maxNorm) : 1.0;
            
            return alpha * sim + (1 - alpha) * radialSim;
        }

        /// <summary>
        /// Post-Hoc Isotropization: Modifies embeddings in-place by ablating the top K principal components.
        /// This forces the embedding space to be more isotropic, drastically improving the utility of Cosine similarity 
        /// by stripping out dominant common-mode directions (e.g., word frequency syntax).
        /// </summary>
        /// <param name="embeddings">The embeddings to isotropize</param>
        /// <param name="topComponents">The top principal components (eigenvectors) to ablate</param>
        public static void IsotropizeInPlace(double[][] embeddings, double[][] topComponents)
        {
            int n = embeddings.Length;
            if (n == 0) return;
            int d = embeddings[0].Length;
            int k = topComponents.Length;
            
            for (int i = 0; i < n; i++)
            {
                for (int c = 0; c < k; c++)
                {
                    // Strip the projection onto component c via the audited SIMD kernels.
                    // Assumes topComponents are orthonormal (PCA eigenvectors), so the
                    // successive ablations commute and the loop order is immaterial.
                    double dot = MatrixOps.Dot(embeddings[i], topComponents[c]);
                    MatrixOps.ScaledSubtract(embeddings[i], topComponents[c], dot, d);
                }
            }
        }
    }
}
