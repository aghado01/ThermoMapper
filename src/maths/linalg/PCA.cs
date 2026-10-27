using System;

namespace Maths.LinAlg
{
    public static class Pca
    {
        /// <summary>
        /// Principal components from the thin SVD of the centred data (<see cref="ThinSvd"/>), without forming the
        /// covariance or the Gram matrix: component i is the i-th right singular vector, signed so that its
        /// largest-magnitude entry is positive, and its eigenvalue is σ_i²/n, the variance with the n denominator.
        /// </summary>
        public static PcaResult Compute(double[][] data, int numComponents = 1, bool center = true, bool whiten = false)
        {
            if (data.Length < 2) throw new ArgumentException("At least 2 samples required.");

            int n = data.Length;
            int d = data[0].Length;

            double[] mean = center ? ComputeMean(data) : new double[d];
            var x = new double[n, d];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < d; j++) x[i, j] = data[i][j] - mean[j];
            ThinSvdResult svd = ThinSvd.Decompose(x);

            double totalVar = 0.0;
            foreach (double s in svd.SingularValues) totalVar += s * s / n;

            int k = Math.Max(0, Math.Min(numComponents, Math.Min(n, d)));
            var components = new double[k][];
            var explained = new double[k];
            var eigenvalues = new double[k];
            for (int i = 0; i < k; i++)
            {
                components[i] = (double[])svd.RightVectors[i].Clone();
                Orient(components[i]);
                eigenvalues[i] = svd.SingularValues[i] * svd.SingularValues[i] / n;
                explained[i] = totalVar > 0 && eigenvalues[i] > 0 ? eigenvalues[i] / totalVar : 0;
            }

            var result = new PcaResult(components, explained, eigenvalues, mean, whiten);

            if (whiten)
                WhitenInPlace(result); // scale to unit variance

            return result;
        }

        /// <summary>
        /// Signs a principal axis in place so that its largest-magnitude entry (the first, on ties) is positive: the
        /// convention <see cref="Compute"/> uses, shared with estimators that return principal axes of their own.
        /// </summary>
        public static void Orient(double[] axis)
        {
            int best = 0;
            for (int j = 1; j < axis.Length; j++) if (Math.Abs(axis[j]) > Math.Abs(axis[best])) best = j;
            if (axis.Length > 0 && axis[best] < 0.0)
                for (int j = 0; j < axis.Length; j++) axis[j] = -axis[j];
        }

        public static Func<double[], double[]> MakeProjector(PcaResult pca)
        {
            return point =>
            {
                var proj = new double[pca.Components.Length];
                for (int i = 0; i < pca.Components.Length; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < point.Length; j++)
                        sum += (point[j] - pca.Mean[j]) * pca.Components[i][j];
                    proj[i] = sum;
                }
                return proj;
            };
        }

        private static double[] ComputeMean(double[][] data)
        {
            int n = data.Length, d = data[0].Length;
            var mean = new double[d];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < d; j++)
                    mean[j] += data[i][j];
            for (int j = 0; j < d; j++) mean[j] /= n;
            return mean;
        }

        private static void WhitenInPlace(PcaResult pca)
        {
            // Scale components by 1/sqrt(eigenvalue)
            for (int i = 0; i < pca.Components.Length; i++)
            {
                double scale = 1.0 / Math.Sqrt(pca.Eigenvalues[i] + 1e-8);
                for (int j = 0; j < pca.Components[i].Length; j++)
                    pca.Components[i][j] *= scale;
            }
        }
    }

    public record PcaResult(
        double[][] Components,           // k x d
        double[] ExplainedVarianceRatio,
        double[] Eigenvalues,
        double[] Mean,
        bool Whitened
    );

}
