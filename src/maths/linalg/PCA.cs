using System;
using System.Linq;

namespace Maths.LinAlg
{
    public static class Pca
    {
        /// <summary>
        /// Compute PCA. Returns projection basis + diagnostics.
        /// </summary>
        public static PcaResult Compute(double[][] data, int numComponents = 1, bool center = true, bool whiten = false)
        {
            if (data.Length < 2) throw new ArgumentException("At least 2 samples required.");

            int n = data.Length;
            int d = data[0].Length;

            double[] mean = center ? ComputeMean(data) : new double[d];
            var centered = center ? Center(data, mean) : data;

            int k = Math.Min(numComponents, Math.Min(n, d));
            var components = new double[k][];
            var explained = new double[k];
            var eigenvalues = new double[k];
            double totalVar = 0.0;

            if (d <= n)
            {
                // Standard PCA: O(n d^2) Covariance matrix
                double[,] cov = ComputeCovariance(centered, n);
                var eig = DenseEigen.DecomposeSymmetric(cov);

                totalVar = eig.Eigenvalues.Where(val => val > 0).Sum();
                for (int i = 0; i < k; i++)
                {
                    components[i] = eig.Eigenvectors[i];
                    eigenvalues[i] = eig.Eigenvalues[i];
                    explained[i] = totalVar > 0 && eigenvalues[i] > 0 ? eigenvalues[i] / totalVar : 0;
                }
            }
            else
            {
                // Fast Dual PCA (SVD equivalent): O(d n^2) Gram matrix
                double[,] gram = ComputeGramMatrix(centered, n);
                var eig = DenseEigen.DecomposeSymmetric(gram);

                totalVar = eig.Eigenvalues.Where(val => val > 0).Sum();
                for (int i = 0; i < k; i++)
                {
                    eigenvalues[i] = eig.Eigenvalues[i];
                    explained[i] = totalVar > 0 && eigenvalues[i] > 0 ? eigenvalues[i] / totalVar : 0;
                    
                    double[] v = new double[d];
                    double lambda = eigenvalues[i];
                    
                    if (lambda > 1e-12)
                    {
                        double scale = 1.0 / Math.Sqrt(n * lambda);
                        for (int row = 0; row < n; row++)
                        {
                            double u_val = eig.Eigenvectors[i][row];
                            for (int col = 0; col < d; col++)
                            {
                                v[col] += centered[row][col] * u_val * scale;
                            }
                        }
                    }
                    components[i] = v;
                }
            }

            var result = new PcaResult(components, explained, eigenvalues, mean, whiten);

            if (whiten)
                WhitenInPlace(result); // scale to unit variance

            return result;
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

        // Helper methods (mean, center, cov) — standard and efficient
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

        private static double[][] Center(double[][] data, double[] mean)
        {
            int n = data.Length, d = data[0].Length;
            var c = new double[n][];
            for (int i = 0; i < n; i++)
            {
                c[i] = new double[d];
                for (int j = 0; j < d; j++)
                    c[i][j] = data[i][j] - mean[j];
            }
            return c;
        }

        private static double[,] ComputeGramMatrix(double[][] centered, int n)
        {
            var gram = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    double sum = MatrixOps.Dot(centered[i], centered[j]);
                    gram[i, j] = gram[j, i] = sum / n;
                }
            }
            return gram;
        }

        private static double[,] ComputeCovariance(double[][] centered, int n)
        {
            int d = centered[0].Length;
            var cov = new double[d, d];
            for (int i = 0; i < d; i++)
            {
                for (int j = i; j < d; j++) // symmetry
                {
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                        sum += centered[k][i] * centered[k][j];
                    cov[i, j] = cov[j, i] = sum / n;
                }
            }
            return cov;
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
