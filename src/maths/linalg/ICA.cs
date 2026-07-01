
using System;
using System.Linq;

namespace Maths.LinAlg
{
    public static class FastIca
    {
        public enum Mode { Deflation, Symmetric }

        /// <summary>
        /// Nonlinearity: most common and robust choice (logcosh).
        /// g(u) = tanh(u), g'(u) = 1 - tanh²(u)
        /// </summary>
        private static double G(double u) => Math.Tanh(u);
        private static double GPrime(double u) => 1.0 - Math.Tanh(u) * Math.Tanh(u);

        public static IcaResult Compute(
            double[][] data,
            int numComponents,
            Mode mode = Mode.Deflation,
            int maxIter = 200,
            double tol = 1e-6,
            int seed = 42)
        {
            if (data.Length < 2) throw new ArgumentException("Need at least 2 samples.");

            int n = data.Length;
            int d = data[0].Length;
            int k = Math.Min(numComponents, d);

            // 1. Whiten via PCA (critical preprocessing step)
            var pca = Pca.Compute(data, numComponents: k, center: true, whiten: true);
            var whitened = ApplyProjection(data, Pca.MakeProjector(pca));

            // 2. Initialize unmixing matrix W (k x k)
            var W = OrthogonalRandomInit(k, k, seed);

            // 3. Run chosen algorithm
            if (mode == Mode.Deflation)
                DeflationIteration(whitened, W, maxIter, tol);
            else
                SymmetricIteration(whitened, W, maxIter, tol);

            // 4. Reconstruct mixing matrix in original space
            var mixing = MatrixTranspose(MatrixMultiply(pca.Components /* k x d */, MatrixTranspose(W))); // W in whitened space

            return new IcaResult(
                UnmixingMatrix: W,
                MixingMatrix: mixing,
                Mean: pca.Mean,
                WhiteningPca: pca
            );
        }

        private static void DeflationIteration(double[][] X, double[][] W, int maxIter, double tol)
        {
            int k = W.Length;
            int n = X.Length;

            for (int p = 0; p < k; p++) // extract one component at a time
            {
                double[] wp = W[p]; // current weight vector

                for (int iter = 0; iter < maxIter; iter++)
                {
                    var wpOld = (double[])wp.Clone();

                    // Compute projections
                    double[] proj = new double[n];
                    for (int i = 0; i < n; i++)
                        proj[i] = Dot(X[i], wp);

                    // Update rule (one-unit FastICA)
                    double[] sum1 = new double[wp.Length];
                    double sumGPrime = 0;

                    for (int i = 0; i < n; i++)
                    {
                        double g = G(proj[i]);
                        double gp = GPrime(proj[i]);
                        for (int j = 0; j < wp.Length; j++)
                            sum1[j] += X[i][j] * g;
                        sumGPrime += gp;
                    }

                    for (int j = 0; j < wp.Length; j++)
                        wp[j] = sum1[j] / n - (sumGPrime / n) * wp[j];

                    // Orthogonalize against previous components (Gram-Schmidt)
                    for (int j = 0; j < p; j++)
                        wp = Subtract(wp, Scale(W[j], Dot(wp, W[j])));

                    // Normalize
                    wp = Normalize(wp);

                    // Convergence check
                    double diff = Math.Abs(Dot(wp, wpOld) - 1.0); // |w^T w_old - 1|
                    if (diff < tol) break;

                    W[p] = wp;
                }
            }
        }

        private static void SymmetricIteration(double[][] X, double[][] W, int maxIter, double tol)
        {
            int n = X.Length;
            int k = W.Length;

            for (int iter = 0; iter < maxIter; iter++)
            {
                var WOld = DeepClone(W);

                // Compute projections: X @ W^T
                double[][] proj = new double[n][];
                for (int i = 0; i < n; i++)
                {
                    proj[i] = new double[k];
                    for (int j = 0; j < k; j++)
                        proj[i][j] = Dot(X[i], W[j]); // assuming W rows
                }

                // Update W
                for (int j = 0; j < k; j++)
                {
                    double[] sum1 = new double[k]; // wait, actually for each row
                    double sumGp = 0;

                    for (int i = 0; i < n; i++)
                    {
                        double g = G(proj[i][j]);
                        double gp = GPrime(proj[i][j]);
                        for (int d = 0; d < k; d++) // wait, correction for matrix
                            sum1[d] += X[i][d] * g; // full vector
                        sumGp += gp;
                    }

                    // Corrected symmetric update per column/row
                    for (int d = 0; d < k; d++)
                        W[j][d] = sum1[d] / n - (sumGp / n) * W[j][d];
                }

                // Symmetric orthogonalization: W = (W W^T)^{-1/2} W
                double[][] decorrelated = SymmetricDecorrelation(W);
                for (int i = 0; i < k; i++)
                    W[i] = decorrelated[i];

                // Convergence
                double maxDiff = 0;
                for (int i = 0; i < k; i++)
                    maxDiff = Math.Max(maxDiff, Math.Abs(Dot(W[i], WOld[i]) - 1.0));

                if (maxDiff < tol) break;
            }
        }

        private static double[][] SymmetricDecorrelation(double[][] W)
        {
            // W = (W W^T)^{-1/2} W  via eigendecomposition
            int k = W.Length;
            var WWt = MatrixMultiply(W, MatrixTranspose(W));
            var eig = DenseEigen.DecomposeSymmetric(ToRectangular(WWt));

            // D^{-1/2} as jagged diagonal
            var sqrtInv = new double[k][];
            for (int i = 0; i < k; i++)
            {
                sqrtInv[i] = new double[k];
                sqrtInv[i][i] = 1.0 / Math.Sqrt(eig.Eigenvalues[i] + 1e-8);
            }

            var U = eig.Eigenvectors;
            // (U D^{-1/2} U^T) W
            return MatrixMultiply(MatrixMultiply(MatrixMultiply(U, sqrtInv), MatrixTranspose(U)), W);
        }

        // === Helper matrix/vector utilities (add to LinearAlgebra if missing) ===
        private static double Dot(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
            return sum;
        }

        private static double[] Scale(double[] v, double s)
        {
            var r = new double[v.Length];
            for (int i = 0; i < v.Length; i++) r[i] = v[i] * s;
            return r;
        }

        private static double[] Subtract(double[] a, double[] b)
        {
            var r = new double[a.Length];
            for (int i = 0; i < a.Length; i++) r[i] = a[i] - b[i];
            return r;
        }

        private static double[] Normalize(double[] v)
        {
            double norm = Math.Sqrt(Dot(v, v));
            return norm < 1e-12 ? v : Scale(v, 1.0 / norm);
        }

        private static double[][] DeepClone(double[][] m)
        {
            var c = new double[m.Length][];
            for (int i = 0; i < m.Length; i++)
                c[i] = (double[])m[i].Clone();
            return c;
        }

        private static double[][] ApplyProjection(double[][] data, Func<double[], double[]> projector)
        {
            return data.Select(projector).ToArray();
        }

        // Basic matrix ops
        private static double[][] MatrixTranspose(double[][] m)
        {
            int rows = m.Length, cols = m[0].Length;
            var t = new double[cols][];
            for (int i = 0; i < cols; i++)
            {
                t[i] = new double[rows];
                for (int j = 0; j < rows; j++) t[i][j] = m[j][i];
            }
            return t;
        }

        private static double[][] MatrixMultiply(double[][] a, double[][] b)
        {
            int m = a.Length, inner = b.Length, n = b[0].Length;
            var c = new double[m][];
            for (int i = 0; i < m; i++)
            {
                c[i] = new double[n];
                for (int l = 0; l < inner; l++)
                    for (int j = 0; j < n; j++)
                        c[i][j] += a[i][l] * b[l][j];
            }
            return c;
        }

        private static double[,] ToRectangular(double[][] matrix)
        {
            int rows = matrix.Length;
            int cols = matrix[0].Length;
            var rectangular = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    rectangular[i, j] = matrix[i][j];
            return rectangular;
        }

        private static double[][] OrthogonalRandomInit(int rows, int cols, int seed)
        {
            var rng = new Random(seed);
            var W = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                W[i] = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    // Box-Muller for standard normal
                    double u1 = 1.0 - rng.NextDouble();
                    double u2 = rng.NextDouble();
                    W[i][j] = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                }
            }
            // Gram-Schmidt orthogonalization
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < i; j++)
                    W[i] = Subtract(W[i], Scale(W[j], Dot(W[i], W[j])));
                W[i] = Normalize(W[i]);
            }
            return W;
        }
    }

    public record IcaResult(
        double[][] UnmixingMatrix,
        double[][] MixingMatrix,
        double[] Mean,
        PcaResult WhiteningPca
    );

}
