using System;
using System.Linq;

namespace Maths.LinAlg
{
    /// <summary>
    /// Symmetric eigenvalue decomposition (Jacobi + refinement via Cholesky if needed).
    /// Suitable for covariance matrices in PCA/ICA.
    /// </summary>
    public static class Eigen
    {
        /// <summary>
        /// A pivot rotates when |a_pq| exceeds <paramref name="tol"/>; a sweep that rotates nothing ends the solve. The result
        /// reports how the solve ended, read from the diagonalization it leaves.
        /// </summary>
        public static EigenResult DecomposeSymmetric(double[,] matrix, int maxSweeps = 256, double tol = 1e-12)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n) throw new ArgumentException("Matrix must be square.");
            DecompositionLog.Record(DecompositionKernel.ReferenceJacobi, n, n);

            var A = (double[,])matrix.Clone(); // working copy
            var V = Identity(n);               // eigenvectors
            int sweeps = 0;
            long rotations = 0;

            for (int sweep = 0; sweep < maxSweeps; sweep++)
            {
                sweeps = sweep + 1;
                long sweepRotations = 0;

                // Classical cyclic Jacobi sweep over all off-diagonal pairs.
                for (int p = 0; p < n - 1; p++)
                {
                    for (int q = p + 1; q < n; q++)
                    {
                        if (Math.Abs(A[p, q]) <= tol) continue;

                        Rotate(A, V, p, q, n);
                        sweepRotations++;
                    }
                }

                rotations += sweepRotations;
                if (sweepRotations == 0) break;
            }
            (double offDiagonalNorm, double largestOffDiagonal) = DenseEigenConvergence.OffDiagonal(A);
            var convergence = new DenseEigenConvergence(sweeps, rotations, offDiagonalNorm, largestOffDiagonal, largestOffDiagonal <= tol);

            // Extract eigenvalues + sort descending
            var eigenvalues = new double[n];
            for (int i = 0; i < n; i++) eigenvalues[i] = A[i, i];

            var idx = Enumerable.Range(0, n).OrderByDescending(i => eigenvalues[i]).ToArray();

            var sortedEig = new double[n];
            var sortedVec = new double[n][];

            for (int i = 0; i < n; i++)
            {
                sortedEig[i] = eigenvalues[idx[i]];
                sortedVec[i] = new double[n];
                for (int j = 0; j < n; j++)
                    sortedVec[i][j] = V[j, idx[i]]; // columns
            }

            return new EigenResult(sortedEig, sortedVec) { Convergence = convergence };
        }

        private static void Rotate(double[,] A, double[,] V, int p, int q, int n)
        {
            double theta = 0.5 * Math.Atan2(2 * A[p, q], A[q, q] - A[p, p]);
            double c = Math.Cos(theta);
            double s = Math.Sin(theta);

            // Update A (symmetric)
            double app = A[p, p], aqq = A[q, q], apq = A[p, q];
            A[p, p] = c * c * app - 2 * s * c * apq + s * s * aqq;
            A[q, q] = s * s * app + 2 * s * c * apq + c * c * aqq;
            A[p, q] = A[q, p] = 0;

            for (int i = 0; i < n; i++)
            {
                if (i == p || i == q) continue;
                double aip = A[i, p], aiq = A[i, q];
                A[i, p] = A[p, i] = c * aip - s * aiq;
                A[i, q] = A[q, i] = s * aip + c * aiq;
            }

            // Update eigenvectors
            for (int i = 0; i < n; i++)
            {
                double vip = V[i, p], viq = V[i, q];
                V[i, p] = c * vip - s * viq;
                V[i, q] = s * vip + c * viq;
            }
        }

        private static double[,] Identity(int n)
        {
            var I = new double[n, n];
            for (int i = 0; i < n; i++) I[i, i] = 1.0;
            return I;
        }
    }

    public record EigenResult(double[] Eigenvalues, double[][] Eigenvectors) // Eigenvectors as rows for convenience in PCA
    {
        /// <summary>How the solve that produced this result ended; the default when its producer reports none.</summary>
        public DenseEigenConvergence Convergence { get; init; }
    }

}
