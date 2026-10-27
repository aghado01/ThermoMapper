using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// A thin singular value decomposition A = U·diag(σ)·Vᵀ, with k = min(m, n) singular values in descending order.
    /// </summary>
    /// <param name="SingularValues">σ₁ ≥ σ₂ ≥ … ≥ σ_k ≥ 0.</param>
    /// <param name="LeftVectors">k orthonormal left singular vectors, each of length m.</param>
    /// <param name="RightVectors">k orthonormal right singular vectors, each of length n.</param>
    /// <param name="NumericalRank">The number of singular values above <paramref name="RankTolerance"/>.</param>
    /// <param name="RankTolerance">max(m, n)·ε·σ₁, ε the machine epsilon: the threshold LAPACK-based tools (R's
    /// <c>Matrix::rankMatrix</c>, NumPy's <c>matrix_rank</c>) use, and the BASS truncation oracle's.</param>
    public sealed record ThinSvdResult(
        double[] SingularValues, double[][] LeftVectors, double[][] RightVectors, int NumericalRank, double RankTolerance);

    /// <summary>
    /// Thin SVD by Householder QR and one-sided Jacobi rotations on the triangular factor (Hestenes 1958; Drmač and
    /// Veselić 2008). It never forms AᵀA, so singular values far below σ₁ keep their absolute accuracy of about
    /// ε·σ₁ rather than the √ε·σ₁ a Gram eigendecomposition leaves them. Left vectors whose singular value is exactly
    /// zero are completed to an orthonormal set. The response decomposition under BASS's PCA basis reads its rank from
    /// <see cref="ThinSvdResult.NumericalRank"/>.
    /// </summary>
    public static class ThinSvd
    {
        private const int MaxSweeps = 60;
        private const double MachineEpsilon = 2.220446049250313e-16;   // 2⁻⁵², R's .Machine$double.eps

        /// <summary>Decomposes the m × n matrix <paramref name="a"/>.</summary>
        public static ThinSvdResult Decompose(double[,] a)
        {
            ArgumentNullException.ThrowIfNull(a);
            int m = a.GetLength(0), n = a.GetLength(1);
            if (m >= n)
            {
                var columns = new double[n][];
                for (int j = 0; j < n; j++)
                {
                    columns[j] = new double[m];
                    for (int i = 0; i < m; i++) columns[j][i] = a[i, j];
                }
                return Tall(columns, m, n, transposed: false);
            }
            else
            {
                var columns = new double[m][];     // the columns of Aᵀ, which are A's rows
                for (int i = 0; i < m; i++)
                {
                    columns[i] = new double[n];
                    for (int j = 0; j < n; j++) columns[i][j] = a[i, j];
                }
                return Tall(columns, n, m, transposed: true);
            }
        }

        // The SVD of the rows × cols matrix whose columns are given (rows ≥ cols); for Aᵀ, U and V swap at the end.
        private static ThinSvdResult Tall(double[][] x, int rows, int cols, bool transposed)
        {
            // Householder QR: x's columns become R's columns, and the reflectors are kept.
            var reflectors = new double[cols][];
            for (int j = 0; j < cols; j++)
            {
                double[] xj = x[j];
                double norm = Math.Sqrt(PairwiseSum.Dot(xj.AsSpan(j), xj.AsSpan(j)));
                if (norm == 0.0) continue;
                double alpha = xj[j] >= 0.0 ? -norm : norm;
                var v = new double[rows - j];
                xj.AsSpan(j).CopyTo(v);
                v[0] -= alpha;                                 // x_j − α, with no cancellation
                double vv = PairwiseSum.Dot(v, v);
                reflectors[j] = v;
                xj[j] = alpha;
                Array.Clear(xj, j + 1, rows - j - 1);
                for (int k = j + 1; k < cols; k++) Reflect(v, vv, x[k].AsSpan(j));
            }

            // One-sided Jacobi on R: rotate column pairs until every pair is orthogonal to working precision.
            var w = new double[cols][];
            var vt = new double[cols][];
            for (int j = 0; j < cols; j++)
            {
                w[j] = new double[cols];
                for (int i = 0; i <= j; i++) w[j][i] = x[j][i];
                vt[j] = new double[cols];
                vt[j][j] = 1.0;
            }
            double orthogonal = cols * MachineEpsilon;      // |cos| of the angle below which a pair counts as orthogonal
            bool converged = false;
            for (int sweep = 0; sweep < MaxSweeps && !converged; sweep++)
            {
                converged = true;
                for (int p = 0; p < cols - 1; p++)
                    for (int q = p + 1; q < cols; q++)
                    {
                        double alpha = PairwiseSum.Dot(w[p], w[p]);
                        double beta = PairwiseSum.Dot(w[q], w[q]);
                        double gamma = PairwiseSum.Dot(w[p], w[q]);
                        if (alpha == 0.0 || beta == 0.0 || Math.Abs(gamma) <= orthogonal * Math.Sqrt(alpha) * Math.Sqrt(beta))
                            continue;
                        converged = false;
                        double zeta = (beta - alpha) / (2.0 * gamma);
                        double t = (zeta >= 0.0 ? 1.0 : -1.0) / (Math.Abs(zeta) + Math.Sqrt(1.0 + zeta * zeta));
                        double c = 1.0 / Math.Sqrt(1.0 + t * t);
                        double s = c * t;
                        Rotate(w[p], w[q], c, s);
                        Rotate(vt[p], vt[q], c, s);
                    }
            }
            if (!converged) throw new InvalidOperationException($"One-sided Jacobi did not converge in {MaxSweeps} sweeps.");

            // σ_j = ‖w_j‖, in descending order; the left vectors of R are w_j/σ_j.
            var sigma = new double[cols];
            for (int j = 0; j < cols; j++) sigma[j] = Math.Sqrt(PairwiseSum.Dot(w[j], w[j]));
            var order = new int[cols];
            for (int j = 0; j < cols; j++) order[j] = j;
            Array.Sort(order, (i, j) => sigma[j].CompareTo(sigma[i]) != 0 ? sigma[j].CompareTo(sigma[i]) : i.CompareTo(j));

            var singularValues = new double[cols];
            var left = new double[cols][];
            var right = new double[cols][];
            int completeFrom = cols;
            for (int k = 0; k < cols; k++)
            {
                int j = order[k];
                singularValues[k] = sigma[j];
                right[k] = vt[j];
                var u = new double[cols];
                if (sigma[j] > 0.0)
                    for (int i = 0; i < cols; i++) u[i] = w[j][i] / sigma[j];
                else if (completeFrom == cols)
                    completeFrom = k;
                left[k] = u;
            }
            for (int k = completeFrom; k < cols; k++) Complete(left, k, cols);

            // U = Q·[U_R; 0], applying the reflectors in reverse.
            for (int k = 0; k < cols; k++)
            {
                var u = new double[rows];
                left[k].CopyTo(u, 0);
                for (int j = cols - 1; j >= 0; j--)
                {
                    double[] v = reflectors[j];
                    if (v is null) continue;
                    Reflect(v, PairwiseSum.Dot(v, v), u.AsSpan(j));
                }
                left[k] = u;
            }

            double tolerance = Math.Max(rows, cols) * MachineEpsilon * (cols > 0 ? singularValues[0] : 0.0);
            int rank = 0;
            while (rank < cols && singularValues[rank] > tolerance) rank++;
            return transposed
                ? new ThinSvdResult(singularValues, right, left, rank, tolerance)
                : new ThinSvdResult(singularValues, left, right, rank, tolerance);
        }

        // y ← (I − 2vvᵀ/vᵀv)·y.
        private static void Reflect(double[] v, double vv, Span<double> y)
        {
            double f = 2.0 * PairwiseSum.Dot(v, y) / vv;
            for (int i = 0; i < v.Length; i++) y[i] -= f * v[i];
        }

        private static void Rotate(double[] p, double[] q, double c, double s)
        {
            for (int i = 0; i < p.Length; i++)
            {
                double a = p[i], b = q[i];
                p[i] = c * a - s * b;
                q[i] = s * a + c * b;
            }
        }

        // Fills left[k] with a unit vector orthogonal to left[0..k), from the standard basis vector that keeps the most
        // after projection; two passes of Gram–Schmidt.
        private static void Complete(double[][] left, int k, int n)
        {
            var best = new double[n];
            double bestNorm = -1.0;
            var candidate = new double[n];
            for (int e = 0; e < n && bestNorm < 0.5; e++)
            {
                Array.Clear(candidate);
                candidate[e] = 1.0;
                for (int pass = 0; pass < 2; pass++)
                    for (int i = 0; i < k; i++)
                    {
                        double projection = PairwiseSum.Dot(candidate, left[i]);
                        for (int r = 0; r < n; r++) candidate[r] -= projection * left[i][r];
                    }
                double norm = Math.Sqrt(PairwiseSum.Dot(candidate, candidate));
                if (norm > bestNorm)
                {
                    bestNorm = norm;
                    candidate.CopyTo(best, 0);
                }
            }
            for (int r = 0; r < n; r++) best[r] /= bestNorm;
            left[k] = best;
        }
    }
}
