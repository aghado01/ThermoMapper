using System;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

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
    /// ε·σ₁ rather than the √ε·σ₁ a Gram eigendecomposition leaves them. A column of R below ε·‖R‖_F is zero to
    /// working precision: rotations leave it alone, and its left vector is completed to an orthonormal set rather than
    /// read from its rounding noise. Norms and rotation angles are computed without underflow or overflow. The response decomposition under BASS's PCA basis reads its rank from
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
            foreach (double value in a)
                if (!double.IsFinite(value)) throw new ArgumentException("The matrix has a non-finite entry.", nameof(a));
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
                double norm = Norm(xj.AsSpan(j));
                if (norm == 0.0) continue;
                double alpha = xj[j] >= 0.0 ? -norm : norm;
                var v = new double[rows - j];
                xj.AsSpan(j).CopyTo(v);
                v[0] -= alpha;                                 // x_j − α, with no cancellation
                double scale = 0.0;                            // the reflector does not see v's scale, so keep vᵀv in range
                foreach (double e in v) scale = Math.Max(scale, Math.Abs(e));
                for (int i = 0; i < v.Length; i++) v[i] /= scale;
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
            // A column below ε·‖R‖_F is zero to working precision: rotating it only chases rounding noise, and its
            // direction carries none of A. ‖R‖_F is invariant under the rotations.
            double frobenius = 0.0;
            foreach (double[] column in w) frobenius = double.Hypot(frobenius, Norm(column));
            double negligible = MachineEpsilon * frobenius;
            var scratchP = new double[cols];
            var scratchQ = new double[cols];
            bool converged = false;
            for (int sweep = 0; sweep < MaxSweeps && !converged; sweep++)
            {
                converged = true;
                for (int p = 0; p < cols - 1; p++)
                    for (int q = p + 1; q < cols; q++)
                    {
                        double np = Norm(w[p]), nq = Norm(w[q]);
                        if (np <= negligible || nq <= negligible) continue;
                        double cos = Cosine(w[p], np, w[q], nq, scratchP, scratchQ);
                        if (Math.Abs(cos) <= orthogonal) continue;
                        // The rotation that diagonalizes the pair's 2 × 2 Gram: t = tan θ solves t² + 2ζt − 1 = 0 with
                        // ζ = (‖w_q‖² − ‖w_p‖²)/(2·w_pᵀw_q). Written through ρ = min/max of the two norms,
                        // t = sign(ζ)·ρ/(|ζρ| + hypot(ρ, ζρ)), nothing overflows however far apart the norms are.
                        double rho = Math.Min(np, nq) / Math.Max(np, nq);
                        double zetaRho = (1.0 - rho * rho) / (2.0 * Math.Abs(cos));
                        double sign = (nq >= np) == (cos > 0.0) ? 1.0 : -1.0;
                        double t = sign * rho / (zetaRho + double.Hypot(rho, zetaRho));
                        double c = 1.0 / Math.Sqrt(1.0 + t * t);
                        double s = c * t;
                        if (s == 0.0) continue;                    // below what a rotation can resolve
                        converged = false;
                        Rotate(w[p], w[q], c, s);
                        Rotate(vt[p], vt[q], c, s);
                    }
            }
            if (!converged) throw new InvalidOperationException($"One-sided Jacobi did not converge in {MaxSweeps} sweeps.");

            // σ_j = ‖w_j‖, in descending order; the left vectors of R are w_j/σ_j, except for negligible columns.
            var sigma = new double[cols];
            for (int j = 0; j < cols; j++) sigma[j] = Norm(w[j]);
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
                if (sigma[j] > negligible)
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

        // ‖x‖ without underflow or overflow: the plain sum of squares when it sits well inside the normal range,
        // otherwise the sum for x/max|x_i|.
        private static double Norm(ReadOnlySpan<double> x)
        {
            double squares = PairwiseSum.Dot(x, x);
            if (squares > 1e-280 && squares < 1e280) return Math.Sqrt(squares);
            double max = 0.0;
            foreach (double e in x) max = Math.Max(max, Math.Abs(e));
            if (max == 0.0) return 0.0;
            var scaled = new double[x.Length];
            for (int i = 0; i < x.Length; i++) scaled[i] = x[i] / max;
            return max * Math.Sqrt(PairwiseSum.Dot(scaled, scaled));
        }

        // The cosine of the angle between x and y, given their norms: the plain dot product when both norms are well
        // inside the normal range, otherwise the dot product of the normalized vectors.
        private static double Cosine(double[] x, double nx, double[] y, double ny, double[] sx, double[] sy)
        {
            if (nx > 1e-140 && nx < 1e140 && ny > 1e-140 && ny < 1e140) return PairwiseSum.Dot(x, y) / nx / ny;
            for (int i = 0; i < x.Length; i++)
            {
                sx[i] = x[i] / nx;
                sy[i] = y[i] / ny;
            }
            return PairwiseSum.Dot(sx, sy);
        }

        // y ← (I − 2vvᵀ/vᵀv)·y.
        private static void Reflect(double[] v, double vv, Span<double> y)
        {
            double f = 2.0 * PairwiseSum.Dot(v, y) / vv;
            MatrixOps.ScaledSubtract(y, v, f, v.Length);
        }

        // (p, q) ← (c·p − s·q, s·p + c·q), elementwise, so vector width does not change the bits.
        private static void Rotate(double[] p, double[] q, double c, double s)
        {
            int n = p.Length, i = 0;
            ref double pr = ref MemoryMarshal.GetArrayDataReference(p);
            ref double qr = ref MemoryMarshal.GetArrayDataReference(q);
            var vc = Vector256.Create(c);
            var vs = Vector256.Create(s);
            for (; i + 3 < n; i += 4)
            {
                var a = Vector256.LoadUnsafe(ref pr, (nuint)i);
                var b = Vector256.LoadUnsafe(ref qr, (nuint)i);
                Vector256.Subtract(Vector256.Multiply(vc, a), Vector256.Multiply(vs, b)).StoreUnsafe(ref pr, (nuint)i);
                Vector256.Add(Vector256.Multiply(vs, a), Vector256.Multiply(vc, b)).StoreUnsafe(ref qr, (nuint)i);
            }
            for (; i < n; i++)
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
