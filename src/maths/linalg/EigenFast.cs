// Maths/LinAlg/EigenFast_v3.cs
#nullable enable
using System;
using System.Buffers;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

namespace Maths.LinAlg
{
    public static unsafe class EigenFast
    {
        /// <summary>
        /// v3: High-Performance Upper-Triangle Algebraic Cyclic Jacobi Solver.
        /// Zero transcendentals, zero symmetry mirroring overhead, fully vectorized.
        /// A pivot rotates when |a_pq| > tol √|a_pp a_qq|, the element-wise relative rule (<see cref="JacobiPivot"/>). A sweep that
        /// rotates nothing ends the solve. An unrepresentable required rotation leaves Converged false on the final working state.
        /// Eigenvalues are written descending; eigenvectors as the columns of the flat column-major <paramref name="outEigenvectors"/>,
        /// or none when it is empty, which skips their rotations, about half the work, and leaves the eigenvalues the same bits. The
        /// working copy and the sort's scratch are pooled, so a call into caller-owned buffers allocates nothing once the pool is warm.
        /// </summary>
        public static DenseEigenConvergence DecomposeSymmetricInPlace(
            ReadOnlySpan<double> flatMatrix,
            int n,
            Span<double> outEigenvalues,
            Span<double> outEigenvectors,
            int maxSweeps = 256,
            double tol = 1e-12)
        {
            JacobiPivot.Validate(tol);
            if (flatMatrix.Length != n * n) throw new ArgumentException("Matrix dimensions must match N x N.");
            if (outEigenvalues.Length != n) throw new ArgumentException("Eigenvalues buffer size must equal N.");
            bool vectors = outEigenvectors.Length != 0;
            if (vectors && outEigenvectors.Length != n * n) throw new ArgumentException("Eigenvectors buffer size must equal N * N, or 0 for eigenvalues alone.");

            DecompositionLog.Record(DecompositionKernel.FastJacobi, n, n);
            int matrixSize = n * n;
            double[] workingArray = ArrayPool<double>.Shared.Rent(matrixSize);
            int sweeps = 0;
            long rotations = 0;
            double offDiagonalNorm, largestOffDiagonal;
            bool converged = true;
            bool stoppedWithoutRotations = false;

            try
            {
                fixed (double* pSrc = flatMatrix)
                fixed (double* A = workingArray)
                fixed (double* V = outEigenvectors)
                {
                    Unsafe.CopyBlock(A, pSrc, (uint)(sizeof(double) * matrixSize));

                    // Initialize Eigenvectors matrix to Identity (flat column-major)
                    if (vectors)
                    {
                        outEigenvectors.Clear();
                        for (int i = 0; i < n; i++) V[i * n + i] = 1.0;
                    }

                    for (int sweep = 0; sweep < maxSweeps; sweep++)
                    {
                        sweeps = sweep + 1;
                        long sweepRotations = 0;

                        // Sweep all unique pairs (p, q) over the strict upper triangle
                        for (int p = 0; p < n - 1; p++)
                        {
                            for (int q = p + 1; q < n; q++)
                            {
                                // In column-major layout, row p, col q is at location [q * n + p]
                                double apq = A[q * n + p];

                                // Compute Algebraic Givens Rotations (No Atan2, Cos, or Sin)
                                int pBase = p * n;
                                int qBase = q * n;
                                double app = A[pBase + p];
                                double aqq = A[qBase + q];
                                if (JacobiPivot.Passes(app, aqq, apq, tol)) continue;
                                if (!JacobiPivot.TryTangent(app, aqq, apq, out double t)) continue;
                                sweepRotations++;

                                double c = 1.0 / Math.Sqrt(1.0 + t * t);
                                double s = t * c;
                                bool canVectorize = Vector256.IsHardwareAccelerated;
                                Vector256<double> vc = canVectorize ? Vector256.Create(c) : default;
                                Vector256<double> vs = canVectorize ? Vector256.Create(s) : default;

                                // Update the tracking basis vector matrix V (Fully Vectorized Column Transformations)
                                int i = 0;
                                if (vectors && canVectorize && n >= 4)
                                {
                                    int simdLen = n & ~3;

                                    for (; i < simdLen; i += 4)
                                    {
                                        var vvip = Vector256.Load(V + pBase + i);
                                        var vviq = Vector256.Load(V + qBase + i);

                                        var vNewVip = Vector256.Subtract(Vector256.Multiply(vc, vvip), Vector256.Multiply(vs, vviq));
                                        var vNewViq = Vector256.Add(Vector256.Multiply(vs, vvip), Vector256.Multiply(vc, vviq));

                                        Vector256.Store(vNewVip, V + pBase + i);
                                        Vector256.Store(vNewViq, V + qBase + i);
                                    }
                                }
                                if (vectors)
                                    for (; i < n; i++)
                                    {
                                        double vip = V[pBase + i];
                                        double viq = V[qBase + i];
                                        V[pBase + i] = c * vip - s * viq;
                                        V[qBase + i] = s * vip + c * viq;
                                    }

                                int zone1 = 0;
                                if (canVectorize && p >= 4)
                                {
                                    int simdLen = p & ~3;

                                    for (; zone1 < simdLen; zone1 += 4)
                                    {
                                        var vaip = Vector256.Load(A + pBase + zone1);
                                        var vaiq = Vector256.Load(A + qBase + zone1);

                                        var vNewAip = Vector256.Subtract(Vector256.Multiply(vc, vaip), Vector256.Multiply(vs, vaiq));
                                        var vNewAiq = Vector256.Add(Vector256.Multiply(vs, vaip), Vector256.Multiply(vc, vaiq));

                                        Vector256.Store(vNewAip, A + pBase + zone1);
                                        Vector256.Store(vNewAiq, A + qBase + zone1);
                                    }
                                }

                                for (; zone1 < p; zone1++)
                                {
                                    int idx_ip = pBase + zone1;
                                    int idx_iq = qBase + zone1;

                                    double aip = A[idx_ip];
                                    double aiq = A[idx_iq];

                                    A[idx_ip] = c * aip - s * aiq;
                                    A[idx_iq] = s * aip + c * aiq;
                                }

                                for (i = p + 1; i < q; i++)
                                {
                                    int idx_ip = i * n + p;
                                    int idx_iq = qBase + i;

                                    double aip = A[idx_ip];
                                    double aiq = A[idx_iq];

                                    A[idx_ip] = c * aip - s * aiq;
                                    A[idx_iq] = s * aip + c * aiq;
                                }

                                for (i = q + 1; i < n; i++)
                                {
                                    int rowBase = i * n;
                                    int idx_ip = rowBase + p;
                                    int idx_iq = rowBase + q;

                                    double aip = A[idx_ip];
                                    double aiq = A[idx_iq];

                                    A[idx_ip] = c * aip - s * aiq;
                                    A[idx_iq] = s * aip + c * aiq;
                                }

                                // Authoritative update of the pivot intersections
                                A[pBase + p] = app - t * apq;
                                A[qBase + q] = aqq + t * apq;
                                A[qBase + p] = 0.0; // Clear the upper-triangle element explicitly
                            }
                        }

                        rotations += sweepRotations;
                        // Every pivot passed the test: the next sweep would rotate nothing either.
                        if (sweepRotations == 0)
                        {
                            stoppedWithoutRotations = true;
                            break;
                        }
                    }

                    // Extract final parsed eigenvalues from the diagonal trace
                    for (int i = 0; i < n; i++) outEigenvalues[i] = A[i * n + i];
                }

                // The returned diagonalization's off-diagonal, read after the last sweep: a sample taken during a sweep precedes
                // that sweep's rotations.
                (offDiagonalNorm, largestOffDiagonal) = DenseEigenConvergence.OffDiagonal(workingArray.AsSpan(0, matrixSize), n);
                for (int q = 1; q < n; q++)
                    for (int p = 0; p < q; p++)
                        converged &= JacobiPivot.Passes(workingArray[p * n + p], workingArray[q * n + q], workingArray[q * n + p], tol);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(workingArray);
            }

            SortDescending(outEigenvalues, outEigenvectors, n);
            return new DenseEigenConvergence(sweeps, rotations, offDiagonalNorm, largestOffDiagonal, converged)
            { Stalled = stoppedWithoutRotations && !converged };
        }

        // Eigenvalues descending with their eigenvector columns. The exchange order is the one this solver always used, so a repeated
        // eigenvalue's vectors keep their order; the exchanges run on an index permutation and the columns move once, following its
        // cycles, O(n²) where swapping whole columns at each exchange cost up to O(n³).
        private static void SortDescending(Span<double> values, Span<double> vectors, int n)
        {
            int[] order = ArrayPool<int>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) order[i] = i;
                for (int i = 0; i < n - 1; i++)
                    for (int j = i + 1; j < n; j++)
                        if (values[j] > values[i])
                        {
                            (values[i], values[j]) = (values[j], values[i]);
                            (order[i], order[j]) = (order[j], order[i]);
                        }
                if (vectors.Length == 0) return;

                // Position i takes the column that started at order[i]; each cycle is walked once, its first column held aside.
                double[] held = ArrayPool<double>.Shared.Rent(n);
                try
                {
                    for (int start = 0; start < n; start++)
                    {
                        if (order[start] < 0 || order[start] == start) continue;
                        vectors.Slice(start * n, n).CopyTo(held);
                        int at = start;
                        while (true)
                        {
                            int from = order[at];
                            order[at] = -1;
                            if (from == start)
                            {
                                held.AsSpan(0, n).CopyTo(vectors.Slice(at * n, n));
                                break;
                            }
                            vectors.Slice(from * n, n).CopyTo(vectors.Slice(at * n, n));
                            at = from;
                        }
                    }
                }
                finally
                {
                    ArrayPool<double>.Shared.Return(held);
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(order);
            }
        }

        /// <summary>
        /// Decomposes a symmetric row-major <c>double[,]</c> matrix and returns eigenvectors as jagged rows.
        /// The input is flattened to the fast core's column-major layout before solving.
        /// </summary>
        public static EigenResult DecomposeSymmetric(
            ReadOnlySpan<double> flatColumnMajorMatrix,
            int n,
            int maxSweeps = 256,
            double tol = 1e-12)
        {
            ArrayPool<double> pool = ArrayPool<double>.Shared;
            double[] outEigVecs = pool.Rent(n * n);
            double[] outEigVals = new double[n];

            try
            {
                DenseEigenConvergence convergence = DecomposeSymmetricInPlace(
                    flatColumnMajorMatrix,
                    n,
                    outEigVals,
                    outEigVecs.AsSpan(0, n * n),
                    maxSweeps,
                    tol);

                var sortedVec = new double[n][];
                for (int i = 0; i < n; i++)
                {
                    sortedVec[i] = new double[n];
                    outEigVecs.AsSpan(i * n, n).CopyTo(sortedVec[i]);
                }

                return new EigenResult(outEigVals, sortedVec) { Convergence = convergence };
            }
            finally
            {
                pool.Return(outEigVecs);
            }
        }

        public static EigenResult DecomposeSymmetric(
            double[,] matrix,
            int maxSweeps = 256,
            double tol = 1e-12)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n)
                throw new ArgumentException("Matrix must be square.", nameof(matrix));

            ArrayPool<double> pool = ArrayPool<double>.Shared;
            double[] flatMatrix = pool.Rent(n * n);

            try
            {
                for (int c = 0; c < n; c++)
                {
                    for (int r = 0; r < n; r++)
                        flatMatrix[c * n + r] = matrix[r, c];
                }

                return DecomposeSymmetric(flatMatrix.AsSpan(0, n * n), n, maxSweeps, tol);
            }
            finally
            {
                pool.Return(flatMatrix);
            }
        }
    }
}
