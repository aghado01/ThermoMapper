// Maths/LinAlg/EigenFast_v3.cs
#nullable enable
using System;
using System.Buffers;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace Maths.LinAlg
{
    public static unsafe class EigenFast
    {
        /// <summary>
        /// v3: High-Performance Upper-Triangle Algebraic Cyclic Jacobi Solver.
        /// Zero transcendentals, zero symmetry mirroring overhead, fully vectorized.
        /// </summary>
        public static void DecomposeSymmetricInPlace(
            ReadOnlySpan<double> flatMatrix,
            int n,
            Span<double> outEigenvalues,
            Span<double> outEigenvectors,
            int maxSweeps = 256,
            double tol = 1e-12,
            DenseEigenFastVariant fastVariant = DenseEigenFastVariant.Default)
        {
            if (!Enum.IsDefined(fastVariant))
                throw new ArgumentOutOfRangeException(nameof(fastVariant));

            if (flatMatrix.Length != n * n) throw new ArgumentException("Matrix dimensions must match N x N.");
            if (outEigenvalues.Length != n) throw new ArgumentException("Eigenvalues buffer size must equal N.");
            if (outEigenvectors.Length != n * n) throw new ArgumentException("Eigenvectors buffer size must equal N * N.");

            int matrixSize = n * n;
            double[] workingArray = ArrayPool<double>.Shared.Rent(matrixSize);
            bool useFma = fastVariant == DenseEigenFastVariant.Fma && Fma.IsSupported;

            try
            {
                fixed (double* pSrc = flatMatrix)
                fixed (double* A = workingArray)
                fixed (double* V = outEigenvectors)
                {
                    Unsafe.CopyBlock(A, pSrc, (uint)(sizeof(double) * matrixSize));

                    // Initialize Eigenvectors matrix to Identity (flat column-major)
                    outEigenvectors.Clear();
                    for (int i = 0; i < n; i++) V[i * n + i] = 1.0;

                    for (int sweep = 0; sweep < maxSweeps; sweep++)
                    {
                        double offDiagNormSq = 0.0;

                        // Sweep all unique pairs (p, q) over the strict upper triangle
                        for (int p = 0; p < n - 1; p++)
                        {
                            for (int q = p + 1; q < n; q++)
                            {
                                // In column-major layout, row p, col q is at location [q * n + p]
                                double apq = A[q * n + p];
                                offDiagNormSq += apq * apq;

                                if (Math.Abs(apq) < tol) continue;

                                // Compute Algebraic Givens Rotations (No Atan2, Cos, or Sin)
                                int pBase = p * n;
                                int qBase = q * n;
                                double app = A[pBase + p];
                                double aqq = A[qBase + q];
                                double theta = (aqq - app) / (2.0 * apq);

                                double t;
                                if (Math.Abs(theta) < 1e12)
                                {
                                    double signTheta = theta >= 0.0 ? 1.0 : -1.0;
                                    t = signTheta / (Math.Abs(theta) + Math.Sqrt(1.0 + theta * theta));
                                }
                                else
                                {
                                    t = 1.0 / (2.0 * theta); // Asymptotic expansion for giant theta to prevent division limits
                                }

                                double c = 1.0 / Math.Sqrt(1.0 + t * t);
                                double s = t * c;
                                bool canVectorize = Vector256.IsHardwareAccelerated;
                                Vector256<double> vc = canVectorize ? Vector256.Create(c) : default;
                                Vector256<double> vs = canVectorize ? Vector256.Create(s) : default;

                                // Update the tracking basis vector matrix V (Fully Vectorized Column Transformations)
                                int i = 0;
                                if (canVectorize && n >= 4)
                                {
                                    int simdLen = n & ~3;

                                    if (useFma)
                                    {
                                        for (; i < simdLen; i += 4)
                                        {
                                            var vvip = Vector256.Load(V + pBase + i);
                                            var vviq = Vector256.Load(V + qBase + i);

                                            var vcvip = Vector256.Multiply(vc, vvip);
                                            var vcviq = Vector256.Multiply(vc, vviq);
                                            var vNewVip = Fma.MultiplyAddNegated(vs, vviq, vcvip);
                                            var vNewViq = Fma.MultiplyAdd(vs, vvip, vcviq);

                                            Vector256.Store(vNewVip, V + pBase + i);
                                            Vector256.Store(vNewViq, V + qBase + i);
                                        }
                                    }
                                    else
                                    {
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
                                }
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

                                    if (useFma)
                                    {
                                        for (; zone1 < simdLen; zone1 += 4)
                                        {
                                            var vaip = Vector256.Load(A + pBase + zone1);
                                            var vaiq = Vector256.Load(A + qBase + zone1);

                                            var vcAip = Vector256.Multiply(vc, vaip);
                                            var vcAiq = Vector256.Multiply(vc, vaiq);
                                            var vNewAip = Fma.MultiplyAddNegated(vs, vaiq, vcAip);
                                            var vNewAiq = Fma.MultiplyAdd(vs, vaip, vcAiq);

                                            Vector256.Store(vNewAip, A + pBase + zone1);
                                            Vector256.Store(vNewAiq, A + qBase + zone1);
                                        }
                                    }
                                    else
                                    {
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

                        if (Math.Sqrt(offDiagNormSq) < tol) break;
                    }

                    // Extract final parsed eigenvalues from the diagonal trace
                    for (int i = 0; i < n; i++) outEigenvalues[i] = A[i * n + i];
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(workingArray);
            }

            // In-place sorting matching structural memory constraints
            for (int i = 0; i < n - 1; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    if (outEigenvalues[j] > outEigenvalues[i])
                    {
                        double tmpLambda = outEigenvalues[i];
                        outEigenvalues[i] = outEigenvalues[j];
                        outEigenvalues[j] = tmpLambda;

                        Span<double> colVi = outEigenvectors.Slice(i * n, n);
                        Span<double> colVj = outEigenvectors.Slice(j * n, n);
                        for (int row = 0; row < n; row++)
                        {
                            double tmpV = colVi[row];
                            colVi[row] = colVj[row];
                            colVj[row] = tmpV;
                        }
                    }
                }
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
            double tol = 1e-12,
            DenseEigenFastVariant fastVariant = DenseEigenFastVariant.Default)
        {
            ArrayPool<double> pool = ArrayPool<double>.Shared;
            double[] outEigVecs = pool.Rent(n * n);
            double[] outEigVals = new double[n];

            try
            {
                DecomposeSymmetricInPlace(
                    flatColumnMajorMatrix,
                    n,
                    outEigVals,
                    outEigVecs.AsSpan(0, n * n),
                    maxSweeps,
                    tol,
                    fastVariant);

                var sortedVec = new double[n][];
                for (int i = 0; i < n; i++)
                {
                    sortedVec[i] = new double[n];
                    outEigVecs.AsSpan(i * n, n).CopyTo(sortedVec[i]);
                }

                return new EigenResult(outEigVals, sortedVec);
            }
            finally
            {
                pool.Return(outEigVecs);
            }
        }

        public static EigenResult DecomposeSymmetric(
            double[,] matrix,
            int maxSweeps = 256,
            double tol = 1e-12,
            DenseEigenFastVariant fastVariant = DenseEigenFastVariant.Default)
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

                return DecomposeSymmetric(flatMatrix.AsSpan(0, n * n), n, maxSweeps, tol, fastVariant);
            }
            finally
            {
                pool.Return(flatMatrix);
            }
        }
    }
}
