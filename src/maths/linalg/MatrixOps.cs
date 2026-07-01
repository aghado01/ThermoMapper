using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Maths.LinAlg
{
    public static class MatrixOps
    {
        /// <summary>
        /// Generates a column-major flat block with <paramref name="k"/> orthonormal
        /// vectors in R^<paramref name="n"/>. Layout is [col0 | col1 | ...].
        /// </summary>
        public static double[] RandomOrthonormal(int n, int k, int seed)
        {
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(n, 0, nameof(n));
            ArgumentOutOfRangeException.ThrowIfLessThan(k, 0, nameof(k));
            if (k > n)
                throw new ArgumentOutOfRangeException(nameof(k), "Column count must not exceed row count.");
            if (seed == 0)
                throw new InvalidOperationException("RandomOrthonormal requires a seed for reproducibility");

            var rng = new Random(seed);
            var block = new double[n * k];
            for (int idx = 0; idx < block.Length; idx++)
                block[idx] = 2.0 * rng.NextDouble() - 1.0;

            Orthonormalize(block, n, k);
            return block;
        }

        /// <summary>
        /// Calculates the dot product of two vectors with native JIT specialization.
        /// Safely tiers down from AVX-512 to AVX2 depending on the executing core's capabilities.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        {
            int len = a.Length;
            if (len != b.Length) throw new ArgumentException("Vector dimensions must match.");

            double sum = 0.0;
            int i = 0;

            // 1. Tier 1: AVX-512 Path (8 doubles per register)
            // If the current core supports it, JIT compiles this and erases the rest.
            // If unsupported, JIT removes this entire block before execution.
            if (Vector512.IsHardwareAccelerated && len >= 8)
            {
                var vecSum = Vector512<double>.Zero;
                int simdLen = len & ~7; // Strip down to multiples of 8

                ref double aRef = ref MemoryMarshal.GetReference(a);
                ref double bRef = ref MemoryMarshal.GetReference(b);

                for (; i < simdLen; i += 8)
                {
                    var va = Vector512.LoadUnsafe(ref aRef, (uint)i);
                    var vb = Vector512.LoadUnsafe(ref bRef, (uint)i);
                    vecSum = Vector512.Add(vecSum, Vector512.Multiply(va, vb));
                }
                sum = Vector512.Sum(vecSum);
            }
            // 2. Tier 2: AVX2 Path (4 doubles per register)
            else if (Vector256.IsHardwareAccelerated && len >= 4)
            {
                var vecSum = Vector256<double>.Zero;
                int simdLen = len & ~3; // Strip down to multiples of 4

                ref double aRef = ref MemoryMarshal.GetReference(a);
                ref double bRef = ref MemoryMarshal.GetReference(b);

                for (; i < simdLen; i += 4)
                {
                    var va = Vector256.LoadUnsafe(ref aRef, (uint)i);
                    var vb = Vector256.LoadUnsafe(ref bRef, (uint)i);
                    vecSum = Vector256.Add(vecSum, Vector256.Multiply(va, vb));
                }
                sum = Vector256.Sum(vecSum);
            }

            // 3. Clean up the scalar remainder (or fallback for non-SIMD environments)
            for (; i < len; i++)
            {
                sum += a[i] * b[i];
            }

            return sum;
        }

        /// <summary>
        /// High-performance in-place vector subtraction mapping: qj -= r * qi.
        /// Uses aggressive tiered JIT optimization matching execution hardware capabilities.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void ScaledSubtract(Span<double> qj, ReadOnlySpan<double> qi, double r, int length)
        {
            int i = 0;
            ref double qjRef = ref MemoryMarshal.GetReference(qj);
            ref double qiRef = ref MemoryMarshal.GetReference(qi);

            // AVX-512 Specialization
            if (Vector512.IsHardwareAccelerated && length >= 8)
            {
                var vr = Vector512.Create(r);
                int simdLen = length & ~7;
                for (; i < simdLen; i += 8)
                {
                    var vqj = Vector512.LoadUnsafe(ref qjRef, (uint)i);
                    var vqi = Vector512.LoadUnsafe(ref qiRef, (uint)i);
                    var vres = Vector512.Subtract(vqj, Vector512.Multiply(vr, vqi));
                    vres.StoreUnsafe(ref qjRef, (uint)i);
                }
            }
            // AVX2 Specialization
            else if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                var vr = Vector256.Create(r);
                int simdLen = length & ~3;
                for (; i < simdLen; i += 4)
                {
                    var vqj = Vector256.LoadUnsafe(ref qjRef, (uint)i);
                    var vqi = Vector256.LoadUnsafe(ref qiRef, (uint)i);
                    var vres = Vector256.Subtract(vqj, Vector256.Multiply(vr, vqi));
                    vres.StoreUnsafe(ref qjRef, (uint)i);
                }
            }

            // Fallback remainder loop
            for (; i < length; i++)
            {
                qj[i] -= r * qi[i];
            }
        }

        /// <summary>
        /// In-place modified Gram-Schmidt over a flat column-major block.
        /// Degenerate columns fall back to a projected standard basis vector.
        /// </summary>
        public static void Orthonormalize(double[] block, int n, int k, double tolerance = 1e-12)
        {
            if (block is null) throw new ArgumentNullException(nameof(block));
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(n, 0, nameof(n));
            ArgumentOutOfRangeException.ThrowIfLessThan(k, 0, nameof(k));
            if (block.Length != n * k)
                throw new ArgumentException("Block length must equal n * k.", nameof(block));

            for (int j = 0; j < k; j++)
            {
                Span<double> qj = block.AsSpan(j * n, n);

                for (int pass = 0; pass < 2; pass++)
                {
                    for (int i = 0; i < j; i++)
                    {
                        ReadOnlySpan<double> qi = block.AsSpan(i * n, n);
                        double projection = Dot(qj, qi);
                        if (projection != 0.0)
                            ScaledSubtract(qj, qi, projection, n);
                    }
                }

                double norm = Math.Sqrt(Dot(qj, qj));
                if (norm < tolerance)
                {
                    qj.Clear();
                    qj[j % n] = 1.0;

                    for (int pass = 0; pass < 2; pass++)
                    {
                        for (int i = 0; i < j; i++)
                        {
                            ReadOnlySpan<double> qi = block.AsSpan(i * n, n);
                            double projection = Dot(qj, qi);
                            if (projection != 0.0)
                                ScaledSubtract(qj, qi, projection, n);
                        }
                    }

                    norm = Math.Sqrt(Dot(qj, qj));
                    if (norm < tolerance)
                        throw new InvalidOperationException("Unable to construct a full-rank orthonormal basis.");
                }

                double invNorm = 1.0 / norm;
                for (int row = 0; row < n; row++)
                    qj[row] *= invNorm;
            }
        }

        /// <summary>
        /// Transposes a row-major jagged array (n rows × p columns) to a column-major jagged array
        /// (p columns × n rows). The result has contiguous column spans suitable for SIMD dot products.
        /// </summary>
        /// <param name="data">Row-major data: n rows, each of length p.</param>
        /// <param name="n">Number of rows.</param>
        /// <param name="p">Number of columns.</param>
        /// <returns>Column-major layout: p columns, each a double[n] of contiguous row values.</returns>
        public static double[][] TransposeToColumnMajor(double[][] data, int n, int p)
        {
            if (data is null) throw new ArgumentNullException(nameof(data));
            if (data.Length != n) throw new ArgumentException("Row count mismatch.", nameof(n));
            if (n == 0 || p == 0) return Array.Empty<double[]>();

            var columns = new double[p][];
            for (int j = 0; j < p; j++)
            {
                var col = new double[n];
                for (int i = 0; i < n; i++)
                    col[i] = data[i][j];
                columns[j] = col;
            }
            return columns;
        }

        /// <summary>
        /// Computes the column Gram matrix G = X^T X (p×p symmetric) using SIMD-accelerated dot products.
        /// Expects column-major input for cache-friendly contiguous access.
        /// </summary>
        /// <param name="columns">Column-major data: p columns, each of length n (contiguous).</param>
        /// <param name="n">Number of rows (column length).</param>
        /// <param name="p">Number of columns (dimension of Gram matrix).</param>
        /// <returns>Symmetric Gram matrix G where G[i,j] = Dot(col_i, col_j).</returns>
        public static double[,] ColumnGramMatrix(double[][] columns, int n, int p)
        {
            if (columns is null) throw new ArgumentNullException(nameof(columns));
            if (columns.Length < p) throw new ArgumentException("Column count mismatch.", nameof(p));

            var gram = new double[p, p];
            for (int i = 0; i < p; i++)
            {
                ReadOnlySpan<double> coli = columns[i].AsSpan(0, n);
                for (int j = i; j < p; j++)
                {
                    ReadOnlySpan<double> colj = columns[j].AsSpan(0, n);
                    double dot = Dot(coli, colj);
                    gram[i, j] = gram[j, i] = dot;
                }
            }
            return gram;
        }

        /// <summary>
        /// Convenience: transpose then compute Gram matrix. Use when the transpose itself is not reusable.
        /// For mxPBF-style workloads where G_X, G_Y, and G_Z = G_X + G_Y are needed, prefer explicit
        /// <see cref="TransposeToColumnMajor"/> to share the transpose cost across all three.
        /// </summary>
        public static double[,] ColumnGramMatrixFromRowMajor(double[][] data, int n, int p)
        {
            var columns = TransposeToColumnMajor(data, n, p);
            return ColumnGramMatrix(columns, n, p);
        }
    }
}
