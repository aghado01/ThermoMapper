// ============================================================================
// viz-core/LocalTangent.cs
// ============================================================================
// Per-point empirical tangent vectors via local PCA over k-NN neighborhoods.
//
// For each point x_i, the local neighborhood N(i) is supplied as a list of
// indices into the global points array. The algorithm:
//
//   1. Compute the arithmetic mean of the neighborhood (centred scatter).
//   2. Accumulate the D×D flat-Euclidean sample scatter matrix.
//   3. Extract the leading eigenvector via power iteration.
//
// This is the Wing-2 (Diagnostic TDA) flow source for VectorFieldLayer.
// No model assumptions are made — input is raw feature coordinates and a
// pre-built k-NN adjacency list. The caller (VizApi) reuses the same
// adjacency list already computed for EdgeLayer construction.
//
// Output: flat double[] of length N×D. Each D-element slice is a unit
// tangent vector. Points whose neighborhood is smaller than minNeighbors
// receive a zero vector (filtered at render time by length < ε guard).
//
// All allocations are rented from ArrayPool; the hot path is zero-alloc
// per-point beyond the two reused scratch buffers.
// ============================================================================
using System;
using System.Buffers;
using Maths.Geometry;

namespace Viz
{
    public static class LocalTangent
    {
        /// <summary>
        /// Minimum number of neighbors required to produce a non-zero tangent.
        /// Points with fewer neighbors receive a zero vector.
        /// </summary>
        public const int DefaultMinNeighbors = 2;

        /// <summary>
        /// Power-iteration count. 20 iterations is more than sufficient for
        /// convergence on the small, well-conditioned D×D scatter matrices
        /// arising from local k-NN neighborhoods (D typically 3–4).
        /// </summary>
        public const int DefaultPowerIterations = 20;

        /// <summary>
        /// Computes a unit tangent vector for every point via local PCA over
        /// its k-NN neighborhood.
        /// </summary>
        /// <param name="points">
        /// N point coordinates. Each element must have the same length D.
        /// </param>
        /// <param name="adjacency">
        /// Ragged array of length N. <c>adjacency[i]</c> holds the indices of
        /// the k-nearest neighbours of point i (not including i itself).
        /// </param>
        /// <param name="minNeighbors">
        /// Minimum neighbor count needed to compute a tangent. Points below
        /// this threshold receive a zero vector in the output.
        /// </param>
        /// <param name="powerIterations">
        /// Number of power-iteration steps used to extract the leading
        /// eigenvector of the local scatter matrix.
        /// </param>
        /// <returns>
        /// Flat double[] of length N×D. Slice [i*D .. i*D+D) is the unit
        /// tangent at point i, or zero if the neighborhood was too small.
        /// </returns>
        public static double[] Compute(
            double[][] points,
            int[][] adjacency,
            int minNeighbors = DefaultMinNeighbors,
            int powerIterations = DefaultPowerIterations)
        {
            return Compute(
                points,
                adjacency,
                manifold: null,
                minNeighbors: minNeighbors,
                powerIterations: powerIterations);
        }

        /// <summary>
        /// Computes a unit tangent vector for every point via local PCA over
        /// its k-NN neighborhood, optionally using tangent-space scatter from
        /// a supplied Riemannian manifold.
        /// </summary>
        public static double[] Compute(
            double[][] points,
            int[][] adjacency,
            IRiemannianManifold? manifold,
            int minNeighbors = DefaultMinNeighbors,
            int powerIterations = DefaultPowerIterations)
        {
            if (points is null) throw new ArgumentNullException(nameof(points));
            if (adjacency is null) throw new ArgumentNullException(nameof(adjacency));
            if (points.Length != adjacency.Length)
                throw new ArgumentException(
                    "points and adjacency must have the same length.", nameof(adjacency));

            int n = points.Length;
            if (n == 0) return Array.Empty<double>();

            int d = points[0].Length;
            if (manifold is not null && manifold.Dimension != d)
                throw new ArgumentException(
                    "manifold.Dimension must match the point dimension.", nameof(manifold));

            var result = new double[n * d];
            bool useEuclideanPath = manifold is null || manifold is EuclideanVectorManifold;

            // Reuse three rented scratch buffers across all points.
            double[] meanBuf = ArrayPool<double>.Shared.Rent(d);
            double[] scatterBuf = ArrayPool<double>.Shared.Rent(d * d);
            double[] vecBuf = ArrayPool<double>.Shared.Rent(d);
            double[] nextBuf = ArrayPool<double>.Shared.Rent(d);
            double[] tangentBuf = ArrayPool<double>.Shared.Rent(d);
            double[] workBuf = ArrayPool<double>.Shared.Rent(d);

            try
            {
                Span<double> mean = meanBuf.AsSpan(0, d);
                Span<double> scatter = scatterBuf.AsSpan(0, d * d);
                Span<double> vec = vecBuf.AsSpan(0, d);
                Span<double> next = nextBuf.AsSpan(0, d);
                Span<double> tangent = tangentBuf.AsSpan(0, d);
                Span<double> work = workBuf.AsSpan(0, d);

                for (int i = 0; i < n; i++)
                {
                    int[] nbrs = adjacency[i];
                    if (nbrs is null || nbrs.Length < minNeighbors) continue;
                    int k = nbrs.Length;

                    scatter.Clear();

                    if (useEuclideanPath)
                    {
                        mean.Clear();
                        for (int j = 0; j < k; j++)
                        {
                            double[] q = points[nbrs[j]];
                            for (int dim = 0; dim < d; dim++)
                                mean[dim] += q[dim];
                        }

                        double invK = 1.0 / k;
                        for (int dim = 0; dim < d; dim++) mean[dim] *= invK;

                        for (int j = 0; j < k; j++)
                        {
                            double[] q = points[nbrs[j]];
                            for (int r = 0; r < d; r++)
                            {
                                double vr = q[r] - mean[r];
                                for (int c = 0; c < d; c++)
                                    scatter[r * d + c] += vr * (q[c] - mean[c]);
                            }
                        }
                    }
                    else
                    {
                        ComputeKarcherMean(points, nbrs, manifold!, mean, tangent, work);

                        for (int j = 0; j < k; j++)
                        {
                            manifold!.LogMap(mean, points[nbrs[j]], tangent);
                            for (int r = 0; r < d; r++)
                            {
                                double vr = tangent[r];
                                for (int c = 0; c < d; c++)
                                    scatter[r * d + c] += vr * tangent[c];
                            }
                        }
                    }
                    // No normalisation — eigenvector direction is scale-invariant.

                    // ── 3. Power iteration → leading eigenvector ──────────────
                    // Seed with a non-degenerate vector. Using (1,0,…,0) works
                    // unless the matrix is degenerate, in which case tangent is
                    // undefined and the zero output is correct.
                    vec.Clear();
                    vec[0] = 1.0;

                    for (int iter = 0; iter < powerIterations; iter++)
                    {
                        next.Clear();
                        for (int r = 0; r < d; r++)
                            for (int c = 0; c < d; c++)
                                next[r] += scatter[r * d + c] * vec[c];

                        double norm = 0;
                        for (int dim = 0; dim < d; dim++) norm += next[dim] * next[dim];
                        norm = Math.Sqrt(norm);

                        if (norm < 1e-12) goto zeroed; // degenerate neighbourhood

                        double invNorm = 1.0 / norm;
                        for (int dim = 0; dim < d; dim++) vec[dim] = next[dim] * invNorm;
                    }

                    // ── 4. Write unit tangent into result ─────────────────────
                    int offset = i * d;
                    for (int dim = 0; dim < d; dim++) result[offset + dim] = vec[dim];
                    continue;

                zeroed:; // tangent stays zero for this point
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(meanBuf);
                ArrayPool<double>.Shared.Return(scatterBuf);
                ArrayPool<double>.Shared.Return(vecBuf);
                ArrayPool<double>.Shared.Return(nextBuf);
                ArrayPool<double>.Shared.Return(tangentBuf);
                ArrayPool<double>.Shared.Return(workBuf);
            }

            return result;
        }

        private static void ComputeKarcherMean(
            double[][] points,
            int[] neighbors,
            IRiemannianManifold manifold,
            Span<double> mean,
            Span<double> tangent,
            Span<double> work)
        {
            points[neighbors[0]].AsSpan().CopyTo(mean);
            double invCount = 1.0 / neighbors.Length;

            for (int iteration = 0; iteration < 16; iteration++)
            {
                work.Clear();

                for (int j = 0; j < neighbors.Length; j++)
                {
                    manifold.LogMap(mean, points[neighbors[j]], tangent);
                    for (int dim = 0; dim < mean.Length; dim++)
                        work[dim] += tangent[dim];
                }

                for (int dim = 0; dim < mean.Length; dim++)
                    work[dim] *= invCount;

                if (manifold.Norm(mean, work) < 1e-8)
                    return;

                manifold.ExpMap(mean, work, tangent);
                tangent.CopyTo(mean);
            }
        }

        /// <summary>
        /// BFS orientation propagation: flips tangent signs in-place so that
        /// adjacent vectors agree in direction (dot product ≥ 0).
        ///
        /// Each connected component in the CSR adjacency is seeded independently,
        /// so disconnected sub-graphs are all propagated. Zero-tangent points
        /// (degenerate neighbourhoods) are skipped and left as-is.
        ///
        /// After this call, coherence values computed over graph edges reflect
        /// genuine manifold continuity rather than arbitrary sign choices from
        /// power iteration.  For a non-orientable manifold (e.g. Möbius tube),
        /// BFS will encounter a sign contradiction when closing the loop — the
        /// resulting seam in the flow field is the correct topological signal.
        /// </summary>
        /// <param name="tangents">Flat N×D tangent array from <see cref="Compute(double[][], int[][], int, int)"/>. Modified in-place.</param>
        /// <param name="n">Number of points.</param>
        /// <param name="d">Embedding dimension.</param>
        /// <param name="rowPointers">CSR row-pointer array of length N+1.</param>
        /// <param name="targets">CSR column-index (target) array of length E.</param>
        public static void PropagateOrientation(
            double[] tangents, int n, int d,
            int[] rowPointers, int[] targets)
        {
            if (n == 0 || d == 0) return;

            bool[] visited = new bool[n];
            int[] queue = new int[n];

            for (int seed = 0; seed < n; seed++)
            {
                if (visited[seed]) continue;

                // Skip zero-tangent seeds — no orientation to propagate from.
                int seedOff = seed * d;
                bool hasNonZero = false;
                for (int dim = 0; dim < d; dim++)
                    if (tangents[seedOff + dim] != 0.0) { hasNonZero = true; break; }
                if (!hasNonZero) { visited[seed] = true; continue; }

                int head = 0, tail = 0;
                queue[tail++] = seed;
                visited[seed] = true;

                while (head < tail)
                {
                    int i = queue[head++];
                    int iOff = i * d;
                    int rowStart = rowPointers[i];
                    int rowEnd = rowPointers[i + 1];

                    for (int idx = rowStart; idx < rowEnd; idx++)
                    {
                        int j = targets[idx];
                        if (visited[j]) continue;
                        visited[j] = true;

                        // Align t_j to t_i: flip if anti-aligned.
                        int jOff = j * d;
                        double dot = 0.0;
                        for (int dim = 0; dim < d; dim++)
                            dot += tangents[iOff + dim] * tangents[jOff + dim];
                        if (dot < 0.0)
                            for (int dim = 0; dim < d; dim++)
                                tangents[jOff + dim] = -tangents[jOff + dim];

                        queue[tail++] = j;
                    }
                }
            }
        }

        /// <summary>
        /// Per-point flow coherence: mean dot-product of each point's tangent
        /// with its graph neighbours' tangents.
        ///
        /// Call after <see cref="PropagateOrientation"/> so that sign flips are
        /// already resolved.  The result lies in [-1, 1] per point:
        /// <list type="bullet">
        ///   <item>≈ 1 — neighbourhood is well-aligned (intrinsic graph edge)</item>
        ///   <item>≈ 0 — neighbourhood is random (degenerate / ambient shortcut)</item>
        ///   <item>&lt; 0 — anti-aligned (should not occur post-propagation on
        ///     orientable components; a negative seam on a Möbius tube is
        ///     topologically correct)</item>
        /// </list>
        /// Points with no neighbours receive coherence 0.
        /// </summary>
        /// <param name="tangents">Flat N×D tangent array (post-propagation).</param>
        /// <param name="n">Number of points.</param>
        /// <param name="d">Embedding dimension.</param>
        /// <param name="rowPointers">CSR row-pointer array of length N+1.</param>
        /// <param name="targets">CSR column-index (target) array of length E.</param>
        /// <returns>double[] of length N with per-point coherence scores.</returns>
        public static double[] ComputeCoherence(
            double[] tangents, int n, int d,
            int[] rowPointers, int[] targets)
        {
            var coherence = new double[n];

            for (int i = 0; i < n; i++)
            {
                int rowStart = rowPointers[i];
                int rowEnd = rowPointers[i + 1];
                int degree = rowEnd - rowStart;
                if (degree == 0) continue;

                int iOff = i * d;
                double sum = 0.0;
                for (int idx = rowStart; idx < rowEnd; idx++)
                {
                    int jOff = targets[idx] * d;
                    double dot = 0.0;
                    for (int dim = 0; dim < d; dim++)
                        dot += tangents[iOff + dim] * tangents[jOff + dim];
                    sum += dot;
                }
                coherence[i] = sum / degree;
            }

            return coherence;
        }
    }
}
