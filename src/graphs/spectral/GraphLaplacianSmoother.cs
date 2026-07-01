// src/graphs/spectral/GraphLaplacianSmoother.cs
#nullable enable
using System;
using System.Runtime.CompilerServices;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Spectral;

/// <summary>
/// Graph-Laplacian smoothing of an N×3 point layout. Solves
/// <c>(I + λL)X = X₀</c> for all three coordinate axes in a single fused
/// conjugate-gradient loop — one CSR traversal per iteration instead of three
/// — where <c>L = D − W</c> is the combinatorial graph Laplacian of
/// <see cref="CsrGraph"/>.
///
/// Matrix-free: the Laplacian is never materialised; only the action
/// <c>(I + λL)·p</c> is computed in <see cref="MatvecFused"/> using the CSR
/// adjacency. Distinct from <see cref="GraphLaplacian"/> (which DOES
/// materialise) and from <see cref="Spectral"/> (which decomposes a
/// materialised Laplacian).
/// </summary>
public static class GraphLaplacianSmoother
{
    /// <summary>
    /// Smooths <paramref name="coords"/> under the graph topology.
    /// </summary>
    /// <param name="graph">Adjacency structure; edge weights used as W.</param>
    /// <param name="coords">N input points, each a length-3 array [x, y, z].</param>
    /// <param name="lambda">Smoothing strength. Use <see cref="AutoLambda"/> to
    /// scale automatically to the point density of the dataset.</param>
    /// <param name="maxIter">CG iteration cap (default 300).</param>
    /// <param name="tol">CG residual tolerance (default 1e-6).</param>
    /// <returns>Smoothed coordinates; same shape as <paramref name="coords"/>.</returns>
    public static double[][] Smooth(
        CsrGraph   graph,
        double[][] coords,
        double     lambda,
        int        maxIter = 300,
        double     tol     = 1e-6)
    {
        int n = graph.NodeCount;

        ReadOnlySpan<int>    rowPtrs = graph.RowPointers.AsSpan();
        ReadOnlySpan<int>    targets = graph.Targets.AsSpan();
        ReadOnlySpan<double> weights = graph.Weights.AsSpan();

        double[] degree = new double[n];
        for (int i = 0; i < n; i++)
            for (int e = rowPtrs[i]; e < rowPtrs[i + 1]; e++)
                degree[i] += weights[e];

        double[] b0 = new double[n], b1 = new double[n], b2 = new double[n];
        for (int i = 0; i < n; i++)
        {
            b0[i] = coords[i][0];
            b1[i] = coords[i][1];
            b2[i] = coords[i][2];
        }

        (double[] x0, double[] x1, double[] x2) =
            FusedCgSolve(graph, degree, lambda, b0, b1, b2, maxIter, tol);

        double[][] result = new double[n][];
        for (int i = 0; i < n; i++)
            result[i] = new double[] { x0[i], x1[i], x2[i] };
        return result;
    }

    /// <summary>
    /// Returns lambda = 0.1 / medianEdgeLength², normalising smoothing strength
    /// to the point density of the dataset.
    /// </summary>
    public static double AutoLambda(CsrGraph graph, double[][] coords)
    {
        int m = graph.Targets.Length;
        double[] lengths = new double[m];
        int idx = 0;

        ReadOnlySpan<int> rowPtrs = graph.RowPointers.AsSpan();
        ReadOnlySpan<int> targets = graph.Targets.AsSpan();

        for (int i = 0; i < graph.NodeCount; i++)
        {
            double[] pi = coords[i];
            for (int e = rowPtrs[i]; e < rowPtrs[i + 1]; e++)
            {
                double[] pj = coords[targets[e]];
                double dx = pi[0] - pj[0], dy = pi[1] - pj[1], dz = pi[2] - pj[2];
                lengths[idx++] = Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }
        }

        Array.Sort(lengths, 0, idx);
        double median = lengths[idx / 2];
        return 0.1 / (median * median);
    }

    // ── Fused CG ─────────────────────────────────────────────────────────────
    // Three independent scalar CG solves sharing one CSR matvec per iteration.
    // Each coordinate maintains its own alpha/beta so the solvers are truly
    // independent; convergence is declared when all three residuals are below tol.

    private static (double[] x0, double[] x1, double[] x2) FusedCgSolve(
        CsrGraph graph, double[] degree, double lambda,
        double[] b0, double[] b1, double[] b2,
        int maxIter, double tol)
    {
        int n = graph.NodeCount;

        double[] x0 = new double[n], x1 = new double[n], x2 = new double[n];
        double[] r0 = new double[n], r1 = new double[n], r2 = new double[n];
        double[] p0 = new double[n], p1 = new double[n], p2 = new double[n];
        double[] ap0 = new double[n], ap1 = new double[n], ap2 = new double[n];

        Array.Copy(b0, r0, n); Array.Copy(b1, r1, n); Array.Copy(b2, r2, n);
        Array.Copy(r0, p0, n); Array.Copy(r1, p1, n); Array.Copy(r2, p2, n);

        double rs0 = Dot(r0, r0), rs1 = Dot(r1, r1), rs2 = Dot(r2, r2);

        for (int iter = 0; iter < maxIter; iter++)
        {
            MatvecFused(graph, degree, lambda, p0, p1, p2, ap0, ap1, ap2);

            double alpha0 = rs0 / Dot(p0, ap0);
            double alpha1 = rs1 / Dot(p1, ap1);
            double alpha2 = rs2 / Dot(p2, ap2);

            for (int i = 0; i < n; i++)
            {
                x0[i] += alpha0 * p0[i]; r0[i] -= alpha0 * ap0[i];
                x1[i] += alpha1 * p1[i]; r1[i] -= alpha1 * ap1[i];
                x2[i] += alpha2 * p2[i]; r2[i] -= alpha2 * ap2[i];
            }

            double rn0 = Dot(r0, r0), rn1 = Dot(r1, r1), rn2 = Dot(r2, r2);
            if (Math.Sqrt(Math.Max(rn0, Math.Max(rn1, rn2))) < tol) break;

            double beta0 = rn0 / rs0, beta1 = rn1 / rs1, beta2 = rn2 / rs2;
            for (int i = 0; i < n; i++)
            {
                p0[i] = r0[i] + beta0 * p0[i];
                p1[i] = r1[i] + beta1 * p1[i];
                p2[i] = r2[i] + beta2 * p2[i];
            }

            rs0 = rn0; rs1 = rn1; rs2 = rn2;
        }

        return (x0, x1, x2);
    }

    // Single CSR traversal writing into all three coordinate output vectors.
    // (I + λL)p[i] = (1 + λ·deg[i])·p[i] − λ·Σ_e w_e·p[target_e]
    private static void MatvecFused(
        CsrGraph graph, double[] degree, double lambda,
        double[] p0, double[] p1, double[] p2,
        double[] ap0, double[] ap1, double[] ap2)
    {
        ReadOnlySpan<int>    rowPtrs = graph.RowPointers.AsSpan();
        ReadOnlySpan<int>    targets = graph.Targets.AsSpan();
        ReadOnlySpan<double> weights = graph.Weights.AsSpan();

        int n = graph.NodeCount;
        for (int i = 0; i < n; i++)
        {
            double scale = 1.0 + lambda * degree[i];
            ap0[i] = scale * p0[i];
            ap1[i] = scale * p1[i];
            ap2[i] = scale * p2[i];

            int start = rowPtrs[i], end = rowPtrs[i + 1];
            for (int e = start; e < end; e++)
            {
                double lw = lambda * weights[e];
                int    t  = targets[e];
                ap0[i] -= lw * p0[t];
                ap1[i] -= lw * p1[t];
                ap2[i] -= lw * p2[t];
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Dot(double[] a, double[] b)
        => MatrixOps.Dot(a.AsSpan(), b.AsSpan());
}
