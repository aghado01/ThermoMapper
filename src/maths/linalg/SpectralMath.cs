// src/maths/linalg/SpectralMath.cs
#nullable enable
using System;
using System.Collections.Generic;

namespace Maths.LinAlg;

/// <summary>
/// Eigenvalue / eigenvector pair returned by spectral routines.
/// <c>Vector</c> is a freshly allocated copy — callers may mutate freely.
/// </summary>
public readonly record struct EigenPair(double Lambda, double[] Vector);

/// <summary>
/// Graph-free dense spectral primitives. Operates on materialised symmetric
/// matrices, dispatches eigendecomposition via <see cref="DenseEigen"/>, and
/// selects the bottom-K eigenpairs (smallest eigenvalues, sorted ascending).
///
/// Pair this with a graph-aware adapter — e.g.
/// <c>Graphs.Spectral.Spectral.ComputeBottomK</c> — when consuming
/// <c>CsrGraph</c> inputs: that adapter materialises the Laplacian and then
/// delegates the spectral work to the helpers here.
/// </summary>
public static class SpectralMath
{
    /// <summary>
    /// Bottom-K eigenpairs of a symmetric dense matrix (idiomatic
    /// <c>double[,]</c> layout). Returns the K smallest eigenvalues with
    /// their eigenvectors, sorted ascending by eigenvalue.
    /// </summary>
    public static IReadOnlyList<EigenPair> BottomK(
        double[,] matrix, int k, DenseEigenOptions options = default)
    {
        if (k <= 0) return Array.Empty<EigenPair>();
        EigenResult eig = DenseEigen.DecomposeSymmetric(matrix, options: options);
        return SelectBottomK(eig, k);
    }

    /// <summary>
    /// Bottom-K eigenpairs of a symmetric matrix supplied as a flat
    /// column-major span (LAPACK / EigenFast layout). <paramref name="n"/>
    /// is the matrix dimension; the span length must equal <c>n × n</c>.
    /// </summary>
    public static IReadOnlyList<EigenPair> BottomK(
        ReadOnlySpan<double> flatColumnMajorMatrix, int n, int k,
        DenseEigenOptions options = default)
    {
        if (k <= 0) return Array.Empty<EigenPair>();
        EigenResult eig = DenseEigen.DecomposeSymmetric(flatColumnMajorMatrix, n, options: options);
        return SelectBottomK(eig, k);
    }

    /// <summary>
    /// Selects the bottom-K eigenpairs from an existing decomposition.
    /// <see cref="DenseEigen"/> returns eigenvalues sorted descending, so this
    /// takes the tail K and re-sorts ascending. Vectors are copied; callers
    /// may mutate them without affecting the source <see cref="EigenResult"/>.
    /// </summary>
    public static IReadOnlyList<EigenPair> SelectBottomK(EigenResult eig, int k)
    {
        int total = eig.Eigenvalues.Length;
        if (total == 0 || k <= 0) return Array.Empty<EigenPair>();

        int n = eig.Eigenvectors.Length > 0 ? eig.Eigenvectors[0].Length : 0;
        var pairs = new List<EigenPair>(Math.Min(k, total));
        for (int i = total - 1; i >= 0 && pairs.Count < k; i--)
        {
            var vector = new double[n];
            Array.Copy(eig.Eigenvectors[i], vector, n);
            pairs.Add(new EigenPair(eig.Eigenvalues[i], vector));
        }
        pairs.Sort(static (left, right) => left.Lambda.CompareTo(right.Lambda));
        return pairs;
    }
}
