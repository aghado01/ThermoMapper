using System;
using System.Collections.Generic;

namespace Maths.LinAlg;

/// <summary>
/// Functions of a symmetric matrix read from its eigensystem A = VΛVᵀ (spectral-solver brief, B6): f(A) = V·f(Λ)·Vᵀ, and the factor
/// forms V·f(Λ) over chosen eigenpairs, as columns or as rows, so a root F with FFᵀ = A or a row root R with RᵀR = A comes from one
/// decomposition. Which eigenpairs count, and what is refused, are the caller's rules: the helpers take the selection, not a cutoff.
/// </summary>
public static class MatrixFunctions
{
    /// <summary>
    /// f(A) = Σ_q f(λ_q) v_q v_qᵀ, symmetrized to (B + Bᵀ)/2, from <paramref name="eigen"/> (eigenvectors as rows of
    /// <see cref="EigenResult.Eigenvectors"/>). Terms with f(λ_q) = 0 and an eigenvector's zero components are skipped: they would add a
    /// signed zero to a sum that starts at +0, which leaves it unchanged, so the result is the same bits and a diagonal A costs O(n²).
    /// </summary>
    public static double[,] Symmetric(EigenResult eigen, Func<double, double> f)
    {
        ArgumentNullException.ThrowIfNull(eigen);
        ArgumentNullException.ThrowIfNull(f);
        int n = eigen.Eigenvalues.Length;
        var result = new double[n, n];
        var nonzero = new int[n];
        for (int q = 0; q < n; q++)
        {
            double weight = f(eigen.Eigenvalues[q]);
            if (weight == 0.0) continue;
            double[] v = eigen.Eigenvectors[q];
            int count = 0;
            for (int i = 0; i < n; i++)
                if (v[i] != 0.0) nonzero[count++] = i;
            for (int a = 0; a < count; a++)
            {
                int i = nonzero[a];
                for (int b = 0; b < count; b++)
                {
                    int j = nonzero[b];
                    result[i, j] += weight * v[i] * v[j];
                }
            }
        }
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) result[i, j] = result[j, i] = 0.5 * (result[i, j] + result[j, i]);
        return result;
    }

    /// <summary>The factor with columns f(λ_q)·v_q for the eigenpairs <paramref name="kept"/> selects, in its order: n × kept.Count,
    /// so with f = √ and the positive eigenpairs F·Fᵀ = A.</summary>
    public static double[,] FactorColumns(EigenResult eigen, IReadOnlyList<int> kept, Func<double, double> f)
    {
        ArgumentNullException.ThrowIfNull(eigen);
        ArgumentNullException.ThrowIfNull(kept);
        ArgumentNullException.ThrowIfNull(f);
        int n = eigen.Eigenvalues.Length;
        var factor = new double[n, kept.Count];
        for (int c = 0; c < kept.Count; c++)
        {
            double weight = f(eigen.Eigenvalues[kept[c]]);
            double[] v = eigen.Eigenvectors[kept[c]];
            for (int a = 0; a < n; a++) factor[a, c] = weight * v[a];
        }
        return factor;
    }

    /// <summary>The factor with rows f(λ_q)·v_qᵀ for the eigenpairs <paramref name="kept"/> selects, in its order: kept.Count × n, so
    /// with f = √ and the positive eigenpairs Rᵀ·R = A.</summary>
    public static double[,] FactorRows(EigenResult eigen, IReadOnlyList<int> kept, Func<double, double> f)
    {
        ArgumentNullException.ThrowIfNull(eigen);
        ArgumentNullException.ThrowIfNull(kept);
        ArgumentNullException.ThrowIfNull(f);
        int n = eigen.Eigenvalues.Length;
        var factor = new double[kept.Count, n];
        for (int r = 0; r < kept.Count; r++)
        {
            double weight = f(eigen.Eigenvalues[kept[r]]);
            double[] v = eigen.Eigenvectors[kept[r]];
            for (int a = 0; a < n; a++) factor[r, a] = weight * v[a];
        }
        return factor;
    }
}
