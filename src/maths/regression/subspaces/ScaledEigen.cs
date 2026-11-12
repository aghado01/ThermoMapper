using System;
using System.Buffers;
using System.Linq;
using Maths.LinAlg;

namespace Maths.Regression.Subspaces;

/// <summary>
/// A symmetric eigendecomposition read at unit scale (E8, slice 4). <see cref="DenseEigen"/>'s Jacobi sweeps stop rotating once every
/// off-diagonal entry is below an absolute threshold, so a matrix of small entries (a covariance in units of 10⁻³, the C of a slowly
/// varying function) would keep off-diagonals as large as its smallest eigenvalues. Dividing by the largest entry first makes the
/// threshold relative, here 10⁻¹⁴ of that entry; the eigenvalues are scaled back.
/// </summary>
internal static class ScaledEigen
{
    public static EigenResult Of(double[,] m)
    {
        ArgumentNullException.ThrowIfNull(m);
        double scale = 0.0;
        foreach (double x in m) scale = Math.Max(scale, Math.Abs(x));
        if (!(scale > 0.0) || !double.IsFinite(scale)) return DenseEigen.DecomposeSymmetric(m);
        int p = m.GetLength(0);
        var unit = new double[p, m.GetLength(1)];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < m.GetLength(1); j++) unit[i, j] = m[i, j] / scale;
        EigenResult e = DenseEigen.DecomposeSymmetric(unit, tol: 1e-14);
        DenseEigenConvergence c = e.Convergence;
        return e with
        {
            Eigenvalues = [.. e.Eigenvalues.Select(x => x * scale)],
            Convergence = c with { OffDiagonalNorm = c.OffDiagonalNorm * scale, LargestOffDiagonal = c.LargestOffDiagonal * scale },
        };
    }

    /// <summary>
    /// The same decomposition into caller-owned outputs (spectral-solver brief, B2): eigenvalues descending into
    /// <paramref name="eigenvalues"/>, and eigenvectors as the columns of <paramref name="columns"/> (p × p), the layout E8 reads, or
    /// none when it is null, which skips their rotations. The matrix is divided into a pooled column-major buffer, so the values and
    /// vectors are the same bits as <see cref="Of"/>'s, with no rows allocated between.
    /// </summary>
    public static DenseEigenConvergence Into(double[,] m, Span<double> eigenvalues, double[,]? columns)
    {
        ArgumentNullException.ThrowIfNull(m);
        int p = m.GetLength(0);
        if (m.GetLength(1) != p) throw new ArgumentException("The matrix is square.", nameof(m));
        double scale = 0.0;
        foreach (double x in m) scale = Math.Max(scale, Math.Abs(x));
        if (!(scale > 0.0) || !double.IsFinite(scale)) return DenseEigen.DecomposeSymmetricInto(m, eigenvalues, columns);
        double[] unit = ArrayPool<double>.Shared.Rent(p * p);
        try
        {
            for (int j = 0; j < p; j++)
                for (int i = 0; i < p; i++) unit[j * p + i] = m[i, j] / scale;
            DenseEigenConvergence c = DenseEigen.DecomposeSymmetricInto(unit.AsSpan(0, p * p), p, eigenvalues, columns, tol: 1e-14);
            for (int i = 0; i < p; i++) eigenvalues[i] *= scale;
            return c with { OffDiagonalNorm = c.OffDiagonalNorm * scale, LargestOffDiagonal = c.LargestOffDiagonal * scale };
        }
        finally
        {
            ArrayPool<double>.Shared.Return(unit);
        }
    }
}
