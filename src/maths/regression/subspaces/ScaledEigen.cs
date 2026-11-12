using System;
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
}
