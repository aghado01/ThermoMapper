using System;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>A flattened tensor-product design with its per-dimension basis counts.</summary>
/// <param name="Design">Row i (point (xᵢ, yᵢ)), column <c>j·NuY + k</c> = <c>B_j(xᵢ)·B_k(yᵢ)</c>.</param>
public sealed record TensorDesign(double[,] Design, int NuX, int NuY);

/// <summary>
/// Builds the flattened bivariate tensor-product B-spline design (He, Yang &amp; Kang 2024). Each row is the outer
/// product of the two 1-D design rows, laid out as <c>j·ν_y + k</c> (x outer, y inner). Each row has only
/// <c>(degree_x+1)·(degree_y+1)</c> non-zeros, so the resulting <c>ZᵀZ</c> is banded (block-banded) and consumes
/// the same <see cref="BandedDesign"/>/<see cref="Maths.LinAlg.BandCholesky"/> path as the 1-D case — paired with
/// a <see cref="TensorPenalty"/> it gives a tensor P-spline. Put the smaller dimension as y (inner) to minimize
/// the flattened bandwidth.
/// </summary>
public static class TensorProductDesign
{
    public static TensorDesign Build(
        SplineBasis basisX, KnotConfig knotsX, double[] xs,
        SplineBasis basisY, KnotConfig knotsY, double[] ys)
    {
        ArgumentNullException.ThrowIfNull(basisX);
        ArgumentNullException.ThrowIfNull(basisY);
        ArgumentNullException.ThrowIfNull(xs);
        ArgumentNullException.ThrowIfNull(ys);
        if (xs.Length != ys.Length)
            throw new ArgumentException("Coordinate arrays must have equal length.", nameof(ys));

        double[,] zx = basisX.Design(knotsX, xs);
        double[,] zy = basisY.Design(knotsY, ys);
        int n = xs.Length;
        int nuX = zx.GetLength(1);
        int nuY = zy.GetLength(1);

        var z = new double[n, nuX * nuY];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < nuX; j++)
            {
                double zxij = zx[i, j];
                if (zxij == 0.0) continue;            // skip the inactive x-block (B-spline sparsity)
                int baseIdx = j * nuY;
                for (int k = 0; k < nuY; k++)
                    z[i, baseIdx + k] = zxij * zy[i, k];
            }
        return new TensorDesign(z, nuX, nuY);
    }
}
