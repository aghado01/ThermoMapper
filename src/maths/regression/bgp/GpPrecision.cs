using System;
using Maths.LinAlg;

namespace Maths.Regression.Bgp;

/// <summary>Reduces a covariance factor to selected precision blocks for GP leave-out scoring.</summary>
internal sealed class GpPrecision
{
    private readonly int _n;
    private readonly double[] _l;

    public GpPrecision(DenseCholesky factor)
    {
        _n = factor.Dimension;
        _l = new double[_n * (_n + 1) / 2];
        for (int i = 0; i < _n; i++)
            for (int j = 0; j <= i; j++)
                _l[Row(i) + j] = factor.Factor(i, j);
    }

    // Only the nonzero tails of L⁻¹'s requested columns are needed. A partition requests each column once;
    // they can be released after its unit, and no off-block entry of R⁻¹ = L⁻ᵀL⁻¹ is computed.
    public double[,] Block(int[] members)
    {
        int m = members.Length;
        var columns = new double[m][];
        for (int a = 0; a < m; a++)
        {
            int j = members[a];
            var column = new double[_n - j];
            column[0] = 1.0 / _l[Row(j) + j];
            for (int i = j + 1; i < _n; i++)
            {
                int row = Row(i);
                double sum = 0.0;
                for (int k = j; k < i; k++) sum -= _l[row + k] * column[k - j];
                column[i - j] = sum / _l[row + i];
            }
            columns[a] = column;
        }
        var block = new double[m, m];
        for (int a = 0; a < m; a++)
            for (int b = 0; b <= a; b++)
            {
                int i = members[a], j = members[b];
                double sum = 0.0;
                for (int k = Math.Max(i, j); k < _n; k++)
                    sum += columns[a][k - i] * columns[b][k - j];
                block[a, b] = block[b, a] = sum;
            }
        return block;
    }

    private static int Row(int i) => i * (i + 1) / 2;
}
