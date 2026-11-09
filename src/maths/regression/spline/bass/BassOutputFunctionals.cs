using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

internal static class BassOutputFunctionals
{
    internal static void Fill(BassDraw draw, IReadOnlyList<BassMonitoredPoint> points, Span<double> values, Func<int, double[]> column)
    {
        double[,] sigma = draw.Sigma!;
        int k = sigma.GetLength(0);
        var root = new double[k, k];
        double logDet = 0;
        for (int i = 0; i < k; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = sigma[i, j];
                for (int r = 0; r < j; r++) sum -= root[i, r] * root[j, r];
                if (i == j) { root[i, i] = Math.Sqrt(sum); logDet += Math.Log(sum); }
                else root[i, j] = sum / root[j, j];
            }
        values[0] = logDet;
        var fitted = values[1..];
        double[,] coefficients = draw.OutputCoefficients!;
        for (int p = 0; p < fitted.Length; p++) fitted[p] = coefficients[0, points[p].Location!.Value];
        for (int j = 0; j < draw.Structure.TermIds.Length; j++)
        {
            double[] term = column(draw.Structure.TermIds[j]);
            for (int p = 0; p < fitted.Length; p++) fitted[p] += coefficients[j + 1, points[p].Location!.Value] * term[p];
        }
    }
}
