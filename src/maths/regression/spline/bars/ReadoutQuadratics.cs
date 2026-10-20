using System;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bars;

// The per-point variance shape every BARS model scales: scale · z(t)ᵀA⁻¹z(t) for each readout row z(t), from a band
// factor of A. Unavailable readouts are NaN.
internal static class ReadoutQuadratics
{
    public static double[] Of(BandCholesky factor, double[,] readoutDesign, double scale)
    {
        int g = readoutDesign.GetLength(0), nu = readoutDesign.GetLength(1);
        if (nu != factor.Dimension)
            throw new ArgumentException("Readout design must have the training design's columns.", nameof(readoutDesign));
        var v = new double[g];
        var row = new double[nu];
        for (int t = 0; t < g; t++)
        {
            for (int j = 0; j < nu; j++) row[j] = readoutDesign[t, j];
            double[] sol = factor.Solve(row);
            double q = 0.0;
            for (int j = 0; j < nu; j++) q += row[j] * sol[j];
            v[t] = scale * q;
        }
        return v;
    }

    public static double[] Unavailable(double[,] readoutDesign)
    {
        var v = new double[readoutDesign.GetLength(0)];
        Array.Fill(v, double.NaN);
        return v;
    }
}
