using System;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bars;

// The per-point variance shape every BARS model scales: scale · z(t)ᵀA⁻¹z(t) for each readout row z(t), from a
// factor of A. Unavailable readouts are NaN.
internal static class ReadoutQuadratics
{
    public static double[] Of(ISpdFactor factor, double[,] readoutDesign, double scale)
        => factor.InverseQuadraticForms(readoutDesign, scale);

    public static double[] Unavailable(double[,] readoutDesign)
    {
        var v = new double[readoutDesign.GetLength(0)];
        Array.Fill(v, double.NaN);
        return v;
    }
}
