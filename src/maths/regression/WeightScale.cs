using System;
using System.Collections.Generic;

namespace Maths.Regression;

/// <summary>
/// Draw weights read over their whole finite range. A reader of linear weights forms ratios of weighted sums: a normalized
/// weight, a weighted mean, a cumulative pick. It reads each weight scaled by 2^−e, the power of two that brings the largest
/// into [1, 2). The scaling is exact wherever the scaled weight is normal, and cancels in the ratio. So ordinary weights give
/// the bits they gave unscaled. A sum of n scaled weights is below 2n, and a weighted sum of values below twice the values'
/// own absolute sum: <c>[1e308, 1e308]</c>, whose sum overflows, and <c>[5e307, 5e307]</c>, whose weighted sums did, read as
/// <c>[1, 1]</c> does. A positive weight below 2^−1074 of the largest scales to 0 and drops out of the sums, where its share
/// was already below rounding.
/// <see cref="PredictiveChecks"/>' rows read log weights by the same scaling.
/// </summary>
internal static class WeightScale
{
    /// <summary>The exponent e = ⌊log₂ of the largest positive weight⌋; 0 when <paramref name="weights"/> is null or holds
    /// none positive and finite.</summary>
    public static int Exponent(IReadOnlyList<double>? weights)
    {
        if (weights is null) return 0;
        double largest = 0.0;
        foreach (double w in weights)
            if (w > largest && double.IsFinite(w)) largest = w;
        return Exponent(largest);
    }

    /// <summary>The exponent of <paramref name="largest"/>; 0 unless it is positive and finite.</summary>
    public static int Exponent(double largest) => largest > 0.0 && double.IsFinite(largest) ? Math.ILogB(largest) : 0;

    /// <summary><paramref name="weight"/>·2^−<paramref name="exponent"/>: exact where the result is normal.</summary>
    public static double Scaled(double weight, int exponent) => Math.ScaleB(weight, -exponent);
}
