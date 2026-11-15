using System;

namespace Maths.Distributions;

public static partial class SpecialFunctions
{
    // Stirling's series for ln Γ past its leading terms: B₂ₖ/(2k(2k − 1)) for k = 1 … 6, the coefficients of y^−(2k−1).
    private static readonly double[] StirlingCoefficients =
    {
        1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0, -1.0 / 1680.0, 1.0 / 1188.0, -691.0 / 360360.0,
    };

    // At and above this argument the six-term series is within 2·10⁻¹⁸ of ln Γ's remainder.
    private const double StirlingFloor = 16.0;

    /// <summary>
    /// <c>ln Γ(x + s) − ln Γ(x)</c>, the log of the rising factorial (x)ₛ, for x &gt; 0 and s ≥ 0, without subtracting two
    /// <see cref="LogGamma"/> values. Below 16 the recurrence Γ(x + s)/Γ(x) = x/(x + s) · Γ(x + 1 + s)/Γ(x + 1) raises x, each step
    /// a −ln(1 + s/x); from there Stirling's series is differenced term by term, (x − ½)·ln(1 + s/x) + s·(ln(x + s) − 1) plus the
    /// series' difference. So the ratio keeps its digits at x near 10⁵, where each ln Γ is near 10⁶ and their difference would
    /// lose about five of them: the power spherical law's normalizer at large concentration.
    /// </summary>
    public static double LogGammaRatio(double x, double s)
    {
        if (!(x > 0.0) || double.IsPositiveInfinity(x)) throw new ArgumentOutOfRangeException(nameof(x), "Γ(x + s)/Γ(x) needs a finite x > 0.");
        if (!(s >= 0.0) || double.IsPositiveInfinity(s)) throw new ArgumentOutOfRangeException(nameof(s), "Γ(x + s)/Γ(x) needs a finite s ≥ 0.");
        if (s == 0.0) return 0.0;

        double steps = 0.0;
        while (x < StirlingFloor)
        {
            double ratio = s / x;   // infinite only at a subnormal x, where x + s is s
            steps -= double.IsPositiveInfinity(ratio) ? Math.Log(s) - Math.Log(x) : Log1p(ratio);
            x += 1.0;
        }
        double y = x + s;
        return steps + (x - 0.5) * Log1p(s / x) + s * (Math.Log(y) - 1.0) + StirlingRemainder(y) - StirlingRemainder(x);
    }

    // Σₖ B₂ₖ/(2k(2k − 1))·y^−(2k−1), by Horner in 1/y².
    private static double StirlingRemainder(double y)
    {
        double inverse = 1.0 / y, square = inverse * inverse, sum = 0.0;
        for (int k = StirlingCoefficients.Length - 1; k >= 0; k--) sum = sum * square + StirlingCoefficients[k];
        return sum * inverse;
    }
}
