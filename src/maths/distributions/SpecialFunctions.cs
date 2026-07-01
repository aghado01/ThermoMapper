using System;

namespace Maths.Distributions;

/// <summary>
/// Special functions backing the continuous distributions — currently the log-gamma function, the normalizer
/// every Beta/Gamma log-density needs. Lanczos approximation (g = 7), accurate to ~15 significant digits over
/// the positive reals, with the reflection formula for <c>x &lt; 0.5</c>.
/// </summary>
public static class SpecialFunctions
{
    private static readonly double[] LanczosG7 =
    {
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
    };

    /// <summary>Natural log of the gamma function, <c>ln Γ(x)</c>, for <c>x &gt; 0</c>.</summary>
    public static double LogGamma(double x)
    {
        if (x < 0.5)
            return Math.Log(Math.PI / Math.Sin(Math.PI * x)) - LogGamma(1.0 - x);

        x -= 1.0;
        double a = LanczosG7[0];
        double t = x + 7.5;
        for (int i = 1; i < LanczosG7.Length; i++)
            a += LanczosG7[i] / (x + i);

        return 0.5 * Math.Log(2.0 * Math.PI) + (x + 0.5) * Math.Log(t) - t + Math.Log(a);
    }
}
