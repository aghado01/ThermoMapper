using System;

namespace LatentRecoveryProbe;

/// <summary>The dilogarithm Li₂(x) = Σ_{k≥1} x^k/k² on [0, 1]: the series below ½, Euler's reflection above.</summary>
public static class Dilogarithm
{
    public static readonly double ZetaTwo = Math.PI * Math.PI / 6.0;

    public static double Li2(double x)
    {
        if (!(x >= 0.0 && x <= 1.0)) throw new ArgumentOutOfRangeException(nameof(x), "Li2 is evaluated on [0, 1] here.");
        if (x == 0.0) return 0.0;
        if (x == 1.0) return ZetaTwo;
        if (x > 0.5) return ZetaTwo - Math.Log(x) * Math.Log(1.0 - x) - Series(1.0 - x);
        return Series(x);
    }

    private static double Series(double x)
    {
        double sum = 0.0, power = x;
        for (int k = 1; k < 2000; k++)
        {
            double term = power / ((double)k * k);
            sum += term;
            if (term < 1e-19) break;
            power *= x;
        }
        return sum;
    }
}
