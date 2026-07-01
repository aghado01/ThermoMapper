using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The gamma distribution as a sampler+density primitive: Marsaglia–Tsang (2000) squeeze sampling for
/// <c>shape ≥ 1</c>, with the Stuart boost <c>G(a) = G(a+1)·U^(1/a)</c> for <c>shape &lt; 1</c>. Shared by the
/// <see cref="Beta"/> sampler (Beta = G_a/(G_a+G_b)) and the scale-mixture robustness weights (inverse-gamma).
/// Shape/scale (not rate) parameterization; static and zero-alloc for inner-loop use.
/// </summary>
public static class Gamma
{
    /// <summary>One draw from <c>Gamma(shape, scale)</c> (mean <c>shape·scale</c>).</summary>
    public static double Sample(Xoshiro256PlusPlus rng, double shape, double scale = 1.0)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(shape > 0.0)) throw new ArgumentOutOfRangeException(nameof(shape), "Shape must be positive.");
        if (!(scale > 0.0)) throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive.");

        if (shape < 1.0)
        {
            double u = 1.0 - rng.NextDouble();   // (0,1]
            return Sample(rng, shape + 1.0, scale) * Math.Pow(u, 1.0 / shape);
        }

        double d = shape - 1.0 / 3.0;
        double c = 1.0 / Math.Sqrt(9.0 * d);
        while (true)
        {
            double x, v;
            do
            {
                x = Normal.Sample(rng);
                v = 1.0 + c * x;
            }
            while (v <= 0.0);

            v = v * v * v;
            double u = 1.0 - rng.NextDouble();
            double x2 = x * x;
            if (u < 1.0 - 0.0331 * x2 * x2)
                return d * v * scale;
            if (Math.Log(u) < 0.5 * x2 + d * (1.0 - v + Math.Log(v)))
                return d * v * scale;
        }
    }

    /// <summary>Log density of <c>Gamma(shape, scale)</c> at <paramref name="x"/>; −∞ for <c>x ≤ 0</c>.</summary>
    public static double LogPdf(double x, double shape, double scale = 1.0)
    {
        if (!(shape > 0.0)) throw new ArgumentOutOfRangeException(nameof(shape));
        if (!(scale > 0.0)) throw new ArgumentOutOfRangeException(nameof(scale));
        if (x <= 0.0) return double.NegativeInfinity;
        return (shape - 1.0) * Math.Log(x) - x / scale
             - shape * Math.Log(scale) - SpecialFunctions.LogGamma(shape);
    }
}
