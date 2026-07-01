using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The normal (Gaussian) distribution as a sampler+density primitive over the project die — the standard
/// normal draws also back the Marsaglia–Tsang <see cref="Gamma"/> sampler. Static (zero-alloc) so it is safe
/// in MCMC inner loops.
/// </summary>
public static class Normal
{
    /// <summary>One standard-normal draw via the Box–Muller transform (cosine branch).</summary>
    public static double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double u1 = 1.0 - rng.NextDouble();   // (0,1] avoids log(0)
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    /// <summary>One draw from <c>N(mean, sd²)</c>.</summary>
    public static double Sample(Xoshiro256PlusPlus rng, double mean, double sd)
    {
        if (!(sd > 0.0)) throw new ArgumentOutOfRangeException(nameof(sd), "Standard deviation must be positive.");
        return mean + sd * Sample(rng);
    }

    /// <summary>Log density of <c>N(mean, sd²)</c> at <paramref name="x"/>.</summary>
    public static double LogPdf(double x, double mean, double sd)
    {
        if (!(sd > 0.0)) throw new ArgumentOutOfRangeException(nameof(sd), "Standard deviation must be positive.");
        double z = (x - mean) / sd;
        return -0.5 * z * z - Math.Log(sd) - 0.5 * Math.Log(2.0 * Math.PI);
    }
}
