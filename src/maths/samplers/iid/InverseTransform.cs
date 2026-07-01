using System;
using Maths.Rng;

namespace Maths.Samplers.Iid;

/// <summary>
/// Inverse-transform Monte Carlo sampling: draw <c>X = Q(U)</c> with
/// <c>U ~ Uniform[0,1)</c> and <c>Q</c> the quantile (inverse CDF) of the target
/// distribution. The generic engine behind any "draw M samples from a
/// distribution you can invert" routine.
/// </summary>
/// <remarks>
/// SPC's Wang sampler inherited the Boltzmann energy draw <c>H = −T·ln(1−r)</c>
/// — exponential inverse-transform sampling — from the Swendsen–Wang framework.
/// The closed-form correlation makes that draw unnecessary <i>there</i>, which
/// is exactly why the machinery belongs here as a general tool rather than bolted
/// onto the sampler: a pocket capability, reusable wherever a known quantile and
/// a stream of uniforms meet.
/// </remarks>
public static class InverseTransform
{
    /// <summary>
    /// Draw <paramref name="count"/> samples from the distribution whose
    /// quantile (inverse CDF) is <paramref name="quantile"/>. The delegate
    /// receives <c>U ∈ [0,1)</c> from <paramref name="rng"/>.
    /// </summary>
    public static double[] Sample(Xoshiro256PlusPlus rng, Func<double, double> quantile, int count)
    {
        ArgumentNullException.ThrowIfNull(rng);
        ArgumentNullException.ThrowIfNull(quantile);
        if (count < 0) throw new ArgumentOutOfRangeException(nameof(count));

        var draws = new double[count];
        for (int i = 0; i < count; i++)
            draws[i] = quantile(rng.NextDouble());
        return draws;
    }

    /// <summary>
    /// Draw <paramref name="count"/> samples from an exponential distribution of
    /// mean <paramref name="scale"/> (rate <c>1/scale</c>) via
    /// <c>Q(u) = −scale·ln(1−u)</c> — the Boltzmann energy budget
    /// <c>H = −T·ln(1−r)</c> at <c>scale = T</c>.
    /// </summary>
    public static double[] Exponential(Xoshiro256PlusPlus rng, double scale, int count)
    {
        if (scale <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(scale), "scale must be positive.");
        return Sample(rng, u => -scale * Math.Log(1.0 - u), count);
    }
}
