using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;

namespace Maths.Regression.Tests;

/// <summary>
/// A BARS posterior with two separated modes of equal mass, for judging tempering where regeneration matters (OQ2).
/// Two equal narrow bumps sit at ¼ and ¾, the noise is mirrored so the data are symmetric about ½, and a cap of
/// <see cref="Cap"/> knots fits both bumps only unequally: three knots on one, two on the other. The target is
/// invariant under ξ ↦ 1 − ξ, which swaps the two allocations, so each mode holds exactly half the mass, a ground
/// truth from outside any estimator. The better-fitted bump carries the curve's maximum, so the peak location's side
/// of ½ names the mode.
/// <para>The death barrier, from a mode's best configuration to its best single deletion, grows with n: about 49 nats
/// at n = 200 and 200 at n = 800 under the default Normal model (calibrated 2026-10-27). Untempered chains did not
/// cross in 20,000 steps at either size.</para>
/// </summary>
internal static class BimodalFixture
{
    /// <summary>The knot cap: enough for three knots on one bump and two on the other, not three on both.</summary>
    public const int Cap = 5;

    private const double Width = 300.0;     // exp(−300 (x − c)²): standard deviation 0.041
    private const double Noise = 0.06;

    /// <summary>The count prior: flat on 0..<see cref="Cap"/>.</summary>
    public static IComplexityPrior Prior() => new UniformPrior(Cap);

    /// <summary>n mirrored design points and responses: x and y at i equal 1 − x and y at n − 1 − i exactly.</summary>
    public static (double[] X, double[] Y) Data(int n)
    {
        if (n % 2 != 0) throw new ArgumentException("n must be even.", nameof(n));
        var rng = new Xoshiro256PlusPlus(seed: 23);
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n / 2; i++)
        {
            x[i] = (i + 0.5) / n;
            x[n - 1 - i] = 1.0 - x[i];
            y[i] = Bump(x[i] - 0.25) + Bump(x[i] - 0.75) + Noise * Gaussian(rng);
            y[n - 1 - i] = y[i];
        }
        return (x, y);
    }

    /// <summary>The reflection ξ ↦ 1 − ξ, which maps each mode onto the other.</summary>
    public static KnotConfig Mirror(KnotConfig config)
    {
        double[] knots = config.InteriorKnots;
        var mirrored = new double[knots.Length];
        for (int j = 0; j < knots.Length; j++) mirrored[j] = 1.0 - knots[knots.Length - 1 - j];
        return new KnotConfig(mirrored);
    }

    private static double Bump(double d) => Math.Exp(-Width * d * d);

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
