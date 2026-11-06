using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The generalized hyperbolic as a normal variance-mean mixture (Barndorff-Nielsen 1977): x = m·v + s·√v·ζ with ζ ~ N(0, 1) and
/// v ~ GIG(p, a, b) (<see cref="GeneralizedInverseGaussian"/>). The Normal-Wald, or normal inverse Gaussian, is p = −½, b = 1; the skewed t
/// a = 0; the variance gamma b = 0.
/// <para><b>Density.</b> Closed: with L(p, a, b) the GIG's log normalizer (<see cref="GeneralizedInverseGaussian.LogNormalizer"/>),
/// ln f(x) = m·x/s² − ½ln(2πs²) + L(p − ½, a + m²/s², b + x²/s²) − L(p, a, b); the integral over v (<see cref="NormalVarianceMeanMixture.IntegratedLogPdf"/>)
/// meets it.</para>
/// <para><b>Tails, quantiles and draws</b> are the mixture's (<see cref="NormalVarianceMeanMixture"/>), over the GIG's density in u = ln v,
/// log-concave there. The variance gamma at p ≤ ½ has the mixture's pole at x = 0.</para>
/// </summary>
public sealed class GeneralizedHyperbolic : NormalVarianceMeanMixture
{
    private readonly double _logNormalizer;
    private readonly GeneralizedInverseGaussian _mixing;

    private GeneralizedHyperbolic(GeneralizedInverseGaussian mixing, double drift, double scale)
        : base(drift, scale, ModeInLog(mixing.P, mixing.A, mixing.B))
    {
        _mixing = mixing;
        _logNormalizer = GeneralizedInverseGaussian.LogNormalizer(mixing.P, mixing.A, mixing.B);
    }

    /// <summary>x = <paramref name="drift"/>·v + <paramref name="scale"/>·√v·ζ, v ~ GIG(<paramref name="p"/>, <paramref name="a"/>,
    /// <paramref name="b"/>), a proper GIG; the drift finite, the scale positive and finite.</summary>
    public static GeneralizedHyperbolic Of(double p, double a, double b, double drift, double scale)
    {
        GeneralizedInverseGaussian mixing = GeneralizedInverseGaussian.Of(p, a, b);
        if (!double.IsFinite(drift)) throw new ArgumentOutOfRangeException(nameof(drift), "The drift must be finite.");
        if (!(scale > 0.0 && double.IsFinite(scale))) throw new ArgumentOutOfRangeException(nameof(scale), "The scale must be positive and finite.");
        return new GeneralizedHyperbolic(mixing, drift, scale);
    }

    /// <summary>v's law.</summary>
    public GeneralizedInverseGaussian Mixing => _mixing;

    /// <inheritdoc/>
    public override double LogPdf(double x)
    {
        double s2 = Scale * Scale, m = Drift;
        return m * x / s2 - 0.5 * Math.Log(2.0 * Math.PI * s2)
               + GeneralizedInverseGaussian.LogNormalizer(_mixing.P - 0.5, _mixing.A + m * m / s2, _mixing.B + x * x / s2) - _logNormalizer;
    }

    // The GIG's density in u = ln v: p·u − (a·eᵘ + b·e⁻ᵘ)/2 − L(p, a, b), each term only where its coefficient is not zero.
    private protected override double LogMixing(double u)
        => _mixing.P * u - 0.5 * ((_mixing.A > 0.0 ? _mixing.A * Math.Exp(u) : 0.0) + (_mixing.B > 0.0 ? _mixing.B * Math.Exp(-u) : 0.0)) - _logNormalizer;

    private protected override double SampleMixing(Xoshiro256PlusPlus rng) => _mixing.Sample(rng);

    // b > 0 drives the density to 0 faster than any power; at b = 0, the gamma's v^(p−1).
    private protected override double LeftExponent => _mixing.B > 0.0 ? double.PositiveInfinity : _mixing.P;

    // The GIG's mode in u = ln v, where p = (a·eᵘ − b·e⁻ᵘ)/2.
    private static double ModeInLog(double p, double a, double b)
    {
        if (a == 0.0) return Math.Log(-0.5 * b / p);
        if (b == 0.0) return Math.Log(2.0 * p / a);
        double root = Math.Sqrt(p * p + a * b);
        return Math.Log(p >= 0.0 ? (p + root) / a : b / (root - p));
    }
}
