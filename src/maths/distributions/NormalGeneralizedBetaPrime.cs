using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The normal variance-mean mixture over the generalized beta prime: x = m·v + s·√v·ζ with ζ ~ N(0, 1) and v ~ GBP(p, a, b)
/// (<see cref="GeneralizedBetaPrime"/>), as the normal inverse Gaussian is named for its mixing law. At p = 1, a = b = ½ and m = 0 it is the
/// horseshoe (Carvalho, Polson and Scott 2010), whose density (2π³)^(−½)·e^(z²/2)·E₁(z²/2)/s at z = x/s is the one closed case; the
/// Strawderman–Berger is p = ½, a = 1, b = ½. Its density, tails, quantiles and draws are the mixture's (<see cref="NormalVarianceMeanMixture"/>),
/// the density by the integral over v, over the GBP's density in u = ln v, log-concave there. With pa ≤ ½, the horseshoe's and the
/// Strawderman–Berger's case, the density has the mixture's pole at x = 0; tails fall as |x|^(−2pb − 1) without a drift.
/// </summary>
public sealed class NormalGeneralizedBetaPrime : NormalVarianceMeanMixture
{
    private readonly GeneralizedBetaPrime _mixing;

    private NormalGeneralizedBetaPrime(GeneralizedBetaPrime mixing, double drift, double scale) : base(drift, scale, mixing.ModeInLog)
        => _mixing = mixing;

    /// <summary>x = <paramref name="drift"/>·v + <paramref name="scale"/>·√v·ζ, v ~ GBP(<paramref name="p"/>, <paramref name="a"/>,
    /// <paramref name="b"/>); the drift finite, the scale positive and finite.</summary>
    public static NormalGeneralizedBetaPrime Of(double p, double a, double b, double drift, double scale)
    {
        GeneralizedBetaPrime mixing = GeneralizedBetaPrime.Of(p, a, b);
        if (!double.IsFinite(drift)) throw new ArgumentOutOfRangeException(nameof(drift), "The drift must be finite.");
        if (!(scale > 0.0 && double.IsFinite(scale))) throw new ArgumentOutOfRangeException(nameof(scale), "The scale must be positive and finite.");
        return new NormalGeneralizedBetaPrime(mixing, drift, scale);
    }

    /// <summary>v's law.</summary>
    public GeneralizedBetaPrime Mixing => _mixing;

    private protected override double LogMixing(double u) => _mixing.LogDensityInLog(u);

    private protected override double SampleMixing(Xoshiro256PlusPlus rng) => _mixing.Sample(rng);

    private protected override double LeftExponent => _mixing.P * _mixing.A;
}
