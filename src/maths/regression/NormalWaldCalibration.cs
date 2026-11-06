using System;
using Maths.Distributions;

namespace Maths.Regression;

/// <summary>
/// The normal inverse Gaussian NIG(α, β, δ, μ) by its moments: X = μ + β·V + √V·Z with V ~ IG(mean δ/γ, shape δ²) and γ = √(α² − β²).
/// In the BASS core's terms (<c>BassNoise.NormalWald</c>, ε = σ(β′v + √v·ζ), v ~ GIG(−½, γ′², 1)) it is σ = δ, β′ = β·δ, γ′ = γ·δ.
/// </summary>
/// <param name="Mu">μ, the location.</param>
/// <param name="Delta">δ, the scale.</param>
/// <param name="Beta">β, the skew.</param>
/// <param name="Alpha">α, the tail, α &gt; |β|.</param>
public readonly record struct NormalInverseGaussian(double Mu, double Delta, double Beta, double Alpha)
{
    /// <summary>γ = √(α² − β²).</summary>
    public double Gamma => Math.Sqrt((Alpha - Beta) * (Alpha + Beta));

    /// <summary>σ² of the core's form, δ².</summary>
    public double Sigma2 => Delta * Delta;

    /// <summary>β′ of the core's form, β·δ.</summary>
    public double CoreSkew => Beta * Delta;

    /// <summary>γ′ of the core's form, γ·δ.</summary>
    public double CoreGamma => Gamma * Delta;

    /// <summary>The steepness ξ = (1 + δγ)^(−½) of the NIG's shape triangle.</summary>
    public double Steepness => 1.0 / Math.Sqrt(1.0 + Delta * Gamma);
}

/// <summary>
/// The Normal-Wald's calibrations, GBASS's <c>nw_est_mom</c> and <c>nw_gamma_prior</c> corrected (E1's slice 0 confirmed both defects).
/// <list type="bullet">
/// <item><b>Moments.</b> GBASS fits the four moments by <c>optim</c> with its kurtosis 3 + 3(1 + 4β²)/(α²δγ); the NIG's is
/// 3 + 3(1 + 4β²/α²)/(δγ). With it the fit is closed: with skewness S and excess kurtosis K, ρ = β/α has ρ² = S²/(3K − 4S²), δγ =
/// 3(1 + 4ρ²)/K, γ² = δγ/((1 − ρ²)·Var), then δ, α = γ/√(1 − ρ²), β = ρα with S's sign and μ = mean − δβ/γ. The moments are an NIG's only
/// where 3K &gt; 4S²; elsewhere it is refused.</item>
/// <item><b>The prior on γ.</b> GBASS's normal on γ reads P(ξ &lt; q₁) and P(ξ &gt; q₂), ξ = (1 + γ)^(−½), with a stray −1 in the second
/// (written 1 − q₂⁻² − 1 for 1 − q₂⁻²), and its default targets (q₁, q₂, p₁, p₂) = (0.1, 0.9, 0.5, 0.05), half the mass above γ = 99 and 5%
/// below 0.235, are out of any normal's reach. Here the prior is the normal's part on γ &gt; 0, the one the sweep draws γ from, and the
/// targets are solved exactly: for each sd, the mean that puts p₁ above g₁ = q₁⁻² − 1, then the sd that puts p₂ below g₂ = q₂⁻² − 1.
/// At P(γ &gt; g₁) = ½ the most such a prior can put below g₂ is 1 − 2^(−g₂/g₁), the exponential limit as the mean goes to −∞; GBASS's
/// defaults ask 0.05 where that is 0.0016. The defaults here, (0.5, 0.9, 0.5, 0.03), half the mass at ξ &lt; ½ (γ &gt; 3) and 3% at
/// ξ &gt; 0.9 (γ &lt; 0.235, the heaviest tails), sit inside the 0.053 that bound allows. A prior on ln γ reaches any targets in closed
/// form but makes γ's conditional non-conjugate, a slice step in place of the truncated normal's exact draw (docs/Maths).</item>
/// </list>
/// </summary>
public static class NormalWaldCalibration
{
    /// <summary>The NIG with the sample's <paramref name="mean"/>, <paramref name="variance"/>, <paramref name="skewness"/> and
    /// <paramref name="kurtosis"/> (the fourth standardized moment, 3 for a normal), in closed form.</summary>
    public static NormalInverseGaussian Moments(double mean, double variance, double skewness, double kurtosis)
    {
        if (!double.IsFinite(mean)) throw new ArgumentOutOfRangeException(nameof(mean), "The mean must be finite.");
        if (!(variance > 0.0 && double.IsFinite(variance))) throw new ArgumentOutOfRangeException(nameof(variance), "The variance must be positive and finite.");
        double s2 = skewness * skewness, excess = kurtosis - 3.0;
        if (!(3.0 * excess > 4.0 * s2) || !double.IsFinite(excess))
            throw new ArgumentException($"Skewness {skewness} and kurtosis {kurtosis} are no NIG's: an NIG has 3(κ − 3) > 4S².", nameof(kurtosis));
        double rho2 = s2 / (3.0 * excess - 4.0 * s2), deltaGamma = 3.0 * (1.0 + 4.0 * rho2) / excess;
        double gamma = Math.Sqrt(deltaGamma / ((1.0 - rho2) * variance)), delta = deltaGamma / gamma;
        double alpha = gamma / Math.Sqrt(1.0 - rho2), beta = Math.CopySign(Math.Sqrt(rho2) * alpha, skewness);
        return new NormalInverseGaussian(mean - delta * beta / gamma, delta, beta, alpha);
    }

    /// <summary>The NIG of <paramref name="data"/>'s sample moments (<see cref="Moments(double, double, double, double)"/>).</summary>
    public static NormalInverseGaussian Moments(ReadOnlySpan<double> data)
    {
        int n = data.Length;
        if (n < 4) throw new ArgumentException("Four moments need four values.", nameof(data));
        double mean = 0.0;
        foreach (double x in data) mean += x;
        mean /= n;
        double m2 = 0.0, m3 = 0.0, m4 = 0.0;
        foreach (double x in data)
        {
            double d = x - mean, d2 = d * d;
            m2 += d2;
            m3 += d2 * d;
            m4 += d2 * d2;
        }
        double variance = m2 / (n - 1), sd = Math.Sqrt(variance);
        return Moments(mean, variance, m3 / n / (sd * sd * sd), m4 / n / (variance * variance));
    }

    /// <summary>
    /// The normal N(m, s²) on γ &gt; 0 that puts <paramref name="p1"/> of its mass at steepness ξ &lt; <paramref name="q1"/> (γ &gt; q₁⁻² − 1)
    /// and <paramref name="p2"/> at ξ &gt; <paramref name="q2"/> (γ &lt; q₂⁻² − 1), 0 &lt; q₁ &lt; q₂ &lt; 1, as (m, s). Targets no such
    /// normal reaches are refused, with the most it reaches.
    /// </summary>
    public static (double Mean, double Sd) GammaPrior(double q1 = 0.5, double q2 = 0.9, double p1 = 0.5, double p2 = 0.03)
    {
        if (!(q1 > 0.0 && q1 < q2 && q2 < 1.0)) throw new ArgumentOutOfRangeException(nameof(q1), "The steepness levels lie in 0 < q₁ < q₂ < 1.");
        if (!(p1 > 0.0 && p2 > 0.0 && p1 + p2 < 1.0)) throw new ArgumentOutOfRangeException(nameof(p1), "The masses are positive and leave room between them.");
        double g1 = 1.0 / (q1 * q1) - 1.0, g2 = 1.0 / (q2 * q2) - 1.0;
        double Below(double s) => LogBelow(MeanAt(s), s, g2) - Math.Log(p2);
        // Scan ln s for a change of sign, then bisect it.
        double previous = double.NaN, best = double.NegativeInfinity;
        for (int k = 0; k <= 400; k++)
        {
            double s = g1 * Math.Exp(-12.0 + 24.0 * k / 400.0), value = Below(s);
            if (double.IsFinite(value)) best = Math.Max(best, value);
            if (k > 0 && double.IsFinite(previous) && double.IsFinite(value) && Math.Sign(previous) != Math.Sign(value))
            {
                double lo = g1 * Math.Exp(-12.0 + 24.0 * (k - 1) / 400.0), hi = s, low = previous;
                for (int i = 0; i < 200; i++)
                {
                    double middle = Math.Sqrt(lo * hi), at = Below(middle);
                    if (Math.Sign(at) == Math.Sign(low)) { lo = middle; low = at; }
                    else hi = middle;
                }
                double sd = Math.Sqrt(lo * hi);
                return (MeanAt(sd), sd);
            }
            previous = value;
        }
        throw new ArgumentException($"No normal on γ > 0 puts {p1} above γ = {g1:G4} and {p2} below γ = {g2:G4}: the most it puts below is about " +
                                    $"{p2 * Math.Exp(best):G3}.", nameof(p2));

        // The mean that puts p1 above g1 at sd s: that mass rises with the mean.
        double MeanAt(double s)
        {
            double lo = -1e8 * s, hi = g1 + 60.0 * s;
            for (int i = 0; i < 300; i++)
            {
                double middle = 0.5 * (lo + hi);
                if (LogAbove(middle, s, g1) < Math.Log(p1)) lo = middle;
                else hi = middle;
            }
            return 0.5 * (lo + hi);
        }
    }

    // ln P(γ > g) and ln P(γ < g) under N(m, s²) on γ > 0.
    private static double LogAbove(double m, double s, double g)
        => Normal.StandardLogCdf((m - g) / s) - Normal.StandardLogCdf(m / s);

    // P(0 < γ < g) as a difference of the two small normal tails, whichever side of the mean holds them.
    private static double LogBelow(double m, double s, double g)
    {
        double mass = m > 0.0
            ? Math.Exp(Normal.StandardLogCdf((g - m) / s)) - Math.Exp(Normal.StandardLogCdf(-m / s))
            : Math.Exp(Normal.StandardLogCdf(m / s)) - Math.Exp(Normal.StandardLogCdf((m - g) / s));
        return Math.Log(mass) - Normal.StandardLogCdf(m / s);
    }
}
