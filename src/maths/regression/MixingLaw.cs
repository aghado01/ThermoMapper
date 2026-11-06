using System;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Blocks;

namespace Maths.Regression;

/// <summary>
/// The mixing law of a normal variance-mean mixture, per observation (E1): εᵢ | vᵢ ~ N(m·vᵢ, s²·vᵢ) with vᵢ drawn from this law, so the
/// noise is Gaussian given its latent. m is the unit drift √σ²·β and s² the unit variance σ²·c of the engine that holds the law
/// (GBASS's √w·β and w·c). It supplies what a sampler and a reader of the noise need of v: its prior draw, its conditional draw given a
/// residual at a likelihood power b, the residual's density with v integrated out, which a pointwise density reads since each vᵢ
/// belongs to one observation (OQ17), and v's mean and variance. BARS's <c>IWeightResampler</c> is this seam at m = 0, its weights
/// the precisions 1/v. Shared by BASS and BARS.
/// </summary>
public interface IMixingLaw
{
    /// <summary>v from its prior: the conditional at b = 0.</summary>
    double SamplePrior(Xoshiro256PlusPlus rng);

    /// <summary>
    /// v from its conditional given <paramref name="residual"/> r under the likelihood N(r | m·v, s²·v) raised to <paramref name="power"/>
    /// b, with s² = <paramref name="unitVariance"/> and m = <paramref name="unitDrift"/>: ∝ p(v)·v^(−b/2)·exp(−b(r − m·v)²/(2s²v)).
    /// </summary>
    double SampleConditional(double residual, double unitVariance, double unitDrift, double power, Xoshiro256PlusPlus rng);

    /// <summary>ln ∫ N(r | m·v, s²·v) p(v) dv, the noise's log density at <paramref name="residual"/> with v integrated out.</summary>
    double MarginalLogDensity(double residual, double unitVariance, double unitDrift);

    /// <summary>E[v]; +∞ where it is infinite.</summary>
    double Mean { get; }

    /// <summary>Var[v]; +∞ where it is infinite.</summary>
    double Variance { get; }
}

/// <summary>
/// Student t mixing: v ~ IG(ν/2, ν/2), GIG(−ν/2, 0, ν), so at m = 0 the noise is s·t_ν. Its conditional is conjugate,
/// v | r ~ IG((ν + b)/2, (ν + b·r²/s²)/2), drawn as the reciprocal of the precision's Gamma by the same call BARS's
/// <c>StudentTWeights</c> makes, so at b = 1 the two draw the same numbers from the same stream. The marginal is t_ν(0, s), with
/// ln Γ((ν + 1)/2) − ln Γ(ν/2) and ln(1 + r²/(νs²)) taken so neither loses its digits as ν grows. A drift m ≠ 0, a skewed t whose
/// conditional is a GIG with a &gt; 0, waits for the GIG sampler (E1's slice 2) and is refused.
/// </summary>
public sealed class StudentTMixing : IMixingLaw
{
    /// <summary>Student t mixing on ν = <paramref name="degreesOfFreedom"/> (positive and finite).</summary>
    public StudentTMixing(double degreesOfFreedom)
    {
        if (!(degreesOfFreedom > 0.0 && double.IsFinite(degreesOfFreedom)))
            throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "ν must be positive and finite.");
        DegreesOfFreedom = degreesOfFreedom;
    }

    /// <summary>ν.</summary>
    public double DegreesOfFreedom { get; }

    /// <summary>v's prior as GIG(−ν/2, 0, ν).</summary>
    public GigParameters Prior => GigParameters.OfInverseGamma(0.5 * DegreesOfFreedom, 0.5 * DegreesOfFreedom);

    /// <inheritdoc/>
    public double SamplePrior(Xoshiro256PlusPlus rng) => SampleConditional(0.0, 1.0, 0.0, 0.0, rng);

    /// <inheritdoc/>
    public double SampleConditional(double residual, double unitVariance, double unitDrift, double power, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(unitVariance > 0.0 && double.IsFinite(unitVariance))) throw new ArgumentOutOfRangeException(nameof(unitVariance), "s² must be positive and finite.");
        if (!(power >= 0.0 && power <= 1.0)) throw new ArgumentOutOfRangeException(nameof(power), "The power lies in [0, 1].");
        if (unitDrift != 0.0)
            throw new NotSupportedException("A drift makes t mixing's conditional a GIG with a > 0, which waits for the GIG sampler (E1's slice 2).");
        if (!double.IsFinite(residual)) throw new ArgumentOutOfRangeException(nameof(residual), "The residual must be finite.");
        // StudentTWeights' step: the precision 1/v ~ Gamma((ν + b)/2, (ν + b·z²)/2), z² = r²/s².
        double z2 = residual * residual / unitVariance;
        double half = 0.5 * DegreesOfFreedom;
        return 1.0 / ConjugateUpdates.NormalPrecision(half, half, power, power * z2).Sample(rng);
    }

    /// <inheritdoc/>
    public double MarginalLogDensity(double residual, double unitVariance, double unitDrift)
    {
        if (!(unitVariance > 0.0 && double.IsFinite(unitVariance))) throw new ArgumentOutOfRangeException(nameof(unitVariance), "s² must be positive and finite.");
        if (unitDrift != 0.0)
            throw new NotSupportedException("A drift makes t mixing's marginal a skewed t, with Bessel K of real order, which waits for E1's slice 3.");
        double nu = DegreesOfFreedom;
        return LogGammaHalfRatio(0.5 * nu) - 0.5 * Math.Log(Math.PI * nu * unitVariance)
               - 0.5 * (nu + 1.0) * SpecialFunctions.Log1p(residual * residual / (nu * unitVariance));
    }

    /// <summary>ν/(ν − 2) for ν &gt; 2; +∞ otherwise.</summary>
    public double Mean => DegreesOfFreedom > 2.0 ? DegreesOfFreedom / (DegreesOfFreedom - 2.0) : double.PositiveInfinity;

    /// <summary>2ν²/((ν − 2)²(ν − 4)) for ν &gt; 4; +∞ otherwise.</summary>
    public double Variance
    {
        get
        {
            double nu = DegreesOfFreedom;
            return nu > 4.0 ? 2.0 * nu * nu / ((nu - 2.0) * (nu - 2.0) * (nu - 4.0)) : double.PositiveInfinity;
        }
    }

    /// <summary>
    /// ln Γ(x + ½) − ln Γ(x), x &gt; 0: the difference of the two log gammas below 10³, where it loses under 10⁻¹² absolute; above, the
    /// difference of Stirling's series taken term by term, ½ln x + x·(ln(1 + h) − h) + (1/(12(x + ½)) − 1/(12x)) − (1/(360(x + ½)³) −
    /// 1/(360x³)) with h = 1/(2x), whose dropped terms are under 10⁻¹⁸ there, so the ratio keeps its digits as x grows.
    /// </summary>
    internal static double LogGammaHalfRatio(double x)
    {
        if (x < 1e3) return SpecialFunctions.LogGamma(x + 0.5) - SpecialFunctions.LogGamma(x);
        double h = 0.5 / x;
        // ln(1 + h) − h = −h²/2 + h³/3 − h⁴/4 + h⁵/5 − h⁶/6, the rest under 10⁻²³ at h ≤ 5·10⁻⁴.
        double tail = h * h * (-0.5 + h * (1.0 / 3.0 + h * (-0.25 + h * (0.2 - h / 6.0))));
        double y = x + 0.5;
        return 0.5 * Math.Log(x) + x * tail + (1.0 / (12.0 * y) - 1.0 / (12.0 * x)) - (1.0 / (360.0 * y * y * y) - 1.0 / (360.0 * x * x * x));
    }
}
