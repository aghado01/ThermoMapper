using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression;

/// <summary>
/// The global conditionals of a normal variance-mean mixture's noise, εᵢ = σ(β·vᵢ + √(c·vᵢ)·ζᵢ), given the latents and a fit f, with
/// rᵢ = yᵢ − fᵢ, at a likelihood power b (E1, GBASS's w and β steps): σ² under a skew, through u = σ⁻¹, and the skew β. Shared by BASS's
/// blocks and BARS's, so the two engines draw the same conditionals that GBASS's oracle checks.
/// </summary>
public static class VarianceMeanConditionals
{
    /// <summary>
    /// u = σ⁻¹ from its conditional under a skew, MHN(bN + K + 2α₀, ½(b·Σᵢ rᵢ²/(c·vᵢ) + Q + 2β₀), b·β·Σᵢ rᵢ/c) (Sun, Kong and Pal 2021):
    /// <paramref name="squares"/> is Σᵢ rᵢ²/(c·vᵢ), <paramref name="sum"/> Σᵢ rᵢ, and K and Q the coefficients inside a prior whose
    /// variance scales with σ² and their energy per unit σ², the coefficients' prior untempered. σ²'s prior is IG(α₀, β₀)
    /// (<paramref name="prior"/>), or Jeffreys at α₀ = β₀ = 0 when null.
    /// </summary>
    public static double SampleInverseScale(int observations, double squares, double sum, double quadratic, int count, InverseGamma? prior,
                                            double skew, double varianceScale, double power, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double shape = prior?.Shape ?? 0.0, scale = prior?.Scale ?? 0.0;
        return ModifiedHalfNormal.Of(power * observations + count + 2.0 * shape, 0.5 * (power * squares + quadratic + 2.0 * scale),
                                     power * skew * sum / varianceScale).Sample(rng);
    }

    /// <summary>
    /// The skew β from its conditional under the prior N(<paramref name="priorMean"/>, <paramref name="priorSd"/>²): rᵢ/σ = β·vᵢ + √(c·vᵢ)·ζᵢ
    /// makes it normal with precision b·Σvᵢ/c + 1/s² and mean (b·Σrᵢ/(σc) + m/s²) over it, <paramref name="sumV"/> = Σvᵢ and
    /// <paramref name="sumR"/> = Σrᵢ.
    /// </summary>
    public static double SampleSkew(double sumV, double sumR, double sigma, double varianceScale, double priorMean, double priorSd,
                                    double power, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double inverse = 1.0 / (priorSd * priorSd), precision = power * sumV / varianceScale + inverse;
        return (power * sumR / (sigma * varianceScale) + priorMean * inverse) / precision + Normal.SampleStandard(rng) / Math.Sqrt(precision);
    }
}
