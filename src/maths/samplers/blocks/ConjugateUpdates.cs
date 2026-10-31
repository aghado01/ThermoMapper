using Maths.Distributions;

namespace Maths.Samplers.Blocks;

/// <summary>
/// The conjugate updates of a normal scale. Given n zero-mean normal terms with sum of squares S, an inverse-gamma
/// prior IG(a, b) on their variance has posterior IG(a + n/2, b + S/2), and a gamma prior Gamma(a, b) (b a rate) on
/// their precision has posterior Gamma(a + n/2, b + S/2). Each returns the posterior, so the caller has the draw and
/// the moments: E[σ² | ·] = (b + S/2)/(a + n/2 − 1) is a Rao–Blackwellized variance estimate.
/// <para><b>Limits of one update (N-X3).</b> Both are limits of the update of a normal scale on the generalized inverse Gaussian
/// family (<see cref="GigParameters"/>, GIG(p, a, b) ∝ x^(p−1)·exp(−(a·x + b/x)/2)). A GIG(p, a, b) prior on a variance takes
/// the posterior GIG(p − n/2, a, b + S), and on a precision GIG(p + n/2, a + S, b). The inverse-gamma update is the first at
/// a = 0 (IG(α, θ) is GIG(−α, 0, 2θ)) and the gamma update the second at b = 0 (Gamma(α, β) is GIG(α, 2β, 0)), with the
/// shapes, scales and rates they have always had, bit for bit. Generalized-hyperbolic noise's latent-scale block updates by the
/// full form, so it generalizes these two rather than sitting beside them. No GIG sampler lives here.</para>
/// </summary>
public static class ConjugateUpdates
{
    /// <summary>The posterior of a normal variance under an inverse-gamma prior.</summary>
    public static InverseGamma NormalVariance(double priorShape, double priorScale, double count, double sumOfSquares)
        => InverseGamma.ShapeScale(priorShape + 0.5 * count, priorScale + 0.5 * sumOfSquares);

    /// <summary>The posterior of a normal precision under a gamma prior with a rate.</summary>
    public static Gamma NormalPrecision(double priorShape, double priorRate, double count, double sumOfSquares)
        => Gamma.ShapeRate(priorShape + 0.5 * count, priorRate + 0.5 * sumOfSquares);

    /// <summary>The posterior of a normal variance under a GIG prior on it: GIG(p − n/2, a, b + S). Its a = 0 limit is
    /// <see cref="NormalVariance(double, double, double, double)"/>.</summary>
    public static GigParameters NormalVariance(GigParameters prior, double count, double sumOfSquares)
        => new(prior.P - 0.5 * count, prior.A, prior.B + sumOfSquares);

    /// <summary>The posterior of a normal precision under a GIG prior on it: GIG(p + n/2, a + S, b). Its b = 0 limit is
    /// <see cref="NormalPrecision(double, double, double, double)"/>.</summary>
    public static GigParameters NormalPrecision(GigParameters prior, double count, double sumOfSquares)
        => new(prior.P + 0.5 * count, prior.A + sumOfSquares, prior.B);
}
