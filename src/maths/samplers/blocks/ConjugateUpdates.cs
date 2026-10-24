using Maths.Distributions;

namespace Maths.Samplers.Blocks;

/// <summary>
/// The conjugate updates of a normal scale. Given n zero-mean normal terms with sum of squares S, an inverse-gamma
/// prior IG(a, b) on their variance has posterior IG(a + n/2, b + S/2), and a gamma prior Gamma(a, b) (b a rate) on
/// their precision has posterior Gamma(a + n/2, b + S/2). Each returns the posterior, so the caller has the draw and
/// the moments: E[σ² | ·] = (b + S/2)/(a + n/2 − 1) is a Rao–Blackwellized variance estimate.
/// </summary>
public static class ConjugateUpdates
{
    /// <summary>The posterior of a normal variance under an inverse-gamma prior.</summary>
    public static InverseGamma NormalVariance(double priorShape, double priorScale, double count, double sumOfSquares)
        => InverseGamma.ShapeScale(priorShape + 0.5 * count, priorScale + 0.5 * sumOfSquares);

    /// <summary>The posterior of a normal precision under a gamma prior with a rate.</summary>
    public static Gamma NormalPrecision(double priorShape, double priorRate, double count, double sumOfSquares)
        => Gamma.ShapeRate(priorShape + 0.5 * count, priorRate + 0.5 * sumOfSquares);
}
