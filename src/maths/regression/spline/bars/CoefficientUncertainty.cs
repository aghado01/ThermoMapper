namespace Maths.Regression.Spline.Bars;

/// <summary>
/// How a BARS observation model's coefficient draws and conditional variances relate to the exact conditional
/// posterior of the coefficients given the knots. Readouts that are nonlinear in the fitted curve (argmax, peak
/// counts, spans, intervals) inherit this label.
/// </summary>
public enum CoefficientUncertainty
{
    /// <summary>Draws and variances are exact (the Normal model).</summary>
    Exact,

    /// <summary>Draws and variances come from a Laplace approximation at the fitted estimate (Whittle; Poisson
    /// with a short Metropolis–Hastings correction).</summary>
    Approximate,

    /// <summary>No coefficient uncertainty: the draw is the conditional mean, so nonlinear readouts reflect knot
    /// uncertainty only (GibbsLoss, whose posterior is a generalized one).</summary>
    KnotOnly,
}
