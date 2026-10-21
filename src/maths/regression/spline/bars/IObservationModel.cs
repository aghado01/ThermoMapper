using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The likelihood component of the free-knot model: given a design matrix it returns the log marginal
/// likelihood with the spline coefficients and scale integrated out — the only quantity the marginalized
/// (k,ξ) chain consumes. The single seam behind which the Normal model (exact), Poisson/Whittle (Laplace-BIC)
/// and the two robustness routes sit. <see cref="EffectiveDimension"/> is model-reported so the complexity
/// prior can charge the right penalty (ν for Normal, a robust trace for bounded-loss) without assuming it
/// equals the spline dimension.
/// </summary>
/// <remarks>
/// The scale-mixture robustness route adds an outer weight-resampling step (a later member of this seam);
/// the Normal model below needs only these two operations.
/// </remarks>
public interface IObservationModel
{
    /// <summary>
    /// Log <c>p(y | design)</c> with coefficients and scale marginalized out — returned up to an additive
    /// constant that depends only on the sample size (it cancels in every marginal-likelihood ratio).
    /// <paramref name="weights"/> null means unit weights.
    /// </summary>
    double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights);

    /// <summary>
    /// The same marginal as a function of the likelihood's inverse temperature b: the coefficients and scale
    /// integrated against <c>p(y | β, σ)^b</c> (the likelihood-tempered path, Q5). The design's work is done here,
    /// once, so each b is then a query. The default serves b = 1 only, from <see cref="LogMarginalLikelihood"/>; a
    /// model with a closed form at every b overrides it.
    /// </summary>
    ITemperedMarginal TemperedMarginal(double[,] design, double[] y, double[]? weights)
        => new UntemperedMarginal(LogMarginalLikelihood(design, y, weights));

    /// <summary>True when <see cref="TemperedMarginal"/> serves every b in (0, 1] at unit weights, so the model can
    /// be tempered on the likelihood path. The default serves b = 1 only.</summary>
    bool CanTemperLikelihood => false;

    /// <summary>Effective number of parameters for the complexity penalty (ν for the Normal model).</summary>
    int EffectiveDimension(double[,] design, double[]? weights);

    /// <summary>
    /// Posterior mean coefficients (the conditional Bayes fit's β̂) for a design; evaluated against a design at
    /// any points this is the per-draw fitted curve the ensemble averages into the Bayes fit.
    /// </summary>
    double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights);

    /// <summary>
    /// A draw of the coefficients from their conditional posterior given the knots, for readouts that are
    /// nonlinear in the fitted curve (argmax, peak counts, spans, intervals). Readouts that are affine in the curve
    /// use <see cref="PosteriorMeanCoefficients"/> instead, the lower-variance Rao-Blackwellized estimate. The
    /// default returns the mean, which carries no coefficient uncertainty; a model that keeps it reports
    /// <see cref="CoefficientUncertainty.KnotOnly"/>. Uses <paramref name="rng"/> only when it draws.
    /// </summary>
    double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => PosteriorMeanCoefficients(design, y, weights);

    /// <summary>
    /// The conditional posterior variance of the fitted curve at each row of <paramref name="readoutDesign"/>,
    /// given the knots, or null when the model has none. Pooled with the spread of the conditional means across
    /// knot draws, it gives the pointwise posterior variance by the law of total variance.
    /// </summary>
    double[]? ConditionalVariance(double[,] design, double[] y, double[]? weights, double[,] readoutDesign) => null;

    /// <summary>How <see cref="SampleCoefficients"/> and <see cref="ConditionalVariance"/> relate to the exact
    /// conditional posterior.</summary>
    CoefficientUncertainty Uncertainty => CoefficientUncertainty.KnotOnly;
}
