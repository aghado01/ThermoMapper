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

    /// <summary>Effective number of parameters for the complexity penalty (ν for the Normal model).</summary>
    int EffectiveDimension(double[,] design, double[]? weights);

    /// <summary>
    /// Posterior mean coefficients (the conditional Bayes fit's β̂) for a design; evaluated against a design at
    /// any points this is the per-draw fitted curve the ensemble averages into the Bayes fit.
    /// </summary>
    double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights);

    /// <summary>
    /// A coefficient vector for the per-draw readout. Defaults to the posterior <i>mean</i> (the Rao-Blackwell
    /// choice — exact for the Normal model, where its credible intervals are already integrated over the chain);
    /// models whose conditional β-posterior is non-Gaussian override this with an actual draw (e.g.
    /// <see cref="PoissonModel"/> corrects the skew at small counts). Uses <paramref name="rng"/> only when it draws.
    /// </summary>
    double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => PosteriorMeanCoefficients(design, y, weights);
}
