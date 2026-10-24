using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// One knot configuration's conditional posterior, fitted once, so the readouts a draw takes from it share the fit: the
/// coefficients' mean (the affine readouts, Rao–Blackwellized), the fitted curve's conditional variance at readout rows,
/// and a joint draw of the coefficients and σ² (the nonlinear readouts and the retained state). Each read returns fresh
/// arrays.
/// </summary>
public interface IConditionalFit
{
    /// <summary>As <see cref="IObservationModel.PosteriorMeanCoefficients"/>.</summary>
    double[] MeanCoefficients();

    /// <summary>As <see cref="IObservationModel.ConditionalVariance"/>.</summary>
    double[]? ConditionalVariance(double[,] readoutDesign);

    /// <summary>As <see cref="IObservationModel.SampleConditional"/>.</summary>
    ConditionalDraw Sample(Xoshiro256PlusPlus rng);
}

/// <summary>The seam's default for a model that keeps no shared fit: each read calls the model's own method, which fits
/// afresh.</summary>
internal sealed class SeparateFits(IObservationModel model, double[,] design, double[] y, double[]? weights) : IConditionalFit
{
    public double[] MeanCoefficients() => model.PosteriorMeanCoefficients(design, y, weights);

    public double[]? ConditionalVariance(double[,] readoutDesign) => model.ConditionalVariance(design, y, weights, readoutDesign);

    public ConditionalDraw Sample(Xoshiro256PlusPlus rng) => model.SampleConditional(design, y, weights, rng);
}

/// <summary>A fit whose draw is its mean, with no conditional variance: the knot-only uncertainty of a generalized
/// posterior (<see cref="CoefficientUncertainty.KnotOnly"/>).</summary>
internal sealed class MeanOnlyFit(double[] mean) : IConditionalFit
{
    public double[] MeanCoefficients() => (double[])mean.Clone();

    public double[]? ConditionalVariance(double[,] readoutDesign) => null;

    public ConditionalDraw Sample(Xoshiro256PlusPlus rng) => new((double[])mean.Clone(), double.NaN);
}
