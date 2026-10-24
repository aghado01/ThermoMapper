namespace Maths.Regression.Spline.Bars;

/// <summary>
/// One draw from an observation model's conditional posterior given the knots: the spline coefficients and the noise
/// variance σ² drawn with them. σ² is NaN for a model without a noise scale (Poisson, Whittle, Gibbs loss) and when
/// no draw is possible. <see cref="CoefficientUncertainty"/> says how exact the coefficients are.
/// </summary>
public readonly record struct ConditionalDraw(double[] Coefficients, double Sigma2);
