using System;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// One design's log marginal likelihood as a function of the likelihood's inverse temperature b, built once by
/// <see cref="IObservationModel.TemperedMarginal"/>. For a collapsed model it integrates the coefficients and scale
/// against <c>p(y | β, σ)^b</c>, which is not b times the untempered marginal. Values omit the same kind of constant
/// as <see cref="IObservationModel.LogMarginalLikelihood"/>: it may depend on b and the data, never on the design.
/// </summary>
public interface ITemperedMarginal
{
    /// <summary>The log marginal likelihood at <paramref name="inverseTemperature"/>; b = 1 is untempered.</summary>
    double LogMarginalLikelihood(double inverseTemperature);
}

// The default for a model without a tempered closed form: its untempered marginal, at b = 1 only.
internal sealed class UntemperedMarginal : ITemperedMarginal
{
    private readonly double _value;
    public UntemperedMarginal(double value) => _value = value;

    public double LogMarginalLikelihood(double inverseTemperature) => inverseTemperature == 1.0
        ? _value
        : throw new NotSupportedException("This observation model has no tempered marginal; it serves b = 1 only.");
}

// A configuration outside the model's support: −∞ at every temperature.
internal sealed class ImpossibleMarginal : ITemperedMarginal
{
    public static readonly ImpossibleMarginal Instance = new();
    public double LogMarginalLikelihood(double inverseTemperature) => double.NegativeInfinity;
}
