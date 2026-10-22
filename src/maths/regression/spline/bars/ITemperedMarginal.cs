using System;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// One design's log marginal likelihood as a function of the likelihood's inverse temperature b, built once by
/// <see cref="IObservationModel.TemperedMarginal"/>. It comes in two forms, both of which run from the prior at b = 0
/// to the posterior at b = 1 with the knot prior untempered:
/// <list type="bullet">
/// <item><b>Integrated:</b> the coefficients and scale integrated against <c>p(y | β, σ)^b</c>, not linear in b. The
/// Normal model's, from its closed form; it is the form a sampler that keeps σ² in its state (BASS) needs.</item>
/// <item><b>Linear:</b> <c>b·ℓ₁</c>, the untempered marginal raised to b. The Laplace-approximated models', where the
/// integrated form would need a fit per rung. Its Occam penalty at a hot rung is b times the untempered one, much
/// weaker than the integrated form's.</item>
/// </list>
/// Values omit the same kind of constant as <see cref="IObservationModel.LogMarginalLikelihood"/>: it may depend on b
/// and the data, never on the design.
/// </summary>
public interface ITemperedMarginal
{
    /// <summary>The log marginal likelihood at <paramref name="inverseTemperature"/>; b = 1 is untempered.</summary>
    double LogMarginalLikelihood(double inverseTemperature);

    /// <summary>Why the design is outside the model's support, or <see cref="DesignRejection.None"/> when it is admitted.</summary>
    DesignRejection Rejection => DesignRejection.None;
}

/// <summary>
/// Why a knot configuration scores −∞: outside the prior's support, or refused by the observation model. Counted by
/// reason (<see cref="SplineTarget"/>), so the truncation the support rules and admission impose stays visible.
/// </summary>
public enum DesignRejection
{
    /// <summary>Admitted: the configuration has a finite score.</summary>
    None,

    /// <summary>The knots are not on the open ordered simplex 0 &lt; ξ₁ &lt; … &lt; ξ_k &lt; 1.</summary>
    OutsideSimplex,

    /// <summary>The count prior gives the knot count probability zero.</summary>
    ZeroCountPrior,

    /// <summary>No more observations than basis functions (m ≤ ν): the marginal is improper.</summary>
    TooFewObservations,

    /// <summary>A factored matrix failed admission: rank-deficient or ill-conditioned.</summary>
    Inadmissible,

    /// <summary>The model's iterative fit (IRLS) did not settle within its iteration cap.</summary>
    NotConverged,

    /// <summary>An identically zero response leaves no posterior for the scale.</summary>
    ZeroResponse,
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

// The linear form, b·ℓ₁, for b in [0, 1]. A marginal of −∞ stays −∞ at every b, b = 0 included, so a configuration
// outside the model's support stays outside it.
internal sealed class LinearMarginal : ITemperedMarginal
{
    private readonly double _value;
    public LinearMarginal(double value) => _value = value;

    // A model refuses a design with ImpossibleMarginal and its reason; the linear form is for admitted ones.
    public static ITemperedMarginal Of(double value, DesignRejection rejection)
        => rejection == DesignRejection.None ? new LinearMarginal(value) : new ImpossibleMarginal(rejection);

    public double LogMarginalLikelihood(double inverseTemperature)
    {
        if (!(inverseTemperature >= 0.0 && inverseTemperature <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature), "The inverse temperature must lie in [0, 1].");
        return double.IsNegativeInfinity(_value) ? double.NegativeInfinity : inverseTemperature * _value;
    }
}

// A configuration outside the model's support: −∞ at every temperature, with the reason.
internal sealed class ImpossibleMarginal : ITemperedMarginal
{
    public ImpossibleMarginal(DesignRejection rejection) => Rejection = rejection;
    public DesignRejection Rejection { get; }
    public double LogMarginalLikelihood(double inverseTemperature) => double.NegativeInfinity;
}
