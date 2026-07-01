using System;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The reversible-jump target for free-knot regression: the marginalized posterior over knot configurations,
/// exposed to the engine as <c>Δ log-marginal-likelihood + Δ log-prior</c> with coefficients integrated out.
/// Holds the fixed run data (design points + response) and composes a <see cref="SplineBasis"/>, an
/// <see cref="IObservationModel"/> and an <see cref="IComplexityPrior"/>.
/// </summary>
public sealed class SplineTarget : IRjTarget<KnotConfig>
{
    private readonly IBasis _basis;
    private readonly IObservationModel _model;
    private readonly IComplexityPrior _prior;
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly double[]? _weights;

    /// <param name="basis">Spline basis (fixed degree).</param>
    /// <param name="model">Observation model supplying the marginal likelihood.</param>
    /// <param name="prior">Prior on the knot count.</param>
    /// <param name="x">Design points, each in [0,1].</param>
    /// <param name="y">Responses, same length as <paramref name="x"/>.</param>
    /// <param name="weights">Optional per-point weights (null = unit).</param>
    public SplineTarget(IBasis basis, IObservationModel model, IComplexityPrior prior,
                        double[] x, double[] y, double[]? weights = null)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(prior);
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        if (x.Length != y.Length)
            throw new ArgumentException("x and y must have the same length.");
        if (weights is not null && weights.Length != x.Length)
            throw new ArgumentException("weights must match the design length.", nameof(weights));

        _basis = basis;
        _model = model;
        _prior = prior;
        _x = x;
        _y = y;
        _weights = weights;
    }

    /// <summary>Log marginal posterior at a configuration (up to the constant shared across configs).</summary>
    public double LogPosterior(KnotConfig config)
    {
        ArgumentNullException.ThrowIfNull(config);
        return _model.LogMarginalLikelihood(_basis.Design(config, _x), _y, _weights)
             + _prior.LogPrior(config.Count);
    }
}
