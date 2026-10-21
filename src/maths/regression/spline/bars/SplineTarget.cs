using System;
using Maths.Distributions;
using Maths.Samplers;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The reversible-jump target for free-knot regression: the marginalized posterior over knot configurations, with
/// coefficients and scale integrated out. Its log prior is the count prior and the knot measure; its log likelihood
/// is the observation model's marginal, at any inverse temperature the model supports. Given k, the knot locations
/// are the order statistics of k uniforms on (0,1), with density <c>k!</c> on the ordered simplex, so the count
/// prior <c>p(k)</c> is the knot-count marginal (R01).
/// Holds the fixed run data (design points + response) and composes a <see cref="SplineBasis"/>, an
/// <see cref="IObservationModel"/> and an <see cref="IComplexityPrior"/>.
/// </summary>
public sealed class SplineTarget : ISamplingTarget<KnotConfig>
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

    /// <summary>
    /// Evaluate a configuration: log prior <c>ln p(k) + ln k!</c>, and the model's marginal built once for its
    /// design (up to the constant shared across configs). Off the support, the open ordered simplex
    /// <c>0 &lt; ξ₁ &lt; … &lt; ξ_k &lt; 1</c>, or where the count prior is zero, no design is built.
    /// </summary>
    public ITargetEvaluation Evaluate(KnotConfig config)
    {
        ArgumentNullException.ThrowIfNull(config);
        double[] knots = config.InteriorKnots;
        double previous = 0.0;
        foreach (double knot in knots)
        {
            if (!(knot > previous && knot < 1.0)) return KnotEvaluation.OffSupport;
            previous = knot;
        }
        int k = knots.Length;
        double logPrior = _prior.LogPrior(k) + SpecialFunctions.LogGamma(k + 1.0);
        if (double.IsNegativeInfinity(logPrior)) return KnotEvaluation.OffSupport;
        return new KnotEvaluation(logPrior, _model.TemperedMarginal(_basis.Design(config, _x), _y, _weights));
    }

    private sealed class KnotEvaluation : ITargetEvaluation
    {
        public static readonly KnotEvaluation OffSupport = new(double.NegativeInfinity, null);
        private readonly ITemperedMarginal? _marginal;

        public KnotEvaluation(double logPrior, ITemperedMarginal? marginal)
        {
            LogPrior = logPrior;
            _marginal = marginal;
        }

        public double LogPrior { get; }

        public double LogLikelihood(double inverseTemperature)
            => _marginal?.LogMarginalLikelihood(inverseTemperature) ?? double.NegativeInfinity;
    }
}
