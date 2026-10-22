using System;
using Maths.Distributions;
using Maths.Samplers;
using Maths.Samplers.Ensemble;

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
/// <para>The support gate (Q6, step 3.7): every basis function needs column mass Σᵢ Bⱼ(xᵢ) ≥
/// <see cref="MinColumnMass"/>, one data point's worth by default. It reads only the inputs and the knots, so it is a
/// declared restriction of the knot prior's support, checked before any design is factored; 0 turns it off.</para>
/// </summary>
public sealed class SplineTarget : ISamplingTarget<KnotConfig>
{
    private readonly IBasis _basis;
    private readonly IObservationModel _model;
    private readonly IComplexityPrior _prior;
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly double[]? _weights;
    private readonly Tally<DesignRejection>? _outcomes;

    /// <summary>The support gate's default minimum column mass: one data point's worth.</summary>
    public const double DefaultMinColumnMass = 1.0;

    /// <param name="basis">Spline basis (fixed degree).</param>
    /// <param name="model">Observation model supplying the marginal likelihood.</param>
    /// <param name="prior">Prior on the knot count.</param>
    /// <param name="x">Design points, each in [0,1].</param>
    /// <param name="y">Responses, same length as <paramref name="x"/>.</param>
    /// <param name="weights">Optional per-point weights (null = unit).</param>
    /// <param name="outcomes">When given, every evaluation's admission outcome is counted into it. A chain owns its
    /// tally, so the target holds no shared mutable state.</param>
    /// <param name="minColumnMass">The support gate's minimum column mass per basis function (0 turns it off).</param>
    public SplineTarget(IBasis basis, IObservationModel model, IComplexityPrior prior,
                        double[] x, double[] y, double[]? weights = null, Tally<DesignRejection>? outcomes = null,
                        double minColumnMass = DefaultMinColumnMass)
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
        if (!(minColumnMass >= 0.0 && double.IsFinite(minColumnMass)))
            throw new ArgumentOutOfRangeException(nameof(minColumnMass), "The minimum column mass must be finite and non-negative.");

        _basis = basis;
        _model = model;
        _prior = prior;
        _x = x;
        _y = y;
        _weights = weights;
        _outcomes = outcomes;
        MinColumnMass = minColumnMass;
    }

    /// <summary>The support gate's minimum column mass Σᵢ Bⱼ(xᵢ) per basis function; 0 when the gate is off.</summary>
    public double MinColumnMass { get; }

    /// <summary>
    /// Evaluate a configuration: log prior <c>ln p(k) + ln k!</c>, and the model's marginal built once for its
    /// design (up to the constant shared across configs). Off the support, the open ordered simplex
    /// <c>0 &lt; ξ₁ &lt; … &lt; ξ_k &lt; 1</c>, or where the count prior is zero, no design is built; a design that
    /// fails the support gate is never factored.
    /// </summary>
    public ITargetEvaluation Evaluate(KnotConfig config)
    {
        ArgumentNullException.ThrowIfNull(config);
        double[] knots = config.InteriorKnots;
        double previous = 0.0;
        foreach (double knot in knots)
        {
            if (!(knot > previous && knot < 1.0)) return Rejected(DesignRejection.OutsideSimplex);
            previous = knot;
        }
        int k = knots.Length;
        double logPrior = _prior.LogPrior(k) + SpecialFunctions.LogGamma(k + 1.0);
        if (double.IsNegativeInfinity(logPrior)) return Rejected(DesignRejection.ZeroCountPrior);
        double[,] design = _basis.Design(config, _x);
        if (!HasSupport(design)) return Rejected(DesignRejection.ThinSupport);
        ITemperedMarginal marginal = _model.TemperedMarginal(design, _y, _weights);
        _outcomes?.Add(marginal.Rejection);
        return new KnotEvaluation(logPrior, marginal);
    }

    // Every column's mass Σᵢ Bⱼ(xᵢ) reaches the minimum.
    private bool HasSupport(double[,] design)
    {
        if (MinColumnMass <= 0.0) return true;
        int m = design.GetLength(0), nu = design.GetLength(1);
        for (int j = 0; j < nu; j++)
        {
            double mass = 0.0;
            for (int i = 0; i < m; i++) mass += design[i, j];
            if (mass < MinColumnMass) return false;
        }
        return true;
    }

    private ITargetEvaluation Rejected(DesignRejection reason)
    {
        _outcomes?.Add(reason);
        return KnotEvaluation.OffSupport;
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
