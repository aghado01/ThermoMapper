using System;
using Maths.Distributions;
using Maths.Rng;
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
    public ITargetEvaluation Evaluate(KnotConfig config) => Evaluate(config, _outcomes);

    /// <summary>
    /// Exact draws from this target's b = 0 member on <paramref name="path"/>, restricted to its support by rejection,
    /// for a ladder's reference rung (OQ2). On both paths the knot locations are the order statistics of k uniforms, as
    /// the target's own knot measure has them; the paths differ in the count law.
    /// <list type="bullet">
    /// <item><b>Powered:</b> k ~ Poisson(1). The powered density's b = 0 limit is the base measure itself, counting on k
    /// times Lebesgue measure on the ordered simplex, whose k-th piece has volume 1/k!: normalized, Poisson(1) with
    /// uniform locations. It needs nothing of the model but its support, so it serves every model.</item>
    /// <item><b>Likelihood:</b> k from the count prior, normalized over the counts that a design can admit: the knot
    /// prior itself. It needs a model that <see cref="IObservationModel.CanTemperToPrior"/>.</item>
    /// </list>
    /// The support includes every refusal the target makes. The declared prior restrictions (the support gate, too few
    /// observations) are the prior's own truncation; numerical refusals restrict it further, and the draw is exact for
    /// the restricted law the ladder's swaps evaluate. Reference draws are not counted in the outcome tally.
    /// </summary>
    public IReferenceSampler<KnotConfig> Reference(TemperingPath path)
    {
        if (path is not (TemperingPath.Powered or TemperingPath.Likelihood)) throw new ArgumentOutOfRangeException(nameof(path));
        if (path == TemperingPath.Likelihood && !_model.CanTemperToPrior)
            throw new InvalidOperationException($"{_model.GetType().Name} does not serve b = 0 on the likelihood path.");
        return new KnotReference(this, path);
    }

    private ITargetEvaluation Evaluate(KnotConfig config, Tally<DesignRejection>? outcomes)
    {
        ArgumentNullException.ThrowIfNull(config);
        double[] knots = config.InteriorKnots;
        double previous = 0.0;
        foreach (double knot in knots)
        {
            if (!(knot > previous && knot < 1.0)) return Rejected(DesignRejection.OutsideSimplex, outcomes);
            previous = knot;
        }
        int k = knots.Length;
        double logPrior = _prior.LogPrior(k) + SpecialFunctions.LogGamma(k + 1.0);
        if (double.IsNegativeInfinity(logPrior)) return Rejected(DesignRejection.ZeroCountPrior, outcomes);
        double[,] design = _basis.Design(config, _x);
        if (!HasSupport(design)) return Rejected(DesignRejection.ThinSupport, outcomes);
        ITemperedMarginal marginal = _model.TemperedMarginal(design, _y, _weights);
        outcomes?.Add(marginal.Rejection);
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

    private static ITargetEvaluation Rejected(DesignRejection reason, Tally<DesignRejection>? outcomes)
    {
        outcomes?.Add(reason);
        return KnotEvaluation.OffSupport;
    }

    // Rejection from the unrestricted b = 0 law: draw k, then k sorted uniforms, until the target admits the draw. The
    // count law is tabulated until its terms fall 60 nats below the largest, past the data size, so its truncation is
    // below e^−60 of the mass; a count law that never decays that far is not proper.
    private sealed class KnotReference : IReferenceSampler<KnotConfig>
    {
        private const int MaxAttempts = 1_000_000;
        private const int MaxCount = 100_000;
        private readonly SplineTarget _target;
        private readonly Temperature _rung;
        private readonly double[] _cumulative;   // the count law's CDF over the tabulated counts

        public KnotReference(SplineTarget target, TemperingPath path)
        {
            _target = target;
            Path = path;
            _rung = new Temperature(0.0, path);
            int m = target._x.Length;
            var logMass = new System.Collections.Generic.List<double>();
            double max = double.NegativeInfinity;
            for (int k = 0; ; k++)
            {
                double lp = path == TemperingPath.Powered ? -SpecialFunctions.LogGamma(k + 1.0) : target._prior.LogPrior(k);
                logMass.Add(lp);
                if (lp > max) max = lp;
                if (k >= m && (double.IsNegativeInfinity(lp) || lp < max - 60.0)) break;
                if (k == MaxCount) throw new InvalidOperationException("The count prior does not decay; it is not proper.");
            }
            if (double.IsNegativeInfinity(max)) throw new InvalidOperationException("The count prior has no mass.");
            _cumulative = new double[logMass.Count];
            double total = 0.0;
            for (int k = 0; k < logMass.Count; k++)
            {
                total += Math.Exp(logMass[k] - max);
                _cumulative[k] = total;
            }
            for (int k = 0; k < _cumulative.Length; k++) _cumulative[k] /= total;
        }

        public TemperingPath Path { get; }

        public long Attempts { get; private set; }

        public (KnotConfig State, ITargetEvaluation Evaluation) Draw(Xoshiro256PlusPlus rng)
        {
            ArgumentNullException.ThrowIfNull(rng);
            for (int attempt = 0; attempt < MaxAttempts; attempt++)
            {
                double u = rng.NextDouble();
                int k = 0;
                while (k < _cumulative.Length - 1 && u >= _cumulative[k]) k++;
                var knots = new double[k];
                for (int j = 0; j < k; j++) knots[j] = rng.NextDouble();
                Array.Sort(knots);
                var config = new KnotConfig(knots);
                Attempts++;
                ITargetEvaluation evaluation = _target.Evaluate(config, null);
                if (!double.IsNegativeInfinity(_rung.LogDensity(evaluation))) return (config, evaluation);
            }
            throw new InvalidOperationException($"No admitted configuration in {MaxAttempts} reference draws.");
        }
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

        public double LogLikelihoodSlope(double inverseTemperature)
            => _marginal is null
                ? throw new NotSupportedException("A configuration off the support has no slope in b.")
                : _marginal.LogMarginalLikelihoodSlope(inverseTemperature);
    }
}
