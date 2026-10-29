using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The cut posterior's conditional model (D09): p(θ, σ², v | y, e) for one emulator draw e at a time, with
/// y = y_m + y_s ∘ K w_e(θ) + δ + ε in physical units (D07). It binds the config to the observations and the emulator and
/// validates them together; <see cref="Target"/> gives the conditional target for draw s. The observations never update
/// the emulator: that is the cut.
/// <list type="bullet">
/// <item><b>Priors.</b> θ uniform on its box (in [0, 1] coordinates) and over each categorical parameter's levels; each
/// variance group's declared prior; v ~ N(0, s_v²·I) for the functional discrepancy. All proper.</item>
/// <item><b>Likelihood.</b> N(y; prediction + μ_δ + D·v, Σ), Σ the declared formulation's
/// (<see cref="CalibrationCovariance"/>), tempered on D05's likelihood path: the prior stays whole.</item>
/// </list>
/// </summary>
public sealed class CalibrationModel
{
    private readonly double[] _y;

    /// <summary>The model of <paramref name="config"/> on <paramref name="observations"/> (m values, the emulator's output
    /// order) and <paramref name="emulator"/>, whose parameters must be the config's.</summary>
    public CalibrationModel(CalibrationConfig config, double[] observations, ICalibrationEmulator emulator)
        : this(config, observations, emulator, null)
    {
    }

    // A model over an explicit measurement model, for the checks on priors no declaration builds.
    internal CalibrationModel(CalibrationConfig config, double[] observations, ICalibrationEmulator emulator, MeasurementModel? measurement)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(observations);
        ArgumentNullException.ThrowIfNull(emulator);
        config.Validate();
        if (!ParameterLayout.Same(config.Parameters, emulator.Parameters))
            throw new ArgumentException("The emulator serves other parameters than the config declares.", nameof(emulator));
        if (observations.Length != emulator.Outputs)
            throw new ArgumentException($"There are {observations.Length} observations and the emulator predicts {emulator.Outputs} outputs.", nameof(observations));
        foreach (double v in observations)
            if (!double.IsFinite(v)) throw new ArgumentException("The observations must be finite.", nameof(observations));
        if (emulator.DrawCount < 1) throw new ArgumentException("The emulator holds no draw.", nameof(emulator));
        Config = config;
        Emulator = emulator;
        Layout = new ParameterLayout(config.Parameters);
        Measurement = measurement ?? new MeasurementModel(config.Measurement, observations.Length);
        Covariance = CalibrationCovariance.Create(config.Covariance, Measurement, emulator);
        _y = (double[])observations.Clone();
        Discrepancy = Covariance.DiscrepancyBasis;
        DiscrepancyScale = Covariance.DiscrepancyScale;
        LogLevelPrior = 0.0;
        foreach (CategoricalParameter p in Layout.Categorical) LogLevelPrior -= Math.Log(p.Levels.Length);
    }

    /// <summary>The config.</summary>
    public CalibrationConfig Config { get; }

    /// <summary>The emulator.</summary>
    public ICalibrationEmulator Emulator { get; }

    /// <summary>The parameters by kind.</summary>
    public ParameterLayout Layout { get; }

    /// <summary>The measurement model.</summary>
    public MeasurementModel Measurement { get; }

    /// <summary>The covariance formulation, built.</summary>
    public CalibrationCovariance Covariance { get; }

    /// <summary>The observations y.</summary>
    public ReadOnlySpan<double> Observations => _y;

    /// <summary>The outputs m.</summary>
    public int Outputs => _y.Length;

    /// <summary>The functional discrepancy's basis D, m × q; null without one.</summary>
    public double[,]? Discrepancy { get; }

    /// <summary>The functional discrepancy's scale s_v; 0 without one.</summary>
    public double DiscrepancyScale { get; }

    /// <summary>q, the discrepancy coefficients sampled: 0 without a functional discrepancy.</summary>
    public int Coefficients => Discrepancy?.GetLength(1) ?? 0;

    /// <summary>Whether every prior is proper, as a b = 0 rung needs (θ's, the levels' and v's always are).</summary>
    public bool PriorsProper => Measurement.IsProper;

    internal double LogLevelPrior { get; }

    /// <summary>The conditional target p(θ, σ², v | y, e_s) for emulator draw <paramref name="draw"/>.</summary>
    public CalibrationTarget Target(int draw)
    {
        if ((uint)draw >= (uint)Emulator.DrawCount) throw new ArgumentOutOfRangeException(nameof(draw), $"The emulator holds {Emulator.DrawCount} draws.");
        return new CalibrationTarget(this, draw);
    }

    /// <summary>A state drawn from the prior: θ and the levels uniform, each variance from its group's prior, v at 0.</summary>
    public CalibrationState Start(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        var unit = new double[Layout.Continuous.Count];
        for (int i = 0; i < unit.Length; i++) unit[i] = rng.NextDouble();
        var levels = new int[Layout.Categorical.Count];
        for (int j = 0; j < levels.Length; j++) levels[j] = rng.NextInt(Layout.Categorical[j].Levels.Length);
        var variances = new double[Measurement.Groups];
        Measurement.SamplePrior(rng, variances);
        return new CalibrationState(unit, levels, variances, new double[Coefficients]);
    }

    /// <summary>
    /// A state drawn exactly from the prior, as an SMC population's first stage needs: θ and the levels uniform, each
    /// variance from its group's prior, v ~ N(0, s_v²·I). <see cref="Start"/> instead puts v at 0, a chain's start. Refused
    /// when a variance group's prior is improper.
    /// </summary>
    public CalibrationState DrawPrior(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!PriorsProper) throw new InvalidOperationException("A variance group's prior is improper: there is no prior to draw from.");
        var unit = new double[Layout.Continuous.Count];
        for (int i = 0; i < unit.Length; i++) unit[i] = rng.NextDouble();
        var levels = new int[Layout.Categorical.Count];
        for (int j = 0; j < levels.Length; j++) levels[j] = rng.NextInt(Layout.Categorical[j].Levels.Length);
        var variances = new double[Measurement.Groups];
        Measurement.SamplePrior(rng, variances);
        var coefficients = new double[Coefficients];
        for (int k = 0; k < coefficients.Length; k++) coefficients[k] = DiscrepancyScale * Normal.SampleStandard(rng);
        return new CalibrationState(unit, levels, variances, coefficients);
    }

    /// <summary>θ's continuous part of <paramref name="state"/> in physical coordinates, through each parameter's recorded
    /// transform (<see cref="ContinuousParameter.ToPhysical"/>).</summary>
    public double[] Physical(CalibrationState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        var physical = new double[Layout.Continuous.Count];
        Layout.ToPhysical(state.Unit, physical);
        return physical;
    }

    /// <summary>The discrepancy δ = μ_δ + D·v of <paramref name="state"/> (m values): μ_δ alone, or 0, without a functional
    /// discrepancy.</summary>
    public double[] Delta(CalibrationState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        return Delta(state.Coefficients);
    }

    /// <summary>The discrepancy δ = μ_δ + D·v at coefficients <paramref name="coefficients"/>, as a kept draw holds them.</summary>
    public double[] Delta(ReadOnlySpan<double> coefficients)
    {
        if (coefficients.Length != Coefficients) throw new ArgumentException($"The discrepancy has {Coefficients} coefficients.", nameof(coefficients));
        var delta = Covariance.Mean is { } mean ? (double[])mean.Clone() : new double[Outputs];
        if (Discrepancy is { } d)
            for (int j = 0; j < Outputs; j++)
                for (int k = 0; k < d.GetLength(1); k++) delta[j] += d[j, k] * coefficients[k];
        return delta;
    }

    /// <summary>The prior's log density at <paramref name="state"/>: −∞ off the support (θ outside [0, 1], a level out of
    /// range, a variance not positive and finite, v not finite).</summary>
    public double LogPrior(CalibrationState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        if (!Shaped(state)) throw new ArgumentException("The state does not have the model's dimensions.", nameof(state));
        foreach (double u in state.Unit)
            if (!(u >= 0.0 && u <= 1.0)) return double.NegativeInfinity;
        for (int j = 0; j < state.Levels.Length; j++)
            if ((uint)state.Levels[j] >= (uint)Layout.Categorical[j].Levels.Length) return double.NegativeInfinity;
        foreach (double v in state.Variances)
            if (!(v > 0.0 && double.IsFinite(v))) return double.NegativeInfinity;
        double prior = LogLevelPrior + Measurement.LogPrior(state.Variances);
        if (Discrepancy is not null)
        {
            double s2 = DiscrepancyScale * DiscrepancyScale, sum = 0.0;
            foreach (double v in state.Coefficients)
            {
                if (!double.IsFinite(v)) return double.NegativeInfinity;
                sum += v * v;
            }
            prior += -0.5 * state.Coefficients.Length * Math.Log(2.0 * Math.PI * s2) - 0.5 * sum / s2;
        }
        return prior;
    }

    internal bool Shaped(CalibrationState state)
        => state.Unit.Length == Layout.Continuous.Count && state.Levels.Length == Layout.Categorical.Count
        && state.Variances.Length == Measurement.Groups && state.Coefficients.Length == Coefficients;
}

/// <summary>
/// A state of the conditional chain, immutable: θ's continuous part in [0, 1] coordinates, its levels, the group
/// variances σ² and the discrepancy coefficients v. A change makes a new state, which a target evaluates afresh, so no
/// likelihood outlives the values it was computed from. The prediction and the covariance factor are kept with the state
/// for the last two targets (emulator draws) that read it, so a bridge between two draws reads both from the cache: a
/// change to v keeps both, to σ² keeps the prediction, to θ keeps the factor, and a state read by a third target
/// recomputes both (C-O1's stale likelihood, corrected by construction).
/// </summary>
public sealed class CalibrationState
{
    private readonly double[] _unit, _variances, _coefficients;
    private readonly int[] _levels;

    public CalibrationState(double[] unit, int[] levels, double[] variances, double[] coefficients)
    {
        ArgumentNullException.ThrowIfNull(unit);
        ArgumentNullException.ThrowIfNull(levels);
        ArgumentNullException.ThrowIfNull(variances);
        ArgumentNullException.ThrowIfNull(coefficients);
        _unit = (double[])unit.Clone();
        _levels = (int[])levels.Clone();
        _variances = (double[])variances.Clone();
        _coefficients = (double[])coefficients.Clone();
    }

    private CalibrationState(double[] unit, int[] levels, double[] variances, double[] coefficients, PredictionCache? prediction, FactorCache? factor)
    {
        _unit = unit;
        _levels = levels;
        _variances = variances;
        _coefficients = coefficients;
        Prediction = prediction;
        Factor = factor;
    }

    /// <summary>θ's continuous part in [0, 1] coordinates.</summary>
    public ReadOnlySpan<double> Unit => _unit;

    /// <summary>θ's levels, an index into each categorical parameter's levels.</summary>
    public ReadOnlySpan<int> Levels => _levels;

    /// <summary>The group variances σ².</summary>
    public ReadOnlySpan<double> Variances => _variances;

    /// <summary>The discrepancy coefficients v.</summary>
    public ReadOnlySpan<double> Coefficients => _coefficients;

    // Written once by a target's evaluation; a racing second write stores the same values.
    internal PredictionCache? Prediction { get; set; }
    internal FactorCache? Factor { get; set; }

    /// <summary>This state with θ's continuous part at <paramref name="unit"/>: the factor kept.</summary>
    public CalibrationState WithUnit(double[] unit) => new((double[])unit.Clone(), _levels, _variances, _coefficients, null, Factor);

    /// <summary>This state with <paramref name="levels"/>: the factor kept.</summary>
    public CalibrationState WithLevels(int[] levels) => new(_unit, (int[])levels.Clone(), _variances, _coefficients, null, Factor);

    /// <summary>This state with <paramref name="variances"/>: the prediction kept.</summary>
    public CalibrationState WithVariances(double[] variances) => new(_unit, _levels, (double[])variances.Clone(), _coefficients, Prediction, null);

    /// <summary>This state with <paramref name="coefficients"/>: the prediction and the factor kept.</summary>
    public CalibrationState WithCoefficients(double[] coefficients) => new(_unit, _levels, _variances, (double[])coefficients.Clone(), Prediction, Factor);

    // With levels whose prediction the caller has computed for `owner` (the categorical enumeration's).
    internal CalibrationState WithLevels(int[] levels, CalibrationTarget owner, double[] prediction)
        => new(_unit, (int[])levels.Clone(), _variances, _coefficients, new PredictionCache(owner, prediction), Factor);

    // With levels whose predictions the caller has computed, for one target or two (a bridge's enumeration).
    internal CalibrationState WithLevels(int[] levels, PredictionCache prediction)
        => new(_unit, (int[])levels.Clone(), _variances, _coefficients, prediction, Factor);

    // A distinct state with the same values and caches: an SMC offspring, whose later cache writes stay its own.
    internal CalibrationState Fork() => new(_unit, _levels, _variances, _coefficients, Prediction, Factor);
}

// A target's prediction for a state, keeping the previous target's beside it (Other has no Other of its own).
internal sealed record PredictionCache(CalibrationTarget Owner, double[] Values, PredictionCache? Other = null)
{
    public double[]? For(CalibrationTarget owner)
        => ReferenceEquals(Owner, owner) ? Values : Other is { } other && ReferenceEquals(other.Owner, owner) ? other.Values : null;

    public static PredictionCache Add(PredictionCache? cache, CalibrationTarget owner, double[] values)
        => new(owner, values, cache is null ? null : cache with { Other = null });
}

// A target's covariance factor for a state, keeping the previous target's beside it.
internal sealed record FactorCache(CalibrationTarget Owner, CovarianceFactor Factor, FactorCache? Other = null)
{
    public CovarianceFactor? For(CalibrationTarget owner)
        => ReferenceEquals(Owner, owner) ? Factor : Other is { } other && ReferenceEquals(other.Owner, owner) ? other.Factor : null;

    public static FactorCache Add(FactorCache? cache, CalibrationTarget owner, CovarianceFactor factor)
        => new(owner, factor, cache is null ? null : cache with { Other = null });
}

/// <summary>
/// The conditional target p(θ, σ², v | y, e_s) for one emulator draw (<see cref="CalibrationModel.Target"/>), as an
/// <see cref="ISamplingTarget{TState}"/> on D05's likelihood path: the samplers, the tempering ladder and stage 8's SMC
/// read it through <see cref="CalibrationEvaluation"/>.
/// </summary>
public sealed class CalibrationTarget : ISamplingTarget<CalibrationState>
{
    internal CalibrationTarget(CalibrationModel model, int draw)
    {
        Model = model;
        Draw = draw;
        Noise = model.Emulator.Noise(draw);
    }

    /// <summary>The model.</summary>
    public CalibrationModel Model { get; }

    /// <summary>The emulator draw s this target conditions on.</summary>
    public int Draw { get; }

    /// <summary>The draw's noise covariance.</summary>
    public EmulatorNoise Noise { get; }

    /// <inheritdoc/>
    public ITargetEvaluation Evaluate(CalibrationState state) => EvaluateState(state);

    /// <summary>The evaluation of <paramref name="state"/>, with the prediction and factor it reads.</summary>
    public CalibrationEvaluation EvaluateState(CalibrationState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        double prior = Model.LogPrior(state);
        if (double.IsNegativeInfinity(prior)) return new CalibrationEvaluation(this, prior, double.NegativeInfinity, null, null);
        double[] prediction = Prediction(state);
        CovarianceFactor factor = Factor(state);
        double[] residual = Residual(prediction, state);
        double logLikelihood = -0.5 * (Model.Outputs * Math.Log(2.0 * Math.PI) + factor.LogDet + factor.InverseQuadraticForm(residual));
        return new CalibrationEvaluation(this, prior, logLikelihood, prediction, factor);
    }

    /// <summary>The prediction at <paramref name="state"/>'s θ under this draw, physical units; kept with the state.</summary>
    public double[] Prediction(CalibrationState state)
    {
        PredictionCache? cache = state.Prediction;
        if (cache?.For(this) is { } cached) return cached;
        var prediction = new double[Model.Outputs];
        Model.Emulator.Predict(Draw, Model.Physical(state), state.Levels, prediction);
        foreach (double v in prediction)
            if (!double.IsFinite(v)) throw new InvalidOperationException("The emulator predicted a value that is not finite.");
        state.Prediction = PredictionCache.Add(cache, this, prediction);
        return prediction;
    }

    /// <summary>Σ at <paramref name="state"/>'s variances under this draw, factored; kept with the state.</summary>
    public CovarianceFactor Factor(CalibrationState state)
    {
        FactorCache? cache = state.Factor;
        if (cache?.For(this) is { } cached) return cached;
        CovarianceFactor factor = Model.Covariance.Factor(state.Variances, Noise);
        state.Factor = FactorCache.Add(cache, this, factor);
        return factor;
    }

    /// <summary>y − prediction − δ at <paramref name="state"/>.</summary>
    public double[] Residual(double[] prediction, CalibrationState state)
    {
        ReadOnlySpan<double> y = Model.Observations;
        double[] delta = Model.Delta(state);
        var residual = new double[y.Length];
        for (int j = 0; j < y.Length; j++) residual[j] = y[j] - prediction[j] - delta[j];
        return residual;
    }
}

/// <summary>
/// One state's evaluation under a <see cref="CalibrationTarget"/>: the untempered log prior and the Gaussian
/// log-likelihood ℓ₁ (complete, so ladder evidence can read it), separable on the likelihood path, ℓ_b = b·ℓ₁. At b = 0
/// the chain samples the prior, so every prior must be proper there (the propriety rule): an improper one is refused.
/// </summary>
public sealed class CalibrationEvaluation : ITargetEvaluation
{
    private readonly double _logLikelihood;

    internal CalibrationEvaluation(CalibrationTarget target, double logPrior, double logLikelihood, double[]? prediction, CovarianceFactor? factor)
    {
        Target = target;
        LogPrior = logPrior;
        _logLikelihood = logLikelihood;
        Prediction = prediction;
        Factor = factor;
    }

    /// <summary>The target that evaluated the state.</summary>
    public CalibrationTarget Target { get; }

    /// <inheritdoc/>
    public double LogPrior { get; }

    /// <summary>ℓ₁, the untempered Gaussian log-likelihood; −∞ off the support.</summary>
    public double LogLikelihoodValue => _logLikelihood;

    /// <summary>The prediction the likelihood read; null off the support.</summary>
    public double[]? Prediction { get; }

    /// <summary>The covariance factor the likelihood read; null off the support.</summary>
    public CovarianceFactor? Factor { get; }

    /// <inheritdoc/>
    public double LogLikelihood(double inverseTemperature)
    {
        if (!(inverseTemperature >= 0.0 && inverseTemperature <= 1.0)) throw new ArgumentOutOfRangeException(nameof(inverseTemperature));
        if (inverseTemperature == 0.0 && !Target.Model.PriorsProper)
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature), "A b = 0 rung samples the prior, and a variance group's prior is improper.");
        return double.IsNegativeInfinity(_logLikelihood) ? double.NegativeInfinity : inverseTemperature * _logLikelihood;
    }

    /// <inheritdoc/>
    public double LogLikelihoodSlope(double inverseTemperature) => _logLikelihood;
}
