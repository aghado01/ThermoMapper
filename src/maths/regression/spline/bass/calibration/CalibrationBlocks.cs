using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// What every calibration block shares: its target, the rung it reads from the chain (D05's likelihood path only, where
/// b tempers the likelihood and the prior stays whole) and the hand-back of a new state with its fresh evaluation, so the
/// chain never holds a likelihood older than its state (C-O1's stale likelihood, corrected).
/// </summary>
public abstract class CalibrationBlock : IBlock<ReversibleJumpChain<CalibrationState>>
{
    private protected CalibrationBlock(CalibrationTarget target)
    {
        ArgumentNullException.ThrowIfNull(target);
        Target = target;
    }

    /// <summary>The target the block updates against.</summary>
    public CalibrationTarget Target { get; }

    /// <inheritdoc/>
    public abstract string Name { get; }

    /// <inheritdoc/>
    public abstract bool Update(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng);

    // The chain's inverse temperature, refused off the likelihood path or when the chain samples another target.
    private protected double Rung(ReversibleJumpChain<CalibrationState> chain)
    {
        ArgumentNullException.ThrowIfNull(chain);
        if (chain.Temperature.Path != TemperingPath.Likelihood)
            throw new InvalidOperationException("Calibration tempers on D05's likelihood path; its blocks draw from that path's conditionals.");
        if (chain.CurrentEvaluation is not CalibrationEvaluation e || !ReferenceEquals(e.Target, Target))
            throw new InvalidOperationException("The chain samples another target than the block's.");
        return chain.Temperature.InverseTemperature;
    }

    private protected void Hand(ReversibleJumpChain<CalibrationState> chain, CalibrationState next, CalibrationEvaluation evaluation)
        => chain.SetState(next, evaluation);
}

/// <summary>
/// The functional discrepancy's coefficients v drawn from their Gaussian conditional at the chain's rung b: precision
/// I/s_v² + b·DᵀΣ⁻¹D and mean its inverse times b·DᵀΣ⁻¹(y − prediction − μ_δ), the likelihood path's (R's covariance,
/// <c>discrep.S / itl</c>, is the powered path's, C-O1). At b = 0 it is v's prior. The prediction and the factor are the
/// state's, kept; the likelihood is recomputed at the new v.
/// </summary>
public sealed class DiscrepancyBlock : CalibrationBlock
{
    private readonly GaussianConditional _draw;

    public DiscrepancyBlock(CalibrationTarget target) : base(target)
    {
        if (target.Model.Coefficients == 0) throw new ArgumentException("The model has no functional discrepancy to draw.", nameof(target));
        _draw = new GaussianConditional(target.Model.Coefficients);
    }

    /// <inheritdoc/>
    public override string Name => "discrepancy";

    /// <inheritdoc/>
    public override bool Update(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng)
    {
        double b = Rung(state);
        CalibrationState current = state.Current;
        (DenseCholesky precision, double[] rhs) = Conditional(current, b);
        var v = new double[Target.Model.Coefficients];
        _draw.Draw(precision, rhs, rng, v);
        CalibrationState next = current.WithCoefficients(v);
        Hand(state, next, Target.EvaluateState(next));
        return true;
    }

    /// <summary>v's conditional at <paramref name="state"/> and rung <paramref name="b"/>: its precision, factored, and
    /// the right-hand side whose solve is its mean.</summary>
    internal (DenseCholesky Precision, double[] Rhs) Conditional(CalibrationState state, double b)
    {
        CalibrationModel model = Target.Model;
        double[,] d = model.Discrepancy!;
        int m = model.Outputs, q = d.GetLength(1);
        CovarianceFactor factor = Target.Factor(state);
        double[] prediction = Target.Prediction(state);
        ReadOnlySpan<double> y = model.Observations;
        var r0 = new double[m];
        for (int j = 0; j < m; j++) r0[j] = y[j] - prediction[j] - (model.Covariance.Mean?[j] ?? 0.0);
        var solved = new double[q][];   // Σ⁻¹D's columns
        var column = new double[m];
        for (int k = 0; k < q; k++)
        {
            for (int j = 0; j < m; j++) column[j] = d[j, k];
            solved[k] = factor.Solve(column);
        }
        double inverseScale2 = 1.0 / (model.DiscrepancyScale * model.DiscrepancyScale);
        var precision = new double[q, q];
        var rhs = new double[q];
        for (int k = 0; k < q; k++)
        {
            double s = 0.0;
            for (int j = 0; j < m; j++) s += solved[k][j] * r0[j];
            rhs[k] = b * s;
            for (int l = 0; l <= k; l++)
            {
                double g = 0.0;
                for (int j = 0; j < m; j++) g += d[j, l] * solved[k][j];
                precision[k, l] = precision[l, k] = b * g + (k == l ? inverseScale2 : 0.0);
            }
        }
        var factored = new DenseCholesky(q);
        factored.Decompose(precision);
        if (!factored.Status.PositiveDefinite) throw new InvalidOperationException("The discrepancy's conditional precision did not factor.");
        return (factored, rhs);
    }
}

/// <summary>
/// Each variance group's σ²_g by a random-walk Metropolis step on u = ln σ²_g, whose target carries the Jacobian of the
/// change of variables, e^u, untempered: p(σ²_g = e^u | ·)·e^u with the likelihood at the chain's rung (R tempers the
/// Jacobian and the prior with it, the powered path's move, C-O1). Every proposal builds Σ afresh. During burn-in each
/// group's step adapts towards a 0.44 acceptance rate (<see cref="Adapting"/>); then it is held.
/// </summary>
public sealed class VarianceBlock : CalibrationBlock
{
    private readonly double[] _logStep;
    private readonly int[] _windowAttempts, _windowAccepts;
    private int _windows;
    private const int Window = 50;
    private const double TargetRate = 0.44;

    /// <param name="target">The target.</param>
    /// <param name="step">Each group's starting step on ln σ².</param>
    public VarianceBlock(CalibrationTarget target, double step = 0.5) : base(target)
    {
        if (!(step > 0.0 && double.IsFinite(step))) throw new ArgumentOutOfRangeException(nameof(step));
        int groups = target.Model.Measurement.Groups;
        _logStep = new double[groups];
        Array.Fill(_logStep, Math.Log(step));
        _windowAttempts = new int[groups];
        _windowAccepts = new int[groups];
    }

    /// <inheritdoc/>
    public override string Name => "variances";

    /// <summary>Whether the steps adapt (burn-in); set false to hold them.</summary>
    public bool Adapting { get; set; }

    /// <summary>Each group's current step on ln σ².</summary>
    public IReadOnlyList<double> Steps => Array.ConvertAll(_logStep, Math.Exp);

    /// <summary>Attempts and acceptances per group since the block was built.</summary>
    public long Attempts { get; private set; }

    /// <inheritdoc cref="Attempts"/>
    public long Accepts { get; private set; }

    /// <inheritdoc/>
    public override bool Update(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng)
    {
        _ = Rung(state);
        bool any = false;
        for (int g = 0; g < _logStep.Length; g++)
        {
            CalibrationState current = state.Current;
            double[] variances = current.Variances.ToArray();
            double u = Math.Log(variances[g]);
            double logTarget = state.CurrentLogDensity + u;
            CalibrationState? candidate = null;
            CalibrationEvaluation? evaluation = null;
            double LogDensity(double proposal)
            {
                double[] next = (double[])variances.Clone();
                next[g] = Math.Exp(proposal);
                candidate = current.WithVariances(next);
                evaluation = Target.EvaluateState(candidate);
                double density = state.Temperature.LogDensity(evaluation);
                return double.IsNegativeInfinity(density) ? double.NegativeInfinity : density + proposal;
            }
            bool accepted = RandomWalkMetropolis.Step(ref u, ref logTarget, Math.Exp(_logStep[g]), LogDensity, rng);
            Attempts++;
            _windowAttempts[g]++;
            if (accepted)
            {
                Hand(state, candidate!, evaluation!);
                Accepts++;
                _windowAccepts[g]++;
                any = true;
            }
            if (Adapting && _windowAttempts[g] == Window) Adapt(g);
        }
        return any;
    }

    // Each window of 50 moves the log step by ±min(0.1, 1/√windows), up when the window accepted above the target rate.
    private void Adapt(int g)
    {
        if (g == 0) _windows++;
        double delta = Math.Min(0.1, 1.0 / Math.Sqrt(Math.Max(1, _windows)));
        double rate = _windowAccepts[g] / (double)_windowAttempts[g];
        _logStep[g] += rate > TargetRate ? delta : -delta;
        _windowAttempts[g] = 0;
        _windowAccepts[g] = 0;
    }
}

/// <summary>
/// The categorical parameters (an extension beyond R, D09). While the level combinations number at most the cap, their
/// joint conditional at the chain's rung is enumerated whole and drawn from (Francom et al. 2019): every combination's
/// prediction in one call (<see cref="ICalibrationEmulator.PredictLevels"/>), its likelihood through the state's factor,
/// the uniform prior cancelling. Above the cap each categorical parameter in turn takes a Metropolis step to a uniformly
/// chosen other level, whose mixing is measured, not assumed.
/// </summary>
public sealed class CategoricalBlock : CalibrationBlock
{
    private readonly int[][]? _combinations;

    /// <param name="target">The target.</param>
    /// <param name="cap">The most combinations the block enumerates.</param>
    public CategoricalBlock(CalibrationTarget target, int cap) : base(target)
    {
        ParameterLayout layout = target.Model.Layout;
        if (layout.Categorical.Count == 0) throw new ArgumentException("The model has no categorical parameter.", nameof(target));
        if (cap < 1) throw new ArgumentOutOfRangeException(nameof(cap), "The cap is at least 1.");
        if (layout.Combinations > cap) return;
        _combinations = new int[layout.Combinations][];
        for (long c = 0; c < layout.Combinations; c++)
        {
            _combinations[c] = new int[layout.Categorical.Count];
            layout.Levels(c, _combinations[c]);
        }
    }

    /// <inheritdoc/>
    public override string Name => Enumerated ? "levels (enumerated)" : "levels (single site)";

    /// <summary>Whether the joint conditional is enumerated (the combinations within the cap).</summary>
    public bool Enumerated => _combinations is not null;

    /// <summary>Single-site attempts and acceptances (above the cap), the mixing's first measure.</summary>
    public long SiteAttempts { get; private set; }

    /// <inheritdoc cref="SiteAttempts"/>
    public long SiteAccepts { get; private set; }

    /// <inheritdoc/>
    public override bool Update(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng)
    {
        double b = Rung(state);
        CalibrationState current = state.Current;
        if (_combinations is null) return SingleSite(state, rng);
        (double[] probabilities, double[,] predictions) = Enumerate(current, b);
        double u = rng.NextDouble(), cumulative = 0.0;
        int chosen = probabilities.Length - 1;
        for (int c = 0; c < probabilities.Length; c++)
        {
            cumulative += probabilities[c];
            if (u < cumulative) { chosen = c; break; }
        }
        var prediction = new double[Target.Model.Outputs];
        for (int j = 0; j < prediction.Length; j++) prediction[j] = predictions[chosen, j];
        CalibrationState next = current.WithLevels(_combinations[chosen], Target, prediction);
        Hand(state, next, Target.EvaluateState(next));
        return true;
    }

    /// <summary>The enumerated conditional at <paramref name="state"/> and rung <paramref name="b"/>: each combination's
    /// probability, in <see cref="ParameterLayout.Levels"/>' order.</summary>
    public double[] Conditional(CalibrationState state, double b) => Enumerate(state, b).Probabilities;

    private (double[] Probabilities, double[,] Predictions) Enumerate(CalibrationState state, double b)
    {
        if (_combinations is null) throw new InvalidOperationException("The combinations exceed the cap: the conditional is not enumerated.");
        CalibrationModel model = Target.Model;
        int m = model.Outputs, count = _combinations.Length;
        var predictions = new double[count, m];
        model.Emulator.PredictLevels(Target.Draw, model.Physical(state), _combinations, predictions);
        CovarianceFactor factor = Target.Factor(state);
        double[] delta = model.Delta(state);
        ReadOnlySpan<double> y = model.Observations;
        var logWeight = new double[count];
        var r = new double[m];
        double max = double.NegativeInfinity;
        for (int c = 0; c < count; c++)
        {
            for (int j = 0; j < m; j++)
            {
                if (!double.IsFinite(predictions[c, j])) throw new InvalidOperationException("The emulator predicted a value that is not finite.");
                r[j] = y[j] - predictions[c, j] - delta[j];
            }
            logWeight[c] = -0.5 * b * factor.InverseQuadraticForm(r);
            max = Math.Max(max, logWeight[c]);
        }
        double total = 0.0;
        for (int c = 0; c < count; c++) total += logWeight[c] = Math.Exp(logWeight[c] - max);
        for (int c = 0; c < count; c++) logWeight[c] /= total;
        return (logWeight, predictions);
    }

    private bool SingleSite(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng)
    {
        ParameterLayout layout = Target.Model.Layout;
        bool any = false;
        for (int j = 0; j < layout.Categorical.Count; j++)
        {
            CalibrationState current = state.Current;
            int levels = layout.Categorical[j].Levels.Length;
            int[] next = current.Levels.ToArray();
            int other = rng.NextInt(levels - 1);
            next[j] = other >= next[j] ? other + 1 : other;
            CalibrationState candidate = current.WithLevels(next);
            CalibrationEvaluation evaluation = Target.EvaluateState(candidate);
            double logAccept = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity;
            SiteAttempts++;
            if (Math.Log(1.0 - rng.NextDouble()) < logAccept)
            {
                Hand(state, candidate, evaluation);
                SiteAccepts++;
                any = true;
            }
        }
        return any;
    }
}
