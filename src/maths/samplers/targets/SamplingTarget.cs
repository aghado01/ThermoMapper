using System;

namespace Maths.Samplers;

/// <summary>
/// The model-neutral target every sampler consumes: reversible-jump chains, tempering ladders and, later, SMC. A state's
/// evaluation splits its log density into the untempered log prior and the log likelihood at any inverse temperature
/// b, with the state's expensive work (a design, a factorization) done once, so another temperature is a query rather
/// than a re-evaluation. That is what a likelihood-tempered swap or an SMC reweighting needs.
/// <para>Values are needed only up to an additive constant that may depend on b and the data but never on the state:
/// it cancels in Metropolis ratios, in swaps and in normalized SMC weights. An evidence estimate needs complete
/// values.</para>
/// </summary>
public interface ISamplingTarget<in TState>
{
    /// <summary>Evaluate <paramref name="state"/> once; temperatures are then queries on the result.</summary>
    ITargetEvaluation Evaluate(TState state);
}

/// <summary>One state's evaluation under an <see cref="ISamplingTarget{TState}"/>.</summary>
public interface ITargetEvaluation
{
    /// <summary>The untempered log prior; −∞ off the support, where no likelihood work is done and the log
    /// likelihood is −∞ too.</summary>
    double LogPrior { get; }

    /// <summary>
    /// ℓ_b, the log likelihood at inverse temperature <paramref name="inverseTemperature"/> in [0, 1]; b = 1 is the
    /// full likelihood. A separable model's is b·ℓ₁. A collapsed model integrates its nuisance parameters against the
    /// tempered likelihood, which is not linear in b. A model may refuse b ≠ 1 (<see cref="NotSupportedException"/>),
    /// or b = 0 when a nuisance prior is improper (<see cref="ArgumentOutOfRangeException"/>).
    /// </summary>
    double LogLikelihood(double inverseTemperature);
}

/// <summary>Which tempered density a replica samples.</summary>
public enum TemperingPath
{
    /// <summary>The powered posterior <c>(prior · likelihood)^b</c>: hot replicas flatten the prior as well.</summary>
    Powered,

    /// <summary>The likelihood-tempered path <c>prior · likelihood^b</c> (the power posterior of Friel and Pettitt
    /// 2008): the prior stays whole, and hot replicas approach it.</summary>
    Likelihood,
}

/// <summary>An inverse temperature and the path it is taken along. At b = 1 both paths are the posterior.</summary>
public readonly struct Temperature : IEquatable<Temperature>
{
    /// <summary>b in [0, 1].</summary>
    public double InverseTemperature { get; }

    /// <summary>The path.</summary>
    public TemperingPath Path { get; }

    public Temperature(double inverseTemperature, TemperingPath path)
    {
        if (!(inverseTemperature >= 0.0 && inverseTemperature <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature), "The inverse temperature must lie in [0, 1].");
        if (path is not (TemperingPath.Powered or TemperingPath.Likelihood))
            throw new ArgumentOutOfRangeException(nameof(path));
        InverseTemperature = inverseTemperature;
        Path = path;
    }

    /// <summary>b = 1: the posterior.</summary>
    public static Temperature Posterior => new(1.0, TemperingPath.Likelihood);

    /// <summary>
    /// The log density of an evaluated state at this temperature: <c>b·(ln prior + ℓ₁)</c> on the powered path,
    /// <c>ln prior + ℓ_b</c> on the likelihood path. A state off the support is −∞ without any likelihood work.
    /// </summary>
    public double LogDensity(ITargetEvaluation evaluation)
    {
        ArgumentNullException.ThrowIfNull(evaluation);
        double prior = evaluation.LogPrior;
        if (double.IsNegativeInfinity(prior)) return double.NegativeInfinity;
        return Path == TemperingPath.Powered
            ? InverseTemperature * (prior + evaluation.LogLikelihood(1.0))
            : prior + evaluation.LogLikelihood(InverseTemperature);
    }

    /// <summary>The untempered log posterior <c>ln prior + ℓ₁</c>.</summary>
    public static double LogPosterior(ITargetEvaluation evaluation) => Posterior.LogDensity(evaluation);

    public bool Equals(Temperature other) => InverseTemperature.Equals(other.InverseTemperature) && Path == other.Path;
    public override bool Equals(object? obj) => obj is Temperature t && Equals(t);
    public override int GetHashCode() => HashCode.Combine(InverseTemperature, Path);
    public static bool operator ==(Temperature a, Temperature b) => a.Equals(b);
    public static bool operator !=(Temperature a, Temperature b) => !a.Equals(b);
    public override string ToString() => $"b = {InverseTemperature} ({Path})";
}

/// <summary>Targets from functions, for samplers and tests that need no model of their own.</summary>
public static class SamplingTarget
{
    /// <summary>
    /// A target given by its log density alone, read as a likelihood under a flat prior, so b tempers all of it and
    /// the two paths agree. A state where the density is −∞ is off the support.
    /// </summary>
    public static ISamplingTarget<TState> FromLogDensity<TState>(Func<TState, double> logDensity)
    {
        ArgumentNullException.ThrowIfNull(logDensity);
        return new FunctionTarget<TState>(state =>
        {
            double d = logDensity(state);
            return new SeparableEvaluation(double.IsNegativeInfinity(d) ? double.NegativeInfinity : 0.0, d);
        });
    }

    /// <summary>A separable target: <c>ℓ_b = b·ℓ₁</c>, with the prior untempered.</summary>
    public static ISamplingTarget<TState> Separable<TState>(Func<TState, double> logPrior, Func<TState, double> logLikelihood)
    {
        ArgumentNullException.ThrowIfNull(logPrior);
        ArgumentNullException.ThrowIfNull(logLikelihood);
        return new FunctionTarget<TState>(state =>
        {
            double prior = logPrior(state);
            return new SeparableEvaluation(prior, double.IsNegativeInfinity(prior) ? double.NegativeInfinity : logLikelihood(state));
        });
    }

    /// <summary>The untempered log posterior of one state, <c>ln prior + ℓ₁</c>.</summary>
    public static double LogPosterior<TState>(this ISamplingTarget<TState> target, TState state)
    {
        ArgumentNullException.ThrowIfNull(target);
        return Temperature.LogPosterior(target.Evaluate(state));
    }

    private sealed class FunctionTarget<TState> : ISamplingTarget<TState>
    {
        private readonly Func<TState, ITargetEvaluation> _evaluate;
        public FunctionTarget(Func<TState, ITargetEvaluation> evaluate) => _evaluate = evaluate;
        public ITargetEvaluation Evaluate(TState state) => _evaluate(state);
    }
}

/// <summary>A separable evaluation: a log prior and an untempered log likelihood ℓ₁, with <c>ℓ_b = b·ℓ₁</c>.</summary>
public sealed class SeparableEvaluation : ITargetEvaluation
{
    private readonly double _logLikelihood;

    public SeparableEvaluation(double logPrior, double logLikelihood)
    {
        LogPrior = logPrior;
        _logLikelihood = logLikelihood;
    }

    public double LogPrior { get; }

    /// <summary><c>b·ℓ₁</c>; a zero likelihood stays −∞ at any b, so the support survives b = 0.</summary>
    public double LogLikelihood(double inverseTemperature)
        => double.IsNegativeInfinity(_logLikelihood) ? double.NegativeInfinity : inverseTemperature * _logLikelihood;
}
