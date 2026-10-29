using System;
using System.Collections.Generic;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Samplers.Evidence;

/// <summary>Which estimator turns a ladder's per-rung draws into a log normalizer ratio.</summary>
public enum LadderEstimator
{
    /// <summary>
    /// Stepping-stone (Xie, Lewis, Fan, Kuo and Chen 2011): the sum over adjacent rungs of <c>ln Ẑ_{k+1}/Z_k</c>, each
    /// ratio importance-sampled from rung k's draws of <c>q_{k+1}/q_k</c>. Unbiased for each ratio at any spacing, and it
    /// reads ℓ_b at the rungs only, which every target supplies. The default.
    /// </summary>
    SteppingStone,

    /// <summary>
    /// Path sampling, or thermodynamic integration (Gelman and Meng 1998; Friel and Pettitt 2008): the trapezoid over b
    /// of each rung's mean slope <c>∂_b ln q_b</c>, with Friel, Hurn and Wyse's (2014) correction as its discretization
    /// diagnostic. On the likelihood path it needs the target's slope.
    /// </summary>
    PathSampling,
}

/// <summary>
/// One rung's draws of a scalar for a ladder estimate: ordered chains from MCMC, whose effective sample size comes from
/// their autocorrelation (<see cref="ChainDiagnostics.MeanEss"/>), or a weighted population from SMC, whose effective
/// sample size comes from its weights (<see cref="SmcWeights.EffectiveSampleSize"/>).
/// </summary>
public sealed class RungDraws
{
    private readonly double[][]? _chains;
    private readonly double[]? _values, _weights;   // a population; the weights normalized

    private RungDraws(double[][]? chains, double[]? values, double[]? weights)
    {
        _chains = chains;
        _values = values;
        _weights = weights;
    }

    /// <summary>Ordered draws, one array per chain, all of one length; the rung's chains across a ladder ensemble.</summary>
    public static RungDraws FromChains(params double[][] chains)
    {
        ArgumentNullException.ThrowIfNull(chains);
        if (chains.Length < 1) throw new ArgumentException("At least one chain is required.", nameof(chains));
        int n = chains[0]?.Length ?? 0;
        if (n < 1) throw new ArgumentException("A chain needs draws.", nameof(chains));
        var copy = new double[chains.Length][];
        for (int c = 0; c < chains.Length; c++)
        {
            ArgumentNullException.ThrowIfNull(chains[c]);
            if (chains[c].Length != n) throw new ArgumentException("Every chain must hold the same number of draws.", nameof(chains));
            copy[c] = (double[])chains[c].Clone();
        }
        return new RungDraws(copy, null, null);
    }

    /// <summary>
    /// A population with its weights (null: equal), as an SMC stage leaves it; draws of zero weight are dropped.
    /// Resampled offspring share ancestors, which weights do not see, so its effective sample size is an upper bound.
    /// </summary>
    public static RungDraws FromPopulation(double[] values, double[]? weights = null)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (values.Length < 1) throw new ArgumentException("A population needs draws.", nameof(values));
        if (weights is not null && weights.Length != values.Length)
            throw new ArgumentException("One weight per draw is required.", nameof(weights));
        var w = new double[values.Length];
        double total = 0.0;
        for (int i = 0; i < w.Length; i++)
        {
            w[i] = weights is null ? 1.0 : weights[i];
            if (!(w[i] >= 0.0) || double.IsPositiveInfinity(w[i]))
                throw new ArgumentException("Weights must be finite and non-negative.", nameof(weights));
            total += w[i];
        }
        if (!(total > 0.0)) throw new ArgumentException("The weights carry no mass.", nameof(weights));
        var kept = new List<double>(w.Length);
        var keptWeights = new List<double>(w.Length);
        for (int i = 0; i < w.Length; i++)
            if (w[i] > 0.0) { kept.Add(values[i]); keptWeights.Add(w[i] / total); }
        return new RungDraws(null, kept.ToArray(), keptWeights.ToArray());
    }

    /// <summary>The draws held.</summary>
    public int Count => _chains is null ? _values!.Length : _chains.Length * _chains[0].Length;

    // Chains of one shape: the per-draw contributions of aligned rungs can be summed draw by draw.
    internal bool SameShape(RungDraws other)
        => _chains is not null && other._chains is not null && _chains.Length == other._chains.Length
           && _chains[0].Length == other._chains[0].Length;

    internal double[][] Chains => _chains!;

    // The largest draw, refusing NaN and +∞ (−∞ is a zero ratio).
    internal double Max(string rung)
    {
        double max = double.NegativeInfinity;
        foreach (double v in All())
        {
            if (double.IsNaN(v) || double.IsPositiveInfinity(v)) throw new ArgumentException($"{rung} has a draw of {v}.");
            if (v > max) max = v;
        }
        return max;
    }

    internal void RequireFinite(string rung)
    {
        foreach (double v in All())
            if (!double.IsFinite(v)) throw new ArgumentException($"{rung} has a draw of {v}.");
    }

    // The mean of map(draw), the draws' variance about it, and the variance of the mean: the variance over the effective
    // sample size, which is zero for constant draws and NaN when the ESS is unavailable.
    internal (double Mean, double Variance, double VarianceOfMean, double Ess) Summarize(Func<double, double> map)
    {
        if (_chains is null)
        {
            double mean = 0.0, variance = 0.0;
            for (int i = 0; i < _values!.Length; i++) mean += _weights![i] * map(_values[i]);
            for (int i = 0; i < _values.Length; i++)
            {
                double d = map(_values[i]) - mean;
                variance += _weights![i] * d * d;
            }
            double ess = SmcWeights.EffectiveSampleSize(_weights);
            return (mean, variance, variance / ess, ess);
        }

        var mapped = new double[_chains.Length][];
        double sum = 0.0;
        int count = 0;
        for (int c = 0; c < _chains.Length; c++)
        {
            mapped[c] = Array.ConvertAll(_chains[c], v => map(v));
            foreach (double v in mapped[c]) sum += v;
            count += mapped[c].Length;
        }
        double m = sum / count, ss = 0.0;
        foreach (double[] chain in mapped)
            foreach (double v in chain) ss += (v - m) * (v - m);
        double var = count > 1 ? ss / (count - 1) : 0.0;
        double effective = Effective(mapped);
        return (m, var, var == 0.0 ? 0.0 : var / effective, effective);
    }

    // The ESS of a mean over chains: estimated, +∞ for constant draws (no error), NaN when unavailable.
    internal static double Effective(double[][] chains)
    {
        EssDiagnostic ess = ChainDiagnostics.MeanEss(chains);
        return ess.State switch
        {
            DiagnosticState.Estimated => ess.Value,
            DiagnosticState.Constant => double.PositiveInfinity,
            _ => double.NaN,
        };
    }

    private IEnumerable<double> All()
    {
        if (_chains is null)
        {
            foreach (double v in _values!) yield return v;
            yield break;
        }
        foreach (double[] chain in _chains)
            foreach (double v in chain) yield return v;
    }
}

/// <summary>One rung's part of a ladder estimate: its b, the value it contributes (the mean slope for path sampling,
/// <c>ln Ẑ_{k+1}/Z_k</c> to the next rung up for stepping-stone), that value's Monte Carlo error and its draws' effective
/// sample size.</summary>
public readonly record struct RungEstimate(double InverseTemperature, double Value, double MonteCarloError, double EffectiveSampleSize);

/// <summary>
/// A ladder's estimate of <c>ln Z(b_max) − ln Z(b_min)</c>: an evidence when the lowest rung is a normalized reference
/// (b = 0 under a proper prior) and the target's densities are complete.
/// </summary>
/// <param name="Method">The estimator.</param>
/// <param name="LogRatio">The estimate.</param>
/// <param name="MonteCarloError">Its Monte Carlo standard error.</param>
/// <param name="DiscretizationError">Path sampling: Friel, Hurn and Wyse's correction, <c>−Σ Δ_k²/12·(V_{k+1} − V_k)</c>,
/// an estimate of the truth minus the trapezoid. It reads each rung's slope variance V_k as the integrand's slope in b,
/// which holds where ln q_b is linear in b (the powered path, a separable likelihood); a collapsed likelihood adds
/// <c>E_b[∂²ℓ_b/∂b²]</c>, which it does not see. Stepping-stone has none: NaN.</param>
/// <param name="Rungs">Per rung, in ascending b: path sampling's every rung, stepping-stone's every rung but the top.</param>
public sealed record LadderEstimate(
    LadderEstimator Method, double LogRatio, double MonteCarloError, double DiscretizationError, IReadOnlyList<RungEstimate> Rungs);

/// <summary>
/// Independent chains' ladder estimates of one log normalizer ratio pooled into one (<see cref="LadderEvidence.Pool"/>),
/// each chain on its own ladder. Each chain's error reads its own draws only, so where its rungs barely communicate it
/// under-covers, and the chains' spread shows it; the reported error is the larger of the two readings.
/// </summary>
/// <param name="Method">The chains' estimator.</param>
/// <param name="LogRatio">The chains' mean, z̄.</param>
/// <param name="MonteCarloError">The reported error: the larger of <paramref name="WithinChainError"/> and
/// <paramref name="BetweenChainError"/>, or the within-chain error alone for one chain; NaN when either is.</param>
/// <param name="WithinChainError">The chains' own errors pooled, √Σ s_c² / C.</param>
/// <param name="BetweenChainError">The standard error of the chains' mean from their spread,
/// √(Σ(ẑ_c − z̄)² / (C(C − 1))); NaN for one chain.</param>
/// <param name="DiscretizationError">The chains' mean discretization diagnostic; NaN for stepping-stone.</param>
/// <param name="Chains">Each chain's estimate.</param>
public sealed record PooledLadderEstimate(
    LadderEstimator Method, double LogRatio, double MonteCarloError, double WithinChainError, double BetweenChainError,
    double DiscretizationError, IReadOnlyList<LadderEstimate> Chains)
{
    /// <summary>The between-chain error over the within-chain one: well above 1 where the chains' own errors
    /// under-cover; NaN for one chain.</summary>
    public double BetweenToWithin => BetweenChainError / WithinChainError;
}

/// <summary>
/// Evidence over a tempering ladder: <c>ln Z(b_max) − ln Z(b_min)</c> from per-rung draws, by stepping-stone or path
/// sampling (<see cref="LadderEstimator"/>). It knows nothing about models: it reduces the ladder's per-rung scalar
/// draws, which a traced <see cref="ParallelTempering{TState}"/> records (<see cref="RungFunctional"/>) and an SMC run's
/// stage populations supply, to a log normalizer ratio with its Monte Carlo error.
/// <para>The Monte Carlo error is the delta method's, over effective sample sizes, read two ways, and the larger is
/// reported. Taken as independent, each rung's variance enters over its own effective sample size, which misses the
/// covariance swaps leave between rungs. When every rung's draws are chains of one shape, as a ladder's traces are, scan
/// t's contributions are also summed across rungs and the error of that sum's mean is read, which counts the covariance,
/// but one autocorrelation truncation over the sum can stop short of a slowly mixing rung's tail.</para>
/// <para>Independent chains, each on its own ladder, pool by <see cref="Pool"/>, whose error also reads the chains'
/// spread.</para>
/// </summary>
public static class LadderEvidence
{
    /// <summary>The estimator <see cref="FromLadders"/> uses unless told otherwise.</summary>
    public const LadderEstimator DefaultEstimator = LadderEstimator.SteppingStone;

    /// <summary>
    /// Stepping-stone from each rung's draws of <c>ln q_c(x) − ln q_b(x)</c>, x drawn at the rung's b and c the next
    /// rung's b toward b_max. The rung at b_max has no draws (its entry may be null). The ladder may be in either order.
    /// </summary>
    public static LadderEstimate SteppingStone(ReadOnlySpan<double> inverseTemperatures, IReadOnlyList<RungDraws?> logRatios)
    {
        ArgumentNullException.ThrowIfNull(logRatios);
        int[] order = Ascending(inverseTemperatures, logRatios.Count, nameof(logRatios));
        int intervals = order.Length - 1;
        var draws = new RungDraws[intervals];
        var maps = new Func<double, double>[intervals];
        var coefficients = new double[intervals];
        var variances = new double[intervals];
        var rungs = new RungEstimate[intervals];
        double logRatio = 0.0;
        for (int k = 0; k < intervals; k++)
        {
            double b = inverseTemperatures[order[k]];
            string name = $"Rung b = {b:R}";
            draws[k] = logRatios[order[k]] ?? throw new ArgumentException($"{name} needs its draws of the log ratio to the next rung.", nameof(logRatios));
            double shift = draws[k].Max(name);
            if (double.IsNegativeInfinity(shift)) throw new ArgumentException($"{name}'s draws give the next rung zero density everywhere.", nameof(logRatios));
            maps[k] = d => Math.Exp(d - shift);
            (double mean, _, double varianceOfMean, double ess) = draws[k].Summarize(maps[k]);
            double rungLogRatio = shift + Math.Log(mean);
            coefficients[k] = 1.0 / mean;
            variances[k] = varianceOfMean / (mean * mean);
            rungs[k] = new RungEstimate(b, rungLogRatio, Math.Sqrt(variances[k]), ess);
            logRatio += rungLogRatio;
        }
        return new LadderEstimate(LadderEstimator.SteppingStone, logRatio, MonteCarloError(draws, maps, coefficients, variances),
                                  double.NaN, rungs);
    }

    /// <summary>
    /// Path sampling from each rung's draws of the slope <c>∂_b ln q_b(x)</c> at its own b: the trapezoid over the ladder,
    /// in either order, with the discretization diagnostic of <see cref="LadderEstimate.DiscretizationError"/>.
    /// </summary>
    public static LadderEstimate PathSampling(ReadOnlySpan<double> inverseTemperatures, IReadOnlyList<RungDraws> slopes)
    {
        ArgumentNullException.ThrowIfNull(slopes);
        int[] order = Ascending(inverseTemperatures, slopes.Count, nameof(slopes));
        int n = order.Length;
        var b = new double[n];
        var draws = new RungDraws[n];
        var maps = new Func<double, double>[n];
        var coefficients = new double[n];
        var variances = new double[n];
        var spread = new double[n];
        var rungs = new RungEstimate[n];
        for (int k = 0; k < n; k++) b[k] = inverseTemperatures[order[k]];
        double logRatio = 0.0;
        for (int k = 0; k < n; k++)
        {
            string name = $"Rung b = {b[k]:R}";
            draws[k] = slopes[order[k]] ?? throw new ArgumentException($"{name} needs its slope draws.", nameof(slopes));
            draws[k].RequireFinite(name);
            maps[k] = u => u;
            (double mean, double variance, double varianceOfMean, double ess) = draws[k].Summarize(maps[k]);
            double weight = 0.5 * ((k > 0 ? b[k] - b[k - 1] : 0.0) + (k < n - 1 ? b[k + 1] - b[k] : 0.0));
            coefficients[k] = weight;
            variances[k] = weight * weight * varianceOfMean;
            spread[k] = variance;
            rungs[k] = new RungEstimate(b[k], mean, Math.Sqrt(varianceOfMean), ess);
            logRatio += weight * mean;
        }
        double correction = 0.0;
        for (int k = 0; k < n - 1; k++)
        {
            double h = b[k + 1] - b[k];
            correction -= h * h / 12.0 * (spread[k + 1] - spread[k]);
        }
        return new LadderEstimate(LadderEstimator.PathSampling, logRatio, MonteCarloError(draws, maps, coefficients, variances),
                                  correction, rungs);
    }

    /// <summary>
    /// The estimate from the rung traces of one or more ladders run on the same rungs and path, each traced
    /// (<see cref="ParallelTempering{TState}.TraceRungs"/>) for the same number of scans; each ladder is one chain per
    /// rung. Path sampling needs traces that hold slopes.
    /// </summary>
    public static LadderEstimate FromLadders<TState>(IReadOnlyList<ParallelTempering<TState>> ladders,
                                                     LadderEstimator method = DefaultEstimator)
    {
        ArgumentNullException.ThrowIfNull(ladders);
        if (ladders.Count < 1) throw new ArgumentException("At least one ladder is required.", nameof(ladders));
        double[] betas = ladders[0].Betas;
        for (int l = 0; l < ladders.Count; l++)
        {
            ParallelTempering<TState> ladder = ladders[l];
            if (ladder.RungTraces is null) throw new InvalidOperationException($"Ladder {l} does not trace its rungs (TraceRungs).");
            if (ladder.Path != ladders[0].Path || !ladder.Betas.AsSpan().SequenceEqual(betas))
                throw new ArgumentException("Every ladder must run on the same rungs and path.", nameof(ladders));
            if (method == LadderEstimator.PathSampling && !ladder.TracesSlopes)
                throw new InvalidOperationException(
                    "The target supplies no slope in b (ITargetEvaluation.LogLikelihoodSlope); stepping-stone needs none.");
        }

        int functional = method == LadderEstimator.PathSampling ? RungFunctional.Slope : RungFunctional.LogRatio;
        var rungs = new RungDraws?[betas.Length];
        for (int i = 0; i < betas.Length; i++)
        {
            if (method == LadderEstimator.SteppingStone && i == 0) continue;   // the cold rung, b_max, has no ratio
            var chains = new double[ladders.Count][];
            for (int l = 0; l < ladders.Count; l++) chains[l] = ladders[l].RungTraces![i].Draws(functional);
            rungs[i] = RungDraws.FromChains(chains);
        }
        return method switch
        {
            LadderEstimator.SteppingStone => SteppingStone(betas, rungs),
            LadderEstimator.PathSampling => PathSampling(betas, rungs!),
            _ => throw new ArgumentOutOfRangeException(nameof(method)),
        };
    }

    /// <summary>
    /// Independent chains' estimates by one estimator, each over its own ladder, pooled (<see cref="PooledLadderEstimate"/>):
    /// their mean, the within-chain error √Σ s_c² / C, the between-chain error √(Σ(ẑ_c − z̄)² / (C(C − 1))) for C ≥ 2, the
    /// larger of the two reported (the within-chain error alone for one chain), and their mean discretization diagnostic.
    /// </summary>
    public static PooledLadderEstimate Pool(IReadOnlyList<LadderEstimate> chains)
    {
        ArgumentNullException.ThrowIfNull(chains);
        if (chains.Count < 1) throw new ArgumentException("At least one chain's estimate is required.", nameof(chains));
        LadderEstimator method = (chains[0] ?? throw new ArgumentException("A chain's estimate is missing.", nameof(chains))).Method;
        double sum = 0.0, variance = 0.0, discretization = 0.0;
        foreach (LadderEstimate chain in chains)
        {
            if (chain is null) throw new ArgumentException("A chain's estimate is missing.", nameof(chains));
            if (chain.Method != method) throw new ArgumentException("Every chain's estimate must come from one estimator.", nameof(chains));
            sum += chain.LogRatio;
            variance += chain.MonteCarloError * chain.MonteCarloError;
            discretization += chain.DiscretizationError;
        }
        int c = chains.Count;
        double mean = sum / c, within = Math.Sqrt(variance) / c, between = double.NaN;
        if (c >= 2)
        {
            double spread = 0.0;
            foreach (LadderEstimate chain in chains) spread += (chain.LogRatio - mean) * (chain.LogRatio - mean);
            between = Math.Sqrt(spread / (c * (c - 1.0)));
        }
        return new PooledLadderEstimate(method, mean, c >= 2 ? Math.Max(within, between) : within, within, between, discretization / c,
                                        [.. chains]);
    }

    /// <summary>
    /// The estimate from a tempered SMC run that traced its stages (<see cref="SmcOptions.TraceStages"/>), on one target
    /// and path: each stage is a rung, its population weighted as the stage left it. Stepping-stone reproduces the run's
    /// <see cref="SmcResult{TState}.LogNormalizerRatio"/>, and its error is the run's genealogy estimate
    /// (<see cref="SmcResult{TState}.LogNormalizerError"/>), which sees shared ancestry and is NaN unless every
    /// resampling was multinomial; its rungs' errors are their own weights'. Path sampling needs slopes, and its error is from the stages' weights, which ignore shared ancestry:
    /// the genealogy estimate covers the normalizer, not a sum of stage means.
    /// </summary>
    public static LadderEstimate FromSmc<TState>(SmcResult<TState> result, LadderEstimator method = DefaultEstimator)
    {
        ArgumentNullException.ThrowIfNull(result);
        IReadOnlyList<SmcStageTrace> stages = result.StageTraces
            ?? throw new InvalidOperationException("The run did not trace its stages (SmcOptions.TraceStages).");
        var b = new double[stages.Count];
        var rungs = new RungDraws?[stages.Count];
        for (int k = 0; k < stages.Count; k++)
        {
            SmcStageTrace stage = stages[k];
            if (stage.NewTarget || stage.Temperature.Path != stages[0].Temperature.Path)
                throw new ArgumentException("A ladder evidence needs a run on one target and path.", nameof(result));
            b[k] = stage.Temperature.InverseTemperature;
            if (method == LadderEstimator.PathSampling)
            {
                if (double.IsNaN(stage.Slopes[0]))
                    throw new InvalidOperationException(
                        "The target supplies no slope in b (ITargetEvaluation.LogLikelihoodSlope); stepping-stone needs none.");
                rungs[k] = RungDraws.FromPopulation(stage.Slopes, stage.Weights);
            }
            else if (k + 1 < stages.Count)
                rungs[k] = RungDraws.FromPopulation(stage.LogRatios!, stage.Weights);
        }
        return method switch
        {
            LadderEstimator.SteppingStone => SteppingStone(b, rungs) with { MonteCarloError = result.LogNormalizerError },
            LadderEstimator.PathSampling => PathSampling(b, rungs!),
            _ => throw new ArgumentOutOfRangeException(nameof(method)),
        };
    }

    // Rung indices in ascending b, from a ladder that is strictly ascending or strictly descending.
    private static int[] Ascending(ReadOnlySpan<double> inverseTemperatures, int draws, string drawsName)
    {
        ReadOnlySpan<double> b = inverseTemperatures;
        if (b.Length < 2) throw new ArgumentException("A ladder needs at least two rungs.", nameof(inverseTemperatures));
        if (draws != b.Length) throw new ArgumentException("One entry per rung is required.", drawsName);
        foreach (double v in b)
            if (!double.IsFinite(v)) throw new ArgumentException("Every rung's b must be finite.", nameof(inverseTemperatures));
        bool up = b[1] > b[0];
        for (int k = 1; k < b.Length; k++)
            if (!(up ? b[k] > b[k - 1] : b[k] < b[k - 1]))
                throw new ArgumentException("The ladder must be strictly monotone.", nameof(inverseTemperatures));
        var order = new int[b.Length];
        for (int k = 0; k < order.Length; k++) order[k] = up ? k : b.Length - 1 - k;
        return order;
    }

    // The rungs' own errors combined as independent; for aligned chains, the larger of that and the error of the mean of
    // each draw's contribution summed over rungs.
    private static double MonteCarloError(RungDraws[] draws, Func<double, double>[] maps, double[] coefficients, double[] variances)
    {
        double total = 0.0;
        foreach (double v in variances) total += v;
        double independent = Math.Sqrt(total);
        bool aligned = true;
        foreach (RungDraws rung in draws) aligned &= rung.SameShape(draws[0]);
        return aligned ? Math.Max(independent, Joint(draws, maps, coefficients)) : independent;
    }

    private static double Joint(RungDraws[] draws, Func<double, double>[] maps, double[] coefficients)
    {
        int chains = draws[0].Chains.Length, length = draws[0].Chains[0].Length;
        var sums = new double[chains][];
        double sum = 0.0;
        for (int c = 0; c < chains; c++)
        {
            sums[c] = new double[length];
            for (int k = 0; k < draws.Length; k++)
            {
                double[] chain = draws[k].Chains[c];
                for (int t = 0; t < length; t++) sums[c][t] += coefficients[k] * maps[k](chain[t]);
            }
            foreach (double v in sums[c]) sum += v;
        }
        int count = chains * length;
        double mean = sum / count, ss = 0.0;
        foreach (double[] chain in sums)
            foreach (double v in chain) ss += (v - mean) * (v - mean);
        double variance = count > 1 ? ss / (count - 1) : 0.0;
        return variance == 0.0 ? 0.0 : Math.Sqrt(variance / RungDraws.Effective(sums));
    }
}
