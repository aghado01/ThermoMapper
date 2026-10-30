using System;
using Maths.Samplers.Evidence;
using Maths.Samplers.Smc;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A multicanonical run's draws as a weighted sample at one b: each draw's weight <c>e^{b·ℓ₁}/w(ℓ₁)</c>, self-normalized
/// over every chain's draws, chain-major in the draws' order, with its Kish effective sample size, pooled and per chain.
/// </summary>
/// <param name="InverseTemperature">b.</param>
/// <param name="Weights">Per draw, chain-major, summing to 1.</param>
/// <param name="EffectiveSampleSize">The pooled weights' Kish ESS, (Σw)²/Σw².</param>
/// <param name="ChainEffectiveSampleSizes">Each chain's own weights' Kish ESS.</param>
public sealed record WeightedDraws(double InverseTemperature, double[] Weights, double EffectiveSampleSize, double[] ChainEffectiveSampleSizes);

/// <summary>
/// A reweighted expectation at one b: the pooled self-normalized estimate, its Monte Carlo error from the production
/// chains' spread, <c>sd(chain estimates)/√C</c> (NaN for one chain), and each chain's own estimate.
/// </summary>
public sealed record ReweightedEstimate(double InverseTemperature, double Mean, double MonteCarloError, double[] ChainMeans);

/// <summary>
/// The evidence curve ln Z(b) = ln E_prior[e^{b·ℓ₁}] over a grid of b from one multicanonical run, with no rungs and
/// no path discretization: the ratio estimator <c>Z(b) = E_w[e^{b·ℓ₁}/w(ℓ₁)] / E_w[1/w(ℓ₁)]</c>, valid when the prior is
/// normalized. Each b's value pools every chain's draws, and its Monte Carlo error is the chains' spread,
/// <c>sd(ln Ẑ_c(b))/√C</c> (NaN for one chain). The record beside <see cref="LadderEstimate"/>'s ladders.
/// </summary>
/// <param name="InverseTemperatures">The grid.</param>
/// <param name="LogEvidence">ln Ẑ(b), pooled.</param>
/// <param name="MonteCarloError">Its error from the chains' spread.</param>
/// <param name="ChainLogEvidence">Per chain, ln Ẑ_c at each b.</param>
/// <param name="EffectiveSampleSize">At each b, the Kish ESS of the pooled weights e^{b·ℓ₁}/w.</param>
/// <param name="NormalizerEffectiveSampleSize">The Kish ESS of the denominator's weights 1/w: low when the window's low
/// extension leaves the prior's tail thinly sampled.</param>
public sealed record LogEvidenceCurve(
    double[] InverseTemperatures, double[] LogEvidence, double[] MonteCarloError, double[][] ChainLogEvidence,
    double[] EffectiveSampleSize, double NormalizerEffectiveSampleSize)
{
    /// <summary>ln Ẑ(b_to) − ln Ẑ(b_from), with its error from the chains' spread of the same difference.</summary>
    public (double LogRatio, double MonteCarloError) LogRatio(int from, int to)
    {
        int chains = ChainLogEvidence.Length;
        var differences = new double[chains];
        for (int c = 0; c < chains; c++) differences[c] = ChainLogEvidence[c][to] - ChainLogEvidence[c][from];
        return (LogEvidence[to] - LogEvidence[from], Reweighting.SpreadError(differences));
    }
}

/// <summary>
/// The reweighting of a multicanonical run's production (<see cref="MulticanonicalDraws"/>) to the likelihood-tempered
/// posterior at any b: draw s weighs <c>e^{b·ℓ₁(x_s)}/w(ℓ₁(x_s))</c>, self-normalized, since π_w ∝ prior·w(ℓ₁) and the
/// tempered posterior ∝ prior·e^{b·ℓ₁}. It reduces one run's draws to the weighted laws at every b and the ln Z(b)
/// curve. Any b in [0, 1] is defined; inside [b_lo, 1] the weights stay bounded where π_w's histogram is flat, and
/// below b_lo the low extension's draws weigh more the further they lie, which the Kish ESS shows.
/// </summary>
public static class Reweighting
{
    /// <summary>ln(e^{b·ℓ₁}/w(ℓ₁)) for one draw: a zero likelihood weighs nothing at b &gt; 0 and 1/w at b = 0.</summary>
    public static double LogImportance(double level, double logWeight, double inverseTemperature)
    {
        double tilt = double.IsNegativeInfinity(level) ? (inverseTemperature > 0.0 ? double.NegativeInfinity : 0.0) : inverseTemperature * level;
        return tilt - logWeight;
    }

    /// <summary>The draws as a weighted sample at <paramref name="inverseTemperature"/>.</summary>
    public static WeightedDraws Weigh(MulticanonicalDraws draws, double inverseTemperature)
    {
        ArgumentNullException.ThrowIfNull(draws);
        CheckTemperature(inverseTemperature);
        int n = draws.Draws, chains = draws.Chains;
        var logs = new double[chains * n];
        var chainEss = new double[chains];
        var scratch = new double[n];
        for (int c = 0; c < chains; c++)
        {
            Span<double> part = logs.AsSpan(c * n, n);
            LogImportances(draws, c, inverseTemperature, part);
            SmcWeights.Normalize(part, scratch);
            chainEss[c] = SmcWeights.EffectiveSampleSize(scratch);
        }
        var weights = new double[logs.Length];
        SmcWeights.Normalize(logs, weights);
        return new WeightedDraws(inverseTemperature, weights, SmcWeights.EffectiveSampleSize(weights), chainEss);
    }

    /// <summary>E_b[functional <paramref name="functional"/>], reweighted.</summary>
    public static ReweightedEstimate Expectation(MulticanonicalDraws draws, double inverseTemperature, int functional)
    {
        ArgumentNullException.ThrowIfNull(draws);
        if ((uint)functional >= (uint)draws.FunctionalDim) throw new ArgumentOutOfRangeException(nameof(functional));
        return Expectation(draws, inverseTemperature, (_, values) => values[functional]);
    }

    /// <summary>E_b[f(ℓ₁, functionals)], reweighted: f of each draw's level and functionals.</summary>
    public static ReweightedEstimate Expectation(MulticanonicalDraws draws, double inverseTemperature, DrawFunction f)
    {
        ArgumentNullException.ThrowIfNull(draws);
        ArgumentNullException.ThrowIfNull(f);
        CheckTemperature(inverseTemperature);
        int n = draws.Draws, chains = draws.Chains;
        var logs = new double[n];
        var weights = new double[n];
        var chainMeans = new double[chains];
        double pooledNumerator = 0.0, pooledDenominator = 0.0, pooledMax = double.NegativeInfinity;
        var chainNumerators = new double[chains];
        var chainDenominators = new double[chains];
        var chainMax = new double[chains];
        for (int c = 0; c < chains; c++)
        {
            LogImportances(draws, c, inverseTemperature, logs);
            double max = double.NegativeInfinity;
            foreach (double l in logs) max = Math.Max(max, l);
            if (double.IsNegativeInfinity(max)) throw new InvalidOperationException($"Every draw of chain {c} weighs nothing at b = {inverseTemperature}.");
            double numerator = 0.0, denominator = 0.0;
            ReadOnlySpan<double> levels = draws.Levels(c);
            for (int s = 0; s < n; s++)
            {
                double w = Math.Exp(logs[s] - max);
                if (w == 0.0) continue;
                numerator += w * f(levels[s], draws.Functionals(c, s));
                denominator += w;
            }
            chainMeans[c] = numerator / denominator;
            (chainNumerators[c], chainDenominators[c], chainMax[c]) = (numerator, denominator, max);
            pooledMax = Math.Max(pooledMax, max);
        }
        for (int c = 0; c < chains; c++)
        {
            double scale = Math.Exp(chainMax[c] - pooledMax);
            pooledNumerator += scale * chainNumerators[c];
            pooledDenominator += scale * chainDenominators[c];
        }
        return new ReweightedEstimate(inverseTemperature, pooledNumerator / pooledDenominator, SpreadError(chainMeans), chainMeans);
    }

    /// <summary>
    /// ln Z(b) over <paramref name="grid"/> by the ratio estimator, which needs the inner target's prior normalized:
    /// otherwise the curve is off by its log mass, a constant in b.
    /// </summary>
    public static LogEvidenceCurve LogEvidence(MulticanonicalDraws draws, ReadOnlySpan<double> grid)
    {
        ArgumentNullException.ThrowIfNull(draws);
        if (grid.Length == 0) throw new ArgumentException("The grid needs at least one b.", nameof(grid));
        foreach (double b in grid) CheckTemperature(b);
        int n = draws.Draws, chains = draws.Chains;
        var logs = new double[chains * n];

        // The denominator, E_w[1/w]: its log, pooled and per chain.
        var chainNormalizer = new double[chains];
        for (int c = 0; c < chains; c++)
        {
            LogImportances(draws, c, 0.0, logs.AsSpan(c * n, n));
            chainNormalizer[c] = LogMean(logs.AsSpan(c * n, n));
        }
        double normalizer = LogMean(logs);
        var normalized = new double[logs.Length];
        SmcWeights.Normalize(logs, normalized);
        double normalizerEss = SmcWeights.EffectiveSampleSize(normalized);

        var value = new double[grid.Length];
        var error = new double[grid.Length];
        var ess = new double[grid.Length];
        var perChain = new double[chains][];
        for (int c = 0; c < chains; c++) perChain[c] = new double[grid.Length];
        var column = new double[chains];
        for (int g = 0; g < grid.Length; g++)
        {
            for (int c = 0; c < chains; c++)
            {
                LogImportances(draws, c, grid[g], logs.AsSpan(c * n, n));
                perChain[c][g] = column[c] = LogMean(logs.AsSpan(c * n, n)) - chainNormalizer[c];
            }
            value[g] = LogMean(logs) - normalizer;
            error[g] = SpreadError(column);
            SmcWeights.Normalize(logs, normalized);
            ess[g] = SmcWeights.EffectiveSampleSize(normalized);
        }
        return new LogEvidenceCurve(grid.ToArray(), value, error, perChain, ess, normalizerEss);
    }

    // sd of the chains' estimates over √C; NaN for one chain.
    internal static double SpreadError(ReadOnlySpan<double> estimates)
    {
        int c = estimates.Length;
        if (c < 2) return double.NaN;
        double mean = 0.0;
        foreach (double e in estimates) mean += e;
        mean /= c;
        double ss = 0.0;
        foreach (double e in estimates) ss += (e - mean) * (e - mean);
        return Math.Sqrt(ss / (c - 1) / c);
    }

    private static void LogImportances(MulticanonicalDraws draws, int chain, double inverseTemperature, Span<double> logs)
    {
        ReadOnlySpan<double> levels = draws.Levels(chain);
        MulticanonicalWeights w = draws.Weights;
        for (int s = 0; s < levels.Length; s++) logs[s] = LogImportance(levels[s], w.LogWeight(levels[s]), inverseTemperature);
    }

    // ln of the mean of e^{logs}.
    private static double LogMean(ReadOnlySpan<double> logs)
    {
        double max = double.NegativeInfinity;
        foreach (double l in logs) max = Math.Max(max, l);
        if (double.IsNegativeInfinity(max)) return max;
        double s = 0.0;
        foreach (double l in logs) s += Math.Exp(l - max);
        return max + Math.Log(s / logs.Length);
    }

    private static void CheckTemperature(double inverseTemperature)
    {
        if (!(inverseTemperature >= 0.0 && inverseTemperature <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature), "The inverse temperature must lie in [0, 1].");
    }
}

/// <summary>A function of one production draw: its untempered log likelihood and its functionals.</summary>
public delegate double DrawFunction(double level, ReadOnlySpan<double> functionals);
