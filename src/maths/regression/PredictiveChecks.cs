using System;
using System.Collections.Generic;
using Maths.Samplers.Ensemble;

namespace Maths.Regression;

/// <summary>An estimate of the expected log pointwise predictive density over units, with its parts per unit.</summary>
/// <param name="Elpd">Σ over units of each unit's estimate.</param>
/// <param name="StandardError">√(U·Var(per-unit estimates)), the usual standard error of the total.</param>
/// <param name="EffectiveParameters">The estimated effective number of parameters: p_WAIC, or lppd − elpd_loo.</param>
/// <param name="PerUnit">Each unit's estimate.</param>
/// <param name="ParetoK">PSIS-LOO only: each unit's Pareto shape k̂; above
/// <see cref="PredictiveChecks.ParetoKThreshold"/> that unit's estimate is unreliable (Vehtari et al. 2024). Null for
/// WAIC.</param>
/// <param name="MonteCarloError">PSIS-LOO only: the Monte Carlo standard error of <see cref="Elpd"/>, √Σ of the units'
/// squared errors; NaN when any unit's k̂ is above the threshold, as in loo, and for WAIC.</param>
/// <param name="MonteCarloErrors">PSIS-LOO only: each unit's Monte Carlo standard error. Null for WAIC.</param>
/// <param name="RelativeEfficiency">PSIS-LOO only: each unit's relative efficiency r_eff = ESS/S, 1 for independent
/// draws. Null for WAIC.</param>
public sealed record PredictiveCheck(double Elpd, double StandardError, double EffectiveParameters, double[] PerUnit,
                                     double[]? ParetoK, double MonteCarloError = double.NaN,
                                     double[]? MonteCarloErrors = null, double[]? RelativeEfficiency = null);

/// <summary>
/// Model-free consumers of a <see cref="PointwiseLogDensity"/> (OQ17): they read the draws × units matrix and know no
/// model. <see cref="Waic"/> is Watanabe's (2010) widely applicable information criterion; <see cref="PsisLoo"/> is
/// leave-one-unit-out cross-validation by Pareto-smoothed importance sampling (Vehtari, Gelman &amp; Gabry 2017; Vehtari,
/// Simpson, Gelman, Yao &amp; Gabry 2024), following the loo package's algorithm. WAIC treats the draws as independent.
/// PSIS-LOO reads the density's chains, when it has them, for each unit's relative efficiency, as loo does.
/// </summary>
public static class PredictiveChecks
{
    /// <summary>
    /// WAIC: per unit, lppd_u = log mean_s p(y_u | θ_s) and p_u = Var_s log p(y_u | θ_s), and elpd_u = lppd_u − p_u.
    /// </summary>
    public static PredictiveCheck Waic(PointwiseLogDensity density)
    {
        ArgumentNullException.ThrowIfNull(density);
        if (density.LeavesUnitsOut)
            throw new ArgumentException(
                "WAIC's terms must condition on each unit's own data; these are leave-unit-out predictives (the collapsed content). Use PsisLoo.",
                nameof(density));
        double[,] ll = density.Values;
        int s = ll.GetLength(0), units = ll.GetLength(1);
        if (s < 2) throw new ArgumentException("WAIC needs at least two draws.", nameof(density));
        var elpd = new double[units];
        double pTotal = 0.0;
        var column = new double[s];
        for (int u = 0; u < units; u++)
        {
            for (int d = 0; d < s; d++) column[d] = ll[d, u];
            double lppd = LogSumExp(column) - Math.Log(s);
            double mean = 0.0;
            foreach (double v in column) mean += v;
            mean /= s;
            double variance = 0.0;
            foreach (double v in column) variance += (v - mean) * (v - mean);
            variance /= s - 1;
            elpd[u] = lppd - variance;
            pTotal += variance;
        }
        return Summarize(elpd, pTotal, null);
    }

    /// <summary>
    /// PSIS-LOO: per unit, the importance log-ratios −log p(y_u | θ_s) are Pareto-smoothed in their upper tail (the largest
    /// M = ⌈min(S/5, 3√(S/r_eff))⌉, at least 5 for a fit), truncated at their raw maximum, and normalized to w̃_s; elpd_u
    /// is the log of the weighted mean of p(y_u | θ_s). k̂ is the tail's fitted shape.
    /// <para>
    /// The relative efficiency r_eff is the unit's ESS over S, of its likelihoods p(y_u | θ_s) within the density's chains:
    /// Stan's estimator on split chains (<see cref="ChainDiagnostics.MeanEss"/>, loo's <c>relative_eff</c>). Independent
    /// draws, or a unit whose ESS cannot be estimated, take 1; chains of unequal length are each cut to the shortest for
    /// it. A slowly mixing unit gets a longer tail. Its Monte Carlo standard error is loo's: the self-normalized
    /// estimator's relative variance Σ_s w̃_s²(p_s/p̂ − 1)²/r_eff (Vehtari et al. 2024, eq. 6), moved to the log scale
    /// by a log-normal's moments, √log(1 + that).
    /// </para>
    /// <para>
    /// With leave-unit-out terms (the collapsed content), the same arithmetic gives the non-factorized model's PSIS-LOO
    /// (Bürkner, Gabry &amp; Vehtari 2021): the ratios are 1/p(y_u | y_rest, structure), and elpd_u estimates log
    /// p(y_u | y_rest). The effective number of parameters, lppd − elpd, then has no in-sample lppd to subtract from and is
    /// NaN. Those ratios vary only with the chain's slow part and clump, which makes k̂ pessimistic (see
    /// <see cref="PointwiseContent.CollapsedPredictive"/>). A tail whose largest ratios tie has k̂ = −∞: bounded in the
    /// sample, it is left unsmoothed.
    /// </para>
    /// </summary>
    public static PredictiveCheck PsisLoo(PointwiseLogDensity density)
    {
        ArgumentNullException.ThrowIfNull(density);
        double[,] ll = density.Values;
        int s = ll.GetLength(0), units = ll.GetLength(1);
        if (s < 2) throw new ArgumentException("PSIS-LOO needs at least two draws.", nameof(density));
        var elpd = new double[units];
        var khat = new double[units];
        var mcse = new double[units];
        var rEff = new double[units];
        double lppdTotal = 0.0;
        var logRatios = new double[s];
        var column = new double[s];
        var weighted = new double[s];
        for (int u = 0; u < units; u++)
        {
            for (int d = 0; d < s; d++) { column[d] = ll[d, u]; logRatios[d] = -ll[d, u]; }
            lppdTotal += LogSumExp(column) - Math.Log(s);
            rEff[u] = RelativeEfficiency(column, density.ChainLengths);
            int tail = (int)Math.Ceiling(Math.Min(0.2 * s, 3.0 * Math.Sqrt(s / rEff[u])));
            khat[u] = SmoothLogWeights(logRatios, tail);   // logRatios now holds the smoothed log weights
            double logNorm = LogSumExp(logRatios);
            for (int d = 0; d < s; d++)
            {
                logRatios[d] -= logNorm;                   // normalized: log w̃_s
                weighted[d] = logRatios[d] + column[d];
            }
            double e = LogSumExp(weighted);
            elpd[u] = e;
            double relativeVariance = 0.0;
            for (int d = 0; d < s; d++)
            {
                double t = Math.Exp(logRatios[d]) * double.ExpM1(column[d] - e);   // w̃_s(p_s/p̂ − 1), within [−1, 1]: no overflow
                relativeVariance += t * t;
            }
            mcse[u] = double.IsNegativeInfinity(e) ? double.NaN : Math.Sqrt(double.LogP1(relativeVariance / rEff[u]));
        }
        double total = 0.0, mcVariance = 0.0, threshold = ParetoKThreshold(s);
        bool reliable = true;
        for (int u = 0; u < units; u++)
        {
            total += elpd[u];
            mcVariance += mcse[u] * mcse[u];
            reliable &= !(khat[u] > threshold);
        }
        PredictiveCheck check = Summarize(elpd, density.LeavesUnitsOut ? double.NaN : lppdTotal - total, khat);
        return check with
        {
            MonteCarloError = reliable ? Math.Sqrt(mcVariance) : double.NaN,
            MonteCarloErrors = mcse,
            RelativeEfficiency = rEff,
        };
    }

    /// <summary>
    /// The k̂ above which a unit's PSIS estimate is unreliable at S draws: min(1 − 1/log₁₀S, 0.7) (Vehtari et al. 2024;
    /// loo's <c>ps_khat_threshold</c>). It is 0.7 from S ≈ 2,154 draws and lower below, as fewer draws tolerate a
    /// lighter tail.
    /// </summary>
    public static double ParetoKThreshold(int draws)
    {
        if (draws < 2) throw new ArgumentOutOfRangeException(nameof(draws), "The threshold needs at least two draws.");
        return Math.Min(1.0 - 1.0 / Math.Log10(draws), 0.7);
    }

    // A unit's r_eff = ESS/S of its likelihoods within the chains: 1 for independent draws or an ESS that cannot be
    // estimated. The likelihoods are scaled by their largest, which the ESS ignores; unequal chains are cut to the shortest.
    private static double RelativeEfficiency(double[] logLik, IReadOnlyList<int>? chainLengths)
    {
        if (chainLengths is null) return 1.0;
        int m = chainLengths.Count, length = int.MaxValue;
        foreach (int c in chainLengths) length = Math.Min(length, c);
        double max = double.NegativeInfinity;
        foreach (double v in logLik) max = Math.Max(max, v);
        var chains = new double[m][];
        int start = 0;
        for (int j = 0; j < m; j++)
        {
            chains[j] = new double[length];
            for (int i = 0; i < length; i++) chains[j][i] = Math.Exp(logLik[start + i] - max);
            start += chainLengths[j];
        }
        EssDiagnostic ess = ChainDiagnostics.MeanEss(chains);
        return ess.State == DiagnosticState.Estimated ? ess.Value / (m * (double)length) : 1.0;
    }

    /// <summary>
    /// The generalized Pareto distribution fitted to exceedances by Zhang &amp; Stephens' (2009) empirical Bayes estimator,
    /// with the loo package's weakly informative adjustment of the shape toward ½: k̂ ← (n·k̂ + 5)/(n + 10). Returns the
    /// shape k (positive for a heavy tail) and scale σ. <paramref name="exceedances"/> must be sorted ascending and
    /// positive.
    /// </summary>
    public static (double Shape, double Scale) FitGeneralizedPareto(ReadOnlySpan<double> exceedances)
    {
        int n = exceedances.Length;
        if (n < 2) throw new ArgumentException("A fit needs at least two exceedances.", nameof(exceedances));
        const double prior = 3.0;
        int m = 30 + (int)Math.Floor(Math.Sqrt(n));
        double xStar = exceedances[(int)Math.Floor(n / 4.0 + 0.5) - 1];   // the first quartile
        double xMax = exceedances[n - 1];
        var theta = new double[m];
        var logLik = new double[m];
        for (int j = 0; j < m; j++)
        {
            theta[j] = 1.0 / xMax + (1.0 - Math.Sqrt(m / (j + 0.5))) / prior / xStar;
            double k = 0.0;
            foreach (double x in exceedances) k += Math.Log(1.0 - theta[j] * x);
            k /= n;
            logLik[j] = n * (Math.Log(-theta[j] / k) - k - 1.0);   // the profile log-likelihood
        }
        double norm = LogSumExp(logLik), thetaHat = 0.0;
        for (int j = 0; j < m; j++) thetaHat += theta[j] * Math.Exp(logLik[j] - norm);
        double shape = 0.0;
        foreach (double x in exceedances) shape += Math.Log(1.0 - thetaHat * x);
        shape /= n;
        double scale = -shape / thetaHat;
        shape = (n * shape + 10.0 * 0.5) / (n + 10.0);
        return (shape, scale);
    }

    // Pareto-smooth the upper tail of one unit's log ratios in place, then truncate at their raw maximum; returns k̂, or
    // +∞ when the tail is too short or flat to fit.
    private static double SmoothLogWeights(double[] logRatios, int tail)
    {
        int s = logRatios.Length;
        double max = double.NegativeInfinity;
        foreach (double v in logRatios) max = Math.Max(max, v);
        for (int d = 0; d < s; d++) logRatios[d] -= max;
        double khat = double.PositiveInfinity;
        if (tail >= 5 && tail < s)
        {
            var order = new int[s];
            for (int d = 0; d < s; d++) order[d] = d;
            var keys = (double[])logRatios.Clone();
            Array.Sort(keys, order);
            double cutoff = keys[s - tail - 1];   // the largest value below the tail
            // A tail flat to rounding cannot be fitted, and needs no smoothing: its largest ratios tie, so they are bounded in
            // the sample and no draw dominates. That is the norm for the collapsed content, whose ratios repeat with the
            // structure, and k̂ is then −∞. Too few draws for a tail leave k̂ at +∞, undiagnosable.
            if (!(Math.Abs(keys[s - 1] - keys[s - tail]) >= 2.220446049250313e-16 / 100.0)) khat = double.NegativeInfinity;
            else
            {
                double expCutoff = Math.Exp(cutoff);
                var exceedances = new double[tail];
                for (int t = 0; t < tail; t++) exceedances[t] = Math.Exp(keys[s - tail + t]) - expCutoff;
                (double shape, double scale) = FitGeneralizedPareto(exceedances);
                if (double.IsFinite(shape))
                {
                    for (int t = 0; t < tail; t++)
                    {
                        double p = (t + 0.5) / tail;
                        double q = Math.Abs(shape) < 1e-12
                            ? -scale * Math.Log(1.0 - p)
                            : scale * double.ExpM1(-shape * Math.Log(1.0 - p)) / shape;
                        logRatios[order[s - tail + t]] = Math.Log(q + expCutoff);
                    }
                    khat = shape;
                }
            }
        }
        for (int d = 0; d < s; d++) logRatios[d] = Math.Min(logRatios[d], 0.0) + max;   // truncate at the raw maximum
        return khat;
    }

    private static PredictiveCheck Summarize(double[] perUnit, double effectiveParameters, double[]? paretoK)
    {
        int units = perUnit.Length;
        double total = 0.0;
        foreach (double e in perUnit) total += e;
        double mean = total / units, ss = 0.0;
        foreach (double e in perUnit) ss += (e - mean) * (e - mean);
        double se = units > 1 ? Math.Sqrt(units * ss / (units - 1)) : double.NaN;
        return new PredictiveCheck(total, se, effectiveParameters, perUnit, paretoK);
    }

    internal static double LogSumExp(ReadOnlySpan<double> values)
    {
        double max = double.NegativeInfinity;
        foreach (double v in values) max = Math.Max(max, v);
        if (double.IsNegativeInfinity(max)) return max;
        double sum = 0.0;
        foreach (double v in values) sum += Math.Exp(v - max);
        return max + Math.Log(sum);
    }
}
