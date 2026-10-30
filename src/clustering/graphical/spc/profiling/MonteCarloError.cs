using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Maths.Samplers.Ensemble;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>
/// Monte Carlo standard errors of the sweep curves, from the per-draw record: sd/√ESS over the draws pooled across the
/// replicas at one temperature, with the effective sample size from <see cref="ChainDiagnostics.MeanEss"/> and the
/// replicas as the chains.
/// <list type="bullet">
/// <item><see cref="OfMean"/>, for a curve that is a mean of a per-draw quantity: χ_FK, so S₂ = χ_FK/N.</item>
/// <item><see cref="OfVariance"/>, for a curve that is a variance: Var(m), and C_v = Var(E)/T². A variance is the mean
/// of y = (x − x̄)², x̄ the pooled mean, so its error is y's sd/√ESS; the autocorrelation that matters is y's.</item>
/// </list>
/// </summary>
/// <remarks>
/// NaN where the error is undefined: no replicas, replicas of unequal length, fewer than four draws, replicas frozen at
/// different values, or an ESS the diagnostics cannot estimate. Zero where the quantity is constant. A variance's x̄ is
/// plugged in, as the curves' own reductions plug it in; the bias that leaves is O(1/n).
/// </remarks>
public static class MonteCarloError
{
    /// <summary>The error of the mean of x, from the replicas' per-draw series, which must be of equal length.</summary>
    public static double OfMean(IReadOnlyList<double[]> replicas)
    {
        if (Check(replicas) is double early)
            return early;
        double[][] chains = Chains(replicas);
        return Error(chains, Mean(chains));
    }

    /// <summary>The error of the variance of x, from the replicas' per-draw series, which must be of equal length.</summary>
    public static double OfVariance(IReadOnlyList<double[]> replicas)
    {
        // x's own spread first. Replicas each constant but not all equal never mixed, and frozen symmetric about x̄
        // they give a constant y that would read as no error; MeanEss reports only the constant case, so check here.
        if (Check(replicas) is double early)
            return early;

        double[][] chains = Chains(replicas);
        double mean = Mean(chains);
        var y = new double[chains.Length][];
        for (int r = 0; r < y.Length; r++)
        {
            double[] x = chains[r];
            var yr = new double[x.Length];
            for (int t = 0; t < x.Length; t++)
            {
                double d = x[t] - mean;
                yr[t] = d * d;
            }
            y[r] = yr;
        }
        return Error(y, Mean(y));
    }

    /// <summary>
    /// Group <paramref name="frames"/> by <see cref="Accumulator.Temperature"/> (ascending, as
    /// <see cref="SweepCurves.ByTemperature"/> groups them) and estimate an error at each temperature from the series
    /// <paramref name="seriesOf"/> reads off each frame. A frame whose series is null (no per-draw record) makes its
    /// temperature NaN.
    /// </summary>
    /// <param name="frames">The per-draw accumulator frames across the sweep.</param>
    /// <param name="seriesOf">The per-draw series a curve is reduced from, or null when the frame has none.</param>
    /// <param name="estimate">The error from the replicas' series at one temperature.</param>
    public static (double[] Temperatures, double[] Values) ByTemperature(
        IReadOnlyList<Accumulator> frames,
        Func<Accumulator, double[]?> seriesOf,
        Func<IReadOnlyList<double[]>, double, double> estimate)
    {
        ArgumentNullException.ThrowIfNull(frames);
        ArgumentNullException.ThrowIfNull(seriesOf);
        ArgumentNullException.ThrowIfNull(estimate);

        var grouped = new SortedDictionary<double, List<double[]>?>();
        foreach (Accumulator frame in frames)
        {
            double[]? series = seriesOf(frame);
            if (!grouped.TryGetValue(frame.Temperature, out List<double[]>? replicas))
                grouped[frame.Temperature] = series is null ? null : [series];
            else if (replicas is not null && series is not null)
                replicas.Add(series);
            else
                grouped[frame.Temperature] = null;
        }

        var temperatures = new double[grouped.Count];
        var values = new double[grouped.Count];
        int index = 0;
        foreach ((double temperature, List<double[]>? replicas) in grouped)
        {
            temperatures[index] = temperature;
            values[index] = replicas is null ? double.NaN : estimate(replicas, temperature);
            index++;
        }
        return (temperatures, values);
    }

    /// <summary>The error of C_v = Var(E)/T² at each temperature; NaN at T ≤ 0.</summary>
    public static (double[] Temperatures, double[] Values) SpecificHeat(IReadOnlyList<Accumulator> frames)
        => ByTemperature(frames, f => Recorded(f, f.EnergyPerDraw),
            (replicas, t) => t > 0.0 ? OfVariance(replicas) / (t * t) : double.NaN);

    /// <summary>The error of Var(m), the magnetization variance, at each temperature.</summary>
    public static (double[] Temperatures, double[] Values) MagnetizationVariance(IReadOnlyList<Accumulator> frames)
        => ByTemperature(frames, f => Recorded(f, f.MagPerDraw), (replicas, _) => OfVariance(replicas));

    /// <summary>
    /// The error of the FK susceptibility χ_FK = ⟨Σ_c |c|²⟩/N at each temperature, from the per-draw Σ|c|². The pair
    /// connectivity S₂ = χ_FK/N carries this error over N.
    /// </summary>
    public static (double[] Temperatures, double[] Values) SusceptibilityFkCluster(IReadOnlyList<Accumulator> frames)
        => ByTemperature(frames, PerSiteSumSqClusterSizes, (replicas, _) => OfMean(replicas));

    private static double[]? Recorded(Accumulator frame, double[] series)
        => series.Length == frame.DrawCount ? series : null;

    private static double[]? PerSiteSumSqClusterSizes(Accumulator frame)
    {
        long[] sums = frame.SumSqClusterSizesPerDraw;
        int n = frame.Spins.Length;
        if (sums.Length != frame.DrawCount || n <= 0)
            return null;
        var series = new double[sums.Length];
        for (int t = 0; t < sums.Length; t++)
            series[t] = (double)sums[t] / n;
        return series;
    }

    // NaN for no replicas, unequal lengths or fewer than four draws; 0 for one value everywhere; NaN for replicas each
    // constant and not all equal (as the diagnostics define Frozen); null to go on and estimate.
    private static double? Check(IReadOnlyList<double[]> replicas)
    {
        ArgumentNullException.ThrowIfNull(replicas);
        if (replicas.Count == 0)
            return double.NaN;
        int n = replicas[0]?.Length ?? 0;
        foreach (double[] series in replicas)
            if (series is null || series.Length != n)
                return double.NaN;
        if (n < 4)
            return double.NaN;

        double first = replicas[0][0];
        bool allEqual = true, eachConstant = true;
        foreach (double[] series in replicas)
            foreach (double x in series)
            {
                if (x != first) allEqual = false;
                if (x != series[0]) eachConstant = false;
            }
        if (allEqual) return 0.0;
        return eachConstant ? double.NaN : null;
    }

    private static double[][] Chains(IReadOnlyList<double[]> replicas)
    {
        if (replicas is double[][] chains)
            return chains;
        chains = new double[replicas.Count][];
        for (int r = 0; r < chains.Length; r++)
            chains[r] = replicas[r];
        return chains;
    }

    private static double Mean(double[][] chains)
    {
        double sum = 0.0;
        long count = 0;
        foreach (double[] chain in chains)
        {
            foreach (double x in chain)
                sum += x;
            count += chain.Length;
        }
        return sum / count;
    }

    // sd/√ESS over the pooled draws, the sd with an n − 1 denominator.
    private static double Error(double[][] chains, double mean)
    {
        double ss = 0.0;
        long count = 0;
        foreach (double[] chain in chains)
        {
            foreach (double x in chain)
                ss += (x - mean) * (x - mean);
            count += chain.Length;
        }
        double variance = ss / (count - 1);

        EssDiagnostic ess = ChainDiagnostics.MeanEss(chains);
        return ess.State switch
        {
            DiagnosticState.Estimated => Math.Sqrt(variance / ess.Value),
            DiagnosticState.Constant  => 0.0,
            _                         => double.NaN,
        };
    }
}
