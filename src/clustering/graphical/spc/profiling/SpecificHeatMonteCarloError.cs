using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Maths.Samplers.Ensemble;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>
/// The Monte Carlo standard error of the specific heat C_v at each temperature, from the per-draw energy record.
/// C_v is the mean of y = (E − Ē)²/T² over the draws pooled across the replicas, Ē being their pooled mean, so its
/// error is sd(y)/√ESS(y), with the effective sample size from <see cref="ChainDiagnostics.MeanEss"/> and the
/// replicas as the chains. The autocorrelation that matters is y's, not E's: C_v's error is the error of a variance.
/// </summary>
/// <remarks>
/// NaN where the error is undefined: T ≤ 0, a frame without its per-draw record, replicas of unequal length, fewer than
/// four draws, replicas frozen at different energies, or an ESS the diagnostics cannot estimate. Zero where E or y is
/// constant. Ē is plugged in, as C_v's own reduction plugs it in; the bias that leaves is O(1/n).
/// </remarks>
public static class SpecificHeatMonteCarloError
{
    /// <summary>
    /// Group <paramref name="frames"/> by <see cref="Accumulator.Temperature"/> (ascending, as
    /// <see cref="SweepCurves.ByTemperature"/> groups them) and return C_v's error at each temperature.
    /// </summary>
    public static (double[] Temperatures, double[] Values) ByTemperature(IReadOnlyList<Accumulator> frames)
    {
        ArgumentNullException.ThrowIfNull(frames);

        var grouped = new SortedDictionary<double, List<double[]>?>();
        foreach (Accumulator frame in frames)
        {
            bool recorded = frame.EnergyPerDraw.Length == frame.DrawCount;
            if (!grouped.TryGetValue(frame.Temperature, out List<double[]>? series))
                grouped[frame.Temperature] = recorded ? [frame.EnergyPerDraw] : null;
            else if (series is not null && recorded)
                series.Add(frame.EnergyPerDraw);
            else
                grouped[frame.Temperature] = null;
        }

        var temperatures = new double[grouped.Count];
        var values = new double[grouped.Count];
        int index = 0;
        foreach ((double temperature, List<double[]>? series) in grouped)
        {
            temperatures[index] = temperature;
            values[index] = series is null ? double.NaN : Of(series, temperature);
            index++;
        }
        return (temperatures, values);
    }

    /// <summary>
    /// C_v's error at one temperature from the replicas' per-draw energy series, which must be of equal length.
    /// </summary>
    /// <param name="energyPerReplica">Each replica's E for every retained draw.</param>
    /// <param name="temperature">T &gt; 0.</param>
    public static double Of(IReadOnlyList<double[]> energyPerReplica, double temperature)
    {
        ArgumentNullException.ThrowIfNull(energyPerReplica);
        if (!(temperature > 0.0) || energyPerReplica.Count == 0)
            return double.NaN;

        int n = energyPerReplica[0].Length;
        foreach (double[] series in energyPerReplica)
            if (series is null || series.Length != n)
                return double.NaN;
        if (n < 4)
            return double.NaN;

        // E's own spread first. Replicas each constant but not all equal never mixed, and frozen symmetric about Ē
        // they give a constant y that would read as no error; MeanEss reports only the constant case, so check here.
        switch (Spread(energyPerReplica))
        {
            case DiagnosticState.Constant: return 0.0;
            case DiagnosticState.Frozen:   return double.NaN;
        }

        double total = (double)energyPerReplica.Count * n;
        double mean = 0.0;
        foreach (double[] series in energyPerReplica)
            foreach (double e in series)
                mean += e;
        mean /= total;

        double invTSq = 1.0 / (temperature * temperature);
        var y = new double[energyPerReplica.Count][];
        double meanY = 0.0;
        for (int r = 0; r < y.Length; r++)
        {
            double[] series = energyPerReplica[r];
            var yr = new double[n];
            for (int t = 0; t < n; t++)
            {
                double d = series[t] - mean;
                yr[t] = d * d * invTSq;
                meanY += yr[t];
            }
            y[r] = yr;
        }
        meanY /= total;

        double varY = 0.0;
        foreach (double[] yr in y)
            foreach (double v in yr)
                varY += (v - meanY) * (v - meanY);
        varY /= total - 1.0;

        EssDiagnostic ess = ChainDiagnostics.MeanEss(y);
        return ess.State switch
        {
            DiagnosticState.Estimated => Math.Sqrt(varY / ess.Value),
            DiagnosticState.Constant  => 0.0,
            _                         => double.NaN,
        };
    }

    // As the diagnostics define it: Constant is one value everywhere; Frozen, each replica constant and not all equal.
    private static DiagnosticState Spread(IReadOnlyList<double[]> replicas)
    {
        double first = replicas[0][0];
        bool allEqual = true, eachConstant = true;
        foreach (double[] series in replicas)
            foreach (double e in series)
            {
                if (e != first) allEqual = false;
                if (e != series[0]) eachConstant = false;
            }
        if (allEqual) return DiagnosticState.Constant;
        return eachConstant ? DiagnosticState.Frozen : DiagnosticState.Estimated;
    }
}
