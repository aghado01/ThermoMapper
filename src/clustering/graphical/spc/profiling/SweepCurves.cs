using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>
/// The sweep-tier half of the observable split: assembles a curve over the
/// temperature axis from a flat list of per-draw <see cref="Accumulator"/>
/// frames. The model reductions
/// (<see cref="Graphs.Models.Potts.Observables.Susceptibility"/> &amp; siblings)
/// say what each frame <i>means</i> as a Gibbs-measure quantity; this says how
/// those per-draw values assemble — group by temperature, average across the
/// replicas at each T. Sweep-native, inference-agnostic.
/// </summary>
public static class SweepCurves
{
    /// <summary>
    /// Group <paramref name="frames"/> by <see cref="Accumulator.Temperature"/>
    /// (ascending) and average <paramref name="perFrame"/> over the replicas at
    /// each temperature. Returns parallel <c>(temperatures, values)</c> arrays.
    /// </summary>
    /// <param name="frames">The per-draw accumulator frames across the sweep.</param>
    /// <param name="perFrame">
    /// The model reduction applied to each frame (e.g. a
    /// <c>Graphs.Models.Potts.Observables.*.Reduce</c> call). Multi-valued
    /// curves (magnetization moments) assemble each component with a separate
    /// call and combine the averaged components afterward — see
    /// <see cref="SweepProfile.From"/>.
    /// </param>
    public static (double[] Temperatures, double[] Values) ByTemperature(
        IReadOnlyList<Accumulator> frames,
        Func<Accumulator, double> perFrame)
    {
        if (frames.Count == 0)
            return (Array.Empty<double>(), Array.Empty<double>());

        var grouped = new SortedDictionary<double, (double Sum, int Count)>();
        for (int i = 0; i < frames.Count; i++)
        {
            Accumulator frame = frames[i];
            double value = perFrame(frame);

            if (grouped.TryGetValue(frame.Temperature, out (double Sum, int Count) acc))
                grouped[frame.Temperature] = (acc.Sum + value, acc.Count + 1);
            else
                grouped[frame.Temperature] = (value, 1);
        }

        var temperatures = new double[grouped.Count];
        var values = new double[grouped.Count];
        int index = 0;
        foreach ((double temperature, (double sum, int count)) in grouped)
        {
            temperatures[index] = temperature;
            values[index] = count > 0 ? sum / count : 0.0;
            index++;
        }

        return (temperatures, values);
    }

    /// <summary>
    /// <see cref="ByTemperature"/> with the replicas' spread kept beside each mean: group
    /// <paramref name="frames"/> by temperature (ascending) and return, at each, the replicas' mean of
    /// <paramref name="perFrame"/> (bit-identical to <see cref="ByTemperature"/>'s value), their sample
    /// standard deviation (n − 1 denominator; NaN for a single replica) and their count.
    /// </summary>
    /// <param name="frames">The per-draw accumulator frames across the sweep.</param>
    /// <param name="perFrame">The model reduction applied to each frame.</param>
    public static SweepCurve ByTemperatureWithSpread(
        IReadOnlyList<Accumulator> frames,
        Func<Accumulator, double> perFrame)
    {
        var grouped = new SortedDictionary<double, List<double>>();
        for (int i = 0; i < frames.Count; i++)
        {
            Accumulator frame = frames[i];
            if (!grouped.TryGetValue(frame.Temperature, out List<double>? values))
                grouped[frame.Temperature] = values = new List<double>();
            values.Add(perFrame(frame));
        }

        var temperatures = new double[grouped.Count];
        var means = new double[grouped.Count];
        var deviations = new double[grouped.Count];
        var replicas = new int[grouped.Count];
        int index = 0;
        foreach ((double temperature, List<double> values) in grouped)
        {
            // Summed from the first value in frame order, as ByTemperature sums, so the means agree bit for bit.
            double sum = values[0];
            for (int r = 1; r < values.Count; r++) sum += values[r];
            double mean = sum / values.Count;

            double squares = 0.0;
            foreach (double v in values) squares += (v - mean) * (v - mean);

            temperatures[index] = temperature;
            means[index] = mean;
            deviations[index] = values.Count > 1 ? Math.Sqrt(squares / (values.Count - 1)) : double.NaN;
            replicas[index] = values.Count;
            index++;
        }

        return new SweepCurve(temperatures, means, deviations, replicas);
    }
}

/// <summary>
/// A curve over the temperature axis with its spread across replicas: at each temperature, the
/// replicas' mean, their sample standard deviation (NaN for a single replica) and their count.
/// Parallel arrays, temperatures ascending.
/// </summary>
public sealed record SweepCurve(double[] Temperatures, double[] Means, double[] StandardDeviations, int[] Replicas);
