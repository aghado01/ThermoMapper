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
}
