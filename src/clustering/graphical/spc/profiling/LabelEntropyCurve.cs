using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Observables;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>
/// The cluster-size entropy curve over the temperature axis, with its sweep-tier assembler:
/// <see cref="From"/> groups SW <see cref="Accumulator"/> frames by temperature, <b>pools</b> the
/// additive <c>ClusterSizeHistogram</c> across replicas at each T, then reduces <b>once</b> via
/// <see cref="ClusterSizeEntropy"/>.
/// </summary>
/// <remarks>
/// <b>Pool-then-reduce</b>, deliberately <i>not</i> <see cref="SweepCurves"/>' reduce-then-average:
/// entropy is nonlinear, so averaging per-replica entropies ≠ entropy of the pooled distribution
/// (Jensen). The histograms are the commutative sufficient-statistic accumulated per draw; the
/// entropy reduction (the <see cref="ClusterSizeEntropy"/> def in <c>graphs/observables</c>) is
/// applied at the end, on the pooled counts. SW keeps only the accumulation.
/// </remarks>
public readonly record struct LabelEntropyCurve(
    double[] Temperatures,
    double[] Entropy)
{
    public static LabelEntropyCurve From(IReadOnlyList<Accumulator> frames)
    {
        if (frames.Count == 0)
            return new LabelEntropyCurve(Array.Empty<double>(), Array.Empty<double>());

        var grouped = new SortedDictionary<double, int[]>();
        for (int i = 0; i < frames.Count; i++)
        {
            Accumulator frame = frames[i];
            if (!grouped.TryGetValue(frame.Temperature, out int[]? histogram))
            {
                histogram = new int[frame.ClusterSizeHistogram.Length];
                grouped[frame.Temperature] = histogram;
            }
            else if (histogram.Length < frame.ClusterSizeHistogram.Length)
            {
                Array.Resize(ref histogram, frame.ClusterSizeHistogram.Length);
                grouped[frame.Temperature] = histogram;
            }

            for (int j = 0; j < frame.ClusterSizeHistogram.Length; j++)
                histogram[j] += frame.ClusterSizeHistogram[j];
        }

        var temperatures = new double[grouped.Count];
        var entropy = new double[grouped.Count];
        int index = 0;
        foreach ((double temperature, int[] histogram) in grouped)
        {
            temperatures[index] = temperature;
            entropy[index] = ClusterSizeEntropy.EntropyNats(histogram);
            index++;
        }

        return new LabelEntropyCurve(temperatures, entropy);
    }
}
