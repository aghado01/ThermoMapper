using System.Collections.Generic;
using Graphs.Primitives;

namespace Graphs.Observables
{
    public readonly record struct EdgeWeightSummary(
        int EdgeCount,
        double MinWeight,
        double MaxWeight,
        double MedianWeight,
        double MeanWeight,
        int NearZeroBridges);

    public static class EdgeWeights
    {
        public static EdgeWeightSummary Summary(CsrGraph graph, double nearZeroEpsilon = 1e-10)
        {
            var weights = new List<double>();
            double minWeight = double.PositiveInfinity;
            double maxWeight = double.NegativeInfinity;
            double weightSum = 0.0;
            int nearZeroCount = 0;

            for (int source = 0; source < graph.NodeCount; source++)
            {
                int rowStart = graph.RowPointers[source];
                int rowEnd = graph.RowPointers[source + 1];
                for (int edge = rowStart; edge < rowEnd; edge++)
                {
                    int target = graph.Targets[edge];
                    if (target <= source)
                        continue;

                    double weight = graph.Weights[edge];
                    weights.Add(weight);
                    if (weight < minWeight) minWeight = weight;
                    if (weight > maxWeight) maxWeight = weight;
                    if (weight <= nearZeroEpsilon) nearZeroCount++;
                    weightSum += weight;
                }
            }

            if (weights.Count == 0)
                return new EdgeWeightSummary(0, 0.0, 0.0, 0.0, 0.0, 0);

            weights.Sort();
            double median = Statistics.MedianOfSorted(System.Runtime.InteropServices.CollectionsMarshal.AsSpan(weights));

            return new EdgeWeightSummary(
                EdgeCount: weights.Count,
                MinWeight: minWeight,
                MaxWeight: maxWeight,
                MedianWeight: median,
                MeanWeight: weightSum / weights.Count,
                NearZeroBridges: nearZeroCount);
        }
    }
}
