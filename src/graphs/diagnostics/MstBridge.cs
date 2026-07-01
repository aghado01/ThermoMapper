using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace Graphs.Diagnostics
{
    public readonly record struct MstBridgeReport(
        int BridgeCount,
        double MinBridgeWeight,
        double MaxBridgeWeight,
        double MedianBridgeWeight,
        double BridgeWeightSkewness);

    public static class MstBridge
    {
        public static MstBridgeReport Compare(CsrGraph withoutMst, CsrGraph withMst)
        {
            if (withoutMst.NodeCount != withMst.NodeCount)
                throw new ArgumentException("Graphs must share the same NodeCount.", nameof(withMst));

            var baselineEdges = new HashSet<long>();
            for (int source = 0; source < withoutMst.NodeCount; source++)
            {
                int rowStart = withoutMst.RowPointers[source];
                int rowEnd = withoutMst.RowPointers[source + 1];
                for (int edge = rowStart; edge < rowEnd; edge++)
                {
                    int target = withoutMst.Targets[edge];
                    if (target <= source)
                        continue;

                    baselineEdges.Add((((long)source) << 32) | (uint)target);
                }
            }

            var bridgeWeights = new List<double>();
            double minWeight = double.PositiveInfinity;
            double maxWeight = double.NegativeInfinity;

            for (int source = 0; source < withMst.NodeCount; source++)
            {
                int rowStart = withMst.RowPointers[source];
                int rowEnd = withMst.RowPointers[source + 1];
                for (int edge = rowStart; edge < rowEnd; edge++)
                {
                    int target = withMst.Targets[edge];
                    if (target <= source)
                        continue;

                    long key = (((long)source) << 32) | (uint)target;
                    if (baselineEdges.Contains(key))
                        continue;

                    double weight = withMst.Weights[edge];
                    bridgeWeights.Add(weight);
                    if (weight < minWeight) minWeight = weight;
                    if (weight > maxWeight) maxWeight = weight;
                }
            }

            if (bridgeWeights.Count == 0)
                return new MstBridgeReport(0, 0.0, 0.0, 0.0, double.NaN);

            bridgeWeights.Sort();
            var sortedSpan = System.Runtime.InteropServices.CollectionsMarshal.AsSpan(bridgeWeights);
            double median   = Statistics.MedianOfSorted(sortedSpan);
            double skewness = Statistics.Skewness(sortedSpan);

            return new MstBridgeReport(
                BridgeCount: bridgeWeights.Count,
                MinBridgeWeight: minWeight,
                MaxBridgeWeight: maxWeight,
                MedianBridgeWeight: median,
                BridgeWeightSkewness: skewness);
        }
    }
}
