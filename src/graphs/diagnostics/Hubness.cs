using System;
using Graphs;
using Graphs.Primitives;

namespace Graphs.Diagnostics
{
    public readonly record struct HubnessReport(
        int K,
        int NodeCount,
        int MaxInDegree,
        double MeanInDegree,
        double InDegreeSkewness,
        int HubCount,
        int AntiHubCount,
        double TopHubCoverage);

    public static class Hubness
    {
        /// <summary>
        /// Analyze hub-spoke pathology in the directed KNN result. Pass the
        /// directed NeighborSelection returned by SelectKnn, not a mutual-KNN result.
        /// </summary>
        public static HubnessReport Analyze(
            NeighborSelection directedKnn,
            int k,
            double hubMultiple = 10.0,
            double antiHubMultiple = 0.1)
        {
            int n = directedKnn.AllNeighbors.Length;
            if (n == 0)
                return new HubnessReport(k, 0, 0, 0.0, double.NaN, 0, 0, 0.0);

            var inDegree = new int[n];
            for (int source = 0; source < n; source++)
            {
                Neighbor[] row = directedKnn.AllNeighbors[source];
                for (int i = 0; i < row.Length; i++)
                    inDegree[row[i].Index]++;
            }

            int maxInDegree = 0;
            long totalInDegree = 0;
            int hubThreshold = (int)Math.Ceiling(hubMultiple * k);
            int antiHubThreshold = (int)Math.Floor(antiHubMultiple * k);
            int hubCount = 0;
            int antiHubCount = 0;

            for (int i = 0; i < n; i++)
            {
                int degree = inDegree[i];
                if (degree > maxInDegree) maxInDegree = degree;
                totalInDegree += degree;
                if (degree >= hubThreshold) hubCount++;
                if (degree <= antiHubThreshold) antiHubCount++;
            }

            double meanInDegree = (double)totalInDegree / n;
            double skewness    = Statistics.Skewness(inDegree.AsSpan());

            int topN = Math.Max(1, (int)Math.Sqrt(n));
            var sorted = (int[])inDegree.Clone();
            Array.Sort(sorted);

            long topMass = 0;
            for (int i = sorted.Length - topN; i < sorted.Length; i++)
                topMass += sorted[i];

            double topHubCoverage = totalInDegree > 0 ? (double)topMass / totalInDegree : 0.0;

            return new HubnessReport(
                K: k,
                NodeCount: n,
                MaxInDegree: maxInDegree,
                MeanInDegree: meanInDegree,
                InDegreeSkewness: skewness,
                HubCount: hubCount,
                AntiHubCount: antiHubCount,
                TopHubCoverage: topHubCoverage);
        }
    }
}
