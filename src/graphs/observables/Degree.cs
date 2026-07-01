using Graphs.Primitives;

namespace Graphs.Observables
{
    public readonly record struct DegreeReport(
        int NodeCount,
        int MinDegree,
        int MaxDegree,
        double MeanDegree,
        int IsolatedCount,
        int UndersampledCount);

    public static class Degree
    {
        public static DegreeReport Distribution(CsrGraph graph)
        {
            int n = graph.NodeCount;
            if (n == 0)
                return new DegreeReport(0, 0, 0, 0.0, 0, 0);

            int minDegree = int.MaxValue;
            int maxDegree = 0;
            long degreeSum = 0;
            int isolatedCount = 0;
            int undersampledCount = 0;

            for (int node = 0; node < n; node++)
            {
                int degree = graph.Degree(node);
                if (degree < minDegree) minDegree = degree;
                if (degree > maxDegree) maxDegree = degree;
                degreeSum += degree;
                if (degree == 0) isolatedCount++;
                if (degree == 1) undersampledCount++;
            }

            return new DegreeReport(
                NodeCount: n,
                MinDegree: minDegree,
                MaxDegree: maxDegree,
                MeanDegree: (double)degreeSum / n,
                IsolatedCount: isolatedCount,
                UndersampledCount: undersampledCount);
        }
    }
}
