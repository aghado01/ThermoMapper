using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Observables
{
    public readonly record struct AlgebraicConnectivityReport(
        double Lambda2,
        bool LikelyWeaklyConnected,
        int NodeCount,
        bool Computed);

    public static class AlgebraicConnectivity
    {
        public static AlgebraicConnectivityReport Compute(
            CsrGraph graph,
            int maxNodesForDense = 2000,
            double weakThreshold = 1e-6)
        {
            int n = graph.NodeCount;
            if (n == 0)
                return new AlgebraicConnectivityReport(0.0, false, 0, false);
            if (n == 1)
                return new AlgebraicConnectivityReport(0.0, false, 1, true);
            if (n > maxNodesForDense)
                return new AlgebraicConnectivityReport(double.NaN, false, n, false);

            var laplacian = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                int rowStart = graph.RowPointers[i];
                int rowEnd = graph.RowPointers[i + 1];
                double weightedDegree = 0.0;

                for (int edge = rowStart; edge < rowEnd; edge++)
                {
                    int j = graph.Targets[edge];
                    double w = graph.Weights[edge];
                    laplacian[i, j] = -w;
                    weightedDegree += w;
                }

                laplacian[i, i] = weightedDegree;
            }

            EigenResult eigen = DenseEigen.DecomposeSymmetric(laplacian);
            double lambda2 = eigen.Eigenvalues[n - 2];

            return new AlgebraicConnectivityReport(
                Lambda2: lambda2,
                LikelyWeaklyConnected: lambda2 < weakThreshold,
                NodeCount: n,
                Computed: true);
        }
    }
}
