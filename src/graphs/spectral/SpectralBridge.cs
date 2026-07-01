using System;
using Graphs.Primitives;

namespace Graphs.Spectral
{
    /// <summary>
    /// Builds a unit line field over graph nodes from a scalar eigenvector (or Fiedler vector)
    /// and ambient positions — graph-spectral analysis, not Riemannian manifold geometry.
    /// </summary>
    public static class SpectralBridge
    {
        public static double[] LineFieldFromEigenvector(
            CsrGraph graph,
            double[][] positions,
            double[] eigenvector)
        {
            if (positions is null) throw new ArgumentNullException(nameof(positions));
            if (eigenvector is null) throw new ArgumentNullException(nameof(eigenvector));
            if (positions.Length != graph.NodeCount)
                throw new ArgumentException("positions length must match graph.NodeCount.", nameof(positions));
            if (eigenvector.Length != graph.NodeCount)
                throw new ArgumentException("eigenvector length must match graph.NodeCount.", nameof(eigenvector));
            if (positions.Length == 0)
                return Array.Empty<double>();

            int ambientDim = positions[0].Length;
            for (int i = 1; i < positions.Length; i++)
            {
                if (positions[i].Length != ambientDim)
                    throw new ArgumentException("All position rows must share the same ambient dimension.", nameof(positions));
            }

            int n = graph.NodeCount;
            var lineField = new double[n * ambientDim];

            for (int i = 0; i < n; i++)
            {
                int rowStart = graph.RowPointers[i];
                int rowEnd = graph.RowPointers[i + 1];
                double normSq = 0.0;

                for (int edge = rowStart; edge < rowEnd; edge++)
                {
                    int j = graph.Targets[edge];
                    double weight = graph.Weights[edge];
                    double dphi = eigenvector[j] - eigenvector[i];
                    int offset = i * ambientDim;

                    for (int dim = 0; dim < ambientDim; dim++)
                    {
                        double contribution = weight * dphi * (positions[j][dim] - positions[i][dim]);
                        lineField[offset + dim] += contribution;
                    }
                }

                int lineOffset = i * ambientDim;
                for (int dim = 0; dim < ambientDim; dim++)
                {
                    double value = lineField[lineOffset + dim];
                    normSq += value * value;
                }

                if (normSq < 1e-20)
                {
                    for (int dim = 0; dim < ambientDim; dim++)
                        lineField[lineOffset + dim] = 0.0;
                    continue;
                }

                double invNorm = 1.0 / Math.Sqrt(normSq);
                for (int dim = 0; dim < ambientDim; dim++)
                    lineField[lineOffset + dim] *= invNorm;
            }

            return lineField;
        }

        public static double[] LineFieldFromFiedler(
            CsrGraph graph,
            double[][] positions,
            int seed = 0,
            LaplacianType lapType = LaplacianType.Combinatorial,
            SolverKind solverKind = SolverKind.Auto)
        {
            var pairs = Spectral.ComputeBottomK(
                graph,
                seed: seed,
                k: 2,
                lapType: lapType,
                solverKind: solverKind);

            if (pairs.Count < 2)
            {
                throw new InvalidOperationException(
                    $"Spectral.ComputeBottomK returned {pairs.Count} eigenpairs; need at least 2 for Fiedler extraction.");
            }

            return LineFieldFromEigenvector(graph, positions, pairs[1].Vector);
        }
    }
}
