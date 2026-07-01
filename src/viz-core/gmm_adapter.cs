using System;
using Clustering.Statistical.GMM;

namespace Viz.Adapters.Gmm
{
    public static class GmmVizAdapter
    {
        /// <summary>
        /// Flattens a fitted GMM into a GaussianLayer. If a merge strategy is supplied,
        /// the layer's ComponentToClusterMap is populated so the renderer can group
        /// components into clusters.
        /// </summary>
        public static GaussianLayer ToGaussianLayer(
            GaussianMixtureModel gmm,
            string name = "Fitted GMM",
            IComponentMergeStrategy? mergeStrategy = null)
        {
            int k = gmm.K;
            int d = gmm.Dimension;
            var means = new double[k * d];
            var covariances = new double[k * d * d];
            var weights = new double[k];

            for (int ki = 0; ki < k; ki++)
            {
                var comp = gmm.Components[ki];
                weights[ki] = comp.Weight;
                for (int dim = 0; dim < d; dim++)
                    means[ki * d + dim] = comp.Mean[dim];
                for (int row = 0; row < d; row++)
                    for (int col = 0; col < d; col++)
                        covariances[ki * d * d + row * d + col] = comp.Covariance[row, col];
            }

            int[]? componentToClusterMap = mergeStrategy?.Merge(
                gmm.Components, gmm.GetFinalResponsibilities());

            return new GaussianLayer(name, means, covariances, weights, k, d, componentToClusterMap);
        }

        /// <summary>Per-point hard component assignment as a LabelLayer.</summary>
        public static LabelLayer ToComponentLabels(
            GaussianMixtureModel gmm, double[][] data,
            string name = "GMM Component")
        {
            return new LabelLayer(name, gmm.Predict(data), LabelLayerKind.GmmComponent);
        }

        /// <summary>
        /// Per-point cluster assignment via a merge strategy as a LabelLayer.
        /// Component predictions are remapped through the strategy's component-to-cluster map.
        /// </summary>
        public static LabelLayer ToClusterLabels(
            GaussianMixtureModel gmm, double[][] data,
            IComponentMergeStrategy strategy,
            string name = "GMM Cluster")
        {
            int[] componentLabels = gmm.Predict(data);
            int[] map = strategy.Merge(gmm.Components, gmm.GetFinalResponsibilities());
            int n = componentLabels.Length;
            var clusterLabels = new int[n];
            for (int i = 0; i < n; i++)
                clusterLabels[i] = map[componentLabels[i]];
            return new LabelLayer(name, clusterLabels, LabelLayerKind.GmmCluster);
        }

        /// <summary>
        /// Slices column k from the N×K responsibility matrix cached by the last Fit call.
        /// Throws if Fit has not been called.
        /// </summary>
        public static NodeSignalLayer ToResponsibilityScalar(
            GaussianMixtureModel gmm, int k, string? name = null)
        {
            double[,] r = gmm.GetFinalResponsibilities()
                ?? throw new InvalidOperationException(
                    "GMM has no cached responsibilities; call Fit first.");
            int n = r.GetLength(0);
            var values = new double[n];
            for (int i = 0; i < n; i++)
                values[i] = r[i, k];
            return new NodeSignalLayer(
                name ?? $"Responsibility (k={k})", values, ScalarSource.Responsibility);
        }

        /// <summary>
        /// Per-point minimum squared Mahalanobis distance over all components,
        /// square-rooted to linearize the scale.
        /// </summary>
        public static NodeSignalLayer ToMahalanobisScalar(
            GaussianMixtureModel gmm, double[][] data,
            string name = "Mahalanobis (min over k)")
        {
            double[,] mahal = gmm.Mahal(data);
            int n = data.Length;
            int kCount = gmm.K;
            var values = new double[n];
            for (int i = 0; i < n; i++)
            {
                double minSq = double.MaxValue;
                for (int ki = 0; ki < kCount; ki++)
                {
                    double v = mahal[i, ki];
                    if (v < minSq) minSq = v;
                }
                values[i] = Math.Sqrt(minSq);
            }
            return new NodeSignalLayer(name, values, ScalarSource.MahalanobisDistance);
        }
    }
}
