namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Maps each fitted component to a cluster index, producing the hierarchical
    /// output layer over a flat EM fit. See docs/gmm.md for the available
    /// strategies and when to use which.
    /// </summary>
    public interface IComponentMergeStrategy
    {
        /// <summary>
        /// Returns an <c>int[]</c> of length <c>components.Length</c> where entry k
        /// is the cluster index for component k. Cluster indices are dense in
        /// <c>[0, clusterCount)</c>; cluster-label ordering is strategy-defined.
        /// </summary>
        /// <param name="components">Fitted components with up-to-date caches.</param>
        /// <param name="responsibilities">
        /// N×K final responsibility matrix. Required by responsibility-based
        /// strategies (e.g. <see cref="EntropyMergeStrategy"/>); ignored by
        /// geometry-only strategies (e.g. <see cref="ModalMergeStrategy"/>).
        /// </param>
        int[] Merge(GaussianComponent[] components, double[,]? responsibilities = null);
    }
}
