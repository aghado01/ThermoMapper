using System;
using Graphs.Distance;

namespace Clustering.Graphical.HdbScan;

/// <summary>
/// Bridges a metric spec string to the struct-generic <see cref="HdbscanRunner.Run{TMetric}"/>
/// via the canonical <see cref="MetricRegistry"/> — no metric-family switch of its
/// own. The registry owns the spec vocabulary + the kind→struct dispatch (so this
/// and <c>DistanceMetricFactory</c> can't drift); HDBSCAN supplies only a visitor
/// that runs the algorithm with the chosen struct, preserving JIT inlining in
/// Prim's loop.
/// </summary>
internal static class HdbscanMetricDispatch
{
    public static HdbscanResult Run(
        HdbscanRunner runner, string metricSpec,
        double[] data, int dim, int minPts,
        int? minClusterSize, bool allowSingleCluster,
        ClusterSelectionMethod selectionMethod, double clusterSelectionEpsilon,
        MstAlgorithm mstAlgorithm, int graphNeighbors)
    {
        MetricSpec spec = MetricRegistry.Parse(metricSpec);
        return MetricRegistry.Invoke(
            spec,
            new RunVisitor(runner, data, dim, minPts, minClusterSize, allowSingleCluster,
                selectionMethod, clusterSelectionEpsilon, mstAlgorithm, graphNeighbors));
    }

    /// <summary>Carries the run parameters so <see cref="Visit{TMetric}"/> can call
    /// the struct-generic runner with the concrete metric the registry selected.
    /// Holds <c>double[]</c> (not a <c>ReadOnlySpan</c>) so it can live as a field;
    /// the array converts to a span at the call.</summary>
    private readonly struct RunVisitor : IMetricVisitor<HdbscanResult>
    {
        private readonly HdbscanRunner          _runner;
        private readonly double[]               _data;
        private readonly int                    _dim;
        private readonly int                    _minPts;
        private readonly int?                   _minClusterSize;
        private readonly bool                   _allowSingleCluster;
        private readonly ClusterSelectionMethod _selectionMethod;
        private readonly double                 _clusterSelectionEpsilon;
        private readonly MstAlgorithm           _mstAlgorithm;
        private readonly int                    _graphNeighbors;

        public RunVisitor(
            HdbscanRunner runner, double[] data, int dim, int minPts,
            int? minClusterSize, bool allowSingleCluster,
            ClusterSelectionMethod selectionMethod, double clusterSelectionEpsilon,
            MstAlgorithm mstAlgorithm, int graphNeighbors)
        {
            _runner                  = runner;
            _data                    = data;
            _dim                     = dim;
            _minPts                  = minPts;
            _minClusterSize          = minClusterSize;
            _allowSingleCluster      = allowSingleCluster;
            _selectionMethod         = selectionMethod;
            _clusterSelectionEpsilon = clusterSelectionEpsilon;
            _mstAlgorithm            = mstAlgorithm;
            _graphNeighbors          = graphNeighbors;
        }

        public HdbscanResult Visit<TMetric>(TMetric metric) where TMetric : struct, IDistanceMetric
            => _mstAlgorithm == MstAlgorithm.SparseKnn
                ? _runner.RunSparse<TMetric>(_data, _dim, _minPts, metric, _graphNeighbors, _minClusterSize, _allowSingleCluster, _selectionMethod, _clusterSelectionEpsilon)
                : _runner.Run<TMetric>(_data, _dim, _minPts, metric, _minClusterSize, _allowSingleCluster, _selectionMethod, _clusterSelectionEpsilon);
    }
}
