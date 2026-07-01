// ============================================================================
// TDA.Mapper — MapperHdbScan.cs
// ============================================================================
// IClusterer adapter for HDBSCAN (Clustering.Graphical.HdbScan.HdbscanRunner).
//
// Defaults are tuned for Mapper preimages (typically 10–200 points, dense):
//   minPts=5, minClusterSize=null (= minPts), allowSingleCluster=true.
// Other callers should have their calling defaults tuned to their needs.
//
// allowSingleCluster=true means a bin that is one dense blob with a few
// outliers returns one cluster (with low-probability outliers) instead of
// all-noise — usually what Mapper wants. Set to false to recover sklearn-
// default semantics.
//
// HDBSCAN labels noise points as -1. Those flow through Mapper as an extra
// LocalClusterId per bin; downstream consumers that want to exclude noise
// from the nerve should filter -1 nodes post-hoc.
// ============================================================================

#nullable enable
using System;
using Clustering.Graphical.HdbScan;
using TDA.Mapper;

namespace TDA.Mapper.Clusterers;

public sealed class HdbscanClusterer : IClusterer
{
    public int  MinPts             { get; }
    public int? MinClusterSize     { get; }
    public bool AllowSingleCluster { get; }

    public string Name => MinClusterSize is int mcs
        ? $"HDBSCAN (minPts={MinPts}, minSz={mcs})"
        : $"HDBSCAN (minPts={MinPts})";

    public HdbscanClusterer(
        int  minPts             = 5,
        int? minClusterSize     = null,
        bool allowSingleCluster = true)
    {
        if (minPts < 2)
            throw new ArgumentOutOfRangeException(nameof(minPts), "minPts must be >= 2.");
        if (minClusterSize is int mcs && mcs < 2)
            throw new ArgumentOutOfRangeException(nameof(minClusterSize), "minClusterSize must be >= 2.");

        MinPts             = minPts;
        MinClusterSize     = minClusterSize;
        AllowSingleCluster = allowSingleCluster;
    }

    public ClusterResult Cluster(double[][] subset)
    {
        if (subset is null || subset.Length == 0)
            return new ClusterResult(Array.Empty<int>(), 0);
        if (subset.Length == 1)
            return new ClusterResult(new[] { 0 }, 1);

        // Flatten, the core-distance clamp, and runner construction all live in
        // HdbscanSession now. Mapper preimages are Euclidean, so the default
        // Metric ("euclidean") applies; no evaluators on this path.
        var result = HdbscanSession.Run(subset, new HdbscanSettings
        {
            MinPts             = MinPts,
            MinClusterSize     = MinClusterSize,
            AllowSingleCluster = AllowSingleCluster,
        });

        return new ClusterResult(result.Labels, result.ClusterCount);
    }
}
