using System.Collections.Generic;

namespace Clustering.Geometric.EMeans;

/// <summary>One epoch of a run.</summary>
/// <param name="Epoch">Zero-based epoch index.</param>
/// <param name="OnlineDistortion">Mean squared distance from each point to the centroid it was assigned to, measured
/// against the centroids as they stood when its batch was assigned.</param>
/// <param name="Revivals">Dead clusters revived during the epoch.</param>
public sealed record EMeansEpochTrace(int Epoch, double OnlineDistortion, int Revivals);

/// <summary>
/// What a run resolved and did: the config as given, the values its auto-picks settled on, the derived seed streams,
/// and the per-epoch trace.
/// </summary>
/// <param name="Config">The config the run was given.</param>
/// <param name="N">Number of points.</param>
/// <param name="Dimension">Coordinates per point.</param>
/// <param name="BatchSize">Resolved batch size <c>B</c>.</param>
/// <param name="PeakStep">Resolved peak step <c>α₀</c>.</param>
/// <param name="SlowStep">The slow EMA's constant step.</param>
/// <param name="Updates">Total EMA updates <c>M = T·⌈n/B⌉</c>.</param>
/// <param name="StreamSeeds">Seeds of the initialization, shuffle and revival streams, in that order.</param>
/// <param name="AutomationLog">Which auto-picks fired, what they chose, and any fallback taken.</param>
/// <param name="Epochs">Per-epoch trace.</param>
public sealed record EMeansManifest(
    EMeansConfig Config,
    int N,
    int Dimension,
    int BatchSize,
    double PeakStep,
    double SlowStep,
    int Updates,
    int[] StreamSeeds,
    IReadOnlyList<string> AutomationLog,
    IReadOnlyList<EMeansEpochTrace> Epochs);

/// <summary>
/// An e-means run's outcome. The decayed statistics are retained as they stood after the last update; the three
/// centroid sets are views of the run, and <see cref="Centroids"/> is the one the config's finalization selected.
/// </summary>
/// <param name="Labels">Nearest centroid in <see cref="Centroids"/> for every point.</param>
/// <param name="Centroids">The selected centroid set.</param>
/// <param name="K">Number of clusters.</param>
/// <param name="Distortion">Mean squared distance from each point to its centroid under <see cref="Labels"/>.</param>
/// <param name="Counts">The training EMA's counts <c>N</c>.</param>
/// <param name="Sums">The training EMA's sums <c>S</c>.</param>
/// <param name="Scatter">The training EMA's decayed sums of squared distances from each assigned point to its
/// centroid as it stood at assignment; <c>Scatter/Counts</c> is a cluster's mean squared radius.</param>
/// <param name="SlowCounts">The slow EMA's counts.</param>
/// <param name="SlowSums">The slow EMA's sums.</param>
/// <param name="EmaCentroids">The training EMA's centroids <c>S/N</c>.</param>
/// <param name="SlowEmaCentroids">The slow EMA's centroids.</param>
/// <param name="LastEpochMeanCentroids">Each cluster's mean over its last-epoch assignments; a cluster with none
/// carries its EMA centroid.</param>
/// <param name="Manifest">Resolved settings and trace.</param>
public sealed record EMeansResult(
    int[] Labels,
    double[][] Centroids,
    int K,
    double Distortion,
    double[] Counts,
    double[][] Sums,
    double[] Scatter,
    double[] SlowCounts,
    double[][] SlowSums,
    double[][] EmaCentroids,
    double[][] SlowEmaCentroids,
    double[][] LastEpochMeanCentroids,
    EMeansManifest Manifest);
