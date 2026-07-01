// ============================================================================
// TDA.Mapper — MapperGMM.cs
// ============================================================================
// GMM-backed IClusterer adapters for Mapper preimage clustering.
//
// Not every GMM paradigm in Clustering.GMM is appropriate here. Mapper calls
// IClusterer once per preimage patch (typically 10–200 points) with no external
// labels and expects a hard partition. Adapters in this file must be:
//   - Self-contained per call (no cross-patch state)
//   - Able to determine cluster count from the patch itself
//   - Meaningful at small N
//
// Current adapters:
//   ModalGmmClusterer  — overspecified EM + mode-topology collapse via ModeAscent
//                        + ModalMergeStrategy. Cluster count emerges from the
//                        density landscape, not from K directly.
//
// Paradigms intentionally absent:
//   Direct EM (single fit, heuristic K) — cluster count is arbitrary; no
//     advantage over KMeansPlusPlusClusterer for Mapper preimages.
//   Semi-supervised — requires external confidence labels; Mapper has none.
//   Entropy merge alone — cut point is still a K-selection problem.
// ============================================================================

#nullable enable
using System;
using System.Linq;
using Clustering.Statistical.GMM;
using TDA.Mapper;

namespace TDA.Mapper.Clusterers;

/// <summary>
/// IClusterer adapter using the GMM modal-topology paradigm.
/// <para>
/// Fits a GMM with an overspecified K, then runs gradient ascent from each
/// component mean to the nearest local mode of the mixture log-density
/// (<see cref="ModeAscent"/>). Components that converge to the same mode are
/// merged into one cluster (<see cref="ModalMergeStrategy"/>). The returned
/// cluster count reflects the number of distinct density modes, not the EM K.
/// </para>
/// </summary>
public sealed class ModalGmmClusterer : IClusterer
{
    public int    KMax          { get; }
    public double ModeTolerance { get; }
    public int    MaxIterations { get; }
    public int    Seed          { get; }
    public double Regularization { get; }

    public string Name => $"GMM-Modal (kMax={KMax})";

    public ModalGmmClusterer(
        int    kMax           = 12,
        double modeTolerance  = double.NaN,
        int    maxIterations  = 100,
        int    seed           = 42,
        double regularization = 1e-6)
    {
        if (kMax < 2)          throw new ArgumentOutOfRangeException(nameof(kMax), "kMax must be ≥ 2.");
        if (maxIterations < 1) throw new ArgumentOutOfRangeException(nameof(maxIterations));
        if (regularization < 0) throw new ArgumentOutOfRangeException(nameof(regularization));

        KMax           = kMax;
        ModeTolerance  = modeTolerance;
        MaxIterations  = maxIterations;
        Seed           = seed;
        Regularization = regularization;
    }

    public ClusterResult Cluster(double[][] subset)
    {
        if (subset is null || subset.Length == 0)
            return new ClusterResult(Array.Empty<int>(), 0);

        if (subset.Length == 1)
            return new ClusterResult(new[] { 0 }, 1);

        int n   = subset.Length;
        int dim = subset[0].Length;

        // Overspecify K: use kMax but ensure at least 2 points per component
        // so EM remains numerically stable. The modal merge will collapse the excess.
        int k = Math.Min(KMax, Math.Max(2, n / 2));

        var rng = new Random(Seed);
        var gmm = new GaussianMixtureModel(k, dim, Regularization);

        if (n >= k * 2)
            gmm.RobustInitialize(subset, rng: rng);
        else
            gmm.RandomInitialize(subset, rng);

        gmm.Fit(subset, MaxIterations);

        // Component → cluster via mode topology.
        // ModalMergeStrategy does not use the responsibility matrix.
        var mergeStrategy = new ModalMergeStrategy(ModeTolerance);
        int[] mergeMap    = mergeStrategy.Merge(gmm.Components);

        // Compose: point → component → cluster
        int[] componentLabels = gmm.Predict(subset);
        int[] clusterLabels   = new int[n];
        for (int i = 0; i < n; i++)
            clusterLabels[i] = mergeMap[componentLabels[i]];

        int nClusters = mergeMap.Max() + 1;
        return new ClusterResult(clusterLabels, nClusters);
    }
}
