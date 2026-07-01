using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;

namespace Clustering.Graphical.HdbScan;

/// <summary>
/// Result of one end-to-end <see cref="HdbscanSession"/> run. Composes the raw
/// <see cref="HdbscanResult"/> (kept as the single source of truth for labels /
/// probabilities / dendrogram) with the evaluator-score dictionary, mirroring
/// <c>SpcSessionResult</c>. Per-cluster stats, CSV, and JSON are presentation
/// and stay in the caller.
/// </summary>
public sealed record HdbscanSessionResult(
    HdbscanResult                       Result,
    IReadOnlyDictionary<string, double> EvaluatorScores)
{
    public int[]      Labels                  => Result.Labels;
    public double[]   MembershipProbabilities => Result.MembershipProbabilities;
    public int        ClusterCount            => Result.ClusterCount;
    public Dendrogram Dendrogram              => Result.Dendrogram;

    /// <summary>Count of points labelled noise (-1).</summary>
    public int NoiseCount
    {
        get
        {
            int c = 0;
            int[] labels = Result.Labels;
            for (int i = 0; i < labels.Length; i++)
                if (labels[i] < 0) c++;
            return c;
        }
    }

    /// <summary>Convenience accessor for the <c>"Purity"</c> entry in
    /// <see cref="EvaluatorScores"/>; <see langword="null"/> when no Purity
    /// evaluator was supplied. Mirrors <c>SpcSessionResult.Purity</c>.</summary>
    public double? Purity =>
        EvaluatorScores.TryGetValue("Purity", out double v) ? v : null;
}

/// <summary>
/// Thin backend orchestration facade for HDBSCAN — the CLI / Mapper / smoke all
/// call into this instead of hand-rolling the same flatten → core-distance
/// clamp → struct-metric dispatch → evaluator pass. Deliberately does <b>not</b>
/// route through <c>GraphCompiler</c>: exact HDBSCAN needs the implicit dense
/// mutual-reachability MST built inside <c>Prim.ComputeMutualReachabilityMst</c>,
/// not a sparse kNN <c>CsrGraph</c>, so there is no graph object to share with
/// SPC here. (The CSR-sparse / approximate overload that consumes a
/// <c>DistanceProjection</c> graph is a deliberate follow-on.)
/// </summary>
public static class HdbscanSession
{
    /// <param name="features">Row-per-point feature matrix; N ≥ 2, uniform
    /// dimensionality.</param>
    /// <param name="settings">Run configuration; defaults to a fresh
    /// <see cref="HdbscanSettings"/> when null.</param>
    /// <param name="externalEvaluators">Optional external cluster evaluators.
    /// Run only when <paramref name="referenceLabels"/> are also provided, so a
    /// caller can pass the standard set unconditionally.</param>
    /// <param name="referenceLabels">Optional ground-truth labels for evaluator
    /// scoring. HDBSCAN noise (-1) is remapped to its own sentinel cluster before
    /// scoring (external evaluators have no noise semantic).</param>
    public static HdbscanSessionResult Run(
        double[][]                               features,
        HdbscanSettings?                         settings           = null,
        IEnumerable<IExternalClusterEvaluator>?  externalEvaluators = null,
        int[]?                                   referenceLabels    = null)
    {
        if (features is null) throw new ArgumentNullException(nameof(features));
        if (features.Length < 2)
            throw new ArgumentException("HDBSCAN requires at least 2 points.", nameof(features));

        var cfg = settings ?? new HdbscanSettings();

        int n   = features.Length;
        int dim = features[0].Length;

        // Floor at 2 (runner invariant), cap at n-1 so core-distance can still
        // find a kth neighbour. Was duplicated across command + mapper.
        int effMinPts = Math.Max(2, Math.Min(cfg.MinPts, n - 1));

        // Sparse-MST kNN graph degree: default max(minPts, 10), floored at minPts,
        // capped at n-1. Ignored on the dense path.
        int effGraphK = Math.Min(n - 1, Math.Max(effMinPts, cfg.GraphNeighbors ?? Math.Max(effMinPts, 10)));

        // Flatten to row-major for the runner. Was duplicated 3×.
        var flat = new double[n * dim];
        for (int i = 0; i < n; i++)
            Array.Copy(features[i], 0, flat, i * dim, dim);

        var runner = new HdbscanRunner(n);
        HdbscanResult result = HdbscanMetricDispatch.Run(
            runner, cfg.Metric,
            flat, dim, effMinPts, cfg.MinClusterSize, cfg.AllowSingleCluster,
            cfg.ClusterSelectionMethod, cfg.ClusterSelectionEpsilon, cfg.MstAlgorithm, effGraphK);

        var scores = new Dictionary<string, double>(StringComparer.Ordinal);
        if (externalEvaluators is not null && referenceLabels is { Length: > 0 })
        {
            int[] labelsForEval = MapNoiseToOwnCluster(result.Labels, result.ClusterCount);
            foreach (var ev in externalEvaluators)
                scores[ev.Name] = ev.Evaluate(labelsForEval, referenceLabels);
        }

        return new HdbscanSessionResult(result, scores);
    }

    /// <summary>
    /// External evaluators (Purity, NMI, ARI, ...) consume dense integer labels
    /// with no special-case noise semantic, so map -1 → <paramref name="clusterCount"/>
    /// (a fresh "noise" cluster) and score that against the reference distribution
    /// rather than silently treating -1 as a regular cluster. HDBSCAN-specific
    /// evaluation convention every caller should apply identically — hence it
    /// lives here, not in the command IO.
    /// </summary>
    internal static int[] MapNoiseToOwnCluster(int[] labels, int clusterCount)
    {
        var mapped = new int[labels.Length];
        for (int i = 0; i < labels.Length; i++)
            mapped[i] = labels[i] < 0 ? clusterCount : labels[i];
        return mapped;
    }
}
