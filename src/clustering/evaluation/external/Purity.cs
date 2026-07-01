using System.Collections.Generic;

namespace Clustering.Evaluation.External;

/// <summary>
/// Purity: <c>(1/N) · Σ_p max_r |p ∩ r|</c> — the fraction of
/// observations whose predicted cluster is dominated by their reference
/// class. Range <c>[0, 1]</c>; 1.0 indicates perfect agreement.
/// </summary>
/// <remarks>
/// <para><b>Reading.</b> Purity rewards predicted clusters that contain
/// a single reference class but does <i>not</i> penalize
/// over-fragmentation — the trivial partition "every observation in its
/// own cluster" achieves purity 1.0 by construction. Read alongside
/// the predicted cluster count, or paired with a complementary index
/// like NMI/ARI that penalizes fragmentation, to detect this failure
/// mode.</para>
///
/// <para><b>Asymmetry.</b> Purity is not symmetric in its arguments:
/// <c>Purity(pred, truth) ≠ Purity(truth, pred)</c> in general. The
/// canonical convention — and the one used here — is to take the max
/// over <i>reference</i> classes within each <i>predicted</i> cluster.</para>
/// </remarks>
public sealed class Purity : IExternalClusterEvaluator
{
    public string Name => "Purity";

    public double Evaluate(int[] predictedLabels, int[] referenceLabels)
    {
        // Score over the assigned subset only — unassigned predictions are
        // dropped (the denominator n becomes the assigned count). The helper
        // also performs the null / equal-length validation.
        (predictedLabels, referenceLabels) =
            EvaluationHelpers.AssignedByPredicted(predictedLabels, referenceLabels);

        int n = predictedLabels.Length;
        if (n == 0) return 0.0;

        // Contingency: (predictedCluster, referenceClass) → count.
        var counts = new Dictionary<(int Pred, int Ref), int>();
        for (int i = 0; i < n; i++)
        {
            var key = (predictedLabels[i], referenceLabels[i]);
            counts[key] = counts.TryGetValue(key, out int c) ? c + 1 : 1;
        }

        // For each predicted cluster, find the max count over reference classes.
        var maxPerCluster = new Dictionary<int, int>();
        foreach (var (key, count) in counts)
        {
            if (!maxPerCluster.TryGetValue(key.Pred, out int current) || count > current)
                maxPerCluster[key.Pred] = count;
        }

        long correctlyAssigned = 0;
        foreach (var v in maxPerCluster.Values) correctlyAssigned += v;
        return (double)correctlyAssigned / n;
    }

    /// <summary>
    /// Convenience helper for one-off purity computation.
    /// </summary>
    public static double Compute(int[] predictedLabels, int[] referenceLabels)
        => new Purity().Evaluate(predictedLabels, referenceLabels);
}
