using System;

namespace Clustering.Evaluation.External;

/// <summary>
/// Adjusted Rand Index (ARI): the Rand Index corrected for chance
/// agreement.
/// <code>ARI = (RI − E[RI]) / (max(RI) − E[RI])</code>
/// computed via the contingency-table identity
/// <code>ARI = (sum_index − expected) / (max_index − expected)</code>
/// where <c>sum_index = Σ C(n_ij, 2)</c>,
/// <c>expected = sum_a · sum_b / C(n, 2)</c>,
/// <c>max_index = (sum_a + sum_b) / 2</c>.
/// </summary>
/// <remarks>
/// <para><b>Range.</b> Up to <c>1.0</c> (perfect agreement);
/// <c>0.0</c> = chance-level agreement; can dip slightly negative for
/// worse-than-random agreement. <b>Higher is better.</b></para>
///
/// <para><b>Why ARI over raw Rand.</b> The raw Rand Index is biased
/// upward — random partitions of similar size achieve high scores by
/// accident. ARI subtracts the expected value under a permutation
/// null model so a score of 0 genuinely means "no better than chance,"
/// which makes scores comparable across datasets of different sizes
/// and class counts.</para>
///
/// <para><b>Symmetry.</b> Symmetric in its two arguments.</para>
///
/// <para><b>Edge cases.</b> Returns 1.0 when both inputs are
/// identical-partition (or trivially when both have a single cluster
/// and N ≤ 1). Returns 0.0 when both are empty.</para>
/// </remarks>
public sealed class AdjustedRandIndex : IExternalClusterEvaluator
{
    public string Name => "ARI";

    public double Evaluate(int[] predictedLabels, int[] referenceLabels)
    {
        // Score over the assigned subset only — unassigned predictions are
        // dropped before the contingency table is built (the denominator n
        // becomes the assigned count).
        (predictedLabels, referenceLabels) =
            EvaluationHelpers.AssignedByPredicted(predictedLabels, referenceLabels);

        int[,] counts = ContingencyTable.Build(
            predictedLabels, referenceLabels,
            out int rows, out int cols,
            out int[] rowSums, out int[] colSums);

        int n = predictedLabels.Length;
        if (n <= 1) return n == 0 ? 0.0 : 1.0;

        // sum_index = Σ C(n_ij, 2)
        // sum_a     = Σ C(rowSums[r], 2)
        // sum_b     = Σ C(colSums[c], 2)
        double sumIndex = 0.0;
        for (int r = 0; r < rows; r++)
        {
            for (int c = 0; c < cols; c++)
            {
                int nrc = counts[r, c];
                if (nrc >= 2) sumIndex += nrc * (nrc - 1) * 0.5;
            }
        }

        double sumA = 0.0;
        for (int r = 0; r < rows; r++)
        {
            int rs = rowSums[r];
            if (rs >= 2) sumA += rs * (rs - 1) * 0.5;
        }

        double sumB = 0.0;
        for (int c = 0; c < cols; c++)
        {
            int cs = colSums[c];
            if (cs >= 2) sumB += cs * (cs - 1) * 0.5;
        }

        double totalPairs = n * (n - 1) * 0.5;          // C(n, 2)
        double expected   = sumA * sumB / totalPairs;
        double maxIndex   = (sumA + sumB) * 0.5;

        double denom = maxIndex - expected;
        if (Math.Abs(denom) < double.Epsilon)
        {
            // Both partitions are effectively trivial (all in one cluster
            // either way) — perfect agreement by convention.
            return 1.0;
        }

        return (sumIndex - expected) / denom;
    }
}
