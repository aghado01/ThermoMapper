using System;

namespace Clustering.Evaluation.External;

/// <summary>
/// Normalized Mutual Information (NMI), arithmetic-mean normalization:
/// <c>NMI = 2·I(U;V) / (H(U) + H(V))</c> where <c>I</c> is mutual
/// information and <c>H</c> is Shannon entropy of the marginal cluster
/// distributions. Range <c>[0, 1]</c>; 1.0 indicates perfect agreement
/// (up to relabeling).
/// </summary>
/// <remarks>
/// <para><b>Normalization choice.</b> Several variants exist (arithmetic
/// mean, geometric mean <c>sqrt(H(U)·H(V))</c>, maximum). The
/// arithmetic mean form is what scikit-learn returns as
/// <c>normalized_mutual_info_score</c> and is the most commonly cited
/// "NMI" in clustering literature; this evaluator returns that form.</para>
///
/// <para><b>Symmetry.</b> Symmetric in the two arguments — swapping
/// predicted and reference labels gives the same value.</para>
///
/// <para><b>Edge cases.</b> Returns 1.0 when both partitions are
/// trivial (single cluster — both entropies zero, undefined ratio
/// resolved to perfect-agreement by convention). Returns 0.0 for
/// empty input.</para>
/// </remarks>
public sealed class NormalizedMutualInformation : IExternalClusterEvaluator
{
    public string Name => "NMI";

    public double Evaluate(int[] predictedLabels, int[] referenceLabels)
    {
        // Score over the assigned subset only — points the clustering declined
        // to place (predicted == Unassigned) are dropped, not densified into a
        // spurious cluster. The denominator n becomes the assigned count.
        (predictedLabels, referenceLabels) =
            EvaluationHelpers.AssignedByPredicted(predictedLabels, referenceLabels);

        int[,] counts = ContingencyTable.Build(
            predictedLabels, referenceLabels,
            out int rows, out int cols,
            out int[] rowSums, out int[] colSums);

        int n = predictedLabels.Length;
        if (n == 0) return 0.0;

        // Mutual information: Σ_{r,c} p(r,c) · log( p(r,c) / (p(r)·p(c)) )
        // computed in nats then ratio-normalized (units cancel).
        double mi = 0.0;
        double invN = 1.0 / n;
        for (int r = 0; r < rows; r++)
        {
            int rs = rowSums[r];
            if (rs == 0) continue;
            for (int c = 0; c < cols; c++)
            {
                int nrc = counts[r, c];
                if (nrc == 0) continue;
                int cs = colSums[c];
                // p(r,c) = nrc / n; p(r) = rs/n; p(c) = cs/n.
                // log( (nrc/n) / ((rs/n)·(cs/n)) ) = log( n·nrc / (rs·cs) )
                mi += (nrc * invN) * Math.Log((double)n * nrc / ((double)rs * cs));
            }
        }

        double hPred = MarginalEntropy(rowSums, n);
        double hRef  = MarginalEntropy(colSums, n);
        double denom = hPred + hRef;

        if (denom <= 0.0)
        {
            // Both partitions are trivial (single cluster). Convention:
            // they "agree" trivially, return 1.
            return 1.0;
        }

        return Math.Clamp(2.0 * mi / denom, 0.0, 1.0);
    }

    private static double MarginalEntropy(int[] sums, int n)
    {
        if (n == 0) return 0.0;
        double invN = 1.0 / n;
        double h = 0.0;
        for (int i = 0; i < sums.Length; i++)
        {
            int s = sums[i];
            if (s == 0) continue;
            double p = s * invN;
            h -= p * Math.Log(p);
        }
        return h;
    }
}
