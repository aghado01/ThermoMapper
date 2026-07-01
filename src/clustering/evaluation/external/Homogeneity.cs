using System;

namespace Clustering.Evaluation.External;

/// <summary>
/// Homogeneity: a partition is homogeneous if each predicted cluster
/// contains observations from only a single reference class.
/// <code>h = 1 − H(reference | predicted) / H(reference)</code>
/// </summary>
/// <remarks>
/// <para><b>Range.</b> <c>[0, 1]</c>; <c>1.0</c> means every predicted
/// cluster is purely from one reference class. <b>Higher is better.</b></para>
///
/// <para><b>Asymmetric.</b> Homogeneity is <i>not</i> symmetric:
/// <c>h(pred, truth) ≠ h(truth, pred)</c>. It measures one direction
/// — "do my clusters preserve class identity?" Pair with
/// <see cref="Completeness"/> ("are class members kept together?") for
/// the symmetric pair, and <see cref="VMeasure"/> for the harmonic
/// mean.</para>
///
/// <para><b>Edge cases.</b> Returns 1.0 when the reference partition
/// is trivial (single class — perfectly homogeneous by definition).
/// Returns 0.0 for empty input.</para>
/// </remarks>
public sealed class Homogeneity : IExternalClusterEvaluator
{
    public string Name => "Homogeneity";

    public double Evaluate(int[] predictedLabels, int[] referenceLabels)
    {
        // Score over the assigned subset only — unassigned predictions are
        // dropped before densification (the denominator n becomes the
        // assigned count).
        (predictedLabels, referenceLabels) =
            EvaluationHelpers.AssignedByPredicted(predictedLabels, referenceLabels);

        int[,] counts = ContingencyTable.Build(
            predictedLabels, referenceLabels,
            out int rows, out int cols,
            out int[] rowSums, out int[] colSums);

        int n = predictedLabels.Length;
        if (n == 0) return 0.0;

        double hRef = MarginalEntropy(colSums, n);
        if (hRef <= 0.0) return 1.0;   // trivial reference: perfect by convention

        // H(reference | predicted) = -Σ_{r,c} (n_rc / n) · log(n_rc / rowSum_r)
        double hRefGivenPred = ConditionalEntropy(counts, rows, cols, rowSums, n);

        return Math.Clamp(1.0 - hRefGivenPred / hRef, 0.0, 1.0);
    }

    internal static double MarginalEntropy(int[] sums, int n)
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

    /// <summary>
    /// H(reference | predicted): conditional entropy of the reference
    /// labels given the predicted clusters. Iterates row-by-row over
    /// the contingency table.
    /// </summary>
    internal static double ConditionalEntropy(
        int[,] counts, int rows, int cols, int[] givenSums, int n)
    {
        if (n == 0) return 0.0;
        double invN = 1.0 / n;
        double h = 0.0;
        for (int r = 0; r < rows; r++)
        {
            int rs = givenSums[r];
            if (rs == 0) continue;
            double invRs = 1.0 / rs;
            for (int c = 0; c < cols; c++)
            {
                int nrc = counts[r, c];
                if (nrc == 0) continue;
                // (n_rc / n) · log(rs / n_rc)
                h -= (nrc * invN) * Math.Log(nrc * invRs);
            }
        }
        return h;
    }
}
