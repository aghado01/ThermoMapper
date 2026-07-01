using System;

namespace Clustering.Evaluation.External;

/// <summary>
/// Completeness: a partition is complete if all observations from a
/// reference class are grouped into a single predicted cluster.
/// <code>c = 1 − H(predicted | reference) / H(predicted)</code>
/// </summary>
/// <remarks>
/// <para><b>Range.</b> <c>[0, 1]</c>; <c>1.0</c> means every reference
/// class is fully contained in a single predicted cluster.
/// <b>Higher is better.</b></para>
///
/// <para><b>Asymmetric.</b> Like <see cref="Homogeneity"/>, this is
/// directional. Completeness asks "do my clusters preserve class
/// togetherness?" — the dual question to homogeneity's "do my clusters
/// preserve class identity?". <see cref="VMeasure"/> combines the two
/// as a harmonic mean.</para>
///
/// <para><b>Symmetry with Homogeneity.</b>
/// <c>Completeness(pred, truth) == Homogeneity(truth, pred)</c>.</para>
///
/// <para><b>Edge cases.</b> Returns 1.0 when the predicted partition
/// is trivial (single cluster — every reference class is trivially
/// fully contained). Returns 0.0 for empty input.</para>
/// </remarks>
public sealed class Completeness : IExternalClusterEvaluator
{
    public string Name => "Completeness";

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

        double hPred = Homogeneity.MarginalEntropy(rowSums, n);
        if (hPred <= 0.0) return 1.0;   // trivial predicted: perfect by convention

        // H(predicted | reference): walk the contingency table
        // column-by-column. Same shape as H(ref|pred) with axes swapped,
        // so we transpose-iterate rather than rebuild the table.
        double hPredGivenRef = ConditionalEntropyColumnMajor(counts, rows, cols, colSums, n);

        return Math.Clamp(1.0 - hPredGivenRef / hPred, 0.0, 1.0);
    }

    /// <summary>
    /// H(predicted | reference): the column-major counterpart of
    /// <see cref="Homogeneity.ConditionalEntropy"/>.
    /// </summary>
    private static double ConditionalEntropyColumnMajor(
        int[,] counts, int rows, int cols, int[] colSums, int n)
    {
        if (n == 0) return 0.0;
        double invN = 1.0 / n;
        double h = 0.0;
        for (int c = 0; c < cols; c++)
        {
            int cs = colSums[c];
            if (cs == 0) continue;
            double invCs = 1.0 / cs;
            for (int r = 0; r < rows; r++)
            {
                int nrc = counts[r, c];
                if (nrc == 0) continue;
                h -= (nrc * invN) * Math.Log(nrc * invCs);
            }
        }
        return h;
    }
}
