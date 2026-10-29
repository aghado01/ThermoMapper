namespace Clustering.Evaluation.External;

/// <summary>
/// Rand Index (RI; Rand 1971): the fraction of point pairs on which two
/// partitions agree — together in both, or apart in both.
/// <code>RI = (C(n,2) + 2·Σ C(n_ij,2) − Σ C(a_i,2) − Σ C(b_j,2)) / C(n,2)</code>
/// computed from the contingency table. It is the Potts link overlap of the two
/// partitions taken over the complete graph; <c>Graphs.Observables.LinkOverlap</c> is the
/// same count restricted to a graph's edges.
/// </summary>
/// <remarks>
/// <para><b>Range.</b> <c>[0, 1]</c>; <c>1.0</c> = identical partitions.
/// <b>Higher is better.</b> Not chance-corrected: random partitions of similar
/// size score high, so compare across datasets with <see cref="AdjustedRandIndex"/>.</para>
///
/// <para><b>Symmetry.</b> Symmetric in its two arguments.</para>
///
/// <para><b>Edge cases.</b> Scores the points whose predicted label is assigned.
/// Returns 0.0 when none are, and 1.0 for a single point (no pairs).</para>
/// </remarks>
public sealed class RandIndex : IExternalClusterEvaluator
{
    public string Name => "RI";

    public double Evaluate(int[] predictedLabels, int[] referenceLabels)
    {
        // Score over the assigned subset only, as every external index does.
        (predictedLabels, referenceLabels) =
            EvaluationHelpers.AssignedByPredicted(predictedLabels, referenceLabels);

        int n = predictedLabels.Length;
        if (n <= 1) return n == 0 ? 0.0 : 1.0;

        int[,] counts = ContingencyTable.Build(
            predictedLabels, referenceLabels,
            out int rows, out int cols,
            out int[] rowSums, out int[] colSums);

        double together = 0.0;                 // Σ C(n_ij, 2): pairs together in both
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
                together += Pairs(counts[r, c]);

        double togetherA = 0.0;                // Σ C(a_i, 2): pairs together in the prediction
        for (int r = 0; r < rows; r++) togetherA += Pairs(rowSums[r]);

        double togetherB = 0.0;                // Σ C(b_j, 2): pairs together in the reference
        for (int c = 0; c < cols; c++) togetherB += Pairs(colSums[c]);

        double totalPairs = Pairs(n);
        // Agreements = together in both + apart in both
        //            = together + (total − togetherA − togetherB + together).
        return (totalPairs + 2.0 * together - togetherA - togetherB) / totalPairs;
    }

    // In double: an int product m·(m − 1) overflows past m ≈ 46,341.
    private static double Pairs(int m) => m * (m - 1.0) * 0.5;
}
