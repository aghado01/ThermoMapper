using System;
using System.Collections.Generic;

namespace Clustering.Evaluation.External;

/// <summary>
/// Internal helper that builds a contingency table from two parallel
/// label arrays. The shape every entropy-based external index (NMI,
/// ARI, Homogeneity, Completeness, V-measure) starts from — densify
/// both label sets once and tally pairwise co-occurrences.
/// </summary>
/// <remarks>
/// <para><b>Output.</b> A dense <c>int[rows, cols]</c> matrix where
/// <c>rows</c> is the predicted-class count and <c>cols</c> is the
/// reference-class count. <c>counts[r, c]</c> is the number of
/// observations whose predicted label maps to index <c>r</c> and whose
/// reference label maps to index <c>c</c>. The dense densified spaces
/// are also returned via out parameters so callers can use them for
/// marginal calculations without re-densifying.</para>
/// </remarks>
internal static class ContingencyTable
{
    /// <summary>
    /// Build the densified contingency table for <paramref name="predictedLabels"/>
    /// and <paramref name="referenceLabels"/>. Both arrays must have the
    /// same length; labels need not be dense — they are mapped to
    /// <c>[0, K)</c> internally.
    /// </summary>
    /// <param name="predictedLabels">Predicted cluster labels.</param>
    /// <param name="referenceLabels">Reference / ground-truth labels.</param>
    /// <param name="rowCount">Number of distinct predicted labels (table rows).</param>
    /// <param name="colCount">Number of distinct reference labels (table cols).</param>
    /// <param name="rowSums">Per-row sums (marginal counts of predicted clusters).</param>
    /// <param name="colSums">Per-col sums (marginal counts of reference classes).</param>
    /// <returns>Dense <c>int[rowCount, colCount]</c> contingency matrix.</returns>
    public static int[,] Build(
        int[] predictedLabels,
        int[] referenceLabels,
        out int rowCount,
        out int colCount,
        out int[] rowSums,
        out int[] colSums)
    {
        ArgumentNullException.ThrowIfNull(predictedLabels);
        ArgumentNullException.ThrowIfNull(referenceLabels);
        if (predictedLabels.Length != referenceLabels.Length)
            throw new ArgumentException(
                $"predictedLabels length ({predictedLabels.Length}) does not match " +
                $"referenceLabels length ({referenceLabels.Length}).");

        int n = predictedLabels.Length;
        if (n == 0)
        {
            rowCount = 0;
            colCount = 0;
            rowSums = Array.Empty<int>();
            colSums = Array.Empty<int>();
            return new int[0, 0];
        }

        // Densify both label spaces in one pass each.
        var rowMap = new Dictionary<int, int>();
        var colMap = new Dictionary<int, int>();
        var predDense = new int[n];
        var refDense = new int[n];
        int nextRow = 0;
        int nextCol = 0;

        for (int i = 0; i < n; i++)
        {
            int pred = predictedLabels[i];
            int refLabel = referenceLabels[i];

            if (!rowMap.TryGetValue(pred, out int r))
            {
                r = nextRow++;
                rowMap[pred] = r;
            }
            predDense[i] = r;

            if (!colMap.TryGetValue(refLabel, out int c))
            {
                c = nextCol++;
                colMap[refLabel] = c;
            }
            refDense[i] = c;
        }

        rowCount = nextRow;
        colCount = nextCol;

        var counts = new int[rowCount, colCount];
        var rs = new int[rowCount];
        var cs = new int[colCount];
        for (int i = 0; i < n; i++)
        {
            int r = predDense[i];
            int c = refDense[i];
            counts[r, c]++;
            rs[r]++;
            cs[c]++;
        }

        rowSums = rs;
        colSums = cs;
        return counts;
    }
}
