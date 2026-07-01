#nullable enable
using System;
using System.Collections.Generic;
using System.Numerics.Tensors;
using System.Threading.Tasks;

namespace Hashish;

/// <summary>
/// Standalone IR-style query helpers over a fitted <see cref="TfIdfModel"/>
/// and an external dense row matrix (flat <c>double[N * Dimension]</c>).
/// Scores use cosine similarity, which collapses to a dot product when the
/// model was fit with <see cref="TfIdfOptions.L2Normalize"/> = true (default).
/// </summary>
public static class TfIdfSearch
{
    /// <summary>
    /// Score a free-text query against every row in <paramref name="denseRows"/>
    /// and return the top-K documents by similarity. Query is sparsified internally,
    /// so cost per row is O(|query nnz|) rather than O(Dimension).
    /// </summary>
    public static (int DocId, double Score)[] ScoreQuery(
        TfIdfModel model,
        double[] denseRows,
        string query,
        int topK = 10)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(denseRows);
        ArgumentNullException.ThrowIfNull(query);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(topK);

        int dim = model.Dimension;
        if (denseRows.Length % dim != 0)
            throw new ArgumentException(
                $"denseRows length {denseRows.Length} is not a multiple of model.Dimension {dim}.");
        int n = denseRows.Length / dim;
        if (n == 0) return Array.Empty<(int, double)>();

        (int[] qIdx, double[] qVal) = model.TransformSparse(query);
        if (qIdx.Length == 0) return Array.Empty<(int, double)>();

        var scores = new double[n];
        if (model.Options.Parallel)
        {
            Parallel.For(0, n, d =>
            {
                scores[d] = SparseDenseDot(qIdx, qVal, denseRows, d * dim);
            });
        }
        else
        {
            for (int d = 0; d < n; d++)
                scores[d] = SparseDenseDot(qIdx, qVal, denseRows, d * dim);
        }

        return TopK(scores, topK);
    }

    /// <summary>
    /// Return the top-K most similar documents to <paramref name="sourceDocId"/>.
    /// Self is excluded. Dot-product cost is O(Dimension) per candidate.
    /// </summary>
    public static (int DocId, double Score)[] NearestDocuments(
        TfIdfModel model,
        double[] denseRows,
        int sourceDocId,
        int topK = 10)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(denseRows);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(topK);

        int dim = model.Dimension;
        if (denseRows.Length % dim != 0)
            throw new ArgumentException(
                $"denseRows length {denseRows.Length} is not a multiple of model.Dimension {dim}.");
        int n = denseRows.Length / dim;
        if (n <= 1) return Array.Empty<(int, double)>();

        ArgumentOutOfRangeException.ThrowIfNegative(sourceDocId);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(sourceDocId, n);

        // Snapshot the source row so parallel workers don't slice the shared buffer
        // for the left operand on every dot. Right operand still slices in-loop.
        var sourceRow = new double[dim];
        Array.Copy(denseRows, sourceDocId * dim, sourceRow, 0, dim);

        var scores = new double[n];
        scores[sourceDocId] = double.NegativeInfinity;

        if (model.Options.Parallel)
        {
            Parallel.For(0, n, d =>
            {
                if (d == sourceDocId) return;
                scores[d] = TensorPrimitives.Dot<double>(
                    sourceRow,
                    denseRows.AsSpan(d * dim, dim));
            });
        }
        else
        {
            for (int d = 0; d < n; d++)
            {
                if (d == sourceDocId) continue;
                scores[d] = TensorPrimitives.Dot<double>(
                    sourceRow,
                    denseRows.AsSpan(d * dim, dim));
            }
        }

        return TopK(scores, topK);
    }

    /// <summary>
    /// Sparse-query dot against a dense row at <paramref name="rowOffset"/> of
    /// <paramref name="denseRows"/>. Walks only the query's nonzero positions.
    /// </summary>
    private static double SparseDenseDot(
        int[] queryIndices,
        double[] queryValues,
        double[] denseRows,
        int rowOffset)
    {
        double s = 0.0;
        for (int k = 0; k < queryIndices.Length; k++)
            s += queryValues[k] * denseRows[rowOffset + queryIndices[k]];
        return s;
    }

    /// <summary>
    /// Top-K largest by score, descending. Uses a bounded min-heap so cost is
    /// O(N log K) — meaningful when K is much smaller than N.
    /// </summary>
    private static (int DocId, double Score)[] TopK(double[] scores, int topK)
    {
        int n = scores.Length;
        int k = Math.Min(topK, n);
        if (k == 0) return Array.Empty<(int, double)>();

        // Min-heap of size K keeps the K largest seen so far.
        // Heap is stored in two parallel arrays; root at index 0.
        var heapIds = new int[k];
        var heapScores = new double[k];
        int heapSize = 0;

        for (int d = 0; d < n; d++)
        {
            double s = scores[d];
            if (double.IsNegativeInfinity(s)) continue;

            if (heapSize < k)
            {
                heapIds[heapSize] = d;
                heapScores[heapSize] = s;
                heapSize++;
                if (heapSize == k) BuildHeap(heapIds, heapScores, heapSize);
            }
            else if (s > heapScores[0])
            {
                heapIds[0] = d;
                heapScores[0] = s;
                SiftDown(heapIds, heapScores, 0, heapSize);
            }
        }

        // Drain into a sorted-descending array.
        var result = new (int, double)[heapSize];
        for (int i = heapSize - 1; i >= 0; i--)
        {
            result[i] = (heapIds[0], heapScores[0]);
            heapSize--;
            if (heapSize > 0)
            {
                heapIds[0] = heapIds[heapSize];
                heapScores[0] = heapScores[heapSize];
                SiftDown(heapIds, heapScores, 0, heapSize);
            }
        }
        return result;
    }

    private static void BuildHeap(int[] ids, double[] scores, int size)
    {
        for (int i = (size - 2) / 2; i >= 0; i--)
            SiftDown(ids, scores, i, size);
    }

    private static void SiftDown(int[] ids, double[] scores, int start, int size)
    {
        int i = start;
        while (true)
        {
            int left = 2 * i + 1;
            int right = left + 1;
            int smallest = i;

            if (left < size && scores[left] < scores[smallest]) smallest = left;
            if (right < size && scores[right] < scores[smallest]) smallest = right;
            if (smallest == i) return;

            (scores[i], scores[smallest]) = (scores[smallest], scores[i]);
            (ids[i], ids[smallest]) = (ids[smallest], ids[i]);
            i = smallest;
        }
    }
}
