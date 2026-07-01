// ============================================================================
// Hashish/CooccurrenceStats.cs
// ============================================================================
// PMI, PPMI, contextual entropy, and conditional probability derived from
// a CooccurrenceModel. All methods are stateless given the model.
//
// PMI reference:
//   Church & Hanks (1990). "Word association norms, mutual information,
//   and lexicography." Computational Linguistics 16(1):22–29.
//
// PPMI (Positive PMI) is the standard variant used in distributional
// semantics — negative PMI values are clamped to 0, which avoids the
// instability of large negative values for rare co-occurrences.
//
// Contextual entropy measures how diverse a token's context distribution
// is. High entropy = appears in many different contexts (semantically
// broad or ambiguous). Low entropy = consistent context (specific role).
// ============================================================================

#nullable enable
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>
/// PMI and distributional statistics derived from a <see cref="CooccurrenceModel"/>.
/// </summary>
public static class CooccurrenceStats
{
    private const double Epsilon = 1e-10;

    // ── PMI ──────────────────────────────────────────────────────────────────

    /// <summary>
    /// Pointwise Mutual Information (natural log):
    /// PMI(a, b) = log( P(a,b) / (P(a) · P(b)) )
    /// </summary>
    /// <returns>
    /// PMI value. Positive = tokens co-occur more than chance.
    /// Negative = tokens co-occur less than chance.
    /// Returns <see cref="double.NegativeInfinity"/> if the pair never co-occurs.
    /// Returns 0 if either token has zero marginal.
    /// </returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Pmi(CooccurrenceModel model, int indexA, int indexB)
    {
        ArgumentNullException.ThrowIfNull(model);
        return ComputePmi(model, indexA, indexB, clampNegative: false);
    }

    /// <summary>String overload of <see cref="Pmi(CooccurrenceModel,int,int)"/>.</summary>
    public static double Pmi(CooccurrenceModel model, string tokenA, string tokenB)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(tokenA, out int ia)) return 0.0;
        if (!model.TokenIndex.TryGetValue(tokenB, out int ib)) return 0.0;
        return ComputePmi(model, ia, ib, clampNegative: false);
    }

    /// <summary>
    /// Positive PMI (PPMI): max(PMI(a,b), 0).
    /// The standard variant for distributional semantic vectors —
    /// negative values are unreliable for sparse co-occurrence data.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Ppmi(CooccurrenceModel model, int indexA, int indexB)
    {
        ArgumentNullException.ThrowIfNull(model);
        return ComputePmi(model, indexA, indexB, clampNegative: true);
    }

    /// <summary>String overload of <see cref="Ppmi(CooccurrenceModel,int,int)"/>.</summary>
    public static double Ppmi(CooccurrenceModel model, string tokenA, string tokenB)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(tokenA, out int ia)) return 0.0;
        if (!model.TokenIndex.TryGetValue(tokenB, out int ib)) return 0.0;
        return ComputePmi(model, ia, ib, clampNegative: true);
    }

    /// <summary>
    /// Builds a PPMI vector for a single token over the full vocabulary.
    /// Useful as a distributional semantic representation for downstream
    /// cosine similarity or SPC distance matrix construction.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double[] PpmiVector(CooccurrenceModel model, int tokenIndex)
    {
        ArgumentNullException.ThrowIfNull(model);
        int v = model.VocabularySize;
        var vec = new double[v];
        for (int j = 0; j < v; j++)
            vec[j] = ComputePmi(model, tokenIndex, j, clampNegative: true);
        return vec;
    }

    /// <summary>String overload of <see cref="PpmiVector(CooccurrenceModel,int)"/>.</summary>
    public static double[] PpmiVector(CooccurrenceModel model, string token)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(token, out int idx))
            return new double[model.VocabularySize];
        return PpmiVector(model, idx);
    }

    /// <summary>
    /// Builds PPMI vectors for the entire vocabulary.
    /// Returns a jagged array where result[i] is the PPMI vector for
    /// vocabulary token i. Suitable for direct use as a distance matrix
    /// input via cosine distance.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double[][] PpmiMatrix(CooccurrenceModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        int v = model.VocabularySize;
        var matrix = new double[v][];
        for (int i = 0; i < v; i++)
            matrix[i] = PpmiVector(model, i);
        return matrix;
    }

    // ── Conditional probability ───────────────────────────────────────────────

    /// <summary>
    /// P(b | a) = count(a,b) / marginal(a).
    /// Probability of seeing token b in the context of token a.
    /// Returns 0 if token a has zero marginal.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double ConditionalProbability(
        CooccurrenceModel model, int indexA, int indexB)
    {
        ArgumentNullException.ThrowIfNull(model);
        int margA = model.Marginals[indexA];
        if (margA == 0) return 0.0;
        return (double)model.Count(indexA, indexB) / margA;
    }

    /// <summary>String overload of
    /// <see cref="ConditionalProbability(CooccurrenceModel,int,int)"/>.</summary>
    public static double ConditionalProbability(
        CooccurrenceModel model, string given, string target)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(given, out int ia)) return 0.0;
        if (!model.TokenIndex.TryGetValue(target, out int ib)) return 0.0;
        return ConditionalProbability(model, ia, ib);
    }

    // ── Contextual entropy ────────────────────────────────────────────────────

    /// <summary>
    /// Shannon entropy of the context distribution for a token (base 2, bits).
    /// H(token) = -Σ_j P(j|token) · log₂ P(j|token)
    ///
    /// High entropy = appears in diverse contexts (semantically broad/ambiguous).
    /// Low entropy  = appears in consistent, specific contexts.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double ContextualEntropy(CooccurrenceModel model, int tokenIndex)
    {
        ArgumentNullException.ThrowIfNull(model);

        int marginal = model.Marginals[tokenIndex];
        if (marginal == 0) return 0.0;

        ReadOnlySpan<int> row = model.Row(tokenIndex);
        double invMarginal = 1.0 / marginal;
        double entropy = 0.0;

        for (int j = 0; j < row.Length; j++)
        {
            int c = row[j];
            if (c == 0) continue;
            double p = c * invMarginal;
            entropy -= p * Math.Log2(p);
        }

        return entropy;
    }

    /// <summary>String overload of
    /// <see cref="ContextualEntropy(CooccurrenceModel,int)"/>.</summary>
    public static double ContextualEntropy(CooccurrenceModel model, string token)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(token, out int idx)) return 0.0;
        return ContextualEntropy(model, idx);
    }

    /// <summary>
    /// Normalized contextual entropy: H(token) / log₂(VocabularySize).
    /// Range [0, 1]. Comparable across corpora with different vocabulary sizes.
    /// Returns 0 if vocabulary size is &lt; 2.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double NormalizedContextualEntropy(CooccurrenceModel model, int tokenIndex)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (model.VocabularySize < 2) return 0.0;
        double maxEntropy = Math.Log2(model.VocabularySize);
        return ContextualEntropy(model, tokenIndex) / maxEntropy;
    }

    /// <summary>String overload of
    /// <see cref="NormalizedContextualEntropy(CooccurrenceModel,int)"/>.</summary>
    public static double NormalizedContextualEntropy(CooccurrenceModel model, string token)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(token, out int idx)) return 0.0;
        return NormalizedContextualEntropy(model, idx);
    }

    // ── Top-N context neighbors ───────────────────────────────────────────────

    /// <summary>
    /// Returns the top-N tokens most strongly associated with
    /// <paramref name="tokenIndex"/> by PPMI, sorted descending.
    /// Useful for sanity-checking the model and for inspecting
    /// cluster-representative tokens.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static (string Token, double Ppmi)[] TopContextNeighbors(
        CooccurrenceModel model, int tokenIndex, int topN = 10)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentOutOfRangeException.ThrowIfLessThan(topN, 1, nameof(topN));

        int v = model.VocabularySize;
        var heap = new List<(string Token, double Ppmi)>(topN + 1);

        for (int j = 0; j < v; j++)
        {
            if (j == tokenIndex) continue;
            double ppmi = ComputePmi(model, tokenIndex, j, clampNegative: true);
            if (ppmi <= 0.0) continue;

            heap.Add((model.Vocabulary[j], ppmi));
            heap.Sort(static (a, b) => a.Ppmi.CompareTo(b.Ppmi)); // ascending
            if (heap.Count > topN) heap.RemoveAt(0);
        }

        heap.Sort(static (a, b) => b.Ppmi.CompareTo(a.Ppmi)); // descending
        return heap.ToArray();
    }

    /// <summary>String overload of
    /// <see cref="TopContextNeighbors(CooccurrenceModel,int,int)"/>.</summary>
    public static (string Token, double Ppmi)[] TopContextNeighbors(
        CooccurrenceModel model, string token, int topN = 10)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.TokenIndex.TryGetValue(token, out int idx))
            return Array.Empty<(string, double)>();
        return TopContextNeighbors(model, idx, topN);
    }

    // ── Private core ──────────────────────────────────────────────────────────

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double ComputePmi(
        CooccurrenceModel model, int ia, int ib, bool clampNegative)
    {
        long total = model.TotalCount;
        if (total == 0) return 0.0;

        int countAB = model.Count(ia, ib);
        if (countAB == 0) return clampNegative ? 0.0 : double.NegativeInfinity;

        int margA = model.Marginals[ia];
        int margB = model.Marginals[ib];
        if (margA == 0 || margB == 0) return 0.0;

        // PMI = log( P(a,b) / (P(a) * P(b)) )
        //     = log( (count_ab / total) / ((marg_a / total) * (marg_b / total)) )
        //     = log( count_ab * total / (marg_a * marg_b) )
        double pmi = Math.Log(
            ((double)countAB * total) /
            ((double)margA * margB + Epsilon));

        return clampNegative ? Math.Max(0.0, pmi) : pmi;
    }
}
