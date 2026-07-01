// ============================================================================
// Hashish/CooccurrenceMatrix.cs
// ============================================================================
// Windowed token co-occurrence matrix builder for text corpora.
//
// Builds an immutable CooccurrenceModel from a corpus of tokenized documents
// using a symmetric sliding window. The model exposes raw co-occurrence counts
// and marginal frequencies needed by CooccurrenceStats (PMI, entropy, etc.).
//
// Design notes:
//   - Vocabulary is built in a single pass; matrix is accumulated in a second
//     pass. Both passes use the same TokenizerPreprocessing path as the rest
//     of Hashish so tokenization is consistent across the module.
//   - Counts are stored as a flat int[] in row-major order (vocab × vocab).
//     For large vocabularies this becomes memory-intensive; use maxVocabSize
//     to cap at the N most frequent tokens (frequency-ranked pruning).
//   - FrozenDictionary for the token→index map: O(1) lookup after build,
//     consistent with IdfModel's pattern.
//   - Thread safety: none on the builder; CooccurrenceModel is immutable.
// ============================================================================

#nullable enable
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Hashish;

/// <summary>
/// Immutable co-occurrence model produced by <see cref="CooccurrenceBuilder"/>.
/// Stores raw symmetric counts and marginal totals.
/// </summary>
public sealed class CooccurrenceModel
{
    private readonly int[] _counts;   // flat row-major vocab×vocab

    internal CooccurrenceModel(
        FrozenDictionary<string, int> tokenIndex,
        string[] vocabulary,
        int[] counts,
        int[] marginals,
        long totalCount,
        int windowSize)
    {
        TokenIndex = tokenIndex;
        Vocabulary = vocabulary;
        Marginals = marginals;
        TotalCount = totalCount;
        WindowSize = windowSize;
        _counts = counts;
    }

    /// <summary>Token → vocabulary index mapping.</summary>
    public FrozenDictionary<string, int> TokenIndex { get; }

    /// <summary>Vocabulary array; index matches <see cref="TokenIndex"/>.</summary>
    public string[] Vocabulary { get; }

    /// <summary>
    /// Per-token marginal co-occurrence count: sum of all co-occurrences
    /// involving token i across all context positions.
    /// </summary>
    public int[] Marginals { get; }

    /// <summary>Sum of all entries in the co-occurrence matrix.</summary>
    public long TotalCount { get; }

    /// <summary>Window radius used during construction.</summary>
    public int WindowSize { get; }

    /// <summary>Number of unique tokens in the vocabulary.</summary>
    public int VocabularySize => Vocabulary.Length;

    /// <summary>
    /// Raw co-occurrence count for the ordered pair (tokenA, tokenB).
    /// The matrix is symmetric: Count(a,b) == Count(b,a).
    /// Returns 0 if either token is outside the vocabulary.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Count(int indexA, int indexB)
        => _counts[indexA * Vocabulary.Length + indexB];

    /// <summary>
    /// Raw co-occurrence count by token string.
    /// Returns 0 if either token is outside the vocabulary.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Count(string tokenA, string tokenB)
    {
        if (!TokenIndex.TryGetValue(tokenA, out int ia)) return 0;
        if (!TokenIndex.TryGetValue(tokenB, out int ib)) return 0;
        return Count(ia, ib);
    }

    /// <summary>
    /// Returns the full co-occurrence row for <paramref name="tokenIndex"/>
    /// as a read-only span over the internal flat array.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<int> Row(int tokenIndex)
        => _counts.AsSpan(tokenIndex * Vocabulary.Length, Vocabulary.Length);

    /// <summary>
    /// Returns the full co-occurrence row for a token string.
    /// Returns an empty span if the token is not in the vocabulary.
    /// </summary>
    public ReadOnlySpan<int> Row(string token)
        => TokenIndex.TryGetValue(token, out int idx) ? Row(idx) : ReadOnlySpan<int>.Empty;
}

/// <summary>
/// Builds a <see cref="CooccurrenceModel"/> from a text corpus using a
/// symmetric sliding window over tokenized documents.
/// </summary>
public static class CooccurrenceBuilder
{
    /// <summary>
    /// Builds a co-occurrence model from a corpus of documents.
    /// </summary>
    /// <param name="documents">Raw text documents. Null entries are skipped.</param>
    /// <param name="windowSize">
    /// Symmetric context window radius. Each token co-occurs with up to
    /// <paramref name="windowSize"/> tokens on either side.
    /// </param>
    /// <param name="maxVocabSize">
    /// If &gt; 0, vocabulary is pruned to the <paramref name="maxVocabSize"/>
    /// most frequent tokens before accumulating counts. 0 = no pruning.
    /// </param>
    /// <param name="ignoreCase">Fold tokens to lowercase before indexing.</param>
    /// <param name="normalizeCompatibility">Unicode compatibility normalization.</param>
    /// <param name="minTokenLength">Minimum token length to include.</param>
    /// <param name="minTokenFrequency">
    /// Tokens appearing in fewer than this many positions across the corpus
    /// are excluded from the vocabulary regardless of <paramref name="maxVocabSize"/>.
    /// </param>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static CooccurrenceModel Build(
        IReadOnlyList<string?> documents,
        int windowSize = 5,
        int maxVocabSize = 0,
        bool ignoreCase = true,
        bool normalizeCompatibility = true,
        int minTokenLength = 1,
        int minTokenFrequency = 2)
    {
        ArgumentNullException.ThrowIfNull(documents);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1, nameof(windowSize));

        // ── Pass 1: tokenize corpus, build raw frequency table ────────────────
        var tokenFreq = new Dictionary<string, int>(4096,
            ignoreCase ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal);

        // Store tokenized documents to avoid re-tokenizing in pass 2.
        var tokenized = new string[documents.Count][];

        for (int d = 0; d < documents.Count; d++)
        {
            string? doc = documents[d];
            if (doc is null) { tokenized[d] = Array.Empty<string>(); continue; }

            string[] tokens = TokenizerPreprocessing.TokenizeWords(
                doc, ignoreCase, normalizeCompatibility, minTokenLength);

            tokenized[d] = tokens;
            foreach (string t in tokens)
            {
                ref int freq = ref CollectionsMarshal.GetValueRefOrAddDefault(tokenFreq, t, out _);
                freq++;
            }
        }

        // ── Build vocabulary: prune by minTokenFrequency, then by maxVocabSize ─
        var vocabList = new List<KeyValuePair<string, int>>(tokenFreq.Count);
        foreach (var pair in tokenFreq)
            if (pair.Value >= minTokenFrequency)
                vocabList.Add(pair);

        // Sort descending by frequency for deterministic pruning.
        vocabList.Sort(static (a, b) => b.Value.CompareTo(a.Value));

        if (maxVocabSize > 0 && vocabList.Count > maxVocabSize)
            vocabList.RemoveRange(maxVocabSize, vocabList.Count - maxVocabSize);

        int vocabSize = vocabList.Count;

        if (vocabSize == 0)
            return new CooccurrenceModel(
                FrozenDictionary<string, int>.Empty,
                Array.Empty<string>(),
                Array.Empty<int>(),
                Array.Empty<int>(),
                0L,
                windowSize);

        // Build token→index map and vocabulary array.
        var comparer = ignoreCase ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
        var indexMap = new Dictionary<string, int>(vocabSize, comparer);
        var vocabArr = new string[vocabSize];

        for (int i = 0; i < vocabSize; i++)
        {
            vocabArr[i] = vocabList[i].Key;
            indexMap[vocabList[i].Key] = i;
        }

        // ── Pass 2: accumulate co-occurrence counts ───────────────────────────
        var counts = new int[vocabSize * vocabSize];
        var marginals = new int[vocabSize];
        long total = 0L;

        for (int d = 0; d < tokenized.Length; d++)
        {
            string[] tokens = tokenized[d];
            int len = tokens.Length;

            for (int center = 0; center < len; center++)
            {
                if (!indexMap.TryGetValue(tokens[center], out int ci)) continue;

                int start = Math.Max(0, center - windowSize);
                int end = Math.Min(len - 1, center + windowSize);

                for (int ctx = start; ctx <= end; ctx++)
                {
                    if (ctx == center) continue;
                    if (!indexMap.TryGetValue(tokens[ctx], out int xi)) continue;

                    // Symmetric increment: both (ci,xi) and (xi,ci).
                    counts[ci * vocabSize + xi]++;
                    counts[xi * vocabSize + ci]++;
                    marginals[ci]++;
                    marginals[xi]++;
                    total += 2;
                }
            }
        }

        return new CooccurrenceModel(
            indexMap.ToFrozenDictionary(comparer),
            vocabArr,
            counts,
            marginals,
            total,
            windowSize);
    }
}
