#nullable enable
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Hashish;

/// <summary>Term-frequency weighting variant applied before IDF multiplication.</summary>
public enum TfVariant
{
    /// <summary>Raw integer term counts.</summary>
    Raw,

    /// <summary>1 + ln(tf) for tf &gt; 0, else 0. Matches sklearn's sublinear_tf.</summary>
    Sublinear,
}

/// <summary>Knobs controlling tokenization, vocabulary pruning, and TF-IDF weighting.</summary>
public sealed class TfIdfOptions
{
    /// <summary>Fold tokens to lowercase before counting.</summary>
    public bool IgnoreCase { get; init; } = true;

    /// <summary>Apply Unicode compatibility (FormKC) normalization to source text.</summary>
    public bool NormalizeCompatibility { get; init; } = true;

    /// <summary>Minimum token length kept after regex match. Tokens shorter than this are dropped.</summary>
    public int MinTokenLength { get; init; } = 1;

    /// <summary>Drop tokens with document frequency below this absolute count. 1 = no floor.</summary>
    public int MinDocFrequency { get; init; } = 1;

    /// <summary>Drop tokens whose document frequency ratio exceeds this. 1.0 = no ceiling.</summary>
    public double MaxDocFrequencyRatio { get; init; } = 1.0;

    /// <summary>TF weighting applied before IDF multiplication.</summary>
    public TfVariant TfVariant { get; init; } = TfVariant.Sublinear;

    /// <summary>IDF formula. Smoothed Lucene/BM25+ variant is the default.</summary>
    public IdfFormula IdfFormula { get; init; } = IdfFormula.Smooth;

    /// <summary>L2-normalize the final row vector. Required for cosine-as-dot.</summary>
    public bool L2Normalize { get; init; } = true;

    /// <summary>Parallelize batch transforms across documents. Default true.</summary>
    public bool Parallel { get; init; } = true;

    internal static TfIdfOptions Effective(TfIdfOptions? options) => options ?? new TfIdfOptions();
}

/// <summary>
/// Tokenized corpus intermediate. Holds the per-document token arrays produced once
/// by <see cref="TfIdf.Tokenize"/> and reused across multiple fits/transforms that
/// share tokenization settings. Skips regex/normalize re-walks when sweeping over
/// different DF cutoffs or TF variants.
/// </summary>
public sealed class TokenizedCorpus
{
    internal TokenizedCorpus(
        string[][] tokens,
        double averageDocumentLength,
        bool ignoreCase,
        bool normalizeCompatibility,
        int minTokenLength)
    {
        Tokens = tokens;
        AverageDocumentLength = averageDocumentLength;
        IgnoreCase = ignoreCase;
        NormalizeCompatibility = normalizeCompatibility;
        MinTokenLength = minTokenLength;
    }

    /// <summary>Per-document normalized token arrays. Read-only after construction.</summary>
    public string[][] Tokens { get; }

    /// <summary>Mean token count across all documents.</summary>
    public double AverageDocumentLength { get; }

    /// <summary>Number of documents.</summary>
    public int DocumentCount => Tokens.Length;

    internal bool IgnoreCase { get; }
    internal bool NormalizeCompatibility { get; }
    internal int MinTokenLength { get; }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void EnsureCompatible(TfIdfOptions options)
    {
        if (options.IgnoreCase != IgnoreCase
            || options.NormalizeCompatibility != NormalizeCompatibility
            || options.MinTokenLength != MinTokenLength)
        {
            throw new InvalidOperationException(
                "TfIdfOptions tokenization knobs disagree with the cached TokenizedCorpus. " +
                "Either re-tokenize via TfIdf.Tokenize or align the options.");
        }
    }
}

/// <summary>Immutable TF-IDF model: vocabulary, IDF weights, and the options used to fit.</summary>
public sealed class TfIdfModel
{
    internal TfIdfModel(
        FrozenDictionary<string, int> vocabulary,
        double[] idf,
        int documentCount,
        double averageDocumentLength,
        TfIdfOptions options)
    {
        Vocabulary = vocabulary;
        InverseDocumentFrequency = idf;
        Dimension = idf.Length;
        DocumentCount = documentCount;
        AverageDocumentLength = averageDocumentLength;
        Options = options;
    }

    /// <summary>Token → vocabulary index. Indices are dense in [0, Dimension).</summary>
    public FrozenDictionary<string, int> Vocabulary { get; }

    /// <summary>IDF weight per vocabulary index. Parallel to <see cref="Vocabulary"/> values.</summary>
    public double[] InverseDocumentFrequency { get; }

    /// <summary>Vocabulary size — the dense-vector dimension produced by Transform.</summary>
    public int Dimension { get; }

    /// <summary>Document count seen at fit time.</summary>
    public int DocumentCount { get; }

    /// <summary>Mean document length (token count) seen at fit time.</summary>
    public double AverageDocumentLength { get; }

    /// <summary>Options used to fit this model. Transform inherits the same knobs.</summary>
    public TfIdfOptions Options { get; }

    /// <summary>Vocab index lookup, or -1 if the token isn't in vocabulary.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int GetIndex(string token)
        => Vocabulary.TryGetValue(token, out int index) ? index : -1;

    /// <summary>
    /// Transform a single document into a dense vector of length <see cref="Dimension"/>.
    /// </summary>
    public double[] Transform(string document)
    {
        ArgumentNullException.ThrowIfNull(document);
        var row = new double[Dimension];
        AccumulateRow(document, row);
        FinalizeRow(row);
        return row;
    }

    /// <summary>
    /// Transform from a pre-tokenized document.
    /// Tokens are looked up against the model vocabulary; out-of-vocabulary tokens are ignored.
    /// </summary>
    public double[] Transform(ReadOnlySpan<string> tokens)
    {
        var row = new double[Dimension];
        AccumulateRowFromTokens(tokens, row);
        FinalizeRow(row);
        return row;
    }

    /// <summary>
    /// Transform a single document into sparse (indices, values) form.
    /// Indices are sorted ascending. Out-of-vocabulary tokens are ignored.
    /// </summary>
    public (int[] Indices, double[] Values) TransformSparse(string document)
    {
        ArgumentNullException.ThrowIfNull(document);
        Dictionary<int, double> counts = AccumulateSparseCounts(document);
        return MaterializeSparse(counts);
    }

    /// <summary>
    /// Batch dense transform of a raw corpus. Returns a flat <c>double[N * Dimension]</c>
    /// where row <c>d</c> lives at offset <c>d * Dimension</c>. Parallelized when
    /// <see cref="TfIdfOptions.Parallel"/> is true.
    /// </summary>
    public double[] TransformAll(IReadOnlyList<string> documents)
    {
        ArgumentNullException.ThrowIfNull(documents);
        int n = documents.Count;
        long totalElements = (long)n * Dimension;
        if (totalElements > Array.MaxLength)
            throw new InvalidOperationException(
                $"Dense matrix {n}×{Dimension} = {totalElements} elements exceeds Array.MaxLength. Use sparse transform.");
        var rows = new double[totalElements];

        if (Options.Parallel)
        {
            Parallel.For(0, n, d =>
            {
                var row = rows.AsSpan(d * Dimension, Dimension);
                AccumulateRow(documents[d], row);
                FinalizeRow(row);
            });
        }
        else
        {
            for (int d = 0; d < n; d++)
            {
                var row = rows.AsSpan(d * Dimension, Dimension);
                AccumulateRow(documents[d], row);
                FinalizeRow(row);
            }
        }
        return rows;
    }

    /// <summary>
    /// Batch dense transform from a cached <see cref="TokenizedCorpus"/>.
    /// Skips re-tokenization — the fastest path for fit-then-transform on the same corpus.
    /// </summary>
    public double[] TransformAll(TokenizedCorpus corpus)
    {
        ArgumentNullException.ThrowIfNull(corpus);
        corpus.EnsureCompatible(Options);

        int n = corpus.DocumentCount;
        long totalElements = (long)n * Dimension;
        if (totalElements > Array.MaxLength)
            throw new InvalidOperationException(
                $"Dense matrix {n}×{Dimension} = {totalElements} elements exceeds Array.MaxLength. Use sparse transform.");
        var rows = new double[totalElements];
        string[][] tokens = corpus.Tokens;

        if (Options.Parallel)
        {
            Parallel.For(0, n, d =>
            {
                var row = rows.AsSpan(d * Dimension, Dimension);
                AccumulateRowFromTokens(tokens[d], row);
                FinalizeRow(row);
            });
        }
        else
        {
            for (int d = 0; d < n; d++)
            {
                var row = rows.AsSpan(d * Dimension, Dimension);
                AccumulateRowFromTokens(tokens[d], row);
                FinalizeRow(row);
            }
        }
        return rows;
    }

    // ── Internal accumulation paths ─────────────────────────────────────────

    private void AccumulateRow(string document, Span<double> row)
    {
        if (string.IsNullOrEmpty(document)) return;

        NormalizationForm form = Options.NormalizeCompatibility
            ? NormalizationForm.FormKC
            : NormalizationForm.FormC;

        string normalized = document.Normalize(form);
        if (Options.IgnoreCase) normalized = normalized.ToLowerInvariant();
        ReadOnlySpan<char> source = normalized.AsSpan();

        var vocabLookup = Vocabulary.GetAlternateLookup<ReadOnlySpan<char>>();
        int minLen = Options.MinTokenLength;

        foreach (ValueMatch m in TokenizerPreprocessing.WordRegex.EnumerateMatches(source))
        {
            if (m.Length < minLen) continue;
            ReadOnlySpan<char> span = source.Slice(m.Index, m.Length);
            if (vocabLookup.TryGetValue(span, out int index))
                row[index] += 1.0;
        }
    }

    private void AccumulateRowFromTokens(ReadOnlySpan<string> tokens, Span<double> row)
    {
        for (int t = 0; t < tokens.Length; t++)
        {
            if (Vocabulary.TryGetValue(tokens[t], out int index))
                row[index] += 1.0;
        }
    }

    private Dictionary<int, double> AccumulateSparseCounts(string document)
    {
        var counts = new Dictionary<int, double>(64);
        if (string.IsNullOrEmpty(document)) return counts;

        NormalizationForm form = Options.NormalizeCompatibility
            ? NormalizationForm.FormKC
            : NormalizationForm.FormC;

        string normalized = document.Normalize(form);
        if (Options.IgnoreCase) normalized = normalized.ToLowerInvariant();
        ReadOnlySpan<char> source = normalized.AsSpan();

        var vocabLookup = Vocabulary.GetAlternateLookup<ReadOnlySpan<char>>();
        int minLen = Options.MinTokenLength;

        foreach (ValueMatch m in TokenizerPreprocessing.WordRegex.EnumerateMatches(source))
        {
            if (m.Length < minLen) continue;
            ReadOnlySpan<char> span = source.Slice(m.Index, m.Length);
            if (!vocabLookup.TryGetValue(span, out int index)) continue;

            ref double slot = ref CollectionsMarshal.GetValueRefOrAddDefault(counts, index, out _);
            slot += 1.0;
        }
        return counts;
    }

    private (int[] Indices, double[] Values) MaterializeSparse(Dictionary<int, double> counts)
    {
        int nnz = counts.Count;
        if (nnz == 0) return (Array.Empty<int>(), Array.Empty<double>());

        var indices = new int[nnz];
        var values = new double[nnz];
        int k = 0;
        foreach (var kvp in counts)
        {
            indices[k] = kvp.Key;
            values[k] = kvp.Value;
            k++;
        }

        Array.Sort(indices, values);

        for (int i = 0; i < nnz; i++)
            values[i] = ApplyTfIdf(values[i], InverseDocumentFrequency[indices[i]]);

        if (Options.L2Normalize)
        {
            double norm = Math.Sqrt(TensorPrimitives.Dot<double>(values, values));
            if (norm > 1e-12)
                TensorPrimitives.Multiply<double>(values, 1.0 / norm, values);
        }

        return (indices, values);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private double ApplyTfIdf(double count, double idf)
        => Options.TfVariant switch
        {
            TfVariant.Raw => count * idf,
            TfVariant.Sublinear => count > 0 ? (1.0 + Math.Log(count)) * idf : 0.0,
            _ => count * idf,
        };

    private void FinalizeRow(Span<double> row)
    {
        switch (Options.TfVariant)
        {
            case TfVariant.Sublinear:
                for (int i = 0; i < row.Length; i++)
                {
                    if (row[i] > 0.0)
                        row[i] = (1.0 + Math.Log(row[i])) * InverseDocumentFrequency[i];
                }
                break;

            case TfVariant.Raw:
            default:
                // counts * idf — vectorized
                TensorPrimitives.Multiply<double>(row, InverseDocumentFrequency, row);
                break;
        }

        if (Options.L2Normalize)
        {
            double norm = Math.Sqrt(TensorPrimitives.Dot<double>(row, row));
            if (norm > 1e-12)
                TensorPrimitives.Multiply<double>(row, 1.0 / norm, row);
        }
    }
}

/// <summary>
/// Tokenization + Fit + FitTransform orchestrator. Separate from
/// <see cref="TfIdfModel"/> so static call sites read like a pipeline.
/// </summary>
public static class TfIdf
{
    /// <summary>
    /// Tokenize a corpus once, materializing a reusable <see cref="TokenizedCorpus"/>.
    /// Subsequent <see cref="Fit(TokenizedCorpus, TfIdfOptions?)"/> and
    /// <see cref="TfIdfModel.TransformAll(TokenizedCorpus)"/> calls skip regex/normalize work.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static TokenizedCorpus Tokenize(
        IReadOnlyList<string> documents,
        TfIdfOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(documents);
        var opts = TfIdfOptions.Effective(options);

        int n = documents.Count;
        var tokens = new string[n][];
        double totalTokens = 0.0;

        for (int d = 0; d < n; d++)
        {
            tokens[d] = TokenizerPreprocessing.TokenizeWords(
                documents[d] ?? string.Empty,
                opts.IgnoreCase,
                opts.NormalizeCompatibility,
                opts.MinTokenLength);
            totalTokens += tokens[d].Length;
        }

        return new TokenizedCorpus(
            tokens,
            n == 0 ? 0.0 : totalTokens / n,
            opts.IgnoreCase,
            opts.NormalizeCompatibility,
            opts.MinTokenLength);
    }

    /// <summary>
    /// Fit a TF-IDF model directly from a raw corpus. Uses the span-based DF builder
    /// in <see cref="InverseDocumentFrequency.Compute"/> — single tokenization pass,
    /// one string allocation per unique corpus token.
    /// </summary>
    public static TfIdfModel Fit(IReadOnlyList<string> documents, TfIdfOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(documents);
        var opts = TfIdfOptions.Effective(options);

        IdfModel idf = InverseDocumentFrequency.Compute(
            documents,
            opts.IdfFormula,
            opts.IgnoreCase,
            opts.NormalizeCompatibility,
            opts.MinTokenLength);

        return BuildModel(idf, opts);
    }

    /// <summary>
    /// Fit a TF-IDF model from a cached tokenized corpus. Reuses already-allocated
    /// token strings — best path when fitting multiple times with different
    /// DF cutoffs or TF variants over the same input.
    /// </summary>
    public static TfIdfModel Fit(TokenizedCorpus corpus, TfIdfOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(corpus);
        var opts = TfIdfOptions.Effective(options);
        corpus.EnsureCompatible(opts);

        int n = corpus.DocumentCount;
        var comparer = opts.IgnoreCase ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
        var df = new Dictionary<string, int>(1024, comparer);
        var seen = new HashSet<string>(256, comparer);

        foreach (string[] docTokens in corpus.Tokens)
        {
            seen.Clear();
            for (int t = 0; t < docTokens.Length; t++)
                seen.Add(docTokens[t]);

            foreach (string token in seen)
            {
                ref int count = ref CollectionsMarshal.GetValueRefOrAddDefault(df, token, out _);
                count++;
            }
        }

        var weights = new Dictionary<string, double>(df.Count, comparer);
        foreach (var kvp in df)
            weights[kvp.Key] = InverseDocumentFrequency.ComputeWeight(n, kvp.Value, opts.IdfFormula);

        var idfModel = new IdfModel(
            n,
            corpus.AverageDocumentLength,
            df.ToFrozenDictionary(comparer),
            weights.ToFrozenDictionary(comparer));

        return BuildModel(idfModel, opts);
    }

    /// <summary>
    /// Fit and produce the dense corpus matrix in one orchestrated pass.
    /// Returns the model plus a flat <c>double[N * Dimension]</c> row matrix
    /// suitable for <see cref="TfIdfSearch"/> queries or downstream pairwise work.
    /// </summary>
    public static (TfIdfModel Model, double[] DenseRows) FitTransform(
        IReadOnlyList<string> documents,
        TfIdfOptions? options = null)
    {
        var opts = TfIdfOptions.Effective(options);
        TokenizedCorpus corpus = Tokenize(documents, opts);
        TfIdfModel model = Fit(corpus, opts);
        double[] rows = model.TransformAll(corpus);
        return (model, rows);
    }

    // ── Vocabulary construction with DF pruning ─────────────────────────────

    private static TfIdfModel BuildModel(IdfModel idf, TfIdfOptions opts)
    {
        int n = idf.DocumentCount;
        int minDf = Math.Max(1, opts.MinDocFrequency);
        int maxDf = opts.MaxDocFrequencyRatio >= 1.0
            ? int.MaxValue
            : Math.Max(minDf, (int)Math.Floor(opts.MaxDocFrequencyRatio * n));

        // Collect surviving tokens deterministically. Sorted by token for cross-run
        // index stability — vocab order doesn't affect math but matters for debug
        // dumps, persisted artifacts, and reproducible neighbor lookups.
        var kept = new List<string>(idf.DocumentFrequency.Count);
        foreach (var kvp in idf.DocumentFrequency)
        {
            if (kvp.Value >= minDf && kvp.Value <= maxDf)
                kept.Add(kvp.Key);
        }

        var comparer = opts.IgnoreCase ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
        kept.Sort(comparer);

        var vocab = new Dictionary<string, int>(kept.Count, comparer);
        var idfArray = new double[kept.Count];
        for (int i = 0; i < kept.Count; i++)
        {
            string token = kept[i];
            vocab[token] = i;
            idfArray[i] = idf.Weights[token];
        }

        return new TfIdfModel(
            vocab.ToFrozenDictionary(comparer),
            idfArray,
            n,
            idf.AverageDocumentLength,
            opts);
    }
}
