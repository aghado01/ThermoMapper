#nullable enable
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;

namespace Hashish;

public enum IdfFormula { Smooth, RobertsonSparckJones, Plain }

/// <summary>Immutable document-frequency and IDF statistics for text corpora.</summary>
public sealed class IdfModel
{
    internal IdfModel(
        int documentCount,
        double averageDocumentLength,
        FrozenDictionary<string, int> documentFrequency,
        FrozenDictionary<string, double> weights)
    {
        DocumentCount = documentCount;
        AverageDocumentLength = averageDocumentLength;
        DocumentFrequency = documentFrequency;
        Weights = weights;
    }

    public int DocumentCount { get; }
    public double AverageDocumentLength { get; }
    public FrozenDictionary<string, int> DocumentFrequency { get; }
    public FrozenDictionary<string, double> Weights { get; }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double GetWeight(string token, double unknownWeight = 0.0)
        => Weights.TryGetValue(token, out double weight) ? weight : unknownWeight;
}

/// <summary>
/// General-purpose IDF primitive. Bm25Stats remains the SimHash-specific helper;
/// this type exposes reusable document-frequency state for vectorizers and sketches.
/// </summary>
public static class InverseDocumentFrequency
{
    /// <summary>
    /// (LastSeenDoc, Count) pair stored per token. LastSeenDoc dedupes within a single
    /// document without a per-doc HashSet — incrementing only when the current doc index
    /// differs from the token's last-seen doc index. One string allocation per unique
    /// corpus token via .NET alternate-lookup; matches re-use spans into the normalized
    /// document text.
    /// </summary>
    private readonly record struct DocFreqEntry(int LastSeenDoc, int Count);

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static IdfModel Compute(
        IReadOnlyList<string> documents,
        IdfFormula formula = IdfFormula.Smooth,
        bool ignoreCase = true,
        bool normalizeCompatibility = true,
        int minTokenLength = 1)
    {
        ArgumentNullException.ThrowIfNull(documents);
        ArgumentOutOfRangeException.ThrowIfLessThan(minTokenLength, 1, nameof(minTokenLength));

        int docCount = documents.Count;
        if (docCount == 0)
        {
            return new IdfModel(
                0,
                0.0,
                FrozenDictionary<string, int>.Empty,
                FrozenDictionary<string, double>.Empty);
        }

        var comparer = ignoreCase ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
        var entries = new Dictionary<string, DocFreqEntry>(1024, comparer);
        var lookup = entries.GetAlternateLookup<ReadOnlySpan<char>>();

        NormalizationForm normForm = normalizeCompatibility
            ? NormalizationForm.FormKC
            : NormalizationForm.FormC;

        double totalTokens = 0.0;

        for (int d = 0; d < docCount; d++)
        {
            string raw = documents[d];
            if (string.IsNullOrEmpty(raw)) continue;

            string normalized = raw.Normalize(normForm);
            if (ignoreCase) normalized = normalized.ToLowerInvariant();
            ReadOnlySpan<char> source = normalized.AsSpan();

            int docTokens = 0;
            foreach (ValueMatch m in TokenizerPreprocessing.WordRegex.EnumerateMatches(source))
            {
                if (m.Length < minTokenLength) continue;
                docTokens++;
                ReadOnlySpan<char> tokenSpan = source.Slice(m.Index, m.Length);

                if (lookup.TryGetValue(tokenSpan, out DocFreqEntry entry))
                {
                    if (entry.LastSeenDoc != d)
                        lookup[tokenSpan] = new DocFreqEntry(d, entry.Count + 1);
                }
                else
                {
                    lookup[tokenSpan] = new DocFreqEntry(d, 1);
                }
            }
            totalTokens += docTokens;
        }

        var df = new Dictionary<string, int>(entries.Count, comparer);
        var weights = new Dictionary<string, double>(entries.Count, comparer);
        foreach (var kvp in entries)
        {
            int count = kvp.Value.Count;
            df[kvp.Key] = count;
            weights[kvp.Key] = ComputeWeight(docCount, count, formula);
        }

        return new IdfModel(
            docCount,
            totalTokens / docCount,
            df.ToFrozenDictionary(comparer),
            weights.ToFrozenDictionary(comparer));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double ComputeWeight(int documentCount, int documentFrequency, IdfFormula formula = IdfFormula.Smooth)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(documentCount, 0, nameof(documentCount));
        ArgumentOutOfRangeException.ThrowIfLessThan(documentFrequency, 0, nameof(documentFrequency));

        if (documentCount == 0 || documentFrequency == 0)
            return 0.0;

        return formula switch
        {
            IdfFormula.Smooth => Math.Log((documentCount + 1.0) / (documentFrequency + 1.0)) + 1.0,
            IdfFormula.RobertsonSparckJones => Math.Log(1.0 + (documentCount - documentFrequency + 0.5) / (documentFrequency + 0.5)),
            IdfFormula.Plain => Math.Log((double)documentCount / documentFrequency),
            _ => throw new NotSupportedException($"IDF formula {formula} is not supported.")
        };
    }
}
