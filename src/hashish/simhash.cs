#nullable enable
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;

namespace Hashish;

/// <summary>
/// BM25-weighted 64-bit SimHash. Stateful instance; build once, call Compute many times.
/// Use <see cref="Bm25Stats.Compute"/> to produce the IDF map from a corpus.
/// </summary>
public sealed class SimHash
{
    private const ulong FnvPrime = 1099511628211UL;
    private const ulong FnvOffsetBasis = 14695981039346656037UL;

    private readonly FrozenDictionary<string, double> _idfMap;
    private readonly double _avgDocLength;
    private readonly double _k1;
    private readonly double _b;
    private readonly double _unknownIdf;
    private readonly double _minWeight;
    private readonly double _maxIdf;

    // Local copy of the word-token regex so this file compiles standalone
    // under Add-Type. Identical options to Bm25Stats.WordRegex.
    private static readonly Regex WordRegex = new(
        @"\w+",
        RegexOptions.Compiled | RegexOptions.CultureInvariant | RegexOptions.NonBacktracking
    );

    /// <param name="idfMap">Frozen IDF weights from <see cref="Bm25Stats.Compute"/>. Null = empty.</param>
    /// <param name="avgDocLength">Mean corpus document length. 0 = use per-document length.</param>
    /// <param name="k1">BM25 term-saturation factor. Must be >= 0.</param>
    /// <param name="b">BM25 field-length normalisation. Must be in [0, 1].</param>
    /// <param name="unknownIdf">IDF assigned to tokens absent from the map.</param>
    /// <param name="minWeight">Tokens with BM25 weight &lt;= this value are skipped.</param>
    /// <param name="maxIdf">IDF values are clamped to this ceiling.</param>
    public SimHash(
        FrozenDictionary<string, double>? idfMap = null,
        double avgDocLength = 0.0,
        double k1 = 1.5,
        double b = 0.75,
        double unknownIdf = 0.0,
        double minWeight = 1e-6,
        double maxIdf = double.PositiveInfinity)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k1, 0.0, nameof(k1));
        ArgumentOutOfRangeException.ThrowIfLessThan(b, 0.0, nameof(b));
        ArgumentOutOfRangeException.ThrowIfGreaterThan(b, 1.0, nameof(b));
        ArgumentOutOfRangeException.ThrowIfLessThan(minWeight, 0.0, nameof(minWeight));

        _idfMap = idfMap ?? FrozenDictionary<string, double>.Empty;
        _avgDocLength = avgDocLength;
        _k1 = k1;
        _b = b;
        _unknownIdf = unknownIdf;
        _minWeight = minWeight;
        _maxIdf = maxIdf;
    }

    /// <summary>Compute the SimHash of <paramref name="text"/>.</summary>
    /// <returns>64-bit SimHash signature; 0 for empty or whitespace-only input.</returns>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public ulong Compute(string text)
    {
        if (string.IsNullOrWhiteSpace(text)) return 0UL;

        var matches = WordRegex.Matches(text);
        if (matches.Count == 0) return 0UL;

        int docLength = matches.Count;
        double avgLen = _avgDocLength <= 0.0 ? docLength : _avgDocLength;

        var tfMap = new Dictionary<string, int>(capacity: docLength, comparer: StringComparer.Ordinal);

        foreach (System.Text.RegularExpressions.Match m in matches)
        {
            string token = m.Value.ToLowerInvariant();
            ref int count = ref CollectionsMarshal.GetValueRefOrAddDefault(tfMap, token, out _);
            count++;
        }

        // 64-double accumulator on the stack — zero GC overhead.
        Span<double> vector = stackalloc double[64];

        // Pre-compute BM25 invariants once per document.
        double k1Plus1 = _k1 + 1.0;
        double bm25DenomConst = _k1 * (1.0 - _b + _b * (docLength / avgLen));

        foreach (var kvp in tfMap)
        {
            double idf = _idfMap.TryGetValue(kvp.Key, out double mapped) ? mapped : _unknownIdf;

            if (idf < 0.0) idf = 0.0;
            else if (idf > _maxIdf) idf = _maxIdf;

            double weight = (idf * k1Plus1 * kvp.Value) / (bm25DenomConst + kvp.Value);
            if (weight <= _minWeight) continue;

            // FNV-1a over chars via ReadOnlySpan — no intermediate string allocation.
            ReadOnlySpan<char> span = kvp.Key.AsSpan();
            ulong tokenHash = FnvOffsetBasis;
            for (int j = 0; j < span.Length; j++)
            {
                tokenHash ^= span[j];
                tokenHash *= FnvPrime;
            }

            // Branchless weight accumulation: add or subtract per bit.
            for (int bit = 0; bit < 64; bit++)
            {
                if ((tokenHash & (1UL << bit)) != 0)
                    vector[bit] += weight;
                else
                    vector[bit] -= weight;
            }
        }

        ulong hash = 0UL;
        for (int bit = 0; bit < 64; bit++)
        {
            if (vector[bit] > 0.0)
                hash |= 1UL << bit;
        }

        return hash;
    }

    /// <summary>
    /// Convenience overload: zero-setup hash with an explicit IDF map.
    /// Equivalent to <c>new SimHash(idfMap, avgDocLength).Compute(text)</c>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong Compute(
        string text,
        FrozenDictionary<string, double> idfMap,
        double avgDocLength,
        double k1 = 1.5,
        double b = 0.75)
        => new SimHash(idfMap, avgDocLength, k1, b).Compute(text);

    /// <summary>
    /// Hamming distance between two SimHash signatures (0 = identical, 64 = opposite).
    /// Uses a single POPCNT hardware instruction on supported CPUs.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int HammingDistance(ulong a, ulong b)
        => System.Numerics.BitOperations.PopCount(a ^ b);
}
