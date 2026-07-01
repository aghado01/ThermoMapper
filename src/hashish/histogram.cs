// ============================================================================
// Hashish/Histogram.cs
// ============================================================================
// Normalized probability mass functions from integer count data, with optional
// Lidstone (add-α) smoothing. Text-native BuildUnigram path maps a token
// sequence through a shared vocabulary into a PMF suitable for KLDivergence
// or FisherRaoSimplex.
//
// Smoothing:
//   alpha = 0.0  — no smoothing; zero counts stay zero
//   alpha = 0.5  — Jeffreys prior (recommended for KL against sparse Q)
//   alpha = 1.0  — Laplace smoothing
//
// The normalize path is zero-alloc: caller supplies the output span.
// Convenience allocating overloads are provided for interactive use.
// ============================================================================

#nullable enable
using System;
using System.Buffers;
using System.Collections.Frozen;
using System.Runtime.CompilerServices;

namespace Hashish;

public static class Histogram
{
    /// <summary>
    /// Converts integer counts to a normalized probability mass function
    /// with optional Lidstone (add-α) smoothing. Writes into <paramref name="output"/>.
    /// <para>
    /// Formula: output[i] = (counts[i] + α) / (Σⱼ counts[j] + α · |support|)
    /// </para>
    /// <para>
    /// With α = 0 and all-zero counts, output is all zeros (degenerate distribution).
    /// With α &gt; 0 and all-zero counts, output is uniform over the support.
    /// </para>
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static void Normalize(
        ReadOnlySpan<int> counts,
        Span<double> output,
        double alpha = 0.0)
    {
        if (output.Length != counts.Length)
            throw new ArgumentException(
                $"Output span length ({output.Length}) must equal counts length ({counts.Length}).");

        int n = counts.Length;
        long total = 0;
        for (int i = 0; i < n; i++)
            total += counts[i];

        double denom = total + alpha * n;
        if (denom <= 0.0)
        {
            output.Clear();
            return;
        }

        double inv = 1.0 / denom;
        for (int i = 0; i < n; i++)
            output[i] = (counts[i] + alpha) * inv;
    }

    /// <summary>Allocating overload of <see cref="Normalize(ReadOnlySpan{int},Span{double},double)"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] Normalize(ReadOnlySpan<int> counts, double alpha = 0.0)
    {
        var output = new double[counts.Length];
        Normalize(counts, output, alpha);
        return output;
    }

    /// <summary>
    /// Builds a unigram PMF from a token sequence against a shared vocabulary.
    /// Tokens not present in <paramref name="vocab"/> are silently skipped.
    /// <para>
    /// The output span must have length equal to <paramref name="vocab"/>.Count.
    /// Uses <see cref="ArrayPool{T}"/> for the intermediate count buffer.
    /// </para>
    /// </summary>
    /// <param name="tokens">Token sequence (e.g. from <see cref="Tokenizer"/>).</param>
    /// <param name="vocab">Token → index map; typically from <see cref="CooccurrenceModel.TokenIndex"/>
    /// or <see cref="TfIdfModel"/>.</param>
    /// <param name="output">Destination PMF span; length must equal vocab.Count.</param>
    /// <param name="alpha">Lidstone smoothing parameter.</param>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static void BuildUnigram(
        ReadOnlySpan<string> tokens,
        FrozenDictionary<string, int> vocab,
        Span<double> output,
        double alpha = 0.0)
    {
        ArgumentNullException.ThrowIfNull(vocab);
        int v = vocab.Count;
        if (output.Length != v)
            throw new ArgumentException(
                $"Output span length ({output.Length}) must equal vocab size ({v}).");

        int[] counts = ArrayPool<int>.Shared.Rent(v);
        counts.AsSpan(0, v).Clear();
        try
        {
            for (int i = 0; i < tokens.Length; i++)
            {
                if (vocab.TryGetValue(tokens[i], out int idx))
                    counts[idx]++;
            }
            Normalize(counts.AsSpan(0, v), output, alpha);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(counts);
        }
    }

    /// <summary>Allocating overload of
    /// <see cref="BuildUnigram(ReadOnlySpan{string},FrozenDictionary{string,int},Span{double},double)"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double[] BuildUnigram(
        ReadOnlySpan<string> tokens,
        FrozenDictionary<string, int> vocab,
        double alpha = 0.0)
    {
        var output = new double[vocab.Count];
        BuildUnigram(tokens, vocab, output, alpha);
        return output;
    }
}
