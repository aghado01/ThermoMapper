#nullable enable
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Hashish;

/// <summary>
/// Context-Triggered Piecewise Hashing (CTPH / ssdeep-style fuzzy hash).
/// Splits content into variable-length chunks via a rolling FNV-1a trigger,
/// then hashes each chunk to produce a comparison-friendly digest.
/// </summary>
public static class ContextTriggeredPiecewiseHash
{
    private const int MinBlockSize = 3;
    private const int MaxBlockSize = 96;
    private const ulong FnvPrime = 1099511628211UL;
    private const ulong FnvOffset = 14695981039346656037UL;

    // Upper bound for stackalloc DP rows in LevenshteinDistance.
    private const int LevenshteinStackLimit = 512;

    /// <summary>
    /// Computes a CTPH digest of the form <c>blockSize:hash1:hash2</c>.
    /// </summary>
    /// <param name="content">Input text.</param>
    /// <param name="blockSize">
    /// Override automatic block-size selection. Must be >= <c>MinBlockSize</c> (3) if provided.
    /// </param>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static string Compute(string content, int blockSize = 0)
    {
        if (string.IsNullOrEmpty(content)) return "0::";

        int bs = blockSize > 0
            ? Math.Max(blockSize, MinBlockSize)
            : SelectBlockSize(content.Length);

        string hash1 = HashSequence(ComputeChunks(content, bs));
        string hash2 = HashSequence(ComputeChunks(content, bs * 2));

        return $"{bs}:{hash1}:{hash2}";
    }

    /// <summary>
    /// Similarity score in [0, 100] between two CTPH digests.
    /// Returns 0 if block sizes differ by more than 2x.
    /// </summary>
    public static double Compare(string hash1, string hash2)
    {
        ReadOnlySpan<char> s1 = hash1.AsSpan();
        ReadOnlySpan<char> s2 = hash2.AsSpan();

        int colon1a = s1.IndexOf(':');
        int colon1b = s1[(colon1a + 1)..].IndexOf(':') + colon1a + 1;
        int colon2a = s2.IndexOf(':');
        int colon2b = s2[(colon2a + 1)..].IndexOf(':') + colon2a + 1;

        if (colon1a < 0 || colon1b <= colon1a || colon2a < 0 || colon2b <= colon2a)
            throw new ArgumentException("Invalid CTPH hash format.");

        int bs1 = int.Parse(s1[..colon1a]);
        int bs2 = int.Parse(s2[..colon2a]);

        double ratio = (double)Math.Max(bs1, bs2) / Math.Min(bs1, bs2);
        if (ratio > 2.0) return 0.0;

        // Choose matching sequence portions based on block-size relationship.
        ReadOnlySpan<char> seq1 = bs1 == bs2
            ? s1[(colon1a + 1)..colon1b]
            : s1[(colon1b + 1)..];
        ReadOnlySpan<char> seq2 = bs1 == bs2
            ? s2[(colon2a + 1)..colon2b]
            : s2[(colon2b + 1)..];

        int distance = LevenshteinDistance(seq1, seq2);
        int maxLen = Math.Max(seq1.Length, seq2.Length);
        if (maxLen == 0) return 100.0;

        return Math.Max(0.0, (1.0 - (double)distance / maxLen) * 100.0);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int SelectBlockSize(int contentLength)
    {
        if (contentLength < 4096) return MinBlockSize;

        int candidate = Math.Max(
            MinBlockSize,
            (int)Math.Ceiling(Math.Log(contentLength / 64.0, 2))
        );
        return Math.Min(candidate, MaxBlockSize);
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static List<ulong> ComputeChunks(string content, int blockSize)
    {
        if (blockSize <= 1) throw new ArgumentOutOfRangeException(nameof(blockSize), "Must be > 1.");

        var chunks = new List<ulong>();
        ulong hash = FnvOffset;
        int windowStart = 0;
        ulong trigger = (ulong)(blockSize - 1);

        // Iterate via ReadOnlySpan — bounds checks eliminated by the JIT for span-length loops.
        ReadOnlySpan<char> span = content.AsSpan();

        for (int i = 0; i < span.Length; i++)
        {
            ulong charCode = span[i];
            hash ^= charCode;
            hash *= FnvPrime;

            if (hash % (ulong)blockSize == trigger)
            {
                chunks.Add(HashChunk(span[windowStart..(i + 1)]));
                windowStart = i + 1;
            }
        }

        if (windowStart < span.Length)
            chunks.Add(HashChunk(span[windowStart..]));

        return chunks;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong HashChunk(ReadOnlySpan<char> chunk)
    {
        ulong hash = FnvOffset;
        for (int i = 0; i < chunk.Length; i++)
        {
            hash ^= chunk[i];
            hash *= FnvPrime;
        }
        return hash;
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static string HashSequence(List<ulong> chunks)
    {
        if (chunks.Count == 0) return string.Empty;

        int take = Math.Min(chunks.Count, 64);

        // MemoryMarshal.Cast: zero-copy reinterpret List<ulong> backing store as byte span.
        Span<ulong> chunkSpan = CollectionsMarshal.AsSpan(chunks)[..take];
        ReadOnlySpan<byte> byteView = MemoryMarshal.Cast<ulong, byte>(chunkSpan);

        string base64 = Convert.ToBase64String(byteView);
        int len = Math.Min(64, base64.Length);
        return base64[..len];
    }

    /// <summary>
    /// Edit distance between two character spans.
    /// Uses stackalloc DP rows for inputs up to <c>LevenshteinStackLimit</c> chars;
    /// falls back to <see cref="ArrayPool{T}"/> for larger inputs.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static int LevenshteinDistance(ReadOnlySpan<char> a, ReadOnlySpan<char> b)
    {
        int n = a.Length, m = b.Length;
        if (n == 0) return m;
        if (m == 0) return n;

        // Stack path: stackalloc must not be inside a try block.
        if (m <= LevenshteinStackLimit)
        {
            Span<int> prev = stackalloc int[m + 1];
            Span<int> curr = stackalloc int[m + 1];
            return LevenshteinCore(a, b, prev, curr);
        }

        // Heap path: ArrayPool for large inputs.
        int[] rentedPrev = ArrayPool<int>.Shared.Rent(m + 1);
        int[] rentedCurr = ArrayPool<int>.Shared.Rent(m + 1);
        try
        {
            return LevenshteinCore(a, b, rentedPrev.AsSpan(0, m + 1), rentedCurr.AsSpan(0, m + 1));
        }
        finally
        {
            ArrayPool<int>.Shared.Return(rentedPrev);
            ArrayPool<int>.Shared.Return(rentedCurr);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static int LevenshteinCore(ReadOnlySpan<char> a, ReadOnlySpan<char> b, Span<int> prev, Span<int> curr)
    {
        int n = a.Length, m = b.Length;
        for (int j = 0; j <= m; j++) prev[j] = j;

        for (int i = 1; i <= n; i++)
        {
            curr[0] = i;
            char aChar = a[i - 1];

            for (int j = 1; j <= m; j++)
            {
                int cost = aChar == b[j - 1] ? 0 : 1;
                curr[j] = Math.Min(
                    Math.Min(prev[j] + 1, curr[j - 1] + 1),
                    prev[j - 1] + cost
                );
            }

            // Swap rows — no copy.
            Span<int> tmp = prev; prev = curr; curr = tmp;
        }

        return prev[m];
    }

    /// <summary>String overload for callers with allocated strings.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int LevenshteinDistance(string a, string b)
        => LevenshteinDistance(a.AsSpan(), b.AsSpan());
}
