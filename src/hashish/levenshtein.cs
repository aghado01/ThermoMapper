// ============================================================================
// Hashish/Levenshtein.cs
// ============================================================================
// Character-level Levenshtein edit distance and normalized similarity.
//
// Uses the two-row DP optimization: only current and previous rows are kept,
// reducing memory from O(m×n) to O(min(m,n)). The shorter string is always
//
// Scratch rows come from ArrayPool for strings longer than the stack threshold,
// avoiding heap allocations in the common short-string case.
// ============================================================================

#nullable enable
using System;
using System.Buffers;
using System.Runtime.CompilerServices;

namespace Hashish;

public static class Levenshtein
{
    // Stack-allocate rows up to this many columns (both rows, so 2× on stack).
    private const int StackThreshold = 256;

    /// <summary>
    /// Character-level edit distance: minimum insertions, deletions, and
    /// substitutions to transform <paramref name="a"/> into <paramref name="b"/>.
    /// Returns 0 for equal strings; returns <c>max(a.Length, b.Length)</c> for
    /// a null/empty string paired with a non-empty one.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static int Distance(ReadOnlySpan<char> a, ReadOnlySpan<char> b)
    {
        // Trim common prefix and suffix — reduces DP table size for near-equal inputs.
        int start = 0;
        while (start < a.Length && start < b.Length && a[start] == b[start])
            start++;
        while (a.Length - start > 0 && b.Length - start > 0
               && a[a.Length - 1] == b[b.Length - 1])
        {
            a = a[..^1];
            b = b[..^1];
        }
        a = a[start..];
        b = b[start..];

        if (a.IsEmpty) return b.Length;
        if (b.IsEmpty) return a.Length;

        // Keep the shorter string on the column axis.
        if (a.Length < b.Length)
        {
            ReadOnlySpan<char> tmp = a;
            a = b;
            b = tmp;
        }

        int rows = a.Length;
        int cols = b.Length; // cols <= rows

        if (cols + 1 <= StackThreshold)
        {
            Span<int> prev = stackalloc int[cols + 1];
            Span<int> curr = stackalloc int[cols + 1];
            return TwoRowDp(a, b, prev, curr, rows, cols);
        }

        int[] rentedPrev = ArrayPool<int>.Shared.Rent(cols + 1);
        int[] rentedCurr = ArrayPool<int>.Shared.Rent(cols + 1);
        try
        {
            return TwoRowDp(a, b,
                rentedPrev.AsSpan(0, cols + 1),
                rentedCurr.AsSpan(0, cols + 1),
                rows, cols);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(rentedPrev);
            ArrayPool<int>.Shared.Return(rentedCurr);
        }
    }

    /// <inheritdoc cref="Distance(ReadOnlySpan{char}, ReadOnlySpan{char})"/>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Distance(string? a, string? b)
        => Distance(a.AsSpan(), b.AsSpan());

    /// <summary>
    /// Normalized edit similarity in [0, 1].
    /// <c>1.0</c> = identical; <c>0.0</c> = completely dissimilar.
    /// Defined as <c>1 - Distance(a, b) / max(a.Length, b.Length)</c>.
    /// Returns 1.0 if both inputs are empty.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Similarity(ReadOnlySpan<char> a, ReadOnlySpan<char> b)
    {
        int maxLen = Math.Max(a.Length, b.Length);
        if (maxLen == 0) return 1.0;
        return 1.0 - (double)Distance(a, b) / maxLen;
    }

    /// <inheritdoc cref="Similarity(ReadOnlySpan{char}, ReadOnlySpan{char})"/>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Similarity(string? a, string? b)
        => Similarity(a.AsSpan(), b.AsSpan());

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static int TwoRowDp(
        ReadOnlySpan<char> a, ReadOnlySpan<char> b,
        Span<int> prev, Span<int> curr,
        int rows, int cols)
    {
        // prev[j] = edit distance for a[0..0] vs b[0..j] (empty a prefix)
        for (int j = 0; j <= cols; j++)
            prev[j] = j;

        for (int i = 1; i <= rows; i++)
        {
            curr[0] = i;
            for (int j = 1; j <= cols; j++)
            {
                int sub = prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
                int del = prev[j] + 1;
                int ins = curr[j - 1] + 1;
                curr[j] = sub < del ? (sub < ins ? sub : ins)
                                    : (del < ins ? del : ins);
            }
            // Swap rows for next iteration.
            var tmp = prev; prev = curr; curr = tmp;
        }

        return prev[cols];
    }
}
