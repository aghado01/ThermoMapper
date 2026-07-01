using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;

namespace Hashish;

/// <summary>
/// Trend Micro Locality Sensitive Hash (TLSH).
/// Produces a fuzzy digest of the form <c>T1{checksum}{lenCode}{body}</c>
/// suitable for near-duplicate detection via distance comparison.
/// </summary>
public static class TrendLocalitySensitiveHash
{
    private const int WindowSize = 5;
    private const int BucketCount = 256;
    private const int MinLength = 50;

    /// <summary>
    /// Computes a TLSH digest. Returns empty string if <paramref name="content"/>
    /// is null, empty, or shorter than <c>50</c> characters.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static string Compute(string content)
    {
        if (string.IsNullOrEmpty(content) || content.Length < MinLength)
            return string.Empty;

        ReadOnlySpan<char> span = content.AsSpan();

        // Bucket array on the stack — 256 × 4 = 1 KB, well within safe stackalloc range.
        Span<int> buckets = stackalloc int[BucketCount];
        ComputeBuckets(span, buckets);

        // Quartile scratch on the stack.
        Span<int> sorted = stackalloc int[BucketCount];
        buckets.CopyTo(sorted);
        MemoryExtensions.Sort(sorted);

        int q1 = sorted[BucketCount / 4];
        int q2 = sorted[BucketCount / 2];
        int q3 = sorted[BucketCount * 3 / 4];

        string body = BuildBody(buckets, q1, q2, q3);
        string checksum = ComputeChecksum(span);
        string lenCode = EncodeLength(content.Length);

        return $"T1{checksum}{lenCode}{body}";
    }

    /// <summary>
    /// Distance between two TLSH digests. Lower = more similar.
    /// Combines length, checksum, and body differences.
    /// </summary>
    public static int Compare(string hash1, string hash2)
    {
        if (!hash1.StartsWith("T1", StringComparison.Ordinal) ||
            !hash2.StartsWith("T1", StringComparison.Ordinal))
            throw new ArgumentException("Invalid TLSH format — expected 'T1' prefix.");

        ReadOnlySpan<char> s1 = hash1.AsSpan(2);
        ReadOnlySpan<char> s2 = hash2.AsSpan(2);

        // Layout: [2-char checksum][1-char len][body...]
        int checksumDist = Math.Abs(
            Convert.ToInt32(s1[..2].ToString(), 16) -
            Convert.ToInt32(s2[..2].ToString(), 16));

        int lenDist = Math.Abs(
            Convert.ToInt32(s1[2..3].ToString(), 16) -
            Convert.ToInt32(s2[2..3].ToString(), 16));

        ReadOnlySpan<char> body1 = s1[3..];
        ReadOnlySpan<char> body2 = s2[3..];

        int bodyDist = 0;
        int minLen = Math.Min(body1.Length, body2.Length);
        for (int i = 0; i < minLen; i++)
            if (body1[i] != body2[i]) bodyDist++;

        return (lenDist * 12) + checksumDist + bodyDist;
    }

    // Fill bucket histogram via sliding window over the content span.
    // ReadOnlySpan<char> slicing avoids Substring allocations; Pearson step is inlined.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void ComputeBuckets(ReadOnlySpan<char> content, Span<int> buckets)
    {
        int limit = content.Length - WindowSize;
        for (int i = 0; i <= limit; i++)
        {
            byte h = PearsonHash(content.Slice(i, WindowSize));
            buckets[h]++;
        }
    }

    // Pearson hash over a character window — iterates ReadOnlySpan<char> directly.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte PearsonHash(ReadOnlySpan<char> window)
    {
        byte hash = 0;
        for (int i = 0; i < window.Length; i++)
        {
            byte b = (byte)((int)window[i] & 0xFF);
            hash = (byte)(((hash ^ b) * 31) & 0xFF);
        }
        return hash;
    }

    // Encode bucket values as 2-bit codes packed into bytes, then hex-encode.
    // Direct bit-packing replaces BitArray — no object allocation, no CopyTo.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static string BuildBody(ReadOnlySpan<int> buckets, int q1, int q2, int q3)
    {
        // 256 buckets × 2 bits = 512 bits = 64 bytes.
        Span<byte> body = stackalloc byte[BucketCount * 2 / 8];

        for (int i = 0; i < BucketCount; i++)
        {
            int code = buckets[i] <= q1 ? 0
                     : buckets[i] <= q2 ? 1
                     : buckets[i] <= q3 ? 2
                     : 3;

            int byteIdx = i * 2 / 8;
            int bitShift = (i * 2) % 8;
            body[byteIdx] |= (byte)(code << bitShift);
        }

        return Convert.ToHexString(body).ToLowerInvariant();
    }

    // Checksum: sum all UTF-8 bytes mod 256, formatted as two hex digits.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static string ComputeChecksum(ReadOnlySpan<char> content)
    {
        int maxBytes = Encoding.UTF8.GetMaxByteCount(content.Length);
        int sum = 0;

        if (maxBytes <= 1024)
        {
            Span<byte> buf = stackalloc byte[1024];
            int written = Encoding.UTF8.GetBytes(content, buf);
            ReadOnlySpan<byte> bytes = buf[..written];
            for (int i = 0; i < bytes.Length; i++) sum += bytes[i];
        }
        else
        {
            byte[] buf = Encoding.UTF8.GetBytes(content.ToString());
            for (int i = 0; i < buf.Length; i++) sum += buf[i];
        }

        return (sum % 256).ToString("x2");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static string EncodeLength(int length)
    {
        if (length <= 0) return "0";
        int code = Math.Min((int)Math.Log(length, 2), 15);
        return code.ToString("x");
    }
}
