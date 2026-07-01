#nullable enable
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>
/// Set-overlap primitives for exact Jaccard similarity and asymmetric containment.
/// Useful before approximate MinHash/LSH indexing, and as a calibration oracle for it.
/// </summary>
public static class JaccardContainment
{
    /// <summary>Jaccard similarity: |A intersect B| / |A union B|.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double Similarity<T>(
        IEnumerable<T> first,
        IEnumerable<T> second,
        IEqualityComparer<T>? comparer = null)
        where T : notnull
    {
        var left = ToSet(first, comparer);
        var right = ToSet(second, comparer);

        if (left.Count == 0 && right.Count == 0)
            return 1.0;

        int intersection = IntersectionCount(left, right);
        int union = left.Count + right.Count - intersection;
        return union == 0 ? 1.0 : (double)intersection / union;
    }

    /// <summary>Jaccard distance: 1 - Jaccard similarity.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Distance<T>(
        IEnumerable<T> first,
        IEnumerable<T> second,
        IEqualityComparer<T>? comparer = null)
        where T : notnull
        => 1.0 - Similarity(first, second, comparer);

    /// <summary>
    /// Asymmetric containment: |query intersect candidate| / |query|.
    /// Returns 1 for an empty query set, matching mathematical subset semantics.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double Containment<T>(
        IEnumerable<T> query,
        IEnumerable<T> candidate,
        IEqualityComparer<T>? comparer = null)
        where T : notnull
    {
        var left = ToSet(query, comparer);
        if (left.Count == 0)
            return 1.0;

        var right = ToSet(candidate, comparer);
        return (double)IntersectionCount(left, right) / left.Count;
    }

    /// <summary>Overlap coefficient: |A intersect B| / min(|A|, |B|).</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double OverlapCoefficient<T>(
        IEnumerable<T> first,
        IEnumerable<T> second,
        IEqualityComparer<T>? comparer = null)
        where T : notnull
    {
        var left = ToSet(first, comparer);
        var right = ToSet(second, comparer);
        int denominator = Math.Min(left.Count, right.Count);
        return denominator == 0 ? 1.0 : (double)IntersectionCount(left, right) / denominator;
    }

    /// <summary>Word-shingle containment from <paramref name="query"/> into <paramref name="candidate"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double WordShingleContainment(
        string query,
        string candidate,
        int width = 2,
        bool ignoreCase = true,
        bool normalizeCompatibility = true)
        => Containment(
            WordShingler.BuildSet(query, width, ignoreCase, normalizeCompatibility),
            WordShingler.BuildSet(candidate, width, ignoreCase, normalizeCompatibility),
            StringComparer.Ordinal);

    /// <summary>
    /// Sørensen–Dice similarity: 2|A ∩ B| / (|A| + |B|).
    /// Gives more weight to shared elements than Jaccard; always ≥ Jaccard similarity.
    /// Returns 1.0 if both sets are empty.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double DiceSimilarity<T>(
        IEnumerable<T> first,
        IEnumerable<T> second,
        IEqualityComparer<T>? comparer = null)
        where T : notnull
    {
        var left = ToSet(first, comparer);
        var right = ToSet(second, comparer);

        int denominator = left.Count + right.Count;
        if (denominator == 0) return 1.0;

        return (2.0 * IntersectionCount(left, right)) / denominator;
    }

    /// <summary>Sørensen–Dice distance: 1 - Dice similarity.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double DiceDistance<T>(
        IEnumerable<T> first,
        IEnumerable<T> second,
        IEqualityComparer<T>? comparer = null)
        where T : notnull
        => 1.0 - DiceSimilarity(first, second, comparer);

    /// <summary>Word-shingle Dice similarity from <paramref name="a"/> and <paramref name="b"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double WordShingleDiceSimilarity(
        string a,
        string b,
        int width = 2,
        bool ignoreCase = true,
        bool normalizeCompatibility = true)
        => DiceSimilarity(
            WordShingler.BuildSet(a, width, ignoreCase, normalizeCompatibility),
            WordShingler.BuildSet(b, width, ignoreCase, normalizeCompatibility),
            StringComparer.Ordinal);

    private static HashSet<T> ToSet<T>(IEnumerable<T> values, IEqualityComparer<T>? comparer)
        where T : notnull
    {
        ArgumentNullException.ThrowIfNull(values);

        if (values is HashSet<T> set && comparer == null)
            return new HashSet<T>(set, set.Comparer);

        return new HashSet<T>(values, comparer);
    }

    private static int IntersectionCount<T>(HashSet<T> left, HashSet<T> right)
        where T : notnull
    {
        if (right.Count < left.Count)
            (left, right) = (right, left);

        int count = 0;
        foreach (T item in left)
        {
            if (right.Contains(item))
                count++;
        }

        return count;
    }
}
