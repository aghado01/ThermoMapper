#nullable enable
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>
/// Word-level n-gram shingling for semantic text overlap. This is distinct from
/// MinHash's character shingles: word shingles preserve token boundaries and are
/// better suited to containment, deduplication, and topic-ish overlap checks.
/// </summary>
public static class WordShingler
{
    /// <summary>
    /// Builds ordered word shingles. Duplicate shingles are preserved.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static string[] Build(
        string text,
        int width = 2,
        bool ignoreCase = true,
        bool normalizeCompatibility = true,
        int minTokenLength = 1,
        string separator = " ")
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(separator);
        ArgumentOutOfRangeException.ThrowIfLessThan(width, 1, nameof(width));

        string[] tokens = TokenizerPreprocessing.TokenizeWords(
            text, ignoreCase, normalizeCompatibility, minTokenLength);

        if (tokens.Length < width)
            return Array.Empty<string>();

        int count = tokens.Length - width + 1;
        var shingles = new string[count];
        for (int i = 0; i < count; i++)
            shingles[i] = string.Join(separator, tokens, i, width);

        return shingles;
    }

    /// <summary>
    /// Builds a deduplicated word-shingle set.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static HashSet<string> BuildSet(
        string text,
        int width = 2,
        bool ignoreCase = true,
        bool normalizeCompatibility = true,
        int minTokenLength = 1,
        string separator = " ")
    {
        string[] shingles = Build(text, width, ignoreCase, normalizeCompatibility, minTokenLength, separator);
        return new HashSet<string>(shingles, StringComparer.Ordinal);
    }
}
